using System;
using System.Reflection;
using System.Runtime.CompilerServices;
using UnityEngine;
using UnityEngine.Rendering;
using UnityEngine.Rendering.Universal;
#if UNITY_6000_0_OR_NEWER
using UnityEngine.Rendering.RenderGraphModule;
#endif

namespace URPSSR
{
#if UNITY_2021
    public class ASPFullScreenRenderPass : ScriptableRenderPass
    {
        public static readonly int blitTexture = Shader.PropertyToID("_BlitTexture");
        public static readonly int blitScaleBias = Shader.PropertyToID("_BlitScaleBias");

        private Material m_Material;
        private int m_PassIndex;
        private bool m_CopyActiveColor;
        private bool m_BindDepthStencilAttachment;
        private RenderTexture m_CopiedColor;

        private static MaterialPropertyBlock s_SharedPropertyBlock = new MaterialPropertyBlock();

        public ASPFullScreenRenderPass(string passName)
        {
            profilingSampler = new ProfilingSampler(passName);
        }

        public void SetupMembers(Material material, int passIndex, bool copyActiveColor,
            bool bindDepthStencilAttachment)
        {
            m_Material = material;
            m_PassIndex = passIndex;
            m_CopyActiveColor = copyActiveColor;
            m_BindDepthStencilAttachment = bindDepthStencilAttachment;
        }

        public override void OnCameraSetup(CommandBuffer cmd, ref RenderingData renderingData)
        {
            RenderTextureDescriptor desc = new RenderTextureDescriptor(
                renderingData.cameraData.cameraTargetDescriptor.width,
                renderingData.cameraData.cameraTargetDescriptor.height);
                
            desc.colorFormat = RenderTextureFormat.ARGB32;
            desc.msaaSamples = 1;
            desc.depthBufferBits = (int)DepthBits.None;
                
            if (m_CopiedColor == null)
            {
                m_CopiedColor = RenderTexture.GetTemporary(desc);
            }
        }

        public void Dispose()
        {
        }
        
        public override void OnCameraCleanup(CommandBuffer cmd)
        {
            if (m_CopiedColor != null)
            {
                RenderTexture.ReleaseTemporary(m_CopiedColor);
                m_CopiedColor = null;
            }
        }
        
        public void ExecuteCopyColorPass(CommandBuffer cmd, RenderTargetIdentifier sourceTexture, RenderTargetIdentifier target)
        {
            Blit(cmd, sourceTexture, target);
        }

        private static void ExecuteMainPass(CommandBuffer cmd, RTHandle sourceTexture, Material material, int passIndex)
        {
            s_SharedPropertyBlock.Clear();
            if (sourceTexture != null)
                s_SharedPropertyBlock.SetTexture(blitTexture, sourceTexture);

            // We need to set the "_BlitScaleBias" uniform for user materials with shaders relying on core Blit.hlsl to work
            s_SharedPropertyBlock.SetVector(blitScaleBias, new Vector4(1, 1, 0, 0));

            cmd.DrawProcedural(Matrix4x4.identity, material, passIndex, MeshTopology.Triangles, 3, 1,
                s_SharedPropertyBlock);
        }
        
        private void DrawQuad(CommandBuffer cmd, Material material, int shaderPass)
        {
            cmd.DrawMesh(RenderingUtils.fullscreenMesh, Matrix4x4.identity, material, 0, shaderPass);   
        }

        public override void Execute(ScriptableRenderContext context, ref RenderingData renderingData)
        {
            if (m_Material == null)
                return;
                
            ref var cameraData = ref renderingData.cameraData;
            CommandBuffer cmd = CommandBufferPool.Get();
            using (new ProfilingScope(null, profilingSampler))
            {
                if (m_CopyActiveColor)
                {
                    CoreUtils.SetRenderTarget(cmd, m_CopiedColor);
                    ExecuteCopyColorPass(cmd, cameraData.renderer.cameraColorTarget, m_CopiedColor);
                }

                if (m_BindDepthStencilAttachment)
                {
                    CoreUtils.SetRenderTarget(cmd, cameraData.renderer.cameraColorTarget,
                        cameraData.renderer.cameraDepthTarget);
                }
                else
                {
                    CoreUtils.SetRenderTarget(cmd, cameraData.renderer.cameraColorTarget);
                }

                m_Material.SetVector("_BlitScaleBias", new Vector4(1,1,0,0));
                m_Material.SetTexture("_BaseMap", m_CopiedColor);
                DrawQuad(cmd, m_Material, 0);
            }

            context.ExecuteCommandBuffer(cmd);
            cmd.Clear();
            CommandBufferPool.Release(cmd);
        }
    }
#endif
#if UNITY_2022_1_OR_NEWER
    public class SSRCompositePass : ScriptableRenderPass
    {
        public static readonly int s_blitTexture = Shader.PropertyToID("_BlitTexture");
        public static readonly int s_prevResolveTexture = Shader.PropertyToID("_PrevSSRColor");
        public static readonly int s_resolveTexture = Shader.PropertyToID("_ResolveSSRColor");
        public static readonly int s_blitScaleBias = Shader.PropertyToID("_BlitScaleBias");
        private int _ssrThinGBufferId = Shader.PropertyToID("_SsrThinGBuffer");
        private Material m_Material;
        private int m_PassIndex;
        private bool m_CopyActiveColor;
        private bool m_BindDepthStencilAttachment;
        private RTHandle m_CopiedColor;
        private RTHandle m_SSRResolvedColor;
        private RTHandle m_PrevSSRColor;

        public bool UsePaddedTexture;
        private float _paddedScreenHeight;
        private float _paddedScreenWidth;
        private float _screenHeight;
        private float _screenWidth;
        private Vector2 _paddedScale;

        private float _ssrIntensity;
        private bool _deferred;
        private MethodInfo _isTemporalAAEnabledMethod;
        
        private static MaterialPropertyBlock s_SharedPropertyBlock = new MaterialPropertyBlock();

        public SSRCompositePass(string passName)
        {
            profilingSampler = new ProfilingSampler(passName);
        }

        public void SetupMembers(RenderingData renderingData, Material material, int passIndex, bool copyActiveColor,
            bool bindDepthStencilAttachment, bool isPadded, bool renderToCustomTarget, float ssrIntensity, bool isDeferred = true)
        {
            m_Material = material;
            m_PassIndex = passIndex;
            m_CopyActiveColor = copyActiveColor;
            m_BindDepthStencilAttachment = bindDepthStencilAttachment;
            UsePaddedTexture = isPadded;
            _ssrIntensity = ssrIntensity;
            _deferred = isDeferred;
            
            var desc = renderingData.cameraData.cameraTargetDescriptor;
            desc.msaaSamples = 1;
            desc.useMipMap = true;
            desc.autoGenerateMips = true;
            desc.mipCount = 8;
            desc.depthBufferBits = (int)DepthBits.None;
            RenderingUtils.ReAllocateIfNeeded(ref m_PrevSSRColor, desc, name: "_PrevSSRColor",
                filterMode: FilterMode.Trilinear);
        }

        public override void OnCameraSetup(CommandBuffer cmd, ref RenderingData renderingData)
        {
            // FullScreenPass manages its own RenderTarget.
            // ResetTarget here so that ScriptableRenderer's active attachement can be invalidated when processing this ScriptableRenderPass.
            ResetTarget();

            if (m_CopyActiveColor)
            {
                var desc = renderingData.cameraData.cameraTargetDescriptor;
                desc.msaaSamples = 1;
                desc.useMipMap = true;
                desc.autoGenerateMips = true;
                desc.mipCount = 8;
                desc.depthBufferBits = (int)DepthBits.None;
                RenderingUtils.ReAllocateIfNeeded(ref m_CopiedColor, desc, name: "_SSRMainColorPaddedMap",
                    filterMode: FilterMode.Trilinear);
            }

            var cameraColorDesc = renderingData.cameraData.cameraTargetDescriptor;
            cameraColorDesc.msaaSamples = 1;
            cameraColorDesc.depthBufferBits = (int)DepthBits.None;
            RenderingUtils.ReAllocateIfNeeded(ref m_SSRResolvedColor, cameraColorDesc, name: "_ResolveSSRColor");
        }

        public void Dispose()
        {
            m_CopiedColor?.Release();
            m_SSRResolvedColor?.Release();
        }

        private void ExecuteCopyColorPass(CommandBuffer cmd, RTHandle sourceTexture)
        {
            Blitter.BlitTexture(cmd, sourceTexture, Vector2.one, 0.0f, false);
        }


        private static void ExecuteMainPass(CommandBuffer cmd, RTHandle sourceTexture, Material material, int passIndex, Camera camera, float intensity)
        {
            SSRPassUtil.PrepareSSRProjectionMatrix(camera);
            material.SetFloat("_SSRIntensity", intensity);
            
            s_SharedPropertyBlock.Clear();
            if (sourceTexture != null)
                s_SharedPropertyBlock.SetTexture(s_blitTexture, sourceTexture);
            s_SharedPropertyBlock.SetVector(s_blitScaleBias, new Vector4(1, 1, 0, 0));
            
            cmd.DrawProcedural(Matrix4x4.identity, material, passIndex, MeshTopology.Triangles, 3, 1,
                s_SharedPropertyBlock);
        }

        public override void Execute(ScriptableRenderContext context, ref RenderingData renderingData)
        {
            ref var cameraData = ref renderingData.cameraData;
            if (cameraData.renderer.cameraColorTargetHandle == null)
                return;
            
            CommandBuffer cmd = CommandBufferPool.Get();
            using (new ProfilingScope(cmd, profilingSampler))
            {
                if (m_CopyActiveColor)
                {
                    CoreUtils.SetRenderTarget(cmd, m_CopiedColor);
                    ExecuteCopyColorPass(cmd, cameraData.renderer.cameraColorTargetHandle);
                }

                CoreUtils.SetRenderTarget(cmd, m_SSRResolvedColor);
                SSRPassUtil.PrepareSSRProjectionMatrix(cameraData.camera);
                m_Material.SetFloat("_SSRIntensity", _ssrIntensity);
                ExecuteMainPass(cmd, m_CopiedColor, m_Material, m_PassIndex, cameraData.camera, _ssrIntensity);
                
                
                CoreUtils.SetRenderTarget(cmd, cameraData.renderer.cameraColorTargetHandle);
                m_Material.SetFloat("_Influence", 1);
                ExecuteMainPass(cmd, m_SSRResolvedColor, m_Material, 4, cameraData.camera, _ssrIntensity);
            }
            
            context.ExecuteCommandBuffer(cmd);
            cmd.Clear();
            CommandBufferPool.Release(cmd);
        }

#if UNITY_6000_0_OR_NEWER
        private static void ExecuteCopyColorPass(RasterCommandBuffer cmd, CopyPassData data)
        {
            Blitter.BlitTexture(cmd, data.inputTexture, Vector2.one, 0.0f, false);
        }
        
        private static void ExecuteCompositePass(RasterCommandBuffer cmd, CompositePassData data)
        {
            data.material.SetFloat("_SSRIntensity", data.ssrIntensity);
            s_SharedPropertyBlock.Clear();
            if (data.inputTexture.IsValid())
                s_SharedPropertyBlock.SetTexture(s_blitTexture, data.inputTexture);
            if (data.resolveTexture.IsValid())
                data.material.SetTexture(s_resolveTexture, data.resolveTexture);

            s_SharedPropertyBlock.SetVector(s_blitScaleBias, new Vector4(1, 1, 0, 0));
            cmd.DrawProcedural(Matrix4x4.identity, data.material, data.passIndex, MeshTopology.Triangles, 3, 1,
                s_SharedPropertyBlock);
        }
        
        private static void ExecuteResolvePass(RasterCommandBuffer cmd, ResolvePassData data)
        {
            SSRPassUtil.PrepareSSRProjectionMatrix(data.camera);
            SSRPassUtil.SetupGbufferTextureToMaterial(data.material, data.gbuffer);
            
            data.material.SetFloat("_SSRIntensity", data.ssrIntensity);
            data.material.SetFloat("_Influence", data.influence);
            s_SharedPropertyBlock.Clear();
            if (data.inputTexture.IsValid())
                s_SharedPropertyBlock.SetTexture(s_blitTexture, data.inputTexture);
            
            if (data.prevResolveTexture.IsValid())
                s_SharedPropertyBlock.SetTexture(s_prevResolveTexture, data.prevResolveTexture);
            s_SharedPropertyBlock.SetVector(s_blitScaleBias, new Vector4(1, 1, 0, 0));
            cmd.DrawProcedural(Matrix4x4.identity, data.material, data.passIndex, MeshTopology.Triangles, 3, 1,
                s_SharedPropertyBlock);
        }
        
        private static void ExecuteMainPass(RasterCommandBuffer cmd, RTHandle sourceTexture, Material material, int passIndex, Camera camera, float intensity, TextureHandle[] gbuffer)
        {
            SSRPassUtil.PrepareSSRProjectionMatrix(camera);
            SSRPassUtil.SetupGbufferTextureToMaterial(material, gbuffer);
            
            material.SetFloat("_SSRIntensity", intensity);
            
            s_SharedPropertyBlock.Clear();
            if (sourceTexture != null)
                s_SharedPropertyBlock.SetTexture(s_blitTexture, sourceTexture);
            s_SharedPropertyBlock.SetVector(s_blitScaleBias, new Vector4(1, 1, 0, 0));
            cmd.DrawProcedural(Matrix4x4.identity, material, passIndex, MeshTopology.Triangles, 3, 1,
                s_SharedPropertyBlock);
        }
        
        private static void ExecuteMainPass(RasterCommandBuffer cmd, RTHandle sourceTexture, Material material, int passIndex)
        {
            s_SharedPropertyBlock.Clear();
            if (sourceTexture != null)
                s_SharedPropertyBlock.SetTexture(s_blitTexture, sourceTexture);
            s_SharedPropertyBlock.SetVector(s_blitScaleBias, new Vector4(1, 1, 0, 0));
            cmd.DrawProcedural(Matrix4x4.identity, material, passIndex, MeshTopology.Triangles, 3, 1,
                s_SharedPropertyBlock);
        }
        
        private void CalculatePaddedScale(UniversalCameraData cameraData)
        {
            _screenWidth = cameraData.cameraTargetDescriptor.width;
            _screenHeight = cameraData.cameraTargetDescriptor.height;

            if (UsePaddedTexture)
            {
                _paddedScreenWidth = Mathf.NextPowerOfTwo((int)_screenWidth);
                _paddedScreenHeight = Mathf.NextPowerOfTwo((int)_screenHeight);
            }
            else
            {
                _paddedScreenWidth = _screenWidth;
                _paddedScreenHeight = _screenHeight;
            }

            Vector2 screenResolution = new Vector2(_screenWidth, _screenHeight);
            Vector2 paddedResolution = new Vector2(_paddedScreenWidth, _paddedScreenHeight);
            _paddedScale = paddedResolution / screenResolution;
        }
        
        public override void RecordRenderGraph(RenderGraph renderGraph, ContextContainer frameData)
        {
            
            UniversalResourceData resourcesData = frameData.Get<UniversalResourceData>();
            UniversalCameraData cameraData = frameData.Get<UniversalCameraData>();

            if (_isTemporalAAEnabledMethod == null)
            {
                var universalCameraDataType = cameraData.GetType();
                _isTemporalAAEnabledMethod = universalCameraDataType?.GetMethod("IsTemporalAAEnabled", 
                    BindingFlags.NonPublic | BindingFlags.Instance);   
            }
            var isTemporalAAEnabled = _isTemporalAAEnabledMethod != null ? (bool)_isTemporalAAEnabledMethod.Invoke(cameraData, null) : false;
            TextureHandle source, destination;
            TextureHandle copyColorTexture;

            Debug.Assert(resourcesData.cameraColor.IsValid());
            CalculatePaddedScale(cameraData);
            var copyDesc = cameraData.cameraTargetDescriptor;
            copyDesc.depthBufferBits = (int)DepthBits.None;
            source = resourcesData.activeColorTexture;
            copyColorTexture = UniversalRenderer.CreateRenderGraphTexture(renderGraph, copyDesc, "_ssrCopyColor", true);
            
            using (var builder = renderGraph.AddRasterRenderPass<CopyPassData>("SSR Copy Color Full Screen", out var passData, profilingSampler))
            {
                passData.inputTexture = source;
                passData.destTexture = copyColorTexture;
  //              builder.AllowPassCulling(false);
                builder.UseTexture(passData.inputTexture, AccessFlags.Read);
                builder.SetRenderAttachment(passData.destTexture, 0, AccessFlags.Write);
                builder.SetRenderFunc((CopyPassData data, RasterGraphContext rgContext) =>
                {
                    ExecuteCopyColorPass(rgContext.cmd, data);
                });
            }
            
            var resolveSSRColorTexture = UniversalRenderer.CreateRenderGraphTexture(renderGraph, copyDesc, "_resolveSSRColor", true);
            var previusFrameSSRColorTexture = renderGraph.ImportTexture(m_PrevSSRColor);
            using (var builder = renderGraph.AddRasterRenderPass<ResolvePassData>(passName, out var passData, profilingSampler))
            {
                passData.material = m_Material;
                passData.passIndex = m_PassIndex;
                passData.camera = cameraData.camera;
                passData.gbuffer = resourcesData.gBuffer;
                passData.prevResolveTexture = previusFrameSSRColorTexture;
                passData.inputTexture = copyColorTexture;
                passData.destination = resolveSSRColorTexture;
                passData.ssrIntensity = _ssrIntensity;
                passData.influence = isTemporalAAEnabled ? 0.1f : 1f;
                SSRPassUtil.SetupRequirementTextureToRead(builder, isTemporalAAEnabled ? ScriptableRenderPassInput.Motion : ScriptableRenderPassInput.Depth, ref resourcesData);
                if (_deferred)
                {
                    SSRPassUtil.SetupGbufferTextureToRead(builder, ref resourcesData);
                }
                else
                {
                    builder.UseGlobalTexture(_ssrThinGBufferId, AccessFlags.Read);
                }
                builder.UseTexture(passData.inputTexture, AccessFlags.Read);
                builder.UseTexture(passData.prevResolveTexture , AccessFlags.Read);
                builder.AllowPassCulling(false);
                builder.SetRenderAttachment(passData.destination, 0, AccessFlags.Write);
                
                builder.SetRenderFunc((ResolvePassData data, RasterGraphContext rgContext) =>
                {
                    ExecuteResolvePass(rgContext.cmd, data);
                });
            }
            
           using (var builder = renderGraph.AddRasterRenderPass<CompositePassData>(passName, out var passData, profilingSampler))
           {
                passData.material = m_Material;
                passData.passIndex = 4;
                passData.inputTexture = copyColorTexture;
                passData.resolveTexture = resolveSSRColorTexture;
                passData.ssrIntensity = _ssrIntensity;
                passData.destTexture = resourcesData.activeColorTexture;
                builder.UseTexture(passData.inputTexture, AccessFlags.Read);
                builder.UseTexture(passData.resolveTexture, AccessFlags.Read);
                builder.SetRenderAttachment(resourcesData.activeColorTexture, 0, AccessFlags.Write);
                builder.SetRenderFunc((CompositePassData data, RasterGraphContext rgContext) =>
                {
                    Blitter.BlitTexture(rgContext.cmd, data.resolveTexture, Vector2.one, passData.material, data.passIndex);
                });
                
            }

           if (isTemporalAAEnabled)
           {
               using (var builder = renderGraph.AddRasterRenderPass<CopyPassData>("SSR Write To History", out var passData, profilingSampler))
               {
                   passData.inputTexture = resourcesData.activeColorTexture;
                   passData.destTexture = previusFrameSSRColorTexture;
                   builder.AllowPassCulling(false);
                   builder.UseTexture(passData.inputTexture, AccessFlags.Read);
                   builder.SetRenderAttachment(passData.destTexture, 0, AccessFlags.Write);
                   builder.SetRenderFunc((CopyPassData data, RasterGraphContext rgContext) =>
                   {
                       ExecuteCopyColorPass(rgContext.cmd, data);
                   });
               }
           }
        }

        private class CopyPassData
        {
            public TextureHandle inputTexture;
            public TextureHandle destTexture;
        }
        
        private class ResolvePassData
        {
            public TextureHandle[] gbuffer;
            public TextureHandle inputTexture;
            public TextureHandle destination;
            public TextureHandle prevResolveTexture;
            public float influence;
            public float ssrIntensity;
            public Camera camera;
            public Material material;
            public int passIndex;
        }

        private class CompositePassData
        {
            public TextureHandle inputTexture;
            public TextureHandle destTexture;
            public TextureHandle resolveTexture;
            public float ssrIntensity;
            public Material material;
            public int passIndex;
        }
#endif
    }
#endif
}