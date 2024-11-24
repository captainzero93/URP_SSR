using System;
using UnityEngine;
using UnityEngine.Rendering;
using UnityEngine.Rendering.Universal;

#if UNITY_6000_0_OR_NEWER
using UnityEngine.Rendering.RenderGraphModule;
#endif

namespace URPSSR
{
#if UNITY_2022_1_OR_NEWER
    public class SSRTracePass : ScriptableRenderPass
    {
        public static readonly int s_blitTexture = Shader.PropertyToID("_BlitTexture");
        public static readonly int s_blitScaleBias = Shader.PropertyToID("_BlitScaleBias");

        private Material _material;
        private int _passIndex;
        
        private String _targetName;
        private RTHandle _customTarget;
        private int _ssrReflectMapID = Shader.PropertyToID("_SSReflectedColorMap");
        public bool DebugReflectMap;

        public bool UsePaddedTexture;
        private float _paddedScreenHeight;
        private float _paddedScreenWidth;
        private float _screenHeight;
        private float _screenWidth;
        private Vector2 _paddedScale;
        private bool _deferred;
        

        private static MaterialPropertyBlock s_SharedPropertyBlock = new MaterialPropertyBlock();

        public SSRTracePass(string passName)
        {
            profilingSampler = new ProfilingSampler(passName);
        }

        public void SetupMembers(Material material, int passIndex, string renderTargetName = "", bool isPadded = false, bool isDeferred = true)
        {
            _material = material;
            _passIndex = passIndex;
            if (renderTargetName != string.Empty)
            {
                _targetName = renderTargetName;
            }
            else
            {
                _targetName = "_SSRFullScreen_Target";
            }

            UsePaddedTexture = isPadded;
            _deferred = isDeferred;
        }

        public override void OnCameraSetup(CommandBuffer cmd, ref RenderingData renderingData)
        {
            if (UsePaddedTexture)
            {
                _screenWidth = renderingData.cameraData.cameraTargetDescriptor.width * URPSSRSettings.GlobalResolutionScale;
                _screenHeight = renderingData.cameraData.cameraTargetDescriptor.height * URPSSRSettings.GlobalResolutionScale;
                _paddedScreenWidth = Mathf.NextPowerOfTwo((int)_screenWidth);
                _paddedScreenHeight = Mathf.NextPowerOfTwo((int)_screenHeight);
            }
            else
            {
                _screenWidth = renderingData.cameraData.cameraTargetDescriptor.width;
                _screenHeight = renderingData.cameraData.cameraTargetDescriptor.height;
                _paddedScreenWidth = _screenWidth * URPSSRSettings.GlobalResolutionScale;
                _paddedScreenHeight = _screenHeight * URPSSRSettings.GlobalResolutionScale;
            }

            Vector2 screenResolution = new Vector2(_screenWidth, _screenHeight);
            Shader.SetGlobalVector("_ScreenResolution", screenResolution);
            if (UsePaddedTexture)
            {
                Vector2 paddedResolution = new Vector2(_paddedScreenWidth, _paddedScreenHeight);
                _paddedScale = paddedResolution / screenResolution;
                Shader.SetGlobalVector("_PaddedResolution", paddedResolution);
                Shader.SetGlobalVector("_PaddedScale", _paddedScale);

                float cx = 1.0f / (512.0f * paddedResolution.x);
                float cy = 1.0f / (512.0f * paddedResolution.y);

                Shader.SetGlobalVector("crossEpsilon", new Vector2(cx, cy));
            }
            else
            {
                _paddedScale = Vector2.one;
                Shader.SetGlobalVector("_PaddedScale", _paddedScale);
            }
            

            var desc = renderingData.cameraData.cameraTargetDescriptor;
            var originalWidth = _screenWidth;
            var originalHeight = _screenHeight;
            desc.msaaSamples = 1;
            desc.colorFormat = RenderTextureFormat.ARGBFloat;
            desc.depthBufferBits = (int)DepthBits.None;
            desc.width = Mathf.CeilToInt(_paddedScreenWidth);
            desc.height = Mathf.CeilToInt(_paddedScreenHeight);
            RenderingUtils.ReAllocateIfNeeded(ref _customTarget, desc, name: "_mySSRTexture", filterMode: FilterMode.Point);
            ConfigureTarget(_customTarget);
            ConfigureClear(ClearFlag.Color, Color.clear);   
        }

        public void Dispose()
        {
            _customTarget?.Release();
        }

        private static void ExecuteMainPass(CommandBuffer cmd, RTHandle sourceTexture, Material material, int passIndex)
        {
            s_SharedPropertyBlock.Clear();
            if (sourceTexture != null)
                s_SharedPropertyBlock.SetTexture(s_blitTexture, sourceTexture);

            // We need to set the "_BlitScaleBias" uniform for user materials with shaders relying on core Blit.hlsl to work
            s_SharedPropertyBlock.SetVector(s_blitScaleBias, new Vector4(1, 1, 0, 0));

            cmd.DrawProcedural(Matrix4x4.identity, material, passIndex, MeshTopology.Triangles, 3, 1,
                s_SharedPropertyBlock);
        }

        private void SetupGlobalVariable(ref CameraData cameraData)
        {
            var SSR_ProjectionMatrix = GL.GetGPUProjectionMatrix(cameraData.camera.projectionMatrix, false);
            var HalfCameraSize = new Vector2(
                (int)(cameraData.camera.pixelWidth * 0.5f),
                (int)(cameraData.camera.pixelHeight * 0.5f));

            Matrix4x4 warpToScreenSpaceMatrix = Matrix4x4.identity;
            warpToScreenSpaceMatrix.m00 = HalfCameraSize.x;
            warpToScreenSpaceMatrix.m03 = HalfCameraSize.x;
            warpToScreenSpaceMatrix.m11 = HalfCameraSize.y;
            warpToScreenSpaceMatrix.m13 = HalfCameraSize.y;

            Matrix4x4 SSR_ProjectToPixelMatrix = warpToScreenSpaceMatrix * SSR_ProjectionMatrix;
            Shader.SetGlobalVector("_WorldSpaceViewDir", cameraData.camera.transform.forward);
            Shader.SetGlobalMatrix("_SSR_ProjectionMatrix", SSR_ProjectToPixelMatrix);
            Shader.SetGlobalVector("_SSR_ScreenSize", new Vector2(cameraData.camera.pixelWidth, cameraData.camera.pixelHeight));
        }
        
        private void SetupGlobalVariable(Matrix4x4 projectionMatrix, int pixelWidth, int pixelHeight, Vector3 camForward)
        {
            var SSR_ProjectionMatrix = GL.GetGPUProjectionMatrix(projectionMatrix, false);
            var HalfCameraSize = new Vector2(
                (int)(pixelWidth * 0.5f),
                (int)(pixelHeight * 0.5f));

            Matrix4x4 warpToScreenSpaceMatrix = Matrix4x4.identity;
            warpToScreenSpaceMatrix.m00 = HalfCameraSize.x;
            warpToScreenSpaceMatrix.m03 = HalfCameraSize.x;
            warpToScreenSpaceMatrix.m11 = HalfCameraSize.y;
            warpToScreenSpaceMatrix.m13 = HalfCameraSize.y;

            Matrix4x4 SSR_ProjectToPixelMatrix = warpToScreenSpaceMatrix * SSR_ProjectionMatrix;
            Shader.SetGlobalVector("_WorldSpaceViewDir", camForward);
            Shader.SetGlobalMatrix("_SSR_ProjectionMatrix", SSR_ProjectToPixelMatrix);
            Shader.SetGlobalVector("_SSR_ScreenSize", new Vector2(pixelWidth, pixelHeight));
        }

        public override void Execute(ScriptableRenderContext context, ref RenderingData renderingData)
        {
            ref var cameraData = ref renderingData.cameraData;
            if (cameraData.renderer.cameraColorTargetHandle == null)
                return;
            CommandBuffer cmd = CommandBufferPool.Get();
            var targetRT = DebugReflectMap ? cameraData.renderer.cameraColorTargetHandle : _customTarget;
            using (new ProfilingScope(cmd, profilingSampler))
            {
                SetupGlobalVariable(ref cameraData);
                ExecuteMainPass(cmd, null, _material, _passIndex);

                if (DebugReflectMap)
                {
                    CoreUtils.SetRenderTarget(cmd, cameraData.renderer.cameraColorTargetHandle);
                    Blitter.BlitCameraTexture(cmd, _customTarget,  cameraData.renderer.cameraColorTargetHandle);
                }
            }

            cmd.SetGlobalTexture(_ssrReflectMapID, targetRT);

            context.ExecuteCommandBuffer(cmd);
            cmd.Clear();
            CommandBufferPool.Release(cmd);
        }

#if UNITY_6000_0_OR_NEWER
        private static void ExecuteCopyColorPass(RasterCommandBuffer cmd, RTHandle sourceTexture)
        {
            Blitter.BlitTexture(cmd, sourceTexture, new Vector4(1, 1, 0, 0), 0.0f, false);
        }
        
        private static void ExecuteMainPass(RasterCommandBuffer cmd, RTHandle sourceTexture, Material material, int passIndex, MainPassData data)
        {
            s_SharedPropertyBlock.Clear();
            if (sourceTexture != null)
                s_SharedPropertyBlock.SetTexture(s_blitTexture, sourceTexture);

            // We need to set the "_BlitScaleBias" uniform for user materials with shaders relying on core Blit.hlsl to work
            s_SharedPropertyBlock.SetVector(s_blitScaleBias, new Vector4(1, 1, 0, 0));
            if (data.IsDeferred)
            { 
                SSRPassUtil.SetupGbufferTextureToMaterial(data.material, data.gBuffer);
            }

            cmd.DrawProcedural(Matrix4x4.identity, material, passIndex, MeshTopology.Triangles, 3, 1,
                s_SharedPropertyBlock);
        }
        
        public override void RecordRenderGraph(RenderGraph renderGraph, ContextContainer frameData)
        {
            UniversalResourceData resourcesData = frameData.Get<UniversalResourceData>();
            UniversalCameraData cameraData = frameData.Get<UniversalCameraData>();

            TextureHandle source, destination;

            Debug.Assert(resourcesData.cameraColor.IsValid());

            source = TextureHandle.nullHandle;
            
           if (UsePaddedTexture)
            {
                _screenWidth = cameraData.cameraTargetDescriptor.width * URPSSRSettings.GlobalResolutionScale;
                _screenHeight = cameraData.cameraTargetDescriptor.height * URPSSRSettings.GlobalResolutionScale;
                _paddedScreenWidth = Mathf.NextPowerOfTwo((int)_screenWidth);
                _paddedScreenHeight = Mathf.NextPowerOfTwo((int)_screenHeight);
            }
            else
            {
                _screenWidth = cameraData.cameraTargetDescriptor.width;
                _screenHeight = cameraData.cameraTargetDescriptor.height;
                _paddedScreenWidth = _screenWidth * URPSSRSettings.GlobalResolutionScale;
                _paddedScreenHeight = _screenHeight * URPSSRSettings.GlobalResolutionScale;
            }

            Vector2 screenResolution = new Vector2(_screenWidth, _screenHeight);
            Shader.SetGlobalVector("_ScreenResolution", screenResolution);
            if (UsePaddedTexture)
            {
                Vector2 paddedResolution = new Vector2(_paddedScreenWidth, _paddedScreenHeight);
                _paddedScale = paddedResolution / screenResolution;
                Shader.SetGlobalVector("_PaddedResolution", paddedResolution);
                Shader.SetGlobalVector("_PaddedScale", _paddedScale);

                float cx = 1.0f / (512.0f * paddedResolution.x);
                float cy = 1.0f / (512.0f * paddedResolution.y);

                Shader.SetGlobalVector("crossEpsilon", new Vector2(cx, cy));
            }
            else
            {
                _paddedScale = Vector2.one;
                Shader.SetGlobalVector("_PaddedScale", _paddedScale);
            }

            var desc = cameraData.cameraTargetDescriptor;
            var originalWidth = _screenWidth;
            var originalHeight = _screenHeight;
            desc.msaaSamples = 1;
            desc.colorFormat = RenderTextureFormat.ARGBFloat;
            desc.depthBufferBits = (int)DepthBits.None;
            desc.width = Mathf.CeilToInt(_paddedScreenWidth);
            desc.height = Mathf.CeilToInt(_paddedScreenHeight);
            destination = UniversalRenderer.CreateRenderGraphTexture(renderGraph, desc, "_SSReflectedColorMap", false);
            var colorCopyDescriptor = cameraData.cameraTargetDescriptor;
            colorCopyDescriptor.msaaSamples = 1;
            colorCopyDescriptor.depthBufferBits = (int)DepthBits.None;
            
            using (var builder =
                   renderGraph.AddRasterRenderPass<MainPassData>(passName, out var passData, profilingSampler))
            {
                passData.CameraForward = cameraData.camera.transform.forward;
                passData.PixelHeight = cameraData.camera.pixelHeight;
                passData.PixelWidth = cameraData.camera.pixelWidth;
                passData.cameraProjectionMatrix = cameraData.camera.projectionMatrix;
                passData.material = _material;
                passData.passIndex = _passIndex;
                passData.gBuffer = resourcesData.gBuffer;
                passData.inputTexture = resourcesData.activeColorTexture;
                passData.IsDeferred = _deferred;

                builder.AllowPassCulling(false);
                builder.UseTexture( passData.inputTexture, AccessFlags.Read);
                SSRPassUtil.SetupRequirementTextureToRead(builder, input, ref resourcesData);
                if (passData.IsDeferred)
                {
                    SSRPassUtil.SetupGbufferTextureToRead(builder, ref resourcesData);
                }

                builder.SetRenderAttachment(destination, 0, AccessFlags.Write);
                builder.SetRenderFunc((MainPassData data, RasterGraphContext rgContext) =>
                {
                    SetupGlobalVariable(data.cameraProjectionMatrix, data.PixelWidth, data.PixelHeight, data.CameraForward);
                    ExecuteMainPass(rgContext.cmd, data.inputTexture, data.material, data.passIndex, data);
                });
                builder.SetGlobalTextureAfterPass(destination, _ssrReflectMapID);
            }
        }

        private class CopyPassData
        {
            internal TextureHandle inputTexture;
        }

        private class MainPassData
        {
            public bool IsDeferred;
            public TextureHandle[] gBuffer;
            public Vector3 CameraForward;
            public int PixelWidth;
            public int PixelHeight;
            public Matrix4x4 cameraProjectionMatrix;
            public Material material;
            public int passIndex;
            public TextureHandle inputTexture;
        }
#endif
    }
#endif
}