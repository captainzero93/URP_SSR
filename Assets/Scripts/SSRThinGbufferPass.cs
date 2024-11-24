using System;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Rendering;
using UnityEngine.Rendering.Universal;

#if UNITY_6000_0_OR_NEWER
using UnityEngine.Rendering.RenderGraphModule;
#endif

namespace URPSSR
{
#if UNITY_2022_1_OR_NEWER
    public class SSRThinGbufferPass : ScriptableRenderPass
    {
        private int _passIndex;

        private String _targetName;
        private RTHandle _customTarget;
        private int _ssrThinGBufferId = Shader.PropertyToID("_SsrThinGBuffer");
        private List<ShaderTagId> m_ShaderTagIdList = new List<ShaderTagId>();
        private FilteringSettings m_filteringSettings;
        private Shader _thinGbufferShader;

        public SSRThinGbufferPass(string passName)
        {
            profilingSampler = new ProfilingSampler(passName);
            m_ShaderTagIdList.Add(new ShaderTagId("SRPDefaultUnlit"));
            m_ShaderTagIdList.Add(new ShaderTagId("UniversalForward"));
            m_ShaderTagIdList.Add(new ShaderTagId("UniversalForwardOnly"));

            m_filteringSettings = new FilteringSettings(RenderQueueRange.opaque);
            if (_thinGbufferShader == null)
            {
                _thinGbufferShader = Shader.Find("SSR/ThinGBuffer");
            }
        }

        public override void OnCameraSetup(CommandBuffer cmd, ref RenderingData renderingData)
        {
        }

        public void Dispose()
        {
            _customTarget?.Release();
        }

        public override void Execute(ScriptableRenderContext context, ref RenderingData renderingData)
        {
            ref var cameraData = ref renderingData.cameraData;
            if (cameraData.renderer.cameraColorTargetHandle == null)
                return;
            CommandBuffer cmd = CommandBufferPool.Get();
            using (new ProfilingScope(cmd, profilingSampler))
            {
            }

            context.ExecuteCommandBuffer(cmd);
            cmd.Clear();
            CommandBufferPool.Release(cmd);
        }

#if UNITY_6000_0_OR_NEWER

        public override void RecordRenderGraph(RenderGraph renderGraph, ContextContainer frameData)
        {
            using (var builder = renderGraph.AddRasterRenderPass<ThinGbufferPass>(passName, out var passData, new ProfilingSampler("Thin GBuffer RenderGraph")))
            {
                // Access the relevant frame data from the Universal Render Pipeline
                UniversalRenderingData universalRenderingData = frameData.Get<UniversalRenderingData>();
                UniversalCameraData cameraData = frameData.Get<UniversalCameraData>();
                UniversalLightData lightData = frameData.Get<UniversalLightData>();

                var sortFlags = SortingCriteria.CommonOpaque;
                DrawingSettings drawSettings = RenderingUtils.CreateDrawingSettings(m_ShaderTagIdList, universalRenderingData, cameraData, lightData, sortFlags);
                drawSettings.overrideShader = _thinGbufferShader;

                var param = new RendererListParams(universalRenderingData.cullResults, drawSettings, m_filteringSettings);
                passData.RendererListHandle = renderGraph.CreateRendererList(param);

                RenderTextureDescriptor desc = new RenderTextureDescriptor(
                    cameraData.cameraTargetDescriptor.width,
                    cameraData.cameraTargetDescriptor.height);
                desc.colorFormat = RenderTextureFormat.R8;
                //use depth pre-pass 
                desc.depthBufferBits = 0;
                TextureHandle destination = UniversalRenderer.CreateRenderGraphTexture(renderGraph, desc, "_SsrThinGBuffer", true);
                passData.Destination = destination;

                builder.UseRendererList(passData.RendererListHandle);
                builder.SetRenderAttachment(passData.Destination, 0);
                builder.AllowPassCulling(false);
                builder.SetRenderFunc((ThinGbufferPass data, RasterGraphContext context) =>
                {
                    context.cmd.DrawRendererList(data.RendererListHandle);
                });

                builder.SetGlobalTextureAfterPass(passData.Destination, _ssrThinGBufferId);
            }
        }

        private class ThinGbufferPass
        {
            public RendererListHandle RendererListHandle;
            public TextureHandle Destination;
        }
#endif
    }
#endif
}