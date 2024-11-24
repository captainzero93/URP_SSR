using System;
using UnityEngine;
using UnityEngine.Rendering;
using UnityEngine.Rendering.Universal;
#if UNITY_6000_0_OR_NEWER
using UnityEngine.Rendering.RenderGraphModule;
#endif

namespace URPSSR
{
    public class HierarchicalZPass : ScriptableRenderPass
    {
        private const int mBUFFERSIZE = 11;
        private const int mTHREADS = 8;

        public static void ReleaseSliceBuffer(ref ComputeBuffer buffer)
        {
            if (buffer != null)
            {
                buffer.Release();
                buffer = null;
            }
        }

        public static void CreateSliceBuffer(ref ComputeBuffer buffer)
        {
            if (buffer == null)
            {
                buffer = new ComputeBuffer(mBUFFERSIZE, sizeof(int) * 2, ComputeBufferType.Default);
            }
        }

        internal struct TargetSlice
        {
            internal int slice;
            internal Vector2Int paddedResolution;
            internal Vector2Int actualResolution;
            internal Vector2 scale;

            public static implicit operator int(TargetSlice target)
            {
                return target.slice;
            }
        }

        private int finalDepthPyramidID;
        private TargetSlice[] tempSlices = new TargetSlice[mBUFFERSIZE];
        private Vector2Int[] sliceResolutions = new Vector2Int[mBUFFERSIZE];
        private ComputeBuffer depthSliceResolutions = null;
        private ComputeShader shader;
        private Vector2 screenSize;

        public HierarchicalZPass(ComputeShader shader, ComputeBuffer depthSliceBuffer)
        {
            profilingSampler = new ProfilingSampler("HierarchicalZPass");
            this.shader = shader;
            this.depthSliceResolutions = depthSliceBuffer;
        }

        public void UpdateBuffer(ComputeBuffer buffer)
        {
            this.depthSliceResolutions = buffer;
        }

        private void PrepareMembers(RenderTextureDescriptor cameraTargetDescriptor)
        {
            int width = (int)(cameraTargetDescriptor.width *
                              URPSSR.URPSSRSettings.GlobalResolutionScale);
            int height = (int)(cameraTargetDescriptor.height *
                               URPSSR.URPSSRSettings.GlobalResolutionScale);

            int paddedWidth = Mathf.NextPowerOfTwo(width);
            int paddedHeight = Mathf.NextPowerOfTwo(height);

            screenSize.x = paddedWidth;
            screenSize.y = paddedHeight;

            for (int i = 0; i < mBUFFERSIZE; i++)
            {
                tempSlices[i].paddedResolution.x = Mathf.Max(paddedWidth >> i, 1);
                tempSlices[i].paddedResolution.y = Mathf.Max(paddedHeight >> i, 1);

                tempSlices[i].actualResolution.x = Mathf.CeilToInt(width / (i + 1));
                tempSlices[i].actualResolution.y = Mathf.CeilToInt(height / (i + 1));

                tempSlices[i].slice = i;

                tempSlices[i].scale.x = tempSlices[i].paddedResolution.x / (float)paddedWidth;
                tempSlices[i].scale.y = tempSlices[i].paddedResolution.y / (float)paddedHeight;
                sliceResolutions[i] = tempSlices[i].paddedResolution;
            }

            finalDepthPyramidID = Shader.PropertyToID("_DepthPyramid");
            depthSliceResolutions.SetData(sliceResolutions);
            Shader.SetGlobalBuffer("_DepthPyramidResolutions", depthSliceResolutions);
        }

        public override void OnCameraSetup(CommandBuffer cmd, ref RenderingData renderingData)
        {
            // Calculate _DepthPyramidResolutions Start //
            //_DepthPyramidResolutions are not calculate inside compute shader, instead, calculated in CPU side
            PrepareMembers(renderingData.cameraData.cameraTargetDescriptor);
            // Calculate _DepthPyramidResolutions End //
        }

        void SetComputeShader(CommandBuffer cmd, RenderTargetIdentifier tArray, int sSlice, int dSlice, int sW, int sH,
            int dW, int dH)
        {
            cmd.SetComputeTextureParam(shader, 0, "_DepthPyramidTextureArray", tArray);
            cmd.SetComputeVectorParam(shader, "_SourceSize", new Vector2(sW, sH));
            cmd.SetComputeVectorParam(shader, "_DestinationSize", new Vector2(dW, dH));
            cmd.SetComputeIntParam(shader, "_SourceSlice", sSlice);
            cmd.SetComputeIntParam(shader, "_DestinationSlice", dSlice);
        }

        public override void Execute(ScriptableRenderContext context, ref RenderingData renderingData)
        {
            if (shader == null)
                return;

            float width = screenSize.x;
            float height = screenSize.y;

            float actualWidth = renderingData.cameraData.cameraTargetDescriptor.width *
                                URPSSRSettings.GlobalResolutionScale;
            float actualHeight = renderingData.cameraData.cameraTargetDescriptor.height *
                                 URPSSRSettings.GlobalResolutionScale;

            {
                var cmd = CommandBufferPool.Get("Init Depth Pyramid");
                cmd.GetTemporaryRTArray(finalDepthPyramidID, (int)width, (int)height, mBUFFERSIZE, 0, FilterMode.Point,
                    RenderTextureFormat.RFloat, RenderTextureReadWrite.Linear, 1, true);
                cmd.SetComputeTextureParam(shader, 1, "_DepthPyramidTextureArray", finalDepthPyramidID);
                cmd.SetComputeVectorParam(shader, "_ScreenSize", new Vector2(actualWidth, actualHeight));
                cmd.DispatchCompute(shader, 1, Mathf.CeilToInt(actualWidth / mTHREADS),
                    Mathf.CeilToInt(actualHeight / mTHREADS), 1);
                context.ExecuteCommandBuffer(cmd);
                CommandBufferPool.Release(cmd);
            }

            {
                var cmd = CommandBufferPool.Get("Calculate Depth Pyramid");
                cmd.SetExecutionFlags(CommandBufferExecutionFlags.AsyncCompute);
                for (int i = 0; i < mBUFFERSIZE - 1; i++)
                {
                    //calculate high z depth for the next scaled down buffer
                    SetComputeShader(cmd,
                        finalDepthPyramidID,
                        tempSlices[i],
                        tempSlices[i + 1],
                        tempSlices[i].paddedResolution.x,
                        tempSlices[i].paddedResolution.y,
                        tempSlices[i + 1].paddedResolution.x,
                        tempSlices[i + 1].paddedResolution.y
                    );

                    int xGroup = Mathf.CeilToInt((float)tempSlices[i + 1].actualResolution.x / mTHREADS);
                    int yGroup = Mathf.CeilToInt((float)tempSlices[i + 1].actualResolution.y / mTHREADS);
                    if (xGroup <= 0 || yGroup <= 0)
                        break;
                    cmd.DispatchCompute(shader, 0, xGroup, yGroup, 1);
                }

                context.ExecuteCommandBufferAsync(cmd, ComputeQueueType.Background);
                CommandBufferPool.Release(cmd);
            }
        }

        public override void OnCameraCleanup(CommandBuffer cmd)
        {
            base.OnCameraCleanup(cmd);
            cmd.ReleaseTemporaryRT(finalDepthPyramidID);
        }
#if UNITY_6000_0_OR_NEWER
        static void SetComputeShader(ComputeCommandBuffer cmd, ComputeShader shader, TextureHandle tArray, int sSlice, int dSlice, int sW, int sH,
            int dW, int dH)
        {
            cmd.SetComputeTextureParam(shader, 0, "_DepthPyramidTextureArray", tArray);
            cmd.SetComputeVectorParam(shader, "_SourceSSize", new Vector2(sW, sH));
            cmd.SetComputeVectorParam(shader, "_DestinationSize", new Vector2(dW, dH));
            cmd.SetComputeIntParam(shader, "_SourceSlice", sSlice);
            cmd.SetComputeIntParam(shader, "_DestinationSlice", dSlice);
        }
        
        public override void RecordRenderGraph(RenderGraph renderGraph, ContextContainer frameData)
        {
            UniversalResourceData resourcesData = frameData.Get<UniversalResourceData>();
            UniversalCameraData cameraData = frameData.Get<UniversalCameraData>();
            // Calculate _DepthPyramidResolutions Start //
            //_DepthPyramidResolutions are not calculate inside compute shader, instead, calculated in CPU side
            PrepareMembers(cameraData.cameraTargetDescriptor);
            // Calculate _DepthPyramidResolutions End //

            float width = screenSize.x;
            float height = screenSize.y;

            float actualWidth = cameraData.cameraTargetDescriptor.width *
                                URPSSRSettings.GlobalResolutionScale;
            float actualHeight = cameraData.cameraTargetDescriptor.height *
                                 URPSSRSettings.GlobalResolutionScale;

            // We need to import buffers when they are created outside of the render graph.
            //  BufferHandle inputHandle = renderGraph.ImportBuffer(inputBuffer);
            // BufferHandle outputHandle = renderGraph.ImportBuffer(outputBuffer);
            var depthDesc = new RenderTextureDescriptor();
            depthDesc.dimension = TextureDimension.Tex2DArray;
            depthDesc.width = (int)width;
            depthDesc.height = (int)height;
            depthDesc.enableRandomWrite = true;
            depthDesc.volumeDepth = mBUFFERSIZE;
            depthDesc.depthBufferBits = 0;
            depthDesc.colorFormat = RenderTextureFormat.RFloat;
            depthDesc.autoGenerateMips = false;
            depthDesc.sRGB = false;
            depthDesc.msaaSamples = 1;
            var depthPyramid = UniversalRenderer.CreateRenderGraphTexture(renderGraph, depthDesc, "_depthPyramid", true, FilterMode.Point, TextureWrapMode.Clamp);
            // Starts the recording of the render graph pass given the name of the pass
            // and outputting the data used to pass data to the execution of the render function.
            // Notice that we use "AddComputePass" when we are working with compute.
            using (var builder = renderGraph.AddComputePass("Init Depth Pyramid", out InitPassData passData))
            {
                // Set the pass data so the data can be transfered from the recording to the execution.
                passData.cs = shader;
                passData.width = width;
                passData.height = height;
                passData.actualWidth = actualWidth;
                passData.actualHeight = actualHeight;
                passData.finalDepthPyramidID = depthPyramid;
                builder.UseTexture(depthPyramid, AccessFlags.Write);
                builder.AllowPassCulling(false);
                // UseBuffer is used to setup render graph dependencies together with read and write flags.
                builder.UseTexture(resourcesData.cameraDepthTexture);
                // The execution function is also call SetRenderfunc for compute passes.
                builder.SetRenderFunc((InitPassData data, ComputeGraphContext cgContext) => ExecuteInitPass(data, cgContext.cmd));
            }


            using (var builder = renderGraph.AddComputePass("Calculate Depth Pyramid", out DownSamplePassData passData))
            {
                // Set the pass data so the data can be transfered from the recording to the execution.
                passData.cs = shader;
                passData.tempSlices = tempSlices;
                passData.finalDepthPyramidID = depthPyramid;
                builder.UseTexture(depthPyramid, AccessFlags.ReadWrite);
                builder.AllowPassCulling(false);
                //builder.EnableAsyncCompute(true);
                // UseBuffer is used to setup render graph dependencies together with read and write flags.
                builder.UseTexture(resourcesData.cameraDepthTexture);
                // The execution function is also call SetRenderfunc for compute passes.
                builder.SetRenderFunc((DownSamplePassData data, ComputeGraphContext cgContext) => ExecuteDownSamplePass(data, cgContext.cmd));
                builder.SetGlobalTextureAfterPass(depthPyramid, finalDepthPyramidID);
            }
        }

        // ExecutePass is the render function set in the render graph recordings.
        // This is good practice to avoid using variables outside of the lambda it is called from.
        // It is static to avoid using member variables which could cause unintended behaviour.
        static void ExecuteInitPass(InitPassData data, ComputeCommandBuffer cmd)
        {
            cmd.SetComputeTextureParam(data.cs, 1, "_DepthPyramidTextureArray", data.finalDepthPyramidID);
            cmd.SetComputeVectorParam(data.cs, "_ScreenSize", new Vector2(data.actualWidth, data.actualHeight));
            cmd.DispatchCompute(data.cs, 1, Mathf.CeilToInt(data.actualWidth / mTHREADS),
                Mathf.CeilToInt(data.actualHeight / mTHREADS), 1);
        }

        static void ExecuteDownSamplePass(DownSamplePassData data, ComputeCommandBuffer cmd)
        {
            for (int i = 0; i < mBUFFERSIZE - 1; i++)
            {
                //calculate high z depth for the next scaled down buffer
                SetComputeShader(cmd,data.cs,
                    data.finalDepthPyramidID,
                    data.tempSlices[i],
                    data.tempSlices[i + 1],
                    data.tempSlices[i].paddedResolution.x,
                    data.tempSlices[i].paddedResolution.y,
                    data.tempSlices[i + 1].paddedResolution.x,
                    data.tempSlices[i + 1].paddedResolution.y
                );

                int xGroup = Mathf.CeilToInt((float)data.tempSlices[i + 1].actualResolution.x / mTHREADS);
                int yGroup = Mathf.CeilToInt((float)data.tempSlices[i + 1].actualResolution.y / mTHREADS);
                if (xGroup <= 0 || yGroup <= 0)
                    break;
                cmd.DispatchCompute(data.cs, 0, xGroup, yGroup, 1);
            }
        }

        private class InitPassData
        {
            // Compute shader.
            public TextureHandle finalDepthPyramidID;
            public ComputeShader cs;
            public float width;
            public float height;
            public float actualWidth;
            public float actualHeight;
        }

        private class DownSamplePassData
        {
            // Compute shader.
            public TextureHandle finalDepthPyramidID;
            public ComputeShader cs;
            public TargetSlice[] tempSlices;
        }
#endif
    }
}