using UnityEngine;
using UnityEngine.Rendering.Universal;
#if UNITY_6000_0_OR_NEWER
using UnityEngine.Rendering.RenderGraphModule;
#endif

namespace URPSSR
{
    public static class SSRPassUtil
    {
        public static readonly int[] s_GBufferShaderPropertyIDs = new int[]
        {
            // Contains Albedo Texture
            Shader.PropertyToID("_GBuffer0"),

            // Contains Specular Metallic Texture
            Shader.PropertyToID("_GBuffer1"),

            // Contains Normals and Smoothness, referenced as _CameraNormalsTexture in other shaders
            Shader.PropertyToID("_GBuffer2"),

            // Contains Lighting texture
            Shader.PropertyToID("_GBuffer3"),

            // Contains Depth texture, referenced as _CameraDepthTexture in other shaders (optional)
            Shader.PropertyToID("_GBuffer4"),

            // Contains Rendering Layers Texture, referenced as _CameraRenderingLayersTexture in other shaders (optional)
            Shader.PropertyToID("_GBuffer5"),

            // Contains ShadowMask texture (optional)
            Shader.PropertyToID("_GBuffer6")
        };

        public static void PrepareSSRProjectionMatrix(Camera camera)
        {
            var SSR_ProjectionMatrix =
                GL.GetGPUProjectionMatrix(camera.projectionMatrix, false);
            var HalfCameraSize = new Vector2((int)(camera.pixelWidth * 0.5f),
                (int)(camera.pixelHeight * 0.5f));

            Matrix4x4 warpToScreenSpaceMatrix = Matrix4x4.identity;
            warpToScreenSpaceMatrix.m00 = HalfCameraSize.x;
            warpToScreenSpaceMatrix.m03 = HalfCameraSize.x;
            warpToScreenSpaceMatrix.m11 = HalfCameraSize.y;
            warpToScreenSpaceMatrix.m13 = HalfCameraSize.y;

            Matrix4x4 SSR_ProjectToPixelMatrix = warpToScreenSpaceMatrix * SSR_ProjectionMatrix;
            Shader.SetGlobalMatrix("_SSR_ProjectionMatrix", SSR_ProjectToPixelMatrix);
            Shader.SetGlobalVector("_SSR_ScreenSize", new Vector2(camera.pixelWidth, camera.pixelHeight));
        }
#if UNITY_6000_0_OR_NEWER
        public static void SetupGbufferTextureToMaterial(Material material, TextureHandle[] gbufferTextures)
        {
            for (int i = 0; i < 3; i++)
            {
                material.SetTexture(s_GBufferShaderPropertyIDs[i], gbufferTextures[i]);
            }
        }
        
        //prepare gbuffer to read in record render graph 
        public static void SetupGbufferTextureToRead(IRasterRenderGraphBuilder builder, ref UniversalResourceData resourceData)
        {
            for (int i = 0; i < 3; i++)
            {
                builder.UseTexture(resourceData.gBuffer[i], AccessFlags.Read);
            }
        }

        public static void SetupRequirementTextureToRead(IRasterRenderGraphBuilder builder, ScriptableRenderPassInput input, ref UniversalResourceData resourcesData)
        {
            bool needsColor = (input & ScriptableRenderPassInput.Color) != ScriptableRenderPassInput.None;
            bool needsDepth = (input & ScriptableRenderPassInput.Depth) != ScriptableRenderPassInput.None;
            bool needsMotion = (input & ScriptableRenderPassInput.Motion) != ScriptableRenderPassInput.None;
            bool needsNormal = (input & ScriptableRenderPassInput.Normal) != ScriptableRenderPassInput.None;

            if (needsColor)
            {
                Debug.Assert(resourcesData.cameraOpaqueTexture.IsValid());
                builder.UseTexture(resourcesData.cameraOpaqueTexture);
                    
            }

            if (needsDepth)
            {
                Debug.Assert(resourcesData.cameraDepthTexture.IsValid());
                builder.UseTexture(resourcesData.cameraDepthTexture);
            }

            if (needsMotion)
            {
                Debug.Assert(resourcesData.motionVectorColor.IsValid());
                builder.UseTexture(resourcesData.motionVectorColor);
                Debug.Assert(resourcesData.motionVectorDepth.IsValid());
                builder.UseTexture(resourcesData.motionVectorDepth);
            }

            if (needsNormal)
            {
                Debug.Assert(resourcesData.cameraNormalsTexture.IsValid());
                builder.UseTexture(resourcesData.cameraNormalsTexture);
            }
        }
#endif
    }
}