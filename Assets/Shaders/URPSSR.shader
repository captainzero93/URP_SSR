Shader "PostProcess/SSR"
{
    Properties
    {
        _StepSize("_StepSize", Range(0, 9)) = 0
        _MaxDistance("_MaxDistance", Range(0,50)) = 1
        _Thickness("_Thickness", Range(0,10)) = 1
        _DitherSize("_DitherSize", Range(1,20)) = 4
        _SSRIntensity("_SSRIntensity", Range(0,1)) = 0
        _FadeDistance("_FadeDistance", Range(0,10)) = 0
        _FadeDistanceFalloff("_FadeDistanceFalloff", Range(0,4)) = 0
        _EdgeDistance("_EdgeDistance", Range(0,100)) = 0
        _EdgeFalloff("EdgeFalloff", Range(0,4)) = 0
        _LodLevel ("_LodLevel", Range(0,8)) = 0
    }
    SubShader
    {
        Tags
        {
            "RenderType"="Transparent" "RenderPipeline" = "UniversalPipeline"
        }
        LOD 100
        ZWrite Off
        Cull Off
        ZTest Always
        Pass
        {
            Name "LinearVS"

            HLSLPROGRAM
            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Core.hlsl"
            #include "Packages/com.unity.render-pipelines.core/Runtime/Utilities/Blit.hlsl"
            #include "SSRAlgo.hlsl"
            #pragma multi_compile_fragment _ _GBUFFER_NORMALS_OCT
            #pragma multi_compile_fragment _ _SSR_DEFERRED
            #pragma vertex Vert
            #pragma fragment FragSSRLinearVS
            ENDHLSL
        }

        Pass
        {
            Name "LinearSS"

            HLSLPROGRAM
            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Core.hlsl"
            #include "Packages/com.unity.render-pipelines.core/Runtime/Utilities/Blit.hlsl"
            #include "SSRAlgo.hlsl"
            #pragma multi_compile_fragment _ _GBUFFER_NORMALS_OCT
            #pragma multi_compile_fragment _ _SSR_DEFERRED
            #pragma vertex Vert
            #pragma fragment FragSSRLinearSS
            ENDHLSL
        }

        Pass
        {
            Name "HizSS"

            HLSLPROGRAM
            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Core.hlsl"
            #include "Packages/com.unity.render-pipelines.core/Runtime/Utilities/Blit.hlsl"
            #include "SSRAlgo.hlsl"
            #pragma multi_compile_fragment _ _GBUFFER_NORMALS_OCT
            #pragma multi_compile_fragment _ _SSR_DEFERRED
            #pragma vertex Vert
            #pragma fragment FragSSRHizSS
            ENDHLSL
        }

        Pass
        {
            Name "Resolve"
            HLSLPROGRAM
            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Core.hlsl"
            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/BRDF.hlsl"
            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/GlobalIllumination.hlsl"
            #include "Packages/com.unity.render-pipelines.core/Runtime/Utilities/Blit.hlsl"
            #include "Common.hlsl"
            #include "SSRAlgo.hlsl"
            #pragma multi_compile_fragment _ _GBUFFER_NORMALS_OCT
            #pragma multi_compile_fragment _ _SSR_DEFERRED
            #pragma vertex Vert
            #pragma fragment frag

            float3 aces_approx(float3 v)
            {
                v *= 0.6f;
                float a = 2.51f;
                float b = 0.03f;
                float c = 2.43f;
                float d = 0.59f;
                float e = 0.14f;
                return clamp((v * (a * v + b)) / (v * (c * v + d) + e), 0.0f, 1.0f);
            }
            
            TEXTURE2D_X_HALF(_PrevSSRColor);
            TEXTURE2D_X_HALF(_SSReflectedColorMap);
            float4 _SSR_FallbackCube_Decodes;
            TEXTURECUBE(_SSR_FallbackCube);
            SAMPLER(sampler_SSR_FallbackCube);
            float _SSRIntensity;
            float _Influence;
            TEXTURE2D_X(_MotionVectorTexture);
            SAMPLER(sampler_MotionVectorTexture);

            float3 CompositeSSRColor(float2 uv, float2 reflectUV, float mask, float2 offset)
            {
                 float3 ssrColor = SAMPLE_TEXTURE2D_X_LOD(_BlitTexture, sampler_LinearClamp, reflectUV + offset, 0).rgb;
                // float3 color = SAMPLE_TEXTURE2D_X_LOD(_BlitTexture, sampler_LinearClamp, uv + offset, 0).rgb;
                 return lerp(0, ssrColor, mask * _SSRIntensity);
            }
            
            half4 frag(Varyings input) : SV_Target
            {
                UNITY_SETUP_STEREO_EYE_INDEX_POST_VERTEX(input);

                //the composite render target is not padded, thererfor when sampling a padded texture, needs to scale the uv. 
                float2 paddedTextureUV = NativeToPaddedUV(input.texcoord);
                float4 reflectBuffer = SAMPLE_TEXTURE2D_X_LOD(_SSReflectedColorMap, sampler_PointClamp, paddedTextureUV, 0);
                float2 reflectUV = (reflectBuffer.rg);
                //other textures are not padded. no need to use padded UV.
                float4 color = SAMPLE_TEXTURE2D_X_LOD(_BlitTexture, sampler_LinearClamp, input.texcoord, 0);
                float rawDepth = SAMPLE_TEXTURE2D_X_LOD(_CameraDepthTexture, sampler_PointClamp, input.texcoord, 0).r;

                UNITY_BRANCH
                if (LinearEyeDepth(rawDepth, _ZBufferParams).r / 300 >= 0.999)
                {
                    return color;
                }
                
                float3 normalWS = 0;
                #if _SSR_DEFERRED
                    half4 gbuffer2 = SAMPLE_TEXTURE2D_X_LOD(_GBuffer2, sampler_PointClamp, input.texcoord, 0);
                    normalWS = UnpackNormal(gbuffer2.xyz);
                #else
                    normalWS = SampleSceneNormals(input.texcoord);
                #endif

                //float4 gbuffer0 = SAMPLE_TEXTURE2D_X_LOD(_GBuffer0, sampler_PointClamp, input.texcoord, 0);
                float4 reflectivity = 0;
                #if _SSR_DEFERRED
                    reflectivity = SAMPLE_TEXTURE2D_X_LOD(_GBuffer1, sampler_PointClamp, input.texcoord, 0).r;
                #else
                     reflectivity = SAMPLE_TEXTURE2D_X_LOD(_SsrThinGBuffer, sampler_PointClamp, input.texcoord, 0).r;
                #endif

                float2 velocity =  SAMPLE_TEXTURE2D_X(_MotionVectorTexture, sampler_MotionVectorTexture, input.texcoord);
                float3 ssrColor =
                    SAMPLE_TEXTURE2D_X_LOD(_BlitTexture, sampler_LinearClamp, reflectUV, 0).rgb;

                float3 prevSsrColor = SAMPLE_TEXTURE2D_X_LOD(_PrevSSRColor, sampler_LinearClamp, reflectUV-velocity, 0).rgb;

                float3 minColor = 9999.0, maxColor = -9999.0;
                
                float mask = saturate(reflectBuffer.b) * reflectBuffer.a;
                            
#if UNITY_VERSION >= 60000101
  // Sample a 3x3 neighborhood to create a box in color space
                for(int x = -1; x <= 1; ++x)
                {
                    for(int y = -1; y <= 1; ++y)
                    {
                        float3 checkColor = CompositeSSRColor(input.texcoord, reflectUV, mask, float2(x,y) * _BlitTexture_TexelSize.xy);
                        minColor = min(minColor, checkColor); // Take min and max
                        maxColor = max(maxColor, checkColor);
                    }
                }
                // Clamp previous color to min/max bounding box
                prevSsrColor = clamp(prevSsrColor, minColor, 2 * maxColor);
#endif
                
                float3 ssrSpecularContribute = ssrColor.rgb;
                ssrColor = lerp(color, (ssrSpecularContribute), reflectivity);
                prevSsrColor = lerp(color, (prevSsrColor), reflectivity);
                float3 blendedColor = lerp(prevSsrColor, ssrColor, _Influence);
                float luminance = dot(blendedColor, float3(0.299, 0.587, 0.114));
                float luminanceWeight = 1.0 / (1.0 + luminance);
                blendedColor =  float4(blendedColor, 1.0) * luminanceWeight;
                float3 finalColor = lerp(color, blendedColor, mask * reflectivity);
                return float4(finalColor.rgb, 1);
            }
            ENDHLSL
        }

        Pass
        {
            ZWrite Off Cull Off
            Blend SrcAlpha OneMinusSrcAlpha
            Name "Composite"
            HLSLPROGRAM
            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Core.hlsl"
            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/BRDF.hlsl"
            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/GlobalIllumination.hlsl"
            #include "Packages/com.unity.render-pipelines.core/Runtime/Utilities/Blit.hlsl"
            #include "Common.hlsl"
            #pragma vertex Vert
            #pragma fragment frag
            float _SSRIntensity;
            TEXTURE2D_X_HALF(_ResolveSSRColor);
            
            
            half4 frag(Varyings input) : SV_Target
            {
                UNITY_SETUP_STEREO_EYE_INDEX_POST_VERTEX(input);
                float4 colorWithSSR = SAMPLE_TEXTURE2D_X_LOD(_BlitTexture, sampler_LinearClamp, input.texcoord, 0);
                return float4(colorWithSSR.rgb, _SSRIntensity);
            }
            ENDHLSL
        }

        Pass
        {
            Name "Write To History"
            HLSLPROGRAM
            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Core.hlsl"
            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/BRDF.hlsl"
            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/GlobalIllumination.hlsl"
            #include "Packages/com.unity.render-pipelines.core/Runtime/Utilities/Blit.hlsl"
            #include "Common.hlsl"
            #pragma vertex Vert
            #pragma fragment frag
            float _SSRIntensity;
            TEXTURE2D_X_HALF(_ResolveSSRColor);
            
            half4 frag(Varyings input) : SV_Target
            {
                UNITY_SETUP_STEREO_EYE_INDEX_POST_VERTEX(input);
                float3 resolveSSRColor = SAMPLE_TEXTURE2D_X_LOD(_ResolveSSRColor, sampler_LinearClamp, input.texcoord, 0).rgb * _SSRIntensity;
                float4 color = SAMPLE_TEXTURE2D_X_LOD(_BlitTexture, sampler_LinearClamp, input.texcoord, 0);
                return float4(resolveSSRColor, 1);
            }
            ENDHLSL
        }
    }
}