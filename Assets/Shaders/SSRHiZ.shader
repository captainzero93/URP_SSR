Shader "Hidden/SSRHiZ"
{
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
            Name "SSR HiZ"

            HLSLPROGRAM
            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Core.hlsl"
            #include "Packages/com.unity.render-pipelines.core/Runtime/Utilities/Blit.hlsl"
            #pragma multi_compile_fragment _ _GBUFFER_NORMALS_OCT
            #pragma vertex Vert
            #pragma fragment frag
            float4 _HiZDepthTexture_TexelSize;
            float _SSR_HiZ_PrevDepthLevel;


            Texture2D _HiZDepthTexture;
            SamplerState sampler_HiZDepthTexture;
            

            float frag(Varyings input) : SV_Target
            {
                UNITY_SETUP_STEREO_EYE_INDEX_POST_VERTEX(input);
    float2 uv = input.texcoord.xy;

    // Sample the hierarchical depth texture at four surrounding offsets
    /*half4 minDepth = half4(
        tex2Dlod(_HiZDepthTexture, float4(uv + float2(-1.0, -1.0) * _HiZDepthTexture_TexelSize.xy, 0, _SSR_HiZ_PrevDepthLevel)).r,
        tex2Dlod(_HiZDepthTexture, float4(uv + float2(-1.0, 1.0) * _HiZDepthTexture_TexelSize.xy, 0, _SSR_HiZ_PrevDepthLevel)).r,
        tex2Dlod(_HiZDepthTexture, float4(uv + float2(1.0, -1.0) * _HiZDepthTexture_TexelSize.xy, 0, _SSR_HiZ_PrevDepthLevel)).r,
        tex2Dlod(_HiZDepthTexture, float4(uv + float2(1.0, 1.0) * _HiZDepthTexture_TexelSize.xy, 0, _SSR_HiZ_PrevDepthLevel)).r
    );*/
                
half4 minDepth = half4(
        _HiZDepthTexture.SampleLevel( sampler_HiZDepthTexture, uv, _SSR_HiZ_PrevDepthLevel, int2(-1.0,-1.0) ).r,
        _HiZDepthTexture.SampleLevel( sampler_HiZDepthTexture, uv, _SSR_HiZ_PrevDepthLevel, int2(-1.0, 1.0) ).r,
        _HiZDepthTexture.SampleLevel( sampler_HiZDepthTexture, uv, _SSR_HiZ_PrevDepthLevel, int2(1.0, -1.0) ).r,
        _HiZDepthTexture.SampleLevel( sampler_HiZDepthTexture, uv, _SSR_HiZ_PrevDepthLevel, int2(1.0, 1.0) ).r
    );
    // Return the maximum depth value from the samples
    return max(max(minDepth.r, minDepth.g), max(minDepth.b, minDepth.a));
            }
            ENDHLSL
        }
    }
}