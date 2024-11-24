#pragma once
#include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Core.hlsl"
#include <HLSLSupport.cginc>
#include "Common.hlsl"

TEXTURE2D_X(_CameraDepthTexture);
UNITY_DECLARE_TEX2D(_HiZDepthTexture);

#if _SSR_DEFERRED
TEXTURE2D_X_HALF(_GBuffer0);
TEXTURE2D_X_HALF(_GBuffer1);
TEXTURE2D_X_HALF(_GBuffer2);
TEXTURE2D_X_HALF(_GBuffer3);
#else
#include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/DeclareNormalsTexture.hlsl"
TEXTURE2D_X_HALF(_SsrThinGBuffer);
#endif

float _Intensity;
float _StepSize = 4;
float _MaxDistance = 1;
float _Thickness;
float _DitherSize = 4;

#define binaryStepCount 5
#define LinearVSSteps 100
#define HiZDepthBias 0.00005
#define LinearSSSteps 100
#define LINEAR_TRACE_DEPTH_BIAS 0.05
#define LINEAR_TRACE_2D_THICKNESS 0.1
#define HizSSSteps 128


static half dither[16] = {
    0.0, 0.5, 0.125, 0.625,
    0.75, 0.25, 0.875, 0.375,
    0.187, 0.687, 0.0625, 0.562,
    0.937, 0.437, 0.812, 0.312
};

float4 TransformViewToHScreen(float3 vpos, float2 screenSize)
{
    float4 cpos = mul(UNITY_MATRIX_P, vpos);
    cpos.xy = float2(cpos.x, cpos.y * _ProjectionParams.x) * 0.5 + 0.5 * cpos.w;
    cpos.xy *= screenSize;
    return cpos;
}


inline half distanceSquared(half2 A, half2 B)
{
    A -= B;
    return dot(A, A);
}

inline half distanceSquared(half3 A, half3 B)
{
    A -= B;
    return dot(A, A);
}

bool intersectsDepthBuffer(half rayZMin, half rayZMax, half sceneZ, half layerThickness)
{
    return (rayZMax >= sceneZ - layerThickness) && (rayZMin <= sceneZ);
}

void rayIterations(Texture2D frontDepth, SamplerState DepthSampler,
                   inout half2 P,
                   inout half stepDirection, inout half end, inout int stepCount, inout int maxSteps,
                   inout bool intersecting,
                   inout half sceneZ, inout half2 dP, inout half3 Q, inout half3 dQ, inout half k, inout half dk,
                   inout half rayZMin, inout half rayZMax, inout half prevZMaxEstimate, inout bool permute,
                   inout half2 hitPixel,
                   half2 invSize, inout half layerThickness)
{
    bool stop = intersecting;
    float2 oldPos = 0;
    float performStitch = 0;
    const float stitchThredhsold = 50;
    UNITY_LOOP
    for (; (P.x * stepDirection) <= end && stepCount < maxSteps && !stop; P += dP, Q.z += dQ.z, k += dk, stepCount += 1)
    {
        rayZMin = prevZMaxEstimate;
        rayZMax = (dQ.z * 0.5 + Q.z) / (dk * 0.5 + k);
        prevZMaxEstimate = rayZMax;

        if (rayZMin > rayZMax)
        {
            Swap(rayZMin, rayZMax);
        }

        hitPixel = permute ? P.yx : P;
        sceneZ = SAMPLE_TEXTURE2D_LOD(frontDepth, DepthSampler, half4(hitPixel * invSize, 0, 0), 0).r;
        sceneZ = -LinearEyeDepth(sceneZ, _ZBufferParams);

        bool isBehind = (rayZMin <= sceneZ - LINEAR_TRACE_DEPTH_BIAS );
        float diff = -(rayZMax - sceneZ);
        if(isBehind)
        {
              if( diff < layerThickness)
              {
                  intersecting = true;
              }
        }
        else if(performStitch < 1 && diff < 10)
        {
            oldPos =  hitPixel;
        }
        stop = intersecting;
    }

    //TODO right now only trace in view space looks correct
    /*if(performStitch && !intersecting)
    {
       intersecting = true;
       hitPixel = lerp(hitPixel, oldPos, 1);
    }*/
   
    P -= dP, Q.z -= dQ.z, k -= dk;
}

bool Linear2D_Trace(Texture2D frontDepth,
                    SamplerState depthSampler,
                    half3 csOrigin,
                    half3 csDirection,
                    half4x4 projectMatrix,
                    half2 csZBufferSize,
                    half jitter,
                    int maxSteps,
                    half layerThickness,
                    half traceDistance,
                    in out half2 hitPixel,
                    int stepSize,
                    in out half3 csHitPoint,
                    in out half stepCount)
{
    half2 invSize = half2(1 / csZBufferSize.x, 1 / csZBufferSize.y);
    hitPixel = half2(-1, -1);

    half nearPlaneZ = -0.01;
    half rayLength = ((csOrigin.z + csDirection.z * traceDistance) > nearPlaneZ)
                         ? ((nearPlaneZ - csOrigin.z) / csDirection.z)
                         : traceDistance;
    half3 csEndPoint = csDirection * rayLength + csOrigin;
    half4 H0 = TransformViewToHScreen(csOrigin, csZBufferSize);
    half4 H1 = TransformViewToHScreen(csEndPoint, csZBufferSize);

    
    half k0 = 1 / H0.w;
    half k1 = 1 / H1.w;
    half2 P0 = H0.xy * k0;
    half2 P1 = H1.xy * k1;
    half3 Q0 = csOrigin * k0;
    half3 Q1 = csEndPoint * k1;
    
    half yMax = csZBufferSize.y - 0.5;
    half yMin = 0.5;
    half xMax = csZBufferSize.x - 0.5;
    half xMin = 0.5;
    half alpha = 0;

    if (P1.y > yMax || P1.y < yMin)
    {
        half yClip = (P1.y > yMax) ? yMax : yMin;
        half yAlpha = (P1.y - yClip) / (P1.y - P0.y);
        alpha = yAlpha;
    }
    if (P1.x > xMax || P1.x < xMin)
    {
        half xClip = (P1.x > xMax) ? xMax : xMin;
        half xAlpha = (P1.x - xClip) / (P1.x - P0.x);
        alpha = max(alpha, xAlpha);
    }

    P1 = lerp(P1, P0, alpha);
    k1 = lerp(k1, k0, alpha);
    Q1 = lerp(Q1, Q0, alpha);

    P1 = (distanceSquared(P0, P1) < 0.0001) ? P0 + half2(0.01, 0.01) : P1;
    half2 delta = P1 - P0;
    bool permute = false;

    if (abs(delta.x) < abs(delta.y))
    {
        permute = true;
        delta = delta.yx;
        P1 = P1.yx;
        P0 = P0.yx;
    }

    half stepDirection = sign(delta.x);
    half invdx = stepDirection / delta.x;
    half2 dP = half2(stepDirection, invdx * delta.y);
    half3 dQ = (Q1 - Q0) * invdx;
    half dk = (k1 - k0) * invdx;

    dP *= stepSize;
    dQ *= stepSize;
    dk *= stepSize;
    P0 += dP * jitter;
    Q0 += dQ * jitter;
    k0 += dk * jitter;

    half3 Q = Q0;
    half k = k0;
    half prevZMaxEstimate = csOrigin.z;
    stepCount = 0;
    half rayZMax = prevZMaxEstimate, rayZMin = prevZMaxEstimate;
    half sceneZ = 100000;
    half end = P1.x * stepDirection;
    bool intersecting = intersectsDepthBuffer(rayZMin, rayZMax, sceneZ, layerThickness);
    half2 P = P0;
    int originalStepCount = 0;

    bool traceBehind_Old = true;
    rayIterations(frontDepth, depthSampler, P, stepDirection, end, originalStepCount,
                  maxSteps,
                  intersecting, sceneZ, dP, Q, dQ, k, dk, rayZMin, rayZMax, prevZMaxEstimate, permute, hitPixel,
                  invSize, layerThickness);

    stepCount = originalStepCount;
    Q.xy += dQ.xy * stepCount;
    csHitPoint = Q * (1 / k);
    return intersecting;
}

inline float ScreenEdgeMask(float2 clipPos)
{
    float yDif = 1 - abs(clipPos.y);
    float xDif = 1 - abs(clipPos.x);
    [flatten]
    if (yDif < 0 || xDif < 0)
    {
        return 0;
    }
    float t1 = smoothstep(0, .2, yDif);
    float t2 = smoothstep(0, .1, xDif);
    return saturate(t2 * t1);
}

bool IsInfinityFar(float rawDepth)
{
    #if UNITY_REVERSED_Z
    // Case for platforms with REVERSED_Z, such as D3D.
    if (rawDepth < 0.00001)
        return true;
    #else
    // Case for platforms without REVERSED_Z, such as OpenGL.
    if(depth > 0.9999)
        return true;
    #endif
    return false;
}

bool BinarySearchVS(float3 rayStep, inout float3 samplePositionVS, inout float2 reflectUV,
                    inout float depthDelta, float oneMinusViewReflectDot)
{
    UNITY_LOOP
    for (int i = 0; i < binaryStepCount; i++)
    {
        rayStep *= 0.5f;
        [flatten]
        if (depthDelta > 0)
        {
            samplePositionVS -= rayStep;
        }
        else if (depthDelta < 0)
        {
            samplePositionVS += rayStep;
        }
        else
        {
            break;
        }

        float4 sampleUV = mul(UNITY_MATRIX_P, float4(samplePositionVS.x, samplePositionVS.y * -1,
                                                     samplePositionVS.z * -1, 1));
        sampleUV /= sampleUV.w;
        sampleUV.x *= 0.5f;
        sampleUV.y *= 0.5f;
        sampleUV.x += 0.5f;
        sampleUV.y += 0.5f;
        reflectUV = sampleUV;
        float sampledDepth = SAMPLE_TEXTURE2D_X_LOD(_CameraDepthTexture, sampler_PointClamp, sampleUV, 0).r;
        depthDelta = samplePositionVS.z - LinearEyeDepth(sampledDepth, _ZBufferParams);

        float minv = 1 / max((oneMinusViewReflectDot * float(i)), 0.001);
        if (abs(depthDelta) > minv)
        {
            return false;
        }
    }
    return true;
}

bool LinearTraceRayVS(Texture2D _DepthTexture, SamplerState sampler_DepthTexture, int NumSteps, float stepSize,
                      float2 BlueNoise, float rawDepth, inout float diff, inout float2 reflectUV,
                      inout float3 rayPos, float3 rayDir)
{
    float2 jitter = BlueNoise + 0.5;
    float StepSize = stepSize;
    //StepSize = StepSize * (jitter.x + jitter.y) + StepSize;

    float oldDepth;
    float oldDiff;
    float2 oldPos;
    reflectUV = 0;
    
    UNITY_LOOP
    for (int i = 0; i < NumSteps; i++)
    {
        rayPos += rayDir * StepSize;
        float4 sampleUV = mul(UNITY_MATRIX_P, float4(rayPos.x, rayPos.y * -1,
                                                     rayPos.z * -1, 1));
        sampleUV /= sampleUV.w;
        sampleUV.x *= 0.5f;
        sampleUV.y *= 0.5f;
        sampleUV.x += 0.5f;
        sampleUV.y += 0.5f;
        [branch]
        if (sampleUV.x >= 1 || sampleUV.x < 0 || sampleUV.y >= 1 || sampleUV.y < 0)
        {
            break;
        }

        float sampledDepth = SAMPLE_TEXTURE2D_X_LOD(_DepthTexture, sampler_DepthTexture, sampleUV, 0).r;
        UNITY_BRANCH
        diff = rayPos.z - LinearEyeDepth(sampledDepth, _ZBufferParams);
        if (diff > 0)
        {
            if( diff < StepSize )
            {
                reflectUV = sampleUV.xy;
                return true;   
            }
            if( (LinearEyeDepth(sampledDepth, _ZBufferParams) - oldDepth) >  -StepSize * 10)
            {
                float blend = (oldDiff - diff)/max(oldDiff, diff)*0.5+0.5;
                reflectUV = lerp(sampleUV, oldPos,blend);
                return true;
            }
        }
        else
        {
            oldDiff = diff;
            oldDepth =  LinearEyeDepth(sampledDepth, _ZBufferParams);
            oldPos.xy = sampleUV;
        }
    }
    return false;
}

half4 FragSSRLinearVS(Varyings input) : SV_Target
{
    UNITY_SETUP_STEREO_EYE_INDEX_POST_VERTEX(input);
    float rawDepth = SAMPLE_TEXTURE2D_X_LOD(_CameraDepthTexture, sampler_PointClamp, input.texcoord, 0).r;

    UNITY_BRANCH
    if (IsInfinityFar(rawDepth))
        return float4(input.texcoord, 0, 1);
    
    float3 normalWS = 0;
    #if _SSR_DEFERRED
        half4 gbuffer2 = SAMPLE_TEXTURE2D_X_LOD(_GBuffer2, sampler_PointClamp, input.texcoord, 0);
        normalWS = normalize(UnpackNormal(gbuffer2.xyz));
    #else
        normalWS = SampleSceneNormals(input.texcoord);
    #endif
    float3 positionWS = ComputeWorldSpacePosition(input.texcoord.xy, rawDepth, UNITY_MATRIX_I_VP);
    float3 reflectRayWS = normalize(reflect((positionWS - _WorldSpaceCameraPos), normalWS));
    float3 reflectRayVS = TransformWorldToViewDir(reflectRayWS);
    float3 rayOrigin = TransformWorldToView(positionWS);
    //handle view space in right-handed coordinates
    reflectRayVS.z *= -1;
    rayOrigin.z *= -1;

    float2 reflectUV = ComputeNormalizedDeviceCoordinates(positionWS, UNITY_MATRIX_VP);
    float3 viewDirWS = normalize(float3(positionWS.xyz) - _WorldSpaceCameraPos);

    float thickness = _StepSize * 2;
    float viewReflectDot = saturate(dot(viewDirWS, reflectRayWS));
    float oneMinusViewReflectDot = sqrt(1 - viewReflectDot);
    _StepSize /= oneMinusViewReflectDot;
    thickness /= oneMinusViewReflectDot;

    float depthDelta = 0;
    float3 samplePositionVS = rayOrigin;
    #if _SSR_DEFERRED
    if (gbuffer2.a <= 0)
       return float4(input.texcoord, 0, 1);
#endif
    float maxRayLength = LinearVSSteps * _StepSize;
    float maxDist = lerp(min(rayOrigin.z, maxRayLength), maxRayLength, viewReflectDot);
    float numSteps_f = maxDist / _StepSize;
    float totalSteps = max(numSteps_f, 0);

    //trace in view space
    bool hasTraceHit = LinearTraceRayVS(_CameraDepthTexture, sampler_PointClamp, totalSteps, _StepSize, float2(0, 0),
                                        rawDepth, depthDelta, reflectUV, samplePositionVS, reflectRayVS);
    float hasValidHit = hasTraceHit ? 1 : 0;

    //binary search refine
    float3 tempRayStep = reflectRayVS * _StepSize;
   // bool refinedHit = BinarySearchVS(tempRayStep, samplePositionVS, reflectUV, depthDelta, oneMinusViewReflectDot);
   // hasValidHit *= refinedHit ? 1 : 0;

    //check back facing
    float3 currentNormal = 0;
    #if _SSR_DEFERRED
        currentNormal = UnpackNormal(SAMPLE_TEXTURE2D_X_LOD(_GBuffer2, sampler_PointClamp, reflectUV, 0).xyz);
    #else
        currentNormal = SampleSceneNormals(reflectUV);
    #endif
 
    float backFaceDot = dot(currentNormal, reflectRayWS);
    hasValidHit *= backFaceDot > 0 ? 0 : 1;

    //fade screen edge
    float screenEdgeFade = ScreenEdgeMask((reflectUV - 0.5) * 2);
    float fadeMask = screenEdgeFade;
    //fade distance

    float3 deltaDir = rayOrigin.xyz - samplePositionVS;
    float distanceFade = dot(deltaDir, deltaDir) / (maxDist * maxDist);
    distanceFade = smoothstep(0, .5, 1 - distanceFade);
    hasValidHit *= distanceFade;
    hasValidHit *= screenEdgeFade;

    return half4(reflectUV, hasValidHit, fadeMask);
}

float4x4 _SSR_ProjectionMatrix;
float2 _SSR_ScreenSize;

half4 FragSSRLinearSS(Varyings input) : SV_Target
{
    UNITY_SETUP_STEREO_EYE_INDEX_POST_VERTEX(input);
    float rawDepth = SAMPLE_TEXTURE2D_X_LOD(_CameraDepthTexture, sampler_PointClamp, input.texcoord, 0).r;

    UNITY_BRANCH
    if (IsInfinityFar(rawDepth))
        return float4(input.texcoord, 0, 1);

    float3 normalWS = 0;
    #if _SSR_DEFERRED
        half4 gbuffer2 = SAMPLE_TEXTURE2D_X_LOD(_GBuffer2, sampler_PointClamp, input.texcoord, 0);
        normalWS = normalize(UnpackNormal(gbuffer2.xyz));
    #else
        normalWS = SampleSceneNormals(input.texcoord);
    #endif
    float3 positionWS = ComputeWorldSpacePosition(input.texcoord.xy, rawDepth, UNITY_MATRIX_I_VP);
    float3 reflectRayWS = normalize(reflect((positionWS - _WorldSpaceCameraPos), normalWS));
    float3 reflectRayVS = TransformWorldToViewDir(reflectRayWS);
    float3 rayOrigin = TransformWorldToView(positionWS);

    float2 ditherUV = input.texcoord.xy * _ScreenParams.xy;
    uint index = (uint(ditherUV.x) % 4) * 4 + uint(ditherUV.y) % 4;
    float jitter = 1 + (1 - dither[index]);
    //64 - 512
    float traceDistance = 200;
    float3 hitPointVS = rayOrigin;
    float2 reflectUV = 0; //ComputeNormalizedDeviceCoordinates(positionWS, UNITY_MATRIX_VP);
    float stepCount = 0;
    float stepSize = _StepSize * 100;
    float hasValidHit = 0;
    bool hasTraceHit = Linear2D_Trace(_CameraDepthTexture, sampler_PointClamp, rayOrigin, reflectRayVS,
                                      _SSR_ProjectionMatrix, _SSR_ScreenSize, jitter, LinearSSSteps, LINEAR_TRACE_2D_THICKNESS,
                                      traceDistance, reflectUV, stepSize, hitPointVS, stepCount);
    reflectUV /= _SSR_ScreenSize;

    hasValidHit = hasTraceHit ? 1 : 0;
    UNITY_BRANCH
    if (hasTraceHit)
    {
        float stepFade = (1 - max(2 * half(stepCount) / half(LinearSSSteps) - 1, 0));
        stepFade = stepFade * stepFade;
        hasValidHit *= stepFade;
    }
    float screenEdgeFade = ScreenEdgeMask((reflectUV - 0.5) * 2);
    float fadeMask = screenEdgeFade;
    hasValidHit *= screenEdgeFade;
    return float4(reflectUV, hasValidHit, fadeMask);
}

/////////////////////////////////////Hierarchical_Z Trace/////////////////////////////////////

float2 _ScreenResolution;
float2 _PaddedResolution;
float2 _PaddedScale;
float2 crossEpsilon;
UNITY_DECLARE_TEX2DARRAY(_DepthPyramid);
Buffer<uint2> _DepthPyramidResolutions;
TEXTURE2D_X(_SSRMainColorPaddedMap);
float3 _WorldSpaceViewDir;

///////////////////////////////////// Padded Texture /////////////////////////////////////////

//converts uv coords of padded texture to uv coords of unpadded texture
inline float2 PaddedToNativeUV(float2 uv)
{
    return uv * _PaddedScale;
}

//converts uv coords of unpadded texture to uv coords of padded texture
inline float2 NativeToPaddedUV(float2 uv)
{
    return uv / _PaddedScale;
}

inline uint2 getScreenResolution()
{
    return _PaddedResolution;
}

inline uint2 getLevelResolution(uint index)
{
    uint2 res = getScreenResolution();
    res.x = res.x >> index;
    res.y = res.y >> index;
    return res;
}

inline float2 scaledUv(float2 uv, uint index)
{
    float2 scaledScreen = getLevelResolution(index);
    float2 realScale = scaledScreen.xy / getScreenResolution();
    uv *= realScale;
    return uv;
}

#define HIZ_MAX_LEVEL 9
#define HIZ_STOP_LEVEL 0
#define HIZ_START_LEVEL 1

/////////////////////////////////////Hierarchical_Z Trace/////////////////////////////////////

inline uint NextPowerOf2(uint value)
{
    uint myNumberPowerOfTwo = 2 << firstbithigh(value - 1);
    return myNumberPowerOfTwo;
}

inline bool floatEqApprox(float a, float b)
{
    const float eps = 0.00001f;
    return abs(a - b) < eps;
}

inline float sampleDepth(float2 uv, uint index)
{
    uv = scaledUv(uv, index);
    return UNITY_SAMPLE_TEX2DARRAY(_DepthPyramid, float3(uv, index));
}

inline float2 cross_epsilon()
{
    return crossEpsilon;
}

inline float2 cell(float2 ray, float2 cell_count)
{
    float2 aligned_uv = floor(ray.xy * cell_count);
return aligned_uv;
    //return floor(ray.xy * cell_count);
}

inline float2 cell_count(float level)
{
    float2 res = getLevelResolution(level);
    return res;
}

inline bool crossed_cell_boundary(float2 cell_id_one, float2 cell_id_two)
{
    return !floatEqApprox(cell_id_one.x, cell_id_two.x) || !floatEqApprox(cell_id_one.y, cell_id_two.y);
}

inline float minimum_depth_plane(float2 ray, float level)
{
    return sampleDepth(ray, level);
}

inline float3 intersectDepthPlane(float3 o, float3 d, float t)
{
    return o + d * t;
}

inline float3 intersectCellBoundary(float3 o, float3 d, float2 cellIndex, float2 cellCount, float2 crossStep,
                                    float2 crossOffset)
{
    float2 cell_size = 1.0 / cellCount;
    float2 planes = cellIndex / cellCount + cell_size * crossStep;
    float2 solutions = (planes - o) / d.xy;
    float3 intersection_pos = o + d * min(solutions.x, solutions.y);

    intersection_pos.xy += (solutions.x < solutions.y) ? float2(crossOffset.x, 0.0) : float2(0.0, crossOffset.y);
    return intersection_pos;
}

inline float3 hiZTrace(float thickness, float3 p, float3 v, float MaxIterations, out float hit, out float iterations,
                       out bool isSky)
{
    const int rootLevel = HIZ_MAX_LEVEL;
    const int endLevel = HIZ_STOP_LEVEL;
    const int startLevel = HIZ_START_LEVEL;
    int level = HIZ_START_LEVEL;

    iterations = 0;
    isSky = false;
    hit = 0;

    [branch]
    if (v.z <= 0)
    {
        return float3(0, 0, 0);
    }

    // scale vector such that z is 1.0f (maximum depth)
    float3 d = v.xyz / v.z;
    // get the cell cross direction and a small offset to enter the next cell when doing cell crossing
    float2 crossStep = float2(d.x >= 0.0f ? 1.0f : -1.0f, d.y >= 0.0f ? 1.0f : -1.0f);
    float2 crossOffset = float2(crossStep.xy * 0.0001 );// float2(crossStep.xy * cross_epsilon() );
    crossStep.xy = saturate(crossStep.xy);

    // set current ray to original screen coordinate and depth
    float3 ray = p.xyz;
    // cross to next cell to avoid immediate self-intersection
    float2 rayCell = cell(ray.xy, cell_count(level));
    ray = intersectCellBoundary(ray, d, rayCell.xy, cell_count(level), crossStep.xy, crossOffset.xy);
    [loop]
    while (level >= endLevel
        && iterations < MaxIterations
        && ray.x >= 0 && ray.x < 1
        && ray.y >= 0 && ray.y < 1
        && ray.z > 0)
    {
        isSky = false;
        // get the cell number of the current ray
        const float2 cellCount = cell_count(level);
        const float2 oldCellIdx = cell(ray.xy, cellCount);

        // get the minimum depth plane in which the current ray resides
        float minZ = minimum_depth_plane(ray.xy, level);

        // intersect only if ray depth is below the minimum depth plane
        float3 tmpRay = ray;
        float min_minus_ray = minZ - ray.z;
        
        tmpRay = min_minus_ray > 0 ? intersectDepthPlane(tmpRay, d, min_minus_ray) : tmpRay;
        
        // get the new cell number as well
        const float2 newCellIdx = cell(tmpRay.xy, cellCount);
        // if the new cell number is different from the old cell number, a cell was crossed
        [branch]
        if (crossed_cell_boundary(oldCellIdx, newCellIdx))
        {
            // intersect the boundary of that cell instead, and go up a level for taking a larger step next iteration
            tmpRay = intersectCellBoundary(ray, d, oldCellIdx, cellCount.xy, crossStep.xy, crossOffset.xy);
            level = min(rootLevel, level + 2.0f);
        }
        else if (level == startLevel)
        {
            float minZOffset = (minZ + (1 - p.z) * thickness);
            isSky = minZ == 1;
            if(minZ >= 1)
                break;
            [flatten]
            if (abs(min_minus_ray) > HiZDepthBias)
            {
                tmpRay = intersectCellBoundary(ray, d, oldCellIdx, cellCount.xy, crossStep.xy, crossOffset.xy);
                level = HIZ_START_LEVEL + 1;
            }
        }
        // go down a level in the hi-z buffer
        --level;
        ray.xyz = tmpRay.xyz;
        ++iterations;
    }
    hit = level < endLevel ? 1 : 0;
    hit = iterations > 0 ? hit : 0;
    return ray;
}

half4 FragSSRHizSS(Varyings input) : SV_Target
{
    UNITY_SETUP_STEREO_EYE_INDEX_POST_VERTEX(input);
    float rawDepth = 1.0-sampleDepth(input.texcoord, 0);

    float2 tempUv = PaddedToNativeUV(input.texcoord);
    float2 originalUV = input.texcoord;
    input.texcoord = tempUv;

    //since we are working with padded textures, we want avoid working on padded pixels
    [branch]
    if (input.texcoord.x > 1.0f || input.texcoord.y > 1.0f)
    {
        return float4(0, 0, 0, 0);
    }

    UNITY_BRANCH
    if (IsInfinityFar(rawDepth))
        return float4(input.texcoord, 0, 0);

    float smoothness = 0;

    float3 normalWS = 0;
    #if _SSR_DEFERRED
        half4 gbuffer2 = SAMPLE_TEXTURE2D_X_LOD(_GBuffer2, sampler_PointClamp, input.texcoord, 0);
        smoothness =  gbuffer2.a;
        normalWS = normalize(UnpackNormal(gbuffer2.xyz));
    #else
        normalWS = SampleSceneNormals(input.texcoord);
    #endif
    float3 positionWS = ComputeWorldSpacePosition(input.texcoord.xy, rawDepth, UNITY_MATRIX_I_VP);
    float3 reflectRayWS = normalize(reflect((positionWS - _WorldSpaceCameraPos), normalWS));
    float3 viewDir = normalize(float3(positionWS.xyz) - _WorldSpaceCameraPos);
    float3 reflectionRay_w = reflect(viewDir, normalWS);
    float3 reflectRayVS = TransformWorldToViewDir(reflectRayWS);
    float3 rayOrigin = TransformWorldToView(positionWS);

    float4 posCS = float4(input.texcoord * 2 - 1, rawDepth, 1);
    posCS.y *= -1;
    float3 vReflectionEndPosInVS = rayOrigin + reflectRayVS * (-rayOrigin.z);
    float4 vReflectionEndPosInCS = mul(UNITY_MATRIX_P, float4(vReflectionEndPosInVS.xyz, 1));
    vReflectionEndPosInCS /= vReflectionEndPosInCS.w;
    vReflectionEndPosInCS.z = 1 - (vReflectionEndPosInCS.z);
    
    float2 ditherUV = (input.texcoord) * _ScreenParams.xy;
    uint index = (uint(ditherUV.x) % 4) * 4 + uint(ditherUV.y) % 4;
    
    float jitter = (1 + (1 - dither[index]));
    float ddd = saturate(dot(_WorldSpaceViewDir, reflectionRay_w));
    float thickness = 0.01 * (1 - ddd);
    float traceDistance = 200;
    float3 hitPointVS = rayOrigin;
    float stepCount = 0;
    float stepSize = _StepSize * 100;
    float hasValidHit = 0;

    
    UNITY_BRANCH
    if(smoothness < 0.5)
    {
        float2 realIntersectUv =  input.texcoord;
        bool hasTraceHit = Linear2D_Trace(_CameraDepthTexture, sampler_PointClamp, rayOrigin, reflectRayVS,
                                          _SSR_ProjectionMatrix, _SSR_ScreenSize, jitter, LinearSSSteps, LINEAR_TRACE_2D_THICKNESS,
                                          traceDistance, realIntersectUv, stepSize, hitPointVS, stepCount);
        realIntersectUv /= _SSR_ScreenSize;

        hasValidHit = hasTraceHit ? 1 : 0;
        UNITY_BRANCH
        if (hasTraceHit)
        {
            float stepFade = (1 - max(2 * half(stepCount) / half(LinearSSSteps) - 1, 0));
            stepFade = stepFade * stepFade;
            hasValidHit *= stepFade;
        }
        float screenEdgeFade = ScreenEdgeMask((realIntersectUv - 0.5) * 2);
        float fadeMask = screenEdgeFade;
        hasValidHit *= screenEdgeFade;
        return float4(realIntersectUv, hasValidHit, fadeMask);
    }
    else
    {
        posCS.z = 1 - (posCS.z);
        
        float3 outReflDirInTS = normalize((vReflectionEndPosInCS - posCS).xyz);
        outReflDirInTS.xy *= float2(0.5f, -0.5f);
        //convert to padded space
        outReflDirInTS.xy = NativeToPaddedUV(outReflDirInTS.xy);
        float3 outSamplePosInTS = float3(NativeToPaddedUV(input.texcoord) + outReflDirInTS.xy * 0, posCS.z + outReflDirInTS.z * -0);

        float hit = 0;
        float mask = smoothstep(0, 0.1f, ddd);
        float iterations = 0;
        bool isSky = 0;
        float3 intersectPoint = hiZTrace(thickness, outSamplePosInTS, outReflDirInTS, HizSSSteps, hit, iterations, isSky);
        float2 realIntersectUv = PaddedToNativeUV(intersectPoint.xy);
        
        if (realIntersectUv.x > 1 || realIntersectUv.x < 0 || realIntersectUv.y > 1 || realIntersectUv.y < 0)
            return 0;
        float edgeMask = 1;
        
        edgeMask = ScreenEdgeMask(realIntersectUv.xy * 2 - 1);
        // TODO use VS trace as fallback to 'stitch' the gap
        
       /* if(edgeMask >= 1 && !hit)
        {
            bool hasTraceHit = Linear2D_Trace(_CameraDepthTexture, sampler_PointClamp, rayOrigin, reflectRayVS,
                                         _SSR_ProjectionMatrix, _SSR_ScreenSize, jitter, LinearSSSteps, LINEAR_TRACE_2D_THICKNESS,
                                         traceDistance, realIntersectUv, stepSize, hitPointVS, stepCount);
            realIntersectUv /= _SSR_ScreenSize;
            hasValidHit = hasTraceHit ? 1 : 0;
            UNITY_BRANCH
            if (hasTraceHit)
            {
                float stepFade = (1 - max(2 * half(stepCount) / half(LinearSSSteps) - 1, 0));
                stepFade = stepFade * stepFade;
                hasValidHit *= stepFade;
            }
            float screenEdgeFade = ScreenEdgeMask((realIntersectUv - 0.5) * 2);
            float fadeMask = screenEdgeFade;
            hasValidHit *= screenEdgeFade;
            return float4(realIntersectUv, hasValidHit, fadeMask);
        }*/
        
        if (rawDepth == 0)
        {
            float3 normalWS = 0;
            float3 tnorm = 0;
            #if _SSR_DEFERRED
            tnorm = UnpackNormal(SAMPLE_TEXTURE2D_X_LOD(_GBuffer2, sampler_PointClamp, realIntersectUv.xy, 0));
            #else
            tnorm = SampleSceneNormals(realIntersectUv.xy);
            #endif
            float d = dot(reflectRayWS, tnorm);
            if (d > 0)
            {
                edgeMask = 0;
            }
        }
        float fadeMask = edgeMask;
        mask *= hit * edgeMask;
      
        return float4(realIntersectUv, mask, fadeMask);   
    }
}
