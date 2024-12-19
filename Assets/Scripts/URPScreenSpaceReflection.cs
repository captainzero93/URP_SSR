using System.Reflection;
using System.Runtime.CompilerServices;
using UnityEngine;
using UnityEngine.Rendering;
using UnityEngine.Rendering.Universal;
using URPSSR;

[assembly: InternalsVisibleTo("Unity.RenderPipelines.Universal.Editor")] 

namespace URPSSR
{
    public static class URPSSRSettings
    {
        const string GlobalInverseScaleShaderString = "_SSRGlobalInvScale";
        private static float s_globalScale = 1.0f;

        public static float GlobalResolutionScale
        {
            get { return s_globalScale; }
            set
            {
                value = Mathf.Clamp(value, 0.1f, 2.0f);
                s_globalScale = value;
                Shader.SetGlobalFloat(GlobalInverseScaleShaderString, 1.0f / s_globalScale);
            }
        }
    }
    
    public class URPScreenSpaceReflection : ScriptableRendererFeature
    {
        public enum SSRMethod
        {
            LinearVS = 0,
            LinearSS,
            HizSS,
            Composite,
        }

        [Header("SSR Setttings")] public SSRMethod SSRPass;

        public bool HizOnly;
        [Range(0.01f, 0.25f)] public float SsrStepSize = 0.1f;
        [Header("SSR Composite Settings")] public bool RenderCompositeResult = true;

        public enum SSRNativeScale
        {
            One = 1,
            Half = 2,
            Quadquarter = 4,
        }
        
        [Header("SSR Composite Settings")] 
        public SSRNativeScale SSRNativeScaleValue = SSRNativeScale.One;
        [Header("SSR Composite Settings")] [Range(0, 2f)]
        public float Intensity = 1;
        private SSRTracePass _ssrTracePass;
        private Material _ssrTraceMat;

        private SSRCompositePass _ssrCompositePass;
        private Material _ssrCompositeMat;

        private SSRThinGbufferPass _ssrThinGbufferPass;

        private ComputeBuffer _hizBuffer;
        public ComputeShader _hizShader;
        private HierarchicalZPass _hizPass;
        

        private GlobalKeyword _ssrRenderingPathKeyword;
        
        public override void Create()
        {
            _ssrRenderingPathKeyword = GlobalKeyword.Create("_SSR_DEFERRED");
            if (_ssrTracePass == null)
                _ssrTracePass = new SSRTracePass("SSR Trace");
            if (_ssrCompositePass == null)
                _ssrCompositePass = new SSRCompositePass("SSR Composite");
            if (_ssrThinGbufferPass == null)
                _ssrThinGbufferPass = new SSRThinGbufferPass("SSR Thin Gbuffer");

            if (_ssrTraceMat == null)
            {
                _ssrTraceMat = new Material(Shader.Find("PostProcess/SSR"));
            }

            if (_ssrCompositeMat == null)
            {
                _ssrCompositeMat = new Material(Shader.Find("PostProcess/SSR"));
            }

            if (_hizShader == null)
            {
                _hizShader = Instantiate(Resources.Load<ComputeShader>("HiZ_shader"));
                _hizShader.hideFlags = HideFlags.DontSave;
            }
            
            if (SSRPass == SSRMethod.HizSS)
            {
                HierarchicalZPass.ReleaseSliceBuffer(ref _hizBuffer);
                HierarchicalZPass.CreateSliceBuffer(ref _hizBuffer);

                if (_hizPass == null)
                {
                    _hizPass = new HierarchicalZPass(_hizShader, _hizBuffer);
                    _hizPass.renderPassEvent = RenderPassEvent.AfterRenderingGbuffer;
                }
                else
                {
                    _hizPass.UpdateBuffer(_hizBuffer);
                }
            }
            else
            {
                HierarchicalZPass.ReleaseSliceBuffer(ref _hizBuffer);
            }
        }

        private PropertyInfo _renderPathReflectionCached;
        
        public override void AddRenderPasses(ScriptableRenderer renderer, ref RenderingData renderingData)
        {
            var isDeferred = false;
            if (_renderPathReflectionCached == null)
            {
                var universalRendererType = renderer.GetType();
                _renderPathReflectionCached = universalRendererType.GetProperty("renderingModeRequested", BindingFlags.NonPublic | BindingFlags.Instance);   
            }
            if((RenderingMode)_renderPathReflectionCached.GetValue(renderer) == RenderingMode.Deferred)
            {
                isDeferred = true;
            }
            else
            {
                _ssrThinGbufferPass.ConfigureInput(ScriptableRenderPassInput.Normal);
                _ssrThinGbufferPass.renderPassEvent = RenderPassEvent.AfterRenderingOpaques;
                renderer.EnqueuePass(_ssrThinGbufferPass);
            }

            if(isDeferred && !Shader.IsKeywordEnabled(_ssrRenderingPathKeyword))
            {
                Shader.EnableKeyword(_ssrRenderingPathKeyword);
            }
            else if (!isDeferred && Shader.IsKeywordEnabled(_ssrRenderingPathKeyword))
            {
                Shader.DisableKeyword(_ssrRenderingPathKeyword);
            }
            
            
            if (_ssrTraceMat == null)
            {
                _ssrTraceMat = new Material(Shader.Find("PostProcess/SSR"));
            }

            if (_ssrCompositeMat == null)
            {
                _ssrCompositeMat = new Material(Shader.Find("PostProcess/SSR"));
            }
            
            URPSSRSettings.GlobalResolutionScale = 1.0f / (int)SSRNativeScaleValue;
            if (SSRPass == SSRMethod.HizSS)
            {
#if UNITY_6000_0_OR_NEWER
                _hizPass.renderPassEvent = RenderPassEvent.AfterRenderingOpaques;
#endif
                _hizPass.ConfigureInput(ScriptableRenderPassInput.Depth);
                renderer.EnqueuePass(_hizPass);
            }

            if (HizOnly)
                return;
#if UNITY_6000_0_OR_NEWER
            _ssrTracePass.renderPassEvent = (RenderPassEvent)RenderPassEvent.AfterRenderingSkybox;
#else
            _ssrTracePass.renderPassEvent = (RenderPassEvent)RenderPassEvent.AfterRenderingOpaques;
#endif
            
            _ssrTracePass.ConfigureInput(isDeferred ? ScriptableRenderPassInput.Depth : ScriptableRenderPassInput.Normal);
            _ssrTraceMat.SetFloat("_StepSize", SsrStepSize);
            _ssrTracePass.SetupMembers(_ssrTraceMat, (int)SSRPass, "_SSR_Linear", SSRPass == SSRMethod.HizSS, isDeferred);
            _ssrTracePass.DebugReflectMap = false;

            renderer.EnqueuePass(_ssrTracePass);

            if (RenderCompositeResult)
            {
                _ssrCompositePass.renderPassEvent = (RenderPassEvent)RenderPassEvent.AfterRenderingTransparents;
                _ssrCompositePass.ConfigureInput(ScriptableRenderPassInput.Depth);
                _ssrCompositePass.SetupMembers(renderingData, _ssrCompositeMat, (int)(SSRMethod.Composite), true,
                    false, SSRPass == SSRMethod.HizSS, false, Intensity, isDeferred);
                renderer.EnqueuePass(_ssrCompositePass);
            }
            
        }

        protected override void Dispose(bool disposing)
        {
            DestroyImmediate(_ssrTraceMat);
            DestroyImmediate(_ssrCompositeMat);
            _ssrTracePass.Dispose();
            _ssrCompositePass.Dispose();
            HierarchicalZPass.ReleaseSliceBuffer(ref _hizBuffer);
        }
    }
}