using UnityEngine;
using UnityEngine.Rendering;
using UnityEngine.Rendering.Universal;
using UnityEditor;
using System.Reflection;
using System.Linq;

public class FixUnity6SSAOResources : EditorWindow
{
    [MenuItem("Tools/Fix SSAO Build Error")]
    static void ShowWindow()
    {
        var window = GetWindow<FixUnity6SSAOResources>("Fix SSAO Build Error");
        window.minSize = new Vector2(400, 300);
        window.Show();
    }

    void OnGUI()
    {
        EditorGUILayout.Space(10);
        EditorGUILayout.LabelField("Unity 6 SSAO Build Fix", EditorStyles.boldLabel);
        EditorGUILayout.Space(5);
        
        EditorGUILayout.HelpBox("This tool helps fix the SSAO-related build errors in Unity 6 URP projects.", MessageType.Info);
        EditorGUILayout.Space(10);

        if (GUILayout.Button("Option 1: Remove SSAO from All Renderers", GUILayout.Height(30)))
        {
            RemoveSSAOFromRenderers();
        }
        
        EditorGUILayout.Space(5);
        
        if (GUILayout.Button("Option 2: Clean Null Renderer References", GUILayout.Height(30)))
        {
            CleanNullRenderers();
        }
        
        EditorGUILayout.Space(5);
        
        if (GUILayout.Button("Option 3: Disable SSAO in Graphics Settings", GUILayout.Height(30)))
        {
            DisableSSAOInGraphicsSettings();
        }
        
        EditorGUILayout.Space(10);
        EditorGUILayout.HelpBox("After using any option, try building again. If the error persists, try the next option.", MessageType.Warning);
    }

    void RemoveSSAOFromRenderers()
    {
        var urpAsset = GraphicsSettings.currentRenderPipeline as UniversalRenderPipelineAsset;
        if (urpAsset == null)
        {
            EditorUtility.DisplayDialog("Error", "No URP Asset found in Graphics Settings!", "OK");
            return;
        }

        // Get renderer list using reflection
        var fieldInfo = urpAsset.GetType().GetField("m_RendererDataList", 
            BindingFlags.NonPublic | BindingFlags.Instance);
        
        if (fieldInfo == null)
        {
            EditorUtility.DisplayDialog("Error", "Could not access renderer list!", "OK");
            return;
        }

        var rendererList = fieldInfo.GetValue(urpAsset) as ScriptableRendererData[];
        if (rendererList == null)
        {
            EditorUtility.DisplayDialog("Error", "Renderer list is null!", "OK");
            return;
        }

        int ssaoFeaturesRemoved = 0;
        
        foreach (var renderer in rendererList)
        {
            if (renderer == null) continue;
            
            var universalRenderer = renderer as UniversalRendererData;
            if (universalRenderer == null) continue;
            
            // Remove SSAO renderer features
            if (universalRenderer.rendererFeatures != null)
            {
                var featuresWithoutSSAO = universalRenderer.rendererFeatures
                    .Where(f => f != null && !f.GetType().Name.Contains("ScreenSpaceAmbientOcclusion"))
                    .ToList();
                
                if (featuresWithoutSSAO.Count < universalRenderer.rendererFeatures.Count)
                {
                    ssaoFeaturesRemoved += universalRenderer.rendererFeatures.Count - featuresWithoutSSAO.Count;
                    
                    // Use reflection to set the renderer features
                    var featuresField = universalRenderer.GetType().GetField("m_RendererFeatures", 
                        BindingFlags.NonPublic | BindingFlags.Instance);
                    if (featuresField != null)
                    {
                        featuresField.SetValue(universalRenderer, featuresWithoutSSAO);
                    }
                    
                    EditorUtility.SetDirty(universalRenderer);
                    Debug.Log($"[SSAO Fix] Removed SSAO features from {universalRenderer.name}");
                }
            }
        }
        
        if (ssaoFeaturesRemoved > 0)
        {
            AssetDatabase.SaveAssets();
            EditorUtility.DisplayDialog("Success", 
                $"Removed {ssaoFeaturesRemoved} SSAO renderer feature(s). Try building again.", "OK");
        }
        else
        {
            EditorUtility.DisplayDialog("Info", 
                "No SSAO renderer features found. Try Option 2 or 3.", "OK");
        }
    }

    void CleanNullRenderers()
    {
        var urpAsset = GraphicsSettings.currentRenderPipeline as UniversalRenderPipelineAsset;
        if (urpAsset == null)
        {
            EditorUtility.DisplayDialog("Error", "No URP Asset found in Graphics Settings!", "OK");
            return;
        }

        var fieldInfo = urpAsset.GetType().GetField("m_RendererDataList", 
            BindingFlags.NonPublic | BindingFlags.Instance);
        
        if (fieldInfo != null)
        {
            var rendererList = fieldInfo.GetValue(urpAsset) as ScriptableRendererData[];
            if (rendererList != null)
            {
                var cleanList = rendererList.Where(r => r != null).ToArray();
                int removed = rendererList.Length - cleanList.Length;
                
                if (removed > 0)
                {
                    fieldInfo.SetValue(urpAsset, cleanList);
                    EditorUtility.SetDirty(urpAsset);
                    AssetDatabase.SaveAssets();
                    
                    EditorUtility.DisplayDialog("Success", 
                        $"Removed {removed} null renderer reference(s). Try building again.", "OK");
                    Debug.Log($"[SSAO Fix] Cleaned {removed} null renderers from URP Asset");
                }
                else
                {
                    EditorUtility.DisplayDialog("Info", 
                        "No null renderers found. Try Option 1 or 3.", "OK");
                }
            }
        }
    }
    
    void DisableSSAOInGraphicsSettings()
    {
        // Try to find and disable SSAO in RenderPipelineGlobalSettings
        var globalSettings = GraphicsSettings.GetSettingsForRenderPipeline<UniversalRenderPipeline>();
        if (globalSettings != null)
        {
            EditorUtility.SetDirty(globalSettings);
            AssetDatabase.SaveAssets();
            Debug.Log("[SSAO Fix] Updated global graphics settings");
        }
        
        // Also try to clear any SSAO settings from Project Settings
        var graphicsSettingsAsset = AssetDatabase.LoadAssetAtPath<Object>("ProjectSettings/GraphicsSettings.asset");
        if (graphicsSettingsAsset != null)
        {
            EditorUtility.SetDirty(graphicsSettingsAsset);
        }
        
        EditorUtility.DisplayDialog("Info", 
            "Attempted to disable SSAO in graphics settings. Try building again.", "OK");
    }
}

// Auto-fix on build
public class SSAOBuildAutoFix : UnityEditor.Build.IPreprocessBuildWithReport
{
    public int callbackOrder => -100;
    
    public void OnPreprocessBuild(UnityEditor.Build.Reporting.BuildReport report)
    {
        Debug.Log("[SSAO AutoFix] Running pre-build SSAO check...");
        
        var urpAsset = GraphicsSettings.currentRenderPipeline as UniversalRenderPipelineAsset;
        if (urpAsset == null) return;
        
        // Auto-clean null renderers
        var fieldInfo = urpAsset.GetType().GetField("m_RendererDataList", 
            BindingFlags.NonPublic | BindingFlags.Instance);
        
        if (fieldInfo != null)
        {
            var rendererList = fieldInfo.GetValue(urpAsset) as ScriptableRendererData[];
            if (rendererList != null)
            {
                var cleanList = rendererList.Where(r => r != null).ToArray();
                if (cleanList.Length < rendererList.Length)
                {
                    fieldInfo.SetValue(urpAsset, cleanList);
                    Debug.Log($"[SSAO AutoFix] Auto-cleaned {rendererList.Length - cleanList.Length} null renderers");
                }
            }
        }
    }
}
