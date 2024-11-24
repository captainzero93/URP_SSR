using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[ExecuteAlways]
public class ShaderVariableController : MonoBehaviour
{
    // Start is called before the first frame update
    void Start()
    {
        
    }

    // Update is called once per frame
    void Update()
    {
        if (Camera.main != null)
        {
            Shader.SetGlobalMatrix("_SSR_Inverse_Projection_Matrix", GL.GetGPUProjectionMatrix(Camera.main.projectionMatrix, false).inverse);
            Shader.SetGlobalMatrix("_SSR_Projection_Matrix", GL.GetGPUProjectionMatrix(Camera.main.projectionMatrix, false));
        }
    }
}
