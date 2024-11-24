using System;
using System.Collections;
using System.Collections.Generic;
using UnityEditor;
using UnityEngine;
using UnityEngine.Rendering;

[ExecuteInEditMode]
public class CharacterReflectionProbe : MonoBehaviour
{
    private ReflectionProbe _probe;

    public RenderTexture _target;

    // Start is called before the first frame update
    private int _renderID;

    private Camera _trackCamera;
    private Vector3 _lastCameraPos;
    IEnumerator Start()
    {
        _trackCamera = Camera.main;
        _probe = GetComponent<ReflectionProbe>();
        while (true)
        {
            yield return new WaitForSeconds(0.016f);

            // render the probe over several frames and blit into TargetTexture once finished.
            _renderID = _probe.RenderProbe(_target);
        }
    }

    // Update is called once per frame
    void Update()
    {
        if (_probe == null)
            return;
        
        if (_trackCamera != null)
        {
            _lastCameraPos =  _trackCamera.transform.position;
        }
        else
        {
            _trackCamera = Camera.main;
        }

        var cameraProjectPoint = new Vector3(
            _lastCameraPos.x,
            _lastCameraPos.y * -1,
            _lastCameraPos.z
        );
        
        _probe.size = Vector3.one;
        _probe.transform.position = cameraProjectPoint;
        //  _probe.intensity = 0;
        if (_probe.IsFinishedRendering(_renderID))
        {
            Shader.SetGlobalTexture("_SSR_FallbackCube", _target);
            Shader.SetGlobalVector("_SSR_FallbackCube_Decodes", Vector4.one);
        }
    }
}