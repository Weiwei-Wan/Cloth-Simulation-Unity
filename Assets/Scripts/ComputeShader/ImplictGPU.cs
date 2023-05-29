using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Rendering;

public class ImplictGPU : MonoBehaviour
{
    float deltaTime = 0.002f;
    private const int THREAD_X = 8;
    private const int THREAD_Y = 8;
    private int nodeNum = 16;
    private bool _initialized = false;

    [SerializeField]
    public GameObject sphere;
    public ComputeShader ImCompute;
    public Shader ClothShader;
    private Material material;

    private float ksStretch = 10000;
    private float ksShear = 10000;
    private float ksBend = 10000;
    private float wind = -5;

    private ComputeBuffer _positionBuffer;
    private ComputeBuffer _normalBuffer;
    private ComputeBuffer _velocitiesBuffer;
    private ComputeBuffer _forceBuffer;
    private ComputeBuffer _ABuffer;
    private ComputeBuffer _dfBuffer;
    private ComputeBuffer _bBuffer;
    private ComputeBuffer _velocityPreBuffer;

    private int _kernelInit;
    private int _kernelStepForce;
    private int _kernelStepJacobian;
    private int _kernelStepVelocity;
    private int _kernelStepPosition;
    private int _groupX;
    private int _groupY;
    
    public void UpdateParameters() {
        ImCompute.SetFloat("ksStretch", ksStretch);
        ImCompute.SetFloat("ksShear", ksShear);
        ImCompute.SetFloat("ksBend", ksBend);
        ImCompute.SetFloat("wind", wind);
    }
 
    public AsyncGPUReadbackRequest Initialize() {
        _kernelInit = ImCompute.FindKernel("Init");
        _kernelStepForce = ImCompute.FindKernel("UpdateF");
        _kernelStepJacobian = ImCompute.FindKernel("UpdateJ");
        _kernelStepVelocity = ImCompute.FindKernel("UpdateV");
        _kernelStepPosition = ImCompute.FindKernel("UpdateP");

        ImCompute.SetInts("nodeNum", nodeNum, nodeNum);
        UpdateParameters();

        _positionBuffer = new ComputeBuffer(nodeNum*nodeNum,16);
        _velocitiesBuffer = new ComputeBuffer(nodeNum*nodeNum,16);
        _normalBuffer = new ComputeBuffer(nodeNum*nodeNum,16);
        _forceBuffer = new ComputeBuffer(nodeNum*nodeNum,16);
        _ABuffer = new ComputeBuffer(9*nodeNum*nodeNum*nodeNum*nodeNum,16);
        _dfBuffer  = new ComputeBuffer(9*nodeNum*nodeNum*nodeNum*nodeNum,16);
        _bBuffer = new ComputeBuffer(nodeNum*nodeNum,16);
        _velocityPreBuffer = new ComputeBuffer(nodeNum*nodeNum,16);

        System.Action<int> setBufferForKernet = (k)=>{
            ImCompute.SetBuffer(k,"vel",_velocitiesBuffer);
            ImCompute.SetBuffer(k,"pos",_positionBuffer);
            ImCompute.SetBuffer(k,"norm",_normalBuffer);
            ImCompute.SetBuffer(k,"force",_forceBuffer);
            ImCompute.SetBuffer(k,"A",_ABuffer);
            ImCompute.SetBuffer(k,"df",_dfBuffer);
            ImCompute.SetBuffer(k,"b",_bBuffer);
            ImCompute.SetBuffer(k,"velocity_pre",_velocityPreBuffer);
        };

        setBufferForKernet(_kernelInit);
        setBufferForKernet(_kernelStepForce);
        setBufferForKernet(_kernelStepJacobian);
        setBufferForKernet(_kernelStepVelocity);
        setBufferForKernet(_kernelStepPosition);

        ImCompute.Dispatch(_kernelInit,_groupX,_groupY,1);

        CreateIndexBuffer();
        material.SetBuffer(ShaderIDs.pos, _positionBuffer);
        material.SetBuffer(ShaderIDs.norm,_normalBuffer);

        return AsyncGPUReadback.Request(_positionBuffer,(req)=>{
            if(req.hasError){
                Debug.LogError("Init error");
            }
            if(req.done && !req.hasError){
                _initialized = true;
            }
        });
    }

    GraphicsBuffer _indexBuffer;
	static class ShaderIDs {
		public static int pos = Shader.PropertyToID( "_positions" );
        public static int norm = Shader.PropertyToID( "_normals" );
	}

    private void CreateIndexBuffer() {
        var quadCount = (nodeNum - 1) * (nodeNum - 1);
        _indexBuffer = new GraphicsBuffer(GraphicsBuffer.Target.Index, quadCount * 6, sizeof(int));
        int[] indicies = new int[_indexBuffer.count];
        for(var x = 0; x < nodeNum - 1; x++){
            for(var y = 0; y < nodeNum - 1; y++){
                var vertexIndex = y * nodeNum + x;
                var quadIndex = y * (nodeNum - 1) + x;
                var upVertexIndex = vertexIndex + nodeNum;
                var offset = quadIndex * 6;
                indicies[offset] = vertexIndex; 
                indicies[offset + 1] = vertexIndex + 1; 
                indicies[offset + 2] = upVertexIndex; 
                indicies[offset + 3] = upVertexIndex; 
                indicies[offset + 4] = vertexIndex + 1; 
                indicies[offset + 5] = upVertexIndex + 1; 
            }
        }
        _indexBuffer.SetData(new List<int>(indicies));
    }

    public IEnumerator StartAsync(){
        yield return Initialize();
        float dt = 0;
        while (true) {
            dt += Time.deltaTime;
            while(dt > deltaTime) {
                ImCompute.SetFloat("dt", deltaTime);
                ImCompute.Dispatch(_kernelStepForce,_groupX,_groupY,1);
                for (int iter=0; iter<10; iter++) {
                    ImCompute.Dispatch(_kernelStepJacobian,_groupX,_groupY,1);
                    ImCompute.Dispatch(_kernelStepVelocity,_groupX,_groupY,1);
                }
                ImCompute.Dispatch(_kernelStepPosition,_groupX,_groupY,1);
                dt -= deltaTime;
            }
            yield return null;
            AsyncGPUReadback.WaitAllRequests();
        }
    }

    void Awake() {
        _groupX = nodeNum / THREAD_X;
        _groupY = nodeNum / THREAD_Y;
        material = new Material(ClothShader);
        UpdateParameters();
        UpdateSpherePos();
        StartCoroutine(StartAsync());
    }

    [ContextMenu("UpdateSetting")]
    private void UpdateSetting() {
        UpdateParameters();
    }

    // collision with sphere
    void UpdateSpherePos() {
        ImCompute.SetVector("sphere", sphere.transform.position);
    }

    void Update() {
        UpdateSpherePos();
    }

    void OnRenderObject() {
        if(!_initialized){
            return;
        }
        material.SetPass(0);
        Graphics.DrawProceduralNow(MeshTopology.Triangles,_indexBuffer,_indexBuffer.count,1);
    }

    void OnDestroy() {
        Debug.Log("Destroy: release buffers");
        if(_positionBuffer != null){
            _positionBuffer.Release();
            _positionBuffer = null;
        }
        if(_velocitiesBuffer != null){
            _velocitiesBuffer.Release();
            _velocitiesBuffer = null;
        }
        if(_forceBuffer != null){
            _forceBuffer.Release();
            _forceBuffer = null;
        }
        if(_ABuffer != null){
            _ABuffer.Release();
            _ABuffer = null;
        }
        if(_dfBuffer != null){
            _dfBuffer.Release();
            _dfBuffer = null;
        }
        if(_bBuffer != null){
            _bBuffer.Release();
            _bBuffer = null;
        }
        if(_velocityPreBuffer != null){
            _velocityPreBuffer.Release();
            _velocityPreBuffer = null;
        }
        if(_indexBuffer != null){
            _indexBuffer.Release();
            _indexBuffer = null;
        }
    }
}