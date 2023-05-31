using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Rendering;

public class ImplictGPU : MonoBehaviour
{
    float deltaTime = 0.02f;
    private const int THREAD_X = 8;
    private const int THREAD_Y = 8;
    private int nodeNum = 16;
     private float ksStretch = 10000;
    private float ksShear = 10000;
    private float ksBend = 10000;
    private float wind = -5;

    private bool _initialized = false;

    [SerializeField]
    public GameObject sphere;
    public ComputeShader ImCompute;
    public Shader ClothShader;
    private Material material;

    private ComputeBuffer _positionBuffer;
    private ComputeBuffer _posPreBuffer;
    private ComputeBuffer _normalBuffer;
    private ComputeBuffer _velocitiesBuffer;
    private ComputeBuffer _forceBuffer;
    private ComputeBuffer _ABuffer;
    private ComputeBuffer _dfBuffer;
    private ComputeBuffer _resBuffer; 
    private ComputeBuffer _rhsBuffer;
    private ComputeBuffer _inertiaBuffer;
    private ComputeBuffer _dBuffer;
    private ComputeBuffer _res2Buffer;
    private ComputeBuffer _AdBuffer;

    private int _kernelInit;
    private int _kernelStepForce;
    private int _kernelStepRhs;
    private int _kernelStepRes2;
    private int _kernelStepD;
    private int _kernelStepAd;
    private int _kernelStepRes2Pos;
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
        _kernelStepPosition = ImCompute.FindKernel("UpdateP");
        _kernelStepRhs = ImCompute.FindKernel("UpdateRhs");
        _kernelStepRes2 = ImCompute.FindKernel("UpdateRes2");
        _kernelStepD = ImCompute.FindKernel("UpdateD");
        _kernelStepAd = ImCompute.FindKernel("UpdateAd");
        _kernelStepRes2Pos = ImCompute.FindKernel("UpdateRes2Pos");

        ImCompute.SetInts("nodeNum", nodeNum, nodeNum);
        UpdateParameters();

        _positionBuffer = new ComputeBuffer(nodeNum*nodeNum,16);
        _posPreBuffer = new ComputeBuffer(nodeNum*nodeNum,16);
        _velocitiesBuffer = new ComputeBuffer(nodeNum*nodeNum,16);
        _normalBuffer = new ComputeBuffer(nodeNum*nodeNum,16);
        _forceBuffer = new ComputeBuffer(nodeNum*nodeNum,16);
        _ABuffer = new ComputeBuffer(9*nodeNum*nodeNum*nodeNum*nodeNum,16);
        _dfBuffer  = new ComputeBuffer(9*nodeNum*nodeNum*nodeNum*nodeNum,16);
        _resBuffer = new ComputeBuffer(nodeNum*nodeNum,16);
        _rhsBuffer = new ComputeBuffer(nodeNum*nodeNum,16);
        _inertiaBuffer = new ComputeBuffer(nodeNum*nodeNum,16);
        _dBuffer = new ComputeBuffer(nodeNum*nodeNum,16);
        _res2Buffer = new ComputeBuffer(nodeNum*nodeNum,16);
        _AdBuffer = new ComputeBuffer(nodeNum*nodeNum,16);

        System.Action<int> setBufferForKernet = (k)=>{
            ImCompute.SetBuffer(k,"vel",_velocitiesBuffer);
            ImCompute.SetBuffer(k,"pos",_positionBuffer);
            ImCompute.SetBuffer(k,"pos_pre",_posPreBuffer);
            ImCompute.SetBuffer(k,"norm",_normalBuffer);
            ImCompute.SetBuffer(k,"force",_forceBuffer);
            ImCompute.SetBuffer(k,"A",_ABuffer);
            ImCompute.SetBuffer(k,"df",_dfBuffer);
            ImCompute.SetBuffer(k,"res",_resBuffer);
            ImCompute.SetBuffer(k,"rhs",_rhsBuffer);
            ImCompute.SetBuffer(k,"inertia",_inertiaBuffer);
            ImCompute.SetBuffer(k,"d",_dBuffer);
            ImCompute.SetBuffer(k,"res2",_res2Buffer);
            ImCompute.SetBuffer(k,"Ad",_AdBuffer);
        };

        setBufferForKernet(_kernelInit);
        setBufferForKernet(_kernelStepForce);
        setBufferForKernet(_kernelStepRhs);
        setBufferForKernet(_kernelStepRes2);
        setBufferForKernet(_kernelStepD);
        setBufferForKernet(_kernelStepAd);
        setBufferForKernet(_kernelStepRes2Pos);
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
                ImCompute.Dispatch(_kernelStepRhs,_groupX,_groupY,1);
                ImCompute.Dispatch(_kernelStepRes2,_groupX,_groupY,1);
                for (int iter=0; iter<20; iter++) {
                    ImCompute.Dispatch(_kernelStepD,_groupX,_groupY,1);
                    ImCompute.Dispatch(_kernelStepAd,_groupX,_groupY,1);
                    ImCompute.Dispatch(_kernelStepRes2Pos,_groupX,_groupY,1);
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
        if(_posPreBuffer != null){
            _posPreBuffer.Release();
            _posPreBuffer = null;
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
        if(_resBuffer != null){
            _resBuffer.Release();
            _resBuffer = null;
        }
        if(_rhsBuffer != null){
            _rhsBuffer.Release();
            _rhsBuffer = null;
        }
        if(_inertiaBuffer != null){
            _inertiaBuffer.Release();
            _inertiaBuffer = null;
        }
        if(_dBuffer != null){
            _dBuffer.Release();
            _dBuffer = null;
        }
        if(_res2Buffer != null){
            _res2Buffer.Release();
            _res2Buffer = null;
        }
        if(_AdBuffer != null){
            _AdBuffer.Release();
            _AdBuffer = null;
        }
        if(_indexBuffer != null){
            _indexBuffer.Release();
            _indexBuffer = null;
        }
    }
}