using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Rendering;

public class FastGPU : MonoBehaviour
{
    float deltaTime = 0.002f;
    private const int THREAD_X = 8;
    private const int THREAD_Y = 8;
    private int nodeNum = 32;
    private bool _initialized = false;

    [SerializeField]
    public GameObject sphere;
    public ComputeShader ExCompute;
    public Shader ClothShader;
    private Material material;

    private float ksStretch = 10000;
    private float ksShear = 10000;
    private float ksBend= 10000;
    private float wind = -5;

    private ComputeBuffer _positionBuffer;
    private ComputeBuffer _normalBuffer;
    private ComputeBuffer _velocitiesBuffer;

    private int _kernelInit;
    private int _kernelStepVelocity;
    private int _kernelStepPosition;
    private int _groupX;
    private int _groupY;
    
    public void UpdateParameters() {
        ExCompute.SetFloat("ksStretch", ksStretch);
        ExCompute.SetFloat("ksShear", ksShear);
        ExCompute.SetFloat("ksBend", ksBend);
        ExCompute.SetFloat("wind", wind);
    }
 
    public AsyncGPUReadbackRequest Initialize() {
        _kernelInit = ExCompute.FindKernel("Init");
        _kernelStepVelocity = ExCompute.FindKernel("UpdateV");
        _kernelStepPosition = ExCompute.FindKernel("UpdateP");
        ExCompute.SetInts("nodeNum", nodeNum, nodeNum);
        UpdateParameters();

        _positionBuffer = new ComputeBuffer(nodeNum*nodeNum,16);
        _velocitiesBuffer = new ComputeBuffer(nodeNum*nodeNum,16);
        _normalBuffer = new ComputeBuffer(nodeNum*nodeNum,16);

        System.Action<int> setBufferForKernet = (k)=>{
            ExCompute.SetBuffer(k,"vel",_velocitiesBuffer);
            ExCompute.SetBuffer(k,"pos",_positionBuffer);
            ExCompute.SetBuffer(k,"norm",_normalBuffer);
        };

        setBufferForKernet(_kernelInit);
        setBufferForKernet(_kernelStepVelocity);
        setBufferForKernet(_kernelStepPosition);

        ExCompute.Dispatch(_kernelInit,_groupX,_groupY,1);

        CreateIndexBuffer();
        material.SetBuffer(ShaderIDs.pos, _positionBuffer );
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
                var vertexIndex = (y * nodeNum + x);
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
        while(true){
            dt += Time.deltaTime;
            while(dt > deltaTime){
                ExCompute.SetFloat("dt", deltaTime);
                ExCompute.Dispatch(_kernelStepVelocity,_groupX,_groupY,1);
                ExCompute.Dispatch(_kernelStepPosition,_groupX,_groupY,1);
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
        ExCompute.SetVector("sphere", sphere.transform.position);
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
        if(_indexBuffer != null){
            _indexBuffer.Release();
            _indexBuffer = null;
        }
    }
}