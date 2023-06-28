using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Rendering;

public class FastGPU : MonoBehaviour
{
    private const int THREAD_X = 8;
    private const int THREAD_Y = 8;

    float deltaTime = 0.002f;
    private int nodeNum = 16;
    private float ksStretch = 1000;
    private float ksShear = 1000;
    private float ksBend= 1000;
    private float wind = -5;
    private bool _initialized = false;

    [SerializeField]
    public GameObject sphere;
    public ComputeShader ExCompute;
    public Shader ClothShader;
    private Material material;

    private ComputeBuffer _positionBuffer;
    private ComputeBuffer _normalBuffer;
    private ComputeBuffer _velocitiesBuffer;
    private ComputeBuffer _LBuffer;
    private ComputeBuffer _JBuffer;
    private ComputeBuffer _QBuffer;
    private ComputeBuffer _ch_LBuffer;
    private ComputeBuffer _ch_vBuffer;
    private ComputeBuffer _rhsBuffer;
    private ComputeBuffer _JdBuffer;
    private ComputeBuffer _forceBuffer;
    private ComputeBuffer _pos_preBuffer;
    private ComputeBuffer _d_valBuffer;
    private ComputeBuffer _resultTempBuffer;
    private ComputeBuffer _resultBuffer;
    private ComputeBuffer _tttBuffer;

    private int _kernelInitJL;
    private int _kernelInitQ;
    private int _kernelInitCHL;
    private int _kernelUpdateDval;
    private int _kernelUpdateJd;
    private int _kernelCholeskySolve1;
    private int _kernelCholeskySolve2;
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
        _kernelInitJL = ExCompute.FindKernel("InitJL");
        _kernelInitQ = ExCompute.FindKernel("InitQ");
        _kernelInitCHL = ExCompute.FindKernel("InitCHL");
        _kernelUpdateDval = ExCompute.FindKernel("UpdateDval");
        _kernelUpdateJd = ExCompute.FindKernel("UpdateJd");
        _kernelCholeskySolve1 = ExCompute.FindKernel("CholeskySolve1");
        _kernelCholeskySolve2 = ExCompute.FindKernel("CholeskySolve2");
        _kernelStepPosition = ExCompute.FindKernel("UpdateP");
        ExCompute.SetInts("nodeNum", nodeNum, nodeNum);
        UpdateParameters();

        _positionBuffer = new ComputeBuffer(nodeNum*nodeNum,16);
        _velocitiesBuffer = new ComputeBuffer(nodeNum*nodeNum,16);
        _normalBuffer = new ComputeBuffer(nodeNum*nodeNum,16);
        _LBuffer = new ComputeBuffer(nodeNum*nodeNum*nodeNum*nodeNum,16);
        _JBuffer = new ComputeBuffer(nodeNum*nodeNum*nodeNum*nodeNum*12,16);
        _QBuffer = new ComputeBuffer(nodeNum*nodeNum*nodeNum*nodeNum,16);
        _ch_LBuffer = new ComputeBuffer(nodeNum*nodeNum*nodeNum*nodeNum,16);
        _ch_vBuffer = new ComputeBuffer(nodeNum*nodeNum,16);
        _rhsBuffer = new ComputeBuffer(nodeNum*nodeNum,16);
        _JdBuffer = new ComputeBuffer(nodeNum*nodeNum,16);
        _forceBuffer = new ComputeBuffer(nodeNum*nodeNum,16);
        _pos_preBuffer = new ComputeBuffer(nodeNum*nodeNum,16);
        _d_valBuffer = new ComputeBuffer(12*nodeNum*nodeNum,16);
        _resultTempBuffer = new ComputeBuffer(nodeNum*nodeNum,16);
        _resultBuffer = new ComputeBuffer(nodeNum*nodeNum,16);
        _tttBuffer = new ComputeBuffer(1,16);

        System.Action<int> setBufferForKernet = (k)=>{
            ExCompute.SetBuffer(k,"vel",_velocitiesBuffer);
            ExCompute.SetBuffer(k,"pos",_positionBuffer);
            ExCompute.SetBuffer(k,"norm",_normalBuffer);
            ExCompute.SetBuffer(k,"L",_LBuffer);
            ExCompute.SetBuffer(k,"J",_JBuffer);
            ExCompute.SetBuffer(k,"Q",_QBuffer);
            ExCompute.SetBuffer(k,"ch_L",_ch_LBuffer);
            ExCompute.SetBuffer(k,"ch_v",_ch_vBuffer);
            ExCompute.SetBuffer(k,"rhs",_rhsBuffer);
            ExCompute.SetBuffer(k,"Jd",_JdBuffer);
            ExCompute.SetBuffer(k,"force",_forceBuffer);
            ExCompute.SetBuffer(k,"pos_pre",_pos_preBuffer);
            ExCompute.SetBuffer(k,"d_val",_d_valBuffer);
            ExCompute.SetBuffer(k,"resultTemp",_resultTempBuffer);
            ExCompute.SetBuffer(k,"result",_resultBuffer);
            ExCompute.SetBuffer(k,"ttt",_tttBuffer);
        };

        setBufferForKernet(_kernelInitJL);
        setBufferForKernet(_kernelInitQ);
        setBufferForKernet(_kernelInitCHL);
        setBufferForKernet(_kernelUpdateDval);
        setBufferForKernet(_kernelUpdateJd);
        setBufferForKernet(_kernelCholeskySolve1);
        setBufferForKernet(_kernelCholeskySolve2); 
        setBufferForKernet(_kernelStepPosition);

        ExCompute.Dispatch(_kernelInitJL,_groupX,_groupY,1);
        ExCompute.Dispatch(_kernelInitQ,_groupX,_groupY,1);
        ExCompute.Dispatch(_kernelInitCHL,_groupX,_groupY,1);

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
                ExCompute.Dispatch(_kernelUpdateDval,_groupX,_groupY,1);
                ExCompute.Dispatch(_kernelUpdateJd,_groupX,_groupY,1);
                ExCompute.Dispatch(_kernelCholeskySolve1,_groupX,_groupY,1);
                ExCompute.Dispatch(_kernelCholeskySolve2,_groupX,_groupY,1);
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
        if (_LBuffer != null){
            _LBuffer.Release();
            _LBuffer= null;
        }
        if (_JBuffer != null){
            _JBuffer.Release();
            _JBuffer= null;
        }
        if (_QBuffer != null){
            _QBuffer.Release();
            _QBuffer= null;
        }
        if (_ch_LBuffer != null){
            _ch_LBuffer.Release();
            _ch_LBuffer= null;
        }
        if (_ch_vBuffer != null){
            _ch_vBuffer.Release();
            _ch_vBuffer= null;
        }
        if (_rhsBuffer != null){
            _rhsBuffer.Release();
            _rhsBuffer= null;
        }
        if (_JdBuffer != null){
            _JdBuffer.Release();
            _JdBuffer= null;
        }
        if (_forceBuffer != null){
            _forceBuffer.Release();
            _forceBuffer= null;
        }
        if (_pos_preBuffer != null){
            _pos_preBuffer.Release();
            _pos_preBuffer= null;
        }
        if (_d_valBuffer != null){
            _d_valBuffer.Release();
            _d_valBuffer= null;
        }
        if (_resultTempBuffer != null){
            _resultTempBuffer.Release();
            _resultTempBuffer= null;
        }
        if (_resultBuffer != null){
            _resultBuffer.Release();
            _resultBuffer= null;
        }
        if(_indexBuffer != null){
            _indexBuffer.Release();
            _indexBuffer = null;
        }
    }
}