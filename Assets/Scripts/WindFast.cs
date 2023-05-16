using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class WindFast : MonoBehaviour
{
    const int row = 16;
    const int col = 16;
    const int node_num = row*col;
    const int element_num = (row-1) * (col-1) * 2;
    int[] element_idx = new int[element_num * 3];
    const float cloth_size1 = 10.0f;
    const float cloth_size2 = 10.0f;

    Mesh mesh;
    MeshFilter meshFilter;
    private GameObject sphere;

    public float gravity;
    public float frontWindForce;
    public bool upLeftFixed;
    public bool upRightFixed;
    public bool downLeftFixed;
    public bool downRightFixed;
    public float ksStretch;
    public float ksShear;
    public float ksBend;

    Vector3[] node_pos = new Vector3[node_num];
    Vector2[] uvs = new Vector2[node_num];

    const int stretch_num = row*(col-1) + col*(row-1);
    const int shear_num = (row - 1) * (col - 1) * 2;
    const int bend_num = row*(col-2) + col*(row-2);
    const int constraint_num = stretch_num + shear_num + bend_num;

    List<int> stable_ind = new List<int>();
    float[,] L = new float[node_num*3, node_num*3];
    float[,] J = new float[node_num*3, constraint_num*3];
    float[,] Q = new float[node_num*3, node_num*3];
    float[,] ch_L = new float[node_num*3, node_num*3];
    Vector3[] inertia = new Vector3[node_num];
    Vector3[] rhs = new Vector3[node_num];
    Vector3[] Jd = new Vector3[node_num];
    Vector3[] My = new Vector3[node_num];
    Vector3[] f_ext = new Vector3[node_num];
    Vector3[] pos = new Vector3[node_num];
    Vector3[] pos_pre = new Vector3[node_num];
    Vector3[] initial_pos = new Vector3[node_num];
    Vector2[] d_springs = new Vector2[constraint_num];
    Vector3[] d_val = new Vector3[constraint_num];

    float dt = 1f;
    float node_mass = 1.0f;
    int cnt = 0;
    
    void Start()
    {
        meshFilter = GetComponent<MeshFilter>();
        mesh = new Mesh();
        meshFilter.mesh = mesh;
        sphere = GameObject.Find("Sphere");
        drawTriangle();
        InitConstraint();
        CholeskyDecomp();
    }

    void drawTriangle()
    {
        float dx_space1 = cloth_size1 / (row - 1);
        float dx_space2 = cloth_size2 / (col - 1);

        for (int j = 0; j < row; j++)
        {
            for (int i = 0; i < col; i++)
            {
                int idx = j * row + i;
                pos[idx] = pos_pre[idx] = initial_pos[idx] = node_pos[idx] = new Vector3((i-row/2)* dx_space1, (j-col/2)* dx_space2, 0);
                uvs[idx] = new Vector2((float)i / (row - 1), (float)j / (col - 1));
            }
        }

        int cnt = 0;
        for (int j = 0; j < row-1; j++)
        {
            for (int i = 0; i < col-1; i++)
            {
                element_idx[cnt++] = j * col + i;
                element_idx[cnt++] = j * col + i + col;
                element_idx[cnt++] = j * col + i + 1;

                element_idx[cnt++] = j * col + i + 1;
                element_idx[cnt++] = j * col + i + col;
                element_idx[cnt++] = j * col + i + 1 + col;
            }
        }
        //add these two triangles to the mesh
        mesh.vertices = node_pos;
        mesh.triangles = element_idx;
        mesh.uv = uvs;
    }

    void ConstructMatrix(int st, int ed, int type)
    {
        d_springs[cnt] = new Vector2(st, ed);
        float stiffness = 1.0f;
        if (type==0) {stiffness = ksStretch;}
        if (type==1) {stiffness = ksBend;}
        if (type==2) {stiffness = ksShear;}
        
        for (int i=0; i<3; i++) {
            L[st*3+i, st *3+i] += stiffness;
            L[ed*3+i, ed *3+i] += stiffness;
            L[ed*3+i, st *3+i] -= stiffness;
            L[st*3+i, ed *3+i] -= stiffness;
            J[st*3+i, cnt*3+i] += stiffness;
            J[ed*3+i, cnt*3+i] -= stiffness;
        }
        cnt += 1;
    }

    void InitConstraint()
    {
        // initialization
        for (int i=0; i<node_num*3; i++) {
            for (int j=0; j<node_num*3; j++) {
                L[i,j] = 0.0f;
            }
            for (int j=0; j<constraint_num*3; j++) {
                J[i,j] = 0.0f;
            }
        }

        for (int j = 0; j < row; j++)
        {
            for (int i = 0; i < col; i++)
            {
                int idx = j * row + i;
                // stretch
                if (i < col - 1) ConstructMatrix(idx, idx + 1, 0);
                if (j < row - 1) ConstructMatrix(idx, idx + col, 0);
                // bend
                if (i < col - 2) ConstructMatrix(idx, idx + 2, 1);
                if (j < row - 2) ConstructMatrix(idx, idx + col * 2, 1);
                // shear
                if ((i < col - 1) && (j < row - 1))
                {
                    ConstructMatrix(idx, idx + 1 + col, 2);
                    ConstructMatrix(idx + 1, idx + col, 2);
                }
            }
        }
    }

    void CholeskyDecomp()
    {
        for (int i=0; i<node_num*3; i++) {
            for (int j=0; j<node_num*3; j++) {
                if(i == j) {Q[i,j] = node_mass + dt * dt * L[i,j];}
                else {Q[i,j] = dt * dt * L[i,j];}
                ch_L[i,j] = 0;
            }
        }

        float[] ch_v = new float[node_num*3];
        for (int j=0; j<node_num*3; j++) {
            for (int i=j; i<node_num*3; i++) {
                ch_v[i] = Q[i,j];
                for(int k = 0; k < j;k++)
                {
                    ch_v[i] -= ch_L[j,k] * ch_L[i,k];
                }
                ch_L[i,j] = ch_v[i] / Mathf.Sqrt(ch_v[j]);
            }
        }
    }

    Vector3[] CholeskySolve(Vector3[] rhs_)
    {
        float[] resultTemp = new float[node_num*3];
        float[] result = new float[node_num*3];
        float[] rhs = new float[node_num*3];
        for (int i=0; i<node_num; i++) {
            rhs[3*i+0] = rhs_[i].x;
            rhs[3*i+1] = rhs_[i].y;
            rhs[3*i+2] = rhs_[i].z;
        }

        for (int i = 0; i < node_num*3; i++) {
            resultTemp[i] = rhs[i] / ch_L[i,i];
            for (int j = 0; j < i; j++) { resultTemp[i] -= ch_L[i,j] / ch_L[i,i] * resultTemp[j];}
        }

        for(int i=node_num*3-1; i>=0; i--) {
            result[i] = resultTemp[i] / ch_L[i,i];
            for (int j = i + 1; j < node_num*3; j++) {
                result[i] -= ch_L[j,i] / ch_L[i,i] * result[j];
            }
        }

        Vector3[] finalPos = new Vector3[node_num];
        for (int i=0; i<node_num; i++) {
            finalPos[i].x = result[3*i+0];
            finalPos[i].y = result[3*i+1];
            finalPos[i].z = result[3*i+2];
            // collision with sphere
            float s_distance = Vector3.Distance(finalPos[i], sphere.transform.position);
            if (s_distance < 2.0f) {
                Vector3 s_pos_dir = (finalPos[i] - sphere.transform.position)/s_distance;
                finalPos[i] = sphere.transform.position + 1.9f*s_pos_dir;
            }
        }
        return finalPos;
    }

    void GravityStep()
    {
        for (int i = 0; i < node_num; i++)
        {
            f_ext[i] = new Vector3(0, -0.001f*gravity, 0.001f*frontWindForce);
        }
    }

    void ComputePos()
    {
        for (int ic = 0; ic < constraint_num; ic++)
        {
            // Hooke's law
            int st = (int) d_springs[ic].x;
            int ed = (int) d_springs[ic].y;
            Vector3 diff = pos[st] - pos[ed];
            float norm = Mathf.Sqrt(diff.x*diff.x + diff.y*diff.y + diff.z*diff.z);
            float distance = Vector3.Distance(initial_pos[st], initial_pos[ed]);
            d_val[ic] = distance * diff / norm;
        }

        // compute Jd
        for (int i = 0; i < node_num; i++) {
            Jd[i] = Vector3.zero;
            for (int j = 0; j < constraint_num * 3; j++) {
                Jd[i].x += J[i*3+0,j] * d_val[j/3].x;
                Jd[i].y += J[i*3+1,j] * d_val[j/3].y;
                Jd[i].z += J[i*3+2,j] * d_val[j/3].z;
            }
        }
        
        for(int i = 0;i < node_num; i++)
        {
            rhs[i] = dt * dt * Jd[i] + node_mass * (2 * pos[i] - pos_pre[i]) + dt * dt * f_ext[i];
            pos_pre[i] = pos[i];
        }
        pos = CholeskySolve(rhs);

        // fixed node
        if (upLeftFixed) {pos[row*(col-1)] = initial_pos[row*(col-1)];}
        if (upRightFixed) {pos[row*col-1] = initial_pos[row*col-1];}
        if (downLeftFixed) {pos[0] = initial_pos[0];}
        if (downRightFixed) {pos[col-1] = initial_pos[col-1];}
    }

    void UpdateMesh()
    {
        for (int i = 0; i < node_num; i++)
        {
            node_pos[i] = pos[i];
        }
        mesh.vertices = node_pos;
        mesh.RecalculateNormals();
    }

    // set up the frame number each second default=60
    // void Awake()
    // {
    //     Application.targetFrameRate = 60;
    // }

    private void Update()
    {
        if (Input.GetMouseButton(0))
        {
            Vector3 screenPos = Camera.main.WorldToScreenPoint(sphere.transform.position);
            Vector3 mousePos = new Vector3(Input.mousePosition.x, Input.mousePosition.y, screenPos.z);
            sphere.transform.position = Camera.main.ScreenToWorldPoint(mousePos);
        }

        if (Input.GetKey("w")) {
            sphere.transform.Translate(0.0f, 0.0f, 0.2f);
        }
        if (Input.GetKey("s")) {
            sphere.transform.Translate(0.0f, 0.0f, -0.2f);
        }

        GravityStep();
        ComputePos();
        UpdateMesh();
    }
}
