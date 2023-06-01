using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class wind2 : MonoBehaviour
{
    float gravity = 9.8f;
    float frontWindForce = 5;
    float damp = 1;
    float ksStretch = 10000;
    float ksShear = 10000;
    float ksBend = 10000;
    float springIniLen = 0.1f;
    float node_mass = 1;
    float dt = 0.02f;
    const int node_row_num = 8;

    bool IsBreak = false;
    bool sum1 = false;
    bool sum2 = false;
    bool ttt = false;
    float rho = 0, rho_old=1, alpha, beta;

    Mesh mesh;
    MeshFilter meshFilter;
    private GameObject sphere;
    const int radius = 1;

    bool upLeftFixed = true;
    bool upRightFixed = true;
    bool downLeftFixed = false;
    bool downRightFixed = false;

    const int node_num = node_row_num * node_row_num;
    const int element_row_num = node_row_num - 1;
    const int element_num = element_row_num * element_row_num * 2;
    int[] element_idx = new int[element_num * 3];

    float stretchSpringLen;
    float shearSpringLen;
    float bendSpringLen;

    // 3 forces in string-nodes system
    int[,] StretchSpringDir = {{1,0},{0,1},{-1,0},{0,-1}};
    int[,] ShearSpringDir = {{2,0},{0,2},{-2,0},{0,-2}};
    int[,] BendSpringDir = {{1,1},{-1,1},{-1,-1},{1,-1}};
    int[,] NormDir = {{1,0},{0,1},{-1,0},{0,-1},};
    
    Vector3[] pos = new Vector3[node_num];
    Vector3[] pos_pre = new Vector3[node_num];
    Vector3[] pos_temp = new Vector3[node_num];
    Vector3[] pos_prepre = new Vector3[node_num];
    Vector3[] force = new Vector3[node_num];
    Vector3[] velocity = new Vector3[node_num];
    Vector3[] norm = new Vector3[node_num];
    Vector3[] velocity_pre = new Vector3[node_num];
    Vector3[] acceleration = new Vector3[node_num];
    Vector3[] acceleration_pre = new Vector3[node_num];

    // for implict Jacobian
    float[] I = {1,0,0,0,1,0,0,0,1};
    float[,,] A = new float[node_num, node_num, 9];
    float[,,] df = new float[node_num, node_num, 9];
    Vector3[] res = new Vector3[node_num];
    Vector3[] rhs = new Vector3[node_num];
    Vector3[] inertia = new Vector3[node_num];
    Vector3[] d = new Vector3[node_num];
    Vector3[] res2 = new Vector3[node_num];
    Vector3[] Ad = new Vector3[node_num];

    void Awake()
    {
        meshFilter = GetComponent<MeshFilter>();
        mesh = new Mesh();
        meshFilter.mesh = mesh;
        sphere = GameObject.Find("Sphere");

        drawTriangle();
        InitConstraint();
        StartCoroutine(StartAsync());
    }

    void drawTriangle()
    {
        float dx_space = 10.0f / (node_row_num - 1);
        for (int j = 0; j < node_row_num; j++)
        {
            for (int i = 0; i < node_row_num; i++)
            {
                int idx = j * node_row_num + i;
                pos[idx] = pos_pre[idx] = pos_temp[idx] = pos_prepre[idx] = new Vector3(i * springIniLen, 0, j * springIniLen);
            }
        }
        stretchSpringLen = Vector3.Distance(pos[1], pos[0]);
        shearSpringLen = Vector3.Distance(pos[2], pos[0]);
        bendSpringLen = Vector3.Distance(pos[1+node_row_num], pos[0]);
        // generate triangle for the mesh
        int cnt = 0;
        for (int j = 0; j < element_row_num; j++)
        {
            for (int i = 0; i < element_row_num; i++)
            {
                element_idx[cnt++] = j * node_row_num + i;
                element_idx[cnt++] = j * node_row_num + i + node_row_num;
                element_idx[cnt++] = j * node_row_num + i + 1;

                element_idx[cnt++] = j * node_row_num + i + 1;
                element_idx[cnt++] = j * node_row_num + i + node_row_num;
                element_idx[cnt++] = j * node_row_num + i + 1 + node_row_num;
            }
        }
        // add these two triangles to the mesh
        mesh.vertices = pos;
        mesh.triangles = element_idx; 
    }

    void InitConstraint()
    {
        for (int i = 0; i < node_num; i++) {
            force[i] = Vector3.zero;
            velocity[i] = Vector3.zero;
            velocity_pre[i] = Vector3.zero;
            for (int j=0; j<node_num; j++) {
                for (int n=0; n<9; n++) {
                    df[i,j,n] = 0;
                    A[i,j,n] = 0;
                }
            }
        }
        acceleration_pre = acceleration;
    }

    bool hasThisPos(int x, int y) {
        return x >= 0 && x < node_row_num && y >=0 && y < node_row_num;
    }

    void updateForce(int a, int b) {
        int index = a + b * node_row_num;
        int nodeAll = node_row_num * node_row_num;
        
        for (int type = 0; type < 3; type++) {
            for (int t = 0; t < 4; t++) {
                // stretch
                int next_idx = a + StretchSpringDir[t,0];
                int next_idy = b + StretchSpringDir[t,1];
                float ks=ksStretch;
                float len = springIniLen;
                // shear
                if (type == 1) {
                    next_idx = a + ShearSpringDir[t,0];
                    next_idy = b + ShearSpringDir[t,1];
                    ks=ksShear;
                    len = 2*springIniLen;
                }
                // bend
                else if (type == 2) {
                    next_idx = a + BendSpringDir[t,0];
                    next_idy = b + BendSpringDir[t,1];
                    ks=ksBend;
                    len = 1.414f*springIniLen;
                }

                if (hasThisPos(next_idx, next_idy)) {
                    int next_index = next_idx + next_idy * node_row_num;
                    Vector3 next_position = pos[next_index];

                    // Hooke's law
                    Vector3 pos_diff = pos[index] - next_position;
                    float distance = Vector3.Distance(pos[index],next_position);
                    Vector3 pos_dir = pos_diff/distance;
                    // calculate spring force
                    force[index] -= ks*(distance-len)*pos_dir;
                    
                    // for implict Jacobian
                    float[] mat = {pos_dir.x*pos_dir.x, pos_dir.x*pos_dir.y, pos_dir.x*pos_dir.z, 
                                    pos_dir.y*pos_dir.x, pos_dir.y*pos_dir.y, pos_dir.y*pos_dir.z, 
                                    pos_dir.z*pos_dir.x, pos_dir.z*pos_dir.y, pos_dir.z*pos_dir.z, }; 
                    for (uint n=0; n<9; n++) {
                        float hes = ks*(I[n] - (len/distance)*(I[n]-mat[n]));
                        df[index,index,n] -= hes;
                        df[index,next_index,n] += hes;
                        A[index,index,n] = I[n]-dt*dt*df[index,index,n];
                        A[index,next_index,n] = -dt*dt*df[index,next_index,n];
                    }
                } 
            }
        }
        //gravity
        force[index] += new Vector3(0,-9.8f*node_mass,0);
        //wind
        //force[index] += 2*dot(normal, float3(0,0,wind)-velocity)*normal;
        //damp
        force[index] += - damp * velocity[index];
        // pin two up corner
        if (b == 0 && (a == 0 || a == node_row_num - 1)) {
            force[index] = Vector3.zero;
        }

    }

    void updateRhs(int a, int b) {
        int index = a + b * node_row_num;
        int nodeAll = node_row_num * node_row_num;
        inertia[index] = node_mass * (2 * pos[index] - pos_pre[index]);
        res[index] = Vector3.zero;
        for (int j = 0; j < nodeAll; j++) {
            res[index].x += df[index,j, 0] * pos[j].x + df[index,j,1] * pos[j].y + df[index,j,2] * pos[j].z;
            res[index].y += df[index,j, 3] * pos[j].x + df[index,j,4] * pos[j].y + df[index,j,5] * pos[j].z;
            res[index].z += df[index,j, 6] * pos[j].x + df[index,j,7] * pos[j].y + df[index,j,8] * pos[j].z;       
        }
        rhs[index] = inertia[index] + dt*dt*force[index] - dt*dt*res[index];
    }

    void updateRes2(int a, int b) {
        int index = a + b * node_row_num;
        int nodeAll = node_row_num * node_row_num;
        d[index] = Vector3.zero;
        Ad[index] = Vector3.zero;
        res2[index] = Vector3.zero;
        Vector3 x = Vector3.zero;
        for (int j = 0; j < nodeAll; j++) {
            x.x += A[index,j,0] * pos[j].x + A[index,j,1] * pos[j].y + A[index,j,2] * pos[j].z;
            x.y += A[index,j,3] * pos[j].x + A[index,j,4] * pos[j].y + A[index,j,5] * pos[j].z;
            x.z += A[index,j,6] * pos[j].x + A[index,j,7] * pos[j].y + A[index,j,8] * pos[j].z; 
        }
        res2[index] = rhs[index] - x;
    }

    void updateD(int a, int b) {
        if (IsBreak) {return;}
        int index = a + b * node_row_num;
        int nodeAll = node_row_num * node_row_num;
        if (!sum1) {
            rho = 0;
            for (int j = 0; j < nodeAll; j++) {
                rho += res2[j].x * res2[j].x + res2[j].y * res2[j].y + res2[j].z * res2[j].z;
            }
            sum1 = true;
        }
        if (rho < 1e-8) {IsBreak=true;}
        d[index] = res2[index] + d[index]*rho/rho_old;
    }

    void updateAd(int a, int b) {
        if (IsBreak) {return;}
        int index = a + b * node_row_num;
        int nodeAll = node_row_num * node_row_num;
        Ad[index] = Vector3.zero;
        for (uint j = 0; j < nodeAll; j++) {
            Ad[index].x += A[index,j,0] * d[j].x + A[index,j,1] * d[j].y + A[index,j,2] * d[j].z;
            Ad[index].y += A[index,j,3] * d[j].x + A[index,j,4] * d[j].y + A[index,j,5] * d[j].z;
            Ad[index].z += A[index,j,6] * d[j].x + A[index,j,7] * d[j].y + A[index,j,8] * d[j].z; 
        }
    }

    void updateAlpha(int a, int b) {
        if (IsBreak) {return;}
        int nodeAll = node_row_num * node_row_num;
        if (!sum2) {
            alpha = 0;
            for (int index = 0; index < nodeAll; index++) {
                alpha += d[index].x * Ad[index].x + d[index].y * Ad[index].y + d[index].z * Ad[index].z;
            }
            sum2 = true;
        }
    }

    void updateRes2Pos(int a, int b) {
        if (IsBreak) {return;}
        int index = a + b * node_row_num;
        res2[index] = res2[index] - Ad[index]*rho/alpha;
        pos_temp[index] = pos_temp[index] + d[index]*rho/alpha;
        // pin two up corner
        if (b == 0 && (a == 0 || a == node_row_num - 1)) {
            pos_temp[index] = pos_pre[index];
        }
        rho_old = rho;
        sum1 = false;
        sum2 = false;
        print(index);
        print(res2[index].x);
    }

    void ComputeJacobianForce() {
        
        for (int i=0; i<node_row_num; i++) {
            for (int j=0; j<node_row_num; j++) {
                updateRhs(i,j);
            }
        }
        for (int i=0; i<node_row_num; i++) {
            for (int j=0; j<node_row_num; j++) {
                updateRes2(i,j);
            }
        }

        int it = 0, it_max = 20;
        
        while (it < it_max) {
            it += 1;

            for (int i=0; i<node_row_num; i++) {
                for (int j=0; j<node_row_num; j++) {
                    updateD(i,j);
                }
            }
            //print(d[5]);
            for (int i=0; i<node_row_num; i++) {
                for (int j=0; j<node_row_num; j++) {
                    updateAd(i,j);
                }
            }

            for (int i=0; i<node_row_num; i++) {
                for (int j=0; j<node_row_num; j++) {
                    updateAlpha(i,j);
                }
            }

            for (int i=0; i<node_row_num; i++) {
                for (int j=0; j<node_row_num; j++) {
                    updateRes2Pos(i,j);
                }
            }
        }
    }

    void Assemble() {
        for (int i=0; i<node_row_num; i++) {
            for (int j=0; j<node_row_num; j++) {
                updateForce(i,j);
            }
        }

        ComputeJacobianForce();

        if (IsBreak) {
            IsBreak = false;
            sum1 = false;
            sum2 = false;
            rho_old = 1;
        }
        
        // if (df[5,5,0]<-49999 && df[5,5,0]>-50001) {
        //     ttt = true;
        //     print("ttttttt");
        // }

        for (int i=0; i<node_num; i++) {
            // only for wind
            if (ttt) pos_temp[i] = pos_pre[i];
            pos[i] = pos_temp[i];
            velocity[i] = (pos[i] - pos_pre[i])/dt;

            // collision with sphere
            float s_distance = Vector3.Distance(pos[i], sphere.transform.position);
            float safe_radius = radius + 0.02f;
            if (s_distance < safe_radius) {
                Vector3 s_pos_dir = Vector3.Normalize(pos[i] - sphere.transform.position);
                pos[i] = pos[i] - (s_distance-safe_radius)*s_pos_dir;
                velocity[i] = velocity[i] - Vector3.Dot(velocity[i],s_pos_dir)*s_pos_dir;
            }

            // save the pre
            pos_prepre[i] = pos_pre[i];
            pos_pre[i] = pos[i];
            velocity_pre[i] = velocity[i];
            acceleration_pre[i] = acceleration[i];
            force[i] = Vector3.zero;
            for (int j=0; j<node_num; j++) {
                for (int n=0; n<9; n++) {
                    df[i,j,n] = 0;
                    A[i,j,n] = 0;
                }
            }
        } 
    }

    void UpdateMesh()
    {
        mesh.vertices = pos;
        mesh.RecalculateNormals();
    }

    public IEnumerator StartAsync(){
        float t = 0;
        while(true){
            t += Time.deltaTime;
            while(t > dt){
                Assemble();
                UpdateMesh();
                t -= dt;
            }
            yield return null;
        }
    }

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
    }
}
