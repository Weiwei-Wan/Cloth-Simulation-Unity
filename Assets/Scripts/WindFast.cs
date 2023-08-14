using System.Collections;
using System.Collections.Generic;
using System.IO;
using UnityEngine;

public class WindFast : MonoBehaviour
{
    StreamWriter writer;
    StreamReader reader;

    float gravity = 0.98f;
    float frontWindForce = 0.5f;
    float Damp = 1;
    float ksStretch = 1000;
    float ksShear = 1000;
    float ksBend = 1000;
    float springIniLen = 0.1f;
    float node_mass = 1f;
    float dt = 0.02f;
    const int node_row_num = 16;
    const float radius = 0.5f;

    const int node_num = node_row_num*node_row_num;
    const int element_row_num = node_row_num - 1;
    const int element_num = (node_row_num-1) * (node_row_num-1) * 2;
    int[] element_idx = new int[element_num * 3];

    Mesh mesh;
    MeshFilter meshFilter;
    private GameObject sphere;

    bool upLeftFixed = true;
    bool upRightFixed = true;
    bool downLeftFixed = false;
    bool downRightFixed = false;

    // 3 forces in string-nodes system
    int[,] StretchSpringDir = {{1,0},{0,1}};
    int[,] ShearSpringDir = {{2,0},{0,2}};
    int[,] BendSpringDir = {{1,1},{-1,1}};
    int[,] NormDir = {{1,0},{0,1},{-1,0},{0,-1},};

    float stretchSpringLen;
    float shearSpringLen;
    float bendSpringLen;

    Vector3[,] L = new Vector3[node_num, node_num];
    Vector3[,] J = new Vector3[node_num, 6*node_num];
    Vector3[,] Q = new Vector3[node_num, node_num];
    Vector3[,] ch_L = new Vector3[node_num, node_num];
    Vector3[] inertia = new Vector3[node_num];
    Vector3[] rhs = new Vector3[node_num];
    Vector3[] Jd = new Vector3[node_num];
    Vector3[] force = new Vector3[node_num];
    Vector3[] pos = new Vector3[node_num];
    Vector3[] pos_pre = new Vector3[node_num];
    Vector3[] initial_pos = new Vector3[node_num];
    Vector3[] velocity = new Vector3[node_num];
    Vector3[] d_val = new Vector3[6*node_num];

    float fast_energy = 0f;

    void WriteData(string message) {
        FileInfo file = new FileInfo(Application.dataPath + "/fast_energy1.txt");
        if (!file.Exists) {
            writer = file.CreateText();
        } else {
            writer = file.AppendText();
        }
        writer.WriteLine(message);
        writer.Flush();
        writer.Dispose();
        writer.Close();
    }
    
    void Start()
    {
        meshFilter = GetComponent<MeshFilter>();
        mesh = new Mesh();
        meshFilter.mesh = mesh;
        sphere = GameObject.Find("Sphere");
        drawTriangle();
        InitConstraint();
        CholeskyDecomp();
        StartCoroutine(StartAsync());
    }

    void drawTriangle() {
        for (int j = 0; j < node_row_num; j++)
        {
            for (int i = 0; i < node_row_num; i++)
            {
                int idx = j * node_row_num + i;
                pos[idx] = pos_pre[idx] = initial_pos[idx] = new Vector3(-i * springIniLen, 0, -j * springIniLen);
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

    bool hasThisPos(int x, int y) {
        return x >= 0 && x < node_row_num && y >=0 && y < node_row_num;
    }

    void InitConstraint() {
        // initialization
        for (int i=0; i<node_num; i++) {
            for (int j=0; j<node_num; j++) {
                L[i,j] = Vector3.zero;
            }
            for (int j=0; j<node_num*6; j++) {
                J[i,j] = Vector3.zero;
            }
        }

        for (int i=0; i<node_row_num; i++) {
            for (int j=0; j<node_row_num; j++) {
                int index = i+j*node_row_num;
                velocity[i] = Vector3.zero;
                for (int type = 0; type < 3; type++) {
                    for (int t = 0; t < 2; t++) {
                        // stretch
                        int next_x = i+StretchSpringDir[t,0];
                        int next_y = j+StretchSpringDir[t,1];
                        float ks=ksStretch;
                        // shear
                        if (type == 1) {
                            next_x = i+ShearSpringDir[t,0];
                            next_y = j+ShearSpringDir[t,1];
                            ks=ksShear;
                        }
                        // bend
                        else if (type == 2) {
                            next_x = i+BendSpringDir[t,0];
                            next_y = j+BendSpringDir[t,1];
                            ks=ksBend;
                        }
                        if (hasThisPos(next_x, next_y)) {
                            int next_index = next_x + next_y*node_row_num;
                            L[index, index] += new Vector3(ks,ks,ks);
                            L[next_index, next_index] += new Vector3(ks,ks,ks);
                            L[next_index, index] = new Vector3(-ks,-ks,-ks);
                            L[index, next_index] = new Vector3(-ks,-ks,-ks);
                            
                            J[index, index*6+(type*2+t)] += new Vector3(ks,ks,ks);
                            J[next_index, index*6+(type*2+t)] -= new Vector3(ks,ks,ks);
                        }
                    }
                }
            }
        }
    }

    void CholeskyDecomp()
    {
        for (int i=0; i<node_num; i++) {
            for (int j=0; j<node_num; j++) {
                if(i == j) {Q[i,j] = new Vector3(node_mass,node_mass,node_mass) + dt * dt * L[i,j];}
                else {Q[i,j] = dt * dt * L[i,j];}
                ch_L[i,j] = Vector3.zero;
            }
        }

        Vector3[] ch_v = new Vector3[node_num];
        for (int j=0; j<node_num; j++) {
            for (int i=j; i<node_num; i++) {
                ch_v[i] = Q[i,j];
                for (int k = 0; k < j;k++) {
                    ch_v[i].x -= ch_L[j,k].x * ch_L[i,k].x;
                    ch_v[i].y -= ch_L[j,k].y * ch_L[i,k].y;
                    ch_v[i].z -= ch_L[j,k].z * ch_L[i,k].z;
                }
                ch_L[i,j].x = ch_v[i].x / Mathf.Sqrt(ch_v[j].x);
                ch_L[i,j].y = ch_v[i].y / Mathf.Sqrt(ch_v[j].y);
                ch_L[i,j].z = ch_v[i].z / Mathf.Sqrt(ch_v[j].z);
            }
        }
        print(ch_v[50]);
    }

    Vector3[] CholeskySolve(Vector3[] rhs_)
    {
        Vector3[] resultTemp = new Vector3[node_num];
        Vector3[] result = new Vector3[node_num];
        Vector3[] rhs = new Vector3[node_num];
        for (int i=0; i<node_num; i++) {
            rhs[i] = rhs_[i];
        }

        for (int i = 0; i < node_num; i++) {
            resultTemp[i].x = rhs[i].x / ch_L[i,i].x;
            resultTemp[i].y = rhs[i].y / ch_L[i,i].y;
            resultTemp[i].z = rhs[i].z / ch_L[i,i].z;
            for (int j = 0; j < i; j++) { 
                resultTemp[i].x -= ch_L[i,j].x / ch_L[i,i].x * resultTemp[j].x;
                resultTemp[i].y -= ch_L[i,j].y / ch_L[i,i].y * resultTemp[j].y;
                resultTemp[i].z -= ch_L[i,j].z / ch_L[i,i].z * resultTemp[j].z;
            }
        }

        for (int i=node_num-1; i>=0; i--) {
            result[i].x = resultTemp[i].x / ch_L[i,i].x;
            result[i].y = resultTemp[i].y / ch_L[i,i].y;
            result[i].z = resultTemp[i].z / ch_L[i,i].z;
            for (int j = i + 1; j < node_num; j++) {
                result[i].x -= ch_L[j,i].x / ch_L[i,i].x * result[j].x;
                result[i].y -= ch_L[j,i].y / ch_L[i,i].y * result[j].y;
                result[i].z -= ch_L[j,i].z / ch_L[i,i].z * result[j].z;
            }
        }

        for (int i=0; i<node_num; i++) {
            pos[i] = result[i];
            velocity[i] = (pos[i] - pos_pre[i])/dt;
            // node_row_numlision with sphere
            float s_distance = Vector3.Distance(pos[i], sphere.transform.position);
            float safe_radius = radius + 0.05f;
            if (s_distance < safe_radius) {
                Vector3 s_pos_dir = Vector3.Normalize(pos[i] - sphere.transform.position);
                pos[i] = pos[i] - (s_distance-safe_radius)*s_pos_dir;
                //velocity[i] = velocity[i] - Vector3.Dot(velocity[i],s_pos_dir)*s_pos_dir;
            }
        }            
        return pos;
    }

    void ComputePos() {
        fast_energy = 0f;
        for (int i=0; i<node_row_num; i++) {
            for (int j=0; j<node_row_num; j++) {
                int index = i+j*node_row_num;
                //gravity
                force[index] = new Vector3(0,-gravity,0) * node_mass;
                //wind
                Vector3 norm = new Vector3(0,0,0);
                for (uint m = 0; m < 4; m++) {
                    uint n = (m + 1) % 4;
                    if (hasThisPos(i+NormDir[m,0],j+NormDir[m,1]) && hasThisPos(i+NormDir[n,0],j+NormDir[n,1])) {
                        Vector3 pos1 = pos[i+NormDir[m,0] + (j+NormDir[m,1])*node_row_num];
                        Vector3 pos2 = pos[i+NormDir[n,0] + (j+NormDir[n,1])*node_row_num];
                        norm += Vector3.Normalize(Vector3.Cross(pos1 - pos[index], pos2 - pos[index]));
                        break;
                    }
                }
                norm = Vector3.Normalize(norm);
                force[index] += 2*Vector3.Dot(norm, new Vector3(0,0,frontWindForce)-velocity[index])*norm;
                force[index] += - Damp * velocity[index];
                
                for (int type = 0; type < 3; type++) {
                    for (int t = 0; t < 2; t++) {
                        // stretch
                        int next_x = i+StretchSpringDir[t,0];
                        int next_y = j+StretchSpringDir[t,1];
                        float length = stretchSpringLen;
                        float ks=ksStretch;
                        // shear
                        if (type == 1) {
                            next_x = i+ShearSpringDir[t,0];
                            next_y = j+ShearSpringDir[t,1];
                            length = shearSpringLen;
                            ks=ksShear;
                        }
                        // bend
                        else if (type == 2) {
                            next_x = i+BendSpringDir[t,0];
                            next_y = j+BendSpringDir[t,1];
                            length = bendSpringLen;
                            ks=ksBend;
                        }
                        if (hasThisPos(next_x, next_y)) {
                            int next_index = next_x + next_y*node_row_num;
                            Vector3 diff = pos[index] - pos[next_index];
                            float _norm = Mathf.Sqrt(diff.x*diff.x + diff.y*diff.y + diff.z*diff.z);
                            d_val[index*6+type*2+t] = length * diff / _norm;
                            //energy
                            fast_energy += ks*(_norm-length)*(_norm-length)*0.5f;
                        }
                    }
                }
            }
        }
        //Debug.Log("   fast_energy:   " + fast_energy);
       // WriteData(" " + fast_energy);

        // compute Jd
        for (int i = 0; i < node_num; i++) {
            Jd[i] = Vector3.zero;
            for (int j = 0; j < 6*node_num; j++) {
                Jd[i].x += J[i,j].x * d_val[j].x;
                Jd[i].y += J[i,j].y * d_val[j].y;
                Jd[i].z += J[i,j].z * d_val[j].z;
            }
            rhs[i] = dt * dt * Jd[i] + node_mass * (2 * pos[i] - pos_pre[i]) + dt * dt * force[i];
            pos_pre[i] = pos[i];
        }
        
        pos = CholeskySolve(rhs);

        // fixed node
        if (upLeftFixed) {pos[0]=initial_pos[0];}
        if (upRightFixed) {pos[node_row_num-1]=initial_pos[node_row_num-1];}
        if (downLeftFixed) {pos[node_row_num*(node_row_num-1)]=initial_pos[node_row_num*(node_row_num-1)];}
        if (downRightFixed) {pos[node_row_num*node_row_num-1]=initial_pos[node_row_num*node_row_num-1];}
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
                ComputePos();
                UpdateMesh();
                t -= dt;
            }
            yield return null;
        }
    }

    private void Update()
    {
        if (Input.GetMouseButton(0)) {
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
