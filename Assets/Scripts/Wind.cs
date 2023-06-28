using System.Collections;
using System.Collections.Generic;
using System.IO;
using UnityEngine;

public class Wind : MonoBehaviour
{
    StreamWriter writer;
    StreamReader reader;

    float gravity = 9.8f;
    float frontWindForce = 5;
    float Damp = 1;
    float ksStretch = 10000;
    float ksShear = 10000;
    float ksBend = 10000;
    float springIniLen = 0.1f;
    float node_mass = 1;
    float dt = 0.05f;
    const int node_row_num = 16;
    const float radius = 0.5f;

    Mesh mesh;
    MeshFilter meshFilter;
    private GameObject sphere;

    public int method;
    public string[] method_options;

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
    int[,] StretchSpringDir = {{1,0},{0,1}};
    int[,] ShearSpringDir = {{2,0},{0,2}};
    int[,] BendSpringDir = {{1,1},{-1,1}};
    int[,] NormDir = {{1,0},{0,1},{-1,0},{0,-1},};
    
    Vector3[] pos = new Vector3[node_num];
    Vector3[] pos_pre = new Vector3[node_num];
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

    float ex_energy = 0f;
    float im_energy = 0f;

    void WriteData(string message) {
        FileInfo file = new FileInfo(Application.dataPath + "/im_energy05.txt");
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

    void drawTriangle() {
        for (int j = 0; j < node_row_num; j++)
        {
            for (int i = 0; i < node_row_num; i++)
            {
                int idx = j * node_row_num + i;
                pos[idx] = pos_pre[idx] = pos_prepre[idx] = new Vector3(-i * springIniLen, 0, -j * springIniLen);
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
        }
        acceleration_pre = acceleration;
    }

    bool hasThisPos(int x, int y) {
        return x >= 0 && x < node_row_num && y >=0 && y < node_row_num;
    }

    void ComputeForce() {
        ex_energy = 0f;
        for (int i=0; i<node_row_num; i++) {
            for (int j=0; j<node_row_num; j++) {
                int index = i+j*node_row_num;
                for (int type = 0; type < 3; type++) {
                    for (int t = 0; t < 2; t++) {
                        // stretch
                        int next_x = i+StretchSpringDir[t,0];
                        int next_y = j+StretchSpringDir[t,1];
                        float ks=ksStretch;
                        float length = stretchSpringLen;
                        // shear
                        if (type == 1) {
                            next_x = i+ShearSpringDir[t,0];
                            next_y = j+ShearSpringDir[t,1];
                            ks=ksShear;
                            length = shearSpringLen;
                        }
                        // bend
                        else if (type == 2) {
                            next_x = i+BendSpringDir[t,0];
                            next_y = j+BendSpringDir[t,1];
                            ks=ksBend;
                            length = bendSpringLen;
                        }
                        if (hasThisPos(next_x, next_y)) {
                            int next_index = next_x + next_y*node_row_num;
                            // Hooke's law
                            Vector3 pos_diff = pos[index]-pos[next_index];
                            float distance = Vector3.Distance(pos[next_index],pos[index]);
                            Vector3 pos_dir = pos_diff/distance;
                            // calculate spring force
                            force[index] -= ks*(distance-length)*pos_dir;
                            force[next_index] += ks*(distance-length)*pos_dir;
                            //energy
                            ex_energy += ks*(distance-length)*(distance-length)*0.5f;
                        }
                    }
                }
               
                //gravity
                force[index] += new Vector3(0,-gravity,0) * node_mass;
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
                //damp
                force[index] += - Damp * velocity[index];
            }
        }
        Debug.Log("   ex_energy:   " + ex_energy);
        WriteData(" " + ex_energy);
    }

    void ComputeJacobianForce()
    {
        // initialization
        for (int i=0; i<node_num; i++) {
            for (int j=0; j<node_num; j++) {
                for (int t=0; t<9; t++) {
                    A[i,j,t] = 0.0f;
                    df[i,j,t] = 0.0f;
                }
            }
        }
        
        im_energy = 0f;
        // internal spring force
        for (int i=0; i<node_row_num; i++) {
            for (int j=0; j<node_row_num; j++) {
                int index = i+j*node_row_num;
                for (int type = 0; type < 3; type++) {
                    for (int t = 0; t < 2; t++) {
                        // stretch
                        int next_x = i+StretchSpringDir[t,0];
                        int next_y = j+StretchSpringDir[t,1];
                        float ks=ksStretch;
                        float length = stretchSpringLen;
                        // shear
                        if (type == 1) {
                            next_x = i+ShearSpringDir[t,0];
                            next_y = j+ShearSpringDir[t,1];
                            ks=ksShear;
                            length = shearSpringLen;
                        }
                        // bend
                        else if (type == 2) {
                            next_x = i+BendSpringDir[t,0];
                            next_y = j+BendSpringDir[t,1];
                            ks=ksBend;
                            length = bendSpringLen;
                        }
                        if (hasThisPos(next_x, next_y)) {
                            int next_index = next_x + next_y*node_row_num;
                            // Hooke's law
                            Vector3 pos_diff = pos[index]-pos[next_index];
                            float distance = Vector3.Distance(pos[next_index],pos[index]);
                            Vector3 pos_dir = pos_diff/distance;
                            // calculate spring force
                            force[index] -= ks*(distance-length)*pos_dir;
                            force[next_index] += ks*(distance-length)*pos_dir;
                            //energy
                            im_energy += ks*(distance-length)*(distance-length)*0.5f;
                            // for implict Jacobian
                            float[] mat = {pos_dir.x*pos_dir.x, pos_dir.x*pos_dir.y, pos_dir.x*pos_dir.z, 
                                           pos_dir.y*pos_dir.x, pos_dir.y*pos_dir.y, pos_dir.y*pos_dir.z, 
                                           pos_dir.z*pos_dir.x, pos_dir.z*pos_dir.y, pos_dir.z*pos_dir.z, }; 
                            for (int n=0; n<9; n++) {
                                float hes = ks*(I[n] - (length/distance)*(I[n]-mat[n]));
                                df[index,index,n] -= hes;
                                df[next_index,next_index,n] -= hes;
                                df[next_index,index,n] += hes;
                                df[index,next_index,n] += hes;
                                // jacobian matrix
                                A[index,index,n] = I[n]-dt*dt*df[index,index,n]; 
                                A[next_index,next_index,n] = I[n]-dt*dt*df[next_index,next_index,n]; 
                                A[index,next_index,n] = -dt*dt*df[index,next_index,n];
                                A[next_index,index,n] = -dt*dt*df[next_index,index,n];
                            }
                        }
                    }
                }

                //gravity
                force[index] += new Vector3(0,-gravity,0) * node_mass;
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
                //damp
                force[index] += - Damp * velocity[index];
            }
        }
        if (upLeftFixed) {force[0]=Vector3.zero;}
        if (upRightFixed) {force[node_row_num-1]=Vector3.zero;}
        if (downLeftFixed) {force[node_row_num*(node_row_num-1)]=Vector3.zero;}
        if (downRightFixed) {force[node_row_num*node_row_num-1]=Vector3.zero;}
        Debug.Log("   im_energy:   " + im_energy);
        WriteData(" " + im_energy);

        for (int i = 0; i < node_num; i++) {
            inertia[i] = node_mass * (2 * pos[i] - pos_pre[i]);
            res[i] = Vector3.zero;
            for (int j = 0; j < node_num; j++) {
                res[i].x += df[i,j,0] * pos[j].x + df[i,j,1] * pos[j].y + df[i,j,2] * pos[j].z;
                res[i].y += df[i,j,3] * pos[j].x + df[i,j,4] * pos[j].y + df[i,j,5] * pos[j].z;
                res[i].z += df[i,j,6] * pos[j].x + df[i,j,7] * pos[j].y + df[i,j,8] * pos[j].z;       
            }
            rhs[i] = inertia[i] + dt*dt*force[i] - dt*dt*res[i];
        }

        for (int i = 0; i < node_num; i++) {
            d[i] = Vector3.zero;
            res2[i] = Vector3.zero;
            Ad[i] = Vector3.zero;
            Vector3 x = Vector3.zero;
            for (int j = 0; j < node_num; j++) {
                x.x += A[i,j,0] * pos[j].x + A[i,j,1] * pos[j].y + A[i,j,2] * pos[j].z;
                x.y += A[i,j,3] * pos[j].x + A[i,j,4] * pos[j].y + A[i,j,5] * pos[j].z;
                x.z += A[i,j,6] * pos[j].x + A[i,j,7] * pos[j].y + A[i,j,8] * pos[j].z; 
            }
            res2[i] = rhs[i] - x;
        }
        int it = 0, it_max = 200;
        float rho = 0, rho_old=1, alpha, beta;
        
        while (it < it_max) {
            it += 1;
            rho = 0;
            alpha = 0;
            beta = 0;
            for (int j=0; j<node_num; j++) {
                rho += res2[j].x * res2[j].x + res2[j].y * res2[j].y + res2[j].z * res2[j].z;
            }
            if (rho < 1e-8) { break;}
            
            beta = rho / rho_old;
            for (int i = 0; i<node_num; i++) {
                d[i] = res2[i] + beta*d[i];
            }

            for (int i = 0; i < node_num; i++) {
                Ad[i] = Vector3.zero;
                for (int j = 0; j < node_num; j++) {
                    Ad[i].x += A[i,j,0] * d[j].x + A[i,j,1] * d[j].y + A[i,j,2] * d[j].z;
                    Ad[i].y += A[i,j,3] * d[j].x + A[i,j,4] * d[j].y + A[i,j,5] * d[j].z;
                    Ad[i].z += A[i,j,6] * d[j].x + A[i,j,7] * d[j].y + A[i,j,8] * d[j].z; 
                }
            }

            for (int j=0; j<node_num; j++) {
                alpha += d[j].x * Ad[j].x + d[j].y * Ad[j].y + d[j].z * Ad[j].z;
            }
            alpha = rho/alpha;

            for (int i=0; i<node_num; i++) {
                res2[i] = res2[i] - alpha * Ad[i];
                pos[i] = pos[i] + alpha * d[i];
            }
            rho_old = rho;
            if (upLeftFixed) {pos[0]=pos_pre[0];}
            if (upRightFixed) {pos[node_row_num-1]=pos_pre[node_row_num-1];}
            if (downLeftFixed) {pos[node_row_num*(node_row_num-1)]=pos_pre[node_row_num*(node_row_num-1)];}
            if (downRightFixed) {pos[node_row_num*node_row_num-1]=pos_pre[node_row_num*node_row_num-1];}
        }
    }

    void Assemble()
    {
        if (method==2) {ComputeJacobianForce();}
        else {ComputeForce();}

        for (int i=0; i<node_num; i++) {
            // pin node
            if (upLeftFixed && i==0) {continue;}
            if (upRightFixed && i==node_row_num-1) {continue;}
            if (downLeftFixed && i==node_row_num*(node_row_num-1)) {continue;}
            if (downRightFixed && i==node_row_num*node_row_num-1) {continue;}

            acceleration[i] = force[i]/node_mass;

            // explict euler method
            if (method==0) {
                velocity[i] = velocity_pre[i] + dt*acceleration[i];
                pos[i] = pos_pre[i] + dt*velocity[i];
            }
            // semi-implict Eluer method
            else if (method==1) {
                velocity[i] = velocity_pre[i] + dt*acceleration[i];
                pos[i] = pos_pre[i] + dt*velocity[i];
            }
            // implict Eluer method
            else if (method==2) {
                // only for wind
                velocity[i] = (pos[i] - pos_pre[i])/dt;
            }
            // Verlet Integration
            else if (method==3) {
                Vector3 x = dt*acceleration[i];
                pos[i] = 2.0f*pos_pre[i]-pos_prepre[i] + dt*x; 
                velocity[i] = velocity_pre[i] + 0.5f*dt*(acceleration[i]+acceleration_pre[i]);
            }

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
