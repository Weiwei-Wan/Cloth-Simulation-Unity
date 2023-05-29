using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Wind : MonoBehaviour
{
    Mesh mesh;
    MeshFilter meshFilter;
    private GameObject sphere;
    const int radius = 1;

    public int method;
    public string[] method_options;

    bool upLeftFixed = true;
    bool upRightFixed = true;
    bool downLeftFixed = false;
    bool downRightFixed = false;

    float gravity = 9.8f;
    float frontWindForce = 5;
    float Damp = 1;

    float ksStretch = 10000;
    float ksShear = 10000;
    float ksBend = 10000;

    const int node_row_num = 16;
    const int node_num = node_row_num * node_row_num;
    const int element_row_num = node_row_num - 1;
    const int element_num = element_row_num * element_row_num * 2;
    int[] element_idx = new int[element_num * 3];

    float springIniLen = 0.1f;
    float stretchSpringLen;
    float shearSpringLen;
    float bendSpringLen;
    float node_mass = 1;
    float dt = 0.002f;

    // 3 forces in string-nodes system
    int[,] StretchSpringDir = {{1,0},{0,1},{-1,0},{0,-1}};
    int[,] ShearSpringDir = {{2,0},{0,2},{-2,0},{0,-2}};
    int[,] BendSpringDir = {{1,1},{-1,1},{-1,-1},{1,-1}};
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
    Vector3[] b = new Vector3[node_num];

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

    void ComputeForce()
    {
        for (int i=0; i<node_row_num; i++) {
            for (int j=0; j<node_row_num; j++) {
                int index = i+j*node_row_num;
                force[index] = new Vector3(0,0,0);
                // stretch
                for (int t = 0; t < 4; t++) {
                    int next_x = i+StretchSpringDir[t,0];
                    int next_y = j+StretchSpringDir[t,1];
                    if (hasThisPos(next_x, next_y)) {
                        int next_index = next_x + next_y*node_row_num;
                        // Hooke's law
                        Vector3 pos_diff = pos[next_index]-pos[index];
                        float distance = Vector3.Distance(pos[next_index],pos[index]);
                        Vector3 pos_dir = pos_diff/distance;
                        // calculate spring force
                        force[index] += ksStretch*(distance-stretchSpringLen)*pos_dir;
                    }
                }
                // shear
                for (int t = 0; t < 4; t++) {
                    int next_x = i+ShearSpringDir[t,0];
                    int next_y = j+ShearSpringDir[t,1];
                    if (hasThisPos(next_x, next_y)) {
                        int next_index = next_x + next_y*node_row_num;
                        // Hooke's law
                        Vector3 pos_diff = pos[next_index]-pos[index];
                        float distance = Vector3.Distance(pos[next_index],pos[index]);
                        Vector3 pos_dir = pos_diff/distance;
                        // calculate spring force
                        force[index] += ksShear*(distance-shearSpringLen)*pos_dir;
                    }
                }
                // bend
                for (int t = 0; t < 4; t++) {
                    int next_x = i+BendSpringDir[t,0];
                    int next_y = j+BendSpringDir[t,1];
                    if (hasThisPos(next_x, next_y)) {
                        int next_index = next_x + next_y*node_row_num;
                        // Hooke's law
                        Vector3 pos_diff = pos[next_index]-pos[index];
                        float distance =  Vector3.Distance(pos[next_index],pos[index]);
                        Vector3 pos_dir = pos_diff/distance;
                        // calculate spring force
                        force[index] += ksBend*(distance-bendSpringLen)*pos_dir;
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
    }
    
    float HlsCalculate(float[] x, int row){
        float hls = 0;
        float[] c = new float[(row-1)*(row-1)];
        if(row == 1){
            return x[0];
        }
        if(row == 2){
            float kk = x[0]*x[3] - x[1]*x[2];
            return kk;
        }
        int i, j, k;
        for(j = 0;j < row ;j++){
            for(i = 0 ; i < row - 1;i++){
                for(k = 0 ; k < row - 1;k++){
                    if(k<j) c[i*(row-1)+k] = x[(i+1)*row+k];
                    if(k>=j) c[i*(row-1)+k] = x[(i+1)*row+k+1];
                }
            }
            hls += Mathf.Pow(-1,j)*x[j]*HlsCalculate(c,row-1);
        }
        return hls;
    }

    float[] adjoint_Matrix(float[] a, int row) {
        int n=0, m=0, nn = 0, mm = 0;
        int i, j;
        float[] company = new float[row*row];
        for( i = 0; i < row ;i++) {
            for(j = 0 ; j < row ; j++) {
                float[] tempArr = new float[(row-1)*(row-1)];
                n = 0;
                m = 0;
                for(int p = 0; p < row ;p++) {
                    for(int q = 0 ;q < row ;q++){
                        if(!(p == i || q ==j)) {
                            tempArr[n*(row-1)+m] = a[p*row+q];
                            m++;
                            if(m == row - 1) {
                                m = 0;
                                n++;
                            }
                        }
                    }
                }	
                float k =  HlsCalculate(tempArr,row-1);
                company[nn*row+mm] = Mathf.Pow(-1,i+j) * HlsCalculate(tempArr,row-1);
                nn++;
                if (nn == row ) {
                    nn = 0;
                    mm++;
                }
            }
        }
        return company;
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

        // internal spring force
        for (int i=0; i<node_row_num; i++) {
            for (int j=0; j<node_row_num; j++) {
                // stretch
                int index = i+j*node_row_num;
                for (int t = 0; t < 4; t++) {
                    int next_x = i+StretchSpringDir[t,0];
                    int next_y = j+StretchSpringDir[t,1];
                    if (hasThisPos(next_x, next_y)) {
                        int next_index = next_x + next_y*node_row_num;
                        // Hooke's law
                        Vector3 pos_diff = pos[next_index]-pos[index];
                        float distance = Vector3.Distance(pos[next_index],pos[index]);
                        Vector3 pos_dir = pos_diff/distance;
                        float ks=ksStretch;
                        float length = stretchSpringLen;
                        // calculate spring force
                        force[index] += ks*(distance-length)*pos_dir;
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
                        }
                    }
                }
                // shear
                for (int t = 0; t < 4; t++) {
                    int next_x = i+ShearSpringDir[t,0];
                    int next_y = j+ShearSpringDir[t,1];
                    if (hasThisPos(next_x, next_y)) {
                        int next_index = next_x + next_y*node_row_num;
                        // Hooke's law
                        Vector3 pos_diff = pos[next_index]-pos[index];
                        float distance = Vector3.Distance(pos[next_index],pos[index]);
                        Vector3 pos_dir = pos_diff/distance;
                        float ks=ksShear;
                        float length = shearSpringLen;
                        // calculate spring force
                        force[index] += ks*(distance-length)*pos_dir;
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
                        }
                    }
                }
                // bend
                for (int t = 0; t < 4; t++) {
                    int next_x = i+BendSpringDir[t,0];
                    int next_y = j+BendSpringDir[t,1];
                    if (hasThisPos(next_x, next_y)) {
                        int next_index = next_x + next_y*node_row_num;
                        // Hooke's law
                        Vector3 pos_diff = pos[next_index]-pos[index];
                        float distance =  Vector3.Distance(pos[next_index],pos[index]);
                        Vector3 pos_dir = pos_diff/distance;
                        float ks=ksBend;
                        float length = bendSpringLen;
                        // calculate spring force
                        force[index] += ks*(distance-length)*pos_dir;
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
        
        // jacobian matrix
        for (int i=0; i<node_num; i++) {
            for (int j=0; j<node_num; j++) {
                for (int t=0; t<9; t++) {
                    if (i != j) {A[i,j,t] = -dt*dt*df[i,j,t];}
                    else {A[i,j,t] = I[t]-dt*dt*df[i,j,t];} 
                }
            }
        }

        // b
        for (int i=0; i<node_num; i++) {
            b[i] = velocity_pre[i] + (dt/node_mass) * force[i];
        }

        // jacobian iteration
        for (int iter=0; iter<10; iter++) {
            for (int i=0; i<node_num; i++) {
                Vector3 r = b[i];
                for (int j=0; j<node_num; j++) {
                    if (i != j) {
                        float x = A[i,j,0]*velocity_pre[j].x + A[i,j,1]*velocity_pre[j].y + A[i,j,2]*velocity_pre[j].z;
                        float y = A[i,j,3]*velocity_pre[j].x + A[i,j,4]*velocity_pre[j].y + A[i,j,5]*velocity_pre[j].z;
                        float z = A[i,j,6]*velocity_pre[j].x + A[i,j,7]*velocity_pre[j].y + A[i,j,8]*velocity_pre[j].z;
                        r = new Vector3(r.x-x, r.y-y, r.z-z);
                    }
                }

                float[] temp = {A[i,i,0], A[i,i,1], A[i,i,2], A[i,i,3], A[i,i,4], A[i,i,5], A[i,i,6], A[i,i,7], A[i,i,8],};
                float[] inverse = adjoint_Matrix(temp, 3);
                float hls = 1.0f/HlsCalculate(temp, 3);
                float xx = hls * (inverse[0]*r.x + inverse[1]*r.y + inverse[2]*r.z);
                float yy = hls * (inverse[3]*r.x + inverse[4]*r.y + inverse[5]*r.z);
                float zz = hls * (inverse[6]*r.x + inverse[7]*r.y + inverse[8]*r.z);
                velocity[i] = new Vector3(xx, yy, zz);
            }

            for (int i=0; i<node_num; i++) {
                velocity_pre[i] = velocity[i];
            }
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
                pos[i] = pos_pre[i] + dt*velocity[i];
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
