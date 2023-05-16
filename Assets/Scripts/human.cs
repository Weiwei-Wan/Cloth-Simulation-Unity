using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class human : MonoBehaviour
{
    Mesh mesh;
    MeshFilter meshFilter;
    private GameObject ObjHuman;

    public float gravity;
    public float frontWindForce;
    public float Damp;
    public float Constrain;
    public float ksStretch;
    public float ksShear;
    public float ksBend;
    public int method;
    public string[] method_options;

    const int col = 40;
    const int row = 8;
    const int node_num = row * col;
    const int element_num = row * col * 2;
    int[] element_idx = new int[element_num * 3];

    const float dress_size1 = 1.5f;
    const float dress_size2 = 1.5f;

    float node_mass = 1.0f;
    const int loopNum = 1;
    const int frameNum = 60;
    float dt = 1/((float) loopNum*frameNum);

    Vector2[] uvs = new Vector2[node_num];

    // 3 forces in string-nodes system
    // stretch: string neighbor 2 nodes
    const int stretch_num = row*col*2-col;
    // shear: string diagonal 2 nodes
    const int shear_num = row*col*2-2*col;
    // blend: next interval 1 nodes
    const int bend_num = row*col*2-2*col;
    const int constraint_num = stretch_num + shear_num + bend_num; 
    Vector3[] d_springs = new Vector3[constraint_num]; // vec3(pos_x, pox_y, spring_type)
    
    Vector3[] pos = new Vector3[node_num];
    Vector3[] initial_pos = new Vector3[node_num];
    Vector3[] pos_pre = new Vector3[node_num];
    Vector3[] pos_prepre = new Vector3[node_num];
    Vector3[] force = new Vector3[node_num];
    Vector3[] velocity = new Vector3[node_num];
    Vector3[] velocity_pre = new Vector3[node_num];
    Vector3[] acceleration = new Vector3[node_num];
    Vector3[] acceleration_pre = new Vector3[node_num];
    Vector3[] offset = new Vector3[node_num];
    Vector3 f_ab = Vector3.zero;

    // for implict Jacobian
    float[] I = {1,0,0,0,1,0,0,0,1};
    float[,,] A = new float[node_num, node_num, 9];
    float[,,] df = new float[node_num, node_num, 9];
    Vector3[] b = new Vector3[node_num];

    void Start()
    {
        Debug.Log(dt);
        meshFilter = GetComponent<MeshFilter>();
        mesh = new Mesh();
        meshFilter.mesh = mesh;
        ObjHuman = GameObject.Find("FinalBaseMesh");
        drawTriangle();
        InitConstraint();
    }

    void drawTriangle()
    {
        for (int j = 0; j < row; j++)
        {
            for (int i = 0; i < col; i++)
            {
                int idx = j * col + i;
                float theta = i*2*Mathf.PI/col;
                pos[idx] = pos_pre[idx] = pos_prepre[idx] = initial_pos[idx] = ObjHuman.transform.position + new Vector3((0.5f*j+dress_size1)*Mathf.Cos(theta), 2.2f-j, (0.5f*j+dress_size2)*Mathf.Sin(theta));
                uvs[idx] = new Vector2((float)-Mathf.Cos(theta), (float)-Mathf.Sin(theta));
            }
        }
        print("drawTriangle");
        // generate triangle for the mesh
        int cnt = 0;
        for (int j = 0; j < row-1; j++)
        {
            for (int i = 0; i < col-1; i++)
            {
                element_idx[cnt++] = j * col + i;
                element_idx[cnt++] = j * col + i + 1;
                element_idx[cnt++] = j * col + i + col;

                element_idx[cnt++] = j * col + i + 1;
                element_idx[cnt++] = j * col + i + 1 + col;
                element_idx[cnt++] = j * col + i + col;
            }
            element_idx[cnt++] = j * col + col-1;
            element_idx[cnt++] = j * col;
            element_idx[cnt++] = j * col + 2*col-1;

            element_idx[cnt++] = j * col;
            element_idx[cnt++] = j * col + col;
            element_idx[cnt++] = j * col + 2*col-1;
        }
        // add these two triangles to the mesh
        mesh.vertices = pos;
        mesh.triangles = element_idx;  // int[]
        mesh.uv = uvs;
    }

    void InitConstraint()
    {
        // the inner spring
        int cnt = 0;
        for (int j = 0; j < col; j++)
        {
            for (int i = 0; i < row; i++)
            {
                int idx = i * col + j;
                // stretch
                if (j == 0) {d_springs[cnt++] = new Vector3(idx+col-1, idx, 0);}
                if (j < col - 1) d_springs[cnt++] = new Vector3(idx, idx + 1, 0);
                if (i < row - 1) d_springs[cnt++] = new Vector3(idx, idx + col, 0);
                // bend
                if (j == 0 || j==1) {d_springs[cnt++] = new Vector3(idx+col-2, idx, 1);}
                if (j < col - 2) d_springs[cnt++] = new Vector3(idx, idx + 2, 1); 
                if (i < row - 2) d_springs[cnt++] = new Vector3(idx, idx + col*2, 1); 
                // shear
                if (i < row - 1) {
                    if (j < col-1) {
                        d_springs[cnt++] = new Vector3(idx, idx + 1 + col, 2);
                        d_springs[cnt++] = new Vector3(idx + 1, idx + col, 2);
                    }
                    else {
                        d_springs[cnt++] = new Vector3(idx + 1, idx, 2);
                        d_springs[cnt++] = new Vector3(idx + col, idx-col+1, 2);
                    }
                }
            }
        }
        for (int i = 0; i < node_num; i++) {
            force[i] = Vector3.zero;
            velocity[i] = Vector3.zero;
            velocity_pre[i] = Vector3.zero;
        }
        ExternalForceStep();
        ComputeForce();
        acceleration_pre = acceleration;
    }

    // add the external forces
    void ExternalForceStep()
    {
        for (int i = 0; i < node_num; i++)
        {
            force[i] = new Vector3(0, -gravity, frontWindForce);
            offset[i] = Vector3.zero;
        }
    }

    void ComputeForce()
    {
        for(int ic = 0; ic < constraint_num;ic++)
        {
            // Hooke's law
            int st = (int) d_springs[ic].x;
            int ed = (int) d_springs[ic].y;
            int type = (int) d_springs[ic].z;
            float length = Vector3.Distance(initial_pos[st], initial_pos[ed]);
            Vector3 pos_diff = pos[st] - pos[ed];
            float distance = Vector3.Distance(pos[st], pos[ed]);
            Vector3 pos_dir = pos_diff/distance;
            // calculate spring force
            if (type==0) {f_ab = ksStretch*(distance - length)*pos_dir;}
            if (type==1) {f_ab = ksShear*(distance - length)*pos_dir;}
            if (type==2) {f_ab = ksBend*(distance - length)*pos_dir;}
            //Debug.Log(f_ab);
            force[st] -= f_ab;
            force[ed] += f_ab;
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

        // external spring force
        for(int ic = 0; ic < constraint_num;ic++)
        {
            // Hooke's law
            int st = (int) d_springs[ic].x;
            int ed = (int) d_springs[ic].y;
            int type = (int) d_springs[ic].z;
            float length = Vector3.Distance(initial_pos[st], initial_pos[ed]);
            Vector3 pos_diff = pos[st] - pos[ed];
            float distance = Vector3.Distance(pos[st], pos[ed]);
            Vector3 pos_dir = pos_diff/distance;
            // calculate spring force
            float ks=1;
            if (type==0) {ks = ksStretch;}
            if (type==1) {ks = ksShear;}
            if (type==2) {ks = ksBend;}
            f_ab = ks*(distance - length)*pos_dir;
            //Debug.Log(f_ab);
            force[st] -= f_ab;
            force[ed] += f_ab;

            // for implict Jacobian
            float[] mat = {pos_dir.x*pos_dir.x, pos_dir.x*pos_dir.y, pos_dir.x*pos_dir.z, 
                         pos_dir.y*pos_dir.x, pos_dir.y*pos_dir.y, pos_dir.y*pos_dir.z, 
                         pos_dir.z*pos_dir.x, pos_dir.z*pos_dir.y, pos_dir.z*pos_dir.z, }; 
            for (int t=0; t<9; t++) {
                float hes = ks*(I[t] - (length/distance)*(I[t]-mat[t]));
                df[st,st,t] -= hes;
                df[ed,ed,t] -= hes;
                df[ed,st,t] += hes;
                df[st,ed,t] += hes;
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
            force[i] -= Damp*velocity_pre[i];
            b[i] = velocity_pre[i] + (dt/node_mass) * force[i];
        }

        // jacobian iteration
        for (int iter=0; iter<5; iter++) {
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
        for (int step=0; step < loopNum; step++) {
            ExternalForceStep();
            if (method==2) {ComputeJacobianForce();}
            else {ComputeForce();}
            
            for (int i=0; i<node_num; i++) {
                // add damp 
                if (method<2) {
                    force[i] -= Damp*velocity_pre[i];
                    acceleration[i] = force[i]/node_mass;
                }

                // explict euler method
                if (method==0) {
                    velocity[i] = velocity_pre[i] + dt*acceleration[i];
                    pos[i] = pos_pre[i] + dt*velocity_pre[i];
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
                // Fast method
                else if (method==4) {
                    
                }

                // collision with human
                Vector3 center = new Vector3(ObjHuman.transform.position.x, pos[i].y, ObjHuman.transform.position.z);
                float bound = ((ObjHuman.transform.position.y+2.2f-pos[i].y)/4)+dress_size1;
                if (Vector3.Distance(pos[i], center) < bound) {
                    Vector3 s_pos_dir = (pos[i] - center)/Vector3.Distance(pos[i], center);
                    pos[i] = center + bound*s_pos_dir;
                }
                // constrain 
                else if (Vector3.Distance(pos[i], pos_pre[i]) > 0.2f) {
                    Vector3 s_pos_dir = (pos[i] - pos_pre[i])/Vector3.Distance(pos[i], pos_pre[i]);
                    pos[i] = pos_pre[i] + 0.2f*s_pos_dir;
                    print("constrain happen");
                }

                if (i<col) {pos[i]=initial_pos[i];}

                // save the pre
                pos_prepre[i] = pos_pre[i];
                pos_pre[i] = pos[i];
                velocity_pre[i] = velocity[i];
                acceleration_pre[i] = acceleration[i];
            } 

            // fixed node
            
        }
    }

    void UpdateMesh()
    {
        mesh.vertices = pos;
        mesh.RecalculateNormals();
    }

    //set up the frame number each second default=60
    void Awake()
    {
        Application.targetFrameRate = frameNum;
    }

    private void Update()
    {
        if (Input.GetKeyDown("a") || Input.GetKeyDown("d")) {
            if (Input.GetKeyDown("a")) {
                ObjHuman.transform.position = new Vector3(ObjHuman.transform.position.x-0.1f, ObjHuman.transform.position.y, ObjHuman.transform.position.z);
            }
            if (Input.GetKeyDown("d")) {
                ObjHuman.transform.position = new Vector3(ObjHuman.transform.position.x+0.1f, ObjHuman.transform.position.y, ObjHuman.transform.position.z);
            }
            for (int j = 0; j < row; j++)
            {
                for (int i = 0; i < col; i++)
                {
                    int idx = j * col + i;
                    float theta = i*2*Mathf.PI/col;
                    initial_pos[idx] = ObjHuman.transform.position + new Vector3((0.5f*j+dress_size1)*Mathf.Cos(theta), 2.2f-j, (0.5f*j+dress_size2)*Mathf.Sin(theta));
                    Vector3 center =  new Vector3(ObjHuman.transform.position.x, pos[idx].y, ObjHuman.transform.position.z);
                    float bound = ((ObjHuman.transform.position.y+2.2f-pos[idx].y)/4)+dress_size1;
                    if (Vector3.Distance(pos[idx], center) < bound) {
                        Vector3 s_pos_dir = (pos[idx] - center)/Vector3.Distance(pos[i], center);
                        pos[idx] = center + bound*s_pos_dir;
                    }

                    if (idx<col) {pos[idx]=initial_pos[idx];}
                }
            }
        }
        Assemble();
        UpdateMesh();
    }
}
