using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class RayTriangleCollision
{
    /// <summary>
    /// 射线与三角形相交检测
    /// </summary>
    /// <param name="source">射线起点坐标</param>
    /// <param name="rayDirection">射线方向以及长度</param>
    /// <param name="p0">三角形一个顶点</param>
    /// <param name="p1">三角形一个顶点</param>
    /// <param name="p2">三角形一个顶点，三角形三个顶点是有顺序要求的，因为输入的顶点次序(顺时针、逆时针)决定三角形法线的朝向</param>
    /// <param name="hitPos">相交情况下交点坐标</param>
    /// <returns>返回值为 true 为相交、false 为不相交</returns>
    public bool IsRayTriangleCollision(Vector3 source, Vector3 rayDirection, Vector3 p0, Vector3 p1, Vector3 p2)
    {
        Vector3 E1 = p1 - p0;
        Vector3 E2 = p2 - p1;
        // 三角形所在平面法向量
        Vector3 triangleNormal = Vector3.Cross(E1, E2);
 
        // 射线起点到三角形一个顶点的向量
        Vector3 PC = p0 - source;
        float dot_rayDir_planeNormal = Dot(rayDirection, triangleNormal);
        float dot_pa_planeNormal = Dot(PC, triangleNormal);
 
        if (   dot_rayDir_planeNormal == 0 
            || (dot_rayDir_planeNormal > 0 && dot_pa_planeNormal < 0) 
            || (dot_rayDir_planeNormal < 0 && dot_pa_planeNormal > 0))
        {
            return false;
        }
 
        // 向量长度不足以到达三角形所在的平面
        if (Mathf.Abs(dot_rayDir_planeNormal) < Mathf.Abs(dot_pa_planeNormal))
        {
            return false;
        }
 
        // 计算射线到平面的距离
        float length = dot_pa_planeNormal / dot_rayDir_planeNormal;
        // 计算射线与平面交点坐标
        Vector3 hitPos = source + rayDirection * length;
 
        float u0, u1, u2;
        float v0, v1, v2;
        // 找到主要的轴，选择投影平面
        if (Mathf.Abs(triangleNormal.x) > Mathf.Abs(triangleNormal.y))
        {
            if (Mathf.Abs(triangleNormal.x) > Mathf.Abs(triangleNormal.z))
            {
                u0 = hitPos.y - p0.y;
                u1 = p1.y - p0.y;
                u2 = p2.y - p0.y;
                v0 = hitPos.z - p0.z;
                v1 = p1.z - p0.z;
                v2 = p2.z - p0.z;
            }
            else
            {
                u0 = hitPos.x - p0.x;
                u1 = p1.x - p0.x;
                u2 = p2.x - p0.x;
                v0 = hitPos.y - p0.y;
                v1 = p1.y - p0.y;
                v2 = p2.y - p0.y;
            }
        }
        else if (Mathf.Abs(triangleNormal.y) > Mathf.Abs(triangleNormal.z))
        {
            u0 = hitPos.x - p0.x;
            u1 = p1.x - p0.x;
            u2 = p2.x - p0.x;
            v0 = hitPos.z - p0.z;
            v1 = p1.z - p0.z;
            v2 = p2.z - p0.z;
        }
        else
        {
            u0 = hitPos.x - p0.x;
            u1 = p1.x - p0.x;
            u2 = p2.x - p0.x;
            v0 = hitPos.y - p0.y;
            v1 = p1.y - p0.y;
            v2 = p2.y - p0.y;
        }
 
        // 计算分母，检查其有效性
        float temp = u1 * v2 - v1 * u2;
        if (temp == 0)
        {
            return false;
        }
        temp = 1.0f / temp;
        //计算重心坐标，每一步都检测边界条件
        float alpha = (u0 * v2 - v0 * u2) * temp;
        if (alpha < 0)
        {
            return false;
        }
 
        float beta = (u1 * v0 - v1 * u0) * temp;
        if (beta < 0)
        {
            return false;
        }
 
        float gamma = 1.0f - alpha - beta;
        if (gamma < 0)
        {
            return false;
        }
 
        return true;
    }

    public Vector3 GetHitPos(Vector3 source, Vector3 rayDirection, Vector3 p0, Vector3 p1, Vector3 p2) {
        Vector3 E1 = p1 - p0;
        Vector3 E2 = p2 - p1;
        // 三角形所在平面法向量
        Vector3 triangleNormal = Vector3.Cross(E1, E2);
 
        // 射线起点到三角形一个顶点的向量
        Vector3 PC = p0 - source;
        float dot_rayDir_planeNormal = Dot(rayDirection, triangleNormal);
        float dot_pa_planeNormal = Dot(PC, triangleNormal);
 
        // 计算射线到平面的距离
        float length = dot_pa_planeNormal / dot_rayDir_planeNormal;
        // 计算射线与平面交点坐标
        Vector3 hitPos = source + rayDirection * (length-0.01f);
        return hitPos;
    }
 
    public float Dot(Vector3 vector1, Vector3 vector2)
    {
        return vector1.x * vector2.x + vector1.y * vector2.y + vector1.z * vector2.z;
    }
}