using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class RayTriangleCollision
{
    /// <summary>
    /// </summary>
    /// <param name="source"> </param>
    /// <param name="rayDirection"> </param>
    /// <param name="p0"> </param>
    /// <param name="p1"> </param>
    /// <param name="p2"> </param>
    /// <param name="hitPos"> </param>
    /// <returns> </returns>
    public bool IsRayTriangleCollision(Vector3 source, Vector3 rayDirection, Vector3 p0, Vector3 p1, Vector3 p2)
    {
        Vector3 E1 = p1 - p0;
        Vector3 E2 = p2 - p1;

        Vector3 triangleNormal = Vector3.Cross(E1, E2);
 
        Vector3 PC = p0 - source;
        float dot_rayDir_planeNormal = Dot(rayDirection, triangleNormal);
        float dot_pa_planeNormal = Dot(PC, triangleNormal);
 
        if (   dot_rayDir_planeNormal == 0 
            || (dot_rayDir_planeNormal > 0 && dot_pa_planeNormal < 0) 
            || (dot_rayDir_planeNormal < 0 && dot_pa_planeNormal > 0))
        {
            return false;
        }
 
        
        if (Mathf.Abs(dot_rayDir_planeNormal) < Mathf.Abs(dot_pa_planeNormal))
        {
            return false;
        }
 
        
        float length = dot_pa_planeNormal / dot_rayDir_planeNormal;
      
        Vector3 hitPos = source + rayDirection * length;
 
        float u0, u1, u2;
        float v0, v1, v2;
        
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
 
        float temp = u1 * v2 - v1 * u2;
        if (temp == 0)
        {
            return false;
        }
        temp = 1.0f / temp;
  
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
       
        Vector3 triangleNormal = Vector3.Cross(E1, E2);
 
        Vector3 PC = p0 - source;
        float dot_rayDir_planeNormal = Dot(rayDirection, triangleNormal);
        float dot_pa_planeNormal = Dot(PC, triangleNormal);
 
        float length = dot_pa_planeNormal / dot_rayDir_planeNormal;

        Vector3 hitPos = source + rayDirection * (length-0.01f);
        return hitPos;
    }
 
    public float Dot(Vector3 vector1, Vector3 vector2)
    {
        return vector1.x * vector2.x + vector1.y * vector2.y + vector1.z * vector2.z;
    }
}