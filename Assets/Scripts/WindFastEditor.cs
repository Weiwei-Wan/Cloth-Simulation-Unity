using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEditor;

[CustomEditor(typeof(WindFast))]
public class WindFastEditor : Editor
{
    public override void OnInspectorGUI()
    {
        WindFast w = (WindFast)target;
        w.upLeftFixed    = EditorGUILayout.Toggle("Up Left Fixed",    w.upLeftFixed);
        w.upRightFixed   = EditorGUILayout.Toggle("Up Right Fixed",   w.upRightFixed);
        w.downLeftFixed  = EditorGUILayout.Toggle("Down Left Fixed",  w.downLeftFixed);
        w.downRightFixed = EditorGUILayout.Toggle("Down Right Fixed", w.downRightFixed);

        w.ksStretch      = EditorGUILayout.FloatField("Ks Stretch: ", w.ksStretch);
        w.ksShear        = EditorGUILayout.FloatField("Ks Shear: ", w.ksShear);
        w.ksBend         = EditorGUILayout.FloatField("Ks Bend: ", w.ksBend);
        
        // force
        w.gravity = EditorGUILayout.FloatField("Gravity", w.gravity);
        w.frontWindForce = EditorGUILayout.FloatField("Front Wind Force", w.frontWindForce);
    }
}
