using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEditor;

[CustomEditor(typeof(Wind))]
public class WindEditor : Editor
{
    public override void OnInspectorGUI()
    {
        Wind w = (Wind) target;

        w.method_options = new string[]{"explict Eluer method", "semi-implict Eluer method", "implict Eluer method", "Verlet Integration"};
        w.method = EditorGUILayout.Popup(w.method, w.method_options);
        
        w.upLeftFixed    = EditorGUILayout.Toggle("Up Left Fixed",    w.upLeftFixed);
        w.upRightFixed   = EditorGUILayout.Toggle("Up Right Fixed",   w.upRightFixed);
        w.downLeftFixed  = EditorGUILayout.Toggle("Down Left Fixed",  w.downLeftFixed);
        w.downRightFixed = EditorGUILayout.Toggle("Down Right Fixed", w.downRightFixed);

        w.ksStretch      = EditorGUILayout.FloatField("Ks Stretch: ", w.ksStretch);
        w.ksShear        = EditorGUILayout.FloatField("Ks Shear: ", w.ksShear);
        w.ksBend         = EditorGUILayout.FloatField("Ks Bend: ", w.ksBend);

        w.Damp           = EditorGUILayout.FloatField("Damp Kd: ", w.Damp);
        //w.Constrain      = EditorGUILayout.FloatField("Constrain value: ", w.Constrain);

        // force   //slider
        w.gravity        = EditorGUILayout.Slider("Gravity", w.gravity, 0, 1);
        w.frontWindForce = EditorGUILayout.Slider("Front Wind Force", w.frontWindForce, 0, 1);
        //w.frontWindForce = EditorGUILayout.FloatField("Front Wind Force", w.frontWindForce);
    }
}
