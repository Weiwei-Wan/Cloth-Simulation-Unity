using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEditor;

[CustomEditor(typeof(human))]
public class humanEditor : Editor
{
    public override void OnInspectorGUI()
    {
        human w = (human) target;

        w.method_options = new string[]{"explict Eluer method", "semi-implict Eluer method", "implict Eluer method", "Verlet Integration", "Fast method"};
        w.method = EditorGUILayout.Popup(w.method, w.method_options);

        // w.ksStretch      = EditorGUILayout.FloatField("Ks Stretch: ", w.ksStretch);
        // w.ksShear        = EditorGUILayout.FloatField("Ks Shear: ", w.ksShear);
        // w.ksBend         = EditorGUILayout.FloatField("Ks Bend: ", w.ksBend);

        // w.Damp           = EditorGUILayout.FloatField("Damp Kd: ", w.Damp);
        
        // // force   //slider
        // w.gravity        = EditorGUILayout.Slider("Gravity", w.gravity, 0, 1);
        // w.frontWindForce = EditorGUILayout.Slider("Front Wind Force", w.frontWindForce, 0, 1);
        //w.fronthumanForce = EditorGUILayout.FloatField("Front human Force", w.fronthumanForce);
    }
}
