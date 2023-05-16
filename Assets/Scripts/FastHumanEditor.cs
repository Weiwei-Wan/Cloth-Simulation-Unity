using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEditor;

[CustomEditor(typeof(FastHuman))]
public class FastHumanEditor : Editor
{
    public override void OnInspectorGUI()
    {
        FastHuman w = (FastHuman) target;

        w.ksStretch      = EditorGUILayout.FloatField("Ks Stretch: ", w.ksStretch);
        w.ksShear        = EditorGUILayout.FloatField("Ks Shear: ", w.ksShear);
        w.ksBend         = EditorGUILayout.FloatField("Ks Bend: ", w.ksBend);

        // force   //slider
        w.gravity        = EditorGUILayout.Slider("Gravity", w.gravity, 0, 1);
        w.frontWindForce = EditorGUILayout.Slider("Front Wind Force", w.frontWindForce, 0, 1);
    }
}
