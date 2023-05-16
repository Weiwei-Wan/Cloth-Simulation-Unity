using System;
using System.Collections;
using System.Collections.Generic;

using UnityEngine;
using UnityEditor;
using UnityEditor.SceneManagement;

[CustomEditor(typeof(CharacterWearableController))]
[CanEditMultipleObjects]
public class WearableInspector : Editor
{
    string wearableName = "";
    string characterFolderName = "";
    WearableType wearableType = WearableType.Top;
    Vector2 scrollPos;

    public override void OnInspectorGUI()
    {
        //base.OnInspectorGUI();

        serializedObject.Update();
        CharacterWearableController controller = (CharacterWearableController)target;
        EditorGUILayout.LabelField("", GUI.skin.horizontalSlider);

        #region Add New Wearable

        GUILayout.Label("Add New Wearable");

        GUILayout.BeginHorizontal();

        GUILayout.Label("Wearable Name: ");
        wearableName = GUILayout.TextField(wearableName, 100);

        GUILayout.EndHorizontal();

        GUILayout.BeginHorizontal();

        GUILayout.Label("Character Folder: ");
        characterFolderName = GUILayout.TextField(characterFolderName, 100);

        GUILayout.EndHorizontal();

        GUILayout.BeginHorizontal();

        //GUILayout.Label("Wearable Type: ");
        wearableType = (WearableType)EditorGUILayout.EnumPopup("Wearable Type: ", wearableType);

        GUILayout.EndHorizontal();

        //Add new wearable button
        if (GUILayout.Button("Add Wearable"))
        {
            Wearable wearable = Wearable.CreateWearable();
            wearable.Name = wearableName;
            wearable.CharacterFolderName = characterFolderName;
            wearable.WearableType = wearableType;
            controller.AddWearable(wearable);
        }

        #endregion Add New Wearable

        EditorGUILayout.LabelField("", GUI.skin.horizontalSlider);

        #region Remove All Wearables

        if (GUILayout.Button("Remove All Wearables"))
        {
            controller.RemoveAllWearables(true, false);
        }
        if (GUILayout.Button("Remove All Wearables(Including Hair)"))
        {
            controller.RemoveAllWearables(true, true);
        }
        EditorGUILayout.LabelField("", GUI.skin.horizontalSlider);

        #endregion Remove All Wearables

        #region Added Wearables List

        GUILayout.Label("Added Wearables");
        EditorGUILayout.Space();

        scrollPos = EditorGUILayout.BeginScrollView(scrollPos, GUILayout.Height(300));

        List<Wearable> wearables = controller.Wearables;
        for(int j=0;j<wearables.Count;j++)
        {
            Wearable wearable = wearables[j];
            GUILayout.Label(String.Format("Name: {0}", wearable.Name));
            GUILayout.Label(String.Format("Character Folder: {0}", wearable.CharacterFolderName));
            GUILayout.Label(String.Format("Type: {0}", wearable.WearableType));
            GUILayout.BeginHorizontal();

            if (GUILayout.Button("Apply"))
            {
                controller.ApplyWearable(wearable, true);
            }

            if (GUILayout.Button("Remove"))
            {
                controller.RemoveWearable(wearable.WearableType, true);
            }

            GUILayout.EndHorizontal();

            if (GUILayout.Button("Delete"))
            {
                controller.DeleteWearable(wearable.Id, true);
            }

            EditorGUILayout.Space();
            EditorGUILayout.Space();
        }

        EditorGUILayout.EndScrollView();

        #endregion Added Wearables List

        if (GUI.changed)
        {
            EditorSceneManager.MarkSceneDirty(controller.gameObject.scene);
        }

        serializedObject.ApplyModifiedProperties();        
    }
}
