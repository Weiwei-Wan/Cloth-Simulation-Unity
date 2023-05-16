using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

#region Wearable Type Enum

public enum WearableType
{
    //Some sample wearables.
    Hair = 1,
    Top = 2,
    Bottom = 3,
    Shoes = 4,
    Jacket = 5,
    Accessory1 = 6,
    Accessory2 = 7,
    Glasses = 8,
    Shoe1 = 9,
    Shoe2 = 10
}

#endregion Wearable Type Enum

#region Wearable Class

[Serializable]
public class Wearable : ScriptableObject
{
    #region Properties and Fields

    public Guid Id { get { return _id; } }

    [SerializeField]
    private Guid _id;

    [SerializeField]
    public GameObject AttachedGameObject { get { return _attachedGameObject; } }

    [SerializeField]
    private GameObject _attachedGameObject;

    public bool IsActive { get { return _isActive; } }

    [SerializeField]
    private bool _isActive;

    [SerializeField]
    public string Name;

    [SerializeField]
    public string CharacterFolderName;

    [SerializeField]
    public WearableType WearableType;

    #endregion Properties and Fields

    #region Constructor

    public static Wearable CreateWearable()
    {
        return ScriptableObject.CreateInstance<Wearable>();
    }

    private Wearable()
    {
        _id = Guid.NewGuid();
        
    }

    #endregion Constructor

    #region Methods

    /// <summary>
    /// Attaches a game object loaded for this wearable.
    /// </summary>
    /// <param name="gameObject"></param>
    public void AttachGameObject(GameObject gameObject)
    {
        _attachedGameObject = gameObject;
        _isActive = true;
    }

    /// <summary>
    /// Removes the attached game object loaded for this wearable.
    /// </summary>
    public void RemoveGameObject()
    {
        _attachedGameObject = null;
        _isActive = false;
    }

    #endregion Methods
}

#endregion Wearable Class