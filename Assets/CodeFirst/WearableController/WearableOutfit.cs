using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[Serializable]
public class WearableOutfit : ScriptableObject
{
    #region Fields

    public Guid Id { get { return _id; } }
    [SerializeField]
    private Guid _id;

    [SerializeField]
    public List<Guid> Wearables;

    [SerializeField]
    public string Name;

    #endregion Fields

    #region Constructor

    public static WearableOutfit CreateWearableOutfit()
    {
        return ScriptableObject.CreateInstance<WearableOutfit>();
    }

    private WearableOutfit()
    {
        _id = Guid.NewGuid();
        Wearables = new List<Guid>();
    }

    #endregion Constructor

    #region Methods

    public void AddWearable(Guid wearableId)
    {
        if (Wearables == null) Wearables = new List<Guid>();
        if(!Wearables.Exists(g=>g == wearableId))
        {
            Wearables.Add(wearableId);
        }        
    }

    public void RemoveWearable(Guid wearableId)
    {
        if (Wearables != null)
        {
            Wearables.RemoveAll(g => g == wearableId);
        }
    }

    #endregion Methods
}
