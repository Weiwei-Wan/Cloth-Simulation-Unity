using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[Serializable]
public class CharacterWearableController : MonoBehaviour
{
    #region Fields and Properties

    private SkinnedMeshRenderer _personMeshRenderer;
    private SkinnedMeshRenderer _wearableMeshRenderer;
    private Dictionary<string, Transform> _personBoneMap;

    public List<Wearable> Wearables { get { return _wearables; } }
    [SerializeField, HideInInspector]
    private List<Wearable> _wearables = new List<Wearable>();

    public List<WearableOutfit> Outfits { get { return _outfits; } }
    [SerializeField, HideInInspector]
    private List<WearableOutfit> _outfits = new List<WearableOutfit>();

    [SerializeField, HideInInspector]
    public List<Wearable> DefaultWearables = new List<Wearable>();

    #endregion Fields and Properties

    #region Mesh Renderer Methods

    /// <summary>
    /// Method to set person mesh renderer.
    /// </summary>
    /// <param name="personMeshRenderer"></param>
    public void SetPersonMeshRenderer(SkinnedMeshRenderer personMeshRenderer)
    {
        _personMeshRenderer = personMeshRenderer;
        _personBoneMap = new Dictionary<string, Transform>();
        foreach (Transform bone in _personMeshRenderer.bones)
        {
            _personBoneMap[bone.gameObject.name] = bone;
        }
    }

    /// <summary>
    /// Loads the person mesh renderer.
    /// </summary>
    private void LoadPersonMeshRenderer()
    {
        Helper helper = new Helper();
        SkinnedMeshRenderer mr = helper.GetSkinMeshRecursive(transform);
        CharacterWearableController wearableController = transform.GetComponent<CharacterWearableController>();

        if (mr != null && mr.ToString() != "null")
        {
            SetPersonMeshRenderer(mr);
        }
    }

    #endregion Mesh Renderer Methods

    #region Outfits

    /// <summary>
    /// Adds a new outfit.
    /// </summary>
    /// <param name="outfit"></param>
    public void AddOutfit(WearableOutfit outfit)
    {
        if (!_outfits.Exists(o => o.Id == outfit.Id))
        {
            _outfits.Add(outfit);
        }
    }

    /// <summary>
    /// Removes an existing outfit.
    /// </summary>
    /// <param name="outfitId">Outfit Id</param>
    public void RemoveOutfit(Guid outfitId)
    {
        _outfits.RemoveAll(o => o.Id == outfitId);
    }

    /// <summary>
    /// Removes an existing outfit.
    /// </summary>
    /// <param name="outfitName">Outfit Name.</param>
    public void RemoveOutfit(string outfitName)
    {
        _outfits.RemoveAll(o => o.Name == outfitName);
    }

    /// <summary>
    /// Applies all the wearables which are contained inside an outfit.
    /// </summary>
    /// <param name="outfitName"></param>
    public void ApplyOutfit(string outfitName)
    {
        WearableOutfit outfit = _outfits.Find(o => o.Name == outfitName);
        if (outfit != null)
        {
            //First remove all existing wearables.
            RemoveAllWearables(false, false);

            foreach (Guid wearableId in outfit.Wearables)
            {
                Wearable wearable = _wearables.Find(w => w.Id == wearableId);
                if (wearable != null)
                {
                    ApplyWearable(wearable);
                }
            }
        }
    }

    #endregion Outfits

    #region Events

    // Start is called before the first frame update
    void Start()
    {
        foreach (Wearable w in DefaultWearables)
        {
            ApplyWearable(w);
        }
    }

    // Update is called once per frame
    void Update()
    {

    }

    #endregion Events

    #region Wearable Methods

    /// <summary>
    /// Sets the wearable list.
    /// </summary>
    /// <param name="wearables"></param>
    public void SetWearable(List<Wearable> wearables)
    {
        _wearables = wearables;
    }

    /// <summary>
    /// Adds a new wearable.
    /// </summary>
    /// <param name="wearable">Wearable object.</param>
    public void AddWearable(Wearable wearable)
    {
        if (wearable == null) return;
        if (_wearables == null) _wearables = new List<Wearable>();

        if(!_wearables.Exists(w=>w.Id == wearable.Id || w.Name == wearable.Name))
        {
            _wearables.Add(wearable);
        }            
    }

    /// <summary>
    /// Deletes a wearable.
    /// </summary>
    /// <param name="id">Wearable Id</param>
    public void DeleteWearable(Guid id)
    {
        DeleteWearable(id, false);
    }

    /// <summary>
    /// Deletes a wearable.
    /// </summary>
    /// <param name="id">Wearable Id</param>
    /// <param name="editorMode">Is in editor mode or not.</param>
    public void DeleteWearable(Guid id, bool editorMode)
    {
        Wearable wearable = _wearables.Find(w => w.Id == id);
        if(wearable != null)
        {
            RemoveWearable(wearable.WearableType, editorMode);
            _wearables.RemoveAll(w => w.Id == id);
        }        
    }

    /// <summary>
    /// Applies the wearable on the character.
    /// </summary>
    /// <param name="wearable">Wearable Object</param>
    public void ApplyWearable(Wearable wearable)
    {
        ApplyWearable(wearable, false);
    }

    /// <summary>
    /// Applies the wearable on the character.
    /// </summary>
    /// <param name="wearable">Wearable Object</param>
    /// <param name="editorMode">Is Editor Mode or not.</param>
    public void ApplyWearable(Wearable wearable, bool editorMode)
    {
        ApplyWearable(wearable.WearableType, wearable.Name, editorMode);
    }

    /// <summary>
    /// Applies the wearable on the character.
    /// </summary>
    /// <param name="wearableType">Type of wearable.</param>
    /// <param name="wearableName">Wearable Name</param>
    public void ApplyWearable(WearableType wearableType, string wearableName)
    {
        ApplyWearable(wearableType, wearableName, false);
    }

    /// <summary>
    /// Applies the wearable on the character.
    /// </summary>
    /// <param name="wearableType">Type of wearable.</param>
    /// <param name="wearableName">Wearable Name</param>
    /// <param name="editorMode">Is Editor Mode or not.</param>
    public void ApplyWearable(WearableType wearableType, string wearableName, bool editorMode)
    {
        //Check if the person mesh renderer is initialized or not.
        //If not initialized then lazy load it.
        if(_personMeshRenderer == null || _personMeshRenderer.ToString() == "null" || _personBoneMap == null)
        {
            LoadPersonMeshRenderer();
        }

        if(_wearables != null && _wearables.Exists(w=>w.Name == wearableName))
        {
            Wearable wearable = _wearables.Find(w => w.Name == wearableName);
            if(wearable != null)
            {
                //Load the wearable gameobject first and cache it.
                GameObject g = Resources.Load(string.Format("Characters/{0}/Wearables/{1}", wearable.CharacterFolderName, wearable.Name), typeof(GameObject)) as GameObject;
                if(g != null)
                {
                    g = Instantiate(g);
                    
                    g.transform.parent = transform;
                    g.transform.localPosition = new Vector3(0, 0, 0);
                    g.transform.localRotation = Quaternion.identity;

                    //Set the wearable mesh renderer
                    Helper helper = new Helper();
                    _wearableMeshRenderer = helper.GetSkinMeshRecursive(g.transform);
                    
                    g.SetActive(false);
                }

                //Remove the existing wearable
                RemoveWearable(wearableType, editorMode);

                //Apply the new wearable
                if(_wearableMeshRenderer != null)
                {
                    g.SetActive(true);
                    Transform[] newBones = new Transform[_wearableMeshRenderer.bones.Length];
                    for (int i = 0; i < _wearableMeshRenderer.bones.Length; ++i)
                    {
                        GameObject bone = _wearableMeshRenderer.bones[i].gameObject;
                        if (_personBoneMap.ContainsKey(bone.name))
                        {
                            newBones[i] = _personBoneMap[bone.name];
                        }
                        else
                        {
                            newBones[i] = _wearableMeshRenderer.bones[i];
                        }
                    }
                    _wearableMeshRenderer.bones = newBones;
                    wearable.AttachGameObject(g);
                }
                else
                {
                    if (editorMode)
                    {
                        DestroyImmediate(g);
                    }
                    else
                    {
                        Destroy(g);
                    }
                }                
            }
        }
    }

    /// <summary>
    /// Removes an existing wearable object.
    /// </summary>
    /// <param name="wearableType"></param>
    public void RemoveWearable(WearableType wearableType)
    {
        RemoveWearable(wearableType, false);
    }

    /// <summary>
    /// Removes an existing wearable object.
    /// </summary>
    /// <param name="wearableType"></param>
    /// <param name="editorMode">Is Editor Mode or not.</param>
    public void RemoveWearable(WearableType wearableType, bool editorMode)
    {
        if (_wearables != null)
        {
            //Iterate over all the associated wearables and remove the active one.
            foreach(Wearable w in _wearables)
            {
                if (w.WearableType == wearableType && w.IsActive)
                {
                    if(editorMode)
                    {
                        DestroyImmediate(w.AttachedGameObject);
                    }
                    else
                    {
                        Destroy(w.AttachedGameObject);
                    }
                    w.RemoveGameObject();
                }                
            }
        }
    }

    /// <summary>
    /// Removes all wearables.
    /// </summary>
    /// <param name="removeHair"></param>
    public void RemoveAllWearables(bool removeHair)
    {
        RemoveAllWearables(false, removeHair);
    }

    /// <summary>
    /// Removes all wearables.
    /// </summary>
    /// <param name="editorMode">Is Editor Mode or not.</param>
    /// <param name="removeHair">Remove hair or not.</param>
    public void RemoveAllWearables(bool editorMode, bool removeHair)
    {
        foreach(Wearable wearable in _wearables)
        {
            if(wearable.WearableType != WearableType.Hair
                || wearable.WearableType == WearableType.Hair && removeHair)
            {
                RemoveWearable(wearable.WearableType, editorMode);
            }            
        }
    }

    #endregion Wearable Methods
}
