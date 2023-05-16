using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Helper
{
    public SkinnedMeshRenderer GetSkinMeshRecursive(Transform t)
    {
        Transform tempTransform;
        SkinnedMeshRenderer smr = null;

        var queue = new Queue<Transform>();
        queue.Enqueue(t);

        while (queue.Count > 0)
        {
            // Take the next node from the front of the queue
            var node = queue.Dequeue();

            // Process the node 'node'
            smr = node.GetComponent<SkinnedMeshRenderer>();
            if (!(smr == null || smr.ToString() == "null"))
            {
                return smr;
            }

            // Add the node’s children to the back of the queue
            int childCount = node.childCount;
            for (int i = 0; i < childCount; i++)
            {
                tempTransform = node.GetChild(i);
                queue.Enqueue(tempTransform);
            }
        }

        // None of the nodes matched the specified predicate.
        return null;
    }

    public Transform FindTransform(Transform parent, string name)
    {
        if (parent.name.Equals(name)) return parent;
        foreach (Transform child in parent)
        {
            Transform result = FindTransform(child, name);
            if (result != null) return result;
        }
        return null;
    }
}
