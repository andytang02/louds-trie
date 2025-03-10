# louds-trie

LOUDS-trie implementation example (C++)

### Merge Explanation

The primary logic of the merge occurs in the recursive function: 

```void TrieImpl::merge_nodes(TrieImpl& trie1, TrieImpl& trie2, TrieImpl* merged_trie, uint64_t p1, bool valid1, uint64_t p2, bool valid2, uint64_t depth)```

At a high level, the function recursively merges the subtrees rooted at parent index p1 for trie1 and parent index p2 for trie2 at layer depth-1, into the corect spot in merged_trie. We will refer to these root nodes as root1 and root2. If valid1 is false, we only merge trie2, and if valid2 is false, we only merge trie1.

The keys in trie1 and trie2 must be merged in lexographical order, which is where most of the implementation detail lies. In particular, if a subtree in trie1 represents a prefix string lexographically less than a subtree in trie2, we should explore the entire subtree in trie1 first and merge that subtree into merged_trie first.

We maintain the invariant that root1 and root2 each correspond to the same prefix string in the respective tries in the recursive call.

The inputs are explained as follows:
- trie1, trie2: the LOUDS tries to be merged
- merged_trie: the merged LOUDS trie we recursively build
- p1: the index in the labels at layer depth-1 in trie1 representing the root1
- p2: the index in the labels at layer depth-1 in trie2 representing the root2
- valid1: whether or not to merge trie1
- valid2: whether or not to merge trie2
- depth: the depth in the trie of subtree being merged

We break down the function into sections:

Initialization:

```
 // Base case: if depth exceeds either trie's levels, stop
 if (depth >= merged_trie->levels_.size()) return;

 // get the level in trie1 and trie2 that we are merging
 Level& level1 = (depth < trie1.levels_.size()) ? trie1.levels_[depth] : trie1.levels_.back();
 Level& level2 = (depth < trie2.levels_.size()) ? trie2.levels_[depth] : trie2.levels_.back();

 // gets the start and end indices of children of root1 in the current level
 // In the special case valid1 is false, sets both i1 and end1 to 0
 uint64_t i1 = (!valid1 || p1 == 0) ? 0 : level1.louds.select(p1 - 1) + 1 - p1;
 uint64_t end1 = (!valid1) ? 0 : level1.louds.select(p1) - p1;

 // gets the start and end indices of children of root2 in the current level
 // In the special case valid2 is false, sets both i1 and end1 to 0
 uint64_t i2 = (!valid2 || p2 == 0) ? 0 : level2.louds.select(p2 - 1) + 1 - p2;
 uint64_t end2 = (!valid2) ? 0 : level2.louds.select(p2) - p2;
```

Loop over children:
```
 // We iterate from i1 to end1 and i2 to end2, which represent the children of root1 and root2
 while (i1 < end1 || i2 < end2) {
```

Case 1, only root2 has children left to merge. Therefore, we explore the entire subtrees at the child nodes of root2:
```
if (i1 >= end1) { // Only root2 has children left to merge
  merged_trie->levels_[depth].louds.add(0); // add a child node in merged_trie 
  merged_trie->levels_[depth].outs.add(level2.outs.get(i2)); // represents whether this child node is a terminal node
  merged_trie->levels_[depth].labels.push_back(level2.labels[i2]); // represents the label of this child node
  if (level2.outs.get(i2)) { // if node is a terminal node
    merged_trie->n_keys_++; // update n_keys 
    merged_trie->levels_[depth + 1].offset++; // update offset
  }
  merged_trie->n_nodes_++; // increment n_nodes

  merge_nodes(trie1, trie2, merged_trie, 0, false, i2, true, depth + 1);  // recurse into subtree at child node, only considering trie2
  i2++; // go to next child node at this level
}
```

Case 2, only root1 has children left to merge. Therefore, we explore the entire subtrees at the child nodes of root1:
```
else if (i2 >= end2) {  // Only root1 has children left to merge
  merged_trie->levels_[depth].louds.add(0); // add a child node in merged_trie 
  merged_trie->levels_[depth].outs.add(level1.outs.get(i1)); // represents whether this child node is a terminal node
  merged_trie->levels_[depth].labels.push_back(level1.labels[i1]); // represents the label of this child node
  if (level1.outs.get(i1)) { // if node is a terminal node
    merged_trie->n_keys_++; // update n_keys 
    merged_trie->levels_[depth + 1].offset++; // update offset
  }
  merged_trie->n_nodes_++; // increment n_nodes

  merge_nodes(trie1, trie2, merged_trie, i1, true, 0, false, depth + 1);  // recurse into subtree at child node, only considering trie1
  i1++;
}
```

Case 3, the current child node in trie1 represents a prefix lexographically less than the current child node in trie2:
```
// note that we maintain that root1 and root2 represent the same prefix string, so this condition means that the current child node in
trie1 represents a prefix lexographically less than the current child node in trie2
else if (level1.labels[i1] < level2.labels[i2]) {
   // logic same as previous case
   merged_trie->levels_[depth].louds.add(0); 
   merged_trie->levels_[depth].outs.add(level1.outs.get(i1));
   merged_trie->levels_[depth].labels.push_back(level1.labels[i1]);
   if (level1.outs.get(i1)) {
     merged_trie->n_keys_++;
     merged_trie->levels_[depth + 1].offset++;
   }
   merged_trie->n_nodes_++;

   merge_nodes(trie1, trie2, merged_trie, i1, true, 0, false, depth + 1);   // recurse into subtree at child node, only considering trie1
   i1++;
 }
```

Case 4, the current child node in trie2 represents a prefix lexographically less than the current child node in trie1:

```
// This is just the reverse of Case 3
else if (level1.labels[i1] > level2.labels[i2]) {
     merged_trie->levels_[depth].louds.add(0); 
     merged_trie->levels_[depth].outs.add(level2.outs.get(i2));
     merged_trie->levels_[depth].labels.push_back(level2.labels[i2]);
     if (level2.outs.get(i2)) {
       merged_trie->n_keys_++;
       merged_trie->levels_[depth + 1].offset++;
     }
     merged_trie->n_nodes_++;

     merge_nodes(trie1, trie2, merged_trie, 0, false, i2, true, depth + 1);  // Recurse into trie2
     i2++;
}
```

Case 5, the current child node in trie1 and current child node in trie2 represent the same prefix string, so we recurse:
```
// level1.labels[i1] == level2.labels[i2], so prefixes match
 else {
      merged_trie->levels_[depth].louds.add(0); // add a child node in merged_trie 
      merged_trie->levels_[depth].outs.add(level1.outs.get(i1) | level2.outs.get(i2)); // is a terminal node if is either terminal node in trie1 or trie2
      merged_trie->levels_[depth].labels.push_back(level1.labels[i1]); // add label
      if (level1.outs.get(i1) | level2.outs.get(i2)) { // update meta data if is a terminal node
        merged_trie->n_keys_++;
        merged_trie->levels_[depth + 1].offset++;
      }
      merged_trie->n_nodes_++;

      merge_nodes(trie1, trie2, merged_trie, i1, true, i2, true, depth + 1);  // recurse into both subtrees simultaneously
      i1++; // go to next child node at this level
      i2++; // go to next child node at this level
    }
  }
}
```

Finally, we mark the end of children at this depth with a 1:
```
merged_trie->levels_[depth].louds.add(1);  // Mark end of children
```

The full implementation is below:


```
void TrieImpl::merge_nodes(TrieImpl& trie1, TrieImpl& trie2, TrieImpl* merged_trie, uint64_t p1, bool valid1, uint64_t p2, bool valid2, uint64_t depth) {
  // Base case: if depth exceeds either trie's levels, stop
  if (depth >= merged_trie->levels_.size()) return;

  Level& level1 = (depth < trie1.levels_.size()) ? trie1.levels_[depth] : trie1.levels_.back();
  Level& level2 = (depth < trie2.levels_.size()) ? trie2.levels_[depth] : trie2.levels_.back();

  uint64_t i1 = (!valid1 || p1 == 0) ? 0 : level1.louds.select(p1 - 1) + 1 - p1;
  uint64_t end1 = (!valid1) ? 0 : level1.louds.select(p1) - p1;

  uint64_t i2 = (!valid2 || p2 == 0) ? 0 : level2.louds.select(p2 - 1) + 1 - p2;
  uint64_t end2 = (!valid2) ? 0 : level2.louds.select(p2) - p2;

  // Merge children
  while (i1 < end1 || i2 < end2) {
    if (i1 >= end1) {  // Only trie2 has children left
      merged_trie->levels_[depth].louds.add(0);
      merged_trie->levels_[depth].outs.add(level2.outs.get(i2));
      merged_trie->levels_[depth].labels.push_back(level2.labels[i2]);
      if (level2.outs.get(i2)) {
        merged_trie->n_keys_++;
        merged_trie->levels_[depth + 1].offset++;
      }
      merged_trie->n_nodes_++;

      merge_nodes(trie1, trie2, merged_trie, 0, false, i2, true, depth + 1);  // Recurse into trie2
      i2++;
    } else if (i2 >= end2) {  // Only trie1 has children left
      merged_trie->levels_[depth].louds.add(0);
      merged_trie->levels_[depth].outs.add(level1.outs.get(i1));
      merged_trie->levels_[depth].labels.push_back(level1.labels[i1]);
      if (level1.outs.get(i1)) {
        merged_trie->n_keys_++;
        merged_trie->levels_[depth + 1].offset++;
      }
      merged_trie->n_nodes_++;

      merge_nodes(trie1, trie2, merged_trie, i1, true, 0, false, depth + 1);  // Recurse into trie1
      i1++;
    } else if (level1.labels[i1] < level2.labels[i2]) {
      merged_trie->levels_[depth].louds.add(0);
      merged_trie->levels_[depth].outs.add(level1.outs.get(i1));
      merged_trie->levels_[depth].labels.push_back(level1.labels[i1]);
      if (level1.outs.get(i1)) {
        merged_trie->n_keys_++;
        merged_trie->levels_[depth + 1].offset++;
      }
      merged_trie->n_nodes_++;

      merge_nodes(trie1, trie2, merged_trie, i1, true, 0, false, depth + 1);  // Recurse into trie1
      i1++;
    } else if (level1.labels[i1] > level2.labels[i2]) {
      merged_trie->levels_[depth].louds.add(0);
      merged_trie->levels_[depth].outs.add(level2.outs.get(i2));
      merged_trie->levels_[depth].labels.push_back(level2.labels[i2]);
      if (level2.outs.get(i2)) {
        merged_trie->n_keys_++;
        merged_trie->levels_[depth + 1].offset++;
      }
      merged_trie->n_nodes_++;

      merge_nodes(trie1, trie2, merged_trie, 0, false, i2, true, depth + 1);  // Recurse into trie2
      i2++;
    } else {  // Labels match
      merged_trie->levels_[depth].louds.add(0);
      merged_trie->levels_[depth].outs.add(level1.outs.get(i1) | level2.outs.get(i2));
      merged_trie->levels_[depth].labels.push_back(level1.labels[i1]);
      if (level1.outs.get(i1) | level2.outs.get(i2)) {
        merged_trie->n_keys_++;
        merged_trie->levels_[depth + 1].offset++;
      }
      merged_trie->n_nodes_++;

      merge_nodes(trie1, trie2, merged_trie, i1, true, i2, true, depth + 1);  // Recurse into both
      i1++;
      i2++;
    }
  }
  merged_trie->levels_[depth].louds.add(1);  // Mark end of children
}
```

The top level function just starts at the root node:
```
TrieImpl* TrieImpl::merge_efficient(TrieImpl& trie1, TrieImpl& trie2) {
  TrieImpl* merged_trie = new TrieImpl(false);
  merged_trie->levels_.resize(max(trie1.levels_.size(), trie2.levels_.size()));

  TrieImpl::merge_nodes(trie1, trie2, merged_trie, 0, true, 0, true, 1);

  return merged_trie;
}
```

Trie calls the TrieImpl version:

```
Trie* Trie::merge_efficient(const Trie& trie1, const Trie& trie2) {
  Trie* merged_trie = new Trie();

  merged_trie->impl_ = TrieImpl::merge_efficient(*trie1.impl_,*trie2.impl_);

  return merged_trie;
}
```
