# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
from dataclasses import dataclass, field
from typing import List, Optional


@dataclass
class Tree:
    _root_node: Optional[int] = None
    _children: List[List[int]] = field(default_factory=list)
    _parent: List[Optional[int]] = field(default_factory=list)
    _visited_children_idx: List[int] = field(default_factory=list)
    _curr_node_idx: Optional[int] = None

    num_nodes: int = 0

    def add_root(self) -> None:
        self._root_node = 0
        self.num_nodes += 1
        self._children.append([])
        self._parent.append(None)
        self._visited_children_idx.append(0)

    def add_child(self, idx: int) -> None:
        self._children[idx].append(self.num_nodes)
        self._parent.append(idx)
        self._children.append([])
        self.num_nodes += 1

    def get_root_idx(self) -> Optional[int]:
        return self._root_node

    def get_children_idx(self, idx: int) -> List[int]:
        """Returns the index of all the children"""
        return self._children[idx]

    def get_parent_idxs(self, idx: int) -> List[int]:
        """Returns all parents of a node"""
        parents = []
        curr_par = self._parent[idx]
        if curr_par is None:
            return parents
        else:
            parents.append(curr_par)
        while curr_par is not None:
            curr_par = self._parent[curr_par]
            if curr_par is not None:
                parents.append(curr_par)
        return parents

    def get_next_idx(self, idx: Optional[int] = None) -> Optional[int]:
        """Returns the index of the next node to visit in a depth first search"""
        if self._curr_node_idx is None:
            # just starting the walk, start at root
            self._curr_node_idx = self._root_node
            return self._curr_node_idx

        if idx is not None:
            self._curr_node_idx = idx

        # try to go deeper on current node
        curr_node_children = self.get_children_idx(self._curr_node_idx)
        nchildren = len(curr_node_children)
        if nchildren and self._visited_children_idx[self._curr_node_idx] != nchildren:
            new_node_idx = curr_node_children[self._visited_children_idx[self._curr_node_idx]]
            self._visited_children_idx[self._curr_node_idx] += 1
            self._curr_node_idx = new_node_idx
            return self._curr_node_idx
        elif self._parent[self._curr_node_idx] is not None:
            # else go up and check for children
            return self.get_next_idx(self._parent[self._curr_node_idx])
        else:
            return self._root_node
