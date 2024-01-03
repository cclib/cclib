# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.
import abc
from abc import ABC


class node:
    _idx = None  # make this a property so it can be getted
    _children = []
    _parent = None
    _last_visited_child_idx = None

    def __init__(self, idx, parent=None):
        _idx = idx
        _parent = None

    def add_child(self, idx):
        if _last_visited_child_idx is None:
            _last_visited_child_idx = 0
        _children.append(node(idx, parent=self))

    def visit_child(self):
        if _last_visited_child_idx is None:
            return None

        if len(_children) == 0:
            return None

        _last_visited_child_idx += 1
        return _children[_last_visited_child_idx - 1]

    def visit_parent(self):
        return _parent


# infinite vertical and horizontal type nodes


# in case we have other types of graphs
class basegraph(ABC):
    @abc.abstractmethod
    def not_sure():
        pass


class graph(basegraph):
    _root_node = None
    num_nodes = 0
    curr_added_idx = 0
    curr_node = None

    def __init__(self):
        add_root(node(num_nodes))
        curr_node = _root_node

    # list init?
    def add_root(self):
        _root_node = node(idx=0)

    def add_child(self):
        # curr_node.add_child
        # increase added_idx?
        pass

    def get_root_idx(self):
        return _root_node.get_idx()

    def get_next_idx(self):
        next_node = curr_node.visit_child()
        if next_node is not None:
            curr_node = next_node
            return curr_node.get_idx()

        while next_node is None or next_node is not _root_node:
            curr_node = curr_node.visit_parent()
            next_node = curr_node.visit_child()

        if next_node is _root_node:
            # we fully explored the tree and there are no unvisited nodes? throw error?
            return -1
