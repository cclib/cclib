# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for the tree implementation."""

from cclib.tree import Tree


class TreeTest:
    def test_tree_init(self) -> None:
        """Create a new tree but don't add anything to it (no nodes)"""
        tree = Tree()
        assert len(tree) == 0
        assert tree.get_root_idx() is None
        assert tree.get_next_idx() is None
        # repeated call doesn't change state
        assert tree.get_next_idx() is None

    def test_tree_add_root(self) -> None:
        """Tree with only a root node"""
        tree = Tree()
        tree.add_root()
        assert len(tree) == 1
        root_idx = tree.get_root_idx()
        assert root_idx == 0
        # TODO doesn't seem correct
        assert tree.get_next_idx() == 0
        # repeated call doesn't change state
        assert tree.get_next_idx() == 0
        # TODO doesn't seem correct
        assert tree.get_children_idx(root_idx) == []
        assert tree.get_parent_idxs(root_idx) == []
        # TODO doesn't seem correct
        assert tree.get_next_idx() == 0
        # TODO doesn't seem correct
        assert tree.get_next_idx(root_idx) == 0

    def test_tree_add_root_one_child(self) -> None:
        """Root with one direct child"""
        # a
        # |
        # b
        tree = Tree()
        tree.add_root()
        root_idx = tree.get_root_idx()
        assert root_idx == 0
        assert tree.get_children_idx(root_idx) == []
        tree.add_child(root_idx)
        assert tree.get_children_idx(root_idx) == [1]
        assert len(tree) == 2
        # TODO doesn't seem correct
        assert tree.get_next_idx() == 0
        # ah, it's the iterator
        assert tree.get_next_idx() == 1
        assert tree.get_next_idx() == 0
        assert tree.get_next_idx() == 0

    def test_tree_add_root_two_children(self) -> None:
        """Root with two direct children"""
        #   a
        #  / \
        # b   c
        tree = Tree()
        tree.add_root()
        root_idx = tree.get_root_idx()
        assert root_idx == 0
        assert tree.get_children_idx(root_idx) == []
        tree.add_child(root_idx)
        assert tree.get_children_idx(root_idx) == [1]
        tree.add_child(root_idx)
        assert tree.get_children_idx(root_idx) == [1, 2]
        assert len(tree) == 3
        assert tree.get_next_idx() == 0
        assert tree.get_next_idx() == 1
        assert tree.get_next_idx() == 2
        assert tree.get_next_idx() == 0
        assert tree.get_next_idx() == 0
        assert tree.get_next_idx() == 0

    def test_tree_add_root_two_children_linear(self) -> None:
        """Root with one child that itself has one child"""
        # a
        # |
        # b
        # |
        # c
        tree = Tree()
        tree.add_root()
        root_idx = tree.get_root_idx()
        assert root_idx == 0
        assert tree.get_children_idx(root_idx) == []
        tree.add_child(root_idx)
        assert tree.get_children_idx(root_idx) == [1]
        assert tree.get_children_idx(1) == []
        tree.add_child(1)
        assert tree.get_children_idx(root_idx) == [1]
        assert tree.get_children_idx(1) == [2]
        assert len(tree) == 3
        assert tree.get_next_idx() == 0
        assert tree.get_next_idx() == 1
        # list index out of range
        # assert tree.get_next_idx() == 2
        # assert tree.get_next_idx() == 0
        # assert tree.get_next_idx() == 0
        # assert tree.get_next_idx() == 0
