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
        # 0
        # |
        # 1
        tree = Tree()
        tree.add_root()
        root_idx = tree.get_root_idx()
        assert root_idx == 0
        assert tree.get_children_idx(root_idx) == []
        tree.add_child(root_idx)
        assert tree.get_children_idx(root_idx) == [1]
        assert tree.get_children_idx(1) == []
        assert tree.get_parent_idxs(1) == [0]
        assert len(tree) == 2
        # TODO doesn't seem correct
        assert tree.get_next_idx() == 0
        # ah, it's the iterator
        assert tree.get_next_idx() == 1
        assert tree.get_next_idx() == 0
        assert tree.get_next_idx() == 0

    def test_tree_add_root_two_children(self) -> None:
        """Root with two direct children"""
        #   0
        #  / \
        # 1   2
        tree = Tree()
        tree.add_root()
        root_idx = tree.get_root_idx()
        assert root_idx == 0
        assert tree.get_children_idx(root_idx) == []
        tree.add_child(root_idx)
        assert tree.get_children_idx(root_idx) == [1]
        assert tree.get_children_idx(1) == []
        tree.add_child(root_idx)
        assert tree.get_children_idx(root_idx) == [1, 2]
        assert tree.get_children_idx(1) == []
        assert tree.get_children_idx(2) == []
        assert tree.get_parent_idxs(root_idx) == []
        assert tree.get_parent_idxs(1) == [0]
        assert tree.get_parent_idxs(2) == [0]
        assert len(tree) == 3
        assert tree.get_next_idx() == 0
        assert tree.get_next_idx() == 1
        assert tree.get_next_idx() == 2
        assert tree.get_next_idx() == 0
        assert tree.get_next_idx() == 0
        assert tree.get_next_idx() == 0

    def test_tree_add_root_two_children_linear(self) -> None:
        """Root with one child that itself has one child"""
        # 0
        # |
        # 1
        # |
        # 2
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
        assert tree.get_children_idx(2) == []
        assert tree.get_parent_idxs(1) == [0]
        assert tree.get_parent_idxs(2) == [1, 0]
        assert len(tree) == 3
        assert tree.get_next_idx() == 0
        assert tree.get_next_idx() == 1
        # list index out of range
        # assert tree.get_next_idx() == 2
        # assert tree.get_next_idx() == 0
        # assert tree.get_next_idx() == 0
        # assert tree.get_next_idx() == 0

    def test_tree_three_layers_left(self) -> None:
        #     0
        #    / \
        #   1   2
        #  /
        # 3
        tree = Tree()
        tree.add_root()
        root_idx = tree.get_root_idx()
        assert root_idx == 0
        tree.add_child(root_idx)
        tree.add_child(root_idx)
        # TODO add_child should return the index it created
        tree.add_child(1)
        assert len(tree) == 4
        assert tree.get_children_idx(root_idx) == [1, 2]
        assert tree.get_children_idx(1) == [3]
        assert tree.get_children_idx(2) == []
        assert tree.get_children_idx(3) == []
        assert tree.get_parent_idxs(1) == [0]
        assert tree.get_parent_idxs(2) == [0]
        assert tree.get_parent_idxs(3) == [1, 0]
        assert tree.get_next_idx() == 0
        assert tree.get_next_idx() == 1
        # list index out of range
        # assert tree.get_next_idx() == 2
        # assert tree.get_next_idx() == 3
        # assert tree.get_next_idx() == 0
        # assert tree.get_next_idx() == 0
        # assert tree.get_next_idx() == 0
        # assert tree.get_next_idx() == 0

    def test_tree_three_layers_right(self) -> None:
        #     0
        #    / \
        #   1   2
        #        \
        #         3
        tree = Tree()
        tree.add_root()
        root_idx = tree.get_root_idx()
        assert root_idx == 0
        tree.add_child(root_idx)
        tree.add_child(root_idx)
        tree.add_child(2)
        assert len(tree) == 4
        assert tree.get_children_idx(root_idx) == [1, 2]
        assert tree.get_children_idx(1) == []
        assert tree.get_children_idx(2) == [3]
        assert tree.get_children_idx(3) == []
        assert tree.get_parent_idxs(1) == [0]
        assert tree.get_parent_idxs(2) == [0]
        assert tree.get_parent_idxs(3) == [2, 0]
        assert tree.get_next_idx() == 0
        assert tree.get_next_idx() == 1
        assert tree.get_next_idx() == 2
        # list index out of range
        # assert tree.get_next_idx() == 3
        # assert tree.get_next_idx() == 0
        # assert tree.get_next_idx() == 0
        # assert tree.get_next_idx() == 0
        # assert tree.get_next_idx() == 0

    def test_tree_three_layers_both(self) -> None:
        #     0
        #    / \
        #   1   2
        #  /     \
        # 3       4
        tree = Tree()
        tree.add_root()
        root_idx = tree.get_root_idx()
        assert root_idx == 0
        tree.add_child(root_idx)
        tree.add_child(root_idx)
        tree.add_child(1)
        tree.add_child(2)
        assert len(tree) == 5
        assert tree.get_children_idx(root_idx) == [1, 2]
        assert tree.get_children_idx(1) == [3]
        assert tree.get_children_idx(2) == [4]
        assert tree.get_children_idx(3) == []
        assert tree.get_children_idx(4) == []
        assert tree.get_parent_idxs(1) == [0]
        assert tree.get_parent_idxs(2) == [0]
        assert tree.get_parent_idxs(3) == [1, 0]
        assert tree.get_parent_idxs(4) == [2, 0]
        assert tree.get_next_idx() == 0
        assert tree.get_next_idx() == 1
        # list index out of range
        # assert tree.get_next_idx() == 2

    def test_tree_three_layers_both_swap(self) -> None:
        """Same as test_tree_three_layers_both but change the order in which leaves are added."""
        #     0
        #    / \
        #   1   2
        #  /     \
        # 4       3
        tree = Tree()
        tree.add_root()
        root_idx = tree.get_root_idx()
        assert root_idx == 0
        tree.add_child(root_idx)
        tree.add_child(root_idx)
        tree.add_child(2)
        tree.add_child(1)
        assert len(tree) == 5
        assert tree.get_children_idx(root_idx) == [1, 2]
        assert tree.get_children_idx(1) == [4]
        assert tree.get_children_idx(2) == [3]
        assert tree.get_children_idx(3) == []
        assert tree.get_children_idx(4) == []
        assert tree.get_parent_idxs(1) == [0]
        assert tree.get_parent_idxs(2) == [0]
        assert tree.get_parent_idxs(3) == [2, 0]
        assert tree.get_parent_idxs(4) == [1, 0]
        assert tree.get_next_idx() == 0
        assert tree.get_next_idx() == 1
        # list index out of range
        # assert tree.get_next_idx() == 2

    def test_tree_three_layers_double_left(self) -> None:
        #     0
        #    / \
        #   1   2
        #  / \
        # 3   4
        tree = Tree()
        tree.add_root()
        root_idx = tree.get_root_idx()
        assert root_idx == 0
        tree.add_child(root_idx)
        tree.add_child(root_idx)
        tree.add_child(1)
        tree.add_child(1)
        assert len(tree) == 5
        assert tree.get_children_idx(root_idx) == [1, 2]
        assert tree.get_children_idx(1) == [3, 4]
        assert tree.get_children_idx(2) == []
        assert tree.get_children_idx(3) == []
        assert tree.get_children_idx(4) == []
        assert tree.get_parent_idxs(1) == [0]
        assert tree.get_parent_idxs(2) == [0]
        assert tree.get_parent_idxs(3) == [1, 0]
        assert tree.get_parent_idxs(4) == [1, 0]
        assert tree.get_next_idx() == 0
        assert tree.get_next_idx() == 1
        # list index out of range
        # assert tree.get_next_idx() == 2
