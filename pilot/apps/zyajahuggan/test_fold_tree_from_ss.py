import pytest
from bootcamp_app import identify_secondary_structure_spans, fold_tree_from_dssp_string, fold_tree_from_ss
from pyrosetta.rosetta.core.kinematics import FoldTree

ss1 =  "   EEEEE   HHHHHHHH  EEEEE   IGNOR EEEEEE   HHHHHHHHHHH  EEEEE  HHHH   "
expected1 = [(4, 8), (12, 19), (22, 26), (36, 41), (45, 55), (58, 62), (65, 68)]

ss2 = "HHHHHHH   HHHHHHHHHHHH      HHHHHHHHHHHHEEEEEEEEEEHHHHHHH EEEEHHH "
expected2 = [(1, 7), (11, 22), (29, 40), (41, 50), (51, 57), (59, 62), (63, 65)]

ss3 = "EEEEEEEEE EEEEEEEE EEEEEEEEE H EEEEE H H H EEEEEEEE"
expected3 = [(1,9), (11, 18), (20, 28), (30, 30), (32, 36), (38, 38), (40, 40), (42, 42), (44, 51)]

@pytest.mark.parametrize("test_input,expected", [(ss1, expected1), (ss2, expected2), (ss3, expected3)])
def test_ss(test_input,expected):
    assert identify_secondary_structure_spans(test_input) == expected


ss = "   EEEEEEE    EEEEEEE         EEEEEEEEE    EEEEEEEEEE   HHHHHH         EEEEEEEEE         EEEEE     "
peptide_edges = [
    (7, 1, -1),
    (7, 10, -1),
    (7, 12, -1),
    (12, 11, -1),
    (12, 14, -1),

    (7, 18, -1),
    (18, 15, -1),
    (18, 21, -1),

    (7, 26, -1),
    (26, 22, -1),
    (26, 30, -1),

    (7, 35, -1),
    (35, 31, -1),
    (35, 39, -1),

    (7, 41, -1),
    (41, 40, -1),
    (41, 43, -1),

    (7, 48, -1),
    (48, 44, -1),
    (48, 53, -1),

    (7, 55, -1),
    (55, 54, -1),
    (55, 56, -1),

    (7, 59, -1),
    (59, 57, -1),
    (59, 62, -1),

    (7, 67, -1),
    (67, 63, -1),
    (67, 71, -1),

    (7, 76, -1),
    (76, 72, -1),
    (76, 80, -1),

    (7, 85, -1),
    (85, 81, -1),
    (85, 89, -1),

    (7, 92, -1),
    (92, 90, -1),
    (92, 99, -1),
]

jump_edges = [
    (7, 12, 1),
    (7, 18, 2),
    (7, 26, 3),
    (7, 35, 4),
    (7, 41, 5),
    (7, 48, 6),
    (7, 55, 7),
    (7, 59, 8),
    (7, 67, 9),
    (7, 76, 10),
    (7, 85, 11),
    (7, 92, 12),
]




@pytest.mark.parametrize("string", [ss])

def test_fold_tree_from_dssp_string(ss_string):
    expected = FoldTree(len(ss))
    expected.clear()

    # peptide edges
    for a, b, label in peptide_edges:
        expected.add_edge(a, b, label)

    # jump edges
    for a, b, label in jump_edges:
        expected.add_edge(a, b, label)

    actual = fold_tree_from_dssp_string(ss_string)

    exp_edges = {
        (expected.edge(i).start(), expected.edge(i).stop(), expected.edge(i).label())
        for i in range(expected.nedges())
    }
    act_edges = {
        (actual.edge(i).start(), actual.edge(i).stop(), actual.edge(i).label())
        for i in range(actual.nedges())
    }

    assert act_edges == exp_edges

