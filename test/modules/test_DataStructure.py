import pytest


def test_simple_structure():
    from mspypeline import DataNode, DataTree
    import numpy as np
    import pandas as pd
    analysis_design = {
        "Ex1": {
            "A": {
                "1": "Ex1_A_1",
                "2": "Ex1_A_2",
            },
            "B": {
                "1": "Ex1_B_1",
                "2": "Ex1_B_2",
            },
        },
        "Ex2": {
            "A": {
                "1": "Ex2_A_1",
                "2": "Ex2_A_2",
            },
            "B": {
                "1": "Ex2_B_1",
                "2": "Ex2_B_2",
            },
        }
    }
    data = pd.DataFrame(np.random.random((10, 8)),
                        columns=["Ex1_A_1", "Ex1_A_2", "Ex1_B_1", "Ex1_B_2",
                                 "Ex2_A_1", "Ex2_A_2", "Ex2_B_1", "Ex2_B_2"])
    should_aggregate_technical_replicates = True
    tree: DataTree = DataTree.from_analysis_design(analysis_design, data, should_aggregate_technical_replicates)
    assert isinstance(tree["Ex1"], DataNode)
    assert isinstance(tree["Ex1"]["A"], DataNode)
    assert isinstance(tree["Ex1_A"], DataNode)
    assert isinstance(tree.root["Ex1"], DataNode)
    assert isinstance(str(tree["Ex1"]), str)
    assert isinstance(repr(tree["Ex1"]), str)
    assert tree.groupby().shape == (10, 2)
    assert "groupcol" in tree.groupby(0, new_col_name="groupcol").columns.names
    assert tree.groupby(0).shape == (10, 2)
    assert tree.groupby(0, method=None).shape == (10, 4)
    assert tree.groupby(1).shape == (10, 4)
    assert tree.groupby("Ex1", method=None).equals(tree["Ex1"].groupby(method=None))
    with pytest.raises(ValueError):
        tree.groupby((1, 1))
    assert tree.aggregate().shape == (10,)
    assert tree.aggregate(index=1).shape == (1,)
    assert tree.aggregate(index=[1, 2]).shape == (2,)
    assert tree.aggregate(method=None).shape == (10, 4)
    assert tree.aggregate("Ex1", method=None).shape == (10, 2)
    with pytest.raises(ValueError):
        tree.aggregate(123)
    assert tree["Ex1"].aggregate(method=None).shape == (10, 2)
    assert len(tree.level_keys_full_name[0]) == 2
    assert len(tree.level_keys_full_name[1]) == 4
    assert tree.root.get_total_number_children() == 4
    assert tree["Ex1"].get_total_number_children() == 2
    assert tree["Ex1_A"].get_total_number_children() == 0
    # #######################
    assert isinstance(tree["Ex1_A_1"], DataNode)
    assert tree.groupby(0, go_max_depth=True).shape == (10, 2)
    assert tree.aggregate(method=None, go_max_depth=True).shape == (10, 8)
    assert tree["Ex1"].aggregate(method=None, go_max_depth=True).shape == (10, 4)
    assert tree.root.get_total_number_children(go_max_depth=True) == 8
    assert tree["Ex1"].get_total_number_children(go_max_depth=True) == 4
    assert tree["Ex1_A"].get_total_number_children(go_max_depth=True) == 2
    assert tree["Ex1_A_1"].get_total_number_children(go_max_depth=True) == 0
