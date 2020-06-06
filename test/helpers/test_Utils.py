def test_get_number_rows_cols_for_fig():
    from mspypeline.helpers import get_number_rows_cols_for_fig
    assert get_number_rows_cols_for_fig([1, 1, 1, 1]) == (2, 2)
    assert get_number_rows_cols_for_fig(4) == (2, 2)


def test_fill_dict():
    from mspypeline.helpers import fill_dict


def test_default_to_regular():
    from mspypeline.helpers import default_to_regular
    from collections import defaultdict
    d = defaultdict(int)
    d["a"] += 1
    assert isinstance(d, defaultdict)
    d = default_to_regular(d)
    assert isinstance(d, dict)
    assert not isinstance(d, defaultdict)


def test_get_analysis_design():
    from mspypeline.helpers import get_analysis_design
    assert get_analysis_design(["A1_1", "A1_2", "A2_1", "A2_2"]) == {
        'A1': {'1': 'A1_1', '2': 'A1_2'},
        'A2': {'1': 'A2_1', '2': 'A2_2'}
    }
    assert get_analysis_design(["A_1_1"]) == {"A": {"1": {"1": "A_1_1"}}}


def test_plot_annotate_line():
    from mspypeline.helpers import plot_annotate_line


def test_venn_names():
    from mspypeline.helpers import venn_names


def test_install_r_dependencies():
    from mspypeline.helpers.Utils import install_r_dependencies


def test_get_number_of_non_na_values():
    from mspypeline.helpers import get_number_of_non_na_values as gna
    assert gna(20) > gna(10) > gna(5) > gna(3)
    assert gna(3) == gna(2) and gna(3) == gna(1)


def test_get_intersection_and_unique():
    from mspypeline.helpers import get_intersection_and_unique
    import pandas as pd
    df1 = pd.DataFrame()
    df2 = pd.DataFrame()
    assert all(map(pd.Series.equals,
                   get_intersection_and_unique(df1, df2),
                   (pd.Series([], dtype=bool), pd.Series([], dtype=bool), pd.Series([], dtype=bool))))
    df1 = pd.DataFrame([[1, 1, 1], [1, 1, 1], [0, 0, 0], [1, 0, 0]])
    df2 = pd.DataFrame([[1, 1, 1], [0, 0, 0], [1, 1, 1], [1, 0, 0]])
    assert all(map(
        pd.Series.equals,
        get_intersection_and_unique(df1, df2),
        (pd.Series([1, 0, 0, 0], dtype=bool), pd.Series([0, 1, 0, 0], dtype=bool), pd.Series([0, 0, 1, 0], dtype=bool))
       ))


def test_dict_depth():
    from mspypeline.helpers import dict_depth
    assert dict_depth({}) == 0
    assert dict_depth({"test": 1}) == 1
    assert dict_depth({"test": {"test": {"test": 1}}}) == 3


def test_get_legend_elements():
    from mspypeline.helpers import get_legend_elements


def test_get_plot_name_suffix():
    from mspypeline.helpers import get_plot_name_suffix
    assert get_plot_name_suffix() == ""
    assert get_plot_name_suffix("test") == "_test"
    assert get_plot_name_suffix(level=1) == "_level_1"
    assert get_plot_name_suffix("test", 1) == "_test_level_1"
