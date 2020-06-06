def test_import():
    import inspect
    # test imports from core
    from mspypeline import MSPParser, MSPInitializer, BasePlotter, MSPGUI, UIHandler, MaxQuantPlotter
    assert inspect.isclass(BasePlotter)
    assert inspect.isclass(MSPInitializer)
    assert inspect.isclass(MSPGUI)
    assert inspect.isclass(MSPParser)
    assert inspect.isclass(UIHandler)
    assert inspect.isclass(MaxQuantPlotter)
    # test imports from file reader
    from mspypeline import BaseReader, MQReader
    assert inspect.isclass(BaseReader)
    assert inspect.isclass(MQReader)
    # test imports from modules
    from mspypeline import DataTree, DataNode, MedianNormalizer, QuantileNormalizer, TailRobustNormalizer
    assert inspect.isclass(DataTree)
    assert inspect.isclass(DataNode)
    assert inspect.isclass(MedianNormalizer)
    assert inspect.isclass(QuantileNormalizer)
    assert inspect.isclass(TailRobustNormalizer)
    from mspypeline import interpolate_data, default_normalizers
    assert inspect.isfunction(interpolate_data)
    assert isinstance(default_normalizers, dict)
    # test imports from helpers
    from mspypeline.helpers import get_intersection_and_unique, get_number_of_non_na_values, default_to_regular,\
        get_number_rows_cols_for_fig, venn_names, plot_annotate_line, get_logger, dict_depth, get_legend_elements,\
        get_plot_name_suffix, get_analysis_design, fill_dict
    assert inspect.isfunction(get_intersection_and_unique)
    assert inspect.isfunction(get_number_of_non_na_values)
    assert inspect.isfunction(default_to_regular)
    assert inspect.isfunction(get_number_rows_cols_for_fig)
    assert inspect.isfunction(venn_names)
    assert inspect.isfunction(plot_annotate_line)
    assert inspect.isfunction(get_logger)
    assert inspect.isfunction(dict_depth)
    assert inspect.isfunction(get_legend_elements)
    assert inspect.isfunction(get_plot_name_suffix)
    assert inspect.isfunction(get_analysis_design)
    assert inspect.isfunction(fill_dict)
    # test plotting backends
    from mspypeline import matplotlib_plots, plotly_plots
    assert inspect.ismodule(matplotlib_plots)
    assert inspect.ismodule(plotly_plots)


def test_config_available():
    from mspypeline import path_package_config
    import os
    assert os.path.isdir(path_package_config)

