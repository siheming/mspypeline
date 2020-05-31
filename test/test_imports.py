def test_import():
    # test classes
    from mspypeline import MSPParser, MSPInitializer, MSPPlots, MSPGUI, UIHandler, MaxQuantPlotter
    import inspect
    assert inspect.isclass(MSPPlots)
    assert inspect.isclass(MSPInitializer)
    assert inspect.isclass(MSPGUI)
    assert inspect.isclass(MSPParser)
    assert inspect.isclass(UIHandler)
    assert inspect.isclass(MaxQuantPlotter)
    # test modules
    from mspypeline import matplotlib_plots, plotly_plots
    assert inspect.ismodule(matplotlib_plots)
    assert inspect.ismodule(plotly_plots)
