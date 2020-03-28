import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go

from mspypeline.helpers import get_number_rows_cols_for_fig


def create_plot():
    import plotly
    import plotly.graph_objs as go

    import pandas as pd
    import numpy as np
    import json

    N = 40
    x = np.linspace(0, 1, N)
    y = np.random.randn(N)
    df = pd.DataFrame({'x': x, 'y': y}) # creating a sample dataframe

    fig = go.Bar(
            x=df['x'], # assign x as the dataframe column 'x'
            y=df['y']
    )
    return go.Figure(fig)


def plot_intensity_histogram(df):
    n_rows, n_cols = get_number_rows_cols_for_fig(df.shape[1])
    fig = make_subplots(rows=n_rows, cols=n_cols)
    for i, col_name in enumerate(df.columns):
        row = (i // n_rows) + 1
        col = (i % n_rows) + 1
        fig.add_trace(go.Histogram(x=df[col_name]), row=row, col=col)
    #

    return fig
