from _plotly_utils.utils import PlotlyJSONEncoder
from flask import (
    Blueprint, render_template
)
import json

from mspypeline.plotting_backend import plotly_plots

bp = Blueprint('analysis_page', __name__)


@bp.route("/analysis/<user_id>", methods=("GET", "POST"))
def index(user_id=0):
    # TODO instead get information from database for user
    """mqinit = MQInitializer(user_id)
    mqinit.init_config()
    # mqinit.configs.update(configs)
    mqinit.read_data()
    # create plotting_backend from initializer
    mqplots = MQPlots.from_MQInitializer(mqinit)
    df = mqplots.histogram_data()
    fig = plotly_plots.plot_intensity_histogram(df)"""

    fig = plotly_plots.create_plot()
    json_data = json.dumps(fig, cls=PlotlyJSONEncoder)
    return render_template("analysis_page.html", user_id=user_id, json=json_data)
