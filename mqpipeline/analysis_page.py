import functools

from flask import (
    Blueprint, flash, g, redirect, render_template, request, session, url_for
)
from werkzeug.security import check_password_hash, generate_password_hash

from mqpipeline.login_db import get_db
from mqpipeline.core.MQPlots import create_plot

bp = Blueprint('analysis_page', __name__)


@bp.route("/<user_id>", methods=("GET", ))
def index(user_id=0):
    json = create_plot()
    return render_template("analysis_page/index.html", user_id=user_id, json=json)
