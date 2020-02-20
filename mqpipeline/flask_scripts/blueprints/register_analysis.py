from flask import Blueprint, request, redirect, url_for, flash, render_template

bp = Blueprint('register_analysis', __name__)


@bp.route("/", methods=("GET", "POST"))
def index():
    return redirect(url_for("register_analysis.create_analysis"))


@bp.route('/start', methods=('GET', 'POST'))
def create_analysis():
    if request.method == "POST":
        path = request.form["path"]
        error = None

        if not path:
            error = "path is required"

        if error is None:

            return redirect(url_for("analysis_page.index", user_id=path))

        flash(error)

    return render_template('register_analysis.html')
