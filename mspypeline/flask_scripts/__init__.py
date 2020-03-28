
__all__ = [
    "create_app_helper",
]


def create_app_helper(test_config):
    from flask import Flask
    import os
    # create and configure the app
    app = Flask(__name__, instance_relative_config=True)
    app.config.from_mapping(
        SECRET_KEY='dev',
        DATABASE=os.path.join(app.instance_path, 'flaskr.sqlite'),
    )

    if test_config is None:
        # load the instance config, if it exists, when not testing
        app.config.from_pyfile('config.py', silent=True)
    else:
        # load the test config if passed in
        app.config.from_mapping(test_config)

    # ensure the instance folder exists
    os.makedirs(app.instance_path, exist_ok=True)

    # a simple page that says hello
    @app.route('/hello')
    def hello():
        return 'Hello, World!'

    from mspypeline.flask_scripts import login_db
    login_db.init_app(app)

    from mspypeline.flask_scripts.blueprints import analysis_page
    app.register_blueprint(analysis_page.bp)

    from mspypeline.flask_scripts.blueprints import register_analysis
    app.register_blueprint(register_analysis.bp)

    app.add_url_rule('/', view_func=register_analysis.create_analysis)

    return app
