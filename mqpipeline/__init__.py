from .version import __version__
import os
# import for "from package import *"
__all__ = [
    "create_app"
]

path_package = script_loc = os.path.dirname(os.path.realpath(__file__))
path_package_config = os.path.join(path_package, "config")


def create_app(test_config=None):
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

    from . import login_db
    login_db.init_app(app)

    from . import analysis_page
    app.register_blueprint(analysis_page.bp)

    return app


if __name__ == "__main__":
    from mqpipeline.core import MQParser, MQUI
    mqparser = MQParser()
    if True:
    #if mqparser.args_dict.get("host_flask"):
        app = create_app()
        app.run()
    gui = MQUI(**mqparser.args_dict)
