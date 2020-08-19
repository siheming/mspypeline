from sphinx.application import Sphinx
from docutils import nodes


def setup(app: Sphinx):
    app.add_config_value("massinclude_dirs", [], "html")
    app.add_node(massinclude_node)
    print(app.config.massinclude_dirs)
    return {"version": "0.1"}


class massinclude_node(nodes.General, nodes.Element):
    pass