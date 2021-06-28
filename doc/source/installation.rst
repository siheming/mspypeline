.. _installation:

Installing mspypeline
=====================

An installation of python is required to use this package. Python can be installed as part of the
`Anaconda <https://www.anaconda.com/products/individual>`__ distribution, which is recommended, or
via `python <https://www.python.org/downloads/>`__
(`python 3.7 <https://www.python.org/downloads/release/python-375/>`__). Then ``mspypeline``
needs to be installed using a terminal (e.g. the terminal provided by anaconda - anaconda prompt).

.. _conda-installation:

Anaconda installation (recommended)
***********************************
When installing ``mspypeline`` via Anaconda, one can choose two possible options to do so.

1. Either a new virtual `conda environment <https://conda.io/projects/conda/en/latest/user-guide/concepts/environments.html>`__
   is created and ``mspypeline`` is installed therein. This can be achieved e.g. by running the following code in a terminal:

.. code-block:: bash

    # create new environment called mspypeline and install mspypeline package
    # + all required packages
    conda create --name mspypeline python=3.7 mspypeline -c conda-forge -c siheming

    # activate the newly created environment
    conda activate mspypeline

    # in the activated environment the python package mspypeline can be called, the GUI
    # can be started or the package can be imported

With a new virtual environment in which python 3.7, ``mspypeline`` and all other required packages are
installed, it is necessary to activate this environment before an analysis from the terminal in order to use
``mspypeline``. Once the installation is performed and the environment is activated an analysis can be started
following :ref:`the next steps <get-started>`.

2. Otherwise, ``mspypeline`` can be installed into the base conda environment. This can be achieved e.g. by running the
   following code in a terminal:

.. code-block:: bash

    conda install -c conda-forge -c siheming mspypeline

With this base installation of ``mspypeline`` and all required packages it is possible to immediately start an analysis
following :ref:`the next steps <get-started>`.


.. _pip-installation:

pip(PyPI) installation
**********************
The ``mspypeline`` python package can optionally also be installed from PyPI and one can choose two possible
installation options.

1. Either a new virtual `python environment <https://docs.python.org/3.7/tutorial/venv.html>`__ is created and
   ``mspypeline`` is installed therein. Since the python version cannot be changed by pip make sure that the correct
   python version ist installed. This can be achieved e.g. by running the following code in a terminal (might differ
   based on your OS system):

.. code-block:: bash

    # create new environment called mspypeline
    python3 -m venv mypypeline

    # activate the newly created environment
    source activate mspypeline

    # install the mspypeline package + all required package within the environment
    pip install mspypeline

With a new virtual environment in which python, ``mspypeline`` and all other required packages are
installed, it is necessary to activate this environment before an analysis from the terminal in order to use
``mspypeline``. Once the installation is performed and the environment is activated an analysis can be started
following :ref:`the next steps <get-started>`.

2. Otherwise, ``mspypeline`` can be installed into the base environment. This can be achieved e.g. by running the
   following code in a terminal:

.. code-block:: bash

    pip install mspypeline

With this base installation of ``mspypeline`` and all required packages it is possible to immediately start an analysis
following :ref:`the next steps <get-started>`.

.. _activate-venv:

Obtaining Sources
*****************
Get the source code by cloning the github project:

.. code-block:: bash

    git clone https://github.com/siheming/mspypeline.git

Download Sources from `PyPI <https://pypi.org/project/mspypeline/>`__.

Download Sources from `conda <https://anaconda.org/siheming/mspypeline>`__.


Python support
**************
Currently supported python versions are 3.7 and 3.8.

Dependencies
************
- `numpy <https://numpy.org/>`__ >= 1.17.4
- `pandas <https://pandas.pydata.org/>`__ >= 1.0.0
- `scipy <https://www.scipy.org/>`__ >= 1.3.1
- `matplotlib <https://matplotlib.org/>`__ >= 3.1.1
- `scikit-learn <https://scikit-learn.org/stable/>`__ >= 0.22.1
- tzlocal >= 2.0.0
- ruamel_yaml >= 3.4.2
- matplotlib-venn >= 0.11.5
- adjusttext >= 0.7.3.1
- plotly >= 4.6.0

Optional Dependencies for R packages
************************************
some plots might require additional R packages to be installed. Because of that additional dependencies are required for
those plots.

- rpy2 >= 3.4.2
