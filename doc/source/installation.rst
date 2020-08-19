.. _installation:

Installing mspypeline
=====================

An installation of python is required to use this package. Python can be installed either
via `python <https://www.python.org/downloads/>`__
(`python 3.7 <https://www.python.org/downloads/release/python-375/>`__), or as part of the
`Anaconda <https://www.anaconda.com/products/individual>`__ distribution, which is recommended. Then ``mspypeline``
needs to be installed using a terminal (e.g. the provided terminal by anaconda).

Python support
**************
Currently supported python versions are 3.7 and 3.8.

Anaconda installation (recommended)
***********************************
Either create a new environment and install mspyeline e.g. by running:

.. code-block:: bash

    conda create --name mspypeline python=3.7 mspypeline -c conda-forge -c sschaum

This will create a new virtual environment called mspypeline with python 3.7 and install the mspypeline package.

Or install mspypeline in the base installation:

.. code-block:: bash

    conda install -c conda-forge -c sschaum mspypeline


pip(PyPI) installation
**********************
Either create a virtual environment and install mspypeline (this might differ based on your OS system).
Since the python version cannot be changed by pip make sure that the correct python version ist installed.

.. code-block:: bash

    python3 -m venv mypypeline
    source activate mspypeline
    pip install mspypeline

Or install into the base installation:

.. code-block:: bash

    pip install mspypeline


.. _activate-venv:

Activating virtual environment
******************************
The commands might differ based on OS systems and versions (e.g. python, pip or conda).

Anaconda
^^^^^^^^
New conda versions allow to run:

.. code-block:: bash

    conda activate [env-name]
    # activate the conda env called mspypeline
    conda activate mspypeline

pip
^^^

.. code-block:: bash

    source activate mspypeline



Obtaining Sources
*****************
Get the source code by cloning the github project:

.. code-block:: bash

    git clone https://github.com/siheming/mspypeline.git

Download Sources from `PyPI <https://pypi.org/project/mspypeline/>`__.

Download Sources from `conda <https://anaconda.org/sschaum/mspypeline>`__.


Dependencies
************
- `numpy <https://numpy.org/>`__ >= 1.17.4
- `pandas <https://pandas.pydata.org/>`__ >= 0.25.3
- `scipy <https://www.scipy.org/>`__ >= 1.3.1
- `matplotlib <https://matplotlib.org/>`__ >= 3.1.1
- `scikit-learn <https://scikit-learn.org/stable/>`__ >= 0.22.1
- tzlocal >= 2.0.0
- ruamel_yaml >= 0.15.46
- matplotlib-venn >= 0.11.5
- adjusttext >= 0.7.3.1
- plotly >= 4.6.0

Optional Dependencies for R packages
************************************
some plots might require additional R packages to be installed. Because of that additional dependencies are required for
those plots.


- rpy2=2.9.4
