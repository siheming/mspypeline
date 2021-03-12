.. _get-started:

.. currentmodule:: mspypeline

.. role:: bash(code)
   :language: bash

Getting started
===============

``mspypeline`` can be used in three ways. If ``mspypeline`` was installed in a virtual environment make sure
:ref:`to activate it first <activate-venv>`.

#. using it in python by importing:

   .. code-block:: python

      import mspypeline

   See the :ref:`Python Quickstart <python-quickstart>`, :ref:`Examples <examples>` or :ref:`API reference <api>`.

#. starting the GUI via the command-line:

   .. code-block:: bash

      python -m mspypeline --gui

   See the :ref:`GUI Qickstart <gui-quickstart>`.

#. using it purely on the command-line:

   .. code-block:: bash

      python -m mspypeline

   See the :ref:`command-line quickstart <cli-quickstart>`.


.. _python-quickstart:

Usage with Python
~~~~~~~~~~~~~~~~~
A simple analysis might look something like this:

   .. code-block:: python

      from mspypeline import MSPInitializer, MaxQuantPlotter

      # create an initializer and read the data
      init = MSPInitializer("path/to/data")
      init.init_config()
      init.read_data()

      plotter = MaxQuantPlotter.from_MSPInitializer(init)
      plotter.create_results()


For more code examples see the :ref:`Examples <examples>`.

.. _`gui-quickstart`:

Usage with GUI
~~~~~~~~~~~~~~
The GUI is divided into different parts to guide you through the general :ref:`workflow`.

* Start by selecting a directory at the top under the "Dir to analyze" option.
.. warning::
    All files to be analyzed need to be in a directory called txt
* The "Yaml file" option selects where the settings for the analysis are from. This can either be the "default"
  settings, or the "file" settings, if they are available.
* If a Quality control report should be generated click the "Create Report" button.
* A Normalizer can be selected, which is further explained in the :ref:`workflow`.
  The plots under "Normalization plots" can help to decide for a normalizer.
* Create a plot by toggling the checkbox for the plot. Then select all intensities and levels
  for which the plot should be created. Currently selected options are indicated by a check mark.
* The "Update" button writes the selected options to a configuration file, the "Start" button additionally creates
  all plots that were selected. Also if the naming convention was not followed a file will be created as described in
  the :ref:`analysis-design`.

.. _cli-quickstart:

Usage with command-line
~~~~~~~~~~~~~~~~~~~~~~~
It is recommended to provide a yaml file, where most settings are already configured. It is **not** possible
at the moment to specify all arguments via the command-line, and will only be implemented if there is any demand.
