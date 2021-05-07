.. _get-started:

.. currentmodule:: mspypeline

.. role:: bash(code)
   :language: bash

Getting started
===============

``mspypeline`` can be used in three ways. If ``mspypeline`` was installed in a virtual environment make sure to first
activate the virtual environment (:ref:`conda env activation <conda-installation>`,
:ref:`pip env activation <pip-installation>`).

#. | starting the GUI via the command-line:
   | (see the :ref:`GUI Qickstart <gui-quickstart>` for next steps)

   .. code-block:: bash

        # either: call the python module mspypeline and start the GUI
        python -m mspypeline --gui

        # or: call the python module mspypeline, start the GUI and receive information from the logger
        python -m mspypeline --gui --loglevel DEBUG


#. | using it in python by importing:
   | (see the :ref:`Python Quickstart <python-quickstart>`, :ref:`Examples <examples>` or :ref:`API reference <api>`
     for next steps.)

   .. code-block:: python

      import mspypeline

#. | using it purely on the command-line:
   | (see the :ref:`command-line quickstart <cli-quickstart>` for next steps.)

   .. code-block:: bash

      python -m mspypeline

.. _`gui-quickstart`:

Usage with GUI
~~~~~~~~~~~~~~
.. warning::
    All files to be analyzed need to be in a directory called txt (as originally created by MaxQuant)
The GUI is divided into different parts to guide you through the general :ref:`workflow`.

#. Start by selecting a directory at the top under the **"Dir to analyze"** option.

#. The **"Yaml file"** option selects where the :ref:`settings` for the analysis are from. This can either be the
   :ref:`"default" settings <default-yaml>`, or the :ref:`"file" settings <default-yaml>`, if they are available.

#. If a :ref:`quality control report <max-quant-report>` should be generated, click the **"Create Report"** button.

#. A :ref:`normalization method <hyperparameter>` can be selected that will be applied to the data. Plot options below
   the **"Normalization plots"** can help to decide for a normalizer.

#. Generally, the GUI is structured in columns and rows like a table. Create a plot by toggling the checkbox in the left
   column for the desired plot (see all available plots in :ref:`analysis options <Plot-Options>`). In the same row,
   select one or more :ref:`intensities <hyperparameter>` and :ref:`levels <analysis-design>` for which the plot should
   be created. Currently selected intensity and level options in the drop down menu are indicated by a check mark.

#. The **"Update"** button writes the selected options to a :ref:`configuration file <default-yaml>` and the **"Start"**
   button creates all plots that were selected.

#. For proper analysis, the :ref:`naming convention <analysis-design>` has to be followed. If the naming convention was
   was violated an auxiliary file will be created as described in
   the :ref:`here <sample-mapping>`.


.. _python-quickstart:


Usage with Python
~~~~~~~~~~~~~~~~~
| For the usage of ``mspypeline`` in a python script, the package or modules of the package have to be first imported.
  An analysis can be created in different ways.
| A simple analysis in which the :class:`~MSPInitializer` is called to read in the data and settings and to create a
  :class:`~MaxQuantPlotter` which enables the calculation and plotting of the data might look something like this:

   .. code-block:: python

      from mspypeline import MSPInitializer, MaxQuantPlotter

      # create an initializer and read the data (txt folder directory) and settings
      init = MSPInitializer("path/to/data")
      init.init_config()
      init.read_data()

      # create a plotter from an initializer instance and create all those results specified in the settings
      plotter = MaxQuantPlotter.from_MSPInitializer(init)
      plotter.create_results()


For more code examples see the :ref:`Examples <examples>`.


.. _cli-quickstart:

Usage with command-line
~~~~~~~~~~~~~~~~~~~~~~~
It is recommended to provide a :ref:`yaml file <settings>`, where most settings are already configured. It is **not**
possible at the moment to specify all arguments via the command-line, and will only be implemented if there is any
demand.

