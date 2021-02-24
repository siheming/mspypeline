.. _file-readers:

File readers
============
File readers are required to translate the format of each specific file into an internal format. File readers scan
available data and create a data dictionary with keys to the data. Data stored on system hardware, is thus only loaded
on demand.


MaxQuant Reader
~~~~~~~~~~~~~~~~
| The MQReader prepares data from MaxQuant files to match the internal data format in order to provide the correct input
  for the plotters.
| Minimum requirement to start and work with the MQReader:

* proteinGroups.txt file.

.. note::
    Currently tested with MaxQuant version: 1.5+
.. warning::
    All files to be analyzed need to be in a directory called txt

The Quality Report will use these files (if available):
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* proteinGroups.txt
* peptides.txt
* parameters.txt
* summary.txt
* msScans.txt
* msmsScans.txt
* evidence.txt

If files are missing some plots of the :ref:`quality control report <plotters>` will not be created or might be empty.

The Plots will use these files:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Only the proteinGroups.txt file is required.