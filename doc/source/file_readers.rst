.. _file-readers:

File readers
============
File readers are required to translate the format of each specific file into an internal format.


Max Quant Reader
~~~~~~~~~~~~~~~~

Minimum requirements: proteinGroups.txt file.

.. note::
    Currently tested with Max Quant version: 1.5+
.. warning::
    All files to be analyzed need to be in a directory called txt

The Report will use these files (if available):
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* proteinGroups.txt
* peptides.txt
* parameters.txt
* summary.txt
* msScans.txt
* msmsScans.txt
* evidence.txt

If files are missing some plots will not be created or might be empty.

The Plots will use these files:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Only the proteinGroups.txt file is required.