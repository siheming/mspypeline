.. MSPypeline documentation master file, created by
   sphinx-quickstart on Fri Apr 10 12:03:31 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MSPypeline's documentation!
======================================

``mspypeline`` is a package for analysing mass spec results.
The source code is available on `github <https://github.com/siheming/mspypeline/>`__.

.. figure:: ./_static/flow_chart_final_docu.pdf
    :width: 300
    :align: right

Usage
=====

| ``mspypeline`` was developed to consolidate the proteomics data analysis of **MaxQuant output tables** by providing a
  tool to analyze complex datasets in a standardized manner with minimal effort and to eliminate the chance of human
  error or obscuring variability during data analysis.
| Designed with an intuitive and concise graphical user interface (GUI), ``mspypeline`` offers researchers,
  unfamiliar with programming or data analysis, the opportunity to explore and visualize their data independently and in
  a time-effective manner. A precisely structured workflow (see figure to the left) is prescribed following logical and
  systematic steps of analysis that starts with a quality control of the data, followed by the assessment and choice of
  data pre-processing operations to finally allow optimal exploratory analysis. By automizing the calculation and
  generation of versatile figures ``mspypeline`` can create a comprehensive and conclusive data analysis within minutes.
  Simultaneously, the more experienced user may interact closer with the ``mspypeline`` package to perform advanced
  analysis leveraging the plethora of customization options.
| At the same time, the software aims to establish a standardized and reproducable analysis procedure, which is
  supported by automated logging of all analysis settings and saving them to a seperate configuration file.
| Thus, ``mspypeline`` provides a platform that supports users with their proteomics data analysis by giving insight
  into the data, offering parameter adaptation when needed, generating custom figures and by reaching biological
  relevant conclusions through differential expression analysis and hypothesis testing.
| Please refer to the other parts of the documentation for :ref:`installation <installation>`,
  :ref:`how to get started <get-started>`, or use the search.


Citation
========
In preparation.


.. toctree::
   :maxdepth: 2
   :caption: Contents


   installation
   get_started
   workflow
   file_readers
   analysis_options
   settings_and_configuration
   package_config
   examples
   gallery
   API <api_reference/index>



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
