.. currentmodule:: mspypeline

Modules
=======
import with:

.. code-block:: python

  from mspypeline import MedianNormalizer, QuantileNormalizer, TailRobustNormalizer, interpolate_data
  from mspypeline.modules.Normalization import BaseNormalizer
  from mspypeline import DataNode, DataTree


Normalization
~~~~~~~~~~~~~

.. autoclass:: mspypeline.modules.Normalization.BaseNormalizer
   :members:

.. autoclass:: MedianNormalizer
   :members:

.. autoclass:: QuantileNormalizer
   :members:

.. autoclass:: TailRobustNormalizer
   :members:

.. autofunction:: interpolate_data

DataStructure
~~~~~~~~~~~~~

.. autoclass:: DataNode
   :members:

.. autoclass:: DataTree
   :members:
