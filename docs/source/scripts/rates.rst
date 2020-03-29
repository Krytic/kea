Event Rate Calculation Scripts
==============================

.. toctree::
   :glob:
   :hidden:


.. py:function:: python plotRatesnoZ.py -i DATA_FOLDER -o OUTPUT_FOLDER

   A script to plot the rates of SNe & Compact mergers without taking metallicity
   into account.

   :param str -i DATA_FOLDER: A string pointing to the data folder containing
    BPASS models.
   :param str -o OUTPUT_FOLDER: A string pointing to a folder to output the images.

   :return: A collection of plots is outputted:

    * The Stellar Formation Rate

      1. Actual rate
      2. Interpolated rate

    * The BPASS event rates per solar mass
    * The Event rates per :math:`Gpc^3`
