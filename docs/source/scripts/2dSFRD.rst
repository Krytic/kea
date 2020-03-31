SFRD over metallicity
=====================


.. toctree::
   :glob:
   :hidden:


.. py:function:: python calculate2dSFRD.py -i DATA_FILE-t TIME_FILE -o OUTPUT_FILE

   A script to calculate the Stellar Formation Rate Density over time per BPASS
   metallicity bin and output a pickle file.

   :param str -i DATA_FILE: The filename of the cosmological simulations.
   :param str -t TIME_FILE: The filename of the time information of the
    cosmological simulation.
   :param str -o OUTPUT_FILE: The filename of the output file.

   :return: A pickle containing a 2D pandas DataFrame with the SFRD over time
    and over metallicity.
