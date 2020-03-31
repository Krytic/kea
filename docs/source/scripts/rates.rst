Event Rate Calculation Scripts
==============================

.. toctree::
   :glob:
   :hidden:


.. py:function:: python calculate2dEventRates.py -s 2DSFRD_FILE \
                                                 -t TIME_FILE \
                                                 -i DATA_FOLDER \
                                                 -o OUTPUT_FILE

   Script to calculate the 2D event rate (using pickles SFRD and BPASS)
   :param str -i 2DSFRD: A string pointing to the 2d SFRD file
   :param str -t TIME_FILE: The time relations file.
   :param str -i DATA_FOLDER: Folder containing the BPASS models
   :param str -o OUTPUT_FILE: Filename of the outputted 2D event rate pickel.
