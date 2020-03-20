import gzip
import os
from hoki import load

def gunzip(source_filepath, dest_filepath, block_size=65536):
    with gzip.open(source_filepath, 'rb') as s_file, \
            open(dest_filepath, 'wb') as d_file:
        while True:
            block = s_file.read(block_size)
            if not block:
                break
            else:
                d_file.write(block)
        d_file.write(block)

def packnload(file):
    gunzip("../bpass_v2.2.1_imf135_300/"+file+".gz", file)
    out = load.model_output(file)
    os.remove(file)
    return out
