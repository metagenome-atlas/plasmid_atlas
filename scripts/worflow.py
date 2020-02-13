import os
import pandas as pd
import logging


def get_all_samples(SampleTableFile):


    if not os.path.exists(SampleTableFile):
        logging.error(f"The configuration says I have to look for SampleTable at {SampleTableFile}\n"
                      "But this file doesn't exists! The SampleTable should be a table in the format: \n"
                      "Samples\tOPTIONAL_HEADER\n"
                      "Sample1\t\n"
                      "Sample2\t\n"
                      )
        raise IOError("sampletable doesn't exist")

    SampleTable = pd.read_csv(SampleTableFile, sep='\t', index_col=0)

    return list(SampleTable.index)
