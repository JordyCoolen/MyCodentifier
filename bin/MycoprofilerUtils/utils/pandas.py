#!/usr/bin/env python

"""
    Functions that use pandas and are generic to use
"""
import pandas as pd

def read_report(inputFile):
    """
        Loads centrifuge report file to pandas

        :param inputFile: Full path to inputFile (*.report)
        :return df: pandas dataframe
    """

    df = pd.read_csv(inputFile, sep='\t', header=0) # file has a header
    return(df)