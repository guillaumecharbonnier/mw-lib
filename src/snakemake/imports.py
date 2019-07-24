# Imports
# I want to be able to print to stderr instead of stdout so --rulegraph and --dag outputs are not 'corrupted' by my printed messages.
from __future__ import print_function
import sys

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

#from Bio import Entrez
#from snakemake.utils import R # To use R inside rules
#from snakemake.utils import report # To write report in reStructuredText compiled into html.
#from snakemake.remote.NCBI import RemoteProvider as NCBIRemoteProvider
#NCBI = NCBIRemoteProvider(email="yourEmail@yourProvider") # email required by NCBI to prevent abuse
import re
import os
import datetime
import os.path
import csv
import pandas
import numpy
import glob # To look for files using UNIX regex.
import getpass # To get user names to send onsuccess email to the correct recipient.


