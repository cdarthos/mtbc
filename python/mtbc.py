import random
from operator import itemgetter

from python.mtbc_package.mtbc_ncbi import MtbcGetRandomSRA

import pandas
import pandas as pd
import requests
import xmltodict
from Bio import Entrez
from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from pandas.api.types import CategoricalDtype

from python.mtbc_package.mtbc_tools import MtbcAcclistToFASTA

if __name__ == "__main__":
    # mtbc = MtbcGetRandomSRA(debug=True, list_length=1000,select_mycobacterium_canettii=True,
    # select_mycobacterium_mungi=True, select_mycobacterium_orygis=True)
    mtbc = MtbcGetRandomSRA(debug=True, list_length=10, id1=2)
    print(mtbc.ncbi_random_acc_list)
    #mtbc_2 = MtbcAcclistToFASTA(mtbc)

