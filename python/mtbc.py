import random
from operator import itemgetter

import pandas as pd
import xmltodict
from Bio import Entrez


def mtbc_random():
    retmax = 1000000
    Entrez.email = "Your.Name.Here@example.org"

    handle = Entrez.esearch(db="sra", term="txid77643[ORGN]", retmax=retmax, retstart=0)
    record = Entrez.read(handle)
    handle.close()

    sample_list = random.choices(record['IdList'], k=10)

    handle_epost = Entrez.epost(db="sra", id=",".join(map(str, sample_list)), retmax=retmax)
    record_epost = Entrez.read(handle_epost)
    handle_epost.close()

    handle_esummary2 = Entrez.esummary(db="sra", webenv=record_epost['WebEnv'], query_key=record_epost['QueryKey'],
                                       retmax=retmax)
    record_esummary2 = Entrez.read(handle_esummary2)
    handle_esummary2.close()

    runs = list(map(itemgetter('Runs'), record_esummary2))

    runs_1 = list(map(lambda run: "<root>" + run + "</root>", runs))

    runs_acc = list(map(xmltodict.parse, runs_1))

    acc_list = pd.json_normalize(runs_acc)['root.Run.@acc']
    return acc_list


if __name__ == "__main__":
    mtbc_random()
