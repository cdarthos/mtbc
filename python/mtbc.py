import random
from operator import itemgetter
import pandas as pd
import xmltodict
from Bio import Entrez


class MtbcRandom:
    def __init__(self):
        self.retmax = 1000000
        self.random_number = 10
        Entrez.email = 'A.N.Other@example.com'
        self.base_taxa = 'txid77643[ORGN]'

    def mtbc_random_search_id(self):
        handle = Entrez.esearch(db="sra",
                                term=self.base_taxa,
                                retmax=self.retmax,
                                retstart=0)
        record = Entrez.read(handle)
        handle.close()

        return record['IdList']

    def mtbc_random_get_accession_number(self):
        id_list = self.mtbc_random_search_id()
        sample_list = random.choices(id_list,
                                     k=self.random_number)

        handle_epost = Entrez.epost(db="sra",
                                    id=",".join(map(str, sample_list)),
                                    retmax=self.retmax)
        record_epost = Entrez.read(handle_epost)
        handle_epost.close()

        handle_esummary = Entrez.esummary(db="sra",
                                          webenv=record_epost['WebEnv'],
                                          query_key=record_epost['QueryKey'],
                                          retmax=self.retmax)
        record_esummary = Entrez.read(handle_esummary)
        handle_esummary.close()

        runs = list(map(itemgetter('Runs'),
                        record_esummary))

        runs_1 = list(map(lambda run: "<root>" + run + "</root>",
                          runs))

        runs_acc = list(map(xmltodict.parse, runs_1))

        acc_list = pd.json_normalize(runs_acc)['root.Run.@acc']
        return acc_list.to_list()


if __name__ == "__main__":
    mtbc = MtbcRandom()
    print(mtbc.mtbc_random_search_id()[:10])
    print()
    print(mtbc.mtbc_random_get_accession_number())
