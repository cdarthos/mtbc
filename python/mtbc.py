import random
from operator import itemgetter
import pandas as pd
import xmltodict
from Bio import Entrez


class MtbcRandom:
    # Mycobacterium tuberculosis complex 77643
    taxa_mycobacterium_tuberculosis_complex = '77643'

    # Mycobacterium canettii 78331
    taxa_mycobacterium_canettii = '78331'

    # Mycobacterium mungi 1844474
    taxa_mycobacterium_mungi = '1844474'

    # Mycobacterium orygis 1305738
    taxa_mycobacterium_orygis = '1305738'

    # Mycobacterium tuberculosis 1773
    taxa_mycobacterium_tuberculosis = '1773'

    def __init__(self,
                 select_mycobacterium_canettii=False,
                 select_mycobacterium_mungi=False,
                 select_mycobacterium_orygis=False,
                 select_mycobacterium_tuberculosis=False,
                 ):
        self.retmax = 1000000
        self.random_number = 10
        Entrez.email = 'A.N.Other@example.com'
        self.select_taxa = {self.taxa_mycobacterium_canettii: select_mycobacterium_canettii,
                            self.taxa_mycobacterium_mungi: select_mycobacterium_mungi,
                            self.taxa_mycobacterium_orygis: select_mycobacterium_orygis,
                            self.taxa_mycobacterium_tuberculosis: select_mycobacterium_tuberculosis}

    def construct_search_request(self):

        if True in self.select_taxa.values():
            request = " OR ".join('txid' + key + '[ORGN]' for key, value in self.select_taxa.items() if value)
            return request

        else:
            request = 'txid' + self.taxa_mycobacterium_tuberculosis_complex + '[ORGN]'
            return request

    def mtbc_random_search_id(self):
        retmax = 1000000
        handle = Entrez.esearch(db="sra",
                                term=self.construct_search_request(),
                                retmax=retmax,
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
    mtbc = MtbcRandom(select_mycobacterium_canettii=True, select_mycobacterium_mungi=True)
    print(mtbc.construct_search_request())
    print()
    print(mtbc.mtbc_random_search_id()[:10])
    print()
    print(mtbc.mtbc_random_get_accession_number())
    print()
