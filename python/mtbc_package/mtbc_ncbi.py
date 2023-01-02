import json
import logging
import random
import uuid
# from json import JSONEncoder
from operator import itemgetter
import pandas as pd
import xmltodict
from Bio import Entrez
from pymongo import MongoClient
from sklearn.utils import shuffle


log = logging.getLogger("mtbc_ncbi")
log.setLevel(logging.DEBUG)
class MtbcData:
    def __init__(self):

        # initial user variable

        self.sample_list = None
        self.all_id_to_acc = None
        self.outgroup = None
        self.list_length = None
        Entrez.email = None
        self.email = None
        self.select_taxa = {self.taxa_mycobacterium_canettii: None,
                            self.taxa_mycobacterium_mungi: None,
                            self.taxa_mycobacterium_orygis: None,
                            self.taxa_mycobacterium_tuberculosis: None}

        # initialize empty variable
        self.ml_tree = None
        self.nj_tree = None
        self.ncbi_random_acc_list = []
        self.ncbi_random_id_list = None
        self.ncbi_all_id = None
        self.ncbi_request_all_id = None
        self.align_with_alignIO = None
        self.df_mutation = None
        self.sequence_dict = {}
        #self.alignement = {}
        self.fasta = ""

        self._id = ""



class MtbcGetRandomSRA:
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
                 outgroup='',
                 list_length=10,
                 email='A.N.Other@example.com',
                 all_id_to_acc=False):

        # initial user variable
        self.sample_list = None
        self.all_id_to_acc = all_id_to_acc
        self.outgroup = outgroup
        self.list_length = list_length
        Entrez.email = email
        self.email = email
        self.select_taxa = {self.taxa_mycobacterium_canettii: select_mycobacterium_canettii,
                            self.taxa_mycobacterium_mungi: select_mycobacterium_mungi,
                            self.taxa_mycobacterium_orygis: select_mycobacterium_orygis,
                            self.taxa_mycobacterium_tuberculosis: select_mycobacterium_tuberculosis}

        # initialize empty variable
        self.ml_tree = None
        self.nj_tree = None
        self.ncbi_random_acc_list = []
        self.ncbi_random_id_list = None
        self.ncbi_all_id = None
        self.ncbi_request_all_id = None
        self.align_with_alignIO = None
        self.df_mutation = None
        self.sequence_dict = {'NC_000962.3': {}}
        self.fasta = None

        self._id = str(uuid.uuid4()) + "_" + str(list_length)

        # main program
        self.construct_search_request()
        log.info("self.ncbi_request_all_id")
        self.get_all_id()
        log.info("self.ncbi_all_id")

        self.select_random_ncbi_id_number()
        log.info("self.sample_list")

        self.acc_number_from_ncbi_id()
        log.info("self.acc_list")

        log.info(len(self.ncbi_random_acc_list))


    def construct_search_request(self):
        log.info("construct_search_request")
        if True in self.select_taxa.values():
            self.ncbi_request_all_id = " OR ".join(
                'txid' + key + '[ORGN]' for key, value in self.select_taxa.items() if value)

        else:
            self.ncbi_request_all_id = 'txid' + self.taxa_mycobacterium_tuberculosis_complex + '[ORGN]'

    def get_all_id(self):
        log.info("get_all_id")
        retmax = 1000000
        handle = Entrez.esearch(db="sra",
                                term=self.ncbi_request_all_id,
                                retmax=retmax,
                                retstart=0)
        record = Entrez.read(handle)
        handle.close()
        self.ncbi_all_id = record['IdList']

    def get_random_acc_list_DEV(self):
        log.info("get_random_acc_list_DEV")
        shuffled_id_list = shuffle(self.ncbi_all_id)

    def select_random_ncbi_id_number(self):
        log.info("select_random_ncbi_id_number")
        self.ncbi_random_id_list = random.choices(self.ncbi_all_id,
                                                  k=self.list_length)

    def acc_number_from_ncbi_id(self):
        if self.all_id_to_acc:
            self.no_limit_acc_number_from_ncbi_id()
        else:
            self.limit_acc_number_from_ncbi_id()

    def no_limit_acc_number_from_ncbi_id(self):
        handle_epost = Entrez.epost(db="sra", id=",".join(map(str, self.ncbi_all_id)))
        record_epost = Entrez.read(handle_epost)
        handle_epost.close()
        log.debug(record_epost)
        log.debug(record_epost['WebEnv'])

        for start in range(0, len(self.ncbi_all_id), 5000):
            handle_esummary = Entrez.esummary(db="sra",
                                              query_key=record_epost['QueryKey'],
                                              WebEnv=record_epost['WebEnv'],
                                              retsart=start,
                                              retmax=5000
                                              )
            record_esummary = Entrez.read(handle_esummary)
            handle_esummary.close()
            runs = list(map(itemgetter('Runs'),
                            record_esummary))
            runs_1 = list(map(lambda run: "<root>" + run + "</root>",
                              runs))
            runs_acc = list(map(xmltodict.parse, runs_1))
            acc_list = pd.json_normalize(runs_acc)['root.Run.@acc']
            self.ncbi_random_acc_list = self.ncbi_random_acc_list + acc_list.to_list()

    def limit_acc_number_from_ncbi_id(self):
        handle_esummary = Entrez.esummary(db="sra",
                                          id=",".join(map(str, self.ncbi_random_id_list))
                                          )
        record_esummary = Entrez.read(handle_esummary)
        logging.info("record_esummary")
        log.info(record_esummary)
        handle_esummary.close()

        runs = list(map(itemgetter('Runs'),
                        record_esummary))

        runs_1 = list(map(lambda run: "<root>" + run + "</root>",
                          runs))

        runs_acc = list(map(xmltodict.parse, runs_1))

        acc_list = pd.json_normalize(runs_acc)['root.Run.@acc']
        self.ncbi_random_acc_list = acc_list.dropna().to_list()

    def add_outgroup(self):
        self.ncbi_random_acc_list.append(self.outgroup)



    def to_json(self):
        self.ncbi_all_id = None
        return self.__dict__

    def to_db(self):
        self.ncbi_all_id = None
        try:
            client = MongoClient('mongodb://localhost:27017/')
            db_mtbc = client.db_mtbc
            request_data = db_mtbc.request_data
        except:
            log.error("error to connect mongo db")
        get_id = request_data.insert_one(self.to_json()).inserted_id
        log.info(get_id)
        client.close()


if __name__ == "__main__":
    mtbc_inst1 = MtbcGetRandomSRA(list_length=10)
    # print(mtbc.to_json())
