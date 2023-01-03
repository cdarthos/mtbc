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

FORMAT = "%(levelname)s:%(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)







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
                 ncbi_list_length=100,
                 email='A.N.Other@example.com',
                 all_id_to_acc=False,
                 target_list_length=10,
                 snp_select=[],
                 snp_reject=[]):

        # initial user variable
        self.target_list_length = target_list_length
        self.snp_select = snp_select
        self.snp_reject = snp_reject
        self.sample_list = None
        self.all_id_to_acc = all_id_to_acc
        self.outgroup = outgroup
        self.ncbi_list_length = ncbi_list_length
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
        self.final_acc_list = None
        self.final_acc_list_length = None
        self.fasta = None

        self._id = str(uuid.uuid4()) + "_" + str(ncbi_list_length)

        # main program
        self.construct_search_request()
        logging.info("self.ncbi_request_all_id")
        self.get_all_id()
        logging.info("self.ncbi_all_id")

        self.select_random_ncbi_id_number()
        logging.info("self.sample_list")

        self.acc_number_from_ncbi_id()
        logging.info("self.acc_list")




    def construct_search_request(self):
        logging.info("construct_search_request")
        if True in self.select_taxa.values():
            self.ncbi_request_all_id = " OR ".join(
                'txid' + key + '[ORGN]' for key, value in self.select_taxa.items() if value)

        else:
            self.ncbi_request_all_id = 'txid' + self.taxa_mycobacterium_tuberculosis_complex + '[ORGN]'

    def get_all_id(self):
        logging.info("get_all_id")
        retmax = 1000000
        handle = Entrez.esearch(db="sra",
                                term=self.ncbi_request_all_id,
                                retmax=retmax,
                                retstart=0)
        record = Entrez.read(handle)
        handle.close()
        self.ncbi_all_id = record['IdList']

    def get_random_acc_list_DEV(self):
        logging.info("get_random_acc_list_DEV")
        shuffled_id_list = shuffle(self.ncbi_all_id)

    def select_random_ncbi_id_number(self):
        logging.info("select_random_ncbi_id_number")
        self.ncbi_random_id_list = random.choices(self.ncbi_all_id,
                                                  k=self.ncbi_list_length)

    def acc_number_from_ncbi_id(self):
        if self.all_id_to_acc:
            self.no_limit_acc_number_from_ncbi_id()
        else:
            self.limit_acc_number_from_ncbi_id()

    def no_limit_acc_number_from_ncbi_id(self):
        handle_epost = Entrez.epost(db="sra", id=",".join(map(str, self.ncbi_all_id)))
        record_epost = Entrez.read(handle_epost)
        handle_epost.close()
        #logging.debug(record_epost)
        logging.debug(record_epost['WebEnv'])

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
        logging.debug(record_esummary)
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