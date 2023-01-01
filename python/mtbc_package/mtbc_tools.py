import io
import json
import os
import pprint

import pandas
import requests
from Bio import Entrez, AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from pandas import CategoricalDtype
import logging
from pymongo import MongoClient
from .mtbc_ncbi import MtbcGetRandomSRA
import dendropy
from dendropy.interop import raxml


log = logging.getLogger("mtbc_tool")


class MtbcAcclistToFASTA:

    def __init__(self,
                 mtbc_get_random_sra: MtbcGetRandomSRA):

        # initial user variable
        self.fasta = None
        self.ml_tree = mtbc_get_random_sra.ml_tree
        self.list_length = mtbc_get_random_sra.list_length
        Entrez.email = mtbc_get_random_sra.email
        self.email = mtbc_get_random_sra.email
        self.select_taxa = vars(mtbc_get_random_sra.select_taxa)

        # initialize empty variable

        self.nj_tree = mtbc_get_random_sra.nj_tree
        self.ncbi_random_acc_list = mtbc_get_random_sra.ncbi_random_acc_list
        self.sample_list = mtbc_get_random_sra.sample_list
        self.ncbi_all_id = None
        self.ncbi_request_all_id = mtbc_get_random_sra.ncbi_request_all_id
        self.align_with_alignIO = mtbc_get_random_sra.align_with_alignIO
        self.df_mutation = mtbc_get_random_sra.df_mutation
        self._id = mtbc_get_random_sra._id
        self.sequence_dict = {'NC_000962.3': {}}
        # self.alignement = {}


        self.mtbc_request()
        log.info("self.sequence")
        log.debug(self.sequence_dict)

        self.reconstruct_sequence_to_fasta_file()

        # self.to_json_file()
        # self.to_db()

    def mtbc_request(self):
        ncbi_random_acc_list_len = len(self.ncbi_random_acc_list)
        index = 1
        log.info("mtbc_request")
        for sra in self.ncbi_random_acc_list:
            logging.info(str(index) + "/" + str(ncbi_random_acc_list_len))
            index += 1

            head = {'Content-Type': 'application/x-www-form-urlencoded',
                    'Host': 'gnksoftware.synology.me:30002'}
            parameters = {
                'search': '{"type":"boolean","quantifier":"allOf","filters":[{"type":"terms","terms":["' +
                          str(sra) +
                          '"],"attribute":"publicId.keyword","quantifier":"anyOf"}]}',
                'alignmentParams': '{"types":["snp"],"includeDataInComment":true}'}
            url = 'http://gnksoftware.synology.me:30002/strains/alignment'
            r = requests.post(url, parameters, head)

            if len(r.text) != 0:
                log.info(r.text[1:].split("\n")[0])
                self.sequence_dict[r.text[1:].split("\n")[0]] = {}
                for diff in r.text.split("\n")[2:]:
                    if len(diff) != 0:
                        ###############
                        # premiere etape, on selectionne juste les mutations (pas d'indels/insertion)
                        ###############
                        if len(diff.split(":")[3]) == 1:
                            self.sequence_dict['NC_000962.3'][diff.split(":")[1]] = diff.split(":")[2]
                            self.sequence_dict[r.text[1:].split("\n")[0]][diff.split(":")[1]] = diff.split(":")[3]

    def reconstruct_sequence_to_fasta_file(self):

        cat_type = CategoricalDtype(categories=list("ATCG"))
        df_mutation = pandas.DataFrame.from_dict(self.sequence_dict, dtype=cat_type)
        log.debug(df_mutation.info(memory_usage="deep"))

        handle = io.StringIO()
        for column in df_mutation.columns:
            df_mutation[column] = df_mutation[column].fillna(df_mutation['NC_000962.3'], axis=0)
            handle.write(">" + column + "\n")
            handle.write("".join(df_mutation[column].to_list()) + "\n")


        self.fasta = handle.getvalue()


    def to_db(self):
        self.ncbi_all_id = None
        try:
            client = MongoClient('mongodb://localhost:27017/')
            db_mtbc = client.db_mtbc
            request_data = db_mtbc.request_data
        except:
            log.error("error to connect mongo db")

        get_id = request_data.update_one(
            {"_id": self._id},
            {"$set": self.to_json()})
        log.info(get_id)
        client.close()

    def to_json(self):
        self.ncbi_all_id = None
        return self.__dict__


class MtbcTree:

    def create_nj_tree_static(id, fasta):
        log.info("create_nj_tree_static")
        calculator = DistanceCalculator('identity')
        handle_fasta = io.StringIO(fasta)
        align = AlignIO.read(handle_fasta, 'fasta')
        #align = AlignIO.read('alignement/{0}'.format(id), 'fasta')
        dist_matrix = calculator.get_distance(align)
        log.info("dist_matrix")
        log.debug(dist_matrix)
        constructor = DistanceTreeConstructor()
        nj_tree = constructor.nj(dist_matrix)
        handle = io.StringIO()
        Phylo.write(nj_tree, handle, "newick", )
        resultat = handle.getvalue()
        handle.close()
        return resultat

    def create_ml_tree_static(id, fasta):
        log.info("create_ml_tree_static")
        data = dendropy.DnaCharacterMatrix.get(
            data = fasta, schema="fasta")
        log.info("data")
        log.info(data)
        rx = raxml.RaxmlRunner()
        ml_tree = rx.estimate_tree(
            char_matrix=data,
            raxml_args=["--no-bfgs"])
        ml_tree_str = ml_tree.as_string(schema="newick")
        log.info(ml_tree_str)
        return ml_tree_str

    def to_json(self):
        self.ncbi_all_id = None
        return self.__dict__
