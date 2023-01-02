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

FORMAT = "%(levelname)s:%(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)

class MtbcAcclistToFASTA:

    def __init__(self,
                 mtbc_get_random_sra: MtbcGetRandomSRA):

        # initial user variable
        self.fasta = None
        self.ncbi_random_acc_list = mtbc_get_random_sra.ncbi_random_acc_list
        self.df_mutation = mtbc_get_random_sra.df_mutation
        self._id = mtbc_get_random_sra._id


        self.sequence_dict = mtbc_get_random_sra.sequence_dict
        if len(self.sequence_dict) < 2:
            self.mtbc_request()
        logging.info("self.sequence")
        logging.debug(self.sequence_dict)
        self.reconstruct_sequence_to_fasta_file()


    def mtbc_request(self):
        ncbi_random_acc_list_len = len(self.ncbi_random_acc_list)
        index = 1
        logging.info("mtbc_request")

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
                logging.info(r.text[1:].split("\n")[0])
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
        logging.debug(df_mutation.info(memory_usage="deep"))

        handle = io.StringIO()
        for column in df_mutation.columns:
            df_mutation[column] = df_mutation[column].fillna(df_mutation['NC_000962.3'], axis=0)
            handle.write(">" + column + "\n")
            handle.write("".join(df_mutation[column].to_list()) + "\n")

        self.fasta = handle.getvalue()

    def to_json(self):
        return self.__dict__


class MtbcTree:

    def create_nj_tree_static(id, fasta):
        logging.info("create_nj_tree_static")
        calculator = DistanceCalculator('identity')
        handle_fasta = io.StringIO(fasta)
        align = AlignIO.read(handle_fasta, 'fasta')
        dist_matrix = calculator.get_distance(align)
        logging.info("dist_matrix")
        logging.debug(dist_matrix)
        constructor = DistanceTreeConstructor()
        nj_tree = constructor.nj(dist_matrix)
        handle = io.StringIO()
        Phylo.write(nj_tree, handle, "newick", )
        resultat = handle.getvalue()
        handle.close()
        return resultat

    def create_ml_tree_static(id, fasta):
        logging.info("create_ml_tree_static")
        data = dendropy.DnaCharacterMatrix.get(
            data=fasta, schema="fasta")
        logging.info("data")
        rx = raxml.RaxmlRunner()
        ml_tree = rx.estimate_tree(
            char_matrix=data,
            raxml_args=["--no-bfgs"])
        ml_tree_str = ml_tree.as_string(schema="newick")
        logging.info(ml_tree_str)
        return ml_tree_str

    def to_json(self):
        return self.__dict__
