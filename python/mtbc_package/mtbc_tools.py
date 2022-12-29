import json
import os

import pandas
import requests
from Bio import Entrez, AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from pandas import CategoricalDtype
import logging
from .custom_encoder import customEncoder
from .mtbc_ncbi import MtbcGetRandomSRA
import dendropy
from dendropy.interop import raxml

log = logging.getLogger("mtbc_tool")


class MtbcAcclistToFASTA:

    def __init__(self,
                 mtbc_get_random_sra: MtbcGetRandomSRA):

        # initial user variable
        self.ml_tree = mtbc_get_random_sra.ml_tree
        self.list_length = mtbc_get_random_sra.list_length
        Entrez.email = mtbc_get_random_sra.email
        self.email = mtbc_get_random_sra.email
        self.select_taxa = mtbc_get_random_sra.select_taxa

        # initialize empty variable
        self.nj_tree = mtbc_get_random_sra.nj_tree
        self.ncbi_random_acc_list = mtbc_get_random_sra.ncbi_random_acc_list
        self.sample_list = mtbc_get_random_sra.sample_list
        self.ncbi_all_id = None
        self.ncbi_request_all_id = mtbc_get_random_sra.ncbi_request_all_id
        self.align_with_alignIO = mtbc_get_random_sra.align_with_alignIO
        self.df_mutation = mtbc_get_random_sra.df_mutation
        self.id = mtbc_get_random_sra.id
        self.sequence_dict = {'NC_000962.3': {}}
        self.alignement = {}

        with open('alignement/{0}'.format(self.id), 'w') as writer:
            writer.write("")
        self.mtbc_request()
        log.info("self.sequence")
        log.debug(self.sequence_dict)
        self.reconstruct_sequence_to_fasta_file()

        self.to_json_file()

    def mtbc_request(self):
        acc_list_len = len(self.acc_list)
        index = 1
        log.info("mtbc_request")
        for sra in self.acc_list:
            logging.info(str(index) + "/" + str(acc_list_len))
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

        log.info(df_mutation.info(memory_usage="deep"))
        with open('alignement/{0}'.format(self.id), 'w') as writer:
            for column in df_mutation.columns:
                df_mutation[column] = df_mutation[column].fillna(df_mutation['NC_000962.3'], axis=0)
                writer.writelines(">" + column + "\n")
                writer.writelines("".join(df_mutation[column].to_list()) + "\n")

    def to_json_file(self):
        self.ncbi_all_id = None
        with open("request/{0}".format(self.id), 'w') as json_request:
            json.dump(self, json_request, cls=customEncoder)


class MtbcTree:
    def __init__(self,
                 mtbc_get_random_sra: MtbcGetRandomSRA):
        # initial user variable
        self.ml_tree = mtbc_get_random_sra.ml_tree
        self.list_length = mtbc_get_random_sra.list_length
        Entrez.email = mtbc_get_random_sra.email
        self.email = mtbc_get_random_sra.email
        self.select_taxa = mtbc_get_random_sra.select_taxa

        # initialize empty variable
        self.nj_tree = mtbc_get_random_sra.nj_tree
        self.ncbi_random_acc_list = mtbc_get_random_sra.ncbi_random_acc_list
        self.sample_list = mtbc_get_random_sra.sample_list
        self.ncbi_all_id = None
        self.ncbi_request_all_id = mtbc_get_random_sra.ncbi_request_all_id
        self.align_with_alignIO = mtbc_get_random_sra.align_with_alignIO
        self.df_mutation = mtbc_get_random_sra.df_mutation
        self.id = mtbc_get_random_sra.id
        self.sequence_dict = mtbc_get_random_sra.sequence_dict
        self.alignement = mtbc_get_random_sra.alignement
        if not os.path.exists('alignement/{0}'.format(self.id)):
            MtbcAcclistToFASTA(self)
        self.align_reconstruct()

    def align_reconstruct(self):
        self.align_with_alignIO = AlignIO.read('alignement/{0}'.format(self.id), 'fasta')

    def create_nj_tree(self):
        with open('nj_tree/{0}'.format(self.id), 'w') as writer:
            writer.write("")
        calculator = DistanceCalculator('identity')
        dist_matrix = calculator.get_distance(self.align_with_alignIO)
        log.info("dist_matrix")
        log.debug(dist_matrix)
        constructor = DistanceTreeConstructor()
        self.nj_tree = constructor.nj(dist_matrix)
        Phylo.write(self.nj_tree, 'nj_tree/{0}'.format(self.id), "newick", )
        self.to_json_file()

    def create_ml_tree(self):
        with open('ml_tree/{0}'.format(self.id), 'w') as writer:
            writer.write("")
        data = dendropy.DnaCharacterMatrix.get(
            path='alignement/{0}'.format(self.id),
            schema="fasta")
        log.info("data")
        log.info(data)
        rx = raxml.RaxmlRunner()
        ml_tree = rx.estimate_tree(
            char_matrix=data,
            raxml_args=["--no-bfgs"])
        self.ml_tree = ml_tree.as_string(schema="newick")
        log.info(self.ml_tree)
        self.to_json_file()
        with open('ml_tree/{0}'.format(self.id), 'w') as writer:
            writer.writelines(self.ml_tree)

    def to_json_file(self):
        self.ncbi_all_id = None
        with open("request/{0}".format(self.id), 'w') as json_request:
            json.dump(self, json_request, cls=customEncoder)
