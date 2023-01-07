import io

import pandas
import requests
from pandas import CategoricalDtype
import logging
from .mtbc_ncbi import MtbcGetRandomSRA

FORMAT = "%(levelname)s:%(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)


class MtbcAcclistToFASTA:

    def __init__(self,
                 mtbc_get_random_sra: MtbcGetRandomSRA,
                 sequence_dict,
                 target_list_length,
                 final_acc_list,
                 snp_select=[],
                 snp_reject=[]):

        # initial user variable
        self.snp_select = snp_select
        self.snp_reject = snp_reject
        self.final_acc_list_length = None
        self.target_list_length: int = int(target_list_length)
        self.fasta = None
        self.ncbi_random_acc_list = mtbc_get_random_sra.ncbi_random_acc_list
        self.df_mutation = None
        self._id = mtbc_get_random_sra._id
        self.sequence_dict = sequence_dict
        self.final_acc_list: list = final_acc_list
        logging.info("len(self.sequence_dict)")
        logging.info(len(self.sequence_dict))

        if len(self.sequence_dict) < 2:
            self.mtbc_request()
        logging.info("self.sequence")
        logging.debug(self.sequence_dict)
        self.reconstruct_sequence_to_fasta_file()

    def mtbc_request(self):
        logging.info("mtbc_request")
        if self.final_acc_list:
            logging.info("if self.final_acc_list:")
            sra_list = self.final_acc_list
        else:
            logging.info("if self.final_acc_list: else :")
            sra_list = self.ncbi_random_acc_list
        index = 1
        reject = 0
        logging.info("mtbc_request")
        self.final_acc_list_length = 0

        ## check if snp reject
        check_snp_reject = False
        check_snp_select = False

        reject = 0
        not_select = 0

        if len(self.snp_reject) > 0:
            check_snp_reject = True

        if len(self.snp_select) > 0:
            check_snp_select = True

        for sra in sra_list:
            if self.final_acc_list_length >= self.target_list_length:
                return
            logging.info("ncbi")
            logging.info(str(index) + "/" + str(len(self.ncbi_random_acc_list)))
            logging.info("final")
            logging.info(str(self.final_acc_list_length) + "/" + str(self.target_list_length))

            head = {'Content-Type': 'application/x-www-form-urlencoded',
                    'Host': 'gnksoftware.synology.me:30002'}
            parameters = {
                'search': '{"type":"boolean","quantifier":"allOf","filters":[{"type":"terms","terms":["' +
                          str(sra) +
                          '"],"attribute":"publicId.keyword","quantifier":"anyOf"}]}',
                'alignmentParams': '{"types":["snp"],"includeDataInComment":true}'}
            url = 'http://gnksoftware.synology.me:30002/strains/alignment'
            r = requests.post(url, parameters, head)

            ## reject
            logging.info(("list of reject snp : " + str(self.snp_reject)))
            if check_snp_reject and any(snp in r.text for snp in self.snp_reject):
                logging.info("#########################################")
                logging.info("SRA : " + str(sra) + " is reject because contains snp_reject ")
                logging.info("Nombre de SNP_reject : " + str(reject))
                logging.info("#########################################")
                reject += 1
                continue

            logging.info(("list of select snp : " + str(self.snp_select)))
            if check_snp_select and not all(snp in r.text for snp in self.snp_select):
                logging.info("#########################################")
                logging.info("SRA : " + str(sra) + " is reject because not contains snp_select ")
                logging.info("Nombre de SNP_not_select : " + str(not_select))
                logging.info("#########################################")
                not_select += 1
                continue

            index += 1
            if len(r.text) != 0:
                logging.info(r.text[1:].split("\n")[0])
                self.final_acc_list.append(r.text[1:].split("\n")[0])
                self.sequence_dict[r.text[1:].split("\n")[0]] = {}
                self.final_acc_list_length += 1
                logging.debug(self.final_acc_list)
                for diff in r.text.split("\n")[2:]:
                    if len(diff) != 0:
                        ###############
                        # premiere etape, on selectionne juste les mutations (pas d'indels/insertion)
                        ###############
                        if len(diff.split(":")[3]) == 1:
                            self.sequence_dict['NC_000962.3'][diff.split(":")[1]] = diff.split(":")[2]
                            self.sequence_dict[r.text[1:].split("\n")[0]][diff.split(":")[1]] = diff.split(":")[3]



        logging.info("#########################################")
        logging.info("#########################################")
        logging.info("#########################################")
        logging.info("#########################################")
        logging.info("Nombre de SNP_reject : " + str(reject))
        logging.info("#########################################")
        logging.info("#########################################")
        logging.info("Nombre de SNP_not_select : " + str(not_select))
        logging.info("#########################################")
        logging.info("#########################################")
        logging.info("#########################################")
        logging.info("#########################################")

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
