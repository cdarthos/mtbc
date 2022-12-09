import random
from operator import itemgetter

import pandas
import pandas as pd
import requests
import xmltodict
from Bio import Entrez
from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor


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

    sequence_dict = {'NC_000962.3': {}}
    alignement = {}

    def __init__(self,
                 select_mycobacterium_canettii=False,
                 select_mycobacterium_mungi=False,
                 select_mycobacterium_orygis=False,
                 select_mycobacterium_tuberculosis=False,
                 retmax=10000000,
                 list_length=10,
                 debug=False,
                 email = 'A.N.Other@example.com'):

        # initial user variable
        self.debug = debug
        self.retmax = retmax
        self.list_length = list_length
        Entrez.email = email
        self.select_taxa = {self.taxa_mycobacterium_canettii: select_mycobacterium_canettii,
                            self.taxa_mycobacterium_mungi: select_mycobacterium_mungi,
                            self.taxa_mycobacterium_orygis: select_mycobacterium_orygis,
                            self.taxa_mycobacterium_tuberculosis: select_mycobacterium_tuberculosis}

        # initialize empty variable
        self.acc_list = None
        self.sample_list = None
        self.ncbi_all_id = None
        self.ncbi_request_all_id = None
        self.align_with_alignIO = None
        self.df_mutation = None

        # main program
        self.construct_search_request()
        if self.debug:
            print("self.ncbi_request_all_id")
            print(self.ncbi_request_all_id)
        self.get_all_id()
        if self.debug:
            print("self.ncbi_all_id")
            print(self.ncbi_all_id)
        self.select_random_ncbi_id_number()
        if self.debug:
            print("self.sample_list")
            print(self.sample_list)
        self.get_acc_number_from_ncbi_id()
        if self.debug:
            print("self.acc_list")
            print(self.acc_list)
        self.mtbc_request()
        if self.debug:
            print("self.sequence")
            print(self.sequence_dict)
        self.reconstruct_sequence_to_fasta_file()
        self.align_reconstruct()
        if self.debug:
            print("self.align_with_alignIO")
            print(self.align_with_alignIO)
        self.create_nj_tree()

    def construct_search_request(self):
        if True in self.select_taxa.values():
            self.ncbi_request_all_id = " OR ".join(
                'txid' + key + '[ORGN]' for key, value in self.select_taxa.items() if value)

        else:
            self.ncbi_request_all_id = 'txid' + self.taxa_mycobacterium_tuberculosis_complex + '[ORGN]'

    def get_all_id(self):
        retmax = 1000000
        handle = Entrez.esearch(db="sra",
                                term=self.ncbi_request_all_id,
                                retmax=retmax,
                                retstart=0)
        record = Entrez.read(handle)
        if self.debug:
            print("Entrez.esearch")
            print(record)
        handle.close()
        self.ncbi_all_id = record['IdList']

    def select_random_ncbi_id_number(self):
        self.sample_list = random.choices(self.ncbi_all_id,
                                          k=self.list_length)

    def get_acc_number_from_ncbi_id(self):
        handle_esummary = Entrez.esummary(db="sra",
                                          id=",".join(map(str, self.sample_list))
                                          )
        record_esummary = Entrez.read(handle_esummary)
        if self.debug:
            print("record_esummary")
            print(record_esummary)
        handle_esummary.close()

        runs = list(map(itemgetter('Runs'),
                        record_esummary))

        runs_1 = list(map(lambda run: "<root>" + run + "</root>",
                          runs))

        runs_acc = list(map(xmltodict.parse, runs_1))

        acc_list = pd.json_normalize(runs_acc)['root.Run.@acc']
        self.acc_list = acc_list.to_list()

    def mtbc_request(self):
        acc_list_len = len(self.acc_list)
        index = 1
        if self.debug:
            print("mtbc_request")
        for sra in self.acc_list:
            if self.debug:
                print(str(index) + "/" + str(acc_list_len))
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

        df_mutation = pandas.DataFrame.from_dict(self.sequence_dict, dtype='category')

        df_ref = pandas.Series(df_mutation['NC_000962.3'])

        if self.debug:
            print(df_mutation.info(memory_usage="deep"))
            print(df_mutation.astype('category').info(memory_usage="deep"))

        with open('alignement/reconstruct_sequence.fasta', 'w') as writer:
            for column in df_mutation.columns:
                df_mutation[column] = df_mutation[column].fillna(df_ref, axis=0)
                writer.writelines(">" + column + "\n")
                writer.writelines("".join(df_mutation[column].to_list()) + "\n")

    def align_reconstruct(self):
        self.align_with_alignIO = AlignIO.read('alignement/reconstruct_sequence.fasta', 'fasta')

    def create_nj_tree(self):
        calculator = DistanceCalculator('identity')
        dist_matrix = calculator.get_distance(self.align_with_alignIO)
        if self.debug:
            print("dist_matrix")
            print(dist_matrix)
        constructor = DistanceTreeConstructor()
        nj_tree = constructor.nj(dist_matrix)
        Phylo.write(nj_tree, "tree/tree1.nwk", "newick")
        if self.debug:
            print("Phylo.draw_ascii(nj_tree)")
            Phylo.draw_ascii(nj_tree)


if __name__ == "__main__":
    mtbc = MtbcRandom(debug=True, list_length=1000)
    # print(mtbc.construct_search_request())
    # print()
    # print(mtbc.id_list[:10])
    # print()
    # print(mtbc.acc_list)
    # print()
    # print(json.dumps(mtbc.alignement, indent=4))
