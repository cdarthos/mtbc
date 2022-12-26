import pandas
import requests
from Bio import Entrez, AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from pandas import CategoricalDtype


class MtbcAcclistToFASTA:

    def __init__(self,
                 mtbc_get_random_sra):

        # initial user variable
        self.debug = mtbc_get_random_sra.debug
        self.retmax = mtbc_get_random_sra.retmax
        self.list_length = mtbc_get_random_sra.list_length
        Entrez.email = mtbc_get_random_sra.email
        self.select_taxa = mtbc_get_random_sra.select_taxa

        # initialize empty variable
        self.nj_tree = None
        self.acc_list = mtbc_get_random_sra.ncbi_random_acc_list
        self.sample_list = None
        self.ncbi_all_id = None
        self.ncbi_request_all_id = None
        self.align_with_alignIO = None
        self.df_mutation = None
        self.id = mtbc_get_random_sra.id
        self.sequence_dict = {'NC_000962.3': {}}
        self.alignement = {}

        self.mtbc_request()
        if self.debug:
            print("self.sequence")
            # print(self.sequence_dict)
        self.reconstruct_sequence_to_fasta_file()
        # self.align_reconstruct()
        # if self.debug:
        #     print("self.align_with_alignIO")
        #     #print(self.align_with_alignIO)
        # self.create_nj_tree()

    def mtbc_request(self):
        acc_list_len = len(self.acc_list)
        index = 1
        if True or self.debug:
            print("mtbc_request")
        for sra in self.acc_list:
            if True or self.debug:
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
                if self.debug:
                    print(r.text[1:].split("\n")[0])
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

        if self.debug:
            print(df_mutation.info(memory_usage="deep"))
        with open('alignement/{0}'.format(self.id), 'w') as writer:
            for column in df_mutation.columns:
                df_mutation[column] = df_mutation[column].fillna(df_mutation['NC_000962.3'], axis=0)
                writer.writelines(">" + column + "\n")
                writer.writelines("".join(df_mutation[column].to_list()) + "\n")

    def align_reconstruct(self):
        self.align_with_alignIO = AlignIO.read('alignement/{0}'.format(self.id), 'fasta')

    def create_nj_tree(self):
        calculator = DistanceCalculator('identity')
        dist_matrix = calculator.get_distance(self.align_with_alignIO)
        if self.debug:
            print("dist_matrix")
            print(dist_matrix)
        constructor = DistanceTreeConstructor()
        nj_tree = constructor.nj(dist_matrix)
        Phylo.write(nj_tree, 'nj_tree/{0}'.format(self.id), "newick")
        if self.debug:
            print("Phylo.draw_ascii(nj_tree)")
            Phylo.draw_ascii(nj_tree)
