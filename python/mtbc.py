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

    sequence = {'NC_000962.3': {}}
    alignement = {}

    def __init__(self,
                 select_mycobacterium_canettii=False,
                 select_mycobacterium_mungi=False,
                 select_mycobacterium_orygis=False,
                 select_mycobacterium_tuberculosis=False,
                 retmax=10000000,
                 list_length=100,
                 debug=False):
        """
        @param retmax:
        @param list_length:
        @param select_mycobacterium_canettii:
        @param select_mycobacterium_mungi:
        @param select_mycobacterium_orygis:
        @param select_mycobacterium_tuberculosis:
        """

        self.align = None
        self.df_mutation = None
        self.debug = debug
        self.retmax = retmax
        self.list_length = list_length
        Entrez.email = 'A.N.Other@example.com'
        self.select_taxa = {self.taxa_mycobacterium_canettii: select_mycobacterium_canettii,
                            self.taxa_mycobacterium_mungi: select_mycobacterium_mungi,
                            self.taxa_mycobacterium_orygis: select_mycobacterium_orygis,
                            self.taxa_mycobacterium_tuberculosis: select_mycobacterium_tuberculosis}
        self.req = self.construct_search_request()
        self.id_list = self.get_all_id()
        self.acc_list = self.mtbc_random_get_accession_number()
        self.request()
        # if debug:
        #    print(json.dumps(self.sequence,indent=4))
        #    print(self.sequence.keys())
        self.reconstruct_sequence()
        self.align_reconstruct()
        self.create_nj_tree()

    def construct_search_request(self):
        if True in self.select_taxa.values():
            request = " OR ".join('txid' + key + '[ORGN]' for key, value in self.select_taxa.items() if value)
            if self.debug:
                print(request)
            return request

        else:
            request = 'txid' + self.taxa_mycobacterium_tuberculosis_complex + '[ORGN]'
            if self.debug:
                print(request)
            return request

    def get_all_id(self):
        retmax = 1000000
        handle = Entrez.esearch(db="sra",
                                term=self.req,
                                retmax=retmax,
                                retstart=0)
        record = Entrez.read(handle)
        if self.debug:
            print(record)
        handle.close()

        return record['IdList']

    def mtbc_random_get_accession_number(self):
        id_list = self.id_list
        sample_list = random.choices(id_list,
                                     k=self.list_length)

        handle_epost = Entrez.epost(db="sra",
                                    id=",".join(map(str, sample_list)),
                                    retmax=self.retmax)
        record_epost = Entrez.read(handle_epost)
        if self.debug:
            print(record_epost)
        handle_epost.close()

        handle_esummary = Entrez.esummary(db="sra",
                                          webenv=record_epost['WebEnv'],
                                          query_key=record_epost['QueryKey'],
                                          retmax=self.retmax)
        record_esummary = Entrez.read(handle_esummary)
        if self.debug:
            print(record_esummary)
        handle_esummary.close()

        runs = list(map(itemgetter('Runs'),
                        record_esummary))

        runs_1 = list(map(lambda run: "<root>" + run + "</root>",
                          runs))

        runs_acc = list(map(xmltodict.parse, runs_1))

        acc_list = pd.json_normalize(runs_acc)['root.Run.@acc']
        return acc_list.to_list()

    def request(self):
        result_alignement = {}
        for sra in self.acc_list:

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
                self.sequence[r.text[1:].split("\n")[0]] = {}
                for diff in r.text.split("\n")[2:]:
                    if len(diff) != 0:
                        ###############
                        # premiere etape, on selectionne juste les mutations (pas d'indels/insertion)
                        ###############
                        if len(diff.split(":")[3]) == 1:
                            self.sequence['NC_000962.3'][diff.split(":")[1]] = diff.split(":")[2]
                            self.sequence[r.text[1:].split("\n")[0]][diff.split(":")[1]] = diff.split(":")[3]
                key = r.text[1:].split("\n")[0]
                value = r.text[1:].split("\n")[1]
                result_alignement[key] = [len(value), value]

        return result_alignement

    def reconstruct_sequence(self):
        df_mutation = pandas.DataFrame.from_dict(self.sequence)
        df_ref = pandas.Series(df_mutation['NC_000962.3'])
        if self.debug:
            print(df_mutation)
        with open('alignement/alignement.fasta', 'w') as writer:
            for column in df_mutation.columns:
                df_mutation[column] = df_mutation[column].fillna(df_ref, axis=0)
                print(">" + column)
                print("".join(df_mutation[column].to_list()))
                writer.writelines(">" + column + "\n")
                writer.writelines("".join(df_mutation[column].to_list()) + "\n")

    def align_reconstruct(self):
        self.align = AlignIO.read('alignement/alignement.fasta', 'fasta')
        if self.debug:
            print(self.align)

    def create_nj_tree(self):
        calculator = DistanceCalculator('identity')
        dist_matrix = calculator.get_distance(self.align)
        constructor = DistanceTreeConstructor()
        nj_tree = constructor.nj(dist_matrix)
        Phylo.write(nj_tree, "tree/tree1.nwk", "newick")


if __name__ == "__main__":
    mtbc = MtbcRandom(debug=True)
    # print(mtbc.construct_search_request())
    # print()
    # print(mtbc.id_list[:10])
    # print()
    # print(mtbc.acc_list)
    # print()
    # print(json.dumps(mtbc.alignement, indent=4))
