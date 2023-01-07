import io
import logging

import dendropy
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from dendropy.interop import raxml


class MtbcTree:

    @staticmethod
    def create_nj_tree_static(fasta):
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
        result = handle.getvalue()
        handle.close()
        return result

    @staticmethod
    def create_ml_tree_static(fasta):
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
