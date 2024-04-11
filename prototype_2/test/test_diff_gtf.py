import unittest
import polars as pl

from prototype_2.diff_gtf import load_gtf, hash_genes_in_gtf, compare_two_gtfs


class TestDiffGtf(unittest.TestCase):

    def test_diff_gtf(self):

        #gtf_file = "../../reference/Homo_sapiens.GRCh38.86.chr1.gtf.gz"
        gtf_file = "../../reference/Homo_sapiens.GRCh38.108.chr1.gtf.gz"
        gtf = load_gtf(gtf_file)
        self.assertIsNotNone(gtf)

        gtf_file_2 = "../../reference/Homo_sapiens.GRCh38.109.chr1.gtf.gz"
        gtf_2 = load_gtf(gtf_file_2)
        self.assertIsNotNone(gtf_2)

        gene_dict = hash_genes_in_gtf(gtf)
        self.assertIsNotNone(gene_dict)

        gene_dict_2 = hash_genes_in_gtf(gtf_2)
        self.assertIsNotNone(gene_dict_2)

        identical_df = compare_two_gtfs(gene_dict, gene_dict_2)
        self.assertIsInstance(identical_df, pl.DataFrame)

