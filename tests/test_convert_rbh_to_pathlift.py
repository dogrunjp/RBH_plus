import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).parents[1] / "src" / "convert_rbh_to_pathlift.py"
SPEC = importlib.util.spec_from_file_location("convert_rbh_to_pathlift", MODULE_PATH)
converter = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = converter
SPEC.loader.exec_module(converter)


class ConvertRbhToPathLiftTests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.fasta = self.root / "source.faa"
        self.mapping = self.root / "map_H_to_T.tsv"
        self.output = self.root / "pathlift.tsv"
        self.fasta.write_text(
            ">ENSP1.1 gene_symbol:G6PD description:first\nAAAA\n"
            ">ENSP2.2 [gene=TKT] second\nBBBB\n"
            ">ENSP3.1 no_symbol\nCCCC\n",
            encoding="utf-8",
        )
        self.mapping.write_text(
            "src_id\ttarget_ids\n"
            "ENSP1.1\tXP_1.1,XP_2.2\n"
            "ENSP2.2\tXP_3.4\n"
            "ENSP3.1\tXP_4.1\n",
            encoding="utf-8",
        )

    def tearDown(self):
        self.tempdir.cleanup()

    def read_output(self):
        with open(self.output, newline="", encoding="utf-8") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_expands_targets_and_skips_sources_without_symbols(self):
        stats = converter.convert_mapping(self.mapping, self.fasta, self.output)

        self.assertEqual(stats.mapping_rows, 3)
        self.assertEqual(stats.source_ids_with_symbol, 2)
        self.assertEqual(stats.source_ids_without_symbol, 1)
        self.assertEqual(stats.output_rows, 3)
        self.assertEqual(
            self.read_output(),
            [
                {"target_id": "XP_1.1", "source_symbol": "G6PD", "source_pid": "ENSP1.1"},
                {"target_id": "XP_2.2", "source_symbol": "G6PD", "source_pid": "ENSP1.1"},
                {"target_id": "XP_3.4", "source_symbol": "TKT", "source_pid": "ENSP2.2"},
            ],
        )

    def test_collapses_proteins_to_genes_and_matches_versionless_ids(self):
        target_map = self.root / "target_ids.tsv"
        target_map.write_text(
            "protein_id\tgene_id\n"
            "XP_1.1\t101\n"
            "XP_2.2\t101\n"
            "XP_3\t202\n",
            encoding="utf-8",
        )

        stats = converter.convert_mapping(
            self.mapping,
            self.fasta,
            self.output,
            target_id_map=target_map,
        )

        self.assertEqual(stats.output_rows, 2)
        self.assertEqual(stats.target_ids_without_mapping, 0)
        self.assertEqual(
            self.read_output(),
            [
                {"target_id": "101", "source_symbol": "G6PD", "source_pid": "ENSP1.1"},
                {"target_id": "202", "source_symbol": "TKT", "source_pid": "ENSP2.2"},
            ],
        )

    def test_missing_target_mapping_fails_without_writing_output(self):
        target_map = self.root / "target_ids.tsv"
        target_map.write_text(
            "protein_id\tgene_id\nXP_1.1\t101\n",
            encoding="utf-8",
        )

        with self.assertRaises(converter.ConversionError):
            converter.convert_mapping(
                self.mapping,
                self.fasta,
                self.output,
                target_id_map=target_map,
            )
        self.assertFalse(self.output.exists())


if __name__ == "__main__":
    unittest.main()
