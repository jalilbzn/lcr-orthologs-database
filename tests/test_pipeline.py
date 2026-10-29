import csv
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from pipeline.anchors import AnchorSettings, interval_between_anchors, rebuild, sw_hit
from pipeline.data import (Family, ONE_BASED, ZERO_BASED, Protein, Region, load_mappings,
                           load_proteins, load_regions, read_crosswalk)
from pipeline.export import fasta_text, filter_groups, metadata_tables
from run_pipeline import check_destination, parser, run


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.sequence = "ACDEFGHIKLMNPQRSTVWY" * 5
        self.proteins = {
            1: Protein(1, "P12345", "", "", "Homo sapiens", self.sequence),
            2: Protein(2, "", "MOUSE00001", "10090", "Mouse", self.sequence),
        }
        self.source = Region(10, 1, 30, 40, self.sequence[30:40], "strict", ZERO_BASED)
        self.target = Region(20, 2, 31, 40, self.sequence[30:40], "anchor", ONE_BASED)
        self.regions = {10: self.source, 20: self.target}
        self.families = {7: Family(7, "P12345", {"Mouse": 2})}

    def csv(self, name, fields, values):
        path = self.root / name
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(fields)
            writer.writerows(values)
        return path

    def test_coordinate_conventions_refer_to_same_residues(self):
        self.assertEqual(self.source.bounds(), self.target.bounds())

    def test_two_anchors_bound_interval_without_off_by_one(self):
        with patch("pipeline.anchors.sw_hit", side_effect=[(0, 30, 1), (40, 100, 1)]):
            interval = interval_between_anchors(self.source, self.sequence, self.sequence,
                                                None, AnchorSettings())
        self.assertEqual(interval, (30, 40))

    def test_single_left_anchor_preserves_legacy_fallback_span(self):
        with patch("pipeline.anchors.sw_hit", side_effect=[(0, 30, 1), None]):
            interval = interval_between_anchors(self.source, self.sequence, self.sequence,
                                                None, AnchorSettings())
        self.assertEqual(interval, (30, 51))

    def test_single_right_anchor_clamps_at_protein_start(self):
        with patch("pipeline.anchors.sw_hit", side_effect=[None, (5, 60, 1)]):
            interval = interval_between_anchors(self.source, self.sequence, self.sequence,
                                                None, AnchorSettings())
        self.assertEqual(interval, (0, 5))

    def test_crossed_or_empty_interval_is_rejected(self):
        for end in (39, 40):
            with self.subTest(end=end), patch("pipeline.anchors.sw_hit", side_effect=[(0, 40, 1), (end, 80, 1)]):
                self.assertIsNone(interval_between_anchors(self.source, self.sequence, self.sequence,
                                                           None, AnchorSettings()))

    def test_insufficient_coverage_is_rejected(self):
        with patch("pipeline.anchors.sw_hit", side_effect=[(0, 30, .49), (40, 90, .49)]):
            self.assertIsNone(interval_between_anchors(self.source, self.sequence, self.sequence,
                                                       None, AnchorSettings()))

    def test_interval_length_limit_is_enforced(self):
        with patch("pipeline.anchors.sw_hit", side_effect=[(0, 30, 1), (50, 90, 1)]):
            self.assertIsNone(interval_between_anchors(self.source, self.sequence, self.sequence,
                                                       None, AnchorSettings(max_interval_length=10)))

    def test_no_alignment_returns_none(self):
        class EmptyAligner:
            def align(self, query, target):
                return []
        self.assertIsNone(sw_hit(self.sequence, self.sequence, EmptyAligner()))

    def test_unsupported_matrix_residue_does_not_abort_rebuild(self):
        self.assertIsNone(sw_hit(self.sequence, self.sequence + "U", None))

    def test_source_is_first_and_fasta_is_not_gapped(self):
        text = fasta_text((7, 10), {20}, self.regions, self.proteins)
        self.assertTrue(text.startswith(">fam=7|src=10|lcr=10|pid=1|species=Homo sapiens|start=30|end=40|seg=strict\n"))
        self.assertEqual(text.count(">"), 2)
        self.assertEqual(text.splitlines()[1], self.sequence[30:40])

    def test_length_filter_keeps_five_aa_and_removes_four_aa_ortholog(self):
        regions = dict(self.regions)
        regions[20] = Region(20, 2, 1, 5, self.sequence[:5], "anchor", ONE_BASED)
        regions[21] = Region(21, 2, 1, 4, self.sequence[:4], "anchor", ONE_BASED)
        original = {(7, 10): {20, 21}}
        groups, counts = filter_groups(original, regions, 5)
        self.assertEqual(groups, {(7, 10): {20}})
        self.assertEqual(original, {(7, 10): {20, 21}})
        self.assertEqual(counts["short_ortholog_entries"], 1)
        self.assertEqual(counts["removed_sequence_entries"], 1)

    def test_short_source_removes_entire_group(self):
        regions = dict(self.regions)
        regions[10] = Region(10, 1, 0, 4, self.sequence[:4], "strict", ZERO_BASED)
        groups, counts = filter_groups({(7, 10): {20}}, regions, 5)
        self.assertEqual(groups, {})
        self.assertEqual(counts["short_source_groups"], 1)
        self.assertEqual(counts["removed_sequence_entries"], 2)

    def test_length_filter_removes_group_without_remaining_orthologs(self):
        regions = dict(self.regions)
        regions[20] = Region(20, 2, 1, 4, self.sequence[:4], "anchor", ONE_BASED)
        groups, counts = filter_groups({(7, 10): {20}}, regions, 5)
        self.assertEqual(groups, {})
        self.assertEqual(counts["groups_without_orthologs"], 1)
        self.assertEqual(counts["removed_sequence_entries"], 2)

    def test_one_aa_threshold_preserves_short_entries(self):
        regions = dict(self.regions)
        regions[20] = Region(20, 2, 1, 1, self.sequence[:1], "anchor", ONE_BASED)
        original = {(7, 10): {20}}
        groups, counts = filter_groups(original, regions, 1)
        self.assertEqual(groups, original)
        self.assertEqual(counts["removed_sequence_entries"], 0)

    def test_invalid_length_threshold_is_rejected(self):
        for threshold in (0, -1):
            with self.subTest(threshold=threshold), self.assertRaisesRegex(ValueError, "at least 1"):
                filter_groups({(7, 10): {20}}, self.regions, threshold)

    def test_filtered_export_and_metadata_include_only_retained_members(self):
        proteins = dict(self.proteins)
        proteins[3] = Protein(3, "", "CHICK00001", "", "Chicken", self.sequence)
        regions = dict(self.regions)
        regions[30] = Region(30, 3, 1, 4, self.sequence[:4], "anchor", ONE_BASED)
        args = parser().parse_args(["--inputs", str(self.root / "inputs"),
                                    "--output", str(self.root / "outputs")])
        captured = {}
        def capture(path, tables):
            captured.update(tables)
            path.write_bytes(b"test workbook")
        with patch("run_pipeline.load_inputs", return_value=(proteins, regions, self.families)), \
             patch("run_pipeline.load_mappings", return_value={(7, 10): {20, 30}}):
            counts = run(args, metadata_writer=capture)
        default = args.output / "anchor_fastas_default" / "7_10.fasta"
        strict = args.output / "anchor_fastas_strict" / "7_10.fasta"
        self.assertEqual(default.read_bytes(), strict.read_bytes())
        self.assertEqual(default.read_text().count(">"), 2)
        self.assertNotIn("lcr=30|", default.read_text())
        self.assertEqual(captured["Groups"]["rows"][0][4:6], [2, 2])
        self.assertEqual({row[0] for row in captured["Proteins"]["rows"]}, {"1", "2"})
        self.assertIn("Minimum LCR length: 5 aa", captured["Groups"]["notes"][0])
        self.assertEqual(counts["length_filter_short_ortholog_entries"], 1)

    def test_filter_removing_all_groups_preserves_existing_output(self):
        output = self.root / "outputs"
        output.mkdir()
        previous = output / "metadata.xlsx"
        previous.write_bytes(b"previous workbook")
        args = parser().parse_args(["--inputs", str(self.root / "inputs"),
                                    "--output", str(output), "--overwrite", "--min-lcr-length", "11"])
        with patch("run_pipeline.load_inputs", return_value=(self.proteins, self.regions, self.families)), \
             patch("run_pipeline.load_mappings", return_value={(7, 10): {20}}):
            with self.assertRaisesRegex(ValueError, "No anchored groups remain"):
                run(args)
        self.assertEqual(previous.read_bytes(), b"previous workbook")

    def test_oma_ids_are_not_written_as_uniprot_accessions(self):
        path = self.csv("protein.csv", ["protein_id", "uniprot_acc", "species", "sequence"],
                        [[2, "MOUSE00001", "Mouse", self.sequence]])
        protein = load_proteins(path)[2]
        self.assertEqual(protein.uniprot_acc, "")
        self.assertEqual(protein.oma_id, "MOUSE00001")

    def test_bad_region_coordinate_is_rejected(self):
        path = self.csv("lcr_seq.csv", ["lcr_id", "protein_id", "start_pos", "end_pos", "sequence", "seg_mode", "coordinate_system"],
                        [[10, 1, 30, 40, self.sequence[30:40], "strict", ONE_BASED]])
        with self.assertRaisesRegex(ValueError, "coordinates do not match"):
            load_regions(path, self.proteins)

    def test_duplicate_protein_id_is_rejected(self):
        path = self.csv("protein.csv", ["protein_id", "uniprot_acc", "species", "sequence"],
                        [[1, "P12345", "Homo sapiens", self.sequence]] * 2)
        with self.assertRaisesRegex(ValueError, "duplicate ID"):
            load_proteins(path)

    def test_orphan_mapping_is_rejected(self):
        path = self.csv("phylogeny_map.csv", ["family_id", "source_lcr_id", "lcr_id"], [[7, 10, 999]])
        with self.assertRaisesRegex(ValueError, "missing referenced record"):
            load_mappings(path, self.proteins, self.regions, self.families)

    def test_wrong_family_member_is_rejected(self):
        path = self.csv("phylogeny_map.csv", ["family_id", "source_lcr_id", "lcr_id"], [[7, 10, 20]])
        with self.assertRaisesRegex(ValueError, "not a family ortholog"):
            load_mappings(path, self.proteins, self.regions, {7: Family(7, "P12345", {})})

    def test_crosswalk_preserves_multiple_matches(self):
        path = self.root / "oma-uniprot.tsv"
        path.write_text("# mapping\nMOUSE00001\tP23456\nMOUSE00001\tQ23456\nHUMAN00001\tP12345\n")
        mappings = read_crosswalk(path, self.proteins)
        tables = metadata_tables({(7, 10): {20}}, self.regions, self.proteins, mappings)
        self.assertEqual(tables["Proteins"]["rows"][1][1], "P23456;Q23456")
        self.assertEqual(tables["Proteins"]["rows"][0][2], "HUMAN00001")

    def test_missing_identifiers_remain_unavailable(self):
        table = metadata_tables({(7, 10): {20}}, self.regions, self.proteins)
        human, mouse = table["Proteins"]["rows"]
        self.assertEqual(human[2], "")
        self.assertIsNone(human[4])
        self.assertEqual(mouse[1], "")
        self.assertEqual(table["Groups"]["rows"][0][4:6], [2, 2])

    def test_rebuild_discards_old_mappings_and_uses_only_human_seeds(self):
        with patch("pipeline.anchors.interval_between_anchors", return_value=(30, 40)):
            regions, groups, counts = rebuild(self.proteins, self.regions, self.families,
                                               AnchorSettings(), aligners={"BLOSUM62": None})
        self.assertEqual(counts, {"attempted": 1, "mapped": 1})
        self.assertNotIn(20, regions)
        self.assertEqual(len(groups[(7, 10)]), 1)
        mapped = regions[next(iter(groups[(7, 10)]))]
        self.assertEqual((mapped.start, mapped.end, mapped.sequence), (31, 40, self.sequence[30:40]))

    def test_input_directory_is_protected(self):
        with self.assertRaisesRegex(ValueError, "must not contain"):
            check_destination(self.root / "inputs", self.root, True)

    def test_unrelated_output_is_protected(self):
        output = self.root / "output"
        output.mkdir()
        (output / "notes.txt").write_text("keep me")
        with self.assertRaisesRegex(ValueError, "unrelated files"):
            check_destination(self.root / "inputs", output, True)

    def test_failed_metadata_write_preserves_previous_successful_output(self):
        inputs, output = self.root / "inputs", self.root / "outputs"
        inputs.mkdir()
        output.mkdir()
        previous = output / "metadata.xlsx"
        previous.write_bytes(b"previous workbook")
        args = parser().parse_args(["--inputs", str(inputs), "--output", str(output), "--overwrite"])
        def fail(path, tables):
            raise RuntimeError("test export failure")
        with patch("run_pipeline.load_inputs", return_value=(self.proteins, self.regions, self.families)), \
             patch("run_pipeline.load_mappings", return_value={(7, 10): {20}}):
            with self.assertRaisesRegex(RuntimeError, "test export failure"):
                run(args, metadata_writer=fail)
        self.assertEqual(previous.read_bytes(), b"previous workbook")


if __name__ == "__main__":
    unittest.main()
