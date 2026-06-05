#!/usr/bin/env python3

import pathlib
import sys
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import tplgtool2


class _FakeGroupedTplg:
    def __init__(self, widget_names, pcm_pairs):
        self.widget_list = [{"widget": {"name": name}} for name in widget_names]
        self.pcm_list = [{"pcm_id": pcm_id, "pcm_name": pcm_name} for pcm_id, pcm_name in pcm_pairs]
        self.graph_list = []


class _FakeGroupedTplgWithGraph:
    def __init__(self, widgets, edges):
        self.widget_list = [
            {"widget": {"name": name, "sname": sname}}
            for name, sname in widgets
        ]
        self.graph_list = [
            {"source": source, "sink": sink}
            for source, sink in edges
        ]


class TestTplgtool2MergeWarnings(unittest.TestCase):
    def test_no_duplicate_warnings(self):
        fake = _FakeGroupedTplg(
            widget_names=["BUF1.0", "PGA1.0"],
            pcm_pairs=[(0, "Speaker"), (1, "Mic")],
        )

        warnings = tplgtool2.get_merge_duplicate_warnings(fake)

        self.assertEqual(warnings, [])

    def test_duplicate_warnings_include_widget_pcm_id_and_name(self):
        fake = _FakeGroupedTplg(
            widget_names=["BUF1.0", "BUF1.0", "PGA1.0"],
            pcm_pairs=[(0, "Speaker"), (0, "Headset"), (2, "Speaker")],
        )

        warnings = tplgtool2.get_merge_duplicate_warnings(fake)

        self.assertEqual(len(warnings), 3)
        self.assertTrue(any("duplicate widget names found" in w and "BUF1.0" in w for w in warnings))
        self.assertTrue(any("duplicate PCM IDs found" in w and "0" in w for w in warnings))
        self.assertTrue(any("duplicate PCM names found" in w and "Speaker" in w for w in warnings))


class TestTplgtool2MergeSources(unittest.TestCase):
    def test_collect_graph_name_sources_collects_widget_sname_and_edges(self):
        fake = _FakeGroupedTplgWithGraph(
            widgets=[("BUF1.0", "BUF_A"), ("PGA1.0", "")],
            edges=[("BUF_A", "PGA1.0"), ("EXTERNAL_SRC", "BUF_A")],
        )

        sources = tplgtool2.collect_graph_name_sources(fake, "topo_a")

        self.assertEqual(sources["BUF1.0"], {"topo_a"})
        self.assertEqual(sources["BUF_A"], {"topo_a"})
        self.assertEqual(sources["PGA1.0"], {"topo_a"})
        self.assertEqual(sources["EXTERNAL_SRC"], {"topo_a"})

    def test_collect_graph_name_sources_accumulates_multiple_files(self):
        topo_a = _FakeGroupedTplgWithGraph(
            widgets=[("BUF1.0", "BUF_SHARED")],
            edges=[("BUF_SHARED", "PGA1.0")],
        )
        topo_b = _FakeGroupedTplgWithGraph(
            widgets=[("BUF2.0", "BUF_SHARED")],
            edges=[("BUF_SHARED", "PGA2.0")],
        )

        sources = tplgtool2.collect_graph_name_sources(topo_a, "topo_a")
        sources = tplgtool2.collect_graph_name_sources(topo_b, "topo_b", sources)

        self.assertEqual(sources["BUF_SHARED"], {"topo_a", "topo_b"})
        self.assertEqual(sources["BUF1.0"], {"topo_a"})
        self.assertEqual(sources["BUF2.0"], {"topo_b"})


if __name__ == "__main__":
    unittest.main()
