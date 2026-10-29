"""Regression tests for compact core construction and MatchMaker API wiring.

ChimeraX objects are lightweight doubles here; actual numerical fitting remains
a cluster integration check, not something these tests simulate or claim.
"""
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace, ModuleType
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pipeline"))
from align_chimerax import compact_core, match_exact_core


class FakeSeq:
    def __init__(self, structure, chain_id, polymer_type):
        self.structure, self.chain_id, self.polymer_type = structure, chain_id, polymer_type

    def bulk_set(self, residues, chars):
        self.residues, self.characters = residues, chars


def structure(site=0):
    seq = "ACDEFGHIKLMNPQRSTVWY" * 18
    if site:
        seq = seq[:site] + "X" * 257 + seq[site:]
    model = SimpleNamespace(name="test", chains=[])
    for c in "AB":
        chain = FakeSeq(model, c, 1)
        chain.bulk_set([SimpleNamespace(number=n, chain_id=c) for n in range(1, len(seq) + 1)], seq)
        model.chains.append(chain)
    return model


class CoreRegressionTests(unittest.TestCase):
    def test_all_sites_keep_every_core_residue(self):
        ref = structure()
        for site in (0, 76, 80, 136):
            model = structure(site)
            rr, rc = compact_core(ref.chains[0], 0, 179, 257)
            mr, mc = compact_core(model.chains[0], site, 179, 257)
            self.assertEqual(rc, mc)
            self.assertEqual(len(mr), 179)
            self.assertEqual(mr[-1].number, 436 if site else 179)
            self.assertIs(mr[-1], model.chains[0].residues[mr[-1].number - 1])

    def call_api(self, site, truncate=False):
        atomic, backend = ModuleType("chimerax.atomic"), ModuleType("chimerax.match_maker.match")
        atomic.StructureSeq = FakeSeq

        def align(session, r, m, *args, **kwargs):
            self.assertEqual(r.characters, m.characters)
            self.assertEqual(len(r.characters), 179)
            self.assertNotIn("?", m.characters)
            self.assertNotIn("X", m.characters)
            return 100, r, m

        def match(session, pairing, pairs, matrix, alg, gap_open, gap_extend, **kwargs):
            self.assertEqual(pairing, "ss")
            self.assertEqual([(r.chain_id, m.chain_id) for r, m in pairs], [("A", "A"), ("B", "B")])
            self.assertEqual(kwargs["cutoff_distance"], 2.0)
            self.assertIs(type(kwargs["ss_fraction"]), float)
            self.assertEqual(kwargs["ss_fraction"], 0.0)
            ra, ma = [], []
            for r, m in pairs:
                _, rs, ms = kwargs["align"](session, r, m, matrix, alg, gap_open, gap_extend, {})
                n = 136 if truncate else 179
                ra += [SimpleNamespace(residue=x, name="CA") for x in rs.residues[:n]]
                ma += [SimpleNamespace(residue=x, name="CA") for x in ms.residues[:n]]
            return [{"full ref atoms": ra, "full match atoms": ma}]

        backend.align, backend.match = align, match
        model = structure(site)
        original = model.chains[0].characters
        with patch.dict(sys.modules, {"chimerax.atomic": atomic, "chimerax.match_maker.match": backend}):
            result = match_exact_core(None, structure(), model, site, 179, 257, 2.0)
        self.assertEqual(model.chains[0].characters, original)
        return result

    def test_joint_api_for_all_conditions(self):
        for site in (0, 76, 80, 136):
            self.assertEqual(len(self.call_api(site)["full match atoms"]), 358)

    def test_original_272_pair_failure_is_still_rejected(self):
        with self.assertRaisesRegex(RuntimeError, "358-site"):
            self.call_api(136, truncate=True)

    def test_missing_residue_rejected(self):
        model = structure(136)
        model.chains[0].residues[393] = None
        with self.assertRaisesRegex(ValueError, "394"):
            compact_core(model.chains[0], 136, 179, 257)

