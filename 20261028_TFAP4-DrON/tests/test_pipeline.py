"""Run: python -m unittest discover -s tests -v"""
import json
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pipeline"))
from common import construct_number, core_spec, classify
from metrics import best_fit_rmsd, summarize
from prepare import prepare


class NumberingTests(unittest.TestCase):
    def test_all_boundaries(self):
        for site, spec in [(76, "1-76,334-436"), (80, "1-80,338-436"), (136, "1-136,394-436")]:
            self.assertEqual(core_spec(site), spec)
            numbers = [construct_number(i, site) for i in range(1, 180)]
            self.assertEqual(len(set(numbers)), 179)
            self.assertEqual(numbers[site - 1], site)
            self.assertEqual(numbers[site], site + 258)
            self.assertEqual(numbers[-1], 436)
            self.assertFalse(any(site < n <= site + 257 for n in numbers))

    def test_sequence_validation(self):
        control = "ACDEFGHIKLMNPQRSTVWY" * 18
        for site in (76, 80, 136):
            sequence = control[:site] + "Z" * 257 + control[site:]
            self.assertEqual(classify(sequence, control, [76, 80, 136], 257), site)
        with self.assertRaises(ValueError):
            classify(control[:76] + "Z" * 256 + control[76:], control, [76, 80, 136], 257)


class MetricTests(unittest.TestCase):
    def setUp(self):
        self.x = np.random.default_rng(100).normal(size=(30, 3)) * 3

    def test_identity_and_rigid_transform(self):
        metrics, d, local = summarize(self.x, self.x)
        self.assertEqual(metrics["core_ca_rmsd"], 0)
        self.assertEqual(metrics["core_ca_lddt"], 1)
        np.testing.assert_allclose(local, 1)
        rotated = self.x @ np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]]) + [50, 20, -10]
        self.assertLess(best_fit_rmsd(rotated, self.x), 1e-12)
        self.assertAlmostEqual(summarize(rotated, self.x)[0]["core_ca_lddt"], 1)

    def test_distortion_and_reflection(self):
        moved = self.x.copy()
        moved[0] += 30
        self.assertLess(summarize(moved, self.x)[0]["core_ca_lddt"], 1)
        reflected = self.x * [-1, 1, 1]
        self.assertGreater(best_fit_rmsd(reflected, self.x), 0.5)


class ManifestTests(unittest.TestCase):
    def test_complete_manifest_and_protect_existing_output(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            cfg = dict(project_root=str(root), input_directory="inputs", output_directory="results",
                       insertion_sites=[76, 80, 136], insertion_length=257, core_end=179,
                       expected_seeds_per_condition=5, expected_models_per_seed=5)
            control = "ACDEFGHIKLMNPQRSTVWY" * 18
            for site in (0, 76, 80, 136):
                sequence = control if not site else control[:site] + "Z" * 257 + control[site:]
                for seed in range(5):
                    path = root / "inputs" / f"site{site}_seed{seed}"
                    path.mkdir(parents=True)
                    stem = "fold_example"
                    request = [dict(name=path.name, modelSeeds=[seed], sequences=[
                        {"proteinChain": dict(sequence=sequence, count=2)},
                        {"dnaSequence": dict(sequence="ATGC", count=2)}])]
                    (path / f"{stem}_job_request.json").write_text(json.dumps(request))
                    for rank in range(5):
                        for suffix in (f"model_{rank}.cif", f"summary_confidences_{rank}.json", f"full_data_{rank}.json"):
                            (path / f"{stem}_{suffix}").write_text("{}")
            config = root / "config.json"
            config.write_text(json.dumps(cfg))
            manifest = prepare(config)
            self.assertEqual(len(manifest["models"]), 100)
            self.assertEqual(manifest["reference_model_id"], "site0_seed0_rank0")
            self.assertEqual(set(manifest["counts"].values()), {25})
            with self.assertRaises(FileExistsError):
                prepare(config)
            # Missing poses must fail before creating a new output directory.
            cfg["output_directory"] = "other"
            config.write_text(json.dumps(cfg))
            next((root / "inputs").rglob("*_model_4.cif")).unlink()
            with self.assertRaises(ValueError):
                prepare(config)
            self.assertFalse((root / "other").exists())


if __name__ == "__main__":
    unittest.main()
