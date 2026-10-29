import json
import sys
import tempfile
import unittest
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pipeline"))
from plot_residue_3d import summarize_residues, ORDER


class Residue3DTests(unittest.TestCase):
    def test_means_exclude_reference_and_do_not_mix_chains(self):
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            rows, models = [], []
            for ci, condition in enumerate(ORDER):
                for i in range(3):
                    mid = f"c{ci}_m{i}"
                    models.append({"model_id": mid})
                    for chain in "AB":
                        for res in (1, 2):
                            rows.append(dict(model_id=mid, condition=condition, chain=chain,
                                             tfap4_residue=res, ca_deviation=i + (10 if chain == "B" else 0),
                                             ca_plddt=70 + i * 10))
            (folder / "manifest.json").write_text(json.dumps(dict(models=models, reference_model_id="c0_m0", config={"core_end": 2})))
            (folder / "alignment_complete.json").write_text(json.dumps({"models": 12}))
            pd.DataFrame(rows).to_csv(folder / "per_residue.csv", index=False)
            summary, _, _ = summarize_residues(folder)
            a = summary[(summary.condition == ORDER[0]) & (summary.chain == "A")].iloc[0]
            b = summary[(summary.condition == ORDER[0]) & (summary.chain == "B")].iloc[0]
            self.assertEqual(a.n_models, 2)
            self.assertEqual(a.deviation_mean, 1.5)
            self.assertEqual(b.deviation_mean, 11.5)
            self.assertEqual(a.plddt_mean, 85)
            self.assertAlmostEqual(a.deviation_sd, 2 ** -0.5)
            pd.DataFrame(rows + [rows[0]]).to_csv(folder / "per_residue.csv", index=False)
            with self.assertRaisesRegex(ValueError, "Duplicate"):
                summarize_residues(folder)
