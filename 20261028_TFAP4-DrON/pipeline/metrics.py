"""Coordinate metrics independent of ChimeraX; distances in angstroms."""
import numpy as np


def distances(x, y):
    return np.linalg.norm(np.asarray(x) - np.asarray(y), axis=-1)


def rmsd(d):
    return float(np.sqrt(np.mean(np.asarray(d) ** 2)))


def best_fit_rmsd(x, y):
    """Proper-rotation Kabsch fit, used only as a symmetry diagnostic."""
    a, b = np.asarray(x), np.asarray(y)
    a, b = a - a.mean(0), b - b.mean(0)
    u, _, vt = np.linalg.svd(a.T @ b)
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(u @ vt)
    return rmsd(distances(a @ u @ correction @ vt, b))


def ca_lddt(x, y, radius=15.0):
    """CA-only distance preservation; includes interchain pairs, excludes self.

    Not full all-atom lDDT and not AlphaFold pLDDT. Reference neighbors <15 A;
    thresholds 0.5,1,2,4 A. All 358 mapped CA sites must be present upstream.
    """
    dx = np.linalg.norm(x[:, None, :] - x[None, :, :], axis=-1)
    dy = np.linalg.norm(y[:, None, :] - y[None, :, :], axis=-1)
    mask = (dy < radius) & ~np.eye(len(y), dtype=bool)
    preservation = np.mean([np.abs(dx - dy) < t for t in (0.5, 1, 2, 4)], axis=0)
    n = mask.sum(axis=1)
    per_residue = np.divide((preservation * mask).sum(1), n,
                            out=np.full(len(y), np.nan), where=n > 0)
    return float(preservation[mask].mean()), per_residue


def summarize(x, y):
    d = distances(x, y)
    result = dict(core_ca_rmsd=rmsd(d), core_ca_mean_deviation=float(d.mean()),
                  core_ca_median_deviation=float(np.median(d)),
                  core_ca_p95_deviation=float(np.percentile(d, 95)),
                  core_ca_max_deviation=float(d.max()))
    for threshold in (1, 2, 4, 8):
        result[f"core_ca_fraction_within_{threshold}A"] = float(np.mean(d <= threshold))
    result["core_ca_lddt"], local = ca_lddt(x, y)
    return result, d, local
