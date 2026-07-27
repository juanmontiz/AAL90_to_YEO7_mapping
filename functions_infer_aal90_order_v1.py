"""
AAL90 parcel-ordering inference.

Discriminates between two conventions for a 90-region AAL matrix / time series:

    'non-symmetric LLRR' : L1..L45, R1..R45     (contiguous hemispheric blocks)
    'symmetric LRLR'     : L1, R1, L2, R2, ...  (interleaved; SPM/AAL native)

NOTE ON TERMINOLOGY
-------------------
'symmetric' / 'non-symmetric' here describe the PARCEL LAYOUT, not matrix
symmetry (SC == SC.T). The latter is reported separately as 'asymmetry'.

Design note
-----------
The discriminating statistic is BLOCK DENSITY, not a homotopic diagonal.
Rationale: AAL indexes regions in spatially contiguous order, so parcels
(2k, 2k+1) are frequently anatomical neighbours within the same hemisphere.
In tractography-derived SC, short-range connectivity dominates, so the
(2k, 2k+1) diagonal is strongly enhanced even in a true LLRR matrix -- a
diagonal-based statistic is therefore maximally misled precisely when the
corpus callosum is under-recovered (crossing-fibre problem in deterministic
DWI). Block-level aggregates cannot be spoofed this way: under a wrong
ordering every 45x45 quadrant mixes hemispheres and the intra/inter density
contrast collapses.
"""

import numpy as np
from scipy.stats import spearmanr


# ─────────────────────────────────────────────────────────────────────────────
# Atlas constants
# ─────────────────────────────────────────────────────────────────────────────
N_AAL90 = 90
N_HEMI = N_AAL90 // 2  # 45 regions per hemisphere

# ─────────────────────────────────────────────────────────────────────────────
# Ordering labels (compare against these constants, never against literals)
# ─────────────────────────────────────────────────────────────────────────────
ORDER_LLRR = "non-symmetric LLRR"
ORDER_LRLR = "symmetric LRLR"
ORDER_UNDET = "undetermined"


# ─────────────────────────────────────────────────────────────────────────────
# Permutations
# ─────────────────────────────────────────────────────────────────────────────
def _deinterleave_perm(n: int = N_AAL90) -> np.ndarray:
    """LRLR -> LLRR. p[i] = index into the original (interleaved) matrix."""
    return np.r_[np.arange(0, n, 2), np.arange(1, n, 2)]


def _interleave_perm(n: int = N_AAL90) -> np.ndarray:
    """LLRR -> LRLR (inverse of _deinterleave_perm)."""
    p = np.empty(n, dtype=int)
    p[0::2] = np.arange(0, n // 2)
    p[1::2] = np.arange(n // 2, n)
    return p


# ─────────────────────────────────────────────────────────────────────────────
# Component statistics (all assume M is in canonical LLRR order)
# ─────────────────────────────────────────────────────────────────────────────
def _block_density_ratio(M: np.ndarray, n_hemi: int = N_HEMI) -> float:
    """
    log10 ratio of intra- to inter-hemispheric mean density.

    True LLRR   -> large positive (~0.9, i.e. ~8x, in typical DWI SC).
    Wrong order -> collapses toward 0 as quadrants mix hemispheres.

    Complexity O(n^2). Immune to index-adjacency confounds.
    """
    intra = np.nanmean(np.r_[M[:n_hemi, :n_hemi].ravel(),
                             M[n_hemi:, n_hemi:].ravel()])
    inter = np.nanmean(M[:n_hemi, n_hemi:])
    if not (np.isfinite(intra) and np.isfinite(inter) and intra > 0 and inter > 0):
        return 0.0
    return float(np.log10(intra / inter))


def _block_balance(M: np.ndarray, n_hemi: int = N_HEMI) -> float:
    """
    |log10(LL / RR)| density imbalance. ~0 when the split at n_hemi is a true
    hemisphere boundary (the hemispheres are near-mirror images); penalty term.
    """
    ll = np.nanmean(M[:n_hemi, :n_hemi])
    rr = np.nanmean(M[n_hemi:, n_hemi:])
    if not (np.isfinite(ll) and np.isfinite(rr) and ll > 0 and rr > 0):
        return 1.0  # bounded penalty; np.inf would poison the margin
    return float(abs(np.log10(ll / rr)))


def _intrahemispheric_homology(M: np.ndarray, n_hemi: int = N_HEMI) -> float:
    """
    Spearman rank correlation between the LL and RR upper-triangular profiles.

    Rank-based, hence independent of density and safe to combine with the
    density ratio. Weighted SC is ~log-normal, so Spearman over Pearson.
    Secondary evidence only: in deterministic DWI this term is weak and can
    even favour the wrong hypothesis (mixed quadrants correlate trivially),
    which is why it is down-weighted in the composite score.
    """
    iu = np.triu_indices(n_hemi, k=1)
    x, y = M[:n_hemi, :n_hemi][iu], M[n_hemi:, n_hemi:][iu]
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3 or np.ptp(x[m]) == 0 or np.ptp(y[m]) == 0:
        return 0.0
    rho = spearmanr(x[m], y[m]).statistic
    return float(rho) if np.isfinite(rho) else 0.0


def _score(M: np.ndarray, n_hemi: int = N_HEMI, w_homology: float = 0.5) -> dict:
    """
    Composite score for the hypothesis that M is in LLRR order.

    total = log10(intra/inter) - |log10(LL/RR)| + w_homology * rho
    """
    dr = _block_density_ratio(M, n_hemi)
    bb = _block_balance(M, n_hemi)
    rho = _intrahemispheric_homology(M, n_hemi)
    return {
        "density_ratio": dr,
        "balance_penalty": bb,
        "homology_rho": rho,
        "total": dr - bb + w_homology * rho,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Main entry point
# ─────────────────────────────────────────────────────────────────────────────
def infer_aal90_order(SC, verbose: bool = True, margin_tol: float = 0.30,
                      w_homology: float = 0.5):
    """
    Infer the parcel ordering of a 90x90 AAL matrix.

    Parameters
    ----------
    SC : array_like, (90, 90)
        Structural or functional connectivity. Symmetrized internally;
        diagonal is masked out.
    verbose : bool
        Print full diagnostics.
    margin_tol : float
        Minimum |total_LLRR - total_LRLR| to declare the decision reliable.
        Scores are log10 density ratios, so 0.30 ~ a factor of 2 in contrast.
    w_homology : float
        Weight of the Spearman LL~RR term. 0.5 for SC; consider 1.0 for FC,
        where interhemispheric homology is genuinely strong (rho ~ 0.6-0.8),
        or 0.0 to rely on block density alone.

    Returns
    -------
    order : {'non-symmetric LLRR', 'symmetric LRLR', 'undetermined'}
        Use the module constants ORDER_LLRR / ORDER_LRLR / ORDER_UNDET for
        comparisons rather than string literals.
    info : dict
        'LLRR', 'LRLR'        : component-score dicts (machine-friendly keys)
        'margin'              : |total_LLRR - total_LRLR|
        'permutation_to_LLRR' : ndarray to canonicalize, or None if undetermined
        'asymmetry'           : relative asymmetry of the input (SC vs SC.T)

    Notes
    -----
    Verification only. Native AAL90 is LRLR by construction (odd index = left);
    the authoritative answer is the atlas label file / centroid MNI x-sign.
    Use this when labels are unavailable or as a cross-check. Note also that
    block density is L/R-symmetric: this cannot tell which block is left.
    """
    SC = np.asarray(SC, dtype=float)
    if SC.ndim != 2 or SC.shape != (N_AAL90, N_AAL90):
        raise ValueError(f"Expected a (90, 90) AAL90 matrix, got {SC.shape}.")

    scale = np.nanmax(np.abs(SC))
    asym = float(np.nanmax(np.abs(SC - SC.T)) / (scale + 1e-12)) if scale > 0 else 0.0

    M = 0.5 * (SC + SC.T)
    np.fill_diagonal(M, np.nan)  # self-connections carry no ordering information

    # H1: already LLRR
    s_llrr = _score(M, w_homology=w_homology)

    # H2: input is LRLR -> de-interleave, then apply the identical statistic
    perm = _deinterleave_perm()
    s_lrlr = _score(M[np.ix_(perm, perm)], w_homology=w_homology)

    margin = abs(s_llrr["total"] - s_lrlr["total"])

    if margin < margin_tol:
        order, perm_out = ORDER_UNDET, None
    elif s_llrr["total"] > s_lrlr["total"]:
        order, perm_out = ORDER_LLRR, np.arange(N_AAL90)
    else:
        order, perm_out = ORDER_LRLR, perm

    if verbose:
        print(f"Input asymmetry (SC vs SC.T) : {asym:.2e}")
        print()
        print(f"{'':26s}{'H1 LLRR':>12s}{'H2 LRLR':>12s}")
        for k, lab in [("density_ratio",   "log10(intra/inter)"),
                       ("balance_penalty", "|log10(LL/RR)|"),
                       ("homology_rho",    "Spearman LL~RR")]:
            print(f"  {lab:24s}{s_llrr[k]:>12.4f}{s_lrlr[k]:>12.4f}")
        print(f"  {'TOTAL':24s}{s_llrr['total']:>12.4f}{s_lrlr['total']:>12.4f}")
        print()
        print(f"Decision margin              : {margin:.4f}  (tol = {margin_tol})")
        print(f"Inferred ordering            : {order}")
        if order == ORDER_UNDET:
            print("  -> Scores too close. Resolve against the atlas label file.")

    return order, {
        "LLRR": s_llrr,
        "LRLR": s_lrlr,
        "margin": margin,
        "permutation_to_LLRR": perm_out,
        "asymmetry": asym,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Canonicalization helper
# ─────────────────────────────────────────────────────────────────────────────
def canonicalize_to_llrr(SC, ts=None, verbose: bool = True):
    """
    Reorder SC (and optionally a (90, T) time-series array) into LLRR order.

    Returns
    -------
    SC_canon : (90, 90) ndarray
    ts_canon : (90, T) ndarray or None
    info     : dict from infer_aal90_order

    Raises
    ------
    RuntimeError if the ordering cannot be determined.
    """
    order, info = infer_aal90_order(SC, verbose=verbose)
    p = info["permutation_to_LLRR"]
    if p is None:
        raise RuntimeError(
            "Ordering undetermined; canonicalization refused. "
            "Check the atlas label file."
        )
    SC = np.asarray(SC, dtype=float)
    SC_canon = SC[np.ix_(p, p)]
    ts_canon = None
    if ts is not None:
        ts = np.asarray(ts, dtype=float)
        if ts.shape[0] != N_AAL90:
            raise ValueError(f"Expected ts with 90 rows, got {ts.shape}.")
        ts_canon = ts[p, :]
    return SC_canon, ts_canon, info


# ─────────────────────────────────────────────────────────────────────────────
# Ground-truth route (use whenever centroids are available)
# ─────────────────────────────────────────────────────────────────────────────
def order_from_centroids(xyz):
    """
    Exact O(n) determination from MNI centroid coordinates.

    Parameters
    ----------
    xyz : (90, 3) array_like
        Parcel centroids in MNI space. Left hemisphere has x < 0.

    Returns
    -------
    order : {'non-symmetric LLRR', 'symmetric LRLR', 'other'}
    is_left : (90,) bool ndarray
        Per-parcel hemisphere assignment. This is the piece the connectivity-
        based test cannot recover.
    """
    xyz = np.asarray(xyz, dtype=float)
    if xyz.shape != (N_AAL90, 3):
        raise ValueError(f"Expected (90, 3) centroids, got {xyz.shape}.")
    is_left = xyz[:, 0] < 0

    # Contiguous blocks (either L-first or R-first)
    if (is_left[:N_HEMI].all() and (~is_left[N_HEMI:]).all()) or \
       ((~is_left[:N_HEMI]).all() and is_left[N_HEMI:].all()):
        return ORDER_LLRR, is_left
    # Interleaved (either parity)
    if (is_left[0::2].all() and (~is_left[1::2]).all()) or \
       ((~is_left[0::2]).all() and is_left[1::2].all()):
        return ORDER_LRLR, is_left
    return "other", is_left


# ─────────────────────────────────────────────────────────────────────────────
# Legacy label mapping (for downstream code keyed to the original strings)
# ─────────────────────────────────────────────────────────────────────────────
LEGACY_LABELS = {
    ORDER_LLRR: "non-symmetric",
    ORDER_LRLR: "symmetric",
    ORDER_UNDET: None,
}


if __name__ == "__main__":
    # Reference behaviour on a deterministic DWI-derived AAL90 SC in LLRR order:
    #   log10(intra/inter) ~ +0.89 (H1) vs ~ -0.13 (H2)
    #   |log10(LL/RR)|     ~  0.01 (both; uninformative unless the split is wrong)
    #   Spearman LL~RR     ~  0.12 (H1) vs ~ 0.21 (H2)  <- weak, hence down-weighted
    #   TOTAL              ~ +0.93 vs -0.03  ->  margin ~ 0.97  ->  ORDER_LLRR
    pass