"""Assign AAL90 regions to Yeo 2011 7-network components by maximal voxel overlap.

Writes aal90_labels_YEO7.csv. Criterion and caveats are documented in README.md and
AAL90_to_YEO7_mapping.pdf.

Row order is L1..L45, R45..R1 (mirrored), matching the data column order
(AAL_labels_data_order in YEO7_AAL90.txt), the order of SC_AAL90.mat.
"""
import csv
from pathlib import Path

import numpy as np
from nilearn import datasets, image

HERE = Path(__file__).resolve().parent
OUT = HERE / "aal90_labels_YEO7.csv"

NAMES = {1: "Visual", 2: "Somatomotor", 3: "Dorsal attention", 4: "Ventral attention",
         5: "Limbic", 6: "Frontoparietal", 7: "Default"}
# Not covered by the cortex-only Yeo atlas: left unassigned (blank yeo7_network),
# so network-level summaries use the other 82 regions.
SUBCORTICAL = {"Caudate", "Putamen", "Pallidum", "Thalamus"}
# Regions largely outside the cortical Yeo mask: assigned by convention
CONVENTION = {"Hippocampus": 7, "Amygdala": 5}
MIN_COVERAGE = 0.10  # min fraction of region voxels carrying a Yeo label


def load_aal_order():
    ns = {}
    exec((HERE / "YEO7_AAL90.txt").read_text(), ns)
    # Data columns are mirrored: L1..L45, then R45..R1
    order = ns["AAL_labels_data_order"]
    assert len(order) == 90
    return order


def main():
    order = load_aal_order()

    aal = datasets.fetch_atlas_aal(version="SPM12")
    yeo = datasets.fetch_atlas_yeo_2011(n_networks=7, thickness="thick")
    yimg = image.load_img(yeo.maps)
    aimg = image.resample_to_img(aal.maps, yimg, interpolation="nearest",
                                 force_resample=True, copy_header=True)
    A = np.asarray(aimg.dataobj).squeeze().astype(int)
    Y = np.asarray(yimg.dataobj).squeeze().astype(int)
    code = dict(zip(aal.labels, map(int, aal.indices)))

    def overlap(base):
        # Pool L and R homologues so both hemispheres share one network
        y = np.concatenate([Y[A == code[base + "_L"]], Y[A == code[base + "_R"]]])
        cnt = np.bincount(y[y > 0], minlength=8)[1:]
        return cnt, (y > 0).mean()

    rows = []
    for i, label in enumerate(order, 1):
        base = label[:-2]
        cnt, cov = overlap(base)
        if base in SUBCORTICAL:
            cls, net, frac, src = "subcortical", None, None, "subcortical"
        elif cov < MIN_COVERAGE:
            assert base in CONVENTION, f"{base}: low Yeo coverage ({cov:.3f}) with no convention"
            cls, net, src = "cortical", CONVENTION[base], "convention (low Yeo coverage)"
            frac = cnt[net - 1] / cnt.sum()
        else:
            cls, net, src = "cortical", int(np.argmax(cnt)) + 1, "overlap"
            frac = cnt[net - 1] / cnt.sum()
        rows.append({"index": i, "label": label, "class": cls,
                     "yeo7_network": "" if net is None else net,
                     "yeo7_network_name": "" if net is None else NAMES[net],
                     "overlap_frac": "" if frac is None else f"{frac:.3f}",
                     "yeo_coverage": f"{cov:.3f}", "source": src})

    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT} ({len(rows)} regions)")


if __name__ == "__main__":
    main()
