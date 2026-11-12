# AAL-90 to Yeo7 Mapping

Resources for assigning the 90 cerebral regions of the Automated Anatomical Labeling atlas (AAL-90) to the seven resting-state networks of Yeo et al. (2011) (Yeo7).
Each region takes the Yeo network it overlaps most in MNI space.

> **Corrected 2026-11-12.** The earlier mapping was hand-made and differed from the atlas overlap for 9 regions; some copies also listed `Frontal_Mid_Orb` twice and omitted `Frontal_Med_Orb`, or placed `Parietal_Inf` in two networks.
> The subcortical "extension" block has been removed.
> Analyses that used the earlier mapping should be rerun.
> See the "Changes" section of `AAL90_to_YEO7_mapping.pdf`.

## Contents

| File | What it is |
|---|---|
| `aal90_labels_YEO7.csv` | **Reference table.** One row per region: `index,label,class,yeo7_network,yeo7_network_name,overlap_frac,yeo_coverage,source`. Rows are in the mirrored data order (see below). |
| `make_aal90_labels_YEO7.py` | Builds the CSV from the atlases (needs `nilearn`, which downloads AAL SPM12 and Yeo 2011). |
| `YEO7_AAL90.txt` | Python-readable dictionary `Yeo7_to_AAL90_labels` (network → AAL labels, built from the CSV) plus the AAL label orders `AAL_labels_non_symmetric` (LLRR) and `AAL_labels_data_order` (mirrored). Load with `exec(open('YEO7_AAL90.txt').read(), ns)`. |
| `example_AAL90_to_YEO7_mapping.ipynb` | Worked example: loads `SC_AAL90.mat`, infers its parcel order, builds network → index lookups and plots network averages of synthetic data. |
| `functions_infer_aal90_order_v1.py` | Infers the parcel order of a 90×90 SC/FC matrix (`infer_aal90_order`) and checks a label list against it (`assert_aal90_order`). |
| `AAL90_to_YEO7_mapping.pdf` / `.tex` | Method, resulting networks, caveats, a suggested Methods sentence and the full region table. |
| `AAL90_atlas.nii` | AAL-90 volume in MNI space, SPM interleaved label values (odd = left, even = right; 1 = Precentral_L, 2 = Precentral_R). |
| `SC_AAL90.mat` | 90×90 structural connectivity matrix in **mirrored** order. |

## Assignment criterion

1. Resample AAL (SPM12) onto the 1 mm Yeo 2011 grid (7 networks, liberal mask) with nearest-neighbour interpolation.
2. Pool the voxels of the left and right homologues of each region.
3. Count Yeo-labelled voxels per network; the network with the most voxels wins.

Exceptions:

- **Caudate, Putamen, Pallidum, Thalamus are unassigned.** Yeo 2011 is cortex-only. Keep them in region-level and whole-brain analyses; network summaries use the other 82 regions.
- **Hippocampus → Default, Amygdala → Limbic** by convention (under 10% of their voxels carry a Yeo label).

| Network | Regions (left and right of each) | n |
|---|---|---|
| Visual | Calcarine, Cuneus, Lingual, Occipital_Sup, Occipital_Mid, Occipital_Inf, Fusiform | 14 |
| Somatomotor | Precentral, Rolandic_Oper, Supp_Motor_Area, Postcentral, Paracentral_Lobule, Heschl, Temporal_Sup | 14 |
| Dorsal attention | Parietal_Sup | 2 |
| Ventral attention / salience | Insula, Cingulum_Mid, SupraMarginal | 6 |
| Limbic | Frontal_Sup_Orb, Olfactory, Rectus, ParaHippocampal, Amygdala, Temporal_Pole_Sup, Temporal_Pole_Mid, Temporal_Inf | 16 |
| Frontoparietal control | Frontal_Mid, Frontal_Mid_Orb, Frontal_Inf_Oper, Frontal_Inf_Tri, Parietal_Inf | 10 |
| Default | Frontal_Sup, Frontal_Inf_Orb, Frontal_Sup_Medial, Frontal_Med_Orb, Cingulum_Ant, Cingulum_Post, Hippocampus, Angular, Precuneus, Temporal_Mid | 20 |
| *(unassigned)* | Caudate, Putamen, Pallidum, Thalamus | 8 |

Assignments that are not clean-cut (overlap fraction below 0.5) are listed in the PDF.

## Parcel order

`SC_AAL90.mat` (and the sleep W/N1/N2/N3 and COMA AAL-90 data it goes with) is in **mirrored** order: L1..L45 in standard AAL order, then R45..R1.
Column `c ≥ 45` holds right region `89 − c`.
`infer_aal90_order` decides among LLRR, mirrored and interleaved by how similar the LL and RR blocks are once homotopic regions are paired (on `SC_AAL90.mat`: mirrored r = 0.97, LLRR 0.19, interleaved 0.12).
Always match regions **by name**, using the label list that matches your data order (`AAL_labels_data_order` for these data).

## Quick use

```python
import csv
net = {r["label"]: r["yeo7_network_name"] for r in csv.DictReader(open("aal90_labels_YEO7.csv"))}
labels = ...  # the 90 AAL labels in your data's column order
idx = {}
for i, lab in enumerate(labels):
    if net[lab]:                      # skip unassigned subcortical regions
        idx.setdefault(net[lab], []).append(i)
```

## Suggested Methods sentence

> Each AAL90 region was assigned to one of the seven Yeo et al. (2011) resting-state networks according to maximal voxel overlap in MNI space (homologous regions pooled across hemispheres). Subcortical nuclei (caudate, putamen, pallidum, thalamus), which the cortical Yeo atlas does not cover, were excluded from network-level summaries, and hippocampus and amygdala, which lie outside the cortical Yeo mask, were assigned to the default and limbic networks, respectively.

## References

- Yeo BTT et al. (2011) *J Neurophysiol* 106:1125–1165. doi:10.1152/jn.00338.2011
- Tzourio-Mazoyer N et al. (2002) *NeuroImage* 15:273–289. doi:10.1006/nimg.2001.0978
- Doucet GE, Lee WH, Frangou S (2019) *Hum Brain Mapp* 40:4577–4587. doi:10.1002/hbm.24722
