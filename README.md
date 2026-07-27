# AAL-90 to Yeo7 Mapping

This folder collects lightweight resources for assigning the 90 parcels of the Automated Anatomical Labeling (AAL-90) atlas to the canonical Yeo 2011 7-network parcellation (Yeo7). It includes the atlas image, the structural connectivity (SC) matrix used in our analyses, the manuscript describing the mapping, and a worked Python notebook that demonstrates how to load the data and aggregate parcel-wise measures at the network level.

## Repository contents

- `AAL90_atlas.nii` – 3D NIfTI volume with the AAL-90 parcellation in MNI space (left-hemisphere parcels first, followed by right-hemisphere parcels). Use this file to extract parcel-wise time-series or images that need to be summarized per network.
- `SC_AAL90.mat` – MATLAB file with a 90×90 structural connectivity matrix aligned to the atlas ordering. The example notebook shows how to infer whether the ordering is L→R blockwise (`non-symmetric`) or interleaved.
- `functions_infer_aal90_order_v1.py` – standalone, importable module used by the notebook to infer parcel ordering from a 90×90 SC/FC matrix via a block-density-ratio statistic (robust to the corpus-callosum crossing-fibre problem that fools a naive diagonal-based check). Also exposes `canonicalize_to_llrr` and an exact centroid-based check (`order_from_centroids`) for when MNI coordinates are available.
- `YEO7_AAL90.txt` – authoritative source for the AAL90→Yeo7 label dictionaries (`Yeo7_to_AAL90_labels`, cortex only, and `Yeo7subcort_to_AAL90_labels`, with a subcortical extension) plus the canonical AAL-90 label order. The notebook mirrors these dictionaries directly.
- `AAL90_to_YEO7_mapping.pdf` (source: `AAL90_to_YEO7_mapping.tex`) – manuscript detailing the rationale for each Yeo7 network definition and reference tables listing the exact parcel names and 1-based indices assigned to every network, including the subcortical extension.
- `example_AAL90_to_YEO7_mapping.ipynb` – Python walkthrough (NumPy, SciPy, Matplotlib) that:
  1. loads `SC_AAL90.mat`, visualizes the connectivity matrix, and infers the parcel ordering using `functions_infer_aal90_order_v1.infer_aal90_order`;
  2. defines `Yeo7_to_AAL90_labels` (7 cortical networks: Visual, Somatomotor, Dorsal Attention, Ventral Attention/Salience, Limbic, Frontal Control, Default Mode) and `Yeo7subcort_to_AAL90_labels` (same 7 plus a Subcortical Extension network covering Caudate/Putamen/Pallidum/Thalamus), matching `YEO7_AAL90.txt`;
  3. builds network→index lookups, and demonstrates how to aggregate parcel-wise vectors (here synthetic values) into group-level bar plots per network.

## Getting started

1. **Environment** – Install Python ≥3.9 with `numpy`, `scipy`, and `matplotlib`. Launch Jupyter Lab/Notebook if you want to run the demo notebook.
2. **Inspect the atlas** – Load `AAL90_atlas.nii` with your preferred neuroimaging library (`nibabel`, `nilearn`, etc.) to extract parcel masks or to verify registrations.
3. **Map parcels to networks** – Either read `AAL90_to_YEO7_mapping.pdf` for human-readable tables, or import `Yeo7_to_AAL90_labels` / `Yeo7subcort_to_AAL90_labels` from `YEO7_AAL90.txt` (or the notebook) to automatically translate AAL parcel-level metrics (BOLD amplitudes, dynamical measures, etc.) into Yeo7 network summaries. Use the subcortical-extension variant if you need all 90 parcels covered; the cortex-only variant leaves Caudate/Putamen/Pallidum/Thalamus unassigned.
4. **Reproducibility** – When using or extending the mapping, cite the references listed in `AAL90_to_YEO7_mapping.pdf`, especially the AAL atlas papers and Yeo et al. 2011 (the source of the 7-network parcellation).

## Notes and suggestions

- The atlas and SC files assume the standard (left-all, then right-all) ordering recommended for AAL-90. `functions_infer_aal90_order_v1.infer_aal90_order` verifies this automatically for any SC/FC matrix you load, and raises rather than guessing if the decision margin is too small — resolve against the atlas label file in that case.
- `Yeo7subcort_to_AAL90_labels` is built by extending `Yeo7_to_AAL90_labels` (`{**Yeo7_to_AAL90_labels, "Subcortical extension": [...]}`) rather than duplicating every entry, so the two dictionaries cannot drift out of sync when the cortical mapping is edited.
- Feel free to export the network tables to CSV/JSON by extending the notebook; the dictionaries already contain the parcel names needed to do so programmatically.
- If you adapt the mapping (e.g., to change limbic or subcortical assignments), keep `YEO7_AAL90.txt`, the notebook, and `AAL90_to_YEO7_mapping.tex` in sync so `AAL90_to_YEO7_mapping.pdf` remains the authoritative reference. Recompile with `pdflatex AAL90_to_YEO7_mapping.tex` (run twice for stable pagination).
