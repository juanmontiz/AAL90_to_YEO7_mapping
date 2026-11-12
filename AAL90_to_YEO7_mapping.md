# AAL90 to Yeo7 mapping

**Juan Manuel Monti** · corrected version, November 12, 2026

> Markdown version of `AAL90_to_YEO7_mapping.tex` / `.pdf`, with the same content.
> The region table at the end is generated from `aal90_labels_YEO7.csv`.

## Summary

Each of the 90 cerebral regions of the Automated Anatomical Labeling atlas (AAL90) is assigned to one of the seven resting-state networks of Yeo et al. (2011) (Yeo7).
A region takes the network it shares the most voxels with in MNI space.
The four subcortical nuclei, which the cortical Yeo atlas does not cover, are left unassigned.
This version replaces an earlier hand-made mapping (see [Changes from the earlier version](#changes-from-the-earlier-version)).

## Atlases

- **Regions:** AAL atlas, SPM12 version, MNI space [Tzourio-Mazoyer 2002]. Only its 90 cerebral labels are used; cerebellar and vermis labels are ignored.
- **Networks:** Yeo et al. 7-network parcellation [Yeo 2011], obtained by clustering resting-state functional connectivity in 1000 healthy subjects. We use the liberal ("thick") cortical mask, MNI152 FreeSurfer-conformed, 1 mm.

The seven networks are Visual, Somatomotor, Dorsal attention, Ventral attention (salience), Limbic, Frontoparietal (control) and Default.

## Assignment criterion

No standard AAL90-to-Yeo7 table exists in the literature, and each study builds its own [Doucet 2019].
We use the most common objective rule, the largest spatial overlap between the two atlases:

1. Resample the AAL atlas onto the 1 mm Yeo grid with nearest-neighbour interpolation.
2. For each AAL region, pool the voxels of the left and right homologues, so both hemispheres of a region receive the same network.
3. Count the voxels of each Yeo network inside the region, ignoring voxels outside the Yeo mask (label 0). The network with the most voxels is the region's network.

Two exceptions:

- **`Caudate`, `Putamen`, `Pallidum` and `Thalamus` are unassigned.** The Yeo 2011 atlas is cortex-only and has no label for them. They stay in every region-level and whole-brain analysis but are excluded from network-level summaries, which therefore use 82 of the 90 regions.
- **`Hippocampus` → Default and `Amygdala` → Limbic, by convention.** Both lie mostly outside the cortical Yeo mask (less than 10% of their voxels carry a Yeo label), so their overlap is not reliable.

The assignment is computed by `make_aal90_labels_YEO7.py` and stored in `aal90_labels_YEO7.csv`. That CSV is the reference for this document, for `YEO7_AAL90.txt` and for the example notebook.

## Resulting networks

Every region is bilateral, so each network holds the left and right copy of each region listed.

| Network | AAL regions | n |
|---|---|--:|
| Visual | `Calcarine`, `Cuneus`, `Lingual`, `Occipital_Sup`, `Occipital_Mid`, `Occipital_Inf`, `Fusiform` | 14 |
| Somatomotor | `Precentral`, `Rolandic_Oper`, `Supp_Motor_Area`, `Postcentral`, `Paracentral_Lobule`, `Heschl`, `Temporal_Sup` | 14 |
| Dorsal attention | `Parietal_Sup` | 2 |
| Ventral attention / salience | `Insula`, `Cingulum_Mid`, `SupraMarginal` | 6 |
| Limbic | `Frontal_Sup_Orb`, `Olfactory`, `Rectus`, `ParaHippocampal`, `Amygdala`*, `Temporal_Pole_Sup`, `Temporal_Pole_Mid`, `Temporal_Inf` | 16 |
| Frontoparietal control | `Frontal_Mid`, `Frontal_Mid_Orb`, `Frontal_Inf_Oper`, `Frontal_Inf_Tri`, `Parietal_Inf` | 10 |
| Default | `Frontal_Sup`, `Frontal_Inf_Orb`, `Frontal_Sup_Medial`, `Frontal_Med_Orb`, `Cingulum_Ant`, `Cingulum_Post`, `Hippocampus`*, `Angular`, `Precuneus`, `Temporal_Mid` | 20 |
| *Unassigned* | `Caudate`, `Putamen`, `Pallidum`, `Thalamus` | 8 |

\* assigned by convention.

As in the original Yeo7 solution, auditory cortex (`Heschl`, `Temporal_Sup`) falls in the somatomotor network [Yeo 2011].

## Caveats

- AAL regions are anatomical, so large regions span several Yeo networks. A single label per region is always an approximation.
- Seven regions have an overlap fraction below 0.5, meaning the chosen network holds less than half of the region's Yeo-labelled voxels:

  | Region | Network | Overlap |
  |---|---|--:|
  | `Hippocampus` | Default (convention) | 0.33 |
  | `Temporal_Inf` | Limbic | 0.36 |
  | `Precuneus` | Default | 0.38 |
  | `Frontal_Inf_Oper` | Frontoparietal | 0.38 |
  | `Supp_Motor_Area` | Somatomotor | 0.40 |
  | `Cingulum_Mid` | Ventral attention | 0.45 |
  | `Parietal_Inf` | Frontoparietal (vs Dorsal attention) | 0.48 |

- Pooling hemispheres hides a few within-hemisphere close calls:
  - `Frontal_Inf_Oper_L`: Ventral attention 0.39 vs Frontoparietal 0.37.
  - `Supp_Motor_Area_L`: Ventral attention 0.36 vs Somatomotor 0.33.
  - `Precuneus_R`: Dorsal attention 0.35 vs Default 0.34.
- Dorsal attention holds only two regions, so network-level averages for it are noisy.
- Results depend slightly on the Yeo mask (thick or thin) and on the resampling direction.

## Suggested Methods sentence

> Each AAL90 region was assigned to one of the seven Yeo et al. (2011) resting-state networks according to maximal voxel overlap in MNI space (homologous regions pooled across hemispheres). Subcortical nuclei (caudate, putamen, pallidum, thalamus), which the cortical Yeo atlas does not cover, were excluded from network-level summaries, and hippocampus and amygdala, which lie outside the cortical Yeo mask, were assigned to the default and limbic networks, respectively.

## Changes from the earlier version

The earlier mapping (December 2025) was built by hand from the literature and differed from the atlas overlap for nine regions:

| Region | Earlier | Now |
|---|---|---|
| `Frontal_Sup` | Frontoparietal | Default |
| `Frontal_Inf_Orb` | Limbic | Default |
| `Frontal_Med_Orb` | Limbic | Default |
| `Cingulum_Ant` | Ventral attention | Default |
| `Frontal_Mid_Orb` | Limbic | Frontoparietal |
| `Temporal_Sup` | Ventral attention | Somatomotor |
| `ParaHippocampal` | Default | Limbic |
| `Temporal_Inf` | Default | Limbic |
| `Parietal_Inf` | Dorsal attention, or both Dorsal attention and Frontoparietal | Frontoparietal |

In addition, the example notebook listed `Frontal_Mid_Orb` twice in Limbic and omitted `Frontal_Med_Orb`.
The four subcortical nuclei formed an extra "Subcortical" block, which is now removed.
**Analyses that used the earlier mapping should be rerun.**

## Region table

- **LLRR idx:** 1-based position in the standard AAL order (L1..L45, then R1..R45).
- **Data idx:** position in the mirrored order of `SC_AAL90.mat` (L1..L45, then R45..R1).
- **Overlap:** fraction of the region's Yeo-labelled voxels (left and right pooled) that fall in the chosen network.
- **Coverage:** fraction of the region's voxels that carry any Yeo label.

| Network | AAL label | LLRR idx | Data idx | Overlap | Coverage | Source |
|---|---|--:|--:|--:|--:|---|
| **Visual** | `Calcarine_L` | 22 | 22 | 0.947 | 0.801 | overlap |
|  | `Calcarine_R` | 67 | 69 | 0.947 | 0.801 | overlap |
|  | `Cuneus_L` | 23 | 23 | 0.894 | 0.842 | overlap |
|  | `Cuneus_R` | 68 | 68 | 0.894 | 0.842 | overlap |
|  | `Lingual_L` | 24 | 24 | 0.998 | 0.837 | overlap |
|  | `Lingual_R` | 69 | 67 | 0.998 | 0.837 | overlap |
|  | `Occipital_Sup_L` | 25 | 25 | 0.877 | 0.753 | overlap |
|  | `Occipital_Sup_R` | 70 | 66 | 0.877 | 0.753 | overlap |
|  | `Occipital_Mid_L` | 26 | 26 | 0.745 | 0.834 | overlap |
|  | `Occipital_Mid_R` | 71 | 65 | 0.745 | 0.834 | overlap |
|  | `Occipital_Inf_L` | 27 | 27 | 0.944 | 0.799 | overlap |
|  | `Occipital_Inf_R` | 72 | 64 | 0.944 | 0.799 | overlap |
|  | `Fusiform_L` | 28 | 28 | 0.621 | 0.867 | overlap |
|  | `Fusiform_R` | 73 | 63 | 0.621 | 0.867 | overlap |
| **Somatomotor** | `Precentral_L` | 1 | 1 | 0.597 | 0.741 | overlap |
|  | `Precentral_R` | 46 | 90 | 0.597 | 0.741 | overlap |
|  | `Rolandic_Oper_L` | 9 | 9 | 0.700 | 0.866 | overlap |
|  | `Rolandic_Oper_R` | 54 | 82 | 0.700 | 0.866 | overlap |
|  | `Supp_Motor_Area_L` | 10 | 10 | 0.402 | 0.766 | overlap |
|  | `Supp_Motor_Area_R` | 55 | 81 | 0.402 | 0.766 | overlap |
|  | `Postcentral_L` | 29 | 29 | 0.783 | 0.789 | overlap |
|  | `Postcentral_R` | 74 | 62 | 0.783 | 0.789 | overlap |
|  | `Paracentral_Lobule_L` | 35 | 35 | 0.933 | 0.681 | overlap |
|  | `Paracentral_Lobule_R` | 80 | 56 | 0.933 | 0.681 | overlap |
|  | `Heschl_L` | 40 | 40 | 0.988 | 0.920 | overlap |
|  | `Heschl_R` | 85 | 51 | 0.988 | 0.920 | overlap |
|  | `Temporal_Sup_L` | 41 | 41 | 0.627 | 0.864 | overlap |
|  | `Temporal_Sup_R` | 86 | 50 | 0.627 | 0.864 | overlap |
| **Dorsal attention** | `Parietal_Sup_L` | 30 | 30 | 0.779 | 0.815 | overlap |
|  | `Parietal_Sup_R` | 75 | 61 | 0.779 | 0.815 | overlap |
| **Ventral attention** | `Insula_L` | 15 | 15 | 0.617 | 0.863 | overlap |
|  | `Insula_R` | 60 | 76 | 0.617 | 0.863 | overlap |
|  | `Cingulum_Mid_L` | 17 | 17 | 0.449 | 0.886 | overlap |
|  | `Cingulum_Mid_R` | 62 | 74 | 0.449 | 0.886 | overlap |
|  | `SupraMarginal_L` | 32 | 32 | 0.599 | 0.854 | overlap |
|  | `SupraMarginal_R` | 77 | 59 | 0.599 | 0.854 | overlap |
| **Limbic** | `Frontal_Sup_Orb_L` | 3 | 3 | 0.716 | 0.804 | overlap |
|  | `Frontal_Sup_Orb_R` | 48 | 88 | 0.716 | 0.804 | overlap |
|  | `Olfactory_L` | 11 | 11 | 0.943 | 0.455 | overlap |
|  | `Olfactory_R` | 56 | 80 | 0.943 | 0.455 | overlap |
|  | `Rectus_L` | 14 | 14 | 0.917 | 0.867 | overlap |
|  | `Rectus_R` | 59 | 77 | 0.917 | 0.867 | overlap |
|  | `ParaHippocampal_L` | 20 | 20 | 0.537 | 0.684 | overlap |
|  | `ParaHippocampal_R` | 65 | 71 | 0.537 | 0.684 | overlap |
|  | `Amygdala_L` | 21 | 21 | 1.000 | 0.037 | convention |
|  | `Amygdala_R` | 66 | 70 | 1.000 | 0.037 | convention |
|  | `Temporal_Pole_Sup_L` | 42 | 42 | 0.550 | 0.733 | overlap |
|  | `Temporal_Pole_Sup_R` | 87 | 49 | 0.550 | 0.733 | overlap |
|  | `Temporal_Pole_Mid_L` | 44 | 44 | 0.647 | 0.741 | overlap |
|  | `Temporal_Pole_Mid_R` | 89 | 47 | 0.647 | 0.741 | overlap |
|  | `Temporal_Inf_L` | 45 | 45 | 0.361 | 0.829 | overlap |
|  | `Temporal_Inf_R` | 90 | 46 | 0.361 | 0.829 | overlap |
| **Frontoparietal** | `Frontal_Mid_L` | 4 | 4 | 0.642 | 0.814 | overlap |
|  | `Frontal_Mid_R` | 49 | 87 | 0.642 | 0.814 | overlap |
|  | `Frontal_Mid_Orb_L` | 5 | 5 | 0.649 | 0.807 | overlap |
|  | `Frontal_Mid_Orb_R` | 50 | 86 | 0.649 | 0.807 | overlap |
|  | `Frontal_Inf_Oper_L` | 6 | 6 | 0.379 | 0.775 | overlap |
|  | `Frontal_Inf_Oper_R` | 51 | 85 | 0.379 | 0.775 | overlap |
|  | `Frontal_Inf_Tri_L` | 7 | 7 | 0.572 | 0.738 | overlap |
|  | `Frontal_Inf_Tri_R` | 52 | 84 | 0.572 | 0.738 | overlap |
|  | `Parietal_Inf_L` | 31 | 31 | 0.483 | 0.926 | overlap |
|  | `Parietal_Inf_R` | 76 | 60 | 0.483 | 0.926 | overlap |
| **Default** | `Frontal_Sup_L` | 2 | 2 | 0.500 | 0.728 | overlap |
|  | `Frontal_Sup_R` | 47 | 89 | 0.500 | 0.728 | overlap |
|  | `Frontal_Inf_Orb_L` | 8 | 8 | 0.597 | 0.755 | overlap |
|  | `Frontal_Inf_Orb_R` | 53 | 83 | 0.597 | 0.755 | overlap |
|  | `Frontal_Sup_Medial_L` | 12 | 12 | 0.841 | 0.789 | overlap |
|  | `Frontal_Sup_Medial_R` | 57 | 79 | 0.841 | 0.789 | overlap |
|  | `Frontal_Med_Orb_L` | 13 | 13 | 0.734 | 0.804 | overlap |
|  | `Frontal_Med_Orb_R` | 58 | 78 | 0.734 | 0.804 | overlap |
|  | `Cingulum_Ant_L` | 16 | 16 | 0.692 | 0.840 | overlap |
|  | `Cingulum_Ant_R` | 61 | 75 | 0.692 | 0.840 | overlap |
|  | `Cingulum_Post_L` | 18 | 18 | 0.967 | 0.611 | overlap |
|  | `Cingulum_Post_R` | 63 | 73 | 0.967 | 0.611 | overlap |
|  | `Hippocampus_L` | 19 | 19 | 0.329 | 0.022 | convention |
|  | `Hippocampus_R` | 64 | 72 | 0.329 | 0.022 | convention |
|  | `Angular_L` | 33 | 33 | 0.779 | 0.834 | overlap |
|  | `Angular_R` | 78 | 58 | 0.779 | 0.834 | overlap |
|  | `Precuneus_L` | 34 | 34 | 0.375 | 0.800 | overlap |
|  | `Precuneus_R` | 79 | 57 | 0.375 | 0.800 | overlap |
|  | `Temporal_Mid_L` | 43 | 43 | 0.568 | 0.832 | overlap |
|  | `Temporal_Mid_R` | 88 | 48 | 0.568 | 0.832 | overlap |
| **Unassigned** | `Caudate_L` | 36 | 36 | — | 0.011 | — |
|  | `Caudate_R` | 81 | 55 | — | 0.011 | — |
|  | `Putamen_L` | 37 | 37 | — | 0.014 | — |
|  | `Putamen_R` | 82 | 54 | — | 0.014 | — |
|  | `Pallidum_L` | 38 | 38 | — | 0.000 | — |
|  | `Pallidum_R` | 83 | 53 | — | 0.000 | — |
|  | `Thalamus_L` | 39 | 39 | — | 0.000 | — |
|  | `Thalamus_R` | 84 | 52 | — | 0.000 | — |

## References

- **[Yeo 2011]** Yeo BTT, Krienen FM, Sepulcre J, Sabuncu MR, Lashkari D, Hollinshead M, Roffman JL, Smoller JW, Zöllei L, Polimeni JR, Fischl B, Liu H, Buckner RL. The organization of the human cerebral cortex estimated by intrinsic functional connectivity. *Journal of Neurophysiology* 106(3):1125–1165, 2011. [doi:10.1152/jn.00338.2011](https://doi.org/10.1152/jn.00338.2011)
- **[Tzourio-Mazoyer 2002]** Tzourio-Mazoyer N, Landeau B, Papathanassiou D, Crivello F, Etard O, Delcroix N, Mazoyer B, Joliot M. Automated anatomical labeling of activations in SPM using a macroscopic anatomical parcellation of the MNI MRI single-subject brain. *NeuroImage* 15(1):273–289, 2002. [doi:10.1006/nimg.2001.0978](https://doi.org/10.1006/nimg.2001.0978)
- **[Doucet 2019]** Doucet GE, Lee WH, Frangou S. Evaluation of the spatial variability in the major resting-state networks across human brain functional atlases. *Human Brain Mapping* 40(15):4577–4587, 2019. [doi:10.1002/hbm.24722](https://doi.org/10.1002/hbm.24722)
