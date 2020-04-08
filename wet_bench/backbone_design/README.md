## Design of CNV Backbone

**TLDR: The `uw_backbone.2020-03-04.bed` is the current CNV backbone design in OPXv7.**

## Background

Original BrocaHR probe design utilized ~3000 SNPs with variable population frequency. Since rare SNPs are unlikely to be helpful in identifying Loss-of-heterozygosity, the panel was redesigned starting with a set of IDT selected probes. Additional analysis of the original panel is in `Analysis of first round probes.ipynb`.

## Initial selection

Initially we asked IDT for the xGen CNV Backbone design. These are described in `xgen_cnv_backbone_panel_probes.bed`. Because these designs had higher density than we envisioned, we asked them to split these designes into 50% subsets, and this got us:
  * IDT_CNV_altA.bed
  * IDT_CNV_Balt.bed
  * IDT_CNV_Cov1.bed
  * IDT_CNV_2cov.bed

For the "cov" variants, an attempt was made to optimize coverage. However, this came at great expense of the spacing and density of the probes (essentially IDT just took the top 50% best performing probes)

## Balanced selection

Using an iterative optimization protocol (see `IDT probe optimizer.ipynb`), several probe designs were generated that had certain miminum spacing + coverage weighting penalties. These are described in `Backbone_selection.pdf`.

Visually examining the spacing and coverage distributions suggested that `selected_600000_14.bed` gave the best balance between coverage and spacing. 

It was also confirmed that the MAFs of the IDT selected probes are common and between 0.3 and 0.7 across EUR and AFR populations. See `Population frequencies.ipynb`.

## Final changes/edits

3 additional ChrY probes were added to the selected_600000_14.bed file:
```
chrY    6670400 6670520
chrY    15472802    15472922
chrY    21730196    21730316
```

## Ordering details

An order was placed for a 16-rxn discovery pool with order #16731366 and internal IDT name 03003092_Krumm_CNV_v2_Disco_SD.



