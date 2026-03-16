## Major Questions/Updates
- Starting manuscript
- Plots *mostly* finalized
- Consider doing ISA??
    - DESeq gives us 112 (rect) and 22 (vag) significant ASVs (should we consider making threshold more stringent???)
    - Core Microbiome gave us only 1 species in CPP-Endo for rectal and 0 species in vaginal ...

# Meeting 8 - March 16th

## Aim 2 - Phylum Composition
### Added new phylum-level taxonomy composition barplots (used to be genera)
![Phylum Barplot](../results/aim2/03-tax_composition/03_taxa_phylum_barplot.png)

## Aim 2 - DESeq2 Phylum

![CPP vs Control](../results/aim2/05-deseq2/05-cpp_control_contrast.png)
![CPP Endo vs Control](../results/aim2/05-deseq2/05-cpp_endo_control_contrast.png)
![CPP Endo vs CPP](../results/aim2/05-deseq2/05-cpp_cpp_endo_contrast.png)

### Updated DESeq2 differential abundance to phylum level (used to be genera)
![Significant ASVs Rectal](../results/aim2/05-deseq2/05-sigASVs_rect.png)
![Significant ASVs Vaginal](../results/aim2/05-deseq2/05-sigASVs_vag.png)

## Aim 3 - Volcano Plots
### Aim 3: Made DESeq2 thresholds more stringent for KO differential abundance (log2FC = 4, padj < 0.01)
![Rectal CPP vs Control](../results/aim3/rectal_CPP_vs_Control_volcano.png)
![Rectal Endo vs Control](../results/aim3/rectal_Endo_vs_Control_volcano.png)
![Rectal Endo vs CPP](../results/aim3/rectal_Endo_vs_CPP_volcano.png)
![Vaginal CPP vs Control](../results/aim3/vaginal_CPP_vs_Control_volcano.png)
![Vaginal Endo vs Control](../results/aim3/vaginal_Endo_vs_Control_volcano.png)
![Vaginal Endo vs CPP](../results/aim3/vaginal_Endo_vs_CPP_volcano.png)

- currently screening through significant KO terms that are differentially abundant in CPP-Endo/CPP-only contrasts (for vaginal (20) and rectal (9) -> 29 total)
