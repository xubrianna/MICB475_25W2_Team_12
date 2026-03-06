## Updates 

# Meeting 7 - March 9th

## New Plots

### Aim 1

#### Faith PD Boxplot
![Faith PD Boxplot](../results/aim1/01_faith_PD_boxplot.png)
> Alpha diversity (Faith's PD) compared across disease groups using boxplots

#### Faith PD Violin
![Faith PD Violin](../results/aim1/01_faith_PD_violin.png)
> Violin plot of Faith's PD showing the distribution shape of alpha diversity within each group.

#### PCoA
![PCoA](../results/aim1/02-pcoa.png)
> PCoA ordination (beta diversity) showing separation of samples by disease group

#### PCoA with Ellipse
![PCoA Ellipse](../results/aim1/02-pcoa_ellipse.png)
> Same PCoA with 95% confidence ellipses per group 

---

### Aim 2

#### Taxa Barplot - Rectal
![Taxa Barplot Rectal](../results/aim2/03-tax_composition/01_taxa_barplot_rectal.png)
> Taxonomic composition of rectal samples at the genus level, grouped by disease status.

#### Taxa Barplot - Vaginal
![Taxa Barplot Vaginal](../results/aim2/03-tax_composition/01_taxa_barplot_vaginal.png)
> Taxonomic composition of vaginal samples at the genus level, grouped by disease status.

#### Taxa Barplot (Combined)
![Taxa Barplot](../results/aim2/03-tax_composition/03_taxa_barplot.png)
> Combined taxonomy barplot across both body sites for an overview of community structure.

#### Venn Diagram - CPP Only vs Healthy (0 detection + 0.5 prevalence)
![Venn all diseases (rectal 0)](../results/aim2/04-core_microbiome/venn_all_diseases_rectal_0.png)
![Venn all diseases (vaginal 0)](../results/aim2/04-core_microbiome/venn_all_diseases_vaginal_0.png)

#### Venn Diagram - CPP Only vs Healthy (0.001 detection + 0.5 prevalence)
![Venn all diseases (rectal 0.001)](../results/aim2/04-core_microbiome/venn_all_diseases_rectal_0.001.png)
![Venn all diseases (vaginal 0.001)](../results/aim2/04-core_microbiome/venn_all_diseases_vaginal_0.001.png)

#### Venn Diagram - CPP Only vs Healthy (0.01 detection + 0.5 prevalence)
![Venn all diseases (rectal 0.01)](../results/aim2/04-core_microbiome/venn_all_diseases_rectal_0.01.png)
![Venn all diseases (vaginal 0.01)](../results/aim2/04-core_microbiome/venn_all_diseases_vaginal_0.01.png)


#### DESeq Volcano - CPP vs Control (sanity check)
![CPP vs Control](../results/aim2/05-deseq2/05-cpp_control_contrast.png)
> Volcano plots (rectal & vaginal) of differentially abundant taxa between CPP and Control. Red points pass padj < 0.05 and |log2FC| > 1.

#### DESeq Volcano - CPP Endo vs Control (sanity check)
![CPP Endo vs Control](../results/aim2/05-deseq2/05-cpp_endo_control_contrast.png)
> Volcano plots (rectal & vaginal) of differentially abundant taxa between CPP Endo and Control.

#### DESeq Volcano - CPP Endo vs CPP
![CPP Endo vs CPP](../results/aim2/05-deseq2/05-cpp_endo_cpp_contrast.png)
> Volcano plots (rectal & vaginal) of differentially abundant taxa between CPP Endo and CPP-only groups.

#### Significant ASVs - Rectal
![SigASVs Rectal](../results/aim2/05-deseq2/05-sigASVs_rect.png)
> Bar plot of significantly differentially abundant genera in rectal samples (CPP Endo vs CPP), with log2 fold change and standard error.

#### Significant ASVs - Vaginal
![SigASVs Vaginal](../results/aim2/05-deseq2/05-sigASVs_vag.png)
> Bar plot of significantly differentially abundant genera in vaginal samples (CPP Endo vs CPP), with log2 fold change and standard error.


