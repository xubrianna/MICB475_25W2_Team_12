## Updates 

# Meeting 7 - March 9th

## Key Takeaways

> **Aim 1 (Diversity):** No significant differences in alpha or beta diversity across disease groups — body site is the primary driver for variation, not disease status.

> **Aim 2 (Taxonomy):** Taxonomic composition is largely similar across disease groups within each body site. Rectal communities show more differentially abundant taxa than vaginal. Core microbiome analysis reveals very few unique taxa in CPP-Endo, especially in vaginal samples.

> **Aim 3 (Functional):** Predicted functional profiles (PICRUSt2 KOs) show more significant differences in vaginal samples than rectal — contrasting with the taxonomic results. 
---

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

![Venn diagram](..results/aim2/04-core_microbiome/venn_all_diseases_vaginal_0.001.png)

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


#### Discussion on results: 
### Aim 1
Alpha diversity (ML):
- No clear differences in Faith’s PD across disease groups (Control, CPP, CPP Endo), suggesting host disease status does not strongly affect alpha diversity.
- Rectal samples show higher Faith’s PD values than vaginal samples, indicating greater phylogenetic diversity in rectal microbiota.
- Patterns are consistent across sites, with variation within groups larger than differences between disease groups.

Beta diversity (ML):
- (Unweighted UniFrac)
- Clear separation by body site, with rectal and vaginal samples clustering apart, indicating distinct microbial community composition between sites.
- Substantial overlap between host disease groups (Control, CPP, CPP Endo) within each body site, suggesting disease status does not strongly drive community differences.
- PCoA Axis 1 explains most variation (24.7%), largely reflecting differences between sample types rather than disease groups.

### Aim 2
Taxa barplot (combined) (ML):
- Mean taxonomic composition (Genus level)
- Clear differences between body sites: vaginal samples are dominated by **Lactobacillus**, while rectal samples show a more diverse mix of genera.
- Similar overall composition across disease groups within each body site, suggesting disease status does not strongly shift dominant taxa.
- Rectal microbiota are more evenly distributed across multiple genera, whereas **vaginal microbiota are largely Lactobacillus-dominated**.
- Body site drives microbial composition more strongly than disease status.

Core microbiome: 
- Rectal samples at 0.5 prevalence and 0.001 detection (remove rare ASVs) has 3% of the ASVs shared that are unique to CPP-Endo -> this has potential 
- Rectal samples at 0.5 prevalence and 0 detection (absence/prevalence) has 4% of the ASVs shared that are unique to CPP-Endo -> this has potential 
- Rectal samples at 0.5 prevalence and 0.01 (keep abundance ASVs) detection has 1% of the ASVs shared that are unique to CPP-Endo -> this has potential 
- There appears to be a lot more that is shared among the healthy individuals than either disease state 
- Vaginal samples at 0.5 prevalence and 0.001 detection has 0% of the ASVs shared that are unique to CPP-Endo -> not much potential for distinguishing CPP-Endo
- Vaginal samples at 0.5 prevalence and 0 detection has 1% of the ASVs shared that are unique to CPP-Endo -> small potential 
- Vaginal samples at 0.5 prevalence and 0.001 detection has 0% of the ASVs shared that are unique to CPP-Endo -> not much potential
- Overall the rectal samples appear to be more shared amongst each disease state than vaginal microbiome which might be more similar between all patients regardless of host disease state 

DESeq2 differential abundance:
- **CPP vs Control (sanity check):** No significantly differentially abundant taxa after FDR correction in either body site — suggests CPP alone does not drive major taxonomic shifts compared to healthy controls
- **CPP Endo vs Control (sanity check):** Similar pattern; no strong differential abundance signal--- reinforces that disease effects are subtle
- **CPP Endo vs CPP (main contrast of interest):**
- Rectal: No taxa reached significance after multiple testing correction (lowest padj ~0.47). Several genera showed significant p-values (e.g., *Faecalitalea*, *Lactobacillus*, *Murdochiella*, *Alistipes*) but did not survive FDR adjustment
- Vaginal: *Fannyhessea* was the only significantly differentially abundant genus (padj = 0.035, log2FC = 4.37), enriched in CPP Endo relative to CPP-only. 
- All other vaginal genera had padj = 1, indicating no other significant differences


DESeq Volcano Plots (ML):
**CPP vs Control**
- Rectal samples show multiple significantly different ASVs, with both positive and negative log2 fold changes, indicating taxa enriched in either CPP or Control groups.
- Vaginal samples show very few significant taxa, suggesting minimal differential abundance between CPP and Control in the vaginal microbiome.
- Significant taxa in rectal samples show larger effect sizes and stronger statistical signals than those in vaginal samples.

**CPP Endo vs Control**
- Rectal microbiota again show many significant ASVs, with clear enrichment patterns in both directions, indicating disease-associated shifts in rectal taxa.
- Vaginal microbiota show fewer significant changes, though several taxa display strong positive fold changes in the CPP Endo group.
- Overall, rectal communities appear more responsive to disease status than vaginal communities.

**CPP Endo vs CPP**
- Rectal samples display numerous significantly different taxa, indicating microbial differences between CPP patients with and without endometriosis.
- Vaginal samples show only a small number of significant taxa, again suggesting weaker disease-associated shifts at this site.
- Some vaginal taxa exhibit large fold changes but relatively few pass significance thresholds, likely reflecting lower statistical power or higher variability.

TLDR;
-**Rectal microbiota show stronger and more consistent differential abundance across disease comparisons than vaginal microbiota.**
-**Vaginal microbiota appear more stable across disease states, with fewer taxa significantly associated with host disease status.**
-**These results suggest that disease-associated microbial shifts may be more detectable in rectal communities than in vaginal communities.**

### Aim 3

#### KO Functional Composition – Rectal

![Rectal KO PCoA](../results/aim3/KO_rectal_pcoa_ellipse.png)

#### KO Functional Composition – Vaginal

![Vaginal KO PCoA](../results/aim3/KO_vaginal_pcoa_ellipse.png)


PCA of KO functional profiles:
- PCA of predicted KEGG KO abundances (PICRUSt2) shows overlapping clusters for Control, CPP, and CPP Endo in both rectal and vaginal samples
- No clear separation by disease group in either body site, suggesting that predicted functional profiles are broadly similar across conditions???

#### Differential KO – Rectal (CPP vs Control)

![Rectal CPP vs Control Volcano](../results/aim3/rectal_CPP_vs_Control_volcano.png)

#### Differential KO – Rectal (Endo vs CPP) - Conparison of Interest

![Rectal Endo vs CPP Volcano](../results/aim3/rectal_Endo_vs_CPP_volcano.png)

#### Differential KO – Rectal (Endo vs Control)

![Rectal Endo vs Control Volcano](../results/aim3/rectal_Endo_vs_Control_volcano.png)

#### Differential KO – Vaginal (CPP vs Control)

![Vaginal CPP vs Control Volcano](../results/aim3/vaginal_CPP_vs_Control_volcano.png)


#### Differential KO – Vaginal (Endo vs CPP) - Comparison of Interest


![Vaginal Endo vs CPP Volcano](../results/aim3/vaginal_Endo_vs_CPP_volcano.png)

#### Differential KO – Vaginal (Endo vs Control)

![Vaginal Endo vs Control Volcano](../results/aim3/vaginal_Endo_vs_Control_volcano.png)


#### AIM3- Discussion

DESeq2 differential KO abundance:
- **Rectal – CPP vs Control:** Some KOs labelled as significant (red points above the padj < 0.05 threshold), indicating some predicted functional differences between CPP and healthy controls in rectal samples. 
- **Rectal – CPP Endo vs Control:** Similar pattern to CPP vs Control, with several significant KOs detected. Suggests that endometriosis-associated CPP shares some functional shifts with CPP-only relative to controls
- **Rectal – CPP Endo vs CPP (Comparison of Interest):** Few significant KOs compared to the case-control comparisons, indicating that the functional profile difference between CPP Endo and CPP-only is more subtle in rectal samples
- **Vaginal – CPP vs Control:** Many significant KOs detected, indicating substantial predicted functional differences between CPP and healthy controls in vaginal samples — more so than in rectal samples
- **Vaginal – CPP Endo vs Control:** Also shows numerous significant KOs, suggesting that CPP Endo is functionally distinct from controls at the vaginal site
- **Vaginal – CPP Endo vs CPP (Comparison of Interest):** Significant KOs present, indicating predicted functional differences between CPP Endo and CPP-only groups in vaginal samples

Overall:
- Interestingly, vaginal samples show more predicted functional differences across disease comparisons than rectal samples — this contrasts with the taxonomic results from Aim 2, where rectal samples had more differentially abundant taxa
- The stronger functional signal in vaginal samples despite limited taxonomic shifts highlights the value of functional profiling (PICRUSt2) as a complement to taxonomic analysis
- Next steps: mapping significant KOs to KEGG pathways to see if they converge on biologically relevant processes (immune, inflammation, etc) 







