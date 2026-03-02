# Meeting 6 — March 2nd

## Aim 1 update: Alpha and Beta Diversity

### Graphs from `R/01-aim1_alpha_diversity.R`

**Rarefaction curve (depth selection = 43,000):**

![Rarefaction curve](../results/aim1/00-rarefaction_curve.png)

**Faith's PD by disease group (faceted by sample site):**

![Faith PD boxplot](../results/aim1/01_faith_PD_boxplot.png)

### Graph from `R/02-aim1_beta_diversity.R`

**Unweighted UniFrac PCoA (faceted by sample site):**

![Unweighted UniFrac PCoA](../results/aim1/02-pcoa.png)

## Summary of findings

- After rarefaction to 43,000 reads/sample, **107 samples** remained for downstream diversity analysis.
- **Alpha diversity (Faith's PD)** did not differ significantly across disease groups within either body site:
	- Vaginal: Kruskal-Wallis $p = 0.856$
	- Rectal: Kruskal-Wallis $p = 0.529$
	- Dunn's post-hoc pairwise tests were all non-significant after BH correction.
- **Beta diversity (Unweighted UniFrac, PERMANOVA)** showed:
	- A significant **body-site effect** (env_medium): $R^2 = 0.212$, $p = 0.001$
	- No significant **disease-group effect** after controlling for site: $R^2 = 0.0146$, $p = 0.455$
	- Within-site disease effects were also non-significant (rectal $p = 0.709$, vaginal $p = 0.457$).

## Interpretation for discussion
- For next steps, disease-associated signals may be more detectable in taxon-level differential abundance** and functional pathway analyses (Aims 2 and 3) rather than broad community diversity alone.

## Questions 
- editing the manuscript - we had few changes, so we will submit?


## Edits from the meeting 
- find a theme to express plot (theme_classic) (theme(___)) -> set a constant theme 

- colour brewer 
- 95% elipses 
- x axis titles (variance) 
- combine the plots 

- try violin plot 

- axis limits higher and lower than extreme limits (make 5 the 0, and have the 20 limit - on the beta plot too 

- geom point vs geom jitter 

- don’t need legend for alpha diversity 
- its okay if not significant - signal can be lost when looking at a top down view (diluted out) 

talking about presentation: 
- have the data for all of the aims by the second meeting before presentation - march 23rd -> start to put together as a story and assemble into a presentation 
- prepare slides for the 30th - go together with Hans (typos, story, etc) 
- presenting another teams 

