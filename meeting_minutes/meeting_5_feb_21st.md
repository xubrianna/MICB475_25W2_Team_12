# Questions

1. Is it overkill to conduct 3 analyses for Aim 2: To identify differentially abundant and core microbial taxa associated with endometriosis among CPP patients, compared to CPP patients without endometriosis and healthy controls
  - Core Microbiome + DESeq2 + ISA? which one is best 
2. For the PICRUSt2 analysis it gives us 3 outputs KO, EC, MetaCyc, do we do all 3 or pick 1 for downstream analysis? 
- individual functions or grouped functions for microbes?
3. Should we optimize the # features or samples for rarefaction?
  - We thought 43,000 since it kept the most ASVs but we lost a lot of samples 
- We decided to choose a sample depth of 43,000 to maximize the number of features (62.2%) in samples for downstream analysis. As a result, we lost 39 samples from our initial pool of 146 which is justified because we have chosen to optimize features.  
4. Ask about rationale and whether it seems justified and biologically valid
5. ask about the healthy controls - are we comparing CPP to healthy first, and then CPP-endo to CPP? is healthy control the baseline? 

discuss overall points/flow: 

- CPP is a broad overarching disease that is a symptom that can stem from many causes
- it is hard to identify what the underlying cause is
- endometriosis is one of the causes, and since in itself is hard to diagnose, CPP can be an effective biomarker
- since endometriosis has links to microbiome, and microbes can cause pain, being able to establish these dysbiostic shifts in the context of endo can help us diagnose endo as the root cause more rapidly instead of waiting 6+ years for better treatment
- not too much known about the particular bodysites, but would be interesting to be able to spot where the dysbiosis is occuring so that physicans can know where to swab from for diagnosis? 

- hypothesis: CPP-endo will have microbial shifts compared to CPP alone, and that it will have upregulation of pathogenic/inflammatory/pain taxa and pathways and reduction of SCFA/essnetial aa pathways (based on previous research)
- main interest: 

- we are first looking at taxonomy, then DEGS and core microbiome to see which taxa are unique to each condition, then trying to establish the functional pathways that could be contributing to pain, then potentially if we find unique taxa trying to train a predictive classifier to differentiate endo-CPP vs CPP alone



