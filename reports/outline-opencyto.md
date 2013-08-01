# Automated Flow-Cytometry Data Analysis with OpenCyto
## John A. Ramey

* Hypothesis: Automated gating can perform as well or better than manual gating to identify responding T-cell subpopulations in vaccine clinical trial data.
* Introduction
    * Goal: Use separation of meaningful cell subpopulations identified via automated gating (pre- and post-vaccination) to identify patients with a vaccine response.
    * We use two intracellular cytokine staining (ICS) data sets
        1. HIV Vaccine Trials Network: HVTN 065 Data Set from Goepfert et al. (2011)
        2. Immune Tolerance Network: ITN507ST Data Set from Newell et al. (2010)
    * Disadvantages of manual gating
        * Time-consuming
        * Inherently subjective
        * Highly variable gate placement from person to person if:
            * An experiment is not well-controlled
            * A marker is not well-resolved
    * Emphasize the pitfalls of manual gating; automation is preferred.
    * OpenCyto yields accurate reproduction of manual gating schemes in an automated manner
    * OpenCyto incorporates prior knowledge through a Bayesian framework.
    * OpenCyto attains fast, robust, accurate gating of rare cell populations
  * OpenCyto
      * Pipeline based on a gating hierarchy defined in a CSV file
      * Data-derived gates for each sample
      * Robust gating with Bayesian mixture models via **flowClust 3.0**
      * Priors are marker-specific, data-driven, and can incorporate expert knowledge
      * Gating parameters can be optimized to to discriminate between patient cohorts (e.g., vaccination status)
      * Custom gating algorithms can be utilized via a plugin system
      * Show that OpenCyto can be used to replicate manual gating with respect to low variability and bias compared to the manual gating statistics.
  * Classification
    * Show that we can improve discrimination by using more features and optimizing the gating, and making this a classification problem (when baseline is available).
    * Rapidly prototype different gating thresholds for the cytokines (i.e., you could design the pipeline for the ENV stimulation and run it on the GAG stimulation for each data set to show that it is robust when data is standardized)
    * Describe how our classifier is constructed
      * We extract all Boolean subsets with associated proportions as features
      * Briefly provide example Boolean subset, similar to FlowCAP 3 talk: (CD4) IL2+ and !IFNg+ and TNFa+
      * We then utilize a LASSO-based classifier using the **glmnet** R package
      * Mention briefly that the shrinkage parameter selected via cross-validation
      * Mention also that **glmnet** employs a variable selection via $L_1$ regularization
      * Rank variables selected using **covTest** package
    * Discuss data sets and results
      * Data sets
        * Data set \#1: HVTN 065
            * 79 patients -- 67 treated patients and 12 placebos
            * Time points of interest: 2 (pre-vaccine) and 12 (post-vaccine)
            * 10 FCS files per patient (5 per time point)
            * 4 (5?) Antigens
            * 3 Cytokines of Interest
            * Goal \#1: Classify the vaccination status of the patient cohorts via the elastic net classifier (the **glmnet** R package)
            * Goal \#2: Identify antigen-specific T-cells responding to the vaccine using the markers selected
        * Data set \#2: ITN (Newell)
        * Demonstrate two gates
            1. Cytokine gate
            2. Transitional gate
      * Results
        * Emphasize that OpenCyto yields similar results to the manual gating
        * Discuss the features selected by **glmnet**
        * Figures:
          * Figure 1: Output from flowClust that demonstrates the fitted mixture model
          * Figure 2: Comparison of automated and manual gates
          * Figure 3: Gated proportions of stimulation groups by features selected for each training patient
          * Figure 4: Comparison of manual gating and automated gating statistics. (show CVs)
        * Tables:
          * Table 1: Response rates by treatment group for each data set (manual and automated).
* Discussion
