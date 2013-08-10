# Automated Gating of Flow-Cytometry Data with the OpenCyto Framework

**Authors**: John A. Ramey, Greg Finak, Mike Jiang, Jafar Taghiyar, Stephen de
   Rosa, Ryan Brinkman, and Raphael Gottardo

## Introduction

Advancements in flow cytometry (FCM) technologies and instrumentation have
enabled rapid quantification of multidimensional attributes for millions of
individual cells to assess cellular heterogeneity and to identify meaningful,
homogeneous cellular subpopulations through a process called gating. However,
the analysis of the resulting large, high-dimensional data sets typically
involves a time-consuming sequential manual gating strategy that is inherently
subjective as the gating depends on the individual performing the gating. This
subjectivity can yield highly variable gate placement from person to person if
an experiment is not well-controlled or if a marker is not well-resolved as
there may be ambiguity as to where a cell-population boundary should be placed
\citep{Aghaeepour:2013dg; Maecker:2012gl; Bashashati:2009em;
Maecker:2005gm}. Furthermore, the manual construction of subjective gates
greatly limits reproducibility and expedited analysis of high-throughput FCM
data \citep{Pyne:2009ec; Hahne:2009hl}. Alternatively, automated, data-driven
algorithms are necessary to expedite the gating of FCM data and to remove the
subjectivity intrinsic to manual gating. This is particularly important in
clinical trials where assays must be extremely well controlled in order to
generate data that is comparable over time.

(**REWORD**) In recent years a large number of automated computational
algorithms have been proposed \citep{Ge:2012hr; Aghaeepour:2010fv;
Naumann:2010jn; Pyne:2009ec; Lo:2008it; Chan:2008gq} and have been shown to
identify cellular populations both reliably and accurately
\citep{Aghaeepour:2013dg}. Despite the level of maturity that data-driven
algorithms have achieved, the software packages available for these algorithms
are often too complex for the average user \citep{McNeil:2013du}, provide
individual solutions for an experiment-specific analysis, and are generally
inadequate in reliably and robustly translating of raw data from a large number
of samples and experiments into quantified accessible results
\citep{Maecker:2010fg}. \cite{Maecker:2012gl} have argued that the resources to
enable a standardized application of all appropriate solutions to most common
FCM immunophenotyping challenges is missing. \cite{McNeil:2013du} and
\cite{Maecker:2012gl} have called for standard, easy-to-use software that can
produce unbiased, reproducible of FCM data, is amenable to differences in
experimental conditions, less computationally intense, and readily available to
the entire flow-cytometry community. Furthermore, \cite{McNeil:2013du} have
called for software and algorithms to identify more precisely and quantify rare
antigen-specific cells of interest. Additionally, \citep{Aghaeepour:2010fv} have
argued that despite its success in recapitulating manual analysis, FCM not
reached its full potential due to the lack of an automated analysis platform to
assist high-throughput data generation.

We have developed the OpenCyto framework, a collection of well-integrated
open-source R packages that delivers robust, reproducible, and data-driven
gating in an automated pipeline by incorporating expert-elicited and data-driven
prior knowledge within a Bayesian model. For a given gating hierarchy, OpenCyto
promotes relatively fast and exhaustive gating that is interpretable in the
context of standard hierarchical, two-dimensional projections of cell
populations, which immunologists and other analysts are used to seeing. Our
automated gating approach allows gating thresholds to be fine-tuned to optimize
detection of informative cell populations in order to discriminate between
patient cohorts based on objective external criteria such as vaccination
status.

The OpenCyto framework provides automated, data-driven gating of
high-dimensional FCM data sets quickly, removing the time-consuming task of
manual gating. By incorporating expert-elicited and data-driven prior knowledge,
OpenCyto attains accurate gating of cell populations, including rare
populations, while controlling variability relative to manual gating, thereby
overcoming the subjectivity in manual gating. Finally, OpenCyto is clearly
valuable in its ability to construct data-driven gates and reproducibly identify
associated biomarkers to distinguish vaccination status within a cohort in
clinical trial data.

The OpenCyto framework includes several built-in flexible and robust gating
algorithms that are applicable to a wide variety of scenarios. Furthermore, we
have included **flowClust 3.0**, which allows expert and data-driven prior
knowledge into the gating process via a Bayesian framework. The OpenCyto
framework, however, is not limited to the gating algorithms included. The
framework allows for custom gating algorithms to be included using a plugin
system. Bioinformaticians will find that custom gating algorithms are easy to
incorporate by simplify inputting and outputting standard **flowCore**
objects. For non-bioinformaticians, OpenCyto is easy to code: OpenCyto can be
run with a small number of lines of code. For standard routines, a user need
only define a gating hierarchy in a CSV template file.

Using the OpenCyto framework, we demonstrate that automated gating can perform
as well or better than manual gating to identify responding T-cell
subpopulations in vaccine clinical trial data. We demonstrate that the OpenCyto
framework can recapitulate manual-gating efforts obtained on an intracellular
cytokine staining (ICS) data set from the HIV Vaccine Trials Network
(HVTN). Specifically, we use the HVTN 065 Data Set from
\cite{Goepfert:2011ci}. From the resulting gates, we can apply multivariate
analysis using subsets of cellular subpopulations. For the manual gating as well
as OpenCyto, we calculated the coefficients of variation of cellular population
proportions across the samples and found that the variability from OpenCyto is
well within the range of that of manual gating, even for rare cellular
subpopulations. Furthermore, based on the gates constructed by OpenCyto, we were
able to discriminate accurately the vaccination status of the patient cohorts as
well as to identify the antigen-specific T-cells responding to the
vaccine. Furthermore, we were able to utilize the separation of meaningful cell
subpopulations identified via automated gating (pre- and post-vaccination) to
identify patients with a vaccine response. Previously, we have applied the
OpenCyto framework to FlowCAP data \cite{Aghaeepour:2013dg} and obtained
excellent results. These data sets are being considered elsewhere.

Additionally, \cite{Bocsi:2008cv} have argued that full automation is often
untenable because simply put, there is no free lunch. However effective in the
majority of cases, an automated method can provide a poor gate in some
instances, such as in identifying rare cytokine subpopulations. Hence, it is
imperative that diagnostic tools be readily available to identify quickly such
failures. Within the OpenCyto framework, we provide several methods to determine
statistically whether a gating algorithm is poor based on either the proportion
of cells gates or based on the position of the gate. These sanity checks are
especially useful in subpopulations of interest where outliers or data artifacts
within a single sample yield a poor automated gate. Effectively, to correct
problem gates, we identify similar samples from which we obtain a gate or
additional data from which we can reapply the gating algorithm.

## The OpenCyto Framework

The OpenCyto framework is a system developed in R for automatically gating
flow-cytometry data.

* Pipeline based on a gating hierarchy defined in a CSV file
    * Simple
    * Reproducible
* Workspaces can be loaded from flowJo (CITE: TODO) via the flowWorkspace
  package (CITE: TODO).
* Data-derived gates for each sample
* Very little R code is necessary for standard analyses. Simply specify a gating
  strategy in the CSV file.
* However, we have created a plugin system for power users to incoporate custom
  gating algorithms.
* Discuss parallel processing
    * Cite \cite{McNeil:2013du}'s statement that algorithms and software need to
      be less computationally intense
    * Parallel processing greatly improves computational performance
    * Several samples can be gated simultaneously to improve computational
      performance
* We have included several gating algorithms for starters
    * The algorithms range from fast, simple algorithms for finding gating
      thresholds for easy-to-gate upstream markers (e.g., CD3)
    * Also, more sophisticated gating algorithms to identify rare,
      antigen-specific cellular subpopulations of interest, e.g., cytokines.
    * We have also provided more sophisticated, robust algorithms for
      identifying higher-dimensional populations of interest, including
      lymphocyte populations
        * Robust gating with Bayesian mixture models via **flowClust 3.0**
        * Useful because priors are marker-specific, data-driven, and can
          incorporate expert knowledge
        * Also useful for automatic detection when several populations
          (clusters) of interest are present
* Gating parameters can be optimized to to discriminate between patient cohorts
  (e.g., vaccination status)
* Standard Boolean gating can be incorporated in the CSV file
    * Stress Polyfunctional gating
* We can quickly process large collections of flow-cytometry data (i.e., FCS files)
    * Also, memory efficient
    * Infrastructure has C++ backend based on HDF5, netCDF, and boost libraries.
    * Also, uses flowWorkspace, flowCore, and ncdfFlow R packages.

The OpenCyto framework allows for custom gating algorithms to be incorporated
with ease. Custom gating algorithms can be utilized via a plugin system. The
user-defined function need only worry about receiving a sample or a subset of
samples from which the function will construct data-derived gates. Furthermore,
the user does not have to worry about the infrastructure for gating a large
collection of samples as OpenCyto takes care of this. Hence, the user can
rapidly prototype different gating approaches and examine the efficacy of a
several candidate tuning parameters.

## Classification Study

* Show that we can improve discrimination by using more features and optimizing
  the gating, and making this a classification problem (when baseline is
  available).
* Rapidly prototype different gating thresholds for the cytokines (i.e., you
  could design the pipeline for the ENV stimulation and run it on the GAG
  stimulation for each data set to show that it is robust when data is
  standardized)

### HVTN065 Data Set

The HVTN065 data set contains 79 patients, 67 of whom were treated patients. The
remaining 12 patients were placed in the placebo cohort. We examined the time
points 2 (pre-vaccine) and 12 (post-vaccine) for each patient. For each time
point, there were three stimulated samples (GAG-1-PTEG, ENV-1-PTEG, and
POL-1-PTEG), two negative controls, and one positive control. Hence, there were
12 FCS files per patient (6 per time point). We applied OpenCyto to gate each of
the samples.

We considered only those FCS files having at least 100,000 cells. This reduced
the number of patients to 73.

For each T-cell subset CD4 and CD8, we extracted all Boolean subsets of the
three cytokines of interest with associated proportions as multivariate features
vectors.

1. TNFa+IFNg+IL2+
2. TNFa+IFNg+IL2-
3. TNFa+IFNg-IL2+
4. TNFa+IFNg-IL2-
5. TNFa-IFNg+IL2+
6. TNFa-IFNg+IL2-
7. TNFa-IFNg-IL2+

We explicity ignored TNFa-IFNg-IL2- because it is redundant given the other
seven features.

**TODO**: Add MFI's to model to see if it helps. (This should improve paper's
  journal quality)

Our goal was to classify the vaccination status of the patient cohorts from the
proportions of the polyfunctional populations of the three cytokines.
Additionally, we aimed to identify antigen-specific T-cells responding to the
vaccine.

Here, we describe our classification study. First, we partition the HVTN065
patients by their treatment status into a treatment group and placebo group. Of
the patients in the treatment group, we randomly partition 60% of the them into
a training data set and the remaining 40% of the patients into a test data set.
We utilize the **glmnet** package using the elastic net with **alpha = 0.5** to
build a classifier from the population proportions for the markers and the
polyfunctional gates obtained using the OpenCyto package. Next, because there
are two visits (i.e., pre- and post-vaccine) for each patient, we pair the
visits in the test data set by patient. For each patient-visit pairing, we
classify the two samples and calculate the difference in their classification
probabilities. Let d = Pr(sample 1 from subject 1 = post-vaccine) - Pr(sample 2
from subject 1 = post-vaccine).  For a given probability threshold, if $d >
threshold$, then we classify sample 1 as post-vaccine and sample 2 as
pre-vaccine. Otherwise, if $d < threshold$, we classify sample 1 as pre-vaccine
and sample 2 as post-vaccine. We calculate the classification accuracy as the
number of correctly classified patients. In the same manner we calculate the
classification accuracy using the placebo patients as a separate test data set.

We selected relevant features using a LASSO-based elastic-net classifier.
Specifically, we used the **glmnet** R package, which employs a variable
selection via $L_1$ regularization. The shrinkage parameter was selected via
10-fold cross-validation. We then ranked the variables selected by $p$-values
obtained from the **covTest** R package.

**TODO**: We should say that 1 stimulation group induces the following subset of
          features ____, while another stimulation group induces the following
          subset of features ____.

Per Greg: 'We also want to do this paired, using the difference in
classification probabilities for two samples from the same subject. i.e.  d =
Pr(sample 1 from subject 1 = post-vaccine) - Pr(sample 2 from subject 1 =
post-vaccine). If d > threshold, then classify sample 1 as post-vaccine and
sample 2 as pre-vaccine, otherwise if d < threshold classify sample 1 as
pre-vaccine and sample 2 as post-vaccine, otherwise mark them as
unclassifiable.'

Next, we calculate ROC curves assuming all vaccinees are true positive and the
placebos are false positive. For each PTID, we compute the absolute value of the
difference in classification probabilties for visits 2 and 12 and then order by
the differences.

Here, we provide the markers that were selected by **glmnet** for each
stimulation group.

### Classification with IFNg+ | IL2+

Here, we compare the automated gates constructed using OpenCyto with the
manually constructed gates. For both cases, we apply logistic regression to
classify patient visits utilizing a single feature, namely IFNg+ | IL2+.

### Automated Gates using OpenCyto

Next, we calculate ROC curves assuming all vaccinees are true positive and the
placebos are false positive. For each PTID, we compute the absolute value of the
difference in classification probabilties for visits 2 and 12 and then order by
the differences.


* Discuss data set and results
    * Data set: HVTN 065
        * **TODO**
            * Update summaries.
            * Question RG wants us to answer: Do features selected on HVTN065
               correlate with antibody response?
            * My clarifying questions:
              * Where are antibody responses recorded? (In CSV file?)
              * How to measure the correlation to answer RG's questions?
  * Results
    * Emphasize that OpenCyto yields similar results to the manual gating
    * Discuss the features selected by **glmnet**
    * Figures:
      * Figure 1: Comparison of automated and manual gates
      * Figure 2: Gated proportions of stimulation groups by features selected
        for each training patient
      * Figure 3: Comparison of manual gating and automated gating statistics. (show CVs)
      * **TODO**: Greg mentioned that we should include boxplots of the
        proportions (the reps are the subjects) and stratify by cohorts. This
        information is available for HVTN065 in the ./data/*.csv file (I
        think). Then facet by the markers (features selected?).
    * Tables:
      * Table 1: Response rates by treatment group for each data set (manual and automated).
    * **TODO**: Scrap emphasis on flowClust. If we decide to include a figure
        (e.g., output of fitted flowClust mixture model) it should go into the
        Supplementary.

## Discussion

## Bibliography

\bibliographystyle{plainnat}
\bibliography{opencyto}
