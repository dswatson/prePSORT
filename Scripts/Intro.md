Transcriptomic Signatures of Treatment Response in Psoriasis
================
David Watson, Amy Foulkes & Mike Barnes

-   [Project Summary](#project-summary)
-   [Data Preparation](#data-preparation)
-   [Baseline Comparison](#baseline-comparison)

Project Summary
===============

Premise
-------

Between 2011 and 2012, participants were recuited in the North West of England to a study entitled "Pharmacogenomic signatures of treatment response in psoriasis" (Research Ethics UK 11/NW/0100) seeking to evaluate gene expression profile from patients with moderate to severe psoriasis.

The overall aim of this pilot project was to evaluate the feasibility of personalised treatments. Biologic therapies (now comprising 5 licensed effective injectable treatments) are expensive (approx £10,000 per patient per year) and have been associated with potentially serious adverse events.

The use of biologic therapies continues to increase. Psoriasis itself is an excellent candidate disease for personalisation of therapies; it is common, chronic and life-limiting with adverse effects on quality of life indices and high impact on healthcare systems. Inflammatory skin disease is an excellent choice for study, with ready access to disease tissue and objective validated scores of disease severity.

Clinical project overview
-------------------------

We prospectively recruited 10 participants, washed out from prior therapies, genetically homogenous (Caucasian with European ancestry to third generation) and collected detailed phenotype information. All patients had disease confirmed by a Consultant Dermatologist and had chronic plaque disease (the most common form, other forms of disease were excluded) with onset at or prior to the age of 40 (Type I or 'early onset' disease). These patients recevied treatment with the tumour necrosis factor inhibitor [Etanercept](https://www.nice.org.uk/guidance/ta199).

Extensive samples were taken, including blood, skin and urine at the 3 clinical time points - baseline (immediately prior to first injection, after washout of 2 weeks from prior therapy or at least 4 x 1/2 life of prior drug); 1 week after treatment; and 12 weeks after treatment. Skin samples were taken from controlled, photoprotected sites on the lower back or upper buttock, from both normal looking skin and from the edge of a plaque of psoriasis. Repeat biopsies were taken at a minimum distance of 2cm from prior biopsies. The severity score Psoriasis Area and Severity Index (PASI) was measured at each visit. All visits were conducted by one doctor and one specialist nurse. Participant journey through treatment was monitored, including adherence to therapy.

Analytic Goals
--------------

**Our aim is to evaluate gene expression profile at baseline, or early in therapy (week 1 data) and ask whether this can predict treatment response at week 12.**

This is fundamental for proof of principle, with a larger Consortium collaborative project Psoriasis Stratification to Optimise Relevant Therapies [PSORT](http://www.PSORT.org.uk) currently recruiting a larger cohort of patients commencing biologic therapy. Furthermore, we have further DNA SNP array, metabolomic and proteomic data for multi-omic data integration following completion of analysis of the RNA-seq data.

The following represents a first draft of our DESeq analysis pipeline for these data. While we ultimately intend to analyse each tissue type both within and across time points, this script is limited to a baseline comparison of gene expression profiles from lesional skin samples as an example to demonstrate our basic approach and raise certain issues we've encountered thus far.

``` r
# Load libraries
library(readr)
library(tximport)
library(ggplot2)
library(sva)
library(DESeq2)
library(IHW)
library(gridExtra)
library(dplyr)
```

Data Preparation
================

Loading
-------

Following RNA-seqencing of all samples, reads were pseudo-aligned using [kallisto](http://www.nature.com/nbt/journal/v34/n5/full/nbt.3519.html). The following code chunk presumes that your working directory contains:

-   the clinical factors file ('Clinical\_factors\_Pre\_PSORT\_Amy\_v6\_26\_7\_16.csv');
-   the isoform mapping file ('Transcripts.tsv'); and
-   a subdirectory with folders corresponding to the kallisto output for each sample ('LesionalSkin').

``` r
pheno <- read_csv('Clinical_factors_Pre_PSORT_Amy_v7_15_8_16.csv')
t2g <- read_tsv('Transcripts.tsv')
dir <- paste(getwd(), 'LesionalSkin', sep='/')
files <- file.path(dir, pheno$sample, 'abundance.tsv')
txi <- tximport(files, type='kallisto', tx2gene=t2g, reader=read_csv)
```

As we discussed at CSAMA, we have several key questions related to this script. We are primarily concerned with the proper categorisation of response, and finding an appropriate design for answering our clinical question using time course data.

Response
--------

Clinically and for the purposes of funding, we measure psoriasis severity using the Psoriasis Area and Severity Index ([PASI](https://www.ncbi.nlm.nih.gov/pubmed/?term=psoriasis+area+severity+index+PASI+Fredriksson)). This provides a score between 0 and 72, where a score of &gt; 10 signifies severe disease. Whilst this score has numerous problems including high inter-individual variability in measurement and poor sensitivity at the lower end of the scale, it is the most widely used measure in psoriasis research and is adopted through clinical guidelines UK & Europe wide for the eligibility for biologic therapies and the assessment of response. It is non-linear. Most commonly, a 75% improvement from baseline score signifies response to treatment, a &lt; 50% improvement signifies non-response, and anything in between is considered an intermediate response. Newer successful therapies are expected to reach PASI 90 as a benchmark for later phases of clinical trials.

During our study, I also recorded additional response measurements including the Physician Global Assessment ([PGA](https://www.ncbi.nlm.nih.gov/pubmed/?term=psoriasis+global+assessment+langley+2004)) and the Dermatology Quality of Life Index ([DLQI](http://www.ncbi.nlm.nih.gov/pubmed/8033378)). The PGA is rarely used in clinical practice but is commonly used in clinical trials. The DLQI is used routinely in clinical practice, but is not specific to psoriasis. It has many flaws, not least that it leads to many NAs, for instance when subjects are embarrassed to answer questions relating to relationships, exercise etc. It is used by funders in the UK as a surrogate marker of response; patients with a less than adequate response (by PASI) to a biologic but with a 5 point improvement in the DLQI can continue their treatment.

### Delta\_PASI Distribution

The distribution of our cohort's percent change in PASI scores between baseline and week 12 (a variable labeled 'Delta\_PASI' in our clinical factors file) suggests that a binary response measure would be completely arbitrary.

``` r
# Build pasi data frame
pasi <- data.frame(Subject    = unique(pheno$subject),
                   Baseline   = unique(pheno$PASI_A), 
                   TwelveWeek = unique(pheno$PASI_C),
                   DeltaPASI  = unique(pheno$Delta_PASI))

# Plot rank vs. DeltaPASI
ggplot(pasi, aes(rank(DeltaPASI), DeltaPASI, label=Subject)) + geom_text() +
  labs(title='Delta_PASI Distribution', x='Rank', y='Delta_PASI') + 
  xlim(-0.25, 10.25) + theme_bw()
```

<p align='center'>
<img src="Intro_files/figure-markdown_github/deltpas-1.png" style="display: block; margin: auto;" />
</p>

This plot, in conjunction with the fact that a single clinician recorded every PASI score in our dataset, led us to conclude that we should analyse response on a continuous scale.

Based on our discussions at CSAMA, we considered transforming our response data to give it a more symmetrical distribution. However, one patient's disease worsened during treatment, resulting in a negative Delta\_PASI score. We considered adding a small value (e.g. ~ .05) to each patient's Delta\_PASI score in order to facilitate a log or Box-Cox transformation, but were unsure if this would be appropriate. For now, we have decided to measure response via simple percentage change.

One simple way to visualise patient response is by plotting baseline vs. 12 week PASI scores. A straight diagonal line across the plot represents total stasis (i.e., 0 change in PASI) while a flat line with slope and intercept 0 would represent 100% improvement in PASI.

``` r
# Build zero-intercept model
mod <- lm(TwelveWeek ~ 0 + Baseline, data=pasi)

# Plot baseline vs. 12 week PASI scores
ggplot(pasi, aes(Baseline, TwelveWeek, label=Subject)) + geom_text() + 
  geom_abline(slope=1, intercept=0) + xlim(0, 40) + ylim(0, 40) + 
  geom_abline(slope=mod$coef, intercept=0, colour='blue') + 
  labs(title='Change in PASI', x='Baseline PASI', y='12 Week PASI') + theme_bw()
```

<p align='center'>
<img src="Intro_files/figure-markdown_github/deltpas2-1.png" style="display: block; margin: auto;" />
</p>

This plot indicates that patient improvement was roughly linear in our study, albeit with a notable outlier in patient AFOX09. The blue regression line has a slope of about .31, suggesting a least squares estimate of 69% expected improvement in PASI. (This number is closer to 72% if subject AFOX09 is removed.)

DeltaPASI does not appear to depend on baseline PASI scores. We determine this both through informal plotting and statistical testing.

``` r
# Plot baseline vs. DeltaPASI
ggplot(pasi, aes(Baseline, DeltaPASI, label=Subject)) + geom_text() + 
  geom_smooth(method='lm', se=FALSE) + 
  labs(title='Baseline vs. Delta_PASI', x='Baseline PASI', y='Delta_PASI') + theme_bw()
```

<p align='center'>
<img src="Intro_files/figure-markdown_github/deltpas3-1.png" style="display: block; margin: auto;" />
</p>

``` r
# Correlation test
cor.test(pasi$DeltaPASI, pasi$Baseline)
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  pasi$DeltaPASI and pasi$Baseline
    ## t = 0.76564, df = 8, p-value = 0.4659
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.4408655  0.7650534
    ## sample estimates:
    ##       cor 
    ## 0.2612898

In fact, what little semblance of a relationship the two variables may appear to have in the plot completely disappears when our outlying observation is removed.

``` r
# Correlation test, take 2
cor.test(pasi$DeltaPASI[-8], pasi$Baseline[-8])
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  pasi$DeltaPASI[-8] and pasi$Baseline[-8]
    ## t = -0.020407, df = 7, p-value = 0.9843
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.6684108  0.6597884
    ## sample estimates:
    ##          cor 
    ## -0.007712975

### Time Course Analysis

Clinically using a baseline within-patient measurement with repeated during-treatment observations makes sense. We understand for our analysis we need to interact time and response. The interaction of a factor *x* and a continuous variable *y* requires the creation of *k* new variables, one for each level of *x*.

``` r
pheno <- pheno %>%
  mutate(wk0.Delta_PASI  = ifelse(time == 'wk0',  Delta_PASI, 0),
         wk1.Delta_PASI  = ifelse(time == 'wk1',  Delta_PASI, 0),
         wk12.Delta_PASI = ifelse(time == 'wk12', Delta_PASI, 0))
```

This solution produces some sensible results, but strikes us as somewhat sloppy and ad hoc. We are currently unaware of any better workaround than this, however we remain very open to any and all suggestions.

Surrogate Variable Analysis
---------------------------

We incorporate surrogate variables to account for unknown batch variation. Because our study involves repeat observations, we include 'subject' as a covariate in our null and full design matrices.

``` r
# Normalise and filter data before running SVA
dds <- DESeqDataSetFromTximport(txi, colData=pheno, design= ~ subject + wk0.Delta_PASI)
dds <- estimateSizeFactors(dds)
dat <- counts(dds, normalized=TRUE)
dat <- dat[rowVars(dat) > 1, ]  # This was a mean filter in the previous script

# Run SVA
mod <- model.matrix(~ subject + wk0.Delta_PASI, data=pheno)
mod0 <- model.matrix(~ subject, data=pheno)
svobj <- svaseq(dat, mod, mod0)
```

    ## Number of significant surrogate variables is:  5 
    ## Iteration (out of 5 ):1  2  3  4  5

``` r
# Add surrogate variables to pheno
pheno <- cbind(pheno, svobj$sv)
colnames(pheno)[51:55] <- paste0('SV', 1:5)
```

Baseline Comparison
===================

Our data is now ready to be run through the DESeq2 pipeline.

``` r
dds <- DESeqDataSetFromTximport(txi, colData=pheno, 
                                design= ~ subject + SV1 + SV2 + SV3 + SV4 + SV5 + wk0.Delta_PASI)
dds <- DESeq(dds)
res <- data.frame(results(dds, filterFun=ihw))
res <- res[order(res$pvalue), ]
```

Let's take a look at the genes declared significantly differentially expressed at 5% FDR.

``` r
na.omit(res[res$padj < .05, ])
```

    ##          baseMean log2FoldChange      lfcSE      stat       pvalue
    ## CCL13    680.3063     -0.7092988 0.12690036 -5.589415 2.278354e-08
    ## ADNP2   1317.1824      0.4962049 0.09365880  5.298006 1.170742e-07
    ## WBP11   1151.4050      0.3716361 0.07350310  5.056060 4.280067e-07
    ## ESAM     304.5407      0.5234042 0.10403981  5.030807 4.884194e-07
    ## IL19     208.3144      0.5290386 0.10583212  4.998847 5.767405e-07
    ## IGLC1    233.7632      0.4839601 0.09975356  4.851557 1.224958e-06
    ## DAPL1   1399.8672     -0.5094046 0.11042474 -4.613138 3.966340e-06
    ## DAPP1   1229.9509     -0.4172611 0.09218295 -4.526445 5.998419e-06
    ## EREG    3401.8610     -0.5191387 0.11611268 -4.470990 7.785828e-06
    ## CTSH    2435.9297     -0.4217540 0.09702756 -4.346745 1.381728e-05
    ## PCDHGC3 1403.9006      0.3575466 0.08626049  4.144964 3.398670e-05
    ## OVOL1   2409.7039     -0.4679856 0.11347592 -4.124096 3.721933e-05
    ## CISD2    848.1727     -0.3360152 0.08275952 -4.060139 4.904347e-05
    ## MTMR1   1851.8285     -0.3319951 0.08396203 -3.954110 7.682013e-05
    ##                 padj    weight
    ## CCL13   0.0002691379 4.4658251
    ## ADNP2   0.0006125959 5.0409525
    ## WBP11   0.0011431164 6.5840664
    ## ESAM    0.0058716712 0.8776404
    ## IL19    0.0015889655 4.7869775
    ## IGLC1   0.0092038582 0.8776404
    ## DAPL1   0.0067541227 5.1632731
    ## DAPP1   0.0072690901 6.2189074
    ## EREG    0.0094124153 4.8486030
    ## CTSH    0.0117209788 6.2189074
    ## PCDHGC3 0.0262094335 6.2189074
    ## OVOL1   0.0316896631 5.1632731
    ## CISD2   0.0444888398 4.4734460
    ## MTMR1   0.0446288526 6.4861455

We find here that using a variance rather than a mean filter during SVA has pushed a few more genes below the 5% FDR threshold. This is line with the findings of [Bourgon et al. (2010)](http://www.pnas.org/content/107/21/9546.long). Independent hypothesis weighting has also helped increase statistical power ([Ignatiadis et al., 2016](http://biorxiv.org/content/early/2016/03/30/034330)).

To visualise the relationship between expression and Delta\_PASI, we first transform the counts. Below are scatterplots for the top four genes.

``` r
# Transform counts
rld <- rlog(dds)

# Prep plot
rld <- as.data.frame(t(assay(rld)))
top <- data.frame(CCL13=rld$CCL13, ADNP2=rld$ADNP2, WBP11=rld$WBP11, ESAM=rld$ESAM,
                  DeltaPASI=pheno$Delta_PASI, Time=pheno$time) 

# Create plots
g1 <- ggplot(top, aes(DeltaPASI, CCL13, colour=Time)) + geom_point() + 
  labs(title='CCL13', y='Expression') + theme_bw()
g2 <- ggplot(top, aes(DeltaPASI, ADNP2, colour=Time)) + geom_point() + 
  labs(title='ADNP2', y='Expression') + theme_bw()
g3 <- ggplot(top, aes(DeltaPASI, WBP11, colour=Time)) + geom_point() + 
  labs(title='WBP11', y='Expression') + theme_bw()
g4 <- ggplot(top, aes(DeltaPASI, ESAM, colour=Time)) + geom_point() + 
  labs(title='ESAM', y='Expression') + theme_bw()

grid.arrange(g1, g2, g3, g4, nrow=2, ncol=2)
```
<p align='center'>
<img src="Intro_files/figure-markdown_github/expr-1.png" style="display: block; margin: auto;" />
</p>