Stringer_2014_sipoglitazar <- function() {
  description <- paste(
    "Simultaneous, cascading indirect-response PK/PD model for fasting plasma",
    "glucose (FPG) and glycosylated hemoglobin (HbA1c) in 780 drug-naive Type 2",
    "diabetes mellitus patients treated for 13 weeks with the PPAR alpha/delta/gamma",
    "agonist sipoglitazar, rosiglitazone 8 mg once daily, or placebo (Stringer 2014).",
    "FPG is a zero-order production / first-order loss turnover pool (KinG, KoutG);",
    "HbA1c is a secondary cascade driven by a power function of FPG (FPG^gamma,",
    "gamma = 0.71) with first-order production and degradation (KinH, KoutH).",
    "The sipoglitazar effect stimulates KoutG through an Emax function of",
    "steady-state exposure, AUC0-24h = total daily dose / CL, where CL is fixed",
    "per UGT2B15 genotype (5.04, 3.35 and 1.53 L/h for UGT2B15*1/*1, *1/*2 and",
    "*2/*2 respectively; Table 1) so that genotype-driven exposure differences",
    "propagate into the glycemic response. Rosiglitazone enters as a fixed 28%",
    "stimulatory step on KoutG (ROTE). Two study-conduct 'lifestyle' effects are",
    "carried: LEFPG, a mean-zero additive random effect on KinG capturing",
    "per-subject diet / exercise improvement or loss of glycemic control, and",
    "LEHB, a step inhibition of KinH that is larger in the placebo arm (3.7%)",
    "than in the sipoglitazar arms (2.0%) and absent on rosiglitazone. This is a",
    "PD-only model: no drug is dosed through the event table -- exposure is",
    "supplied through the DOSE_SIPOGLITAZAR_MGD and UGT2B15 genotype covariates.",
    "No demographic covariate (age, sex, weight, duration of disease) was retained."
  )
  reference <- paste(
    "Stringer F, DeJongh J, Scott G, Danhof M.",
    "A Model-Based Approach to Analyze the Influence of UGT2B15 Polymorphism",
    "Driven Pharmacokinetic Differences on the Pharmacodynamic Response of the",
    "PPAR Agonist Sipoglitazar.",
    "J Clin Pharmacol. 2014;54(4):453-461. doi:10.1002/jcph.227.",
    "Per-genotype sipoglitazar clearance values reproduced in Table 1 originate",
    "from the companion population PK analysis, Stringer F, Ploeger BA, DeJongh J,",
    "et al. J Clin Pharmacol. 2013;53(3):256-263.",
    sep = " "
  )
  vignette <- "Stringer_2014_sipoglitazar"
  units <- list(
    time          = "day",
    dosing        = "mg sipoglitazar (oral total daily dose, carried as the DOSE_SIPOGLITAZAR_MGD covariate; the model has no event-table dose records)",
    concentration = "FPG in mmol/L; HbA1c in % (NGSP); AUC0-24h at steady state in mg*day/L"
  )

  paper_specific_compartments <- c("hba1c")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Both states are turnover pools carried directly on the
  # measured (concentration / percentage) scale rather than as amounts, which
  # is the standard indirect-response parameterisation used by the source paper
  # (Equations 3 and 4).
  compartmentData <- list(
    glucose = list(analyte = "glucose", units = "mmol/L", specimen = "plasma", verified = TRUE),
    hba1c   = list(analyte = "HbA1c",   units = "%",      specimen = "blood cell", verified = TRUE)
  )

  covariateData <- list(
    TRT = list(
      description        = "Per-subject treatment-cohort indicator. 0 = placebo; 1 = sipoglitazar at any regimen OTHER than 32 mg twice daily (8 mg QD, 16 mg QD, 16 mg BID, 32 mg QD, 64 mg QD); 2 = sipoglitazar 32 mg twice daily; 3 = rosiglitazone 8 mg once daily.",
      units              = "(categorical / integer-coded)",
      type               = "categorical",
      reference_category = "0 (placebo)",
      notes              = "Stringer 2014 Methods 'Subjects and Data Collection': 'Patients were treated with sipoglitazar, rosiglitazone, or placebo: sipoglitazar 8 mg once daily (QD), 16 mg QD, 16 mg twice daily (BID), 32 mg QD, 32 mg BID, and 64 mg QD, placebo or rosiglitazone 8 mg QD.' Disposition: sipoglitazar n = 572, rosiglitazone n = 72, placebo n = 136 (total 780). The cohort indicator carries three distinct pieces of model structure. (a) Level 2 is split out from the other sipoglitazar regimens because the paper estimated a SEPARATE FPG baseline for that arm: Methods 'A lower FPG baseline value was observed in the 32 mg BID sipoglitazar group compared to all other treatment groups. The addition of a separate FPG baseline for this group was included in the model' (Table 2 rows 'BSL FPG 9.41' and 'BSL FPG(a) 9.02', footnote a 'BSL FPG value for 32 mg BID group'). (b) Level 3 gates the rosiglitazone treatment effect ROTE on KoutG, because no rosiglitazone plasma concentrations were collected and the effect is a fixed step rather than exposure-driven (Methods: 'In the rosiglitazone group, no plasma concentration data were collected during the treatment period and as such the treatment effect for rosiglitazone (ROTE) was included using a stimulatory step function on KoutG'). (c) Levels 0 vs 1/2 vs 3 select the HbA1c lifestyle effect LEHB (0.037 placebo, 0.037 - 0.017 = 0.020 sipoglitazar, 0 rosiglitazone per Results 'Lifestyle Effect Model': 'No significant lifestyle effect on HbA1c could be indentified on the rosiglitazone group'). The sipoglitazar DOSE LEVEL is NOT carried by TRT -- it is carried separately by DOSE_SIPOGLITAZAR_MGD, so a level-1 or level-2 subject must also have a non-zero daily dose set.",
      source_name        = "treatment group"
    ),
    DOSE_SIPOGLITAZAR_MGD = list(
      description        = "Subject's TOTAL daily sipoglitazar dose in mg/day, summed across both administrations for the twice-daily regimens. 0 mg/day for placebo and rosiglitazone subjects.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Stringer 2014 Methods: studied regimens are 8 mg QD, 16 mg QD, 16 mg BID, 32 mg QD, 32 mg BID and 64 mg QD, i.e. total daily doses of 8, 16, 32, 32, 64 and 64 mg/day. Enters the steady-state exposure computation AUC = DOSE_SIPOGLITAZAR_MGD / (CL * 24), which drives the Emax stimulation of KoutG (Equation 1). The total-daily-dose basis (rather than a per-administration dose) is fixed by the paper's own definition of the driver: Methods 'AUC50 is the AUC0-24h at steady state achieving half the maximal response' and 'individual exposure, AUC (AUC = dose/CL) over the dose interval at steady state'. It is confirmed numerically by the paper's own simulations -- Table 1 assigns 64 mg to UGT2B15*1/*1 (CL 5.04 L/h) giving AUC = 64 / 120.96 = 0.529 mg*day/L, DEF = 0.487 * 0.529 / (1.15 + 0.529) = 0.153, and a steady-state FPG of 9.41 / 1.153 = 8.16 mmol/L, matching the -1.2 mmol/L change from baseline reported in Results and the ~8.2 mmol/L asymptote of Figure 2a. The paper tested but did NOT retain a separate AUC50 between QD and BID regimens ('To explore any potential differences between daily dosing regimens, a different AUC50 value was tested between the BID and QD groups'; only one AUC50 appears in Table 2), so a 32 mg BID subject and a 64 mg QD subject share the same drug effect at the same clearance.",
      source_name        = "dose"
    ),
    UGT2B15_STAR2_HET = list(
      description        = "Binary germline genotype indicator for the UGT2B15*1/*2 heterozygote group: 1 = subject carries exactly one UGT2B15*2 allele; 0 = otherwise (the union of *1/*1 wild-type homozygotes and *2/*2 homozygous carriers, which the paired indicator UGT2B15_STAR2_HOM flags). Time-fixed per subject.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (UGT2B15*1/*1 wild-type homozygote, when UGT2B15_STAR2_HOM and UGT2B15_MISSING are also 0)",
      notes              = "Stringer 2014 Table 1 fixes the sipoglitazar clearance used to compute exposure for this stratum at 3.35 L/h. Genotype counts UGT2B15*1/*1 : *1/*2 : *2/*2 = 149 : 357 : 194 (Table S1 of the supplement), i.e. 21% : 51% : 28% of genotyped subjects (Methods 'Model Qualification'). UGT2B15*2 is the reduced-glucuronidation variant, so *1/*2 subjects have intermediate clearance and *2/*2 subjects the lowest -- Introduction: 'Higher plasma exposure of sipoglitazar was observed in the UGT2B15*2/*2 genotype than subjects homozygous for the wild-type allele UGT2B15*1/*1 (3.3-fold higher) or heterozygous allele UGT2B15*1/*2 (2.2-fold higher).'",
      source_name        = "UGT2B15*1/*2"
    ),
    UGT2B15_STAR2_HOM = list(
      description        = "Binary germline genotype indicator for the UGT2B15*2/*2 homozygous-variant group: 1 = subject carries two UGT2B15*2 alleles; 0 = otherwise (the union of *1/*1 wild-type homozygotes and *1/*2 heterozygotes, which the paired indicator UGT2B15_STAR2_HET flags). Time-fixed per subject.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (UGT2B15*1/*1 wild-type homozygote, when UGT2B15_STAR2_HET and UGT2B15_MISSING are also 0)",
      notes              = "Stringer 2014 Table 1 fixes the sipoglitazar clearance used to compute exposure for this stratum at 1.53 L/h -- 3.3-fold lower than the UGT2B15*1/*1 value of 5.04 L/h, matching the '3.3-fold higher' exposure statement in the Introduction. n = 194 of 700 genotyped subjects (supplement Table S1).",
      source_name        = "UGT2B15*2/*2"
    ),
    UGT2B15_MISSING = list(
      description        = "Binary indicator for a subject whose UGT2B15 genotype was not collected. 1 = genotype missing; 0 = genotype known (in which case UGT2B15_STAR2_HET and UGT2B15_STAR2_HOM together identify the diplotype, and both being 0 means *1/*1). When UGT2B15_MISSING is 1 both star-allele indicators must be 0.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (UGT2B15 genotype known)",
      notes              = "Stringer 2014 Methods 'Intra-Individual Variability and Residual Error': 'Genotype information was not collected in 10% of the population, however, these subjects were included in the analysis using an average clearance value for the population.' Supplement Table S1 footnote: 'genotype information not collected in 80 subjects' (80 of 780 = 10.3%). The paper does not print the numeric population-average clearance it used; the value carried in this model file (3.21 L/h) is DERIVED from the paper's own numbers as the genotype-frequency-weighted mean of the Table 1 clearances using the Table S1 genotype counts -- (149 * 5.04 + 357 * 3.35 + 194 * 1.53) / 700 = 3.205 L/h. See the vignette Assumptions and deviations section.",
      source_name        = "genotype not collected"
    )
  )

  # Screened during forward inclusion on the IIV of BSLG, BSLH, LEFPG and LEHB
  # but NOT retained in the final model, so these never appear in model().
  # Results 'Covariate Analysis': "However neither sex, weight, duration of
  # disease, or age had a significant effect on any of these parameters and no
  # covariates were retained in the final model."
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age at baseline.",
      units       = "years",
      type        = "continuous",
      notes       = "Screened by forward inclusion / backward elimination (Methods 'Covariate Analysis') on the IIV of BSLG, BSLH, LEFPG and LEHB; not retained. Cohort median 56 years (range 34-75), supplement Table S1."
    ),
    SEXF = list(
      description = "Sex indicator: 1 = female, 0 = male.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened and not retained (Results 'Covariate Analysis'). Cohort 388 male : 392 female, supplement Table S1."
    ),
    WT = list(
      description = "Body weight at baseline.",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened and not retained (Results 'Covariate Analysis'). Cohort median 88.8 kg (range 55-160), supplement Table S1."
    ),
    T_DIAG_DIAB = list(
      description = "Time since diagnosis of type 2 diabetes mellitus.",
      units       = "years",
      type        = "continuous",
      notes       = "Screened and not retained (Results 'Covariate Analysis'). Cohort median 1.0 year (range 0-30.9), supplement Table S1. All subjects were drug-naive at entry."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 780L,
    n_studies      = 2L,
    n_observations = "FPG and HbA1c sampled at -1, 0, 2, 4, 6, 8, 10 and 12 weeks in every subject (Methods 'Subjects and Data Collection'); the paper does not report an observation count.",
    age_range      = "34-75 years (inclusion criterion age > 35 and < 75 years; supplement Table S1 realised range 34-75)",
    age_median     = "56 years (supplement Table S1)",
    weight_range   = "55-160 kg (supplement Table S1)",
    weight_median  = "88.8 kg (supplement Table S1)",
    sex_female_pct = 50.3,
    race_ethnicity = "Not reported in the source paper. The Discussion contrasts the observed UGT2B15*1/*1 frequency (21%) against a literature Caucasian frequency of 19-22% and Asian American / Japanese American frequencies of 47% / 100%, implying a predominantly Caucasian cohort, but no race table is given.",
    disease_state  = "Drug-naive adults with type 2 diabetes mellitus. Inclusion criteria: diagnosis of type 2 diabetes with no prior exposure to anti-diabetic medication, screening HbA1c > 7.0% and < 10.0%, age > 35 and < 75 years. Baseline FPG median 9.3 mmol/L (range 2.9-20.8); baseline HbA1c median 7.9% (range 6.9-9.9); duration of disease median 1.0 year (range 0-30.9). Supplement Table S1.",
    dose_range     = "Sipoglitazar 8 mg QD, 16 mg QD, 16 mg BID, 32 mg QD, 32 mg BID or 64 mg QD (total daily dose 8-64 mg/day), rosiglitazone 8 mg QD, or placebo, each for 13 weeks. Disposition: sipoglitazar n = 572, rosiglitazone n = 72, placebo n = 136.",
    regions        = "Not reported in the source paper; the central laboratory was Medical Research Laboratories International, Brussels, Belgium.",
    genotype       = "UGT2B15*1/*1 : *1/*2 : *2/*2 = 149 : 357 : 194 among the 700 genotyped subjects (21% : 51% : 28%); genotype was not collected in the remaining 80 subjects (supplement Table S1 footnote).",
    notes          = "Two 13-week Phase II randomized, double-blind trials pooled for the analysis. All subjects received dietary advice for the entire trial duration, which is what the LEFPG / LEHB 'lifestyle effect' parameters absorb. Estimation in NONMEM 7.1 with FOCE-I and ADVAN6; simulations in Berkeley Madonna 8.3.13. Model stability assessed with 500 bootstrap replicates (92.6% minimised successfully); no eta shrinkage above 12%."
  )

  ini({
    # ------------------------------------------------------------------
    # All parameter estimates are from Stringer 2014 Table 2, "Summary of
    # Parameter Estimates for the Final Model Including Bootstrap
    # Estimates" (Model Estimate column), except the per-genotype
    # clearances which are from Table 1. The parenthesised numbers in
    # Table 2's Model Estimate column are CV% (relative standard errors),
    # not variability.
    #
    # TIME UNIT IS DAYS throughout: KoutG and KoutH are reported in
    # days^-1 and AUC50 in mg*day/L.
    # ------------------------------------------------------------------

    # ---- FPG turnover (Equation 3) ----
    lrbase_fpg       <- log(9.41);  label("Typical baseline FPG for all treatment groups except sipoglitazar 32 mg BID (mmol/L)")  # Table 2 'BSL FPG (mmol/L)' = 9.41 (CV 0.83%)
    lrbase_fpg_bid32 <- log(9.02);  label("Typical baseline FPG for the sipoglitazar 32 mg BID group (mmol/L)")                    # Table 2 'BSL FPG (mmol/L)a' = 9.02 (CV 1.45%), footnote a 'BSL FPG value for 32 mg BID group'
    lkout_fpg        <- log(0.027); label("First-order rate constant for removal of FPG, KoutG (1/day)")                           # Table 2 'KoutG (days-1)' = 0.027 (CV 8.1%)

    # ---- HbA1c turnover (Equation 4) ----
    lrbase_hba1c <- log(7.96);  label("Typical baseline HbA1c (%, NGSP)")                                       # Table 2 'BSL HbA1c (%)' = 7.96 (CV 0.46%)
    lkout_hba1c  <- log(0.031); label("First-order HbA1c degradation rate constant, KoutH (1/day)")             # Table 2 'KoutH (days-1)' = 0.031 (CV 6.1%)
    lgamma       <- log(0.71);  label("Power exponent on FPG driving HbA1c production, FPG^gamma (unitless)")   # Table 2 'Gamma' = 0.71 (CV 3.9%); printed as FPG^gamma in the Methods text and as FPG^lambda in Equation (4) -- the same parameter

    # ---- Sipoglitazar drug effect on KoutG (Equation 1) ----
    # DEF = Emax * AUC / (AUC50 + AUC).  Table 2 reports Emax on a PERCENT
    # scale ('Emax (%) 48.7'; Results: 'The Emax of this effect was estimated
    # at 49%'), while ROTE below is already a fraction ('ROTE 0.28';
    # Results: 'a population mean value of 28%').  Emax is therefore carried
    # here as the fraction 0.487.  The fraction scale is confirmed by the
    # paper's own simulations: it is the only scale under which the Figure 2a
    # steady-state FPG asymptotes (8.16 / 7.85 / 7.28 mmol/L for
    # UGT2B15*1/*1, *1/*2, *2/*2 at 64 mg) and the reported 6-month FPG
    # changes (-1.2 / -1.6 / -2.1 mmol/L) are reproduced.
    lemax  <- log(0.487); label("Maximal fractional stimulation of KoutG by sipoglitazar, Emax (fraction; Table 2 reports 48.7%)")  # Table 2 'Emax (%)' = 48.7 (CV 16.2%), rescaled from percent to fraction
    lauc50 <- log(1.15);  label("Steady-state AUC0-24h giving half the maximal sipoglitazar effect, AUC50 (mg*day/L)")             # Table 2 'AUC50 (mg day/L)' = 1.15 (CV 28%)

    # ---- Rosiglitazone comparator effect on KoutG ----
    e_rosi_kout_fpg <- 0.28; label("Rosiglitazone 8 mg QD treatment effect ROTE: fractional stimulation of KoutG (unitless)")  # Table 2 'ROTE' = 0.28 (CV 8.76%); a step function, not exposure-driven, because no rosiglitazone concentrations were collected

    # ---- 'Lifestyle' effects ----
    # LEFPG: a mean-zero ADDITIVE random effect on KinG.  The structural
    # (population) value is fixed at exactly zero -- Table 2 row 'LEFPGb'
    # reports 0 with footnote b 'LEFPG structural parameter was fixed at 0',
    # and the text below Equation (3) states 'LEFPG is included as an
    # additive random effect with the structural parameter fixed at zero'.
    e_lifestyle_kin_fpg <- fixed(0); label("Population mean lifestyle effect on KinG, LEFPG (unitless); all of the signal is in the between-subject random effect")  # Table 2 'LEFPGb' = 0, footnote b

    # LEHB: a step inhibition of KinH.  Table 2's second lifestyle row is
    # labelled 'LEHBactive' but its footnote c prints Equation (2),
    # LEHBactive = LEHBplacebo - LEHBfactor.  The tabulated 0.017 is
    # therefore the SUBTRACTED FACTOR, not the active value itself: the row
    # carries a CV% (32.2), every other CV% in Table 2 is the relative
    # standard error of an estimated NONMEM parameter taken from the
    # covariance matrix, and footnote c says LEHBactive is COMPUTED from two
    # other quantities -- so LEHBactive is not a THETA and would carry no
    # RSE, while LEHBfactor is.  Reading it as the factor gives LEHBactive = 0.037 -
    # 0.017 = 0.020 exactly, which is what Results 'Lifestyle Effect Model'
    # prints: 'a population mean decrease of 3.7% was identified for the
    # lifestyle effect in the placebo group, while in the actively treated
    # groups the reduction was slightly lower, 2%.'  (Had the row been the
    # active value, that sentence would read '1.7%' to match the one-decimal
    # precision of '3.7%'.)  The two readings differ by only 0.02% HbA1c in
    # every published output, i.e. below the paper's own reporting
    # precision; see the vignette Assumptions and deviations section.
    e_placebo_kin_hba1c <- 0.037; label("Lifestyle effect LEHB in the placebo arm: fractional inhibition of KinH (unitless)")                                                     # Table 2 'LEHBplacebo' = 0.037 (CV 12.5%)
    e_active_kin_hba1c  <- 0.017; label("LEHBfactor: amount subtracted from the placebo lifestyle effect in actively treated arms, so LEHBactive = 0.037 - 0.017 = 0.020")        # Table 2 'LEHBactivec' = 0.017 (CV 32.2%) read through footnote c / Equation (2)

    # ---- Sipoglitazar clearance by UGT2B15 genotype (Table 1) ----
    # These are NOT estimated in this PD analysis: they are the population
    # clearance values carried over from the companion population PK
    # analysis (reference 2) and reproduced in Table 1 of this paper, used
    # only to turn a daily dose into a steady-state AUC.  Reported in L/h;
    # converted to L/day inside model().
    lcl_ugt2b15_s1s1    <- fixed(log(5.04)); label("Sipoglitazar clearance in UGT2B15*1/*1 subjects (L/h)")   # Table 1 column header 'UGT2B15*1/*1 (CL = 5.04 L/h)'
    lcl_ugt2b15_s1s2    <- fixed(log(3.35)); label("Sipoglitazar clearance in UGT2B15*1/*2 subjects (L/h)")   # Table 1 column header 'UGT2B15*1/*2 (CL = 3.35 L/h)'
    lcl_ugt2b15_s2s2    <- fixed(log(1.53)); label("Sipoglitazar clearance in UGT2B15*2/*2 subjects (L/h)")   # Table 1 column header 'UGT2B15*2/*2 (CL = 1.53 L/h)'
    lcl_ugt2b15_missing <- fixed(log(3.21)); label("Sipoglitazar clearance for subjects without a collected genotype (L/h) -- DERIVED, see comment")  # DERIVED, not printed in the paper: Methods states the 10% ungenotyped subjects 'were included in the analysis using an average clearance value for the population'; the genotype-frequency-weighted mean of the Table 1 clearances using the supplement Table S1 counts is (149*5.04 + 357*3.35 + 194*1.53)/700 = 3.205 L/h

    # ------------------------------------------------------------------
    # Inter-individual variability.  Table 2 'Random effects: inter-
    # individual variability (IIV)' reports these on the VARIANCE scale
    # (the row labels are omega^2).  The variance reading is confirmed by
    # the bootstrap column: omega^2 BSL FPG = 0.05 with a bootstrap 95% CI
    # of 0.043-0.055, a relative half-width of ~12%, which matches the
    # ~10% expected for a variance estimated from 780 subjects
    # (sqrt(2/780) = 5.1% RSE) and is roughly twice what an SD-scale
    # estimate would show.
    #
    # NOTE ON TABLE 2's BOOTSTRAP COLUMN: the bootstrap 95% CIs printed
    # against 'omega^2 BSL FPG' (0.008-0.01) and 'omega^2 BSL HbA1c'
    # (0.043-0.055) are TRANSPOSED in the published table -- neither CI
    # contains its own point estimate, and swapping them makes both
    # contain theirs.  Only the bootstrap column is affected; the Model
    # Estimate values carried here are self-consistent.  Recorded in the
    # vignette Errata.
    #
    # The paper states 'The correlation between IIV on baselines was
    # included using the OMEGA BLOCK option' but does NOT report the
    # off-diagonal covariance between the two baseline etas, so it is
    # carried as zero here (see vignette Assumptions and deviations).
    # ------------------------------------------------------------------
    etalrbase_fpg   ~ 0.05  # Table 2 'omega^2 BSL FPG'   = 0.05 (CV 4.2%);  log-normal IIV on the FPG baseline
    etalrbase_hba1c ~ 0.01  # Table 2 'omega^2 BSL HbA1c' = 0.01 (CV 6.1%);  enters through the Box-Cox transformation below, not directly as exp(eta)

    # Both lifestyle effects carry ADDITIVE (not log-normal) between-subject
    # random effects.  This is forced for LEFPG, whose population value is
    # exactly zero -- a log-normal eta cannot be placed on a zero-valued
    # parameter -- and the paper says so explicitly ('implemented in the
    # model using an additive random effect').  LEHB is carried the same way
    # for consistency with its structurally parallel partner; the paper does
    # not state the LEHB eta's distribution (vignette Errata).
    # The LEHB eta is named for `e_placebo_kin_hba1c` because that is the base
    # fixed effect it attaches to: LEHB = LEHBplacebo - LEHBfactor * (active) +
    # eta, so the placebo theta is the parameter carrying the random effect and
    # the active-arm factor is a fixed offset from it.
    etae_lifestyle_kin_fpg ~ 0.04  # Table 2 'omega^2 LEFPG' = 0.04 (CV 9.3%);  additive random effect, SD 0.2
    etae_placebo_kin_hba1c ~ 0.01  # Table 2 'omega^2 LEHB'  = 0.01 (CV 10.6%); additive random effect, SD 0.1

    # Box-Cox shape parameter for the HbA1c baseline eta.  Methods:
    # 'IIV on the baseline for HbA1c was evaluated using a Box-Cox
    # transformation model, which was applied to account for the skewness
    # observed in the individual data' (reference 18 = Petersson KJ, Hanze E,
    # Savic RM, Karlsson MO, Pharm Res 2009;26:2174-2185).  Results: 'The
    # inclusion of the Box-Cox transformation on IIV for the baseline on
    # HbA1c (BSLH) resulted in a significant (P <= 0.001) change in the
    # objective function.'  A positive lambda skews the baseline
    # distribution to the RIGHT, consistent with the supplement Table S1
    # baseline HbA1c median 7.9% sitting closer to the lower bound (6.9)
    # than the upper (9.9).
    boxcox_rbase_hba1c <- fixed(2.39); label("Box-Cox shape parameter lambda for the IIV transformation on the HbA1c baseline (dimensionless)")  # Table 2 'Box-Cox' = 2.39 (CV 17.8%); carried as fixed() because it is a shape constant of the eta distribution rather than a structural PD parameter

    # ---- Residual error (Table 2 'Random effects: residual error') ----
    # Methods 'Intra-Individual Variability and Residual Error':
    # 'Residual variability was included using a proportional model.'
    propSd_glucose <- 0.106; label("Proportional residual error SD for FPG (fraction; Table 2 reports 10.6%)")    # Table 2 'Residual error FPG (%)'   = 10.6 (CV 4.8%)
    propSd_hba1c   <- 0.020; label("Proportional residual error SD for HbA1c (fraction; Table 2 reports 2.0%)")   # Table 2 'Residual error HbA1c (%)' = 2.0  (CV 6.8%)
  })

  model({
    # ------------------------------------------------------------------
    # Constants
    # ------------------------------------------------------------------
    hours_per_day <- 24   # Table 1 reports CL in L/h; the PD model runs in days

    # ------------------------------------------------------------------
    # 1. Treatment-cohort indicators derived from TRT.  Written as sums of
    #    equality tests (rather than logical || ) so each indicator is an
    #    unambiguous 0/1 numeric; the four TRT levels are mutually
    #    exclusive, so the sum below can only be 0 or 1.
    # ------------------------------------------------------------------
    is_placebo  <- (TRT == 0)
    is_sipo32bid <- (TRT == 2)
    is_sipo     <- (TRT == 1) + (TRT == 2)
    is_rosi     <- (TRT == 3)

    # ------------------------------------------------------------------
    # 2. Baselines.  The 32 mg BID sipoglitazar arm gets its own typical
    #    FPG baseline (Table 2 footnote a); every other arm shares 9.41.
    #    Log-normal IIV on the FPG baseline; Box-Cox-transformed eta on the
    #    HbA1c baseline (Petersson 2009 form, matching the existing
    #    GonzalezSales_2015_testosterone implementation).
    # ------------------------------------------------------------------
    rbase_fpg_typ <- exp(lrbase_fpg) * (1 - is_sipo32bid) +
                     exp(lrbase_fpg_bid32) * is_sipo32bid
    rbase_fpg     <- rbase_fpg_typ * exp(etalrbase_fpg)

    # Box-Cox transform of the HbA1c baseline eta: eta_bc = (exp(lambda *
    # eta) - 1) / lambda, then applied on the log scale.  eta = 0 maps to
    # eta_bc = 0, so the population median baseline is preserved at 7.96%.
    eta_rbase_hba1c_bc <- (exp(boxcox_rbase_hba1c * etalrbase_hba1c) - 1) /
                          boxcox_rbase_hba1c
    rbase_hba1c        <- exp(lrbase_hba1c) * exp(eta_rbase_hba1c_bc)

    # ------------------------------------------------------------------
    # 3. Sipoglitazar exposure.  CL is a fixed per-genotype constant
    #    (Table 1); AUC is the steady-state AUC0-24h = total daily dose /
    #    CL.  The wild-type UGT2B15*1/*1 stratum is the implicit reference
    #    (all three indicators zero).
    # ------------------------------------------------------------------
    cl <- exp(lcl_ugt2b15_s1s1) *
            (1 - UGT2B15_STAR2_HET - UGT2B15_STAR2_HOM - UGT2B15_MISSING) +
          exp(lcl_ugt2b15_s1s2)    * UGT2B15_STAR2_HET +
          exp(lcl_ugt2b15_s2s2)    * UGT2B15_STAR2_HOM +
          exp(lcl_ugt2b15_missing) * UGT2B15_MISSING

    auc <- DOSE_SIPOGLITAZAR_MGD / (cl * hours_per_day)   # mg*day/L

    # ------------------------------------------------------------------
    # 4. Drug and treatment effects on KoutG.
    #    Equation (1): DEF = Emax * AUC / (AUC50 + AUC).
    #    Placebo and rosiglitazone subjects carry a zero sipoglitazar daily
    #    dose, so auc = 0 and def = 0 for them without further gating.
    # ------------------------------------------------------------------
    emax  <- exp(lemax)
    auc50 <- exp(lauc50)
    def   <- emax * auc / (auc50 + auc)
    rote  <- e_rosi_kout_fpg * is_rosi

    # ------------------------------------------------------------------
    # 5. Lifestyle effects.
    #    LEFPG: population value fixed at 0 plus an additive random effect.
    #    LEHB : 0.037 on placebo, 0.037 - 0.017 = 0.020 on the sipoglitazar
    #           arms, and structurally absent on rosiglitazone (Results:
    #           'No significant lifestyle effect on HbA1c could be
    #           indentified on the rosiglitazone group'), so both the
    #           typical value and the random effect are switched off there.
    # ------------------------------------------------------------------
    lefpg <- e_lifestyle_kin_fpg + etae_lifestyle_kin_fpg

    lehb_typ <- e_placebo_kin_hba1c * is_placebo +
                (e_placebo_kin_hba1c - e_active_kin_hba1c) * is_sipo
    lehb     <- (lehb_typ + etae_placebo_kin_hba1c) * (1 - is_rosi)

    # ------------------------------------------------------------------
    # 6. Turnover rate constants.  KinG and KinH are not estimated: they
    #    are pinned by the requirement that each pool sit exactly at its
    #    baseline at t = 0 with no drug and no lifestyle effect --
    #      0 = KinG - KoutG * BSL_FPG                      => KinG = KoutG * BSL_FPG
    #      0 = KinH * BSL_FPG^gamma - KoutH * BSL_HbA1c    => KinH = KoutH * BSL_HbA1c / BSL_FPG^gamma
    #    Both are formed per subject from that subject's own baselines, so
    #    the baseline IIV propagates into the production rates.
    # ------------------------------------------------------------------
    gamma      <- exp(lgamma)
    kout_fpg   <- exp(lkout_fpg)
    kout_hba1c <- exp(lkout_hba1c)
    kin_fpg    <- kout_fpg * rbase_fpg
    kin_hba1c  <- kout_hba1c * rbase_hba1c / rbase_fpg^gamma

    # ------------------------------------------------------------------
    # 7. ODE system -- Stringer 2014 Equations (3) and (4).
    #      dFPG/dt   = KinG * (1 + LEFPG) - KoutG * (1 + DEF + ROTE) * FPG
    #      dHbA1c/dt = KinH * (1 - LEHB) * FPG^gamma - KoutH * HbA1c
    #    Both states are carried on the measured scale (mmol/L and %),
    #    not as amounts.
    #
    #    Equation (3) is printed in the paper as 'KinG * (1 +/- LEFPG)' -- a
    #    literal typeset plus-minus sign on p. 455, not an extraction
    #    artefact (it is the only such glyph in the article).  The
    #    plus/minus is immaterial: LEFPG is a mean-zero, symmetrically
    #    distributed additive random effect, so (1 + LEFPG) and (1 - LEFPG)
    #    generate identical distributions.  The '+' form is used here because
    #    it is the one that makes the sign of an individual LEFPG readable
    #    directly (positive LEFPG = raised glucose production = loss of
    #    glycemic control).
    # ------------------------------------------------------------------
    d/dt(glucose) <- kin_fpg * (1 + lefpg) -
                     kout_fpg * (1 + def + rote) * glucose
    d/dt(hba1c)   <- kin_hba1c * (1 - lehb) * glucose^gamma -
                     kout_hba1c * hba1c

    glucose(0) <- rbase_fpg
    hba1c(0)   <- rbase_hba1c

    # ------------------------------------------------------------------
    # 8. Observations and residual error (proportional on both endpoints).
    # ------------------------------------------------------------------
    glucose ~ prop(propSd_glucose)
    hba1c   ~ prop(propSd_hba1c)
  })
}
