Keunecke_2020_regorafenib_phase3 <- function() {
  description <- "Final joint parent + metabolite population PK model for oral regorafenib and its two pharmacologically active metabolites M-2 (N-oxide, BAY 75-7495) and M-5 (N-oxide N-desmethyl, BAY 81-8752), fitted to sparsely sampled data from four phase 3 studies in metastatic colorectal cancer, gastrointestinal stromal tumour and hepatocellular carcinoma (Keunecke 2020 Table 3). Structure is identical to the phase 1 model (Keunecke_2020_regorafenib_phase1): each analyte has a two-compartment disposition plus a gallbladder reservoir giving enterohepatic circulation, released during three 30-minute post-prandial windows per day; M-2 is formed pre-systemically through its own three-transit absorption chain and M-5 systemically from parent clearance. This version adds the retained covariate effects - sex and body mass index on parent clearance, and sex on both metabolite clearances - with women having lower clearance of all three analytes. Because the metabolite clearance random effects were almost perfectly correlated, the authors reduced the model to a single shared metabolite eta rescaled by a factor, so M-5 uses the M-2 eta multiplied by sd_ratio_cl_m5. All volumes and clearances are apparent, relative to oral bioavailability and to the parent fraction (1 - fm_m2)."
  reference <- paste(
    "Keunecke A, Hoefman S, Drenth H-J, Zisowsky J, Cleton A, Ploeger BA.",
    "Population pharmacokinetics of regorafenib in solid tumours:",
    "Exposure in clinical practice considering enterohepatic circulation",
    "and food intake.",
    "Br J Clin Pharmacol. 2020;86(12):2362-2376.",
    "doi:10.1111/bcp.14334.",
    sep = " "
  )
  vignette <- "Keunecke_2020_regorafenib"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    SEXF = list(
      description        = "Sex (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "male (SEXF = 0)",
      notes              = "Retained on parent, M-2 and M-5 clearance. Entered as a multiplicative factor (1 - theta * SEXF) with male as the reference, which is how Appendix Section 3.1 reports it: 'women had on average a 17% lower CLP than men (2.53 L h-1 and 3.05 L h-1, respectively)', and 3.05 * (1 - 0.169) = 2.53. Sex was the covariate giving the largest drop in objective function on CL_P (-20.7 points) and the only one carried through for the metabolites, because adding a second metabolite covariate produced numerical instability (Appendix Section 3.3). Direction confirmed against Figure 4: the male:female ratio of median parent AUC is ~0.86 (within the 0.8-1.25 bioequivalence bounds) and of protein-binding-corrected total AUC is ~0.71 (below the lower bound) - i.e. female exposure is HIGHER for all three analytes. Cohort: 27.9% female (Appendix Table A2).",
      source_name        = "Sex"
    ),
    BMI = list(
      description        = "Body mass index at baseline",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Retained on parent clearance only. Appendix Section 3.1 states the selected continuous covariates 'were included in the model ... with exponential function ... acting on the typical value of apparent clearance (CL)', with coefficient -0.363, and that this 'implicates an expected 0.15% higher CLP for a patient with a 0.1 kg m-2 reduced BMI'. That statement fixes the parameterisation as exp(theta * (BMI - BMIref) / BMIref) with BMIref = 24.5 kg/m^2, the phase 3 median in Appendix Table A2: 0.363 / 24.5 * 0.1 = 0.148%. An un-normalised exp(theta * (BMI - BMIref)) would give 3.7% instead, and would predict an 18-fold AUC ratio between the BMI >= 30 and BMI < 30 groups where Figure 4A shows ~1.10. Higher BMI lowers clearance, which the authors attribute to reduced CYP3A4 expression at high BMI. Cohort median 24.5, 5-95th percentile 18.5-32.4 kg/m^2 (Appendix Table A2).",
      source_name        = "BMI"
    )
  )

  # Screened in the univariate GAM but NOT retained in the final model
  # (Appendix Section 3.1 and 3.3, and section 3.3.2 / 3.3.4 of the paper).
  covariatesDataExcluded <- list(
    RACE_ASIAN = list(
      description = "Asian race (vs rest of world)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Statistically significant on CL_P in the univariate GAM, but excluded from the model-based analysis because the category contained <15% of the total (Appendix Section 3.1)."
    ),
    CHILDPUGH = list(
      description = "Child-Pugh liver score",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Significant in the univariate GAM but <15% of the total in one category, so not carried into the model-based analysis (Appendix Section 3.1). 63.1% missing overall (Appendix Table A2)."
    ),
    CONMED_UGT1A9_INDUCER = list(
      description = "UGT1A9 inducer received during the treatment period",
      units       = "(binary)",
      type        = "binary",
      notes       = "Significant in the univariate GAM but <15% of the total in one category (0.4% yes), so not carried into the model-based analysis (Appendix Section 3.1)."
    ),
    TUMTP = list(
      description = "Tumour type (colorectal carcinoma / gastrointestinal stromal tumour / hepatocellular carcinoma)",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Slope not significantly different from zero in the univariate GAM, so not taken forward (Appendix Section 3.1). Cohort 54% CRC, 8.8% GIST, 37.2% HCC."
    ),
    HGB = list(
      description = "Baseline haemoglobin",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Slope not significantly different from zero in the univariate GAM for CL_P (Appendix Section 3.1). Significant for CL_M-2 but not pursued (Appendix Section 3.3)."
    ),
    TPROT = list(
      description = "Baseline total protein",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Slope not significantly different from zero in the univariate GAM for CL_P (Appendix Section 3.1)."
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Slope not significantly different from zero in the univariate GAM for CL_P (Appendix Section 3.1). It was the second most significant covariate for CL_M-5, but adding it after sex produced numerical issues and poorly estimated parameters (Appendix Section 3.3). Cohort median 70 kg."
    ),
    BILI = list(
      description = "Baseline total bilirubin",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Slope not significantly different from zero in the univariate GAM for CL_P (Appendix Section 3.1); significant for CL_M-5 but not pursued."
    ),
    ALT = list(
      description = "Baseline alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Carried into the model-based covariate analysis for CL_P but gave a non-significant drop in objective function once sex and BMI were in the model (Appendix Section 3.1)."
    ),
    EGFR = list(
      description = "Baseline estimated glomerular filtration rate",
      units       = "mL/min/1.73m^2",
      type        = "continuous",
      notes       = "Carried into the model-based covariate analysis for CL_P but gave a non-significant drop in objective function (Appendix Section 3.1)."
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Significant in the univariate GAM for CL_M-2, but adding it after sex led to numerical issues and poorly estimated parameters, so metabolite covariate evaluation stopped at sex (Appendix Section 3.3)."
    ),
    AST = list(
      description = "Baseline aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Significant in the univariate GAM for CL_M-2 but not pursued (Appendix Section 3.3)."
    )
  )

  compartmentData <- list(
    depot            = list(analyte = "regorafenib", units = "mg", specimen = "administration site", verified = TRUE),
    transit1         = list(analyte = "regorafenib", units = "mg", specimen = "administration site", verified = TRUE),
    transit2         = list(analyte = "regorafenib", units = "mg", specimen = "administration site", verified = TRUE),
    central          = list(analyte = "regorafenib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1      = list(analyte = "regorafenib", units = "mg", specimen = "plasma", verified = TRUE),
    gallbladder      = list(analyte = "regorafenib", units = "mg", specimen = "bile", verified = TRUE),
    transit1_m2      = list(analyte = "regorafenib M-2 (N-oxide, BAY 75-7495)", units = "mg", specimen = "administration site", verified = TRUE),
    transit2_m2      = list(analyte = "regorafenib M-2 (N-oxide, BAY 75-7495)", units = "mg", specimen = "administration site", verified = TRUE),
    transit3_m2      = list(analyte = "regorafenib M-2 (N-oxide, BAY 75-7495)", units = "mg", specimen = "administration site", verified = TRUE),
    central_m2       = list(analyte = "regorafenib M-2 (N-oxide, BAY 75-7495)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_m2   = list(analyte = "regorafenib M-2 (N-oxide, BAY 75-7495)", units = "mg", specimen = "plasma", verified = TRUE),
    gallbladder_m2   = list(analyte = "regorafenib M-2 (N-oxide, BAY 75-7495)", units = "mg", specimen = "bile", verified = TRUE),
    central_m5       = list(analyte = "regorafenib M-5 (N-oxide N-desmethyl, BAY 81-8752)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_m5   = list(analyte = "regorafenib M-5 (N-oxide N-desmethyl, BAY 81-8752)", units = "mg", specimen = "plasma", verified = TRUE),
    gallbladder_m5   = list(analyte = "regorafenib M-5 (N-oxide N-desmethyl, BAY 81-8752)", units = "mg", specimen = "bile", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 906L,
    n_studies      = 4L,
    age_range      = "40-78 years (5th-95th percentile)",
    age_median     = "61 years",
    weight_range   = "48-99 kg (5th-95th percentile)",
    weight_median  = "70 kg",
    bmi_range      = "18.5-32.4 kg/m^2 (5th-95th percentile)",
    bmi_median     = "24.5 kg/m^2",
    sex_female_pct = 27.9,
    race_ethnicity = c(Caucasian = 54.4, Asian = 33.8, Black = 1.1, AmericanIndian = 0.1, Mixed = 0.2, Missing = 10.4),
    disease_state  = "Adults with treatment-refractory metastatic colorectal carcinoma (54%), advanced gastrointestinal stromal tumour (8.8%) or hepatocellular carcinoma after sorafenib (37.2%).",
    dose_range     = "Regorafenib 160 mg once daily (4 x 40 mg tablets), 3 weeks on / 1 week off, versus matching placebo.",
    regions        = "North America, Europe, Israel, Australia, South America and Asia.",
    hepatic_function = "Category A (normal) 56.4%, B1 28.7%, B2 11.7%, C 3.2% (Appendix Table A2 footnote, defined on total bilirubin and AST/ALT).",
    notes          = "Keunecke 2020 Appendix Tables A1 and A2. The four phase 3 studies are CORRECT (14387, NCT01103323, n = 388, 4059 PK observations), GRID (14874, NCT01271712, n = 81, 345 observations), CONCUR (15808, NCT01584830, n = 98, 1065 observations) and RESORCE (15982, NCT01774344, n = 339, 3210 observations). Sampling was sparse: predose trough at cycle 1 day 15 and cycle 2 day 15 in most patients, with 2-4 h and 5-10 h post-dose samples in subsets. Combined with the 62 phase 1 patients the full dataset is 968 patients and 10019 observations. Mealtimes were not recorded in the phase 3 studies; a 1 h interval between breakfast and dosing gave the lowest objective function and was imputed for all phase 3 patients (section 3.3.1), with lunch 4 h and dinner 10 h after breakfast per section 2.2. Median eGFR 100 mL/min/1.73m^2, median albumin 4.0 g/dL, median baseline ALT 25 U/L. Baseline distributions per study are in Appendix Table A2."
  )

  ini({
    # ---------------------------------------------------------------------
    # Absorption. Unchanged from the phase 1 model and fixed there.
    # ---------------------------------------------------------------------
    lka <- fixed(log(0.482)); label("First-order absorption / transit rate constant (1/h)")  # Keunecke 2020 Table 3: ka = 0.482 1/h, Fixed

    # ---------------------------------------------------------------------
    # Formation fractions, on the logit scale, both fixed to the phase 1
    # parent-metabolite estimates (Appendix Section 3.2). Table 3 footnote:
    # "FRM2 and FRM5 on logit scale corresponds with 41% and 25%".
    # ---------------------------------------------------------------------
    logitfm_m2 <- fixed(-0.355); label("Logit of the pre-systemic fraction of dose forming M-2 (unitless)")  # Keunecke 2020 Table 3: FRM2 = -0.355, Fixed (estimated in Table 2)
    logitfm_m5 <- fixed(-1.09);  label("Logit of the fraction of parent clearance forming M-5 (unitless)")   # Keunecke 2020 Table 3: FRM5 = -1.09, Fixed (estimated in Table 2)

    # ---------------------------------------------------------------------
    # Parent disposition. All apparent, all fixed to the phase 3 covariate
    # regorafenib model. CL_P dropped from 4.02 L/h in phase 1 to 3.05 L/h
    # here; section 3.3.1 notes the phase 1 model over-predicted phase 3
    # clearance, so CL_P and its IIV were re-estimated on the phase 3 data
    # (the base-model typical value before covariates was 2.9 L/h).
    # ---------------------------------------------------------------------
    lcl <- fixed(log(3.05)); label("Apparent parent clearance CL_P/(1-FRM2)/F for a male at reference BMI (L/h)")  # Keunecke 2020 Table 3: CL_P/(1-FRM2) = 3.05 L/h, Fixed
    lvc <- fixed(log(10.7)); label("Apparent parent central volume VC_P/(1-FRM2)/F (L)")                           # Keunecke 2020 Table 3: VC_P/(1-FRM2) = 10.7 L, Fixed
    lq  <- fixed(log(11.0)); label("Apparent parent inter-compartmental clearance Q_P/(1-FRM2)/F (L/h)")           # Keunecke 2020 Table 3: Q_P/(1-FRM2) = 11.0 L/h, Fixed
    lvp <- fixed(log(162));  label("Apparent parent peripheral volume VP_P/(1-FRM2)/F (L)")                        # Keunecke 2020 Table 3: VP_P/(1-FRM2) = 162 L, Fixed

    # ---------------------------------------------------------------------
    # Enterohepatic circulation. Identical to phase 1.
    # ---------------------------------------------------------------------
    lkbm  <- fixed(log(0.141)); label("Central-to-gallbladder transfer rate constant k_CG (1/h)")  # Keunecke 2020 Table 3: k_CG = 0.141 1/h, Fixed
    lkehc <- fixed(log(100));   label("Gallbladder emptying rate constant k_GE (1/h)")             # Keunecke 2020 Table 3: k_GE = 100 1/h, Fixed

    tmeal1 <- fixed(8);   label("Clock time of breakfast, the first gallbladder-emptying trigger (h since midnight)")  # Keunecke 2020 Table 3: Breakfast = 08:00, Fixed
    tmeal2 <- fixed(12);  label("Clock time of lunch, the second gallbladder-emptying trigger (h since midnight)")     # Keunecke 2020 Table 3: Lunch = 12:00, Fixed
    tmeal3 <- fixed(18);  label("Clock time of dinner, the third gallbladder-emptying trigger (h since midnight)")     # Keunecke 2020 Table 3: Dinner = 18:00, Fixed
    dge    <- fixed(0.5); label("Duration of gallbladder emptying after each meal (h)")                                # Keunecke 2020 Table 3: DGE = 0.5 h, Fixed

    # ---------------------------------------------------------------------
    # Metabolite clearances, re-estimated on the phase 3 data.
    # ---------------------------------------------------------------------
    lcl_m2 <- log(1.99); label("Apparent M-2 clearance CL_M-2/F for a male (L/h)")  # Keunecke 2020 Table 3: CL_M-2 = 1.99 L/h, SE 0.00508, RSE 0.255%, 95% CI 1.98 to 2.00
    lcl_m5 <- log(1.42); label("Apparent M-5 clearance CL_M-5/F for a male (L/h)")  # Keunecke 2020 Table 3: CL_M-5 = 1.42 L/h, SE 0.0532, RSE 3.76%, 95% CI 1.31 to 1.52

    # ---------------------------------------------------------------------
    # Covariate effects. Sex enters as a multiplicative factor with male as
    # the reference category, BMI as a median-normalised exponential.
    # See the covariateData notes for the two independent checks that fix
    # the BMI parameterisation and the sign of the sex effect.
    # ---------------------------------------------------------------------
    e_sexf_cl    <- fixed(0.169); label("Fractional reduction in parent clearance for female sex (unitless)")  # Keunecke 2020 Table 3: Sex on regorafenib clearance = 0.169, Fixed; Appendix Section 3.1 "women had on average a 17% lower CLP than men (2.53 L h-1 and 3.05 L h-1)"
    e_bmi_cl     <- fixed(-0.363); label("Exponential coefficient for median-normalised BMI on parent clearance (unitless)")  # Keunecke 2020 Table 3: BMI on regorafenib clearance = -0.363, Fixed; Appendix Section 3.1 "The coefficient of the exponential term was estimated at -0.363 ... an expected 0.15% higher CLP for a patient with a 0.1 kg m-2 reduced BMI"
    e_sexf_cl_m2 <- 0.380; label("Fractional reduction in M-2 clearance for female sex (unitless)")  # Keunecke 2020 Table 3: Sex on M-2 clearance = 0.380, SE 0.0192, RSE 5.06%, 95% CI 0.342 to 0.418
    e_sexf_cl_m5 <- 0.761; label("Fractional reduction in M-5 clearance for female sex (unitless)")  # Keunecke 2020 Table 3: Sex on M-5 clearance = 0.761, SE 0.0387, RSE 5.08%, 95% CI 0.685 to 0.837

    # Ratio of the M-5 to the M-2 clearance random-effect SD. Appendix
    # Section 3.3: the variance of IIV on CL_M-5 could not be estimated
    # reliably (zero gradients), so the authors reduced the model "by using
    # the same ETA for CLM-2 as for CLM-5 but allowing for different
    # variances by estimating a factor between these variances". The M-5 eta
    # is therefore deterministic given the M-2 eta, which is what keeps the
    # phase 3 OMEGA well conditioned where the phase 1 3x3 block is not.
    sd_ratio_cl_m5 <- 2.21; label("Ratio of the M-5 to the M-2 clearance random-effect SD (unitless)")  # Keunecke 2020 Table 3: Factor IIV CL_M-5 = 2.21, SE 0.0752, RSE 3.41%, 95% CI 2.06 to 2.35; footnote "The estimated factor for IIV on CLM-5 corresponds with a variance of 0.385 * (2.21^2) = 1.88"

    # ---------------------------------------------------------------------
    # Inter-individual variability (exponential model, Appendix Section 2,
    # so the tabulated omega^2 are log-scale variances).
    # ---------------------------------------------------------------------
    # Keunecke 2020 Table 3, row "omega^2 k a" = 0.127, held constant. Eta
    # traces are kept on their own line above the declaration: rxode2 rewrites
    # a trailing comment on an eta line into a label(), which then duplicates
    # fixed().
    etalka ~ fixed(0.127)

    #
    # Source trace for the block below, kept on its own lines above the
    # declaration: a trailing comment INSIDE a multi-line omega c(...) is
    # rewritten into a bare `;` and makes the model unparseable.
    #   Keunecke 2020 Table 3, row "omega^2 CL_P"            = 0.189 (SE 0.0123, RSE 6.53%)
    #   Keunecke 2020 Table 3, row "omega CL_P,CL_M-2/M-5"   = 0.206 (RSE 8.16%)
    #   Keunecke 2020 Table 3, row "omega^2 CL_M-2"          = 0.385 (RSE 5.68%)
    etalcl + etalcl_m2 ~ c(
      0.189,
      0.206, 0.385
    )

    # Keunecke 2020 Table 3, row "omega^2 FRM2" = 0.156, held constant at the
    # Table 2 (phase 1 parent-metabolite) estimate.
    etalogitfm_m2 ~ fixed(0.156)
    # Keunecke 2020 Table 3, row "omega^2 FRM5" = 0.841, held constant at the
    # Table 2 (phase 1 parent-metabolite) estimate.
    etalogitfm_m5 ~ fixed(0.841)

    # ---------------------------------------------------------------------
    # Residual error. Larger than in phase 1 for all three analytes, as
    # expected for sparse sampling with imputed rather than recorded meals.
    # ---------------------------------------------------------------------
    propSd    <- 0.543;        label("Parent proportional residual error (fraction)")  # Keunecke 2020 Table 3: Parent Prop. error SD = 0.543, SE 0.00124, RSE 0.229%
    addSd_m2  <- fixed(0.001); label("M-2 additive residual error (mg/L)")             # Keunecke 2020 Table 3: M-2 Add. error SD = 0.001, Fixed
    propSd_m2 <- 0.485;        label("M-2 proportional residual error (fraction)")     # Keunecke 2020 Table 3: M-2 Prop. error SD = 0.485, SE 0.00938, RSE 1.93%
    addSd_m5  <- fixed(0.001); label("M-5 additive residual error (mg/L)")             # Keunecke 2020 Table 3: M-5 Add. error SD = 0.001, Fixed
    propSd_m5 <- 1.14;         label("M-5 proportional residual error (fraction)")     # Keunecke 2020 Table 3: M-5 Prop. error SD = 1.14, SE 0.113, RSE 9.93%
  })

  model({
    # 1. Covariate factors.
    #    Sex: multiplicative, male reference (Appendix Section 3.1).
    #    BMI: exponential in the median-normalised deviation from the phase 3
    #    median of 24.5 kg/m^2 (Appendix Table A2). The normalisation is what
    #    reproduces the Appendix's "0.15% higher CL_P per 0.1 kg/m^2 lower
    #    BMI": 0.363 / 24.5 * 0.1 = 0.148%.
    sex_cl    <- 1 - e_sexf_cl    * SEXF
    sex_cl_m2 <- 1 - e_sexf_cl_m2 * SEXF
    sex_cl_m5 <- 1 - e_sexf_cl_m5 * SEXF
    bmi_cl    <- exp(e_bmi_cl * (BMI - 24.5) / 24.5)

    # 2. Individual parameters.
    ka    <- exp(lka + etalka)
    cl    <- exp(lcl + etalcl) * sex_cl * bmi_cl
    vc    <- exp(lvc)
    q     <- exp(lq)
    vp    <- exp(lvp)
    kbm   <- exp(lkbm)
    kehc  <- exp(lkehc)

    cl_m2 <- exp(lcl_m2 + etalcl_m2) * sex_cl_m2
    # One shared metabolite eta, rescaled for M-5 (Appendix Section 3.3).
    cl_m5 <- exp(lcl_m5 + sd_ratio_cl_m5 * etalcl_m2) * sex_cl_m5

    # Individual logit-scale values collected on their own lines so the etas
    # stay mu-referenced, then inverse-logit back onto (0, 1).
    logitfm_m2_ind <- logitfm_m2 + etalogitfm_m2
    logitfm_m5_ind <- logitfm_m5 + etalogitfm_m5
    fm_m2 <- 1 / (1 + exp(-logitfm_m2_ind))
    fm_m5 <- 1 / (1 + exp(-logitfm_m5_ind))

    # 3. Food-intake switch. `t` is absolute time, so the event table origin
    #    must be midnight; tod is the resulting time of day. In the phase 3
    #    studies mealtimes were not recorded and breakfast was imputed 1 h
    #    before the dose (section 3.3.1), i.e. dosing at 09:00.
    tod   <- t - 24 * floor(t / 24)
    meal  <- (tod >= tmeal1) * (tod < tmeal1 + dge) +
             (tod >= tmeal2) * (tod < tmeal2 + dge) +
             (tod >= tmeal3) * (tod < tmeal3 + dge)

    # 4. Parent regorafenib (Figure 1 caption ODEs for A3, A4, A5).
    d/dt(depot)       <- -ka * depot
    d/dt(transit1)    <-  ka * (1 - fm_m2) * depot - ka * transit1
    d/dt(transit2)    <-  ka * transit1 - ka * transit2
    d/dt(central)     <-  ka * transit2 -
                            kbm * central -
                            (cl / vc) * central -
                            (q / vc) * central + (q / vp) * peripheral1 +
                            meal * kehc * gallbladder
    d/dt(peripheral1) <-  (q / vc) * central - (q / vp) * peripheral1
    d/dt(gallbladder) <-  kbm * central - meal * kehc * gallbladder

    # 5. M-2, formed pre-systemically through its own three-transit chain.
    d/dt(transit1_m2)    <-  ka * fm_m2 * depot - ka * transit1_m2
    d/dt(transit2_m2)    <-  ka * transit1_m2 - ka * transit2_m2
    d/dt(transit3_m2)    <-  ka * transit2_m2 - ka * transit3_m2
    d/dt(central_m2)     <-  ka * transit3_m2 -
                               kbm * central_m2 -
                               (cl_m2 / vc) * central_m2 -
                               (q / vc) * central_m2 + (q / vp) * peripheral1_m2 +
                               meal * kehc * gallbladder_m2
    d/dt(peripheral1_m2) <-  (q / vc) * central_m2 - (q / vp) * peripheral1_m2
    d/dt(gallbladder_m2) <-  kbm * central_m2 - meal * kehc * gallbladder_m2

    # 6. M-5, formed systemically as the fraction fm_m5 of parent clearance.
    d/dt(central_m5)     <-  fm_m5 * (cl / vc) * central -
                               kbm * central_m5 -
                               (cl_m5 / vc) * central_m5 -
                               (q / vc) * central_m5 + (q / vp) * peripheral1_m5 +
                               meal * kehc * gallbladder_m5
    d/dt(peripheral1_m5) <-  (q / vc) * central_m5 - (q / vp) * peripheral1_m5
    d/dt(gallbladder_m5) <-  kbm * central_m5 - meal * kehc * gallbladder_m5

    # 7. Bioavailability: the tabulated parent parameters are divided by
    #    (1 - FRM2), so the parent chain must receive the full nominal dose.
    f(depot) <- 1 / (1 - fm_m2)

    # 8. Observations. Metabolite volumes are the parent's by assumption.
    Cc    <- central / vc
    Cc_m2 <- central_m2 / vc
    Cc_m5 <- central_m5 / vc

    Cc    ~ prop(propSd)
    Cc_m2 ~ add(addSd_m2) + prop(propSd_m2)
    Cc_m5 ~ add(addSd_m5) + prop(propSd_m5)
  })
}
