Nguyen_2025_valbenazine <- function() {
  description <- paste(
    "Joint parent-metabolite population pharmacokinetic model of valbenazine",
    "and its active metabolite [+]-alpha-dihydrotetrabenazine",
    "([+]-alpha-HTBZ, NBI-98782) in healthy adults and in patients with",
    "tardive dyskinesia or Huntington's disease chorea (Nguyen 2025).",
    "Valbenazine is absorbed through a chain of four transit compartments",
    "(transit rate constant KTR = KA) into a two-compartment linear",
    "disposition system; a fraction FM = 0.207 of the valbenazine",
    "elimination flux forms [+]-alpha-HTBZ, which itself follows",
    "two-compartment linear disposition. The two analytes SHARE the apparent",
    "central volume VC/F -- the identifiability constraint that makes FM",
    "estimable (previously fixed at 0.35 in the predecessor model).",
    "Relative bioavailability rises saturably with dose (Emax form anchored",
    "at the 1 mg dose, FMAX = 1.36, ED50 = 90.1 mg); the effect is small at",
    "clinically relevant doses (1.09- and 1.16-fold at 60 and 80 mg relative",
    "to 40 mg). Covariates: allometric body weight on CLP/F and VC/F",
    "(70 kg reference), oral-solution formulation and fed state on KTR, fed",
    "state on F1, and CYP2D6 poor- and intermediate-metabolizer status on the",
    "metabolite clearance CLM (51.6% and 28.1% reductions relative to",
    "extensive or ultra-rapid metabolizers). Residual error is log-additive",
    "and stratified by analyte and by sampling design (intensive vs sparse),",
    "switched per observation record through SAMPLE_INTENSIVE.",
    sep = " "
  )
  reference <- paste(
    "Nguyen HQ, Crass RL, Chapel S, Kuan HYS, Loewen G, Brar S.",
    "Population pharmacokinetic and exposure-efficacy analyses of",
    "valbenazine in patients with Huntington's disease: supporting dose",
    "selection for chorea management.",
    "The Journal of Clinical Pharmacology. 2025;65(12):1777-1788.",
    "doi:10.1002/jcph.70092.",
    "Parameter estimates from Table 1; covariate-parameter equations from",
    "Supplemental Table S2; individual empirical Bayes parameter",
    "distributions from Supplemental Table S3.",
    "Companion exposure-efficacy model from the same paper:",
    "modellib('Nguyen_2025_valbenazine_tmc').",
    sep = " "
  )
  vignette <- "Nguyen_2025_valbenazine_huntington_chorea"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte and specimen confirmed against the Nguyen 2025
  # Results ("two-compartment, linear disposition for valbenazine and
  # [+]-alpha-HTBZ and transit compartment absorption for valbenazine") and
  # Figure S3 (final joint parent-metabolite model structure).
  compartmentData <- list(
    depot                = list(analyte = "valbenazine",                        units = "mg", specimen = "administration site", verified = TRUE),
    transit1             = list(analyte = "valbenazine",                        units = "mg", specimen = "administration site", verified = TRUE),
    transit2             = list(analyte = "valbenazine",                        units = "mg", specimen = "administration site", verified = TRUE),
    transit3             = list(analyte = "valbenazine",                        units = "mg", specimen = "administration site", verified = TRUE),
    transit4             = list(analyte = "valbenazine",                        units = "mg", specimen = "administration site", verified = TRUE),
    central              = list(analyte = "valbenazine",                        units = "mg", specimen = "plasma",               verified = TRUE),
    peripheral1          = list(analyte = "valbenazine",                        units = "mg", specimen = "plasma",               verified = TRUE),
    central_htbz         = list(analyte = "[+]-alpha-dihydrotetrabenazine",     units = "mg", specimen = "plasma",               verified = TRUE),
    peripheral1_htbz     = list(analyte = "[+]-alpha-dihydrotetrabenazine",     units = "mg", specimen = "plasma",               verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight. Power (allometric) effect on the valbenazine apparent clearance CLP/F and on the shared apparent central volume VC/F.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline weight. Reference weight is 70 kg, stated explicitly in Nguyen 2025 Supplemental Table S2 ('centered on a reference weight of 70 kg') for both equations (1) and (2) -- NOT the cohort mean of 77.8 kg. CLP/F = 23.2 * (WT/70)^0.602; VC/F = 226 * (WT/70)^1.04. Because VC/F is shared between valbenazine and [+]-alpha-HTBZ, the weight effect on VC/F propagates to the metabolite central volume as well. No weight effect on QP/F, VPP/F, CLM, QM or VPM. Nguyen 2025 Results: the relative difference in predicted valbenazine AUC across the observed weight range versus the median weight was below 30%.",
      source_name        = "WT"
    ),
    DOSE = list(
      description        = "Administered valbenazine dose for the current dose record, in mg. Drives the saturable (Emax) increase in relative oral bioavailability F1 above the 1 mg anchor dose.",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-record administered dose (use case (a) of the DOSE canonical: dose drives a dose-dependent relative bioavailability applied through f(depot); the same shape as Wada 2023 sparsentan). Nguyen 2025 Supplemental Table S2 equation (4): F1 = [1 + FMAX * (DOSE - 1 mg) / (ED50 + DOSE - 1 mg)] * (1 + theta10 * FED), with F1 defined as 1 for a 1 mg dose in the fasted state. Studied dose range 1-150 mg. The model clamps (DOSE - 1) at 0 so a sub-1-mg dose cannot drive F1 below the anchor; the paper never dosed below 1 mg. At clinically relevant doses the effect is modest: 1.09-fold at 60 mg and 1.16-fold at 80 mg relative to 40 mg (Nguyen 2025 Discussion).",
      source_name        = "DOSE"
    ),
    FED = list(
      description        = "Fed-versus-fasted status of the dose record. Proportional reductions in both the transit absorption rate constant KTR and the relative bioavailability F1.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = "Per-dose-record indicator. Nguyen 2025 Supplemental Table S2 equations (3) and (4), with proportional-shift coefficients theta9 = -0.685 on KTR (68.5% slower absorption when fed) and theta10 = -0.0646 on F1 (6.5% lower bioavailability when fed) from Table 1. Nguyen 2025 Results describes the re-evaluated clinical impact of food status as minimal. In KINECT-HD the capsule was taken 'without regard to food' (Supplemental Table S1), so a KINECT-HD simulation may reasonably use either level; the fasted reference (FED = 0) reproduces the typical-value parameters of Table 1 directly.",
      source_name        = "FED"
    ),
    FORM_SOLUTION = list(
      description        = "Oral-solution formulation indicator for the dose record. Proportional increase in the transit absorption rate constant KTR relative to the capsule reference.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (capsule; the formulation used in both phase 3 studies, KINECT-HD and KINECT 3)",
      notes              = "Per-dose-record indicator. Nguyen 2025 Supplemental Table S2 equation (3): theta8 = 1.51 is the PROPORTIONAL shift for the solution, entering as KTR * (1 + 1.51 * FORM_SOLUTION), i.e. a 2.51-fold higher absorption rate constant for the oral solution than for the capsule -- not a 1.51-fold multiplier. Table S2 states this explicitly ('the proportional shift in KA for the solution formulation (SOLN = 1) relative to capsule formulations (SOLN = 0)'). The oral solution appeared only in the six phase 1 studies; both phase 3 studies used capsules (Supplemental Table S1), so FORM_SOLUTION = 0 for any KINECT-HD or KINECT 3 simulation.",
      source_name        = "SOLN"
    ),
    CYP2D6_PM = list(
      description        = "CYP2D6 poor-metabolizer phenotype indicator. Proportional reduction in the [+]-alpha-HTBZ systemic clearance CLM.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP2D6 extensive or ultra-rapid metabolizer, when CYP2D6_IM is also 0)",
      notes              = "Time-fixed germline-genotype-derived phenotype, determined by PCR-based assay on a baseline whole-blood sample (Nguyen 2025 Bioanalytical Methods). Nguyen 2025 Supplemental Table S2 equation (5): CLM = 31.0 * (1 + theta17 * PM) * (1 + theta18 * IM), theta17 = -0.516, so poor metabolizers have 51.6% lower metabolite clearance and a 2.06-fold higher steady-state [+]-alpha-HTBZ AUC (1 / (1 - 0.516) = 2.07; Nguyen 2025 Results reports AUC ratio 2.06 [90% CI 1.7-2.62] and Cmax ratio 1.83 [90% CI 1.55-2.25]). Paired with CYP2D6_IM so that both indicators = 0 selects the pooled extensive-or-ultra-rapid reference. Cohort: 5.6% poor metabolizers among the 425 pooled subjects.",
      source_name        = "PM"
    ),
    CYP2D6_IM = list(
      description        = "CYP2D6 intermediate-metabolizer phenotype indicator. Proportional reduction in the [+]-alpha-HTBZ systemic clearance CLM.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP2D6 extensive or ultra-rapid metabolizer, when CYP2D6_PM is also 0)",
      notes              = "Time-fixed germline-genotype-derived phenotype. Nguyen 2025 Supplemental Table S2 equation (5): theta18 = -0.281, so intermediate metabolizers have 28.1% lower metabolite clearance and a 1.39-fold higher steady-state [+]-alpha-HTBZ AUC. Paired with CYP2D6_PM; the reference (both indicators = 0) is the POOLED extensive-plus-ultra-rapid group, which is why an intermediate-metabolizer indicator is required rather than reusing CYP2D6_EM (whose 0 level would also capture ultra-rapid metabolizers and would therefore apply the intermediate-metabolizer reduction to them). Table 1 footnote a: subjects whose CYP2D6 status was inconclusive and reported as either intermediate or extensive were assigned to the intermediate-metabolizer category. Cohort: 30.5% intermediate metabolizers among the 425 pooled subjects (58.2% extensive, 3.2% ultra-rapid).",
      source_name        = "IM"
    ),
    SAMPLE_INTENSIVE = list(
      description        = "Per-observation sampling-intensity indicator selecting between the intensive (rich) and sparse log-additive residual-error magnitudes, separately for each analyte.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (sparse sampling)",
      notes              = "Per-observation-record indicator. Nguyen 2025 Methods: 'Residual unexplained variability was investigated for refinement by distinguishing errors between parent and metabolite, as well as between sparse and rich sampling data', and Results reports that stratifying the residual error by sampling design reduced the OFV from -7341.67 to -10154.80. Table 1 thetas 11/12 (parent rich/sparse: 0.402 / 0.678) and 19/20 (metabolite rich/sparse: 0.264 / 0.720). The six phase 1 studies contributed rich profiles (SAMPLE_INTENSIVE = 1); the two phase 3 studies (KINECT-HD, KINECT 3) contributed sparse samples (SAMPLE_INTENSIVE = 0) -- Nguyen 2025 attributes the wider phase 3 residual to 'high residual error from sparse sampling and limited dose timing information'. A single subject may carry both levels.",
      source_name        = "sampling design (rich vs sparse)"
    )
  )

  # Covariates that Nguyen 2025 SCREENED but did not retain in the final PK
  # model. Documented here for provenance; they are deliberately absent from
  # model() and from covariateData.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex. Screened as an intrinsic factor on valbenazine and [+]-alpha-HTBZ PK.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Nguyen 2025 Results: 'valbenazine and [+]-alpha-HTBZ exposure remained unaffected by other intrinsic patient factors, including sex, age, race, and kidney and liver function markers, confirming that dose adjustment is not necessary for these subgroups.' No point estimate reported; not retained."
    ),
    AGE = list(
      description = "Age. Screened as an intrinsic factor on valbenazine and [+]-alpha-HTBZ PK.",
      units       = "years",
      type        = "continuous",
      notes       = "Nguyen 2025 Results: not retained (see SEXF note). Cohort mean (SD) age 48.3 (14.3) years."
    ),
    RACE_WHITE = list(
      description = "Race. Screened as an intrinsic factor on valbenazine and [+]-alpha-HTBZ PK.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Nguyen 2025 Results: not retained (see SEXF note). The KINECT-HD efficacy cohort was 96.0% White (Supplemental Table S5), so race was poorly identifiable in the HD subgroup."
    ),
    CRCL = list(
      description = "Kidney function marker. Screened as an intrinsic factor on valbenazine and [+]-alpha-HTBZ PK.",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Nguyen 2025 Results: 'kidney and liver function markers' were screened and not retained. The paper does not name the specific marker or report a point estimate."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 425L,
    n_studies      = 8L,
    n_observations = 14371L,
    age_range      = "mean (SD) 48.3 (14.3) years across the pooled PK dataset; KINECT-HD enrolled adults 18-75 years",
    weight_range   = "mean (SD) 77.8 (15.7) kg across the pooled PK dataset; KINECT-HD 41.0-134 kg (median 75.8)",
    weight_median  = "75.8 kg (KINECT-HD full analysis set, Supplemental Table S6)",
    sex_female_pct = 43.3,
    race_ethnicity = c(White = 96.0, Black = 0.8, Asian = 0.8, Other = 2.4),
    disease_state  = "Pooled: healthy adults (six phase 1 studies), patients with tardive dyskinesia (KINECT 3, NCT02274558), and patients with Huntington's-disease-associated chorea (KINECT-HD, protocol HD3005)",
    dose_range     = "1-150 mg valbenazine orally in the phase 1 studies; 40-80 mg once daily in KINECT-HD; 40 or 80 mg once daily in KINECT 3",
    cyp2d6_status  = c(extensive = 58.2, intermediate = 30.5, poor = 5.6, ultrarapid = 3.2),
    notes          = paste(
      "Pooled PK analysis dataset: 7279 valbenazine and 7092",
      "[+]-alpha-HTBZ plasma concentration records from 425 patients and",
      "healthy participants (Nguyen 2025 Results). Demographics from the",
      "Results 'Population PK Model Development' paragraph; KINECT-HD",
      "baseline characteristics from Supplemental Tables S5 and S6. The",
      "model was built stepwise: phase 1 parent data -> joint",
      "parent-metabolite -> plus KINECT 3 (tardive dyskinesia) -> plus",
      "KINECT-HD (Huntington's chorea); this file encodes the FINAL model",
      "(Table 1). Individual empirical Bayes estimates differ by cohort",
      "(Supplemental Table S3): median CLP/F 22.5 / 26.0 / 22.7 L/h and",
      "median CLM 28.5 / 28.2 / 22.2 L/h in healthy subjects / tardive",
      "dyskinesia / Huntington's chorea respectively.",
      sep = " "
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Valbenazine (parent) structural parameters -- Nguyen 2025 Table 1.
    # Apparent (/F) parameters: the model has no absolute bioavailability, so
    # F1 is a RELATIVE bioavailability anchored at 1 for a 1 mg fasted dose.
    # ---------------------------------------------------------------------
    lcl  <- log(23.2); label("Valbenazine apparent systemic clearance CLP/F at 70 kg (L/h)")                  # Table 1 theta 1: 23.2 (%RSE 3.95; 95% CI 21.4, 25.0)
    lvc  <- log(226);  label("Apparent central volume VC/F at 70 kg, SHARED by valbenazine and [+]-alpha-HTBZ (L)") # Table 1 theta 2: 226 (%RSE 4.49; 95% CI 206, 246)
    lktr <- log(8.43); label("Transit absorption rate constant KTR = KA, capsule in the fasted state (1/h)")  # Table 1 theta 3: 8.43 (%RSE 3.84; 95% CI 7.79, 9.06)
    lq   <- log(22.0); label("Valbenazine apparent inter-compartmental clearance QP/F (L/h)")                 # Table 1 theta 4: 22.0 (%RSE 4.35; 95% CI 20.1, 23.9)
    lvp  <- log(198);  label("Valbenazine apparent peripheral volume VPP/F (L)")                              # Table 1 theta 5: 198 (%RSE 3.44; 95% CI 185, 211)

    # ---------------------------------------------------------------------
    # Saturable dose-dependent relative bioavailability -- Table S2 eq (4):
    #   F1 = [1 + FMAX * (DOSE - 1) / (ED50 + DOSE - 1)] * (1 + theta10 * FED)
    # Anchored so that F1 = 1 for a 1 mg dose in the fasted state.
    # ---------------------------------------------------------------------
    lfmax <- log(1.36); label("Maximum increase in relative oral bioavailability above the 1 mg anchor dose, FMAX (unitless)") # Table 1 theta 6: 1.36 (%RSE 10.9; 95% CI 1.07, 1.65)
    led50 <- log(90.1); label("Valbenazine dose achieving 50% of FMAX, ED50 (mg)")                                            # Table 1 theta 7: 90.1 (%RSE 22.9; 95% CI 49.6, 131)

    # ---------------------------------------------------------------------
    # Extrinsic covariate effects on absorption -- Table S2 eqs (3) and (4).
    # All three are PROPORTIONAL shifts, entering as (1 + theta * COV).
    # ---------------------------------------------------------------------
    e_soln_ktr   <-  1.51;    label("Proportional shift in KTR for the oral solution vs capsule (unitless)")   # Table 1 theta 8:  1.51    (%RSE 5.59;  95% CI 1.35, 1.68)
    e_fed_ktr    <- -0.685;   label("Proportional shift in KTR for the fed vs fasted state (unitless)")        # Table 1 theta 9:  -0.685  (%RSE 0.816; 95% CI -0.696, -0.674)
    e_fed_fdepot <- -0.0646;  label("Proportional shift in F1 for the fed vs fasted state (unitless)")         # Table 1 theta 10: -0.0646 (%RSE 20.2;  95% CI -0.0901, -0.0391)

    # ---------------------------------------------------------------------
    # [+]-alpha-HTBZ (metabolite) structural parameters -- Table 1.
    # CLM / QM / VPM carry no /F because the parent's F1 already scales the
    # amount of valbenazine that reaches the systemic circulation.
    # ---------------------------------------------------------------------
    lcl_htbz <- log(31.0); label("[+]-alpha-HTBZ systemic clearance CLM in extensive or ultra-rapid CYP2D6 metabolizers (L/h)") # Table 1 theta 13: 31.0 (%RSE 4.99; 95% CI 28.0, 34.1)
    lq_htbz  <- log(1.27); label("[+]-alpha-HTBZ inter-compartmental clearance QM (L/h)")                                       # Table 1 theta 15: 1.27 (%RSE 7.60; 95% CI 1.08, 1.45)
    lvp_htbz <- log(97.8); label("[+]-alpha-HTBZ peripheral volume VPM (L)")                                                    # Table 1 theta 16: 97.8 (%RSE 7.12; 95% CI 84.2, 111)

    # Fraction of valbenazine metabolised to [+]-alpha-HTBZ. Estimated on the
    # natural scale under the shared-VC/F constraint that makes it
    # identifiable; kept on the natural scale because the paper reports it
    # that way and the 95% CI (0.199, 0.214) sits far inside (0, 1).
    fm <- 0.207; label("Fraction of valbenazine metabolised to [+]-alpha-HTBZ, FM (molar fraction)")  # Table 1 theta 14: 0.207 (%RSE 1.91; 95% CI 0.199, 0.214). Discussion: previously FIXED at 0.35 in the predecessor model; re-estimated here.

    # CYP2D6 phenotype effects on CLM -- Table S2 eq (5), proportional shifts
    # relative to the POOLED extensive-plus-ultra-rapid reference.
    e_cyp2d6_pm_cl_htbz <- -0.516; label("Proportional shift in CLM for CYP2D6 poor metabolizers (unitless)")         # Table 1 theta 17: -0.516 (%RSE 12.3; 95% CI -0.641, -0.391) -> 51.6% lower CLM
    e_cyp2d6_im_cl_htbz <- -0.281; label("Proportional shift in CLM for CYP2D6 intermediate metabolizers (unitless)") # Table 1 theta 18: -0.281 (%RSE 15.3; 95% CI -0.366, -0.197) -> 28.1% lower CLM

    # ---------------------------------------------------------------------
    # Allometric body-weight exponents -- Table S2 eqs (1) and (2),
    # centered on a 70 kg reference weight (stated in Table S2).
    # ---------------------------------------------------------------------
    e_wt_cl <- 0.602; label("Power exponent on (WT/70) for valbenazine CLP/F (unitless)") # Table 1 theta 21: 0.602 (%RSE 21.2; 95% CI 0.352, 0.852)
    e_wt_vc <- 1.04;  label("Power exponent on (WT/70) for the shared VC/F (unitless)")   # Table 1 theta 22: 1.04  (%RSE 17.0; 95% CI 0.695, 1.39)

    # ---------------------------------------------------------------------
    # Inter-individual variability -- Table 1, reported as omega^2 on the
    # log scale. Table 1 reports no off-diagonal covariances, so the omega
    # matrix is encoded as diagonal.
    # ---------------------------------------------------------------------
    etalcl      ~ 0.150  # Table 1 IIV CLP/F: omega^2 = 0.150 (shrinkage 16.5%; 95% CI 0.116, 0.184) -> 38.7% CV = sqrt(exp(0.150)-1)
    etalvc      ~ 0.194  # Table 1 IIV VC/F:  omega^2 = 0.194 (shrinkage 35.8%; 95% CI 0.0849, 0.303) -> 44.0% CV
    etalktr     ~ 0.187  # Table 1 IIV KA:    omega^2 = 0.187 (shrinkage 40.4%; 95% CI 0.138, 0.236) -> 43.2% CV
    etalcl_htbz ~ 0.225  # Table 1 IIV CLM:   omega^2 = 0.225 (shrinkage 15.2%; 95% CI 0.126, 0.325) -> 47.5% CV

    # ---------------------------------------------------------------------
    # Residual error -- Table 1 thetas 11, 12, 19, 20. NONMEM additive error
    # on the natural-log-transformed molar concentration, i.e. log-normal in
    # linear space, so these map onto expSd / lnorm(). The %-transformed
    # estimates Table 1 reports (40.2%, 67.8%, 26.4%, 72.0%) are the same
    # numbers read as approximate CVs.
    #
    # Four magnitudes: {parent, metabolite} x {intensive, sparse} sampling.
    # A single rxode2 endpoint takes one error parameter, so the per-record
    # magnitude is assembled in model() from SAMPLE_INTENSIVE and handed to
    # lnorm() as an expression -- the Friberg_2012_voriconazole.R pattern.
    # ---------------------------------------------------------------------
    expSdIntensive      <- 0.402; label("Log-scale residual SD for valbenazine, intensive (rich) sampling")     # Table 1 theta 11: 0.402 (%RSE 0.994; 95% CI 0.395, 0.410) -> 40.2%
    expSdSparse         <- 0.678; label("Log-scale residual SD for valbenazine, sparse sampling")               # Table 1 theta 12: 0.678 (%RSE 2.22;  95% CI 0.649, 0.708) -> 67.8%
    expSdIntensive_htbz <- 0.264; label("Log-scale residual SD for [+]-alpha-HTBZ, intensive (rich) sampling")  # Table 1 theta 19: 0.264 (%RSE 1.04;  95% CI 0.259, 0.270) -> 26.4%
    expSdSparse_htbz    <- 0.720; label("Log-scale residual SD for [+]-alpha-HTBZ, sparse sampling")            # Table 1 theta 20: 0.720 (%RSE 2.28;  95% CI 0.688, 0.752) -> 72.0%
  })

  model({
    # ---------------------------------------------------------------------
    # Constants. Nguyen 2025 fitted natural-log-transformed concentrations in
    # MOLAR units, so FM = 0.207 is a molar fraction. This file carries drug
    # amounts in mg and reports concentrations in ng/mL (the units of the
    # published exposure summaries in Table 2 and of the assay LLOQs), so the
    # metabolite-formation flux needs the molar-to-mass conversion
    # MW([+]-alpha-HTBZ) / MW(valbenazine). The two molecular weights are
    # standard chemical constants and are NOT reported in Nguyen 2025:
    #   valbenazine free base   C24H38N2O4, MW 418.57 g/mol (PubChem CID 24756910)
    #   [+]-alpha-HTBZ          C19H29NO3,  MW 319.44 g/mol (PubChem CID 92987)
    # Valbenazine doses are expressed as free base (the marketed capsule
    # contains valbenazine tosylate equivalent to the labelled free-base
    # amount), so 418.57 is the right denominator for a mg dose.
    # ---------------------------------------------------------------------
    mwValbenazine <- 418.57   # g/mol
    mwHtbz        <- 319.44   # g/mol
    refWt         <- 70       # kg -- Table S2 reference weight

    # ---------------------------------------------------------------------
    # Individual valbenazine parameters -- Table S2 eqs (1), (2), (3).
    # ---------------------------------------------------------------------
    cl  <- exp(lcl  + etalcl)  * (WT / refWt)^e_wt_cl
    vc  <- exp(lvc  + etalvc)  * (WT / refWt)^e_wt_vc
    ktr <- exp(lktr + etalktr) * (1 + e_soln_ktr * FORM_SOLUTION) * (1 + e_fed_ktr * FED)
    q   <- exp(lq)
    vp  <- exp(lvp)

    # ---------------------------------------------------------------------
    # Individual [+]-alpha-HTBZ parameters -- Table S2 eq (5). The metabolite
    # central volume is the SAME vc as the parent's (the identifiability
    # constraint of Nguyen 2025 Methods), so it inherits both the allometric
    # weight effect and etalvc.
    # ---------------------------------------------------------------------
    cl_htbz <- exp(lcl_htbz + etalcl_htbz) *
               (1 + e_cyp2d6_pm_cl_htbz * CYP2D6_PM) *
               (1 + e_cyp2d6_im_cl_htbz * CYP2D6_IM)
    q_htbz  <- exp(lq_htbz)
    vp_htbz <- exp(lvp_htbz)

    # ---------------------------------------------------------------------
    # Saturable dose-dependent relative bioavailability -- Table S2 eq (4).
    # max(DOSE - 1, 0) is a defensive clamp: the anchor is the 1 mg dose and
    # the studied range starts there, so a sub-1-mg dose must not push F1
    # below the anchor value of 1.
    # ---------------------------------------------------------------------
    fmax   <- exp(lfmax)
    ed50   <- exp(led50)
    dosexs <- max(DOSE - 1, 0)
    fdepot <- (1 + fmax * dosexs / (ed50 + dosexs)) * (1 + e_fed_fdepot * FED)

    # ---------------------------------------------------------------------
    # ODE system.
    # Absorption: a chain of four transit compartments between the depot and
    # the central compartment, every transfer at the single rate constant
    # KTR = KA (Nguyen 2025 Results: 'A total of four transit absorption
    # compartments were identified'; Table S2 eq (3): KA_i = KTR_i).
    # Disposition: two-compartment linear for each analyte.
    # Metabolism: a fraction FM of the valbenazine elimination flux appears
    # as [+]-alpha-HTBZ in the metabolite central compartment, converted from
    # valbenazine mass to metabolite mass by the molecular-weight ratio.
    # Nguyen 2025 Methods: the biotransformation 'was assumed to follow a
    # first-order process, occurring after valbenazine became available
    # systemically'.
    # ---------------------------------------------------------------------
    d/dt(depot)    <- -ktr * depot
    d/dt(transit1) <-  ktr * depot    - ktr * transit1
    d/dt(transit2) <-  ktr * transit1 - ktr * transit2
    d/dt(transit3) <-  ktr * transit2 - ktr * transit3
    d/dt(transit4) <-  ktr * transit3 - ktr * transit4

    d/dt(central)     <-  ktr * transit4 -
                          cl * central / vc -
                          q  * central / vc + q * peripheral1 / vp
    d/dt(peripheral1) <-  q  * central / vc - q * peripheral1 / vp

    d/dt(central_htbz)     <-  fm * (mwHtbz / mwValbenazine) * cl * central / vc -
                               cl_htbz * central_htbz / vc -
                               q_htbz  * central_htbz / vc + q_htbz * peripheral1_htbz / vp_htbz
    d/dt(peripheral1_htbz) <-  q_htbz  * central_htbz / vc - q_htbz * peripheral1_htbz / vp_htbz

    f(depot) <- fdepot

    # ---------------------------------------------------------------------
    # Observations. Amounts are in mg and volumes in L, so amount/volume is
    # mg/L; x1000 converts to ng/mL. Both analytes divide by the SHARED vc.
    # ---------------------------------------------------------------------
    Cc      <- 1000 * central      / vc
    Cc_htbz <- 1000 * central_htbz / vc

    # Per-record residual magnitude selected by sampling design.
    expSdNow      <- expSdIntensive      * SAMPLE_INTENSIVE + expSdSparse      * (1 - SAMPLE_INTENSIVE)
    expSdNow_htbz <- expSdIntensive_htbz * SAMPLE_INTENSIVE + expSdSparse_htbz * (1 - SAMPLE_INTENSIVE)

    Cc      ~ lnorm(expSdNow)
    Cc_htbz ~ lnorm(expSdNow_htbz)
  })
}
