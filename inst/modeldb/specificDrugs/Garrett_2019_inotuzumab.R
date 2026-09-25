Garrett_2019_inotuzumab <- function() {
  description <- "Two-compartment population PK model for inotuzumab ozogamicin in adults with relapsed/refractory B-cell acute lymphoblastic leukemia (ALL) or B-cell non-Hodgkin lymphoma (NHL); linear plus empirical time-dependent (target-mediated) clearance with baseline body surface area on CL1, CL2 and V1, baseline percentage of peripheral-blood blasts on the time-dependent decay coefficient, and concomitant rituximab on CL1 (Garrett 2019, 11 pooled adult studies)."
  reference <- "Garrett M, Ruiz-Garcia A, Parivar K, Hee B, Boni J. Population pharmacokinetics of inotuzumab ozogamicin in relapsed/refractory acute lymphoblastic leukemia and non-Hodgkin lymphoma. J Pharmacokinet Pharmacodyn. 2019;46(3):211-222. doi:10.1007/s10928-018-9614-9"
  vignette <- "Garrett_2019_inotuzumab"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Garrett 2019 Methods ('Pharmacokinetic
  # sampling and bioanalytical methods': serum concentrations of InO measured by
  # HPLC/MS/MS in ALL studies and by ELISA in NHL studies) and Online Resource 3
  # ($MODEL NCOMP=2 COMP=(CENTRAL) COMP=(PHP), amounts in mg with S1 = V1 in L).
  compartmentData <- list(
    central = list(analyte = "inotuzumab ozogamicin", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "inotuzumab ozogamicin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    BSA_BASE = list(
      description = "Baseline body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed baseline. Enters CL1 and CL2 as power models centered on the population median 1.84 m^2 (exponents 1.54 and 1.64) and Vc as a LINEAR term in the centered deviation, (1 + 0.774 * (BSA_BASE - 1.84)) -- not a power model. The linear form is what Online Resource 3 implements (V1BBSA = (1 + THETA(15) * (BBSA - 1.84))) and is confirmed by the paper's own reported effect sizes: at BSA 1.55 m^2 (10th percentile) and 2.21 m^2 (90th percentile) Garrett 2019 reports Vc changes of -23% and +29%, which the linear form reproduces (-22.4%, +28.6%) and a power form does not (-12.4%, +15.2%). Source column 'BBSA'.",
      source_name = "BBSA"
    ),
    BLSTPB = list(
      description = "Baseline percentage of blasts in peripheral blood",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed baseline. Power effect on cl_exp_kdes (exponent -0.0401) referenced to 5.25%, applied to ALL patients only. Percentage of blasts was not collected in the NHL studies (Garrett 2019 Table 1 shows NA for NHL); Online Resource 3 encodes this as a missing-data sentinel, IF(BLSTPB.EQ.-99) KDESBLSTPB = 1. This model reproduces that branch without a sentinel by multiplying the exponent by DIS_BCPALL, so the power term collapses to exactly 1 for NHL subjects regardless of the BLSTPB value supplied. Distinct from BLSTABL (absolute counts, 10^9/L), which is the column the successor Wu 2024 model uses. Source column 'BLSTPB'.",
      source_name = "BLSTPB"
    ),
    DIS_BCPALL = list(
      description = "B-cell precursor acute lymphoblastic leukemia disease-state indicator (1 = B-cell ALL, 0 = B-cell non-Hodgkin lymphoma)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (B-cell non-Hodgkin lymphoma)",
      notes = "Garrett 2019 calls this the 'ALL effect' and states that it confounds disease type (B-cell ALL vs B-cell NHL) with bioanalytical assay method (HPLC/MS/MS for ALL studies, ELISA for NHL studies), because each method was used exclusively in one tumour type (Table 2 footnote d). Fractional-change (dummy-variable) effects on CL1 (-0.745) and cl_exp_kdes (-0.860), and gates the BLSTPB effect on cl_exp_kdes. Source column 'PTST' (coded 1 = NHL, 2 = ALL in Online Resource 3; recoded here to the canonical 0/1 indicator, which leaves both coefficients verbatim because NHL is the reference in both codings). The main-text equations call the same indicator 'PTSTALL'.",
      source_name = "PTST"
    ),
    CONMED_RITUX = list(
      description = "Concomitant rituximab combination-therapy indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (concomitant rituximab; the most common category in the pooled dataset)",
      notes = "Reference category is WITH rituximab: Online Resource 3 sets CL1RITX = 1 when RITX = 1 (commented '; Most common') and CL1RITX = (1 + THETA(12)) when RITX = 0, so the +0.155 fractional increase in CL1 applies to subjects NOT receiving rituximab. Garrett 2019 Results states this directly: 'the absence of concomitant rituximab use resulted in an estimated increase of CL1 by 16% (95% CI 5-26%) in patients receiving single-agent InO versus those receiving rituximab plus InO'. The source column 'RITX' therefore carries the OPPOSITE polarity to the canonical CONMED_RITUX (1 = receiving rituximab), and the model applies the effect as (1 + 0.155 * (1 - CONMED_RITUX)) so that the published coefficient is preserved verbatim. NOTE the main-text equation block is self-contradictory on this point: it glosses 'RITX' as '(without rituximab) ... 1 if applicable', the reverse of the control stream. The semantics are nonetheless unambiguous -- all three sources agree that no-rituximab raises CL1 by 15.5% -- and the ALL-arm equation confirms it numerically (0.113 * (1 - 0.745) * (1 + 0.155) = 0.0333 L/h, the published typical CL1 for single-agent ALL patients). Source column 'RITX'.",
      source_name = "RITX"
    )
  )

  # Covariates screened during stepwise covariate modelling but NOT retained in
  # the final model (Garrett 2019 'Covariate model development'). Documented for
  # provenance only; none is referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Baseline age",
      units = "years",
      type = "continuous",
      notes = "Screened; Garrett 2019 reports 'Age, race, and baseline creatinine clearance did not have a significant effect on InO PK parameters'. Retained by the successor Wu 2024 pediatric reanalysis on cl_exp_kdes, but not by this adult model."
    ),
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened; superseded by BSA_BASE. Garrett 2019: 'Including BBSA as a covariate also corrected any previous trends associated with body weight (Online Resource 7)'."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Not tested. Garrett 2019: 'All potential covariates were tested in the M3 base model with the exception of sex because it was highly correlated with BBSA'."
    ),
    CRCL = list(
      description = "Baseline creatinine clearance",
      units = "mL/min",
      type = "continuous",
      notes = "Screened and not significant (Garrett 2019 'Covariate model development')."
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/dL",
      type = "continuous",
      notes = "Selected by generalized additive modelling for CL1 and CL2 but eliminated in backward selection at P < 0.001; not in the final model."
    ),
    ALT = list(
      description = "Baseline alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened; not retained."
    ),
    AST = list(
      description = "Baseline aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened; not retained."
    ),
    BILI = list(
      description = "Baseline total bilirubin",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened; not retained. Hepatic impairment by NCI ODWG category was also tested and not retained (Online Resource 4 tabulates the categories)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 765L,
    n_studies = 11L,
    n_observations = 8361L,
    cohorts = list(
      ALL = list(
        n = 234L,
        n_studies = 2L,
        age_median = "46.0 years (range 20.0-79.0)",
        bsa_median = "1.86 m^2 (range 1.27-2.81)",
        weight_median = "74.0 kg (range 30.9-154)",
        regimen = "Single-agent InO in all ALL studies"
      ),
      NHL = list(
        n = 531L,
        n_studies = 9L,
        age_median = "65.0 years (range 18.0-92.0)",
        bsa_median = "1.83 m^2 (range 1.13-2.56)",
        weight_median = "73.2 kg (range 33.5-148)",
        regimen = "3 studies single-agent InO, 5 studies InO plus rituximab, 1 study InO plus rituximab and chemotherapy"
      )
    ),
    age_range = "18-92 years",
    age_median = "61.0 years (Garrett 2019 Table 1, total column)",
    weight_range = "30.9-154 kg",
    weight_median = "73.3 kg (Garrett 2019 Table 1)",
    bsa_range = "1.13-2.81 m^2",
    bsa_median = "1.84 m^2 (Garrett 2019 Table 1; the reference value used in every BSA covariate term)",
    sex_female_pct = 40.1,
    race_ethnicity = c(White = 69.8, Black = 2.61, Asian = 7.06, Japanese = 13.2, Other = 6.93, Unknown = 0.392),
    disease_state = "Relapsed or refractory CD22+ B-cell acute lymphoblastic leukemia (n = 234) or relapsed or refractory B-cell non-Hodgkin lymphoma (n = 531).",
    dose_range = "ALL: 1.2-1.8 mg/m^2 per 21- to 28-day cycle given as 2 or 3 fractionated weekly doses (the pivotal regimen is 1.8 mg/m^2/cycle as 0.8 mg/m^2 on day 1 and 0.5 mg/m^2 on days 8 and 15). NHL: single doses of 0.4-2.4 mg/m^2 on day 1 or 2 of each cycle. All doses by intravenous infusion.",
    regions = "Multi-regional (global phase 3 study 1022 plus Japanese and other regional studies; 13.2% of the pooled population was Japanese).",
    renal_function = "Baseline creatinine clearance median 93.1 mL/min (range 18.2-368); ALL 122 vs NHL 81.8 mL/min.",
    hepatic_function = "NCI ODWG category A (normal) in 79.7%, B1/B2 (mild) in 19.6%, C (moderate) in 0.392%, D (severe) in 0.131% (Online Resource 4).",
    bioanalytic_methods = "ALL studies: validated HPLC with tandem mass spectrometry indirectly measuring N-acetyl-gamma-calicheamicin dimethyl hydrazide conjugated to the InO antibody, LLOQ 1.0 ng/mL. NHL studies: validated ELISA directly measuring conjugated N-acetyl-gamma-calicheamicin, LLOQ 50-667 ng/mL depending on study. Two residual-error magnitudes were estimated to absorb the assay difference (Garrett 2019 Table 2 footnote f).",
    notes = "Final M3 population PK model (data below LLOQ treated as censored via the F_FLAG likelihood, NONMEM 7 levels 2.0/3.0, FOCE-I with Laplacian; Perl-speaks-NONMEM 4.2.0 for the 1000-sample nonparametric bootstrap and pvcVPCs). Baseline demographics are Garrett 2019 Table 1 (continuous) and Online Resource 4 (categorical); the study list is Online Resource 1. Note a paper-internal inconsistency in the observation count: the Results text reports 8361 serum PK samples in total but also 2978 from ALL and 6272 from NHL patients, which sum to 9250; n_observations records the 8361 total actually stated for the analysis dataset. Predecessor of the Wu 2024 pooled adult + pediatric reanalysis; see modellib('Wu_2024_inotuzumab')."
  )

  ini({
    # Structural parameters: final-model typical values, Garrett 2019 Table 2
    # ('NONMEM results OFV = 1450.357', Estimate column). Time unit is hour;
    # clearances L/h, volumes L, cl_exp_kdes 1/h. Reference covariate values are
    # BSA_BASE = 1.84 m^2, BLSTPB = 5.25%, DIS_BCPALL = 0 (NHL),
    # CONMED_RITUX = 1 (with rituximab).
    #
    # NOTE the Table 2 / Online Resource 3 values are the finals; the $THETA
    # initial estimates in Online Resource 3 (0.118248, 6.63858, 0.367849,
    # 0.0321833, 0.0398226, 5.46741) are seeded near but not at the finals and
    # are NOT used here.
    #
    # Garrett 2019 names CL1 the 'linear clearance' and CL2 the 'clearance
    # associated with time-dependent clearance'; total clearance is
    # CL = CL1 + CL2 * exp(-kdes * time). Mapped onto the nlmixr2lib
    # exponential-decay-to-a-constant family: CL1 -> cl_exp_inf (the
    # non-decaying asymptote), CL2 -> cl_exp_component, kdes -> cl_exp_kdes.
    lcl_exp_inf <- log(0.113); label("Linear (asymptotic) clearance for an NHL adult receiving rituximab (L/h)") # Garrett 2019 Table 2 (CL1)
    lvc <- log(6.70); label("Central volume of distribution (L)") # Garrett 2019 Table 2 (V1)
    lcl_exp_component <- log(0.369); label("Initial value of the time-dependent clearance component (L/h)") # Garrett 2019 Table 2 (CL2)
    lcl_exp_kdes <- log(0.0337); label("Decay coefficient of the time-dependent clearance for an NHL adult (1/h)") # Garrett 2019 Table 2 (kdes)
    lq <- log(0.0405); label("Intercompartmental clearance (L/h)") # Garrett 2019 Table 2 (Q)
    lvp <- log(5.10); label("Peripheral volume of distribution (L)") # Garrett 2019 Table 2 (V2)

    # Covariate effects, Garrett 2019 Table 2. Continuous covariates on CL1,
    # CL2 and cl_exp_kdes are power models centered on the reference value;
    # the BSA effect on Vc is LINEAR in the centered deviation; the categorical
    # effects are additive fractional-change dummy variables (1 + theta * I).
    e_bsa_cl_exp_inf <- 1.54; label("Power exponent of baseline BSA on cl_exp_inf (unitless)") # Garrett 2019 Table 2 ('CL1 / BBSA effect')
    e_bsa_cl_exp_component <- 1.64; label("Power exponent of baseline BSA on cl_exp_component (unitless)") # Garrett 2019 Table 2 ('CL2 / BBSA effect')
    e_bsa_vc <- 0.774; label("Linear coefficient of centered baseline BSA on Vc (1/m^2)") # Garrett 2019 Table 2 ('V1 / BBSA effect'); linear form per Online Resource 3 V1BBSA
    e_all_cl_exp_inf <- -0.745; label("Fractional change in cl_exp_inf for B-cell ALL vs NHL (unitless)") # Garrett 2019 Table 2 ('CL1 / ALL effect')
    e_all_cl_exp_kdes <- -0.860; label("Fractional change in cl_exp_kdes for B-cell ALL vs NHL (unitless)") # Garrett 2019 Table 2 ('kdes / ALL effect')
    e_blstpb_cl_exp_kdes <- -0.0401; label("Power exponent of baseline peripheral-blood blast percentage on cl_exp_kdes, ALL only (unitless)") # Garrett 2019 Table 2 ('kdes / BLSTPB effect')
    e_ritux_cl_exp_inf <- 0.155; label("Fractional change in cl_exp_inf in the ABSENCE of concomitant rituximab (unitless)") # Garrett 2019 Table 2 ('CL1 / RITX + InO effect')

    # Inter-individual variability. Garrett 2019 Table 2 reports the IIV
    # 'Estimate' on the percent-CV scale (footnote c: 'interindividual
    # variability of parameter estimates has been reported on the %CV scale,
    # i.e. sqrt(x^2) as the parameters follow a log-normal distribution'), while
    # the accompanying 95% CI columns are on the VARIANCE scale. The variances
    # below are CV^2 and each is confirmed by the midpoint of its own published
    # CI: CL1 0.423^2 = 0.1789 vs CI 0.147-0.211 (midpoint 0.179); V1
    # 0.412^2 = 0.1697 vs 0.150-0.190 (0.170); CL2 0.672^2 = 0.4516 vs
    # 0.370-0.533 (0.4515); kdes 0.455^2 = 0.2070 vs 0.137-0.277 (0.207).
    # The three off-diagonal covariances are printed directly in Table 2 and
    # independently reproduce footnote e's published correlations from these
    # variances: 0.156/sqrt(0.1789*0.1697) = 89.5% ('89.4% for CL1-V1'),
    # 0.213/sqrt(0.1789*0.4516) = 74.9% ('75.0% for CL1-CL2'), and
    # 0.222/sqrt(0.4516*0.1697) = 80.2% ('80.2% for CL2-V1').
    # Block order follows Online Resource 3 $OMEGA BLOCK(3): ETA(1) CL1,
    # ETA(2) V1, ETA(3) CL2; kdes is a separate BLOCK(1). Q and V2 carry no IIV
    # (Garrett 2019: 'the time-dependent model was further improved by removing
    # the random effect on peripheral compartment parameters').
    etalcl_exp_inf + etalvc + etalcl_exp_component ~ c(
      0.178929,
      0.156, 0.169744,
      0.213, 0.222, 0.451584
    ) # Garrett 2019 Table 2 (CV% 42.3 / 41.2 / 67.2; covariances 0.156 / 0.213 / 0.222)
    etalcl_exp_kdes ~ 0.207025 # Garrett 2019 Table 2 (CV% 45.5 for kdes; CI 0.137-0.277 confirms the variance)

    # Residual error. Online Resource 3 fixes $SIGMA 1 FIX and estimates the
    # residual magnitude as a THETA used directly as the standard deviation
    # ('SIG = THETA(7)', 'W = SIG', 'Y = IPRED + W*ERR(1)' with
    # 'IPRED = LOG(F)'). The reported values are therefore residual SDs on the
    # log scale, NOT variances, despite Table 2 labelling the rows 'sigma^2
    # prop' / 'Variance of the ... population' -- that label is a
    # mis-transcription contradicted by the deposited control stream, and the
    # SD reading is what the successor Wu 2024 refit also reports (0.444 NHL /
    # 0.612 adult ALL). Two magnitudes were estimated for the two
    # disease/assay strata (Table 2 footnote f): NHL (ELISA) 0.453, ALL
    # (HPLC/MS/MS) 0.619. The packaged value is the ALL stratum, which is the
    # licensed indication and the arm Garrett 2019 simulates; to simulate the
    # NHL stratum override with `ini(model, expSd = 0.453)` (see vignette).
    expSd <- 0.619; label("Log-scale residual SD for the B-cell ALL / HPLC-MS-MS stratum (unitless)") # Garrett 2019 Table 2 ('sigma^2 prop | ALL' = 0.619, read as an SD per Online Resource 3 $SIGMA 1 FIX); NHL stratum 0.453
  })

  model({
    # Individual PK parameters. Garrett 2019 'Final model results' prints the
    # typical-value equations for an NHL patient as:
    #   CL1  = 0.113 L/h * (1 - 0.745 * PTSTALL) * (BBSA/1.84)^1.54
    #                    * (1 + 0.155 * RITX)
    #   CL2  = 0.369 L/h * (BBSA/1.84)^1.64
    #   V1   = 6.70 L    * (1 + 0.774 * [BBSA - 1.84])
    #   kdes = 0.0337 /h * (1 - 0.860 * PTSTALL) * (BLSTPB/5.25)^-0.0401
    # and the corresponding ALL equations (single-agent, so no rituximab) as
    #   CL1 = 0.0333 L/h * (BBSA/1.84)^1.54  and  kdes = 0.00472 /h *
    #   (BLSTPB/5.25)^-0.0401,
    # which this parameterization reproduces exactly:
    # 0.113 * (1 - 0.745) * (1 + 0.155) = 0.03328 and
    # 0.0337 * (1 - 0.860) = 0.004718.
    #
    # The rituximab term is written on the canonical CONMED_RITUX polarity
    # (1 = receiving rituximab, the reference) as (1 - CONMED_RITUX) so that
    # the published +0.155 coefficient is carried verbatim; see covariateData.
    cl_exp_inf <- exp(lcl_exp_inf + etalcl_exp_inf) *
      (1 + e_all_cl_exp_inf * DIS_BCPALL) *
      (BSA_BASE / 1.84)^e_bsa_cl_exp_inf *
      (1 + e_ritux_cl_exp_inf * (1 - CONMED_RITUX))

    # Vc takes the BSA effect as a LINEAR term in the centered deviation, not a
    # power term (Online Resource 3: V1BBSA = (1 + THETA(15) * (BBSA - 1.84))).
    vc <- exp(lvc + etalvc) * (1 + e_bsa_vc * (BSA_BASE - 1.84))

    cl_exp_component <- exp(lcl_exp_component + etalcl_exp_component) *
      (BSA_BASE / 1.84)^e_bsa_cl_exp_component

    # The BLSTPB effect applies to ALL patients only: peripheral-blood blast
    # percentage was not collected in the NHL studies, and Online Resource 3
    # neutralizes the term for those subjects with a missing-data branch
    # (IF(BLSTPB.EQ.-99) KDESBLSTPB = 1). Multiplying the exponent by
    # DIS_BCPALL reproduces that branch without requiring a -99 sentinel: for
    # an NHL subject the power term is BLSTPB^0 = 1 for any supplied value.
    cl_exp_kdes <- exp(lcl_exp_kdes + etalcl_exp_kdes) *
      (1 + e_all_cl_exp_kdes * DIS_BCPALL) *
      (BLSTPB / 5.25)^(e_blstpb_cl_exp_kdes * DIS_BCPALL)

    q <- exp(lq)
    vp <- exp(lvp)

    # Empirical time-dependent clearance (Garrett 2019 Methods, 'Base
    # pharmacokinetic model development'): CLt = CL2 * exp(-kdes * Time) and
    # CL = CL1 + CLt. In Online Resource 3 the decay is driven by the $DES
    # integration time T on a dataset whose TIME column is time after first
    # dose (the $DATA file is '..._tafd.csv'), so rxode2's `time` is the
    # matching driver: event tables for this model must carry time from the
    # first dose.
    cl_exp_component_t <- cl_exp_component * exp(-cl_exp_kdes * time)
    cl <- cl_exp_inf + cl_exp_component_t

    # Two-compartment model with intravenous input; no depot. Online Resource 3
    # uses ADVAN6 with K10 = CL/V1, K12 = Q/V1, K21 = Q/V2 and S1 = V1.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg and vc in L give central/vc in mg/L = ug/mL; multiply by 1000
    # to report Cc in ng/mL, the scale of the paper's LLOQ (1.0 ng/mL for the
    # ALL assay) and of its simulated AUCtau (29,800 ng*h/mL).
    Cc <- (central / vc) * 1000

    # Residual error is additive on log-transformed data (Garrett 2019: 'the
    # residual error ... for observations was modeled additively based on
    # log-transformed data'; Online Resource 3 Y = IPRED + W*ERR(1) with
    # IPRED = LOG(F)), i.e. log-normal in nlmixr2's linear space.
    Cc ~ lnorm(expSd)
  })
}
