Kirubakaran_2022_tacrolimus <- function() {
  description <- paste0(
    "Two-compartment population pharmacokinetic model for oral immediate-release ",
    "tacrolimus (Prograf) in adult heart transplant recipients (Kirubakaran 2022): ",
    "first-order absorption; FFM-allometric scaling on CL/F and Q/F (exponent 0.75) ",
    "and on V2/F and V3/F (exponent 1.0); haematocrit power effect on CL/F; and a ",
    "state-dependent typical CL/F (without vs with concomitant azole antifungal, ",
    "21.1 vs 4.2 L/h) with a state-dependent CL/F BSV magnitude (61% vs 89.5% CV). ",
    "Structural PK was estimated with NONMEM PRIOR (NWPRI) support from the ",
    "published Sikma 2017 thoracic-transplant tacrolimus popPK model."
  )
  reference <- paste0(
    "Kirubakaran R, Uster DW, Hennig S, Carland JE, Day RO, Wicha SG, Stocker SL. ",
    "Adaptation of a population pharmacokinetic model to inform tacrolimus therapy ",
    "in heart transplant recipients. Br J Clin Pharmacol. 2023;89(4):1162-1175. ",
    "doi:10.1111/bcp.15566. PK structure adapted via the NONMEM PRIOR (NWPRI) ",
    "subroutine from Sikma MA, Hunault CC, Van Maarseveen EM, et al. High ",
    "variability of whole-blood tacrolimus pharmacokinetics early after thoracic ",
    "organ transplantation. Eur J Drug Metab Pharmacokinet. 2020;45(1):123-134. ",
    "doi:10.1007/s13318-019-00591-7."
  )
  vignette <- "Kirubakaran_2022_tacrolimus"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    FFM = list(
      description        = "Fat-free mass",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Power-form allometric scaling on CL/F and Q/F (exponent 0.75) and on V2/F ",
        "and V3/F (exponent 1.0); reference 57 kg (Kirubakaran 2022 Methods 2.5.1 / ",
        "Section 3.3.3; the reference value 57 kg is the median FFM reported in ",
        "Sikma 2017, retained for parameter coherence with the prior model). ",
        "Time-varying (updated at every observation alongside HCT, serum creatinine, ",
        "albumin, AST, ALP per Methods 2.1). Derived per subject from height, ",
        "weight, and sex via the Janmahasatian 2005 formula."
      ),
      source_name        = "FFM"
    ),
    HCT = list(
      description        = "Haematocrit (volume fraction of red blood cells)",
      units              = "fraction (e.g., 0.26)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Power-form effect on CL/F with exponent -0.84 and reference 0.34 ",
        "(Kirubakaran 2022 Table 3 and Section 3.3.3). Reference 0.34 is the ",
        "median HCT in the model-building dataset (Table 1 reports median 0.26 ",
        "in the cohort, but the model-internal centring constant is 0.34). ",
        "Time-varying (carried forward up to one month if missing). Source ",
        "reports HCT as a fraction; do not pass percent values."
      ),
      source_name        = "HCT"
    ),
    CONMED_AZOLE = list(
      description        = "Concomitant azole antifungal therapy indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant azole antifungal)",
      notes              = paste0(
        "Time-varying. Selects the typical CL/F (21.1 L/h without azole, 4.2 L/h ",
        "with azole) and the BSV magnitude on CL/F (61% CV without azole, 89.5% ",
        "CV with azole) per Kirubakaran 2022 Table 3. Pooled across itraconazole ",
        "(n = 715 of 1099 model-building concentrations), voriconazole (n = 23), ",
        "and fluconazole (n = 21) (Section 3.2). All heart transplant recipients ",
        "received itraconazole 200 mg q12h immediately post-transplant as ",
        "Aspergillus prophylaxis, continued for up to 6 months unless ",
        "contraindicated. The Kirubakaran 2022 dataset carries CONMED_AZOLE = 1 ",
        "for one week after azole discontinuation to allow tacrolimus apparent ",
        "clearance to stabilize given itraconazole's long half-life ",
        "(Section 3.3.3); preserve this 1-week post-cessation lag when building ",
        "input data. The 80% CL/F reduction with azole reflects mechanism-based ",
        "CYP3A4 / P-glycoprotein inhibition."
      ),
      source_name        = "AZOLE"
    )
  )

  population <- list(
    n_subjects     = 87L,
    n_studies      = 1L,
    age_range      = "16-70 years (model building 16-70; external evaluation 16-69)",
    age_median     = "53 years (model building); 56 years (external evaluation)",
    weight_range   = "40-111 kg (model building 40-107; external evaluation 45-111)",
    weight_median  = "77 kg (model building); 75 kg (external evaluation)",
    height_range   = "154-195 cm (model building 154-190; external evaluation 150-195)",
    height_median  = "175 cm (both subsets)",
    sex_female_pct = 32.2,
    race_ethnicity = c(White_Caucasian = 48.3, Asian = 11.5, Unknown = 40.2),
    disease_state  = paste0(
      "Heart transplant recipients followed from transplantation to ",
      "approximately 1 year post-transplant; immunosuppressive maintenance with ",
      "oral immediate-release tacrolimus (Prograf), mycophenolate mofetil, and ",
      "tapered prednisolone; basiliximab IV induction in 71/87 (82%); all ",
      "received itraconazole 200 mg q12h prophylaxis until at least the cessation ",
      "of antifungal therapy 6 months post-transplant."
    ),
    dose_range     = paste0(
      "Oral immediate-release tacrolimus q12h, individualized to trough target. ",
      "Median (range) 0.50 mg q12h (0.05-8.00) under concomitant azole antifungal ",
      "and 3.00 mg q12h (0.25-12.00) without concomitant azole (Section 3.2)."
    ),
    regions        = "Single centre, St Vincent's Hospital Sydney, Australia",
    n_concentrations_modelbuild  = 1099L,
    n_concentrations_external    = 348L,
    haematocrit_baseline = "median 0.26 (range 0.21-0.38) (Table 1)",
    albumin_baseline     = "median 34 g/L (range 20-43) (Table 1)",
    creatinine_baseline  = "median 131 umol/L (range 39-276) (Table 1)",
    creatinine_clearance_baseline = "median 63 mL/min (range 25-172) (Cockcroft-Gault, Table 1)",
    diabetes_pct         = 33.3,
    cyp3a5_genotype_known_pct = 51.0,
    notes          = paste0(
      "Retrospective routine-care monitoring data; almost entirely pre-dose ",
      "(trough) concentrations (Section 3.2). Bioanalytical assay: whole-blood ",
      "tacrolimus by UHPLC-MS/MS, quantification range 2-50 ug/L. Model-building ",
      "subset (n = 47, transplanted in 2018) and external-evaluation subset ",
      "(n = 40, transplanted in 2017) reported in Table 1. Cohort enrolled ",
      "1 January 2017 to 31 December 2018."
    )
  )

  ini({
    # ----- Structural PK (Kirubakaran 2022 Table 3 final-model column) -----
    # Two-compartment first-order oral absorption; F = 1 fixed (Section 2.5.1).
    # Reference covariates: FFM = 57 kg (Sikma 2017 median, retained), HCT = 0.34.
    # The model exposes two distinct typical CL/F intercepts -- one without azole
    # and one with concomitant azole antifungal therapy -- because Kirubakaran
    # 2022 estimates separate typical values for these two states (Table 3) and
    # also estimates a state-dependent BSV magnitude (61% CV without azole vs
    # 89.5% CV with azole). Each lcl_<state> intercept therefore has its own
    # eta<lcl_<state>> IIV term.
    lka       <- log(0.508);  label("Absorption rate (Ka, 1/h)")  # Kirubakaran 2022 Table 3 final-model Ka = 0.508 h^-1 (RSE 20%)
    lcl       <- log(21.1);   label("Apparent CL/F at reference covariates without concomitant azole antifungal (CL/F, L/h)")  # Kirubakaran 2022 Table 3 final-model CL/F (without azole) = 21.1 L/h (RSE 11%)
    lcl_azole <- log(4.2);    label("Apparent CL/F at reference covariates with concomitant azole antifungal (CL/F, L/h)")  # Kirubakaran 2022 Table 3 final-model CL/F (with azole) = 4.2 L/h (RSE 12%)
    lvc       <- log(197);    label("Apparent central volume V2/F at reference FFM (V2/F, L)")  # Kirubakaran 2022 Table 3 final-model V2/F = 197 L (RSE 9%)
    lq        <- log(55.0);   label("Apparent intercompartmental clearance Q/F at reference FFM (Q/F, L/h)")  # Kirubakaran 2022 Table 3 final-model Q/F = 55.0 L/h (RSE 10%)
    lvp       <- log(297);    label("Apparent peripheral volume V3/F at reference FFM (V3/F, L)")  # Kirubakaran 2022 Table 3 final-model V3/F = 297 L (RSE 9%)

    # FFM-allometric exponents (fixed at 0.75 / 1.0; Section 3.3.3 final-model
    # equations cite the reported allometry scale of Anderson and Holford 2009 /
    # carried from Sikma 2017). Shared across CL/F + Q/F and V2/F + V3/F
    # respectively per the e_<cov>_<param1>_<param2> convention.
    e_ffm_cl_q  <- 0.75;  label("Power exponent of FFM on CL/F and Q/F (unitless)")  # Kirubakaran 2022 Section 2.5.1 / Section 3.3.3 (fixed allometry)
    e_ffm_vc_vp <- 1.00;  label("Power exponent of FFM on V2/F and V3/F (unitless)")  # Kirubakaran 2022 Section 2.5.1 / Section 3.3.3 (fixed allometry)

    # Haematocrit power effect on CL/F (estimated; ref HCT = 0.34).
    e_hct_cl <- -0.84;  label("Power exponent of haematocrit on CL/F (unitless)")  # Kirubakaran 2022 Table 3 final-model HCT effect on CL/F = -0.84 (RSE 15%)

    # ----- IIV (Kirubakaran 2022 Table 3 final-model column) -----
    # State-dependent BSV magnitude on log(CL/F): 61.0% CV without azole, 89.5%
    # CV with azole. omega^2 = log(CV^2 + 1). Two diagonal IIV terms; one applies
    # per CONMED_AZOLE state at any observation (NONMEM "IF (AZOLE.EQ.1) THEN
    # CL = TVCL * EXP(ETA(2)) ELSE CL = TVCL * EXP(ETA(1))" pattern). Diagonal
    # because the source paper does not report a covariance between the two BSV
    # terms; all between-occasion variabilities (BOV) are fixed to zero in the
    # final model (Methods 2.5.1, Table 3 column 3).
    etalcl       ~ 0.31633   # log(0.610^2 + 1) = 0.31633; 61.0% CV without azole; Kirubakaran 2022 Table 3 BSV CL/F (without azole) (RSE 16%; shrinkage 27%)
    etalcl_azole ~ 0.58836   # log(0.895^2 + 1) = 0.58836; 89.5% CV with azole; Kirubakaran 2022 Table 3 BSV CL/F (with azole) (RSE 11%; shrinkage 3%)

    # ----- Residual error -----
    # Log-transform-both-sides proportional residual error (Section 2.5.1 /
    # Section 3.3.3). NONMEM "additive on log scale" maps to nlmixr2 prop()
    # (linear-space proportional). The reported 41% CV maps to propSd = 0.41.
    propSd <- 0.41;  label("Proportional residual error (fraction)")  # Kirubakaran 2022 Table 3 final-model proportional RUV = 41% (RSE 2%)
  })

  model({
    # ----- 1. Derived covariate terms -----
    # FFM-allometric scaling: reference 57 kg (Kirubakaran 2022 Methods 2.5.1 /
    # Sikma 2017 cohort median). HCT power effect: reference 0.34.
    ffm_cl_q  <- (FFM / 57)^e_ffm_cl_q
    ffm_vc_vp <- (FFM / 57)^e_ffm_vc_vp
    hct_cl    <- (HCT / 0.34)^e_hct_cl

    # ----- 2. Individual PK parameters -----
    # State-dependent CL/F: typical value 21.1 L/h without azole and 4.2 L/h
    # with azole; the IIV term is also state-dependent (etalcl applies when
    # CONMED_AZOLE = 0, etalcl_azole applies when CONMED_AZOLE = 1). Each
    # subject carries one realization of each eta; the AZOLE indicator selects
    # which one drives CL/F at any given observation, reproducing the
    # NONMEM "IF (AZOLE.EQ.1) THEN ... ELSE ..." pattern of the source paper.
    ka          <- exp(lka)
    cl_no_azole <- exp(lcl + etalcl)             * ffm_cl_q * hct_cl
    cl_azole    <- exp(lcl_azole + etalcl_azole) * ffm_cl_q * hct_cl
    cl          <- cl_no_azole * (1 - CONMED_AZOLE) + cl_azole * CONMED_AZOLE
    vc          <- exp(lvc) * ffm_vc_vp
    q           <- exp(lq)  * ffm_cl_q
    vp          <- exp(lvp) * ffm_vc_vp

    # ----- 3. Micro-constants -----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ----- 4. ODE system (2-cmt with first-order oral absorption) -----
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ----- 5. Bioavailability -----
    # F = 1 fixed (Methods 2.5.1; "the population's average oral bioavailability
    # was fixed to one"). No f(depot) statement needed; rxode2 default is 1.

    # ----- 6. Observation and error -----
    # Dose in mg, central amount in mg, vc in L -> mg/L; multiply by 1000 to
    # report tacrolimus whole-blood concentration in ug/L (= ng/mL), which is
    # the bioanalytical assay output and the published concentration unit
    # (Methods 2.3, Figure 2).
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
