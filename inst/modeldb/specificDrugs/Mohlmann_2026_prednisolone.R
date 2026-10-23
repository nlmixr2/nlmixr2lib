Mohlmann_2026_prednisolone <- function() {
  description <- paste(
    "Two-compartment population PK model for total and protein-unbound prednisolone in",
    "children with immune-mediated and systemic inflammatory diseases, oral or IV",
    "(Mohlmann 2026). First-order absorption with F fixed to 1; the disposition is",
    "parameterised on UNBOUND prednisolone with allometric scaling to 70 kg (exponents",
    "fixed at 0.75 and 1). Total prednisolone is recovered algebraically from linear",
    "binding to albumin plus saturable binding to corticosteroid-binding globulin",
    "(transcortin), the latter supplied as a time-varying covariate that follows its",
    "circadian rhythm. Clearance differs by cohort (HCT prophylaxis, HCT GvHD therapy,",
    "autoimmune disease as reference)."
  )
  reference <- paste(
    "Mohlmann JE, Lindemans CA, Jansen MHA, van Luin M, Punt AM, Ezzafzafi S,",
    "Nierkens S, Huitema ADR, de Vries Schultink AHM. Population Pharmacokinetics of",
    "Total and Protein-Unbound Prednisolone in Children with Immune-Mediated and",
    "Systemic Inflammatory Diseases. Clin Pharmacokinet. 2026;65(3):277-292.",
    "doi:10.1007/s40262-025-01601-5"
  )
  vignette <- "Mohlmann_2026_prednisolone"
  units <- list(time = "h", dosing = "mg", concentration = "umol/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling of CL and Q (exponent fixed at 0.75) and of Vc and Vp (exponent",
        "fixed at 1), reference 70 kg (Methods 2.4.1; Online Resource 3 $PK '(WGT/70)**0.75',",
        "'(WGT/70)**1'). The exponents were estimated at 0.67 and 0.84 but the fixed",
        "theoretical values were kept (Results 3.3). Cohort median 36.5 kg, range 5.0-109",
        "(Table 1)."
      ),
      source_name = "WGT"
    ),
    ALB = list(
      description = "Serum albumin concentration",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Structural input to the linear albumin-binding arm, with no reference value and no",
        "covariate-effect coefficient. Converted to umol/L inside model() as ALB / 0.0665",
        "(albumin 66.5 kDa), exactly as Online Resource 3 'ALBX = ALBM/0.0665'. The value",
        "measured on (or closest to) the PK day was used (Methods 2.4.1); missing values",
        "were imputed with the cohort median (Methods 2.3). Cohort median 31.1 g/L, range",
        "13.9-44.3 (Table 1)."
      ),
      source_name = "ALBM"
    ),
    TRANSCORTIN = list(
      description = "Plasma corticosteroid-binding globulin (transcortin) concentration, time-varying",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Structural binding capacity of the saturable CBG arm, with no reference value and no",
        "covariate-effect coefficient. Converted to umol/L inside model() as TRANSCORTIN / 52",
        "(CBG 52 kDa, Discussion 4.3). Time-varying: the authors measured CBG once per",
        "occasion and used the Melin et al. circadian model (Eqs 7-8: two cosines, AMP24",
        "0.055, SHIFT24 1.8 h, AMP12 0.029, SHIFT12 15.7 h, TIME in hours after a 3 PM",
        "baseline) to obtain per-individual estimates of the full 24-h curve, which were",
        "supplied as the data column CBGpm (umol/L; CBGp in mg/L). Supply the curve at every",
        "record. Baseline median 1.04 umol/L (about 54 mg/L), range 0.60-2.48 umol/L",
        "(Results 3.2)."
      ),
      source_name = "CBGp"
    ),
    TX_HCT = list(
      description = "Haematopoietic cell transplant recipient (PIKACHU cohort)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "1 = prednisolone given after allogeneic HCT, as GvHD prophylaxis or as GvHD /",
        "engraftment-syndrome therapy; 0 = autoimmune or systemic inflammatory disease",
        "(CARPE DIEM cohort, the reference). Together with AGVHD_ANY it decomposes the",
        "control-stream column PorT (0 = HCT prophylaxis, 1 = HCT therapy, 2 = autoimmune):",
        "PorT 0 is TX_HCT = 1, AGVHD_ANY = 0; PorT 1 is TX_HCT = 1, AGVHD_ANY = 1; PorT 2 is",
        "TX_HCT = 0, AGVHD_ANY = 0."
      ),
      source_name = "PorT"
    ),
    AGVHD_ANY = list(
      description = "Acute graft-versus-host disease (any organ) being treated with prednisolone",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "1 = HCT recipient receiving the intensified (2 mg/kg/day) prednisolone regimen as",
        "therapy of GvHD (or of engraftment syndrome, which the paper pools into the same",
        "HCT-therapy level); 0 = otherwise. Only meaningful with TX_HCT = 1. Time-varying per",
        "occasion: eight patients were sampled first on prophylaxis and later on GvHD",
        "therapy and were kept as one subject with two occasions (Methods 2.2.1)."
      ),
      source_name = "PorT"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Tested on CL and V, not retained (Methods 2.4.3, Results 3.4)."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested on CL as continuous and as categories split at 0.5, 1, 2 and 12 years; no effect beyond allometric scaling (Methods 2.4.3, Results 3.4)."
    ),
    PAGE = list(
      description = "Postmenstrual age (age in weeks plus 40 weeks) for a Hill maturation function on CL",
      units = "weeks",
      type = "continuous",
      notes = "Maturation (Eq. 6, Hill 3.92 and TM50 54.2 weeks fixed, or estimated) tested on CL in addition to and instead of allometric scaling; not retained (Methods 2.4.3, Results 3.4)."
    ),
    CONMED_CYP3A4_INH = list(
      description = "Concomitant CYP3A4 inhibitor",
      units = "(binary)",
      type = "binary",
      notes = "Tested on CL; a trend toward lower CL but not significant (Results 3.4). 55.9% of occasions (Table 1)."
    ),
    CONMED_CSA = list(
      description = "Concomitant cyclosporin A",
      units = "(binary)",
      type = "binary",
      notes = "Tested on F and CL, not retained (Methods 2.4.3, Results 3.4). 69.1% of occasions (Table 1)."
    ),
    HCT_TCD_ATG = list(
      description = "Antithymocyte globulin serotherapy",
      units = "(binary)",
      type = "binary",
      notes = "Tested on CL and V, not retained (Methods 2.4.3; control-stream column SERO)."
    ),
    TBILI = list(
      description = "Total bilirubin (liver function)",
      units = "umol/L",
      type = "continuous",
      notes = "Liver-function markers (bilirubin, AST, ALT, GGT, alkaline phosphatase) tested on CL and V, not retained (Methods 2.4.3, Results 3.4)."
    ),
    ROUTE_ORAL = list(
      description = "Formulation and oral route (tablet, suspension, feeding tube)",
      units = "(binary)",
      type = "binary",
      notes = "Formulation and oral route tested on ka, and GvHD diagnosis tested on ka; not retained (Methods 2.4.3)."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "prednisolone", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "prednisolone unbound", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "prednisolone unbound", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 60,
    n_studies = 2,
    n_observations = "305 total and 304 unbound prednisolone serum concentrations from 68 PK occasions (1-5 samples per occasion)",
    age_range = "0.2-19 years",
    age_median = "10 years",
    weight_range = "5.0-109 kg",
    weight_median = "36.5 kg",
    sex_female_pct = 43.3,
    race_ethnicity = "not reported",
    disease_state = paste(
      "44 children after allogeneic haematopoietic cell transplantation (prednisolone as",
      "GvHD prophylaxis or as GvHD / engraftment-syndrome therapy; ALL, AML, Hurler",
      "syndrome and other indications) and 16 children with autoimmune or systemic",
      "inflammatory disease (SLE, uveitis, IBD and others)"
    ),
    dose_range = "Prednisolone 0.25-2 mg/kg/day (4-120 mg/day), PK dose 0.5 mg/kg (0.25-1); 66.2% IV, 33.8% oral",
    regions = "The Netherlands (University Medical Centre Utrecht / Princess Maxima Centre)",
    co_medication = "CYP3A4 inhibitors 55.9%, cyclosporin A 69.1%, no CYP3A4 inducers (Table 1)",
    notes = paste(
      "Prospective PIKACHU (HCT, from October 2020) and CARPE DIEM (autoimmune, from April",
      "2022) studies, data to March 2024. Total and unbound (equilibrium dialysis)",
      "prednisolone by LC-MS/MS; CBG by radioimmunoassay. Median albumin 31.1 g/L, median",
      "measured CBG 1.09 umol/L (Table 1). NONMEM 7.5.1, FOCE-I; precision by SIR."
    )
  )

  ini({
    # Final control stream (Online Resource 3) $THETA values; Table 2 prints
    # the same estimates rounded.
    lka <- log(0.68); label("Absorption rate constant (1/h)") # Online Resource 3 THETA(1) = 0.68; Table 2 ka 0.68
    lcl <- log(155); label("Clearance of unbound prednisolone, autoimmune cohort, 70 kg (L/h)") # Online Resource 3 THETA(2) = 155; Table 2 CL70kg 155 L/h
    lvc <- log(338); label("Central volume of unbound prednisolone, 70 kg (L)") # Online Resource 3 THETA(3) = 338; Table 2 Vc,70kg 338 L
    lq <- log(105); label("Intercompartmental clearance, 70 kg (L/h)") # Online Resource 3 THETA(4) = 105; Table 2 Q70kg 105 L/h
    lvp <- log(170); label("Peripheral volume, 70 kg (L)") # Online Resource 3 THETA(5) = 170; Table 2 Vp,70kg 170 L

    # Allometric exponents fixed a priori (Methods 2.4.1; kept fixed in
    # Results 3.3).
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL and Q (unitless)") # Online Resource 3 '(WGT/70)**0.75'
    e_wt_vc <- fixed(1); label("Allometric exponent on Vc and Vp (unitless)") # Online Resource 3 '(WGT/70)**1'

    # Equilibrium dissociation constants of the binding model (Eq. 2), umol/L.
    lkalb <- log(200); label("Dissociation constant of prednisolone for albumin, Kalb (umol/L)") # Online Resource 3 THETA(6) = 200; Table 2 Kalb 200 uM
    lkcbg <- log(0.0476); label("Dissociation constant of prednisolone for corticosteroid-binding globulin, KCBG (umol/L)") # Online Resource 3 THETA(7) = 0.0476; Table 2 KCBG 0.048 uM (rounded); Results 3.4 '0.0476'

    # Cohort effects on CL, multiplicative (THETA^indicator); reference is
    # the autoimmune cohort.
    e_hct_proph_cl <- 1.1; label("CL multiplier, HCT prophylaxis cohort vs autoimmune (unitless)") # Online Resource 3 THETA(8) = 1.1; Table 2 HCT-P 1.10
    e_hct_ther_cl <- 0.975; label("CL multiplier, HCT GvHD-therapy cohort vs autoimmune (unitless)") # Online Resource 3 THETA(9) = 0.975; Table 2 HCT-T 0.97 (rounded)

    # IIV: Online Resource 3 $OMEGA BLOCK(3) (CL, Vc, Vp) with the CL-Vp
    # covariance fixed to 0, plus a diagonal $OMEGA on ka. Variances;
    # Table 2 prints sqrt(omega) as %.
    etalcl + etalvc + etalvp ~ c(
      0.0993,
      0.097, 0.365,
      0, -0.14, 0.153
    ) # Online Resource 3 $OMEGA BLOCK(3); Table 2 IIV CL 31.5%, Vc 60.4%, Vp 39.1%, corr CL-Vc 51.0%, corr Vc-Vp -59.0%
    etalka ~ 0.268 # Online Resource 3 $OMEGA 0.268; Table 2 IIV ka 51.8%

    # Residual error: proportional on each output. $SIGMA holds variances.
    propSd <- sqrt(0.114); label("Proportional residual error, unbound prednisolone (fraction)") # Online Resource 3 $SIGMA EPS(1) 0.114 (CPU); Table 2 33.8%
    propSd_Ctot <- sqrt(0.0703); label("Proportional residual error, total prednisolone (fraction)") # Online Resource 3 $SIGMA EPS(2) 0.0703 (CP); Table 2 26.5%
  })

  model({
    # Cohort indicators (control-stream PP and PT): autoimmune = reference.
    hct_proph <- TX_HCT * (1 - AGVHD_ANY)
    hct_ther <- TX_HCT * AGVHD_ANY

    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl *
      e_hct_proph_cl^hct_proph * e_hct_ther_cl^hct_ther
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q <- exp(lq) * (WT / 70)^e_wt_cl
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc
    kalb <- exp(lkalb)
    kcbg <- exp(lkcbg)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Both states hold UNBOUND prednisolone: CL and Q act on the unbound
    # pool and the bound drug is an output transformation, not a
    # reservoir. Oral doses enter depot (F fixed to 1, Results 3.3); IV
    # doses (the salt-corrected prednisolone base amount) enter central as
    # a bolus.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Unit conversions from Online Resource 3: mg/L prednisolone to umol/L
    # ('S2 = V2/2.77', MW 360.4); albumin g/L to umol/L ('ALBX = ALBM/0.0665',
    # 66.5 kDa); CBG mg/L to umol/L (52 kDa, Discussion 4.3).
    alb_um <- ALB / 0.0665
    cbg_um <- TRANSCORTIN / 52

    # Unbound concentration (umol/L), CPU in the control stream. By library
    # convention the central drug concentration is Cc even when the pool is
    # unbound.
    Cc <- central / vc * 2.77
    # Total concentration (Eqs 1-2; control stream CP): unbound plus
    # saturable CBG binding plus linear albumin binding.
    Ctot <- Cc + Cc * cbg_um / (Cc + kcbg) + Cc * alb_um / kalb

    Cc ~ prop(propSd)
    Ctot ~ prop(propSd_Ctot)
  })
}
