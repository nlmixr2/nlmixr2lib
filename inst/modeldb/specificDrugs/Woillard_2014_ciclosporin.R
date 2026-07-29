Woillard_2014_ciclosporin <- function() {
  description <- "Two-compartment population PK model with Erlang-distributed transit absorption (5 sequential delay compartments) and first-order elimination for oral ciclosporin (CsA) in adult haematopoietic stem cell transplant (HSCT) recipients on graft-versus-host disease prophylaxis (Woillard 2014, NONMEM final model). The apparent peripheral volume of distribution Vp/F is fixed at 500 L; no covariate effects were retained in the final model. Combined additive plus proportional residual error."
  reference <- paste(
    "Woillard JB, Lebreton V, Neely M, Turlure P, Girault S, Debord J,",
    "Marquet P, Saint-Marcoux F. Pharmacokinetic tools for the dose",
    "adjustment of ciclosporin in haematopoietic stem cell transplant",
    "patients. Br J Clin Pharmacol 2014; 78(4):836-846.",
    "doi:10.1111/bcp.12394.",
    sep = " "
  )
  vignette <- "Woillard_2014_ciclosporin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariatesDataExcluded <- list(
    HGB = list(
      description = "Whole-blood haemoglobin concentration.",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Tested but not retained. Median 10.0 g/dL (range 8.2-14.8) in the development cohort (Table 1). Spearman screening followed by NONMEM OFV test (P < 0.001, 1 d.f., i.e. decrease of 10.83 in OFV) did not retain a covariate effect."
    ),
    HCT = list(
      description = "Haematocrit, expressed as a percentage.",
      units       = "%",
      type        = "continuous",
      notes       = "Tested and reached the NONMEM screening threshold on apparent clearance (OFV decreased by 17 points) but was not retained in the final model because its addition worsened the precision of inter-patient variability estimates (Results, 'Covariate analysis'). Median 29 % (range 23-43) in the development cohort (Table 1). The authors note clinical plausibility through CsA's red-blood-cell binding."
    ),
    TBILI = list(
      description = "Total serum bilirubin.",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Tested but not retained. Median 11 umol/L (range 4-74) in the development cohort (Table 1)."
    ),
    CREAT = list(
      description = "Serum creatinine.",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Tested but not retained. Median 77 umol/L (range 34-198) in the development cohort (Table 1)."
    ),
    ALB = list(
      description = "Serum albumin.",
      units       = "g/L",
      type        = "continuous",
      notes       = "Tested but not retained. Median 31.4 g/L (range 18.8-55.1) in the development cohort (Table 1; n = 62 with albumin available out of 72 profiles)."
    ),
    ALT = list(
      description = "Alanine aminotransferase (ALAT in the source paper).",
      units       = "U/L",
      type        = "continuous",
      notes       = "Tested but not retained. Median 35 U/L (range 8-195) in the development cohort (Table 1; source paper reports IU/L, value-identical to the SI canonical U/L). The source paper uses the European 'ALAT' abbreviation; the canonical nlmixr2lib name is ALT and the canonical unit string is U/L per the 2026-06-19 SI register."
    ),
    AST = list(
      description = "Aspartate aminotransferase (ASAT in the source paper).",
      units       = "U/L",
      type        = "continuous",
      notes       = "Tested but not retained. Median 27 U/L (range 10-146) in the development cohort (Table 1; source paper reports IU/L, value-identical to the SI canonical U/L). The source paper uses the European 'ASAT' abbreviation; the canonical nlmixr2lib name is AST and the canonical unit string is U/L per the 2026-06-19 SI register."
    ),
    WT = list(
      description = "Body weight at the sampling occasion.",
      units       = "kg",
      type        = "continuous",
      notes       = "Tested but not retained. Median 71 kg (range 47-101) in the development cohort (Table 1; n = 71 of 72 profiles with weight recorded). The authors did not investigate allometric scaling on the structural parameters."
    )
  )

  population <- list(
    species               = "human",
    n_subjects            = 45L,
    n_studies             = 1L,
    n_profiles            = 87L,
    n_dev_subjects        = 40L,
    n_dev_profiles        = 72L,
    n_validation_profiles = 15L,
    age_range             = "24-70 years (overall range across development and validation)",
    age_median            = "Development cohort 59 years (24-67); validation cohort 63 years (48-70)",
    weight_range          = "47-101 kg",
    weight_median         = "Development cohort 71 kg (47-101); validation cohort 72 kg (66-100)",
    sex_female_pct        = 37.5,
    race_ethnicity        = "Not reported in source paper (single-centre French cohort).",
    disease_state         = "Adult HSCT (haematopoietic stem cell transplant) recipients on ciclosporin for graft-versus-host disease (GVHD) prophylaxis. All patients received a reduced-intensity conditioning regimen prior to HSCT (most commonly fludarabine + busulfan + anti-T lymphocyte globulin). 96 % received peripheral blood stem cells, 4 % cord blood. 89 % received T-cell depletion with anti-T lymphocyte globulin. Co-immunosuppression: CsA alone for HLA-sibling donors (n = 11) or CsA + mycophenolate mofetil for matched unrelated donors (n = 29) or HLA-sibling donors without ATG (n = 5). Primary disease was a mix of myeloid acute leukemia, chronic lymphoid leukemia, non-Hodgkin's and Hodgkin's lymphoma, myelodysplastic syndrome, acute lymphoblastic leukemia, and myelofibrosis.",
    dose_range            = "Ciclosporin oral 3 mg/kg twice daily, starting 3 days before engraftment (day -3) and titrated by blood trough concentrations monitored twice weekly. Target inter-dose AUC(0,12h) = 4.3 mg/L*h, taken from prior renal-transplant literature. If no GVHD was noted, CsA was tapered progressively over 4 weeks starting day 100.",
    regions               = "France (single centre: Limoges University Hospital).",
    sampling_schedule     = "Per profile: 10 samples at pre-dose and 0.33, 0.66, 1, 2, 3, 4, 6, 8 and 12 hours post-dose. One patient had 4 sampling periods, 8 patients had 3, 21 had 2, and 17 had 1. Sampling periods spanned day 0 to day 100 post-transplant.",
    bioanalytical         = "Whole-blood ciclosporin quantified by turbulent-flow LC-MS/MS (Cyclone P online extraction, Propel MS C18 analytical column, TSQ Quantum Discovery MS/MS) with calibration over 10-2000 ug/L, LOD = 10 ug/L and LOQ = 20 ug/L. Inter-assay precision (RSD) -3.1 to 11.8 %, mean relative error 4.0 to 11.7 %.",
    baseline_demographics = "Median (range) in the development cohort (Table 1): age 59 (24-67) years, weight 71 (47-101) kg, sex M/F 25/15, haematocrit 29 (23-43) %, haemoglobin 10.0 (8.2-14.8) g/dL, serum creatinine 77 (34-198) umol/L, total bilirubin 11 (4-74) umol/L, albumin 31.4 (18.8-55.1) g/L, ALAT 35 (8-195) U/L, ASAT 27 (10-146) U/L. Sampling time post-transplant median 4 days (range 0-99).",
    notes                 = "Three independent modelling approaches (NONMEM, iterative two-stage ITS, non-parametric Pmetrics) were fit in parallel to compare Bayesian estimators of CsA AUC(0,12h) under a three-sample limited sampling strategy. The packaged model file encodes the NONMEM final model (Table 2) parameterised in standard CL/V/Q form. The ITS and Pmetrics fits used a gamma-law absorption with macro-constant disposition (FAIV, FBIV, alpha, beta in Table 3) and are not packaged here; see the vignette's Assumptions and deviations section for rationale."
  )

  ini({
    # NONMEM final-model estimates from Woillard 2014 Table 2. Structural
    # model: two-compartment with first-order elimination and Erlang
    # absorption with five sequential delay compartments (Results,
    # "Pharmacokinetic models"). FOCE INTER method via Wings for NONMEM
    # (Methods, "Non-linear mixed effects modelling"). No covariates were
    # retained in the final model (Results, "Covariate analysis").

    # Erlang transit absorption rate constant.
    lktr <- log(5.72);  label("Erlang transit absorption rate constant Ktr (1/h)")           # Table 2, Ktr = 5.72 1/h

    # Apparent inter-compartmental clearance.
    lq   <- log(33.7);  label("Apparent inter-compartmental clearance Q/F (L/h)")            # Table 2, Q/F = 33.7 L/h

    # Apparent central volume of distribution.
    lvc  <- log(222);   label("Apparent central volume of distribution Vc/F (L)")            # Table 2, Vc/F = 222 L

    # Apparent peripheral volume of distribution. The paper notes the
    # peripheral volume was "arbitrarily fixed to 500 L" (Results,
    # "Pharmacokinetic models"), so the value is encoded with fixed().
    lvp  <- fixed(log(500)); label("Apparent peripheral volume of distribution Vp/F (L), fixed") # Table 2, Vp/F = 500 L (fixed; NA SE/CI)

    # Apparent oral clearance.
    lcl  <- log(41.2);  label("Apparent oral clearance CL/F (L/h)")                          # Table 2, CL/F = 41.2 L/h

    # Inter-patient variability (Table 2 IPV column). Conversion from
    # reported %CV to log-scale variance via omega^2 = log(1 + CV^2):
    #   Ktr   CV 50 % -> log(1 + 0.50^2) = 0.22314
    #   Q/F   CV 72 % -> log(1 + 0.72^2) = 0.41759
    #   Vc/F  CV 57 % -> log(1 + 0.57^2) = 0.28140
    #   CL/F  CV 40 % -> log(1 + 0.40^2) = 0.14842
    # No IIV is reported on Vp/F (Vp/F was fixed; IPV is "NA" in Table 2).
    etalktr ~ 0.22314                                                                         # Table 2, IPV Ktr = 50 % CV
    etalq   ~ 0.41759                                                                         # Table 2, IPV Q/F = 72 % CV
    etalvc  ~ 0.28140                                                                         # Table 2, IPV Vc/F = 57 % CV
    etalcl  ~ 0.14842                                                                         # Table 2, IPV CL/F = 40 % CV

    # Combined additive + proportional residual error. The paper's Table 2
    # footer reports "Proportional error = 12.53 % and additive error =
    # 0.023 mg ml-1". The "mg ml-1" units in the footer are inconsistent
    # with the rest of the paper, which expresses ciclosporin
    # concentrations in mg/L (e.g. assay calibration 10-2000 ug/L = 0.010-2
    # mg/L; target AUC(0,12h) = 4.3 mg/L*h). 0.023 mg/mL = 23 mg/L is well
    # above the highest calibrator, whereas 0.023 mg/L = 23 ug/L sits near
    # the LOQ (20 ug/L) and is the physically reasonable value. The
    # additive error is therefore interpreted as 0.023 mg/L; the apparent
    # unit-string typo is documented in the vignette's Assumptions and
    # deviations section.
    propSd <- 0.1253; label("Proportional residual error (fraction)")                         # Table 2 footer, proportional error = 12.53 %
    addSd  <- 0.023;  label("Additive residual error (mg/L)")                                 # Table 2 footer, additive error = 0.023 (paper units "mg ml-1" read as mg/L)
  })

  model({
    # Individual PK parameters (exponential IIV; no covariate effects in
    # the final NONMEM model per Results, "Covariate analysis").
    ktr <- exp(lktr + etalktr)
    q   <- exp(lq   + etalq)
    vc  <- exp(lvc  + etalvc)
    vp  <- exp(lvp)
    cl  <- exp(lcl  + etalcl)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment disposition with Erlang absorption (n = 5 sequential
    # delay compartments). Dose enters depot; depot + transit1 + transit2
    # + transit3 + transit4 form the five sequential delay compartments
    # connected by the common rate constant ktr; transit4 empties into
    # central at rate ktr. Bioavailability F is folded into the apparent
    # parameters CL/F, Q/F, Vc/F, Vp/F.
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot     - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1  - ktr * transit2
    d/dt(transit3)    <-  ktr * transit2  - ktr * transit3
    d/dt(transit4)    <-  ktr * transit3  - ktr * transit4
    d/dt(central)     <-  ktr * transit4  - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central   - k21 * peripheral1

    # Concentration in central compartment. Dose units mg, volume units L
    # give central / vc in mg/L, matching the source paper's reporting
    # convention (target AUC(0,12h) = 4.3 mg/L*h; assay calibration
    # 0.010-2 mg/L).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
