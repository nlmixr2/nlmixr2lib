Wojciechowski_2025_fidanacogeneElaparvovec <- function() {
  description <- paste(
    "Population model of factor IX (FIX) activity after the AAV gene therapy",
    "fidanacogene elaparvovec and/or FIX replacement therapy in hemophilia B",
    "(Wojciechowski 2025). A Hargrove-Schmidt gene/protein-expression model",
    "(productive transgene amount translated into FIX at a synthesis site,",
    "first-order transgene degradation, and for a transient-production",
    "mixture subpopulation a reversible transition to a non-productive",
    "transgene state) feeds a three-compartment FIX disposition model with",
    "an endogenous baseline. Standard half-life FIX products dose into the",
    "same central compartment; extended half-life products use the fixed",
    "three-compartment rFIXFc model of Diao 2014. Allometric weight scaling",
    "referenced to 70 kg, age effects on synthesis, transgene degradation",
    "and relative FIX activity, manufacturing-process effects on synthesis,",
    "Box-Cox-transformed IIV on the translation constant, study effects on",
    "baseline, relative activity and the log-scale residual error.",
    sep = " "
  )
  reference <- paste(
    "Wojciechowski J, Gaitonde P, Hughes JH, Ravva P.",
    "Population Modeling of Factor IX Activity Following Administration of",
    "Fidanacogene Elaparvovec Gene Therapy in Participants with Hemophilia B.",
    "Clin Pharmacokinet. 2025;64(10):1531-1548. doi:10.1007/s40262-025-01535-y.",
    "Parameter values from Table 2; model structure from the ESM NONMEM",
    "control stream (run57.mod). Extended half-life FIX disposition fixed",
    "from Diao L, et al. Clin Pharmacokinet. 2014;53(5):467-477",
    "(doi:10.1007/s40262-013-0129-7), as hard-coded in run57.mod.",
    sep = " "
  )
  vignette <- "Wojciechowski_2025_fidanacogeneElaparvovec"
  units <- list(
    time = "h",
    dosing = "vg (fidanacogene elaparvovec, into transgene); IU (FIX replacement, into central or central_ehl)",
    concentration = "IU/dL (FIX activity)"
  )
  dosing <- c("transgene", "central")

  # Gene-therapy expression states and the extended-half-life FIX product
  # states are specific to this model.
  paper_specific_compartments <- c(
    "transgene",
    "transgene_np",
    "fix_cell",
    "central_ehl",
    "peripheral1_ehl",
    "peripheral2_ehl"
  )

  compartmentData <- list(
    transgene = list(
      analyte = "productive fidanacogene elaparvovec transgene (DEPOT in run57.mod)",
      units = "vg",
      specimen = "tissue",
      verified = TRUE
    ),
    transgene_np = list(
      analyte = "non-productive transgene, transient-production subpopulation only (PERIG)",
      units = "vg",
      specimen = "tissue",
      verified = TRUE
    ),
    fix_cell = list(
      analyte = "transgene-derived FIX at the site of synthesis (CELL)",
      units = "IU",
      specimen = "tissue",
      verified = TRUE
    ),
    central = list(
      analyte = "FIX activity from endogenous, transgene-derived and standard half-life FIX (CENT)",
      units = "IU",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "FIX in the peripheral compartment (PERI)",
      units = "IU",
      specimen = "tissue",
      verified = TRUE
    ),
    peripheral2 = list(
      analyte = "FIX in the extravascular compartment (EXTR)",
      units = "IU",
      specimen = "tissue",
      verified = TRUE
    ),
    central_ehl = list(
      analyte = "extended half-life FIX product activity (CENTA)",
      units = "IU",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1_ehl = list(
      analyte = "extended half-life FIX product, peripheral compartment (PERIA)",
      units = "IU",
      specimen = "tissue",
      verified = TRUE
    ),
    peripheral2_ehl = list(
      analyte = "extended half-life FIX product, extravascular compartment (EXTRA)",
      units = "IU",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Fixed allometric scaling referenced to 70 kg: exponent 0.75 on CL, Q2, Q3;",
        "1 on V1, V2, V3; -0.25 on ksyn, kdeg, ktr,d, ktr,p. Estimated power",
        "exponent 1.08 on the translation constant ktau (referenced to 70 kg).",
        "The fixed Diao 2014 extended-half-life sub-model uses its own (WT/73)",
        "powers 0.436 on CL and 0.396 on V1."
      ),
      source_name = "WT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effects referenced to 35 years on ksyn (-0.704), kdeg (-1.57)",
        "and relative FIX activity FIXA (0.0395). The control stream applies",
        "each age effect only when AGE > 0 (factor 1 otherwise)."
      ),
      source_name = "AGE"
    ),
    FORM_FIDANACOGENE_PROCESS1 = list(
      description = "Fidanacogene elaparvovec manufacturing process 1 indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (process 3, the commercial process)",
      notes = "ksyn multiplied by (1 + 0.356). Mutually exclusive with FORM_FIDANACOGENE_PROCESS2. Source CSGRPN = 1.",
      source_name = "CSGRPN"
    ),
    FORM_FIDANACOGENE_PROCESS2 = list(
      description = "Fidanacogene elaparvovec manufacturing process 2 indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (process 3, the commercial process)",
      notes = "ksyn multiplied by (1 + 0.658). Mutually exclusive with FORM_FIDANACOGENE_PROCESS1. Source CSGRPN = 2.",
      source_name = "CSGRPN"
    ),
    MIX_TRANSIENT_EXPR = list(
      description = "Mixture-model class indicator for the transient FIX production phenotype (subpopulation 2)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (subpopulation 1, constant FIX production)",
      notes = paste(
        "Latent NONMEM $MIX class (MIXNUM = 2). When 1, the productive",
        "transgene transitions to a non-productive state (half-life 4.90 weeks)",
        "and back (62.6 weeks); when 0 both transition rates are 0.",
        "Population probability of subpopulation 1 = exp(0.945)/(1 + exp(0.945))",
        "= 0.720 (Table 2 log-odds 0.945), so for population simulation draw",
        "MIX_TRANSIENT_EXPR ~ Bernoulli(0.280) per subject."
      ),
      source_name = "MIXNUM"
    ),
    STUDY_B1821002 = list(
      description = "Nonacog alfa study B1821002 (NCT00364182) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (other studies)",
      notes = "BASE multiplied by (1 + 0.788); residual SD multiplied by (1 + 1.32).",
      source_name = "PROT"
    ),
    STUDY_B1821010 = list(
      description = "Nonacog alfa study B1821010 (NCT01335061) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (other studies)",
      notes = "FIXA multiplied by (1 + 0.252); residual SD multiplied by (1 + 0.545).",
      source_name = "PROT"
    ),
    STUDY_B1821034 = list(
      description = "Nonacog alfa study B1821034 (previously untreated patients) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (other studies)",
      notes = "BASE multiplied by (1 + 0.994); residual SD multiplied by (1 + 0.843).",
      source_name = "PROT"
    ),
    STUDY_B1821036 = list(
      description = "Nonacog alfa study B1821036 (NCT00093171) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (other studies)",
      notes = "FIXA multiplied by (1 + 0.173); residual SD multiplied by (1 - 0.245).",
      source_name = "PROT"
    ),
    STUDY_C037 = list(
      description = "Fidanacogene elaparvovec studies C0371002, C0371003 and C0371005 (pooled) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (other studies)",
      notes = "Residual SD multiplied by (1 + 0.305). Set to 1 to simulate observations from the gene-therapy studies.",
      source_name = "PROT"
    )
  )

  covariatesDataExcluded <- list(
    CONMED_STEROID = list(
      description = "Concomitant corticosteroid use",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Univariately significant on ktr,d (4.97-fold faster transition to the",
        "non-productive state, Table S4 run 12) but judged descriptive and not",
        "carried forward to the multivariate model (Results, Covariate Model)."
      )
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened on ktau (Table S4 run 9, p = 0.0544); not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 337,
    n_studies = 11,
    age_range = "< 0.1-64.7 years (median 21.9); gene-therapy studies 18.0-62.0 years (median 29.0)",
    weight_range = "1.30-172 kg (median 67.2); gene-therapy studies 55.0-139 kg (median 86.1)",
    sex_female_pct = 0,
    race_ethnicity = c(White = 82.8, Asian = 5.6, Black = 5.0, Other = 5.3, Missing = 1.2),
    disease_state = "Hemophilia B (moderately severe to severe in the gene-therapy studies, FIX activity <= 2%)",
    dose_range = paste(
      "Fidanacogene elaparvovec single 60-min IV infusion, mostly 5e11 vg/kg",
      "(n = 63, three studies); nonacog alfa 50-100 IU/kg on demand or",
      "prophylaxis (n = 274, eight studies)"
    ),
    regions = "Multinational",
    notes = paste(
      "8931 FIX activity observations by one-stage clotting assay (Actin FSL",
      "reagent in the gene-therapy studies); 8.4% BLQ observations excluded.",
      "Manufacturing process 1/2/3 in 10/8/45 gene-therapy participants;",
      "31 of 63 received concomitant corticosteroids (Table 1)."
    )
  )

  ini({
    # FIX disposition (nonacog alfa stage; Table 2). Volumes and clearances in L
    # as reported; converted to dL in model() as in run57.mod (10 * THETA).
    lrbase <- log(1.83); label("Baseline FIX activity BASE (IU/dL)") # Table 2 'Baseline (BASE; IU/dL)' = 1.83
    lcl <- log(0.617); label("Clearance CL for a 70 kg individual (L/h)") # Table 2 'Clearance (CL; L/h)' = 0.617
    lvc <- log(10.3); label("Central volume V1 for a 70 kg individual (L)") # Table 2 'Volume of the central compartment (V1; L)' = 10.3
    lq <- log(1.46); label("Central-peripheral intercompartmental clearance Q2 (L/h)") # Table 2 'Q2; L/h' = 1.46
    lvp <- log(3.29); label("Peripheral volume V2 (L)") # Table 2 'Volume of the peripheral compartment (V2; L)' = 3.29
    lq2 <- log(0.145); label("Central-extravascular intercompartmental clearance Q3 (L/h)") # Table 2 'Q3; L/h' = 0.145
    lvp2 <- log(14.0); label("Extravascular volume V3 (L)") # Table 2 'Volume of the extravascular compartment (V3; L)' = 14.0

    # Gene therapy (Table 2; run57.mod converts each half-life to a rate constant)
    lthalf_ksyn <- log(34.9); label("Half-life of FIX synthesis, ln(2)/ksyn (day)") # Table 2 'Half-life of FIX synthesis (ksyn; d)' = 34.9
    lktau <- log(24.0); label("Translation constant ktau (IU/vg/h, x 1e-13)") # Table 2 'Translation constant (ktau; IU/vg x 10^-13 h^-1)' = 24.0
    lthalf_kdeg <- log(423); label("Half-life of transgene degradation, ln(2)/kdeg (week)") # Table 2 'Half-life of transgene degradation (kdeg; weeks)' = 423
    lthalf_ktrd <- log(4.90); label("Half-life of productive-to-non-productive transgene transition (week)") # Table 2 'ktr,d; weeks' = 4.90
    lthalf_ktrp <- log(62.6); label("Half-life of non-productive-to-productive transgene transition (week)") # Table 2 'ktr,p; weeks' = 62.6
    boxcox_lktau <- -0.631; label("Box-Cox shape parameter on the ktau random effect (unitless)") # Table 2 'Box-Cox transformation parameter on ktau (BXCXK)' = -0.631

    # Fixed allometric exponents, referenced to 70 kg (Methods; Table 2 footnote)
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL, Q2, Q3 (unitless)") # Methods: fixed exponent 0.75 for clearances
    e_wt_vc <- fixed(1); label("Allometric exponent on V1, V2, V3 (unitless)") # Methods: fixed exponent 1 for volumes
    e_wt_krate <- fixed(-0.25); label("Allometric exponent on ksyn, kdeg, ktr,d, ktr,p (unitless)") # Methods: fixed exponent -0.25 for first-order rate constants

    # Covariate effects (Table 2)
    e_age_ksyn <- -0.704; label("Power exponent of age/35 on ksyn (unitless)") # Table 2 'Effect of age on ksyn (referenced to 35 years)' = -0.704
    e_age_kdeg <- -1.57; label("Power exponent of age/35 on kdeg (unitless)") # Table 2 'Effect of age on kdeg (referenced to 35 years)' = -1.57
    e_age_fixa <- 0.0395; label("Power exponent of age/35 on relative FIX activity FIXA (unitless)") # Table 2 'Effect of age on FIXA (referenced to 35 years)' = 0.0395
    e_form_p1_ksyn <- 0.356; label("Fractional change in ksyn for manufacturing process 1 (fraction)") # Table 2 'Effect of manufacturing process 1 on ksyn' = 0.356
    e_form_p2_ksyn <- 0.658; label("Fractional change in ksyn for manufacturing process 2 (fraction)") # Table 2 'Effect of manufacturing process 2 on ksyn' = 0.658
    e_wt_ktau <- 1.08; label("Power exponent of WT/70 on ktau (unitless)") # Table 2 'Effect of body weight on ktau (referenced to 70 kg)' = 1.08

    # Study effects (Table 2; reference = all other studies)
    e_study1002_rbase <- 0.788; label("Fractional change in BASE for study B1821002 (fraction)") # Table 2 'Effect of B1821002 study on BASE' = 0.788
    e_study1034_rbase <- 0.994; label("Fractional change in BASE for study B1821034 (fraction)") # Table 2 'Effect of B1821034 study on BASE' = 0.994
    e_study1010_fixa <- 0.252; label("Fractional change in FIXA for study B1821010 (fraction)") # Table 2 'Effect of B1821010 study on FIXA' = 0.252
    e_study1036_fixa <- 0.173; label("Fractional change in FIXA for study B1821036 (fraction)") # Table 2 'Effect of B1821036 study on FIXA' = 0.173
    e_study1002_expSd <- 1.32; label("Fractional change in residual SD for study B1821002 (fraction)") # Table 2 'Effect of B1821002 study on RUV' = 1.32
    e_study1010_expSd <- 0.545; label("Fractional change in residual SD for study B1821010 (fraction)") # Table 2 'Effect of B1821010 study on RUV' = 0.545
    e_study1034_expSd <- 0.843; label("Fractional change in residual SD for study B1821034 (fraction)") # Table 2 'Effect of B1821034 study on RUV' = 0.843
    e_study1036_expSd <- -0.245; label("Fractional change in residual SD for study B1821036 (fraction)") # Table 2 'Effect of B1821036 study on RUV' = -0.245
    e_studyc037_expSd <- 0.305; label("Fractional change in residual SD for the C037 gene-therapy studies (fraction)") # Table 2 'Effect of C037 studies on RUVADD' = 0.305

    # Extended half-life FIX products: Diao 2014 rFIXFc model, hard-coded
    # (no uncertainty, no IIV) in run57.mod and Table S2.
    lcl_ehl <- fixed(log(2.39)); label("EHL product clearance for a 73 kg individual (dL/h)") # run57.mod POPCLA = 2.39*(WT/73)**0.436 (Diao 2014)
    lvc_ehl <- fixed(log(71.4)); label("EHL product central volume for a 73 kg individual (dL)") # run57.mod POPV1A = 71.4*(WT/73)**0.396 (Diao 2014)
    lq_ehl <- fixed(log(1.67)); label("EHL product central-peripheral clearance (dL/h)") # run57.mod POPQ2A = 1.67 (Diao 2014)
    lvp_ehl <- fixed(log(87.0)); label("EHL product peripheral volume (dL)") # run57.mod POPV2A = 87.0 (Diao 2014)
    lq2_ehl <- fixed(log(39.3)); label("EHL product second intercompartmental clearance (dL/h)") # run57.mod POPQ3A = 39.3 (Diao 2014)
    lvp2_ehl <- fixed(log(39.9)); label("EHL product second peripheral volume (dL)") # run57.mod POPV3A = 39.9 (Diao 2014)
    e_wt_cl_ehl <- fixed(0.436); label("Power exponent of WT/73 on EHL product clearance (unitless)") # run57.mod POPCLA exponent 0.436 (Diao 2014)
    e_wt_vc_ehl <- fixed(0.396); label("Power exponent of WT/73 on EHL product central volume (unitless)") # run57.mod POPV1A exponent 0.396 (Diao 2014)

    # Inter-individual variability (Table 2, variances)
    etalrbase ~ 0.419 # Table 2 omega2 BASE = 0.419
    etalktau ~ 0.657 # Table 2 omega2 ktau = 0.657 (Box-Cox transformed in model)

    # Residual error: additive on log scale; THETA is the SD, SIGMA fixed to 1
    expSd <- 0.276; label("Log-scale residual SD for the reference studies (log IU/dL)") # Table 2 'Standard deviation of residuals (RUV; SD)' = 0.276
  })

  model({
    # Covariate factors (run57.mod $PK). Age effects default to 1 when AGE <= 0.
    cov_age_ksyn <- 1
    cov_age_kdeg <- 1
    cov_age_fixa <- 1
    if (AGE > 0) {
      cov_age_ksyn <- (AGE / 35)^e_age_ksyn
      cov_age_kdeg <- (AGE / 35)^e_age_kdeg
      cov_age_fixa <- (AGE / 35)^e_age_fixa
    }
    cov_form_ksyn <- 1 + e_form_p1_ksyn * FORM_FIDANACOGENE_PROCESS1 + e_form_p2_ksyn * FORM_FIDANACOGENE_PROCESS2
    cov_study_rbase <- 1 + e_study1002_rbase * STUDY_B1821002 + e_study1034_rbase * STUDY_B1821034
    cov_study_fixa <- 1 + e_study1010_fixa * STUDY_B1821010 + e_study1036_fixa * STUDY_B1821036

    # Relative FIX activity between data sources (POPFIXA = 1, no IIV)
    fixa <- cov_study_fixa * cov_age_fixa

    # FIX disposition: L -> dL (x 10), allometric scaling referenced to 70 kg
    rbase <- exp(lrbase + etalrbase) * cov_study_rbase
    cl <- 10 * exp(lcl) * (WT / 70)^e_wt_cl
    vc <- 10 * exp(lvc) * (WT / 70)^e_wt_vc
    q <- 10 * exp(lq) * (WT / 70)^e_wt_cl
    vp <- 10 * exp(lvp) * (WT / 70)^e_wt_vc
    q2 <- 10 * exp(lq2) * (WT / 70)^e_wt_cl
    vp2 <- 10 * exp(lvp2) * (WT / 70)^e_wt_vc

    # Gene therapy rate constants (1/h): days x 24, weeks x 168
    ksyn <- log(2) / (24 * exp(lthalf_ksyn)) * (WT / 70)^e_wt_krate * cov_age_ksyn * cov_form_ksyn
    kdeg <- log(2) / (168 * exp(lthalf_kdeg)) * (WT / 70)^e_wt_krate * cov_age_kdeg
    ktrd <- MIX_TRANSIENT_EXPR * log(2) / (168 * exp(lthalf_ktrd)) * (WT / 70)^e_wt_krate
    ktrp <- MIX_TRANSIENT_EXPR * log(2) / (168 * exp(lthalf_ktrp)) * (WT / 70)^e_wt_krate
    # Box-Cox transformed random effect on the translation constant (IU/vg/h)
    phi_lktau <- (exp(etalktau)^boxcox_lktau - 1) / boxcox_lktau
    ktau <- exp(lktau + phi_lktau) / 1e13 * (WT / 70)^e_wt_ktau

    # Extended half-life FIX product (Diao 2014), referenced to 73 kg
    cl_ehl <- exp(lcl_ehl) * (WT / 73)^e_wt_cl_ehl
    vc_ehl <- exp(lvc_ehl) * (WT / 73)^e_wt_vc_ehl
    q_ehl <- exp(lq_ehl)
    vp_ehl <- exp(lvp_ehl)
    q2_ehl <- exp(lq2_ehl)
    vp2_ehl <- exp(lvp2_ehl)

    # Transgene (vg): productive and non-productive states; translation into
    # FIX does not consume transgene. The non-productive state also degrades
    # at kdeg (run57.mod DADT(9); Table S3 run 56).
    d/dt(transgene) <- -kdeg * transgene - ktrd * transgene + ktrp * transgene_np
    d/dt(transgene_np) <- ktrd * transgene - ktrp * transgene_np - kdeg * transgene_np
    # FIX at the site of synthesis (IU)
    d/dt(fix_cell) <- ktau * transgene - ksyn * fix_cell
    # FIX disposition (IU): endogenous production cl * rbase, transgene-derived
    # FIX, and standard half-life FIX replacement doses into central
    d/dt(central) <- cl * rbase + ksyn * fix_cell -
      q * (central / vc - peripheral1 / vp) -
      q2 * (central / vc - peripheral2 / vp2) -
      cl * central / vc
    d/dt(peripheral1) <- q * (central / vc - peripheral1 / vp)
    d/dt(peripheral2) <- q2 * (central / vc - peripheral2 / vp2)
    # Extended half-life FIX replacement doses (IU)
    d/dt(central_ehl) <- -q_ehl * (central_ehl / vc_ehl - peripheral1_ehl / vp_ehl) -
      q2_ehl * (central_ehl / vc_ehl - peripheral2_ehl / vp2_ehl) -
      cl_ehl * central_ehl / vc_ehl
    d/dt(peripheral1_ehl) <- q_ehl * (central_ehl / vc_ehl - peripheral1_ehl / vp_ehl)
    d/dt(peripheral2_ehl) <- q2_ehl * (central_ehl / vc_ehl - peripheral2_ehl / vp2_ehl)

    # Endogenous steady state at time 0
    central(0) <- vc * rbase
    peripheral1(0) <- vp * rbase
    peripheral2(0) <- vp2 * rbase

    # FIX activity (IU/dL): FIXA-scaled central activity plus the EHL product
    Cc <- fixa * central / vc + central_ehl / vc_ehl
    expSd_study <- expSd * (1 + e_study1002_expSd * STUDY_B1821002 +
      e_study1010_expSd * STUDY_B1821010 +
      e_study1034_expSd * STUDY_B1821034 +
      e_study1036_expSd * STUDY_B1821036 +
      e_studyc037_expSd * STUDY_C037)
    Cc ~ lnorm(expSd_study)
  })
}
