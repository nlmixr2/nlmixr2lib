MedellinGaribay_2015_gentamicin <- function() {
  description <- "Two-compartment IV population PK model for gentamicin in infants 1-24 months (Medellin-Garibay 2015) with linear body-weight scaling on CL and central volume Vc and an additive (CLCR/75)-driven term on CL; intercompartmental clearance Q and peripheral volume Vp are not weight-scaled in the published parameterisation."
  reference <- paste(
    "Medellin-Garibay SE, Rueda-Naharro A, Pena-Cabia S, Garcia B, Romano-Moreno S, Barcia E.",
    "Population pharmacokinetics of gentamicin and dosing optimization for infants.",
    "Antimicrob Agents Chemother. 2015;59(1):482-489.",
    "doi:10.1128/AAC.03464-14.",
    sep = " "
  )
  vignette <- "MedellinGaribay_2015_gentamicin"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight (baseline / time-varying as recorded in the medical record).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Linear (not allometric power) scaling on both CL and Vc per Table 3",
        "footnote b: CL = theta1 * BW + theta5 * (CLCR/75); Vc = theta2 * BW.",
        "Q (theta3) and Vp (theta4) are absolute (no BW scaling) in the published",
        "parameterisation. Population mean 6.4 +/- 2.2 kg (Table 1, n = 208)."
      ),
      source_name        = "BW"
    ),
    CRCL = list(
      description        = paste(
        "Creatinine clearance estimated by the Schwartz equation",
        "CLCR = K * length(cm) / SCr(mg/dL), with K = 0.45 for term infants with",
        "appropriate weight for age, 0.33 for low-weight infants, and 0.55 for",
        "infants > 1 year, applicable to alkaline-picrate (Jaffe) creatinine",
        "assays (including IDMS-traceable forms, as used here). Reported in",
        "mL/min/1.73 m^2 (BSA-normalised)."
      ),
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Used in an additive linear effect on CL with divisive normalisation to",
        "the reference value 75 mL/min/1.73 m^2 (close to the population mean",
        "76.7 +/- 36.9 reported in Table 1). Effect form CL contribution =",
        "theta5 * (CLCR / 75) with theta5 = 0.06 L/h (Table 3). The reference",
        "75 is taken directly from the published equation (Table 3 footnote b)",
        "and is not separately rounded."
      ),
      source_name        = "CLCR"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Postnatal age in months. Screened during covariate model building; not retained in the final two-compartment model. The age range 1-24 months defines the population scope rather than entering as an effect.",
      units       = "months",
      type        = "continuous",
      notes       = "Screened; not retained. See Results 'Base pharmacokinetic models and covariate additions' and Discussion."
    ),
    HT = list(
      description = "Body length (height) in centimetres. Screened during covariate model building; not retained.",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened; not retained. Height enters the model only indirectly via the Schwartz CLCR derivation."
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male). Screened during covariate model building; not retained.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened; not retained (Results 'Base pharmacokinetic models and",
        "covariate additions'; Discussion: 'the covariate of sex did not have",
        "any influence on the performance of the model'). The development",
        "cohort (n = 208) was 44% male, 56% female."
      )
    ),
    BMI = list(
      description = "Body mass index (kg/m^2). Reported as a baseline demographic; not retained as a covariate in the final model.",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Tabulated in Table 1 (population mean 15.8 +/- 2.1) but not retained as a model covariate. Captured here for population description only."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 208L,
    n_studies      = 1L,
    age_range      = "1-24 months",
    age_mean       = "5.8 +/- 4.8 months",
    weight_range   = "Mean 6.4 +/- 2.2 kg (Table 1); individual range not tabulated.",
    weight_mean    = "6.4 +/- 2.2 kg",
    height_mean    = "62.6 +/- 9.1 cm",
    bmi_mean       = "15.8 +/- 2.1 kg/m^2",
    sex_female_pct = 56,
    race_ethnicity = "Not reported (single-centre Spanish cohort).",
    crcl_mean      = "76.7 +/- 36.9 mL/min/1.73 m^2 (Schwartz, Jaffe / IDMS-traceable assay)",
    creat_mean     = "Serum creatinine 0.42 +/- 0.1 mg/dL",
    disease_state  = paste(
      "Hospitalised infants 1-24 months receiving intravenous gentamicin with",
      "routine therapeutic drug monitoring. Indications included infections by",
      "Gram-negative and Gram-positive bacilli, central nervous system /",
      "respiratory / abdominal / urinary / bone / soft-tissue infections,",
      "endocarditis, and septicaemia; gentamicin was used in combination with",
      "ampicillin as empirical therapy for sepsis in newborns and infants."
    ),
    dose_range     = paste(
      "Mean total daily dose 6.0 +/- 1.5 mg/kg/day. Regimens recorded in the",
      "development cohort (n = 208): 4 infants on 2.3 mg/kg q6h; 98 on 2.3",
      "mg/kg q8h; 3 on 3.05 mg/kg q12h; 103 on a single 4.8 mg/kg q24h dose.",
      "Each dose given as a 20-min IV infusion via IVAC syringe pump. The",
      "paper proposes 7 mg/kg q24h with TDM as the optimised regimen for",
      "future practice."
    ),
    regions        = "Spain (Hospital Universitario Severo Ochoa, Leganes, 1990-2011).",
    age_classes    = paste(
      "Pooled-cohort age/weight classes used during covariate exploration:",
      "INF 0 = 1-11 months with low BW (12.1%); INF 1 = 1-11 months with",
      "normal BW (72.9%); INF 2 = > 1 year with normal BW (15%), where 'low'",
      "vs 'normal' BW follows the WHO child growth standards. The final model",
      "does not use these classes as covariates; only continuous BW and CLCR",
      "are retained."
    ),
    n_observations = 335L,
    validation_cohort = paste(
      "Separate external validation cohort of 55 infants (mean age 5.5 +/- 4.4",
      "months, mean BW 6.7 +/- 2.3 kg, mean CLCR 76.7 +/- 44.3 mL/min/1.73 m^2)",
      "with 86 gentamicin observations. Mean prediction error (final model)",
      "-0.2 +/- 1.5 (95% CI -0.6, 0.1) confirms no systematic bias."
    ),
    notes          = paste(
      "Retrospective single-centre cohort study, 1990-2011. The 208-subject",
      "development cohort had 335 serum gentamicin concentrations (mean 4.8",
      "mg/L, range 0.5-15.9 mg/L). Most subjects (96.2%) had 1-2 concentrations",
      "drawn. Baseline demographics from Table 1; dosing detail from Results",
      "'Demographics' paragraph 3."
    )
  )

  ini({
    # Structural population PK parameters (Medellin-Garibay 2015 Table 3 final
    # two-compartment open model with covariates).
    #
    # Table 3 footnote b parameterises CL and Vc as additive / linear functions
    # of body weight rather than the usual allometric (WT/ref)^exponent form:
    #     CL = theta1 * BW + theta5 * (CLCR / 75)
    #     Vc = theta2 * BW
    # Q (theta3) and Vp (theta4) are absolute (L/h and L respectively) and are
    # not weight-scaled in the published parameterisation.
    lcl <- log(0.12);   label("Body-weight coefficient on CL (theta1, L/h per kg BW)")     # Table 3 theta1 = 0.12 +/- 0.01 L/h/kg
    lvc <- log(0.35);   label("Body-weight coefficient on Vc (theta2, L per kg BW)")       # Table 3 theta2 = 0.35 +/- 0.02 L/kg
    lq  <- log(0.23);   label("Intercompartmental clearance Q (theta3, L/h)")              # Table 3 theta3 = 0.23 +/- 0.05 L/h
    lvp <- log(2.3);    label("Peripheral volume of distribution Vp (theta4, L)")          # Table 3 theta4 = 2.3 +/- 1.0 L

    # Additive CLCR effect on CL (Table 3 footnote b, theta5).
    e_crcl_cl <- 0.06;  label("Additive CLCR contribution to CL (theta5, L/h per CLCR/75)") # Table 3 theta5 = 0.06 +/- 0.11 L/h

    # Inter-individual variability (Medellin-Garibay 2015 Table 3 IIV block;
    # reported as CV% in the source). Convert CV% to log-scale variance for
    # log-normal etas: omega^2 = log(1 + CV^2).
    #   IIV CL : 26.4% CV  -> log(1 + 0.264^2) = 0.067374
    #   IIV Vc : 26.5% CV  -> log(1 + 0.265^2) = 0.067869
    # The paper does not report a CL-Vc correlation for the final model, so the
    # block is diagonal.
    etalcl ~ 0.067374   # Table 3 IIV CL 26.4% CV
    etalvc ~ 0.067869   # Table 3 IIV Vc 26.5% CV

    # Residual variability (Medellin-Garibay 2015 Table 3 sigma row; modelled
    # as heteroscedastic / proportional error per Methods 'Pharmacokinetic
    # analysis'). 11.2% CV maps directly to a proportional SD.
    propSd <- 0.112;   label("Proportional residual error (fraction)")    # Table 3 sigma = 11.2 +/- 2.9 %
  })

  model({
    # Individual PK parameters. The structural model from Table 3 footnote b
    # is preserved: CL has BW-linear and CLCR-additive components; Vc is
    # purely BW-linear; Q and Vp are absolute. Log-normal IIV is applied
    # multiplicatively to the typical values for CL and Vc.
    cl <- (exp(lcl) * WT + e_crcl_cl * (CRCL / 75)) * exp(etalcl)
    vc <- exp(lvc) * WT * exp(etalvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment IV PK. Gentamicin is delivered as a 20-min IV infusion
    # in the source paper; the library model does not hard-code the infusion
    # duration, so users can specify rate / dur per dose in their event table.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central                 - k21 * peripheral1

    # Plasma concentration: dose in mg, vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
