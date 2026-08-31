Kim_2026_midazolam_postecmo <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for midazolam and its",
    "primary active metabolite 1-hydroxymidazolam (1-OH-MDZ) in 11 adults",
    "receiving continuous midazolam sedation AFTER venoarterial",
    "extracorporeal membrane oxygenation (VA-ECMO) was discontinued",
    "(Kim 2026, post-ECMO model). Structurally identical to the on-ECMO",
    "model: midazolam two-compartment disposition whose only elimination",
    "pathway is metabolic clearance to 1-OH-MDZ (CL_MF), feeding a",
    "two-compartment 1-OH-MDZ disposition. No covariate is retained,",
    "because the ECMO circuit is no longer present. Parent-to-metabolite",
    "transfer is encoded 1:1 in midazolam-equivalent mass (the paper",
    "reports no molar conversion factor). The companion on-ECMO model is",
    "modellib('Kim_2026_midazolam_ecmo'); the two were fit independently."
  )
  reference <- paste(
    "Kim H, Jin BH, Yang S, Hahn J, Kang S, Kim D, Lee H, Kwack H,",
    "Chae SU, Bae SK, Wi J, Chang MJ.",
    "Effect of extracorporeal membrane oxygenation flow rate on",
    "midazolam clearance: a population pharmacokinetic study.",
    "Anesthesiology. 2026;144(3):485-488.",
    "doi:10.1097/ALN.0000000000005811.",
    sep = " "
  )
  vignette <- "Kim_2026_midazolam"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  # No covariate is retained in the post-ECMO model (Table 1 shows '-'
  # for the only covariate coefficient in the paper), so there is no
  # covariateData. The demographic and laboratory entries below are the
  # screened-covariate list of Supplementary Methods, 'Covariate model
  # development'; Q_ECMO is listed separately because it is structurally
  # inapplicable after decannulation rather than screened-and-dropped.
  # No effect estimate is reported for any of them, so they are
  # documentation only and are deliberately absent from model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as a continuous covariate; not retained. Median 56 (range 22-87), Supplementary Table S1."
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened as a continuous covariate; not retained. Median 70 (range 52.9-81.4), Supplementary Table S1."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Listed among the screened continuous covariates in Supplementary Methods; not retained, and not tabulated in Supplementary Table S1."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as a continuous covariate; not retained. Median 63.9 (range 21-356), Supplementary Table S1."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as a continuous covariate; not retained. Median 44 (range 12-93), Supplementary Table S1."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened as a continuous covariate; not retained. Median 1.0 (range 0.6-1.7), Supplementary Table S1."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as an indicator variable; not retained. 3 of 11 female, Supplementary Table S1."
    ),
    CONMED_SUFENTANIL = list(
      description = "Concomitant sufentanil sedation",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as an indicator variable (Supplementary Methods); not retained. Counts not reported."
    ),
    CONMED_REMIFENTANIL = list(
      description = "Concomitant remifentanil sedation",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as an indicator variable (Supplementary Methods); not retained. Counts not reported."
    ),
    Q_ECMO = list(
      description = "ECMO circuit blood flow rate",
      units       = "L/min",
      type        = "continuous",
      notes       = "Retained in the companion on-ECMO model but structurally inapplicable here: these patients have been decannulated, so there is no circuit flow. Table 1 shows '-' for theta ECMO flow rate on CL1-OH MDZ in the 'After ECMO Support' column."
    )
  )

  compartmentData <- list(
    central           = list(analyte = "midazolam",          units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1       = list(analyte = "midazolam",          units = "mg", specimen = "plasma", verified = TRUE),
    central_1ohm      = list(analyte = "1-OH-midazolam",     units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_1ohm  = list(analyte = "1-OH-midazolam",     units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 11,
    n_studies      = 1,
    age_range      = "22-87 years",
    age_median     = "56 years",
    weight_range   = "52.9-81.4 kg",
    weight_median  = "70 kg",
    sex_female_pct = 27.3,
    disease_state  = paste(
      "Critically ill adults in a cardiovascular intensive care unit",
      "after decannulation from venoarterial ECMO, still receiving",
      "midazolam sedation. Indications for the preceding VA-ECMO run:",
      "ST-elevation myocardial infarction (4), atrial fibrillation (3),",
      "myocarditis (2), cardiac arrest (2). Patients with severe hepatic",
      "or renal dysfunction, or taking strong CYP3A4 inhibitors, were",
      "excluded."
    ),
    dose_range     = paste(
      "Continuous IV midazolam infusion, rate set by the treating",
      "physician to meet each patient's sedation need, with additional",
      "boluses as required; the study did not fix a dose."
    ),
    regions        = "Republic of Korea (Seoul)",
    lab_medians    = paste(
      "AST 63.9 U/L (21-356); ALT 44 U/L (12-93);",
      "serum creatinine 1.0 mg/dL (0.6-1.7)"
    ),
    notes          = paste(
      "Prospective cohort, ClinicalTrials.gov NCT02581280; Yonsei",
      "University IRB 4-2014-0919. Baseline demographics are",
      "Supplementary Table S1 (post-ECMO support column). These 11",
      "patients are the subset of the 19-patient on-ECMO cohort with",
      "post-decannulation sampling; Supplementary Methods states that",
      "'two independent population pharmacokinetic models were",
      "developed', so this model does not share parameters with the",
      "on-ECMO model."
    )
  )

  ini({
    # ================================================================
    # Midazolam (parent) disposition -- Table 1, 'After ECMO Support'
    # column. Two compartments. CL_MF is the ONLY parent elimination
    # pathway in this model (Supplementary Fig. S1): all midazolam
    # leaving the central compartment becomes 1-OH-midazolam.
    # ================================================================
    lvc <- log(8.25)
    label("Midazolam central volume Vc,MDZ (L)")                    # Table 1 After ECMO: Vc,MDZ = 8.25 L (RSE 46.6%; bootstrap 0.47-15.54)
    lvp <- log(37.4)
    label("Midazolam peripheral volume Vp,MDZ (L)")                 # Table 1 After ECMO: Vp,MDZ = 37.4 L (RSE 41.0%; bootstrap 14.25-74.28)
    lq <- log(7.88)
    label("Midazolam intercompartmental clearance Q,MDZ (L/h)")     # Table 1 After ECMO: Q,MDZ = 7.88 L/h (RSE 108.9%; bootstrap 1.28-34.91; the paper flags this RSE as 'somewhat higher')
    lcl_form_1ohm <- log(5.44)
    label("Midazolam-to-1-OH-midazolam formation clearance CL_MF (L/h)")  # Table 1 After ECMO: CL_MF = 5.44 L/h (RSE 17.3%; bootstrap 2.33-6.01)

    # ================================================================
    # 1-OH-midazolam (metabolite) disposition -- Table 1, 'After ECMO
    # Support' column. Two compartments; CL_1-OH is the terminal
    # elimination. No covariate acts on it in this model.
    # ================================================================
    lcl_1ohm <- log(56.6)
    label("1-OH-midazolam clearance (L/h)")                         # Table 1 After ECMO: CL_1-OH MDZ = 56.6 L/h (RSE 18.8%; bootstrap 25.33-67.09)
    lvc_1ohm <- log(2.71)
    label("1-OH-midazolam central volume Vc,1-OH (L)")              # Table 1 After ECMO: Vc,1-OH MDZ = 2.71 L (RSE 26.1%; bootstrap 0.56-3.33)
    lvp_1ohm <- log(98)
    label("1-OH-midazolam peripheral volume Vp,1-OH (L)")           # Table 1 After ECMO: Vp,1-OH MDZ = 98 L (RSE 37.0%; bootstrap 27.96-170.01)
    lq_1ohm <- log(1.14)
    label("1-OH-midazolam intercompartmental clearance Q,1-OH (L/h)")  # Table 1 After ECMO: Q,1-OH MDZ = 1.14 L/h (RSE 43.0%; bootstrap 0.54-2.46)

    # ================================================================
    # Inter-individual variability. As in the on-ECMO model, Table 1's
    # 'Random effects' block is on the VARIANCE scale: the
    # parenthetical is 100 * sqrt(omega^2) for every entry.
    # Check: 0.574 -> 75.8, 0.732 -> 85.6, 1.08 -> 103.9,
    #        0.118 -> 34.4 (all four match the table).
    # ================================================================
    etalcl_1ohm ~ 0.574       # Table 1 After ECMO: omega^2 CL,1-OH MDZ = 0.574 (75.8; bootstrap 0.25-0.94)
    etalcl_form_1ohm ~ 0.732  # Table 1 After ECMO: omega^2 CLmf = 0.732 (85.6; bootstrap 0.40-1.50)
    etalvp ~ 1.08             # Table 1 After ECMO: omega^2 Vp,MDZ = 1.08 (103.9; bootstrap 0.31-4.13)
    etalvc_1ohm ~ 0.118       # Table 1 After ECMO: omega^2 Vc,1-OH MDZ = 0.118 (34.4; bootstrap 0.06-0.53)

    # Table 1 reports omega^2 Vc,MDZ as '0 FIX' in the 'After ECMO
    # Support' column: the IIV on midazolam central volume was
    # estimated in the on-ECMO model but held at zero here. It is
    # omitted rather than written as `etalvc ~ fixed(0)` because a
    # zero-variance diagonal makes OMEGA singular and breaks the
    # Cholesky sampler used by rxSolve (same handling as
    # Wattanakul_2024_primaquine_motherinfant). Mathematically
    # identical: a zero-variance eta contributes nothing.

    # The on-ECMO model's 'ETA correlation' row has no counterpart
    # here -- Table 1 shows '-' for the 'After ECMO Support' column.

    # ================================================================
    # Residual variability. Supplementary Methods gives the error model
    # as Cij = Cpred,ij * (1 + eps_ij) with Var(eps) = sigma^2, and
    # Table 1 tabulates sigma^2 -- so propSd = sqrt(sigma^2).
    # ================================================================
    propSd <- sqrt(0.525)
    label("Proportional residual error, midazolam (fraction)")           # Table 1 After ECMO: sigma^2 Parent = 0.525 (RSE 17.0%; bootstrap 0.30-0.65)
    propSd_1ohm <- sqrt(0.359)
    label("Proportional residual error, 1-OH-midazolam (fraction)")      # Table 1 After ECMO: sigma^2 Metabolite = 0.359 (RSE 11.4%; bootstrap 0.27-0.43)
  })

  model({
    # ------------------------------------------------------------
    # Individual parameters, midazolam (parent). No IIV on vc in this
    # model (omega^2 Vc,MDZ = 0 FIX; see ini()).
    # ------------------------------------------------------------
    vc <- exp(lvc)
    vp <- exp(lvp + etalvp)
    q  <- exp(lq)
    cl_form_1ohm <- exp(lcl_form_1ohm + etalcl_form_1ohm)

    # ------------------------------------------------------------
    # Individual parameters, 1-OH-midazolam (metabolite).
    # ------------------------------------------------------------
    cl_1ohm <- exp(lcl_1ohm + etalcl_1ohm)
    vc_1ohm <- exp(lvc_1ohm + etalvc_1ohm)
    vp_1ohm <- exp(lvp_1ohm)
    q_1ohm  <- exp(lq_1ohm)

    # ------------------------------------------------------------
    # ODE system (Supplementary Fig. S1). Midazolam is dosed by IV
    # infusion into `central`. Its whole elimination flux becomes the
    # 1-OH-midazolam formation flux, transferred 1:1 in
    # midazolam-equivalent mass -- the paper states no molar
    # conversion factor, the same convention used by
    # Franken_2017_midazolam. See the vignette Errata.
    # ------------------------------------------------------------
    d/dt(central) <- -cl_form_1ohm * central / vc -
                      q * (central / vc - peripheral1 / vp)
    d/dt(peripheral1) <- q * (central / vc - peripheral1 / vp)

    d/dt(central_1ohm) <- cl_form_1ohm * central / vc -
                          cl_1ohm * central_1ohm / vc_1ohm -
                          q_1ohm * (central_1ohm / vc_1ohm - peripheral1_1ohm / vp_1ohm)
    d/dt(peripheral1_1ohm) <- q_1ohm * (central_1ohm / vc_1ohm - peripheral1_1ohm / vp_1ohm)

    # ------------------------------------------------------------
    # Observations. States hold mg (of midazolam equivalents) and
    # volumes are L, so amount/volume is mg/L; the assay reports
    # ng/mL, and 1 mg/L = 1000 ng/mL.
    # ------------------------------------------------------------
    Cc <- central / vc * 1000
    Cc_1ohm <- central_1ohm / vc_1ohm * 1000

    Cc ~ prop(propSd)
    Cc_1ohm ~ prop(propSd_1ohm)
  })
}
