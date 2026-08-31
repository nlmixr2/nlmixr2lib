Kim_2026_midazolam_ecmo <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for midazolam and its",
    "primary active metabolite 1-hydroxymidazolam (1-OH-MDZ) in 19 adults",
    "receiving continuous midazolam sedation DURING venoarterial",
    "extracorporeal membrane oxygenation (VA-ECMO) support (Kim 2026,",
    "on-ECMO model). Midazolam: two-compartment disposition whose only",
    "elimination pathway is metabolic clearance to 1-OH-MDZ (CL_MF).",
    "1-OH-MDZ: two-compartment disposition with its own clearance.",
    "Clearance of 1-OH-MDZ increases with ECMO circuit blood flow rate",
    "via a median-normalized proportional effect. Parent-to-metabolite",
    "transfer is encoded 1:1 in midazolam-equivalent mass (the paper",
    "reports no molar conversion factor). The companion post-ECMO model",
    "fit to the same cohort after decannulation is",
    "modellib('Kim_2026_midazolam_postecmo')."
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

  covariateData <- list(
    Q_ECMO = list(
      description        = paste(
        "Blood flow rate delivered through the VA-ECMO circuit."
      ),
      units              = "L/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters 1-OH-midazolam clearance as a median-normalized",
        "proportional effect, CL_1-OH = 49.7 * (1 + 0.336 *",
        "(Q_ECMO - 2.7) / 2.7). NEITHER the functional form NOR the",
        "centering median is stated anywhere in the paper or its",
        "supplements. Supplementary Methods says only that continuous",
        "covariates were 'centered on median values and evaluated using",
        "exponential, power, and proportional models', without saying",
        "which was retained, and Supplementary Table S1 tabulates every",
        "other screened covariate but omits ECMO flow rate, so no median",
        "is published. Both were recovered by fitting the paper's own",
        "Figure 1A, whose five points carry exact numeric labels (1.00,",
        "1.16, 1.34, 1.47, 1.63 at 1-5 L/min); those labels are the",
        "CL_1-OH ratios because the 24-h metabolite concentration is",
        "inversely proportional to CL_1-OH to within about 1%. Of the",
        "three candidate forms only the median-normalized proportional",
        "form fits (RMSE 0.011, vs 0.094 for power; the plain",
        "proportional form would require a negative median of -2.3 L/min",
        "and the exponential form is rejected outright), and solving for",
        "the median at each of the four non-trivial points gives 2.66,",
        "2.47, 2.72, 2.71 L/min. Corroborated by the main text calling",
        "the relationship a 'linear rise' (power and exponential are not",
        "linear in flow) and by the Methods' 'centered on median values'.",
        "See the vignette Errata. Naturally time-varying (flow is",
        "titrated during support and stepped down during weaning); the",
        "paper does not state the time resolution used in the fit.",
        "Operator-ratified canonical (sidecar oare_PMC12777615 q1)."
      ),
      source_name        = "ECMO flow rate"
    )
  )

  # Recorded and evaluated but NOT retained in the final model. All
  # except ECMO_PUMP_SPEED are named in the screened-covariate list of
  # Supplementary Methods, 'Covariate model development'; pump speed is
  # documented only as collected ('Clinical data collection' and 'ECMO
  # system'), and the paper never says it was entered into the stepwise
  # search. No effect estimates are reported for any of these, so they
  # are documentation only and are deliberately absent from model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as a continuous covariate; not retained. Median 53.5 (range 22-87), Supplementary Table S1."
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened as a continuous covariate; not retained. Median 70 (range 52.9-92), Supplementary Table S1."
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
      notes       = "Screened as a continuous covariate; not retained. Median 234 (range 3.8-2883), Supplementary Table S1."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as a continuous covariate; not retained. Median 81 (range 19-1538), Supplementary Table S1."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened as a continuous covariate; not retained. Median 1.2 (range 0.6-4.5), Supplementary Table S1."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as an indicator variable; not retained. 5 of 19 female, Supplementary Table S1."
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
    ECMO_PUMP_SPEED = list(
      description = "ECMO centrifugal-pump rotational speed",
      units       = "RPM",
      type        = "continuous",
      notes       = "Recorded alongside the flow rate (Supplementary Methods, 'Clinical data collection' and 'ECMO system') but absent from the screened-covariate list in 'Covariate model development' and absent from Table 1; only Q_ECMO entered the final model. The reverse of Yang 2017 remifentanil, where pump speed was retained and flow rate was not."
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
    n_subjects     = 19,
    n_studies      = 1,
    age_range      = "22-87 years",
    age_median     = "53.5 years",
    weight_range   = "52.9-92 kg",
    weight_median  = "70 kg",
    sex_female_pct = 26.3,
    disease_state  = paste(
      "Critically ill adults on venoarterial ECMO in a cardiovascular",
      "intensive care unit. Indications: ST-elevation myocardial",
      "infarction (10), atrial fibrillation (3), myocarditis (2),",
      "cardiac arrest (2), non-ST-elevation myocardial infarction (1),",
      "pulmonary embolism (1). Patients with severe hepatic or renal",
      "dysfunction, or taking strong CYP3A4 inhibitors, were excluded."
    ),
    dose_range     = paste(
      "Continuous IV midazolam infusion, rate set by the treating",
      "physician to meet each patient's sedation need, with additional",
      "boluses as required; the study did not fix a dose. Simulations",
      "in the paper span 1-10 mg/h."
    ),
    regions        = "Republic of Korea (Seoul)",
    ecmo_duration  = "median 133.7 h (range 10.6-269)",
    lab_medians    = paste(
      "AST 234 U/L (3.8-2883); ALT 81 U/L (19-1538);",
      "serum creatinine 1.2 mg/dL (0.6-4.5)"
    ),
    notes          = paste(
      "Prospective cohort, ClinicalTrials.gov NCT02581280; Yonsei",
      "University IRB 4-2014-0919. Baseline demographics are",
      "Supplementary Table S1 (on-ECMO support column). The same",
      "patients contribute the post-ECMO model where they survived to",
      "decannulation with continued sampling (11 of 19); the two models",
      "were fit independently. ECMO flow rate itself is NOT tabulated",
      "in Supplementary Table S1."
    )
  )

  ini({
    # ================================================================
    # Midazolam (parent) disposition -- Table 1, 'On ECMO Support'
    # column. Two compartments. CL_MF is the ONLY parent elimination
    # pathway in this model (Supplementary Fig. S1): all midazolam
    # leaving the central compartment becomes 1-OH-midazolam. The
    # paper's limitations paragraph notes that minor routes
    # (4-hydroxylation, glucuronidation) were deliberately excluded.
    # ================================================================
    lvc <- log(36.7)
    label("Midazolam central volume Vc,MDZ (L)")                    # Table 1 On ECMO: Vc,MDZ = 36.7 L (RSE 48.6%; bootstrap 2.15-72.06)
    lvp <- log(57.1)
    label("Midazolam peripheral volume Vp,MDZ (L)")                 # Table 1 On ECMO: Vp,MDZ = 57.1 L (RSE 43.5%; bootstrap 30.54-127.80)
    lq <- log(3.1)
    label("Midazolam intercompartmental clearance Q,MDZ (L/h)")     # Table 1 On ECMO: Q,MDZ = 3.1 L/h (RSE 149.9%; bootstrap 1.67-19.89; the paper flags this RSE as 'somewhat higher')
    lcl_form_1ohm <- log(3.26)
    label("Midazolam-to-1-OH-midazolam formation clearance CL_MF (L/h)")  # Table 1 On ECMO: CL_MF = 3.26 L/h (RSE 18.5%; bootstrap 2.25-4.61)

    # ================================================================
    # 1-OH-midazolam (metabolite) disposition -- Table 1, 'On ECMO
    # Support' column. Two compartments; CL_1-OH is the terminal
    # elimination. Its typical value applies at the reference ECMO
    # flow rate (see e_q_ecmo_cl_1ohm below).
    # ================================================================
    lcl_1ohm <- log(49.7)
    label("1-OH-midazolam clearance at the reference ECMO flow rate (L/h)")  # Table 1 On ECMO: CL_1-OH MDZ = 49.7 L/h (RSE 14.7%; bootstrap 35.06-63.77)
    lvc_1ohm <- log(6.42)
    label("1-OH-midazolam central volume Vc,1-OH (L)")              # Table 1 On ECMO: Vc,1-OH MDZ = 6.42 L (RSE 31.8%; bootstrap 0.16-8.17)
    lvp_1ohm <- log(127)
    label("1-OH-midazolam peripheral volume Vp,1-OH (L)")           # Table 1 On ECMO: Vp,1-OH MDZ = 127 L (RSE 146.1%; bootstrap 34.54-761.86)
    lq_1ohm <- log(0.443)
    label("1-OH-midazolam intercompartmental clearance Q,1-OH (L/h)")  # Table 1 On ECMO: Q,1-OH MDZ = 0.443 L/h (RSE 65.7%; bootstrap 0.1-1.24)

    # ================================================================
    # ECMO flow-rate effect on 1-OH-midazolam clearance.
    #
    # The COEFFICIENT is published; the FORM and the CENTERING MEDIAN
    # are not, and were recovered by fitting the paper's own Fig. 1A
    # (figure-derived, not paper-reported -- see covariateData notes
    # and the vignette Errata). The retained encoding is the
    # median-normalized proportional form
    #
    #   CL_1-OH = 49.7 * (1 + 0.336 * (Q_ECMO - 2.7) / 2.7)
    #
    # which reproduces all five Fig. 1A labels to RMSE 0.011:
    #   flow 1 2 3 4 5 L/min -> ratio vs 1 L/min
    #   predicted 1.000 1.158 1.316 1.474 1.631
    #   published 1.00  1.16  1.34  1.47  1.63
    # ================================================================
    e_q_ecmo_cl_1ohm <- 0.336
    label("ECMO flow-rate proportional coefficient on 1-OH-midazolam CL (unitless)")  # Table 1 On ECMO: theta ECMO flow rate on CL1-OH MDZ = 0.336 (RSE 93.5%; bootstrap 0.30-0.51)

    # ================================================================
    # Inter-individual variability. Table 1's 'Random effects' block is
    # on the VARIANCE scale: the parenthetical for every one of the
    # nine random-effect entries in the table is exactly
    # 100 * sqrt(omega^2), i.e. the footnote's 'CV %' is the
    # approximate log-normal CV, not a relative standard error.
    # Check: 0.246 -> 49.6, 0.497 -> 70.5, 0.693 -> 83.3,
    #        0.176 -> 42.0, 0.137 -> 37.0 (all five match the table).
    # So the tabulated numbers are used directly as variances.
    # ================================================================
    etalcl_1ohm ~ 0.246       # Table 1 On ECMO: omega^2 CL,1-OH MDZ = 0.246 (49.6; bootstrap 0.09-0.52)
    etalcl_form_1ohm ~ 0.497  # Table 1 On ECMO: omega^2 CLmf = 0.497 (70.5; bootstrap 0.20-0.94)
    etalvc ~ 0.693            # Table 1 On ECMO: omega^2 Vc,MDZ = 0.693 (83.3; bootstrap 0.02-3.4)
    etalvp ~ 0.176            # Table 1 On ECMO: omega^2 Vp,MDZ = 0.176 (42.0; bootstrap 0-2.38)
    etalvc_1ohm ~ 0.137       # Table 1 On ECMO: omega^2 Vc,1-OH MDZ = 0.137 (37.0; bootstrap 0.01-5.65)

    # Table 1 carries an 'ETA correlation' row labelled
    # 'omega^2 CL1-OH MDZ, omega^2 Q1-OH MDZ' = 0.228 (bootstrap
    # 0-0.56). It is NOT encoded here: it names an eta on Q,1-OH MDZ,
    # and no variance for that eta is reported anywhere in the paper,
    # so the correlation cannot be built without inventing one.
    # (Operator sidecar oare_PMC12777615 q3.) Recorded verbatim in the
    # vignette Errata so a reader holding the control stream can
    # restore it.

    # ================================================================
    # Residual variability. Supplementary Methods gives the error model
    # as Cij = Cpred,ij * (1 + eps_ij) with Var(eps) = sigma^2, and
    # Table 1 tabulates sigma^2 -- so propSd = sqrt(sigma^2).
    # Cross-check that these are variances and not SDs: the RSE of an
    # estimated variance is roughly sqrt(2/N), and 13.3% implies about
    # 113 observations, which is plausible for 19 serially sampled
    # patients; an SD scale would imply about 28 observations, too few.
    # ================================================================
    propSd <- sqrt(0.404)
    label("Proportional residual error, midazolam (fraction)")           # Table 1 On ECMO: sigma^2 Parent = 0.404 (RSE 13.3%; bootstrap 0.30-0.51)
    propSd_1ohm <- sqrt(0.362)
    label("Proportional residual error, 1-OH-midazolam (fraction)")      # Table 1 On ECMO: sigma^2 Metabolite = 0.362 (RSE 8.5%; bootstrap 0.29-0.41)
  })

  model({
    # ------------------------------------------------------------
    # Reference ECMO flow rate for the covariate centering.
    # FIGURE-DERIVED, not paper-reported: back-solved from Fig. 1A
    # (see ini() comment on e_q_ecmo_cl_1ohm).
    # ------------------------------------------------------------
    q_ecmo_ref <- 2.7   # L/min

    # ------------------------------------------------------------
    # Individual parameters, midazolam (parent).
    # ------------------------------------------------------------
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q  <- exp(lq)
    cl_form_1ohm <- exp(lcl_form_1ohm + etalcl_form_1ohm)

    # ------------------------------------------------------------
    # Individual parameters, 1-OH-midazolam (metabolite). ECMO flow
    # rate acts only on the metabolite's clearance; nothing in the
    # parent's disposition depends on flow.
    # ------------------------------------------------------------
    ecmo_factor <- 1 + e_q_ecmo_cl_1ohm * (Q_ECMO - q_ecmo_ref) / q_ecmo_ref
    cl_1ohm <- exp(lcl_1ohm + etalcl_1ohm) * ecmo_factor
    vc_1ohm <- exp(lvc_1ohm + etalvc_1ohm)
    vp_1ohm <- exp(lvp_1ohm)
    q_1ohm  <- exp(lq_1ohm)

    # ------------------------------------------------------------
    # ODE system (Supplementary Fig. S1). Midazolam is dosed by IV
    # infusion into `central`. Its whole elimination flux becomes the
    # 1-OH-midazolam formation flux, transferred 1:1 in
    # midazolam-equivalent mass -- the paper states no molar
    # conversion factor, the same convention used by
    # Franken_2017_midazolam. A true molar correction (MW 341.8 /
    # 325.8 = 1.049) would rescale CL_1-OH and the 1-OH volumes by
    # under 5%; see the vignette Errata.
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
