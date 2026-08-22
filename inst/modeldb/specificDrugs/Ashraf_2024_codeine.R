Ashraf_2024_codeine <- function() {
  description <- paste(
    "Joint parent + three-metabolite population PK model for oral codeine",
    "and its metabolites morphine, codeine-6-glucuronide (C6G) and",
    "morphine-3-glucuronide (M3G) in 997 ambulatory surgical patients given",
    "a single preoperative 60 mg oral codeine dose (Ashraf 2024). Codeine is",
    "described by a one-compartment disposition with a first-order",
    "absorption depot and an estimated bioavailability. Total codeine",
    "elimination is split into a CYP2D6-mediated branch (ke * fm_morphine) that",
    "forms morphine and a non-CYP2D6 (glucuronidation) branch",
    "(ke * (1 - fm_morphine)) that forms C6G; morphine elimination is in turn",
    "split into an M3G-forming branch (ke_morphine * fm_m3g, fm_m3g fixed at 0.60",
    "from the literature) and other pathways. Morphine, C6G and M3G each",
    "have a one-compartment disposition. The CYP2D6 activity score (CPIC sum",
    "of allele activity values) enters as an ordinal continuous covariate on",
    "the codeine-to-morphine metabolic fraction through an exponential model",
    "referenced to an activity score of 2, followed by the authors'",
    "f / (1 + f) rescaling that constrains the fraction to (0, 1). Codeine",
    "and morphine clearance and central volume are allometrically scaled to",
    "the population median body weight of 80 kg with fixed exponents 0.75",
    "and 1; the glucuronide parameters are not weight-scaled."
  )
  reference <- paste(
    "Ashraf MW, Poikola S, Neuvonen M, Kiiski JI, Kontinen VK, Olkkola KT,",
    "Backman JT, Niemi M, Saari TI. Population Pharmacokinetic",
    "Quantification of CYP2D6 Activity in Codeine Metabolism in Ambulatory",
    "Surgical Patients for Model-Informed Precision Dosing.",
    "Clin Pharmacokinet. 2024;63(11):1547-1560.",
    "doi:10.1007/s40262-024-01433-9.",
    sep = " "
  )
  vignette <- "Ashraf_2024_codeine"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "mg/L"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. The model carries no metabolite-to-parent molecular
  # weight ratio (Fig. 1 draws none on any transfer arrow and the paper
  # reports no molecular weights), so every state is an amount in codeine
  # mass equivalents and each metabolite volume is an apparent volume that
  # absorbs the corresponding molecular-weight ratio. Concentrations
  # nonetheless come out on the measured (assay) scale because the
  # metabolite clearances and volumes were estimated against the measured
  # metabolite concentrations.
  compartmentData <- list(
    depot       = list(analyte = "codeine", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "codeine", units = "mg", specimen = "plasma", verified = TRUE),
    central_morphine = list(analyte = "morphine", units = "mg", specimen = "plasma", verified = TRUE),
    central_c6g = list(analyte = "codeine-6-glucuronide", units = "mg", specimen = "plasma", verified = TRUE),
    central_m3g = list(analyte = "morphine-3-glucuronide", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Enters codeine and morphine clearance and",
        "central volume a priori via CL_i = CL_80kg * (WT / 80)^0.75 and",
        "V_i = V_80kg * (WT / 80)^1 (Online Resource 1 Section 1.6.1,",
        "Eqs. 3-6). The reference weight is the study population median,",
        "80 kg. Allometric exponents are fixed at the canonical 0.75 and 1",
        "and were not estimated. The C6G and M3G clearance and volume",
        "parameters are NOT weight-scaled: the paper restricts allometric",
        "scaling to 'both codeine and morphine' (paper Section 2.5;",
        "Online Resource 1 Section 1.6.1)."
      ),
      source_name        = "WTKG (Online Resource 1 Eqs. 3-4)"
    ),
    CYP2D6 = list(
      description        = paste(
        "CYP2D6 activity score (AS): the CPIC consensus sum of per-allele",
        "activity values for the subject's diplotype."
      ),
      units              = "(activity score, 0-4)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject (germline genotype). Scoring per paper",
        "Table 1 footnote a: 1 for each normal-function allele (*1, *2),",
        "0.5 for each decreased-function allele (*9, *17, *29, *41), 0.25",
        "for each *10 allele, and 0 for each no-function allele (*3, *4,",
        "*5, *6, *40); duplicated normal-function alleles combined with",
        "*10 score 2.25, and other duplications score 3 or 4. Observed",
        "values in the study population were 0, 0.25, 0.5, 0.75, 1, 1.25,",
        "1.5, 2, 2.25, 3 and 4. The model reference value is AS = 2 (the",
        "modal normal-metaboliser score, paper Section 3.2). Used as an",
        "ordinal continuous covariate on the codeine-to-morphine metabolic",
        "fraction; it is NOT dichotomised into the PM / IM / NM / UM",
        "phenotype classes (the paper reports that the activity-score",
        "parameterisation reduced unexplained variability in the metabolic",
        "ratio from 45% to 33% relative to the phenotype-class model)."
      ),
      source_name        = "AS (paper Section 3.2 and Table 1 footnote a)"
    )
  )

  # Screened by the PsN stepwise-covariate-modelling protocol but not
  # retained in the final model (paper Section 3.2; Online Resource 1
  # Section 2.4: no improvement in predictive performance). Documented here
  # so the paper's covariate screen is preserved without declaring unused
  # entries in covariateData. The American Society of Anesthesiologists
  # (ASA) physical status class and the number of cigarettes smoked per day
  # were also screened and rejected; neither has a canonical register entry
  # and neither is used by this model, so they are recorded in the vignette
  # rather than here.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened on codeine clearance with power and exponential models in",
        "the PsN SCM run; not retained (Online Resource 1 Section 2.4)."
      )
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened in the PsN SCM run; not retained."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = paste(
        "Reported in paper Table 1 baseline demographics; not tested as a",
        "model covariate (body weight was used as the size descriptor)."
      )
    ),
    SMOKE_CURRENT = list(
      description = "Current cigarette smoking",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened on codeine clearance (increased glucuronidation) and on",
        "morphine clearance in the PsN SCM run; not retained",
        "(Online Resource 1 Section 2.4)."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 997,
    n_studies      = 1,
    age_mean_sd    = "47.8 (12.9) years",
    weight_mean_sd = "79.8 (13.8) kg",
    weight_median  = "80 kg (allometric reference)",
    bmi_mean_sd    = "25.9 (3.60) kg/m^2",
    sex_female_pct = NULL,
    race_ethnicity = NULL,
    disease_state  = paste(
      "Adults scheduled for elective ambulatory (day-case) surgery.",
      "American Society of Anesthesiologists physical status class 1",
      "(n = 513, 51%), 2 (n = 424, 43%) and 3 (n = 60, 6%). Use of strong",
      "CYP2D6 inhibitors was an exclusion criterion."
    ),
    dose_range     = paste(
      "Single preoperative oral dose of codeine 60 mg given as a",
      "paracetamol 1000 mg + codeine 60 mg combination tablet. Model",
      "simulations in the paper additionally cover 30 mg and 60 mg codeine",
      "given once to four times daily for up to four days."
    ),
    regions        = "Finland (Jorvi Hospital day-surgery unit, Espoo)",
    genotype       = paste(
      "CYP2D6 activity score distribution (n, % of 997): AS 0, 37 (3.7);",
      "0.25, 5 (0.5); 0.5, 21 (2.1); 0.75, 2 (0.2); 1, 240 (24); 1.25,",
      "23 (2.3); 1.5, 67 (6.7); 2, 537 (54); 2.25, 2 (0.2); 3, 61 (6.1);",
      "4, 2 (0.2). Predicted metaboliser status: 37 poor, 268 intermediate,",
      "629 normal and 65 ultrarapid. A total of 64 distinct CYP2D6",
      "genotypes were observed."
    ),
    notes          = paste(
      "Demographics from paper Table 1. Prospective clinical trial",
      "(EudraCT 2015-005561-23), 1000 patients recruited, 997 with",
      "concentration data for one or more analytes. Sampling was sparse:",
      "exactly two blood samples per patient, the first 20-60 min and the",
      "second 180-360 min after the paracetamol-codeine tablet. Assay was",
      "LC-MS/MS with a 0.05 ng/mL lower limit of quantification (0.1 ng/mL",
      "for morphine-6-glucuronide) and day-to-day CV below 11%. NONMEM",
      "7.4.3+ with ADVAN13, FOCE-I and MU-referencing in the log domain;",
      "parameter precision from a sampling-importance-resampling run."
    )
  )

  ini({
    # ============================================================
    # Allometric exponents on body weight -- fixed a priori at the
    # canonical 0.75 (CL) and 1 (V) values. Online Resource 1
    # Section 1.6.1, Eqs. 3-4: WT_CL = (WTKG_i / 80)^0.75 and
    # WT_V = (WTKG_i / 80). Applied to codeine and morphine CL and Vc
    # only, not to the glucuronides.
    # ============================================================
    e_wt_cl_q <- fixed(0.75)
    label("Allometric exponent on codeine and morphine CL (unitless)")
    # Online Resource 1 Eq. 3; reported without uncertainty.

    e_wt_vc_vp <- fixed(1)
    label("Allometric exponent on codeine and morphine Vc (unitless)")
    # Online Resource 1 Eq. 4; reported without uncertainty.

    # ============================================================
    # Codeine disposition -- paper Table 2, 'Parameter estimate'
    # median column. Values are for an 80 kg subject.
    # ============================================================
    lka <- log(8.74)
    label("Codeine absorption rate constant (1/h)")
    # Table 2: k_a,cod = 8.74 1/h (95% CI 6.49-10.98).

    lcl <- log(110)
    label("Codeine total systemic clearance at 80 kg (L/h)")
    # Table 2: CL_cod = 110 L/h (95% CI 59.6-160.4).

    lvc <- log(427.5)
    label("Codeine central volume at 80 kg (L)")
    # Table 2: V_c,cod = 427.5 L (95% CI 231.2-623.8).

    lfdepot <- log(0.84)
    label("Codeine oral bioavailability (fraction)")
    # Table 2: F_cod = 0.84 (95% CI 0.4556-1.00). Estimated, not fixed.

    # ============================================================
    # Morphine disposition -- paper Table 2.
    # ============================================================
    lcl_morphine <- log(357.5)
    label("Morphine total systemic clearance at 80 kg (L/h)")
    # Table 2: CL_mor = 357.5 L/h (95% CI 178.2-536.9). Apparent value in
    # codeine mass equivalents (no molecular-weight ratio in the model).

    lvc_morphine <- log(22.8)
    label("Morphine central volume at 80 kg (L)")
    # Table 2: V_c,mor = 22.8 L (95% CI 9.447-36.2).

    # ============================================================
    # Codeine-6-glucuronide disposition -- paper Table 2. Not
    # allometrically scaled (paper Section 2.5 restricts weight scaling
    # to codeine and morphine).
    # ============================================================
    lcl_c6g <- log(7.96)
    label("C6G total systemic clearance (L/h)")
    # Table 2: CL_C6G = 7.96 L/h (95% CI 4.28-11.66).

    lvc_c6g <- log(5.36)
    label("C6G central volume (L)")
    # Table 2: V_c,C6G = 5.36 L (95% CI 2.88-7.86).

    # ============================================================
    # Morphine-3-glucuronide disposition -- paper Table 2.
    # ============================================================
    lcl_m3g <- log(9.45)
    label("M3G total systemic clearance (L/h)")
    # Table 2: CL_M3G = 9.45 L/h (95% CI 4.67-14.21).

    lvc_m3g <- log(8.47)
    label("M3G central volume (L)")
    # Table 2: V_c,M3G = 8.47 L (95% CI 4.17-12.78).

    # ============================================================
    # Metabolic fractions.
    # ============================================================
    lfm_morphine <- log(0.16)
    label("Codeine fraction metabolised to morphine at CYP2D6 AS 2, before the f/(1+f) rescaling (fraction)")
    # Table 2: f_morphine = 0.16 (95% CI 0.108-0.209), footnote b 'Scaled to
    # activity score 2'. This is theta_A of Online Resource 1 Eq. 11; the
    # fraction actually entering the ODEs is f/(1+f) = 0.138 at AS 2
    # (Online Resource 1 Eq. 14).

    lfm_m3g <- fixed(log(0.60))
    label("Fraction of morphine elimination forming M3G, from the literature (fraction)")
    # Paper Section 3.1 and Online Resource 1 Section 2.1: 'the metabolic
    # fraction of M3G from morphine elimination was fixed at 60%,
    # according to a recent report'. Not estimated; no uncertainty
    # reported.

    # ============================================================
    # CYP2D6 activity-score covariate effect on the codeine-to-morphine
    # metabolic fraction. Online Resource 1 Eq. 11:
    #   f_morphine,pop = theta_BASE * exp(theta_eff * (AS_i - AS_ref))
    # with AS_ref = 2.
    # ============================================================
    e_cyp2d6_fm_morphine <- 1.00
    label("Exponential CYP2D6 activity-score effect on the codeine-to-morphine fraction (per activity-score unit)")
    # Table 2: GEN_eff = 1.00 (95% CI 0.96-1.05). Natural (untransformed)
    # scale: the confidence interval is symmetric about 1.

    # ============================================================
    # Inter-individual variability. Paper Table 2 reports the eta rows as
    # standard deviations (omega), not variances: the final eta on the
    # codeine-to-morphine ratio is 0.334 and the text (Section 3.2)
    # reports the corresponding unexplained variability as 33%. The
    # variances below are therefore the squares of the tabulated values.
    # No covariances were reported. IIV was estimated on k_a, V_c and F
    # for codeine, on morphine clearance, and on the metabolic ratio only
    # (paper Section 3.1; Online Resource 1 Section 2.1).
    # ============================================================
    etalka ~ 3.75^2
    # Table 2: eta k_a,cod = 3.75 (95% CI 2.96-4.55). Very large on the
    # log scale; the design has only two samples per subject, the first
    # 20-60 min post-dose, so the absorption rate is barely identifiable
    # per individual. The paper's own codeine VPC (Fig. 2A) spans roughly
    # two orders of magnitude at the earliest sampling times.

    etalvc ~ 0.181^2
    # Table 2: eta V_c,cod = 0.181 (95% CI 0.158-0.204).

    etalfdepot ~ 0.024^2
    # Table 2: eta F_cod = 0.024 (95% CI 0.019-0.029). Small because F and
    # V_c are confounded in an oral-only design; the codeine exposure
    # variability is carried by V_c.

    etalcl_morphine ~ 0.062^2
    # Table 2: eta CL_mor = 0.062 (95% CI 0.046-0.077).

    etalfm_morphine ~ 0.334^2
    # Table 2: eta R_morphine = 0.334 (95% CI 0.306-0.361), described in the
    # text as 33% unexplained variability in the metabolic ratio after
    # the activity-score covariate was added. Applied to f_morphine BEFORE the
    # f/(1+f) rescaling (Online Resource 1 Eqs. 13-14).

    # ============================================================
    # Residual variability. The paper describes an additive error model
    # (paper Section 2.4; Online Resource 1 Eq. 16, Y = IPRED + EPS)
    # applied to log-transformed observations, i.e. a log-normal
    # (exponential) residual in linear concentration space. A literally
    # additive residual in ng/mL is ruled out by magnitude: 0.103 ng/mL
    # for C6G would be about 0.01% of the observed C6G concentrations
    # (of order 1000 ng/mL, Fig. 2C), whereas the log-scale reading gives
    # 10.3% -- consistent with the reported assay day-to-day CV of below
    # 11%. See vignette Assumptions and deviations.
    # ============================================================
    expSd <- 0.259
    label("Log-scale residual SD for codeine (fraction)")
    # Table 2: epsilon_cod = 0.259 (95% CI 0.242-0.275).

    expSd_morphine <- 0.149
    label("Log-scale residual SD for morphine (fraction)")
    # Table 2: epsilon_mor = 0.149 (95% CI 0.136-0.162).

    expSd_c6g <- 0.103
    label("Log-scale residual SD for C6G (fraction)")
    # Table 2: epsilon_C6G = 0.103 (95% CI 0.094-0.112).

    expSd_m3g <- 0.096
    label("Log-scale residual SD for M3G (fraction)")
    # Table 2: epsilon_M3G = 0.096 (95% CI 0.087-0.105).
  })

  model({
    # ------------------------------------------------------------
    # Reference covariate values.
    # ------------------------------------------------------------
    wt_ref <- 80   # population median body weight (kg), Online Resource 1 Eqs. 3-4
    as_ref <- 2    # reference CYP2D6 activity score, paper Section 3.2

    allo_cl_factor <- (WT / wt_ref)^e_wt_cl_q
    allo_v  <- (WT / wt_ref)^e_wt_vc_vp

    # ------------------------------------------------------------
    # Individual codeine parameters. IIV on k_a, V_c and F only;
    # codeine clearance carries no IIV (paper Table 2).
    # ------------------------------------------------------------
    ka     <- exp(lka + etalka)
    cl     <- exp(lcl) * allo_cl_factor
    vc     <- exp(lvc + etalvc) * allo_v
    fdepot <- exp(lfdepot + etalfdepot)

    # ------------------------------------------------------------
    # Individual metabolite parameters. Only morphine clearance carries
    # IIV. The glucuronide parameters are not weight-scaled.
    # ------------------------------------------------------------
    cl_morphine <- exp(lcl_morphine + etalcl_morphine) * allo_cl_factor
    vc_morphine <- exp(lvc_morphine) * allo_v
    cl_c6g      <- exp(lcl_c6g)
    vc_c6g      <- exp(lvc_c6g)
    cl_m3g      <- exp(lcl_m3g)
    vc_m3g      <- exp(lvc_m3g)

    # ------------------------------------------------------------
    # CYP2D6 activity-score effect on the codeine-to-morphine metabolic
    # fraction. Online Resource 1 Eq. 11 (exponential model, retained as
    # the final covariate model over the linear and power alternatives),
    # Eq. 13 (individual value via an exponential IIV) and Eq. 14 (the
    # f / (1 + f) rescaling that keeps the fraction inside (0, 1)).
    #
    # Note the sign: the main-text equation in Section 3.2 prints
    # exp(theta_B * (AS_ref - AS_i)), which would make the fraction FALL
    # as CYP2D6 activity rises. Online Resource 1 Eq. 11 prints
    # exp(theta_eff * (AS_i - AS_ref)), which is the direction the
    # paper's own results require (Table 3 morphine AUC/F rises
    # monotonically with activity score). The supplement form is used
    # here; see vignette Assumptions and deviations.
    # ------------------------------------------------------------
    fm_morphine_raw <- exp(lfm_morphine + etalfm_morphine + e_cyp2d6_fm_morphine * (CYP2D6 - as_ref))
    fm_morphine     <- fm_morphine_raw / (1 + fm_morphine_raw)
    fm_m3g          <- exp(lfm_m3g)

    # ------------------------------------------------------------
    # Micro-constants. Every state is an amount in codeine mass
    # equivalents; Fig. 1 applies no metabolite-to-parent molecular
    # weight ratio to any transfer, so the metabolite volumes are
    # apparent values that absorb the molecular-weight conversion.
    # ------------------------------------------------------------
    kel          <- cl / vc                    # codeine total elimination
    kel_morphine <- cl_morphine / vc_morphine  # morphine total elimination
    kel_c6g      <- cl_c6g / vc_c6g
    kel_m3g      <- cl_m3g / vc_m3g

    # ------------------------------------------------------------
    # ODE system, paper Fig. 1. Codeine elimination splits into the
    # CYP2D6 branch (ke * fm_morphine -> morphine) and the non-CYP2D6
    # glucuronidation branch (ke * (1 - fm_morphine) -> C6G). Morphine
    # elimination splits into the M3G branch (ke_morphine * fm_m3g) and other
    # pathways (ke_morphine * (1 - fm_m3g)), which leave the system.
    # ------------------------------------------------------------
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    d/dt(central_morphine) <-
      kel * fm_morphine * central - kel_morphine * central_morphine

    d/dt(central_c6g) <-
      kel * (1 - fm_morphine) * central - kel_c6g * central_c6g

    d/dt(central_m3g) <-
      kel_morphine * fm_m3g * central_morphine - kel_m3g * central_m3g

    f(depot) <- fdepot

    # ------------------------------------------------------------
    # Observations. Concentrations are in mg/L (equivalently ug/mL);
    # multiply by 1000 to compare with the paper's ng/mL figures.
    # ------------------------------------------------------------
    Cc          <- central          / vc
    Cc_morphine <- central_morphine / vc_morphine
    Cc_c6g      <- central_c6g      / vc_c6g
    Cc_m3g      <- central_m3g      / vc_m3g

    Cc          ~ lnorm(expSd)
    Cc_morphine ~ lnorm(expSd_morphine)
    Cc_c6g      ~ lnorm(expSd_c6g)
    Cc_m3g      ~ lnorm(expSd_m3g)
  })
}
