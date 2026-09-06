Kim_2023_camostat <- function() {
  description <- paste(
    "Joint parent + metabolite population PK model for oral camostat mesylate",
    "in 15 healthy Korean adult men given a single 100, 200 or 300 mg dose",
    "(Kim 2023). Camostat itself is never measured -- it is hydrolysed by",
    "plasma carboxylesterase too rapidly to quantify -- so the model's parent",
    "analyte is the active metabolite GBPA (4-(4-guanidinobenzoyloxy)phenyl",
    "acetic acid, FOY-251), delivered to a one-compartment central",
    "compartment by first-order absorption with a lag time. GBPA elimination",
    "splits into a linear branch that leaves the system, (1 - fm) * CL, and a",
    "Michaelis-Menten branch, fm * Vmax * C / (Km + C), that forms the",
    "inactive metabolite GBA (4-guanidinobenzoic acid); GBA has its own",
    "one-compartment disposition with linear clearance. Between-subject",
    "variability is estimated on Ka, lag time, GBPA volume, GBPA clearance,",
    "GBA clearance and Vmax, with a -0.89 correlation between the GBPA volume",
    "and Vmax random effects, and residual error is proportional on each",
    "analyte. No covariate effects were retained. IMPORTANT: Vmax is NOT the",
    "value printed in the paper's Table S1 -- that row carries no unit and no",
    "unit assignment reproduces the authors' own published output. The value",
    "here was derived from the paper's Table 2 and Figure 6; see the lvmax",
    "comment and the vignette Errata."
  )
  reference <- paste(
    "Kim G, Moon H-k, Kim T, Yun S-h, Yun H-y, Hong JH, Kim D-D.",
    "Safety Evaluation and Population Pharmacokinetics of Camostat Mesylate",
    "and Its Major Metabolites Using a Phase I Study.",
    "Pharmaceutics. 2023;15(9):2357. doi:10.3390/pharmaceutics15092357.",
    "All parameter estimates are in Supplementary Table S1; no estimate",
    "appears in the main text.",
    sep = " "
  )
  vignette <- "Kim_2023_camostat"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Doses are camostat mesylate in mg. Figure 4 draws no
  # molecular-weight ratio on either metabolic transfer, so every state is an
  # amount in camostat mass equivalents and each volume and clearance is an
  # apparent value that absorbs bioavailability and the molecular-weight
  # ratios. Concentrations nonetheless come out on the measured (assay) scale
  # because the volumes and clearances were estimated against the measured
  # GBPA and GBA plasma concentrations.
  compartmentData <- list(
    depot = list(
      analyte = "camostat mesylate", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "GBPA (FOY-251)", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    central_gba = list(
      analyte = "GBA (4-guanidinobenzoic acid)", units = "mg",
      specimen = "plasma", verified = TRUE
    )
  )

  # No covariate was retained. Paper Section 3.4: 'There are no covariate
  # effects in the PK parameters (height, weight, etc.)'.
  covariateData <- list()

  # Screened during stepwise covariate model building (paper Section 2.6) and
  # not retained (Section 3.4). Recorded here so the covariate screen is
  # preserved without declaring entries that model() never references. The
  # baseline values are paper Table 1.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Table 1 means (SD) by dose group: 70.1 (10.9), 68.5 (10.9) and",
        "74.5 (16.0) kg; overall median 66.5 kg. Enrolment required 55.0-90.0",
        "kg. Screened and not retained; the model carries no allometric",
        "scaling."
      )
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = paste(
        "Table 1 means (SD) by dose group: 173.1 (7.19), 172.8 (5.53) and",
        "175.4 (6.97) cm; overall mean 174.2 cm. Screened and not retained."
      )
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = paste(
        "Table 1 means (SD) by dose group: 23.3 (2.60), 22.9 (2.93) and",
        "24.0 (3.96) kg/m^2. Enrolment required 18.0-29.9 kg/m^2. Reported",
        "as a baseline demographic; not retained."
      )
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Table 1 means (SD) by dose group: 26.8 (6.34), 29.0 (5.56) and",
        "28.0 (5.34) years; overall median 26 years. Enrolment required",
        "19-55 years. Screened and not retained."
      )
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = paste(
        "Table 1 means (SD) by dose group: 17.4 (4.56), 16.6 (4.72) and",
        "13.6 (6.27) IU/L. An AST above the upper limit of normal was an",
        "exclusion criterion. Screened and not retained."
      )
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = paste(
        "Table 1 means (SD) by dose group: 24.2 (8.76), 24.8 (12.5) and",
        "21.8 (12.3) IU/L. An ALT above the upper limit of normal was an",
        "exclusion criterion. Screened and not retained."
      )
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (BSA-normalised)",
      units       = "mL/min/1.73 m^2",
      type        = "continuous",
      notes       = paste(
        "Table 1 means (SD) by dose group: 99.2 (10.7), 113.8 (31.8) and",
        "93.8 (15.3) mL/min/1.73 m^2. An eGFR below 60 was an exclusion",
        "criterion. Screened and not retained."
      )
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = paste(
        "Table 1 means (SD) by dose group: 4.42 (0.21), 4.50 (0.12) and",
        "4.58 (0.29) g/dL, i.e. 44.2, 45.0 and 45.8 g/L. Screened and not",
        "retained."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 15,
    n_studies      = 1,
    n_observations = 360,
    age_range      = "19-55 years by protocol; group means 26.8-29.0, overall median 26 years",
    weight_range   = "55.0-90.0 kg by protocol; group means 68.5-74.5 kg, overall median 66.5 kg",
    height_mean    = "174.2 cm",
    sex_female_pct = 0,
    race_ethnicity = c(Asian = 100),
    disease_state  = "healthy volunteers",
    dose_range     = "single oral 100, 200 or 300 mg camostat mesylate tablet under fasting conditions (n = 5 per dose)",
    regions        = "Republic of Korea (Chungnam National University Hospital)",
    renal_function = "eGFR at or above 60 mL/min/1.73 m^2 required for enrolment",
    notes          = paste(
      "Parallel-group, open-label single-dose Phase 1 study",
      "(ClinicalTrials.gov NCT04782505; IRB CNUH2021-02-018-016). Plasma",
      "sampled at 0, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 5 and 6 h after",
      "dosing; 180 samples assayed for each of GBPA and GBA, giving the 360",
      "observations used in the joint fit. Assay was LC-MS/MS with a",
      "0.5-1000 ng/mL calibration range for both analytes. Estimation was",
      "SAEM in Monolix 2021R2 with a bootstrap for parameter precision.",
      "All subjects were male, so the model carries no information about",
      "female pharmacokinetics; the authors state the simulations should be",
      "used for predictions in males (paper Discussion)."
    )
  )

  ini({
    # ============================================================
    # Absorption. Supplementary Table S1, 'Fixed effect' block.
    # ============================================================
    lka <- log(0.757)
    label("GBPA absorption rate constant (1/h)")
    # Table S1: Ka = 0.757 1/h, RSE 3.71%, bootstrap 0.754 (5-95% 0.697-0.833).

    ltlag <- log(0.204)
    label("Absorption lag time (h)")
    # Table S1: Tlag = 0.204 h, RSE 2.75%, bootstrap 0.203 (5-95% 0.189-0.213).

    # ============================================================
    # GBPA (parent analyte) disposition. Table S1.
    # ============================================================
    lvc <- log(237.9)
    label("Apparent central volume of distribution of GBPA (L)")
    # Table S1 Vd1 = 237.9 L, RSE 22.8%, bootstrap 248.6 (5-95% 175.7-341.0).

    lcl <- log(1269)
    label("Apparent clearance term for GBPA, split between the linear and metabolic branches by fm (L/h)")
    # Table S1 CL1 = 1269, RSE 19.7%, bootstrap 1685 (5-95% 1112-2273). The
    # unit column of Table S1 is blank for this row; L/h is the only reading
    # that is dimensionally admissible for a clearance multiplying a plasma
    # concentration in Figure 4, and 0.184 * 1269 = 233 L/h is a plausible
    # non-metabolic apparent clearance beside the 620-640 L/h total apparent
    # clearance implied by Table 2 (Dose / AUCinf). Retained as printed.

    # ============================================================
    # GBA (metabolite) disposition. Table S1.
    # ============================================================
    lvc_gba <- log(174.6)
    label("Apparent central volume of distribution of GBA (L)")
    # Table S1 Vd2 = 174.6 L, RSE 3.58%, bootstrap 160.0 (5-95% 122.9-213.9).

    lcl_gba <- log(103.7)
    label("Apparent clearance of GBA (L/h)")
    # Table S1 CL2 = 103.7 L/h, RSE 7.41%, bootstrap 94.95 (5-95% 71.13-130.0).
    # Figure 4 names this arm CLmet.

    # ============================================================
    # Metabolic split and Michaelis-Menten conversion of GBPA to GBA.
    # ============================================================
    fm <- 0.816
    label("Fraction of GBPA elimination routed to GBA formation (unitless)")
    # Table S1 estimates the complement: FRAC = 0.184, RSE 0.17%, bootstrap
    # 0.168 (5-95% 0.117-0.199), described as the fraction of cleared as
    # GBPA. Figure 4 multiplies the linear branch by FRAC and the
    # Michaelis-Menten branch by (1 - FRAC), so fm = 1 - 0.184 = 0.816. The
    # Discussion quotes the same number back as about 82% of the GBPA
    # metabolized to GBA in our model.

    lvmax <- fixed(log(650))
    label("Maximum rate of GBPA to GBA metabolic conversion (mg/h; not paper-derived, see comment)")
    # NOT the printed estimate. Table S1 prints Vmax = 2.215 (RSE 23%,
    # bootstrap 1.97, 5-95% 1.334-2.723) with the unit column left blank, and
    # no unit assignment reconciles it: encoded as 2.215 mg/h the metabolic
    # branch carries almost no flux and simulated GBA AUCinf at 100 mg is
    # 5.8 h*ng/mL against the 762.8 in the authors' own Table 2, i.e. under
    # 1%. Neither x60 (1/min to 1/h) nor x1000 (g to mg, or mg to ug) closes
    # the gap. 650 mg/h is the value that minimises the sum of squared log
    # ratios between this model and the paper's own published exposures --
    # Table 2 AUCinf for both analytes at 100, 200 and 300 mg -- and it also
    # reproduces the Figure 6 steady-state peaks. Residual deviations are
    # -4.5%, -3.4% and -1.5% on GBPA AUCinf and -17.7%, -2.3% and +7.3% on
    # GBA AUCinf. Operator-ratified 2026-09-05; full derivation in the
    # vignette Errata.

    lkm <- log(1192)
    label("Michaelis constant for GBPA to GBA conversion (ng/mL)")
    # Table S1 Km = 1.192, RSE 0.19%, bootstrap 1.89, unit column blank.
    # Read as 1.192 mg/L = 1192 ng/mL. mg/L is the only admissible reading:
    # observed GBPA concentrations are 72-274 ng/mL (Table 2 Cmax), so
    # Km = 1192 ng/mL leaves the enzyme only 5-15% saturated at the peak,
    # which is what the paper's own dose-proportionality result requires
    # (log-transformed power-model slopes 1.0038-1.0066 with 95% CIs inside
    # 0.8-1.25 across 100-300 mg). Reading Km as 1.192 ng/mL would leave the
    # enzyme deeply saturated and the model grossly dose-nonlinear. The mild
    # saturation this reading does produce has the sign the paper's own
    # dose-normalised exposures show: GBPA AUCinf per 100 mg rises across
    # doses (156.5, 158.5, 159.2) while GBA AUCinf per 100 mg falls.

    # ============================================================
    # Between-subject variability. Table S1 'Between-subject Variability'
    # block. Monolix reports omega as the standard deviation of the random
    # effect, so the variances below are the squares of the tabulated
    # values. Note a text/table conflict: paper Section 3.4 lists BSV on
    # 'lag time, Vmax, Km, V1, clearance of GBPA, and clearance of GBA' --
    # omitting Ka and including Km -- while Table S1 lists Ka and has no Km
    # row. The table carries the numbers and is followed here.
    # ============================================================
    etalka ~ 0.011881
    # Table S1 BSV Ka = 0.109 (RSE 31.0%, bootstrap 0.095, 5-95%
    # 0.053-0.134); variance = 0.109^2 = 0.011881.

    etaltlag ~ 0.007921
    # Table S1 BSV Tlag = 0.089 (RSE 26.2%, bootstrap 0.093, 5-95%
    # 0.047-0.158); variance = 0.089^2 = 0.007921.

    etalcl ~ 0.395641
    # Table S1 BSV CL1 = 0.629 (RSE 26.5%, bootstrap 0.489, 5-95%
    # 0.286-0.714); variance = 0.629^2 = 0.395641.

    etalcl_gba ~ 0.044100
    # Table S1 BSV CL2 = 0.210 (RSE 26.4%, bootstrap 0.180, 5-95%
    # 0.051-0.281); variance = 0.210^2 = 0.044100.

    etalvc + etalvmax ~ c(0.660969, -0.593327, 0.672400)
    # Table S1 BSV V1 = 0.813 (RSE 24.1%) and BSV Vmax = 0.820 (RSE 19.7%),
    # with Corr_Vmax_V1 = -0.89 (RSE 8.62%, bootstrap -0.868).
    # var_vc = 0.813^2 = 0.660969; var_vmax = 0.820^2 = 0.672400;
    # cov = -0.89 * 0.813 * 0.820 = -0.593327.

    # ============================================================
    # Residual variability. Table S1 'Error model' block; the paper
    # describes a proportional error model for both analytes
    # (Section 3.4 and Figure 4 discussion).
    # ============================================================
    propSd <- 0.25
    label("Proportional residual error SD for GBPA (fraction)")
    # Table S1 Pro_GBPA = 0.25, RSE 6.66%, bootstrap 0.247 (5-95%
    # 0.209-0.288).

    propSd_gba <- 0.23
    label("Proportional residual error SD for GBA (fraction)")
    # Table S1 Pro_GBA = 0.23, RSE 6.74%, bootstrap 0.227 (5-95%
    # 0.171-0.283).
  })

  model({
    # ------------------------------------------------------------
    # Individual parameters. No covariate effects were retained.
    # ------------------------------------------------------------
    ka     <- exp(lka + etalka)
    tlag   <- exp(ltlag + etaltlag)
    vc     <- exp(lvc + etalvc)
    cl     <- exp(lcl + etalcl)
    vc_gba <- exp(lvc_gba)
    cl_gba <- exp(lcl_gba + etalcl_gba)
    vmax   <- exp(lvmax + etalvmax)
    km     <- exp(lkm)

    # ------------------------------------------------------------
    # Observables. Amounts are in mg and volumes in L, so amount / volume
    # is mg/L; the factor 1000 puts both concentrations on the ng/mL scale
    # the paper reports and on which Km is expressed.
    # ------------------------------------------------------------
    Cc     <- central / vc * 1000
    Cc_gba <- central_gba / vc_gba * 1000

    # ------------------------------------------------------------
    # Metabolic formation rate, Figure 4:
    #   (1 - FRAC) * Vmax * C_GBPA / (Km + C_GBPA), with fm = 1 - FRAC.
    # The concentration ratio is dimensionless, so the rate is in mg/h.
    # ------------------------------------------------------------
    rmet <- fm * vmax * Cc / (km + Cc)

    # ------------------------------------------------------------
    # ODE system, Figure 4. GBPA leaves the central compartment by a
    # linear branch, (FRAC) * CL, and by the Michaelis-Menten branch that
    # forms GBA; GBA leaves by the linear CLmet.
    # ------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - (1 - fm) * cl * central / vc - rmet
    d/dt(central_gba) <-  rmet - cl_gba * central_gba / vc_gba

    alag(depot) <- tlag

    Cc     ~ prop(propSd)
    Cc_gba ~ prop(propSd_gba)
  })
}
