Fokkink_2022_immunoglobulin <- function() {
  description <- "Two-compartment population PK model for high-dose intravenous immunoglobulin (IVIg) in Guillain-Barre syndrome, fitted to total serum IgG (IVIg-derived plus an estimated endogenous baseline) with allometric weight scaling and disease-severity (GBS disability score) and concomitant methylprednisolone effects on clearance (Fokkink 2022)"
  reference <- paste(
    "Fokkink WJR, van Tilburg SJ, de Winter BCM, Sassen SDT, van Doorn PA,",
    "Koch BCP, Jacobs BC. Population pharmacokinetic modelling of intravenous",
    "immunoglobulin treatment in patients with Guillain-Barre syndrome.",
    "Clin Pharmacokinet. 2022;61(9):1285-1296. doi:10.1007/s40262-022-01136-z",
    sep = " "
  )
  vignette <- "Fokkink_2022_immunoglobulin"
  units <- list(time = "day", dosing = "g", concentration = "g/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling (WT/70) with fixed exponents: 0.75 on CL and Q, 1 on V1 and V2 (Results 3.1.1; ESM control stream THETA(6) = 0.75 FIX, THETA(7) = 1 FIX).",
      source_name = "WGHT"
    ),
    SCORE_GBS_DISABILITY = list(
      description = "Guillain-Barre syndrome disability score (Hughes GBS-DS), ordinal 0 (healthy) to 6 (dead)",
      units = "(GBS-DS points, 0-6)",
      type = "continuous",
      reference_category = "0-2 (combined reference group)",
      notes = "Grouped inside model() into the paper's three categories: 0-2 (reference), 3-4 (CL x 1.22) and 5 (mechanical ventilation, CL x 1.83). A score of 6 (death) was disregarded by the authors and falls in neither effect group here. The ESM control stream codes the score per record in column GBS (CIDP patients -99); the paper does not state whether the score was time-varying or held at a single visit, so supply the per-record score as recorded. Results 3.1.2 and Table 2.",
      source_name = "GBS"
    ),
    CONMED_METHYLPREDNISOLONE = list(
      description = "Concomitant intravenous methylprednisolone (500 mg/day for 5 days, started within 48 h of the first IVIg dose)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (IVIg only)",
      notes = "CL x 1.31 when 1 (Table 2; ESM control stream THETA(14)^MP).",
      source_name = "MP"
    ),
    DIS_GBS = list(
      description = "Disease type: 1 = Guillain-Barre syndrome, 0 = chronic inflammatory demyelinating polyneuropathy (CIDP) on IVIg maintenance",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CIDP)",
      notes = "Multiplies the endogenous baseline IgG by 0.65 for GBS (Methods 2.3.1, Table 2 'Disease type (if GBS)'; ESM THETA(9)^TYPE). The CIDP group had a higher steady-state baseline because of prolonged IVIg maintenance. The control stream also estimates a CIDP-specific CL factor (THETA(10)) that the paper does not report; it is not encoded, so CIDP subjects take the reference CL.",
      source_name = "TYPE"
    )
  )

  covariatesDataExcluded <- list(
    VNTR_FCGRT = list(
      description = "FcRn (FCGRT) variable-number-tandem-repeat genotype",
      units = "(categorical)",
      type = "categorical",
      notes = "Screened in the univariate analysis (available in 71 GBS patients) but NOT retained: 'Sex, age and VNTR were not significantly correlated with any of the PK parameters' (Results 3.1.2). No point estimate is reported.",
      source_name = "VNTR"
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Screened but not retained (Results 3.1.2).",
      source_name = "SEX"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened but not retained (Results 3.1.2).",
      source_name = "AGE"
    )
  )

  compartmentData <- list(
    central = list(analyte = "IVIg-derived immunoglobulin G", units = "g", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "IVIg-derived immunoglobulin G", units = "g", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 192L,
    n_studies = 3L,
    age_range = "5-89 years (GBS model-building cohort)",
    age_mean = "53 years",
    weight_range = "18-122 kg",
    weight_mean = "72 kg",
    sex_female_pct = 43,
    race_ethnicity = "Not reported (Dutch multicentre trials)",
    disease_state = "Guillain-Barre syndrome (177 patients, unable to walk independently at entry) plus 15 patients with stable chronic inflammatory demyelinating polyneuropathy (CIDP) on IVIg maintenance",
    dose_range = "IVIg 0.4 g/kg/day for 5 consecutive days (GBS); fixed maintenance regimens (CIDP)",
    regions = "Europe (the Netherlands)",
    notes = "Model-building cohort: 177 GBS patients (589 total-IgG samples) from two Dutch GBS Study Group trials (IVIg vs plasma exchange; IVIg + methylprednisolone vs placebo) plus 15 CIDP patients (222 samples) contributing early post-infusion data. Demographics are for the GBS cohort (Table 1): GBS-DS at nadir mean 4.1 (SD 0.7), 25% mechanically ventilated, 34% received methylprednisolone; male 57%. External validation in 177 SID-GBS good-prognosis patients (689 samples). All model-building GBS patients received Gammagard S/D."
  )

  ini({
    # Structural parameters: Table 2, final model. Time unit day (Table 2
    # reports CL and Q in L/day). Reference weight 70 kg.
    lcl <- log(0.28); label("Clearance for a 70 kg patient with GBS-DS 0-2, no methylprednisolone (L/day)") # Table 2 'CL (L/day/70 kg)' final = 0.28 (RSE 9%)
    lvc <- log(2.87); label("Central volume of distribution for a 70 kg patient (L)") # Table 2 'V1 (L)' final = 2.87 (RSE 7%)
    lq <- log(0.25); label("Intercompartmental clearance for a 70 kg patient (L/day)") # Table 2 'Q (L/day)' final = 0.25 (RSE 17%)
    lvp <- log(2.65); label("Peripheral volume of distribution for a 70 kg patient (L)") # Table 2 'V2 (L)' final = 2.65 (RSE 13%)

    # Allometric exponents fixed (Results 3.1.1; ESM THETA(6), THETA(7) FIX).
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent on CL and Q (unitless)") # Results 3.1.1 and ESM THETA(6) = 0.75 FIX
    e_wt_vc_vp <- fixed(1); label("Allometric exponent on V1 and V2 (unitless)") # Results 3.1.1 and ESM THETA(7) = 1 FIX

    # Endogenous baseline IgG (CBAS): Methods 2.3.1 equation, CBAS = theta_CBAS * theta_T^TYPE.
    lbl_igg <- log(15.90); label("Endogenous baseline serum IgG in CIDP patients (g/L)") # Table 2 'CBAS (g/L)' final = 15.90 (RSE 8%)
    e_dis_gbs_bl_igg <- 0.65; label("Multiplicative factor on baseline IgG for GBS vs CIDP (unitless)") # Table 2 'Disease type (if GBS)' final = 0.65 (RSE 9%)

    # Clearance covariates, multiplicative (Results 3.1.2 final CL equation).
    e_conmed_methylprednisolone_cl <- 1.31; label("Multiplicative factor on CL with concomitant methylprednisolone (unitless)") # Table 2 'Methylprednisolone' final = 1.31 (RSE 9%)
    e_score_gbs_disability_34_cl <- 1.22; label("Multiplicative factor on CL for GBS-DS 3 or 4 vs 0-2 (unitless)") # Table 2 'GBS disability score of 3-4' final = 1.22 (RSE 10%)
    e_score_gbs_disability_5_cl <- 1.83; label("Multiplicative factor on CL for GBS-DS 5 vs 0-2 (unitless)") # Table 2 'GBS disability score of 5' final = 1.83 (RSE 15%)

    # IIV reported as %CV (Table 2 legend lists CV); omega^2 = log(1 + CV^2).
    # The ESM control stream estimates an OMEGA BLOCK(2) whose covariance is
    # not reported, so the two etas are encoded as independent.
    etalcl ~ 0.221429 # Table 2 IIV 'CL (%)' final = 49.9 -> log(1 + 0.499^2)
    etalbl_igg ~ 0.052000 # Table 2 IIV 'CBAS (%)' final = 23.1 -> log(1 + 0.231^2)

    # Residual error: proportional on total IgG (ESM $ERROR: Y = FTOT + FTOT * THETA(5) * EPS(1), SIGMA 1 FIX).
    propSd <- 0.12; label("Proportional residual error (fraction)") # Table 2 'Proportional' final = 0.12 (RSE 7%)
  })
  model({
    # GBS disability score grouping (Results 3.1.2): 0-2 reference, 3-4, 5.
    gbsds_34 <- 0
    if (SCORE_GBS_DISABILITY >= 3 & SCORE_GBS_DISABILITY <= 4) gbsds_34 <- 1
    gbsds_5 <- 0
    if (SCORE_GBS_DISABILITY == 5) gbsds_5 <- 1

    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q *
      e_conmed_methylprednisolone_cl^CONMED_METHYLPREDNISOLONE *
      e_score_gbs_disability_34_cl^gbsds_34 *
      e_score_gbs_disability_5_cl^gbsds_5
    vc <- exp(lvc) * (WT / 70)^e_wt_vc_vp
    q <- exp(lq) * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp) * (WT / 70)^e_wt_vc_vp

    bl_igg <- exp(lbl_igg + etalbl_igg) * e_dis_gbs_bl_igg^DIS_GBS

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # IVIg is infused directly into the central compartment; the states hold
    # IVIg-derived IgG only.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Observation is TOTAL serum IgG = IVIg-derived + endogenous baseline
    # (ESM $ERROR: FTOT = F + CBAS).
    Cc <- central / vc + bl_igg
    Cc ~ prop(propSd)
  })
}
