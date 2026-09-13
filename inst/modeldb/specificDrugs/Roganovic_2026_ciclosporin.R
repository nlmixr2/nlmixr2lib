Roganovic_2026_ciclosporin <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption and",
    "elimination for ciclosporin (cyclosporin A) in paediatric and young-adult",
    "renal transplant recipients followed by therapeutic drug monitoring.",
    "Allometric body weight (centred at 40 kg) acts on CL/F with an estimated",
    "exponent of 0.89 and on V/F with an exponent fixed at 1; haemoglobin",
    "enters CL/F as a linear effect centred at 120 g/L, so CL/F falls as",
    "haemoglobin rises. Correlated interindividual variability on CL/F and",
    "V/F, plus interoccasion variability on CL/F where an occasion is one",
    "therapeutic-drug-monitoring day."
  )
  reference <- paste(
    "Roganovic M, Cvetkovic M, Gojkovic I, Spasojevic B, Jovanovic M,",
    "Miljkovic B, Vucicevic K. Population Pharmacokinetics Model of",
    "Cyclosporin A in Children and Young Adult Renal Transplant Patients:",
    "Focus on Haemoglobin Contribution to Exposure Variability.",
    "Pharmaceutics. 2026;18(1):99. doi:10.3390/pharmaceutics18010099"
  )
  vignette <- "Roganovic_2026_ciclosporin"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling on CL/F and V/F, centred on the rounded cohort",
        "median of 40 kg (Roganovic 2026 Section 3.3; cohort median 39.65 kg,",
        "range 9.8-103 kg, Table 1). The CL/F exponent was estimated (0.89)",
        "rather than fixed at 0.75; the V/F exponent was fixed at 1.",
        "Treated as a time-fixed baseline weight in this implementation; the",
        "source is a retrospective chart review over roughly one year, so",
        "growth within a paediatric subject was present in the fitted data."
      ),
      source_name        = "WT"
    ),
    HGB = list(
      description        = "Blood haemoglobin concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Linear (not power) effect on CL/F centred at 120 g/L, the cohort",
        "median (Roganovic 2026 Table 1; range 73-164 g/L). Roganovic 2026",
        "Eq. 4 writes the effect as (1 - 0.00279 * (HGB - 120)), so CL/F",
        "DECREASES as haemoglobin rises. The paper attributes this to",
        "erythrocyte binding: ciclosporin is measured in whole blood and is",
        "extensively bound to red cells, so a lower haemoglobin leaves less",
        "drug sequestered and raises apparent clearance. Time-varying in the",
        "source data (repeat laboratory values across TDM occasions). The",
        "linear form goes negative for HGB above about 478 g/L, far outside",
        "any physiological range, but users extrapolating beyond the fitted",
        "73-164 g/L window should clamp the covariate."
      ),
      source_name        = "HGB"
    ),
    OCC = list(
      description        = "Occasion index for interoccasion variability on CL/F",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Roganovic 2026 Section 2.3 defines an occasion as one",
        "therapeutic-drug-monitoring day: 'Every pair of concentrations",
        "measured on the same day, or a single concentration measured on a",
        "given day, was treated as a separate occasion.' The fitted dataset",
        "reached a maximum of 33 occasions per patient over one year",
        "(Section 3.3), and the paper's own simulations (Table 4) use two",
        "occasions per patient. rxode2 has no native occasion level, so the",
        "single shared IOV variance is expanded here into six",
        "indicator-multiplexed etas -- one per TDM day over a one-week",
        "monitoring window, a superset of the paper's two-occasion",
        "simulation. Occasion 1 carries the estimated variance and occasions",
        "2-6 are fixed to the same value, reproducing NONMEM",
        "$OMEGA BLOCK(1) SAME. The occasion count is an extraction-side",
        "construct; adding further occasions is mechanical (one more",
        "etaiov_cl_k line plus one more term in iov_cl). Records outside",
        "1-6 carry no IOV."
      ),
      source_name        = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at transplantation",
      units       = "years",
      type        = "continuous",
      notes       = "Screened in the stepwise covariate modelling procedure (Roganovic 2026 Section 2.3) but not retained in the final model. No maturation function was applied because only 3 of 58 patients were aged 2 years or younger and ciclosporin metabolism is near-adult by age 2 (Discussion)."
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened in the stepwise covariate modelling procedure but not retained. Height is highly correlated with body weight (Roganovic 2026 Figure 2) and already enters the Schwartz creatinine-clearance estimate."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened in the stepwise covariate modelling procedure but not retained."
    ),
    HCT = list(
      description = "Haematocrit",
      units       = "percent",
      type        = "continuous",
      notes       = "Explicitly tested and NOT selected in either the forward-inclusion or backward-elimination step (Roganovic 2026 Section 3.3). Notable because haematocrit is the erythrocyte marker most other ciclosporin population models retain; the paper's central claim is that haemoglobin is the more sensitive marker in this cohort."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Selected by the stepwise procedure but dropped from the final model: the effect was statistically significant yet changed CL/F by less than 20 percent across the usual covariate range, the paper's pre-specified clinical-relevance threshold, and its relative standard error was 58.8 percent (Roganovic 2026 Section 3.3)."
    ),
    CRCL = list(
      description = "Creatinine clearance estimated with the revised (bedside) Schwartz formula, CRCL = 0.413 * HT / CREAT",
      units       = "mL/min/1.73m^2",
      type        = "continuous",
      notes       = "Selected by the stepwise procedure but dropped: in the model carrying creatinine clearance on CL/F the relative standard error of the effect was 572.7 percent (Roganovic 2026 Section 3.3). Highly correlated with serum creatinine (Figure 2), so the two were tested separately."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as the categorical covariate GEND but not retained. The cohort was 34 of 58 male (58.62 percent), so 41.38 percent female."
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "ciclosporin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ciclosporin", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species          = "human",
    n_subjects       = 58,
    n_studies        = 1,
    n_observations   = 974,
    age_range        = "1-25 years (at transplantation)",
    age_median       = "12 years",
    weight_range     = "9.8-103 kg",
    weight_median    = "39.65 kg",
    sex_female_pct   = 41.38,
    disease_state    = "renal transplant recipients on ciclosporin, corticosteroid and mycophenolic acid",
    dose_range       = "5 mg/kg/day orally at initiation, divided into two or three daily doses, then adjusted by therapeutic drug monitoring",
    regions          = "Serbia (single centre: University Children's Hospital, Belgrade)",
    renal_function   = "creatinine clearance 8.93-131.58 mL/min/1.73m^2 (median 55.21) by the revised Schwartz formula",
    notes            = paste(
      "Retrospective chart review; 47 of 58 patients (81.03 percent) were",
      "children under 18 and 11 (18.97 percent) were young adults aged 18-25.",
      "30 patients received a living-donor and 28 a cadaveric graft. Baseline",
      "characteristics are Roganovic 2026 Table 1. Whole-blood ciclosporin was",
      "measured by Abbott chemiluminescent microparticle immunoassay, lower",
      "limit of quantification 30 ng/mL and upper limit 1500 ng/mL. The 974",
      "steady-state samples were 471 pre-dose troughs (C0), 501 two-hour",
      "post-dose samples (C2), one 4 h and one 6 h sample; samples were drawn",
      "after at least three consecutive days on an unchanged regimen."
    )
  )

  ini({
    # Structural parameters, referenced to a 40 kg patient with haemoglobin 120 g/L
    lka <- fixed(log(1.15))
    label("Absorption rate constant (1/h)")  # Section 3.3 and Discussion: Ka fixed at 1.15 1/h because sparse absorption-phase sampling made it unidentifiable; a sensitivity analysis over a plausible Ka range changed the other estimates minimally

    lcl <- log(15)
    label("Apparent clearance CL/F for a 40 kg patient at haemoglobin 120 g/L (L/h)")  # Table 3: CL/F = 15 L/h/40 kg (RSE 4.7%; bootstrap median 14.92, 95% CI 13.63-16.35); Eq. 4

    lvc <- log(71.1)
    label("Apparent central volume V/F for a 40 kg patient (L)")  # Table 3: V/F = 71.1 L/40 kg (RSE 5.8%; bootstrap median 70.76, 95% CI 63.034-79.71); Eq. 5

    # Allometric body-weight exponents, weight centred at the rounded cohort median of 40 kg
    e_wt_cl <- 0.89
    label("Allometric exponent of (WT/40 kg) on CL/F (unitless)")  # Table 3: theta_ALL = 0.89 (RSE 5.7%; bootstrap median 0.89, 95% CI 0.78-0.98); Eq. 4. Estimated rather than fixed at 0.75 because the cohort spans children and young adults

    e_wt_vc <- fixed(1)
    label("Allometric exponent of (WT/40 kg) on V/F (unitless)")  # Section 3.3: 'the exponent for V/F was kept fixed at 1'; Eq. 5 prints the exponent as 1

    # Covariate effect
    e_hgb_cl <- -0.00279
    label("Linear haemoglobin effect on CL/F, centred at 120 g/L (per g/L)")  # Table 3: theta_HGB = -0.00279 (RSE 24%; bootstrap median -0.0028, 95% CI -0.0041 to -0.0016); Eq. 4 applies it as (1 - 0.00279 * (HGB - 120))

    # Correlated interindividual variability on CL/F and V/F.
    # Table 3 reports IIV on the SD scale as percentages: IIV CL = 34.91% and
    # IIV V = 43.05%, so the variances are 0.3491^2 and 0.4305^2. The
    # off-diagonal is Table 3 row 'omega (CL-V)' = 0.136, which the table
    # footnote defines as the covariance between CL and V. That implies a
    # correlation of 0.136 / (0.3491 * 0.4305) = 0.905 -- high, but expected
    # for apparent parameters, since CL/F and V/F both carry the same
    # unmeasured 1/F factor in an oral sparse-sampling design.
    etalcl + etalvc ~ c(0.3491^2, 0.136, 0.4305^2)  # Table 3 rows 'IIV CL (%)' 34.91 (RSE 15.5), 'omega (CL-V)' 0.136 (RSE 27.3), 'IIV V (%)' 43.05 (RSE 11.8)

    # Interoccasion variability on CL/F, one shared variance across occasions
    # (NONMEM $OMEGA BLOCK(1) SAME). Table 3 row 'IOV CL (%)' = 12.25 on the
    # SD scale, so the variance is 0.1225^2.
    etaiov_cl_1 ~ 0.1225^2         # Table 3 row 'IOV CL (%)' = 12.25 (RSE 12.5; bootstrap median 12.14, 95% CI 8.7-15.87)
    etaiov_cl_2 ~ fixed(0.1225^2)  # shared IOV variance, occasion 2
    etaiov_cl_3 ~ fixed(0.1225^2)  # shared IOV variance, occasion 3
    etaiov_cl_4 ~ fixed(0.1225^2)  # shared IOV variance, occasion 4
    etaiov_cl_5 ~ fixed(0.1225^2)  # shared IOV variance, occasion 5
    etaiov_cl_6 ~ fixed(0.1225^2)  # shared IOV variance, occasion 6

    # Residual error
    propSd <- 0.258
    label("Proportional residual error (fraction)")  # Table 3 row 'Wp (%)' = 0.258 (RSE 4.3; bootstrap median 0.26, 95% CI 0.23-0.28). Read on the SD scale, i.e. 25.8% proportional error: the alternative variance reading would give 50.8% CV, which exceeds the total observed spread of the C2 samples beyond day 90 (mean 716.69, SD 243.68 ng/mL, i.e. 34% CV, Table 2) and so is arithmetically impossible for a residual term
  })

  model({
    # Occasion indicators. Roganovic 2026 defines an occasion as one TDM day;
    # six are carried here (see covariateData$OCC). Records with OCC outside
    # 1-6 zero every indicator and therefore carry no IOV.
    occ1 <- (OCC == 1)
    occ2 <- (OCC == 2)
    occ3 <- (OCC == 3)
    occ4 <- (OCC == 4)
    occ5 <- (OCC == 5)
    occ6 <- (OCC == 6)
    iov_cl <-
      occ1 * etaiov_cl_1 + occ2 * etaiov_cl_2 + occ3 * etaiov_cl_3 +
      occ4 * etaiov_cl_4 + occ5 * etaiov_cl_5 + occ6 * etaiov_cl_6

    ka <- exp(lka)

    # Eq. 4: CL/F = 15 * (WT/40)^0.89 * (1 - 0.00279 * (HGB - 120))
    cl <- exp(lcl + etalcl + iov_cl) * (WT / 40)^e_wt_cl *
      (1 + e_hgb_cl * (HGB - 120))

    # Eq. 5: V/F = 71.1 * (WT/40)^1
    vc <- exp(lvc + etalvc) * (WT / 40)^e_wt_vc

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Doses are in mg and vc is in L, so central/vc is mg/L; the assay and every
    # reported concentration in Roganovic 2026 are in ng/mL (1 mg/L = 1000 ng/mL).
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
