Li_2026_timolol <- function() {
  description <- "One-compartment population PK model with first-order absorption and elimination for topical 0.5% timolol maleate gel in healthy adults and infants with proliferative superficial infantile hemangioma, with allometric weight scaling and study-specific relative bioavailability in the two infant trials (Li 2026). Doses are the mass of gel applied (mg of gel), not mg of timolol."
  reference <- paste(
    "Li L, Li Y, Hu F, Wang D, Sun N, Meng Q, Yu Q, Xu L, Zhao L, Ma L.",
    "Population Pharmacokinetics of Topical Timolol Maleate Gel in Healthy",
    "Volunteers and Infants With Superficial Infantile Hemangioma.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(10):e70334.",
    "doi:10.1002/psp4.70334.",
    sep = " "
  )
  vignette <- "Li_2026_timolol"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "timolol", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "timolol", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling of CL/F (exponent 0.75, fixed) and Vc/F (exponent 1, fixed) with reference weight 60 kg (Li 2026 Table 1 footnote: CL = CLTV x (WT/60)^0.75; Vc = Vc,TV x (WT/60)).",
      source_name = "WT"
    ),
    STUDY_CTR20202613 = list(
      description = "1 = record from the multicentre phase Ib/II study of timolol maleate gel in infants with proliferative superficial infantile hemangioma (CTR20202613); 0 = any other study in the pooled analysis",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (phase I healthy-adult study CTR20191402, F1 = 1)",
      notes = "Multiplies the relative bioavailability by 252 (Li 2026 Table 1 'F1Ib/II'). Set STUDY_CTR20202613 = 0 and STUDY_CTR20221992 = 0 for healthy adults; set exactly one of them to 1 for infants with infantile hemangioma.",
      source_name = "F1Ib/II (study stratum)"
    ),
    STUDY_CTR20221992 = list(
      description = "1 = record from the randomized, double-blind, placebo-controlled phase III study of timolol maleate gel in infants with proliferative superficial infantile hemangioma (CTR20221992); 0 = any other study in the pooled analysis",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (phase I healthy-adult study CTR20191402, F1 = 1)",
      notes = "Multiplies the relative bioavailability by 202 (Li 2026 Table 1 'F1III'). This is the three-times-daily commercial regimen and the default infant stratum for a new simulation.",
      source_name = "F1III (study stratum)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years (adults) / months (infants)",
      type = "continuous",
      notes = "Screened against the CL/F random effect in infants (Li 2026 Figure S2); not retained (p > 0.01)."
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Screened against the CL/F random effect in infants (Li 2026 Figure S2); not retained."
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened (Li 2026 Figure S2); not retained."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "IU/L",
      type = "continuous",
      notes = "Screened (Li 2026 Figure S2); not retained."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "IU/L",
      type = "continuous",
      notes = "Screened (Li 2026 Figure S2); not retained."
    ),
    TBIL = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Screened (Li 2026 Figure S2); not retained."
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units = "mL/min",
      type = "continuous",
      notes = "Screened (Li 2026 Figure S2); not retained."
    ),
    BHW = list(
      description = "Baseline hemangioma thickness",
      units = "mm",
      type = "continuous",
      notes = "Screened (Li 2026 Figure S2, labelled BHT there); not retained."
    ),
    BHS = list(
      description = "Baseline hemangioma surface area",
      units = "cm^2",
      type = "continuous",
      notes = "Screened (Li 2026 Figure S2); not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 138L,
    n_studies = 3L,
    n_observations = 810L,
    age_range = "Healthy adults 21-38 years; infants 1.2-6.0 months",
    age_median = "Healthy adults 29.0 years; infants 2.40 months",
    weight_range = "Healthy adults 49.4-76.2 kg; infants 3.50-9.57 kg",
    weight_median = "Healthy adults 58.4 kg; infants 6.10 kg",
    sex_female_pct = 56.5,
    race_ethnicity = "Not reported (Chinese trials)",
    disease_state = "24 healthy adult volunteers (phase I) and 114 infants with proliferative superficial infantile hemangioma (35 phase Ib/II, 79 phase III).",
    dose_range = "0.5% timolol maleate gel at a nominal ~3 mg of gel per cm^2 per application. Adults: 10% of body surface area (~4.8 g gel per application), single dose then 2, 4 or 6 times daily for 6 days. Infants: applied to the lesion (~9 mg gel for a 3 cm^2 lesion), once to five times daily for 7 days in phase Ib/II and three times daily for 24 weeks in phase III.",
    regions = "China",
    notes = "Li 2026 Results 'Demographics and Data Summary' and Tables S3-S4. 478 concentrations from adults and 332 from infants, measured in dried blood spots by LC-MS/MS (LLOQ 0.0500 ng/mL). All infants were older than 90 days at first dose. The paper does not report the infant dose amounts in its analysis dataset; with the nominal 3 mg of gel per cm^2 of lesion, the model predicts infant concentrations one to two orders of magnitude below the published infant data while reproducing the adult data (see the vignette)."
  )

  ini({
    lka <- log(0.01); label("First-order absorption rate constant (1/h)") # Table 1 'Ka, h-1' = 0.01 (RSE 4.3%)
    lcl <- log(3030e3); label("Apparent clearance CL/F for a 60 kg adult, per mg of gel applied (L/h)") # Table 1 'CL/F, x10^3 L/h' = 3030 (RSE 39.4%)
    lvc <- log(22100e3); label("Apparent central volume Vc/F for a 60 kg adult, per mg of gel applied (L)") # Table 1 'Vc/F, x10^3 L' = 22,100 (RSE 21.3%)
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL/F (unitless)") # Table 1 'WT on CL/F' = 0.75 FIX
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on Vc/F (unitless)") # Table 1 'WT on Vc/F' = 1 FIX
    lfdepot <- fixed(log(1)); label("Relative bioavailability in healthy adults (unitless)") # Table 1 footnote: F1 in healthy adults fixed to 1
    e_study_ctr20202613_fdepot <- 252; label("Relative bioavailability in the phase Ib/II infant study vs healthy adults (ratio)") # Table 1 'F1Ib/II' = 252 (RSE 29.4%)
    e_study_ctr20221992_fdepot <- 202; label("Relative bioavailability in the phase III infant study vs healthy adults (ratio)") # Table 1 'F1III' = 202 (RSE 27.9%)

    etalcl ~ 1.948816 # Table 1 'omega (CL), %' = 139.6 read as the log-scale SD 1.396; variance = 1.396^2

    propSd <- 0.49; label("Proportional residual error (fraction)") # Table 1 'sigma (Prop), %' = 49.0 (RSE 6.0%)
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 60)^e_wt_cl
    vc <- exp(lvc) * (WT / 60)^e_wt_vc
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    f(depot) <- exp(lfdepot) *
      e_study_ctr20202613_fdepot^STUDY_CTR20202613 *
      e_study_ctr20221992_fdepot^STUDY_CTR20221992

    # Amount in mg over volume in L gives mg/L; x 1000 converts to ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
