Grimsley_1999_vancomycin <- function() {
  description <- "One-compartment IV-infusion population PK model for vancomycin in neonates and young infants (Grimsley 1999). Developed from routine therapeutic-drug-monitoring data in 59 neonates (347 concentrations). Clearance scales linearly with body weight and inversely with serum creatinine concentration (CL = 3.56 * WT / CREAT, L/h, WT in kg, CREAT in umol/L); central volume scales linearly with body weight (V = 0.669 * WT, L/kg). The covariate-coupled CL form (no separately estimated exponents) is reported by the paper as the entire structural model."
  reference <- "Grimsley C, Thomson AH. Pharmacokinetics and dose requirements of vancomycin in neonates. Arch Dis Child Fetal Neonatal Ed. 1999;81(3):F221-F227. doi:10.1136/fn.81.3.f221"
  vignette <- "Grimsley_1999_vancomycin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Grimsley 1999 Table 2: median 1.52 kg (range 0.57-4.23). Time-varying within the routine TDM dataset; the structural CL and V terms apply linearly to WT with no separately estimated exponent (CL = 3.56 * WT / CREAT, V = 0.669 * WT per Table 4 final-model column). The covariate-coupled form was preferred over a postconceptual-age (PCA) covariate because PCA and WT were highly correlated (r = 0.89) and adding PCA on top of weight + creatinine did not improve OFV.",
      source_name        = "weight"
    ),
    CREAT = list(
      description        = "Serum creatinine concentration",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Grimsley 1999 Table 2: median 49 umol/L (range 18-172). Stored under the canonical CREAT column. Enters CL as a simple 1/CREAT term (Table 4 final model: CL = 3.56 * WT / CREAT). The paper compared a power-law CREAT^theta form against the 1/CREAT form and found that fixing the exponent at -1 (i.e., theta1 / CREAT) had a negligible OFV penalty (-4.9) relative to the estimated-exponent form, so the simpler 1/CREAT form was retained. One outlier with rapidly changing creatinine (459 -> 316 umol/L over 3 days) was removed before model building; downstream users should treat CREAT as a slowly-varying covariate within the model's domain of applicability.",
      source_name        = "creatinine"
    )
  )

  covariatesDataExcluded <- list(
    PAGE = list(
      description = "Postconceptual age (PCA) - screened but not retained",
      units       = "weeks",
      type        = "continuous",
      notes       = "Grimsley 1999 Table 3 model 8: PCA-on-CL was tested on top of the weight+creatinine model and showed no significant improvement (delta-OFV = 0.0). Highly correlated with weight (r = 0.89) within the cohort, so the paper omitted it from the final model."
    ),
    GA = list(
      description = "Gestational age at birth - screened but not retained",
      units       = "weeks",
      type        = "continuous",
      notes       = "Grimsley 1999 Table 3 model 11: a categorical cutoff (GA < 35 weeks) on CL achieved borderline significance (delta-OFV = 7.9) but the coefficient suggested an unexpected 27% increase in CL for GA < 35 weeks. Authors judged this to be a spurious result driven by limited GA < 35 numbers and omitted it from the final model."
    ),
    PNA = list(
      description = "Postnatal age - screened but not retained",
      units       = "days",
      type        = "continuous",
      notes       = "Grimsley 1999 Results: tested and produced no significant effect on CL or V."
    ),
    SEXF = list(
      description = "Biological sex (1 = female) - screened but not retained",
      units       = "(binary)",
      type        = "binary",
      notes       = "Grimsley 1999 Results: tested and produced no significant effect."
    ),
    NUT_IV = list(
      description = "Intravenous nutrition indicator - screened but not retained",
      units       = "(binary)",
      type        = "binary",
      notes       = "Grimsley 1999 Table 3 model 10: a multiplicative IV-nutrition effect on CL was tested and achieved borderline significance (delta-OFV = 8.2, ~21% higher CL with IV feeding) but had little influence on scatterplots and IIV. Authors concluded the effect was likely confounded with weight and omitted it from the final model."
    ),
    DOPAMINE = list(
      description = "Concurrent dopamine use - screened but not retained",
      units       = "(binary)",
      type        = "binary",
      notes       = "Grimsley 1999 Results: only 2 of 59 patients received dopamine; no significant effect detected (the Seay 1994 popPK in a larger cohort did retain a dopamine multiplicative effect on CL = 0.455 per administration, but Grimsley 1999's smaller dopamine subgroup could not confirm)."
    ),
    APGAR5 = list(
      description = "5-minute Apgar score - screened but not retained",
      units       = "score (0-10)",
      type        = "continuous",
      notes       = "Grimsley 1999 Results: tested and produced no significant effect."
    ),
    CARDIAC = list(
      description = "Cardiac defect indicator - screened but not retained",
      units       = "(binary)",
      type        = "binary",
      notes       = "Grimsley 1999 Results: tested and produced no significant effect."
    ),
    INFECTED = list(
      description = "Confirmed infection indicator - screened but not retained",
      units       = "(binary)",
      type        = "binary",
      notes       = "Grimsley 1999 Results: tested and produced no significant effect."
    ),
    VENT_STATUS = list(
      description = "Ventilation status (air / oxygen / active ventilation) - screened but not retained",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Grimsley 1999 Results: tested and produced no significant effect."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 59L,
    n_studies        = 1L,
    age_range        = "postconceptual age 26-45 weeks; postnatal age 2-76 days",
    age_median       = "postconceptual age 32 weeks; postnatal age 19 days",
    weight_range     = "0.57-4.23 kg",
    weight_median    = "1.52 kg",
    sex_female_pct   = 36,
    race_ethnicity   = "Not reported (single-center cohort at Yorkhill NHS Trust, Glasgow)",
    disease_state    = "Neonates and young infants under 3 months of age requiring intensive care and receiving vancomycin during routine therapeutic drug monitoring; 85% had confirmed infection, 20% had a cardiac defect, 61% required oxygen or active ventilation",
    dose_range       = "Intravenous vancomycin administered as 1 h infusions, starting doses 15 mg/kg q24h to 15 mg/kg q8h per the Guy's, St Thomas's and Lewisham Hospitals Paediatric Formulary; 64 courses of treatment in 59 patients",
    regions          = "United Kingdom (Yorkhill Hospitals neonatal intensive care, Glasgow)",
    gestational_age_range = "25-41 weeks (median 29)",
    renal_function   = "Serum creatinine median 49 umol/L (range 18-172)",
    n_concentrations = 347L,
    notes            = "Patient characteristics from Grimsley 1999 Table 2. 70 patients initially eligible; 11 excluded (1 vancomycin discontinued before sampling, 6 with incomplete dose details, 3 over 3 months old, 1 outlier with rapidly resolving acute renal injury). 347 concentrations comprised 153 peaks (drawn 1 h after end of infusion), 183 troughs (end of dose interval), and 11 mid-dose samples. Modeling done in NONMEM Version IV; one-compartment final model accepted after the two-compartment full-covariate model showed only marginal advantage on residual plots."
  )

  ini({
    # Structural parameters (Grimsley 1999 Table 4 final-model column).
    # The paper publishes the model directly as
    #   CL (L/h) = 3.56 * WT (kg) / CREAT (umol/L)
    #   V  (L)   = 0.669 * WT (kg)
    # with no separately estimated covariate exponents (WT enters CL and V
    # linearly; CREAT enters CL as a 1/CREAT term that the paper compared
    # against an estimated power form and judged equivalent). The typical-
    # value parameters below carry the published coefficients verbatim;
    # WT and CREAT are applied multiplicatively in model() to reproduce
    # Table 4's structural equations.
    lcl <- log(3.56);  label("CL coefficient on WT/CREAT (L*umol/L/h/kg)")  # Grimsley 1999 Table 4: CL = 3.56 * weight / creatinine (SE 3.8%)
    lvc <- log(0.669); label("V per kg (L/kg)")                              # Grimsley 1999 Table 4: V  = 0.669 * weight (SE 4.7%)

    # Inter-individual variability (Grimsley 1999 Table 4 final-model column,
    # reported as percent coefficient of variation; omega^2 = log(CV^2 + 1)
    # for log-normal etas).
    etalcl ~ 0.04729 # log(0.22^2 + 1); 22% CV on CL (SE 28%)
    etalvc ~ 0.03188 # log(0.18^2 + 1); 18% CV on V  (SE 46%)

    # Additive residual error (Grimsley 1999 Table 4 final-model column).
    # The paper notes the residual structure was best described by an
    # additive form rather than proportional or combined.
    addSd <- 4.53; label("Additive residual error (mg/L)") # Grimsley 1999 Table 4: residual error = 4.53 mg/L (SE 13%)
  })
  model({
    # Individual PK parameters. CL is the published coefficient (3.56)
    # multiplied by the covariate combination WT / CREAT; V is the
    # published per-kg coefficient (0.669) multiplied by WT. There are
    # no separately estimated covariate exponents in this model.
    cl <- exp(lcl + etalcl) * WT / CREAT
    vc <- exp(lvc + etalvc) * WT

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, vc in L -> central/vc has units mg/L.
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
