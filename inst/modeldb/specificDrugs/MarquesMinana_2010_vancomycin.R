MarquesMinana_2010_vancomycin <- function() {
  description <- "One-compartment IV-infusion population PK model for vancomycin in neonates (Marques-Minana 2010). Developed from 70 NICU neonates (postmenstrual age 25.1-48.1 weeks; weight 0.7-3.7 kg). Weight-normalized clearance is linear in postmenstrual age and increased by concomitant amoxicillin-clavulanic acid; weight-normalized volume of distribution is decreased by concomitant spironolactone. Additive interindividual variability on CL and V per the paper's Step 4 error-model selection; additive residual error."
  reference <- paste(
    "Marques-Minana M-R, Saadeddin A, Peris J-E.",
    "Population pharmacokinetic analysis of vancomycin in neonates.",
    "A new proposal of initial dosage guideline.",
    "Br J Clin Pharmacol. 2010;70(5):713-722.",
    "doi:10.1111/j.1365-2125.2010.03736.x",
    sep = " "
  )
  vignette <- "MarquesMinana_2010_vancomycin"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Marques-Minana 2010 Table 2 original-group: mean 1.7 kg",
        "(range 0.7-3.7, SD 0.8). Enters CL and V_d linearly with the published",
        "per-kg coefficients (no separately estimated allometric exponent);",
        "the structural equations TVCL = (theta1 * PMA * (1 + theta2 * AMX)) * wt",
        "and TVV = (theta3 * (1 - theta4 * SPI)) * wt come directly from",
        "Table 3 legend."
      ),
      source_name        = "weight"
    ),
    PAGE = list(
      description        = "Postmenstrual age (gestational age + postnatal age)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Marques-Minana 2010 Table 2 original-group: mean 34.6",
        "weeks (range 25.1-48.1, SD 5.3). The source paper carries PMA in weeks",
        "and the published TVCL formula uses PMA in weeks directly; this model",
        "converts the canonical PAGE (months) to weeks via pma_wk = PAGE * 4.35",
        "so the published theta1 = 1.92e-3 L/h/kg/(week of PMA) coefficient",
        "applies unchanged. Reference category n/a (linear effect, no centring).",
        "Drives CL maturation through the linear PMA term in Table 3."
      ),
      source_name        = "PMA (weeks)"
    ),
    CONMED_AMOXCLAV = list(
      description        = "Concomitant amoxicillin-clavulanic acid coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant amoxicillin-clavulanic acid)",
      notes              = paste(
        "Time-fixed per subject's vancomycin course. Marques-Minana 2010",
        "Table 3 theta2 = 0.650: multiplicative deviation on CL when",
        "amoxicillin-clavulanic acid is coadministered (+65% CL vs no",
        "coadministration). Mechanistic rationale (Discussion): amoxicillin",
        "is a substrate of the renal peptide transporters hPepT1 / hPepT2",
        "(Li et al. 2006); since vancomycin undergoes tubular reabsorption,",
        "amoxicillin may competitively inhibit that reabsorption and thereby",
        "increase vancomycin renal clearance. Population prevalence not",
        "tabulated in the paper (one of 21 co-administered drugs evaluated)."
      ),
      source_name        = "AMX"
    ),
    CONMED_SPIRON = list(
      description        = "Concomitant spironolactone coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant spironolactone)",
      notes              = paste(
        "Time-fixed per subject's vancomycin course. Marques-Minana 2010",
        "Table 3 theta4 = 0.344: multiplicative fractional decrease in",
        "weight-normalized V_d when spironolactone is coadministered",
        "(V_d *= (1 - 0.344 * CONMED_SPIRON), i.e. -34.4% V_d). Mechanistic",
        "rationale (Discussion): spironolactone reduces total body water",
        "/ extracellular fluid volume, and hydrophilic antibiotics like",
        "vancomycin distribute primarily into ECF, so the V_d shrinks.",
        "Population prevalence not tabulated in the paper."
      ),
      source_name        = "SPI"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 70L,
    n_studies        = 1L,
    age_range        = "postmenstrual age 25.1-48.1 weeks; gestational age 24.0-42.0 weeks; postnatal age 4.0-63.0 days",
    age_median       = "postmenstrual age 34.6 weeks; gestational age 32.2 weeks; postnatal age 16.9 days",
    weight_range     = "0.7-3.7 kg",
    weight_median    = "1.7 kg",
    sex_female_pct   = 53,
    race_ethnicity   = "Not reported (single-center cohort at 'Hospital La Fe', Valencia, Spain)",
    disease_state    = "Neonates in the neonatal intensive care unit with proven or suspected gram-positive infections treated with intravenous vancomycin",
    dose_range       = "15 mg/kg per dose IV (constant-rate infusion over 60 min); dosing interval 8, 12, 18, or 24 h depending on weight and postnatal age (Table 1)",
    regions          = "Spain ('Hospital La Fe', Valencia)",
    gestational_age_range = "24.0-42.0 weeks (mean 32.2, SD 5.0)",
    postnatal_age_range   = "4.0-63.0 days (mean 16.9, SD 10.9)",
    bsa_range        = "0.1-0.3 m^2 (mean 0.15, SD 0.04)",
    notes            = paste(
      "Patient characteristics from Marques-Minana 2010 Table 2 (original",
      "group n = 70; a separate validation group of 41 neonates with similar",
      "demographics was used for external evaluation but not for model",
      "estimation). Vancomycin assayed by FPIA (TDxFLx, Abbott; CV < 6%).",
      "Two serum samples per patient around the third dose: peak (3 h after",
      "infusion completion) and trough (immediately before next dose). NONMEM",
      "version VI with PREDPP ADVAN1 / TRANS2 subroutines (one-compartment,",
      "CL / V parameterization). The basic model was a one-compartment IV",
      "infusion fit to CL and V with proportional IIV; covariates were then",
      "screened (Step 2) and added sequentially (Step 3) before evaluating",
      "alternative IIV / residual error structures (Step 4); the additive",
      "interindividual error model was selected in Step 4 and survived",
      "backward elimination unchanged in Step 5. Bootstrap of 1000 resampled",
      "datasets (92% success) confirmed the final estimates (Table 3)."
    )
  )

  ini({
    # ===== Structural parameters (Marques-Minana 2010 Table 3 final model) =====
    # The published TVCL formula is TVCL = (theta1 * PMA * (1 + theta2 * AMX)) * wt
    # with PMA in weeks, wt in kg, theta1 in L/h/kg/(week of PMA), theta2 dimensionless.
    # The published TVV formula is TVV = (theta3 * (1 - theta4 * SPI)) * wt with
    # theta3 in L/kg, theta4 dimensionless. Storing the structural coefficients on
    # log scale keeps the canonical 'lcl' / 'lvc' naming; the actual covariate
    # multiplication happens in model() below.
    lcl <- log(1.92e-3); label("Log CL coefficient on PMA*wt (log(L/h/kg/week))") # Marques-Minana 2010 Table 3: theta1 = 1.92e-3 (95% CI 1.72e-3, 2.12e-3)
    lvc <- log(0.572);   label("Log V_d per kg (log(L/kg))")                      # Marques-Minana 2010 Table 3: theta3 = 0.572  (95% CI 0.505, 0.639)

    # ===== Covariate effects (Marques-Minana 2010 Table 3 final model) =====
    e_amx_cl <- 0.650;   label("Multiplicative increase in CL with amoxicillin-clavulanic acid (unitless)") # Table 3: theta2 = 0.650 (95% CI 0.391, 0.909)
    e_spi_vc <- 0.344;   label("Fractional decrease in V_d with spironolactone (unitless)")                 # Table 3: theta4 = 0.344 (95% CI 0.242, 0.446)

    # ===== Inter-individual variability (Marques-Minana 2010 Table 3 final model) =====
    # The paper explicitly selected an additive interindividual error model in
    # Step 4 (Methods 'Pharmacokinetic analysis' / Results paragraph 2); the
    # reported w_CL (CV%) = 35.6 and w_V_d (CV%) = 19.3 are SD / TVP * 100, i.e.
    # the additive ETA standard deviation expressed as %CV relative to the
    # typical-individual TVCL / TVV (PMA = 34.6 wk, wt = 1.7 kg, no AMX, no SPI).
    # Typical CL = 1.92e-3 * 34.6 * 1.7 = 0.1129 L/h -> omega_CL  = 0.356 * 0.1129 = 0.0402 L/h -> omega_CL^2 = 1.616e-3 (L/h)^2.
    # Typical V  = 0.572 * 1.7         = 0.9724 L  -> omega_Vd  = 0.193 * 0.9724 = 0.1877 L   -> omega_Vd^2 = 3.522e-2 L^2.
    # The eta names follow the canonical 'etal<param>' convention to match
    # 'lcl' / 'lvc' (the lint expects matching prefixes). The application in
    # model() is ADDITIVE on the linear-scale cl / vc (cl = exp(lcl) * ... + etalcl),
    # not multiplicative on the log scale -- see model() comment and the
    # vignette's Assumptions & deviations section. Precedent: Sheng_2013_ABPM.R
    # encodes additive IIV on linear-scale parameters; here the same scheme is
    # used while keeping the canonical 'l' prefix on the parameter and eta
    # names so the IIV-name convention is satisfied.
    etalcl ~ 1.616e-3   # SD = 0.0402 L/h, additive on linear cl (35.6% CV of typical CL = 0.1129 L/h)
    etalvc ~ 3.522e-2   # SD = 0.1877 L,  additive on linear vc (19.3% CV of typical V  = 0.9724 L)

    # ===== Residual error (Marques-Minana 2010 Table 3 final model) =====
    # The published sigma (s) = 2.69 microgram/mL carries concentration units,
    # indicating an additive residual error model on the linear concentration
    # scale (consistent with the paper's vancomycin assay quantifying serum
    # concentrations in microgram/mL). 95% CI 0.73; 3.73 from Table 3.
    addSd <- 2.69; label("Additive residual error (mg/L)") # Marques-Minana 2010 Table 3: sigma = 2.69 (95% CI 0.73, 3.73)
  })

  model({
    # ----- Derived covariate terms -----
    # Convert canonical PAGE (months) back to source-paper PMA (weeks) so the
    # published theta1 coefficient and PMA term apply unchanged.
    pma_wk <- PAGE * 4.35

    # ----- Individual PK parameters -----
    # Additive interindividual variability on the linear-scale parameter
    # (paper Step 4): cl_i = TVCL + ETA_CL, vc_i = TVV + ETA_V. The exp(lcl)
    # / exp(lvc) terms recover the linear-scale structural coefficient from
    # its log-scale storage. The eta names (etalcl / etalvc) carry the 'l'
    # prefix to match the parameter names per the IIV-naming convention,
    # but the etas are added linearly OUTSIDE the exp() to preserve the
    # paper's additive form -- not the more common multiplicative form
    # exp(lcl + etalcl) that maps log-scale etas back to lognormal IIV.
    cl <- exp(lcl) * pma_wk * WT * (1 + e_amx_cl * CONMED_AMOXCLAV) + etalcl
    vc <- exp(lvc) * WT             * (1 - e_spi_vc * CONMED_SPIRON) + etalvc

    # ----- Micro-constant -----
    kel <- cl / vc

    # ----- ODE system (one-compartment IV) -----
    # Vancomycin is administered as a constant-rate IV infusion over 60 min
    # (paper Methods: dose delivered by syringe pump). Dosing event rate
    # information is supplied by the simulated event table; NONMEM PREDPP
    # ADVAN1 / TRANS2 corresponds to this one-compartment central-only ODE.
    d/dt(central) <- -kel * central

    # ----- Observation and error -----
    # Dose in mg, vc in L -> central/vc has units mg/L (equivalent to
    # microgram/mL, the published assay unit).
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
