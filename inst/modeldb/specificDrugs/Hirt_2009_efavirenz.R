Hirt_2009_efavirenz <- function() {
  description <- "One-compartment population PK model with first-order absorption and elimination for once-daily oral efavirenz (EFV) in treatment-naive HIV-1-infected West African children (Hirt 2009). CL/F and V/F scale linearly with body weight (shared allometric exponent fixed at 1) and CL/F additionally varies with postnatal age via a power covariate centred at the cohort median 6.35 years (signed exponent -0.535, so apparent clearance decreases with age); the inter-individual variability of V/F is forced to perfect correlation with the eta of CL/F and is constructed as vc_eta_scale * etalcl (the K parameter in Hirt 2009 Table 2); multiplicative residual error."
  reference <- "Hirt D, Urien S, Olivier M, Peyriere H, Nacro B, Diagbouga S, Zoure E, Rouet F, Hien H, Msellati P, Van De Perre P, Treluyer JM. Is the recommended dose of efavirenz optimal in young West African human immunodeficiency virus-infected children? Antimicrob Agents Chemother. 2009;53(10):4407-4413. doi:10.1128/AAC.01594-08"
  vignette <- "Hirt_2009_efavirenz"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear scaling on CL/F and V/F (allometric exponent fixed at 1) with reference weight 16.4 kg (cohort median). The body-weight exponent on CL/F was estimated at 1.13 in a sensitivity analysis but the final model was fitted with the exponent fixed to 1 because dosing is given linearly in mg/kg (Hirt 2009 Discussion paragraph 3 on page 4410).",
      source_name        = "BW"
    ),
    AGE = list(
      description        = "Postnatal age at study entry",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at study entry. Power covariate on CL/F centred at the cohort median 6.35 years: CL/F scales by (AGE / 6.35)^e_age_cl with e_age_cl negative. Apparent per-kg clearance therefore decreases as age increases above the median, consistent with the maturation of hepatic CYP2B6 activity declining from its early-childhood peak toward adult levels (Hirt 2009 Results paragraph 'Population pharmacokinetics' and Discussion paragraph iii on page 4410).",
      source_name        = "age"
    )
  )

  covariatesDataExcluded <- list(
    HT = list(
      description        = "Body height (size)",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested as a covariate on CL/F and V/F via generalized additive modelling in the basic model but did not meet the inclusion criteria (objective-function decrease >= 6.63 and a reduction in inter-subject variability) in the upward model-building procedure (Hirt 2009 Methods 'Modeling strategy and population pharmacokinetic model').",
      source_name        = "size"
    ),
    CREAT = list(
      description        = "Basal serum creatinine",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested as a covariate on CL/F and V/F but not retained (Hirt 2009 Methods). Median 63 umol/L, range 6.6-126.4 (Table 1).",
      source_name        = "serum creatinine concn"
    ),
    ALT = list(
      description        = "Alanine aminotransferase",
      units              = "IU/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested as a covariate on CL/F and V/F but not retained (Hirt 2009 Methods). Median 26 IU/L, range 7-92 (Table 1).",
      source_name        = "alanine aminotransferase concn"
    ),
    AST = list(
      description        = "Aspartate aminotransferase",
      units              = "IU/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested as a covariate on CL/F and V/F but not retained (Hirt 2009 Methods). Median 44 IU/L, range 12-200 (Table 1).",
      source_name        = "aspartate aminotransferase concn"
    ),
    TBILI = list(
      description        = "Total bilirubin",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested as a covariate on CL/F and V/F but not retained (Hirt 2009 Methods). Median 6.5 umol/L, range 0-42.8 (Table 1).",
      source_name        = "total bilirubin concn"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 48L,
    n_studies      = 1L,
    n_profiles     = 200L,
    age_range      = "2.77-14.70 years",
    age_median     = "6.35 years",
    weight_range   = "11-37 kg",
    weight_median  = "16.4 kg",
    sex_female_pct = NA_real_,
    race_ethnicity = c(Black = 100),
    disease_state  = "Treatment-naive HIV-1 infection; CDC clinical category C and/or low CD4 (children <=5 years) or CD4 <=200/uL (>5 years); or category B/A/N with intermediate CD4 plus viral load >100,000 copies/mL.",
    co_medication  = "Once-daily didanosine (ddI; 240 mg/m^2 BSA) plus lamivudine (3TC; 8 mg/kg) given alongside once-daily efavirenz (the ddI-3TC-EFV combination).",
    dose_range     = "Once-daily oral EFV as 200 mg capsules administered by body-weight-band per pediatric label. 46 of 48 children received the label-recommended dose; cohort median dose 250 mg (14.4 mg/kg) once daily. EFV administered at 18:00.",
    regions        = "Burkina Faso (Bobo-Dioulasso); the BURKINAME-ANRS 12103 phase II open-label trial.",
    sampling_window = "Week 2 (n = 48): pre-dose, 1 h, and 3 h post-dose. Months 2 to 5 (n = 9): pre-dose, 1, 2, 3, 6, 12, and 24 h post-dose. 200 EFV plasma concentrations total measured by HPLC-UV (LLOQ 0.5 mg/L); no observations below LLOQ.",
    notes          = "NIH ClinicalTrials.gov NCT00122538. Eligible age range 30 months to 15 years; baseline weight >= 10 kg required. All-Black African cohort; ethnicity could not be tested as a covariate (single-stratum design). Sex distribution was not tabulated in the publication."
  )

  ini({
    # Structural fixed-effect estimates from Hirt 2009 Table 2 (final model).
    # Typical-value parameters are reported in the paper as per-kg quantities;
    # encoded here at the cohort-median reference body weight (16.4 kg) with
    # an explicit linear (exponent = 1) weight covariate to keep the registry
    # convention (typical CL/F and V/F at a reference WT, allometric power on
    # WT in model()). 0.21 L/h/kg * 16.4 kg = 3.444 L/h; 4.48 L/kg * 16.4 kg
    # = 73.472 L.
    lka       <- log(0.45)    ; label("Absorption rate constant ka (1/h)")                                              # Table 2 final ka = 0.45 1/h (RSE 23%)
    lcl       <- log(3.444)   ; label("Apparent clearance CL/F at WT = 16.4 kg, AGE = 6.35 years (L/h)")                # Table 2 final CL/F = 0.21 L/h/kg (RSE 9%) x 16.4 kg
    lvc       <- log(73.472)  ; label("Apparent central volume V/F at WT = 16.4 kg (L)")                                # Table 2 final V/F = 4.48 L/kg (RSE 14%) x 16.4 kg

    # Covariate effects.
    e_wt_cl_vc <- fixed(1)    ; label("Shared linear (power = 1) WT exponent on CL/F and V/F (unitless; fixed)")        # Hirt 2009 Discussion paragraph 3 on page 4410 (BW exponent fixed to 1; estimated value was 1.13)
    # Power exponent for the AGE covariate on CL/F. The paper writes the
    # equation as CL/F = theta_CL/F * BW / (AGE/6.35)^0.535 (the age term in
    # the denominator); equivalent to CL/F * (AGE/6.35)^(-0.535), so the
    # canonical signed exponent is -0.535. Apparent per-kg clearance therefore
    # decreases as age increases above the cohort median.
    e_age_cl   <- -0.535      ; label("Age power exponent on CL/F: CL/F * (AGE / 6.35)^e_age_cl (unitless)")            # Table 2 final theta_age = 0.54 (RSE 28%); paper text page 4410 reports 0.535 with the age term in the denominator

    # Eta-coupling constant: eta_V/F is forced to perfect correlation with
    # eta_CL/F and constructed as vc_eta_scale * etalcl. The single estimated
    # scaling parameter is "K" in Hirt 2009 Table 2 (the observed raw
    # correlation between eta_CL/F and eta_V/F before reparameterisation was
    # r = 0.99; the paper fixed it to exactly 1 and estimated the SD ratio).
    # Stored as a structural theta rather than as a separate etalvc slot,
    # because for a correlation of exactly 1 the two forms are mathematically
    # equivalent and the deterministic form avoids the singular-covariance
    # numerical issue that a 1-correlation block would otherwise trigger
    # (same encoding choice as Prytula 2016 tacrolimus).
    vc_eta_scale <- 0.64     ; label("Scaling factor K relating eta_V/F to eta_CL/F (correlation fixed to 1; eta_V/F = vc_eta_scale * etalcl)")  # Table 2 final K = 0.64 (RSE 10%)

    # IIV on CL/F. The paper reports omega(CL/F) as 61% CV; converting to
    # log-scale variance: omega^2 = log(1 + CV^2) = log(1 + 0.61^2) = 0.31634.
    etalcl ~ 0.31634                                                                                                    # Table 2 final omega(CL/F) = 61% CV; omega^2 = log(1 + 0.61^2)

    # Residual error. Hirt 2009 Methods 'Population pharmacokinetics' reports
    # the residual variability was best described by a multiplicative
    # (proportional) error model with sigma = 35% (RSE 13%).
    propSd <- 0.35           ; label("Proportional residual error (fraction)")                                          # Table 2 final sigma = 35% (RSE 13%)
  })

  model({
    # Allometric scaling: linear in body weight (e_wt_cl_vc fixed to 1)
    # referenced at the cohort median 16.4 kg.
    wt_ref <- WT / 16.4

    # Individual PK parameters.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * wt_ref^e_wt_cl_vc * (AGE / 6.35)^e_age_cl
    vc <- exp(lvc + vc_eta_scale * etalcl) * wt_ref^e_wt_cl_vc

    kel <- cl / vc

    # One-compartment with first-order absorption. Dose lands in `depot`;
    # bioavailability F is implicit in the apparent CL/F and V/F.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Concentration: dose in mg, vc in L -> mg/L (the EFV concentration unit
    # used throughout Hirt 2009).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
