Bouazza_2010_lamivudine <- function() {
  description <- "Two-compartment population PK model for once-daily oral lamivudine in HIV-infected West African children (Bouazza 2010); allometric weight scaling on CL/F, Q/F, Vc/F, and Vp/F with reference body weight 16.8 kg, and absorption rate constant Ka structurally fixed to the disposition distribution-phase eigenvalue (Ka = alpha = 0.71 1/h) from the literature"
  reference <- "Bouazza N, Hirt D, Bardin C, Diagbouga S, Nacro B, Hien H, Zoure E, Rouet F, Ouiminga A, Blanche S, Van De Perre P, Treluyer J-M, Msellati P, Urien S. Is the recommended once-daily dose of lamivudine optimal in West African HIV-infected children? Antimicrob Agents Chemother. 2010;54(9):3938-3943. doi:10.1128/AAC.00306-10"
  vignette <- "Bouazza_2010_lamivudine"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline (single measurement per child); used for allometric scaling on CL/F and Q/F (exponent 0.75, fixed) and on Vc/F and Vp/F (exponent 1, fixed) with reference weight 16.8 kg (cohort median, Methods 'Modeling strategy').",
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 45,
    n_studies      = 1,
    age_range      = "2.5-14 years (median 6.75)",
    age_median     = "6.75 years",
    weight_range   = "11-37 kg",
    weight_median  = "16.8 kg",
    sex_female_pct = 38,
    race_ethnicity = "West African (sub-Saharan); ethnic stratification not reported in source.",
    disease_state  = "Antiretroviral-naive HIV-1 infected children, CDC clinical stage A-C, with severe immunological / virological status meeting Burkina Faso HAART initiation criteria.",
    dose_range     = "Lamivudine oral, 8 mg/kg once daily as 150 mg tablet or 10 mg/ml oral solution; median dose 150 mg (range 90-300 mg). Co-administered with didanosine 240 mg/m^2 q.d. and weight-band efavirenz q.d.",
    regions        = "Burkina Faso (Bobo-Dioulasso); BURKINAME-ANRS 12103 trial (ClinicalTrials.gov NCT00122538).",
    n_observations = 148,
    notes          = "Forty-nine children enrolled, 45 evaluable for PK (17 girls / 28 boys). Sampling schedule: pre-dose, 1 h, 3 h post-dose (39 children) or pre-dose, 1, 2, 3, 6, 12, 24 h post-dose (10 children). Sampling began on day 15 of treatment for 38 children and between months 2-5 of treatment for 11 children -- assumed to be at steady state for the model. Demographics in Table 1 of the source."
  )

  ini({
    # Structural PK parameters -- typical values at reference body weight 16.8 kg.
    # Two-compartment model with the structural assumption Ka = alpha (distribution
    # eigenvalue), which is required by the sparse early-time sampling design (no
    # samples during the absorption phase). See Methods 'Modeling strategy'.
    lka <- fixed(log(0.71)); label("Absorption rate constant (1/h); fixed to the published value from reference 17 (= the model's alpha disposition eigenvalue)") # Table 2 K_a; Methods 'Modeling strategy' fixes K_a per ref. 17
    lcl <- log(16.9);        label("Apparent clearance at WT=16.8 kg (CL/F, L/h)") # Table 2 CL/F
    lvc <- log(30.8);        label("Apparent central volume at WT=16.8 kg (Vc/F, L)") # Table 2 Vc/F
    lvp <- log(58.6);        label("Apparent peripheral volume at WT=16.8 kg (Vp/F, L)") # Table 2 Vp/F
    lq  <- log(4.48);        label("Apparent intercompartmental clearance at WT=16.8 kg (Q/F, L/h)") # Table 2 Q/F

    # Allometric exponents -- fixed at canonical theoretical values per Methods
    # ("from allometric scaling theory, these are typically 0.75 for clearance
    # parameters and 1 for volumes of distribution (2)").
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL/F and Q/F (unitless; fixed at 0.75)") # Methods 'Modeling strategy'
    e_wt_vc <- fixed(1);    label("Allometric exponent on Vc/F and Vp/F (unitless; fixed at 1)")   # Methods 'Modeling strategy'

    # IIV on CL/F only (exponential model; only IIV retained in final model).
    # omega CL/F = 0.30 in Table 2 is reported as the SD on the log scale
    # (approximate CV). Same-author group precedent (Bouazza 2011 lamivudine,
    # PMID 21576437, Table 1 omega CL/F = 0.32 matches the abstract's "32%"
    # IIV statement) supports this reading. nlmixr2's `~` syntax expects the
    # variance, so the encoded value is 0.30^2 = 0.09.
    etalcl ~ 0.09 # Table 2 omega(CL/F) = 0.30 (SD); variance = 0.30^2 = 0.09

    # Residual error -- proportional ("multiplicative error model" per Results).
    # sigma = 0.60 in Table 2 is the proportional SD (fractional 60% CV).
    propSd <- 0.60; label("Proportional residual error (fraction)") # Table 2 sigma
  })

  model({
    # Individual PK parameters (allometric weight scaling around 16.8 kg).
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 16.8)^e_wt_cl
    vc <- exp(lvc)          * (WT / 16.8)^e_wt_vc
    vp <- exp(lvp)          * (WT / 16.8)^e_wt_vc
    q  <- exp(lq)           * (WT / 16.8)^e_wt_cl

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # Concentration: dose in mg, Vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
