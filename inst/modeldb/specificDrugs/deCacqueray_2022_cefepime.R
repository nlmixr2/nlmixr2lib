deCacqueray_2022_cefepime <- function() {
  description <- "One-compartment IV population PK model for cefepime in 59 critically ill infants and children aged 1.1 months to 17.6 years (de Cacqueray 2022); body-weight allometric scaling referenced to 9 kg (fixed exponents 0.75 on CL, 1 on Vc), a power effect of Schwartz-estimated glomerular filtration rate on CL referenced to 153 mL/min/1.73 m2, and correlated log-normal between-subject variability on CL and Vc. Parameters transcribed from the secondary source Gotta 2025, which reprints the final-model equations and deposits the literal Monolix/Simulx source used to simulate this model; the de Cacqueray 2022 primary is not open access."
  reference <- paste(
    "de Cacqueray N, Hirt D, Zheng Y, Bille E, Leger PL, Rambaud J,",
    "Toubiana J, Chosidow A, Vimont S, Callot D, et al. Cefepime",
    "population pharmacokinetics and dosing regimen optimization in",
    "critically ill children with different renal function.",
    "Clin Microbiol Infect. 2022;28(10):1389.e1-1389.e7.",
    "doi:10.1016/j.cmi.2022.05.007.",
    "Parameter values transcribed from the secondary source",
    "Gotta V, Csajka C, Glauser A, Berger C, Pfister M, Paioni P.",
    "Risk of potentially neurotoxic exposure in infants under high-dose",
    "cefepime treatment - a pharmacometric simulation study.",
    "Pharmaceutics. 2025;17(5):544. doi:10.3390/pharmaceutics17050544",
    "(Equations 3 and 4, Section 2.2.2, and Supplemental Data S1).",
    sep = " "
  )
  vignette <- "Gotta_2025_cefepime_infant_neurotoxicity"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "cefepime", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling referenced to 9 kg on both CL (exponent 0.75)",
        "and Vc (exponent 1). Gotta 2025 Eq. 3 gives",
        "CL = 1.21 * (weight/9)^0.75 and Eq. 4 gives V = 4.8 * (weight/9);",
        "the Supplemental Data S1 Monolix model file spells the same two",
        "expressions out as Cltyp and Vtyp. The de Cacqueray 2022 cohort",
        "had an inclusion criterion of weight > 3 kg, so extrapolation",
        "below 3 kg is outside the model-development range.",
        "Source column weight_kg in the deposited covariate file."
      ),
      source_name        = "weight_kg"
    ),
    CRCL = list(
      description        = paste(
        "Estimated glomerular filtration rate, BSA-normalized, computed with",
        "the Schwartz (1976) paediatric equation from body length and plasma",
        "creatinine"
      ),
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power effect on CL: (CRCL / 153)^0.37 per Gotta 2025 Eq. 3.",
        "The reference value 153 mL/min/1.73 m^2 is supranormal because the",
        "de Cacqueray 2022 cohort were critically ill children, in whom",
        "augmented renal clearance is common; it is the cohort-centring",
        "value carried through in the deposited Monolix source",
        "(Gotta 2025 Supplemental Data S1, `Cltyp = Clpop *",
        "(weight_kg/9)^0.75 * (eGFR_Schwartz1976/153)^theGFR`).",
        "Gotta 2025 states that the estimating equation is Schwartz 1976",
        "(their reference 30); the Schwartz 1984 full-term-infant variant",
        "was used only in their sensitivity analysis, not with this model.",
        "Source column eGFR_Schwartz1976 in the deposited covariate file."
      ),
      source_name        = "eGFR_Schwartz1976"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 59L,
    n_studies      = 1L,
    age_range      = "1.1 months - 17.6 years",
    weight_range   = "Inclusion criterion weight > 3 kg; individual weights not reported in the secondary source",
    disease_state  = paste(
      "Critically ill paediatric patients with different renal function,",
      "mainly with lung or bloodstream infections"
    ),
    n_observations = 129L,
    notes          = paste(
      "Population description transcribed from Gotta 2025 Section 2.2.2,",
      "which summarises de Cacqueray 2022. 129 plasma concentrations were",
      "collected from 59 patients; the sampling times are not reported in",
      "the secondary source, and Gotta 2025 characterise the underlying",
      "design as sparse. The de Cacqueray 2022 primary publication is not",
      "open access and was not available on disk when this model was",
      "extracted, so no baseline demographic table, no covariate ranges,",
      "and no parameter uncertainty (RSE / confidence intervals) could be",
      "recorded. Neonates were NOT included in model development, which is",
      "why Gotta 2025 do not apply this model below 1 month of age",
      "(their Table 2 footnote). Re-extract from the primary when it",
      "becomes obtainable."
    )
  )

  ini({
    # Structural PK -- Gotta 2025 Eq. 3 and Eq. 4 (Section 2.2.2), confirmed
    # against the population parameter vector in the deposited Simulx project
    # (Supplemental Data S1): 'PopParameters' = {{{values={1.21, 4.8, 0.39,
    # 0.35, 0.5, 0.37, 0.39}}}} for parameters {Clpop, Vpop, sdCl, sdV,
    # corr_Cl_V, theGFR, b}.
    lcl <- log(1.21); label("Typical CL at WT = 9 kg and CRCL = 153 mL/min/1.73 m2 (L/h)") # Gotta 2025 Eq. 3; Suppl. Data S1 Clpop = 1.21
    lvc <- log(4.8);  label("Typical Vc at WT = 9 kg (L)")                                 # Gotta 2025 Eq. 4; Suppl. Data S1 Vpop = 4.8

    # Allometric exponents. Gotta 2025 Section 2.2.2 describes these as
    # "standard weight-based allometric scaling", and both Eq. 3 and the
    # deposited Monolix source hard-code 0.75 on CL and 1 (implicit, the
    # ratio appears un-exponentiated) on V, i.e. neither was estimated.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL (unitless)") # Gotta 2025 Eq. 3 (hard-coded); Suppl. Data S1 `(weight_kg/9)^0.75`
    e_wt_vc <- fixed(1.00); label("Allometric exponent on Vc (unitless)") # Gotta 2025 Eq. 4 (linear in weight); Suppl. Data S1 `Vpop * (weight_kg/9)`

    # Renal-function covariate effect on CL.
    e_crcl_cl <- 0.37; label("Power exponent on (CRCL / 153) for CL (unitless)") # Gotta 2025 Eq. 3; Suppl. Data S1 theGFR = 0.37

    # Reference (centring) covariate values -- structural constants of the
    # published equations, not estimated quantities.
    wt_ref   <- fixed(9);   label("Reference body weight (kg)")                          # Gotta 2025 Eq. 3 and Eq. 4 denominators
    crcl_ref <- fixed(153); label("Reference estimated GFR (mL/min/1.73 m^2)")           # Gotta 2025 Eq. 3 denominator

    # Inter-individual variability. The deposited Monolix source declares
    #   Cl = {distribution=lognormal, typical=Cltyp, sd=sdCl}
    #   V  = {distribution=lognormal, typical=Vtyp,  sd=sdV}
    #   correlation = {level=id, r(Cl, V)=corr_Cl_V}
    # so sdCl / sdV are standard deviations of the Gaussian random effects
    # (i.e. omega on the log scale), NOT lognormal CVs. That settles the
    # otherwise-ambiguous "39% in CL and 35% in V" of Section 2.2.2:
    #   var(etalcl) = 0.39^2 = 0.1521
    #   var(etalvc) = 0.35^2 = 0.1225
    #   cov         = 0.5 * 0.39 * 0.35 = 0.06825   (r = corr_Cl_V = 0.5)
    # The 2x2 block is positive definite (determinant 0.013974).
    #
    # PROVENANCE CAVEAT: the correlation r(Cl,V) = 0.5 appears only in the
    # Simulx project file deposited with Gotta 2025; neither Gotta 2025's
    # main text nor its supplemental methods attribute it to de Cacqueray
    # 2022. It may therefore be an assumption made by Gotta et al. rather
    # than a de Cacqueray estimate. It is retained here because it is the
    # value on disk and the value that produced Gotta 2025's published
    # Table 2 percentages; verify against the primary when obtainable.
    etalcl + etalvc ~ c(0.1521,
                        0.06825, 0.1225) # Gotta 2025 Section 2.2.2 (39% CL, 35% Vc) + Suppl. Data S1 sdCl = 0.39, sdV = 0.35, corr_Cl_V = 0.5

    # Residual error -- Gotta 2025 Section 2.2.2 "Residual intra-individual
    # variability was 39%", encoded in the deposited Monolix source as
    #   y = {distribution = normal, prediction = Cc,
    #        errorModel = proportional(b)}   with b = 0.39.
    propSd <- 0.39; label("Proportional residual error (fraction)") # Gotta 2025 Section 2.2.2; Suppl. Data S1 b = 0.39
  })

  model({
    # ----- Derived covariate terms -----
    # Allometric size on CL and Vc, and the renal-function power term on CL.
    wt_cl_factor   <- (WT / wt_ref)^e_wt_cl
    wt_vc_factor   <- (WT / wt_ref)^e_wt_vc
    crcl_cl_factor <- (CRCL / crcl_ref)^e_crcl_cl

    # ----- Individual PK parameters -----
    cl <- exp(lcl + etalcl) * wt_cl_factor * crcl_cl_factor
    vc <- exp(lvc + etalvc) * wt_vc_factor

    # ----- Micro-constants -----
    kel <- cl / vc

    # ----- ODE system (one-compartment IV) -----
    d/dt(central) <- -kel * central

    # ----- Observation -----
    # Dose in mg, vc in L -> central / vc has units mg/L (= ug/mL).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
