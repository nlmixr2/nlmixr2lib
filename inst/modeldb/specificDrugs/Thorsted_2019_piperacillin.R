Thorsted_2019_piperacillin <- function() {
  description <- paste(
    "Two-compartment IV population PK model for UNBOUND (free, ultrafiltrate)",
    "piperacillin in febrile children aged 1-18 years receiving cancer",
    "chemotherapy (Thorsted 2019). Body weight enters every clearance and",
    "volume through fixed-exponent allometry (0.75 on CL and Q, 1.0 on Vc and",
    "Vp) referenced to 70 kg. Unexplained variability in clearance is carried",
    "entirely as inter-occasion variability between febrile episodes (up to",
    "four episodes per child); inter-individual variability became",
    "insignificant once the between-episode term was added and is not part",
    "of the final model. Proportional residual error."
  )
  reference <- paste(
    "Thorsted A, Kristoffersson AN, Maarbjerg SF, Schroder H, Wang M,",
    "Brock B, Nielsen EI, Friberg LE. (2019).",
    "Population pharmacokinetics of piperacillin in febrile children",
    "receiving cancer chemotherapy: the impact of body weight and target on",
    "an optimal dosing regimen.",
    "J Antimicrob Chemother 74(10):2984-2993.",
    "doi:10.1093/jac/dkz270.",
    "Erratum: J Antimicrob Chemother 2020;75(1):254-255, doi:10.1093/jac/dkz429",
    "(corrects a Results sentence and Figure 5 only; no model parameter changed)."
  )
  vignette <- "Thorsted_2019_piperacillin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "piperacillin (unbound)", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "piperacillin (unbound)", units = "mg", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Registered at each febrile episode (Patients and methods, 'Patient",
        "population'); episode-level median 39.4 kg (IQR 24.5-51.8, range",
        "9.5-107; Table 1). Enters as fixed-exponent allometry referenced to",
        "70 kg: (WT/70)^0.75 on CL and Q, (WT/70)^1.0 on Vc and Vp (Table 2",
        "footnote b)."
      ),
      source_name = "WT"
    ),
    OCC = list(
      description = "Febrile-episode (occasion) index for the inter-occasion variability in clearance",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "An occasion is a febrile episode (Patients and methods,",
        "'Pharmacokinetic modelling': 'an occasion defined as a febrile",
        "episode'); children contributed 1-4 episodes each (Results,",
        "'Patient characteristics'), so four occasion slots are encoded.",
        "Decomposed inside model() into mutually exclusive indicators",
        "oc1..oc4 that select the per-occasion eta, because rxode2 parses",
        "but cannot simulate the `eta ~ var | occ` syntax; occasions 2-4 are",
        "fixed() to the occasion-1 variance (NONMEM $OMEGA BLOCK(1) SAME).",
        "Records carrying any other OCC value receive no IOV. Because the",
        "final model has no IIV, simulating each episode as a separate ID",
        "with OCC = 1 is statistically equivalent."
      ),
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened by SCM; linear effects of age on CL and Vc survived",
        "backwards deletion but were removed after randomization testing",
        "(Results, 'Pharmacokinetic modelling': 'only body weight was",
        "included in the final model'). Median 12 years (range 1-18; Table 1)."
      )
    ),
    CRCL = list(
      description = "Glomerular filtration rate estimated from serum creatinine by the Schwartz formula for children",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Tested on CL with several parameterisations; largest drop dOFV =",
        "-2.54, not retained. Episode-level median 172.4 (IQR 139.8-210.8,",
        "range 87.0-425.8; Table 1), i.e. a hyperfiltrating population."
      )
    ),
    PAGE = list(
      description = "Postmenstrual age",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Kidney-maturation (sigmoid PMA) function on CL tested but not",
        "significant (dOFV = -0.286; Results)."
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "male",
      notes = "Screened by SCM; not retained. 27/43 (63%) male (Table 1)."
    ),
    DIS_BACTEREMIA = list(
      description = "Bacteraemia during the febrile episode",
      units = "(binary)",
      type = "binary",
      reference_category = "no bacteraemia",
      notes = paste(
        "Screened by SCM together with neutropenia; neither retained.",
        "Bacteraemia in 10/89 episodes, neutropenia in 77/89 (Table 1)."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 43L,
    n_episodes = 89L,
    n_observations = 482L,
    n_studies = 1L,
    age_range = "1-18 years (median 12, IQR 7-14)",
    age_median = "12 years",
    weight_range = "9.5-107 kg (episode-level median 39.4, IQR 24.5-51.8)",
    weight_median = "39.4 kg",
    sex_female_pct = 37,
    race_ethnicity = "Not reported (single-centre Danish paediatric oncology department)",
    disease_state = paste(
      "Children with cancer and chemotherapy-induced fever (> 38.0 C)",
      "starting empirical piperacillin/tazobactam; neutropenia in 87% of",
      "episodes, bacteraemia in 11%."
    ),
    dose_range = paste(
      "Piperacillin/tazobactam (8:1) 100 mg/kg piperacillin as a 5-min IV",
      "infusion approximately every 8 h (300 mg/kg/day, capped at 16000",
      "mg/day)."
    ),
    renal_function = "Schwartz eGFR median 172.4 mL/min/1.73 m^2 (range 87.0-425.8)",
    regions = "Denmark (Aarhus University Hospital)",
    notes = paste(
      "482 free piperacillin serum samples (19 below the 0.5 mg/L LLOQ,",
      "handled with M3) over 89 febrile episodes (1-4 per child), April",
      "2016 - January 2018 (EudraCT 2016-00466-33). NONMEM 7.4.3, Laplacian",
      "with interaction. Concentrations are UNBOUND (UPLC after",
      "ultrafiltration), so all parameters are apparent unbound-drug values."
    )
  )

  ini({
    # Structural parameters: typical values for a 70-kg patient (Table 2,
    # footnote b).
    lcl <- log(15.4); label("Clearance at 70 kg (L/h)") # Table 2 CL = 15.4 L/h (RSE 3.7%)
    lvc <- log(16.0); label("Central volume at 70 kg (L)") # Table 2 Vc = 16.0 L (RSE 4.8%)
    lq <- log(0.237); label("Intercompartmental clearance at 70 kg (L/h)") # Table 2 Q = 0.237 L/h (RSE 9.6%)
    lvp <- log(3.40); label("Peripheral volume at 70 kg (L)") # Table 2 Vp = 3.40 L (RSE 34%)

    # Fixed allometric exponents (Methods 'Pharmacokinetic modelling'; Table 2 footnote b)
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL (unitless)") # Table 2 footnote b, fixed 0.75
    e_wt_q <- fixed(0.75); label("Allometric exponent of body weight on Q (unitless)") # Methods: fixed 0.75 for CLs
    e_wt_vc <- fixed(1.0); label("Allometric exponent of body weight on Vc (unitless)") # Table 2 footnote b, fixed 1.0
    e_wt_vp <- fixed(1.0); label("Allometric exponent of body weight on Vp (unitless)") # Methods: fixed 1.0 for Vs

    # Inter-fever-episode (occasion) variability in CL; no IIV in the final model.
    # omega^2 = log(0.166^2 + 1) = 0.027183
    etaiov_cl_1 ~ 0.027183; label("IOV (between febrile episodes) on clearance, occasion 1 (log-scale variance)") # Table 2 'CV%CL' 16.6% (RSE 11%, shrinkage 5.2%)
    etaiov_cl_2 ~ fixed(0.027183); label("IOV on clearance, occasion 2 (log-scale variance)") # same magnitude as occasion 1
    etaiov_cl_3 ~ fixed(0.027183); label("IOV on clearance, occasion 3 (log-scale variance)") # same magnitude as occasion 1
    etaiov_cl_4 ~ fixed(0.027183); label("IOV on clearance, occasion 4 (log-scale variance)") # same magnitude as occasion 1

    propSd <- 0.332; label("Proportional residual error (fraction)") # Table 2 'CV%ERR' 33.2% (RSE 4.7%)
  })

  model({
    # Occasion (febrile-episode) indicators; any OCC outside 1-4 gets no IOV.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 + oc4 * etaiov_cl_4

    # Individual parameters with fixed-exponent allometry referenced to 70 kg
    cl <- exp(lcl + iov_cl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc) * (WT / 70)^e_wt_vc
    q <- exp(lq) * (WT / 70)^e_wt_q
    vp <- exp(lvp) * (WT / 70)^e_wt_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Unbound piperacillin concentration (mg/L)
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
