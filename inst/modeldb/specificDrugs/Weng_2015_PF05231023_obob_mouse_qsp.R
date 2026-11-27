Weng_2015_PF05231023_obob_mouse_qsp <- function() {
  description <- paste(
    "Preclinical (mouse, ob/ob). QSP. Terminal-clipping PK and oral",
    "glucose tolerance test (OGTT) PD model for PF-05231023, a",
    "long-acting FGF21 analogue made of two FGF21 molecules conjugated",
    "to the CovX-2000 antibody scaffold. Each molecule carries four",
    "termini [C1 N1 C2 N2] (C- and N-terminus of the FGF21 on arm 1 and",
    "arm 2), each clipped independently by first-order proteolysis",
    "(C-terminal half-life 2.5 h, N-terminal 13 h) in the SC depot and",
    "in a single systemic compartment (4% of body weight), giving 16",
    "species per compartment (32 ODE states). SC absorption is",
    "first-order (half-life 48 h) for every species; the CovX-2000",
    "backbone is cleared from the systemic compartment with a 12-day",
    "half-life. The intact pool (species with at least one FGF21",
    "carrying both termini, counted once) and the N-only pool (no intact",
    "FGF21 but at least one intact N-terminus) drive a hyperbolic",
    "instantaneous signal S(t) = (I/KI + N/KN) / (1 + I/KI + N/KN), and",
    "the OGTT glucose-AUC reduction is A times the convolution of S(t)",
    "with exp(-t/gamma), carried as one ODE state. Deterministic: the",
    "paper fits mean responses by least squares and reports no",
    "between-animal or residual variability."
  )
  reference <- paste(
    "Weng Y, Chabot JR, Bernardo B, Yan Q, Zhu Y, Brenner MB, Vage C,",
    "Logan A, Calle R, Talukdar S. Pharmacokinetics (PK),",
    "pharmacodynamics (PD) and integrated PK/PD modeling of a novel long",
    "acting FGF21 clinical candidate PF-05231023 in diet-induced obese",
    "and leptin-deficient obese mice. PLoS One. 2015;10(3):e0119104.",
    "doi:10.1371/journal.pone.0119104."
  )
  vignette <- "Weng_2015_PF05231023_mouse_qsp"
  dosing <- c("depot_1111", "central_1111")
  units <- list(
    time = "h",
    dosing = "mg/kg",
    concentration = "nM"
  )

  # The 16 clipping species are named by the paper's own binary
  # [C1 N1 C2 N2] code (1 = terminus intact, 0 = clipped), in the SC
  # depot and in the systemic (central) compartment.
  paper_specific_compartment_pattern <- "^(depot|central)_[01]{4}$"
  paper_specific_compartments <- c("ogtt_memory")

  covariateData <- list()

  compartmentData <- list(
    depot_1111 = list(
      analyte = "PF-05231023 species [1111] (arm 1 C+/N+, arm 2 C+/N+)",
      units = "mg/kg",
      specimen = "administration site",
      verified = TRUE
    ),
    depot_1110 = list(
      analyte = "PF-05231023 species [1110] (arm 1 C+/N+, arm 2 C+/N+)",
      units = "mg/kg",
      specimen = "administration site",
      verified = TRUE
    ),
    depot_1101 = list(
      analyte = "PF-05231023 species [1101] (arm 1 C+/N+, arm 2 C+/N+)",
      units = "mg/kg",
      specimen = "administration site",
      verified = TRUE
    ),
    depot_1100 = list(
      analyte = "PF-05231023 species [1100] (arm 1 C+/N+, arm 2 C+/N+)",
      units = "mg/kg",
      specimen = "administration site",
      verified = TRUE
    ),
    depot_1011 = list(
      analyte = "PF-05231023 species [1011] (arm 1 C+/N+, arm 2 C+/N+)",
      units = "mg/kg",
      specimen = "administration site",
      verified = TRUE
    ),
    depot_1010 = list(
      analyte = "PF-05231023 species [1010] (arm 1 C+/N+, arm 2 C+/N+)",
      units = "mg/kg",
      specimen = "administration site",
      verified = TRUE
    ),
    depot_1001 = list(
      analyte = "PF-05231023 species [1001] (arm 1 C+/N+, arm 2 C+/N+)",
      units = "mg/kg",
      specimen = "administration site",
      verified = TRUE
    ),
    depot_1000 = list(
      analyte = "PF-05231023 species [1000] (arm 1 C+/N+, arm 2 C+/N+)",
      units = "mg/kg",
      specimen = "administration site",
      verified = TRUE
    ),
    depot_0111 = list(
      analyte = "PF-05231023 species [0111] (arm 1 C+/N+, arm 2 C+/N+)",
      units = "mg/kg",
      specimen = "administration site",
      verified = TRUE
    ),
    depot_0110 = list(
      analyte = "PF-05231023 species [0110] (arm 1 C+/N+, arm 2 C+/N+)",
      units = "mg/kg",
      specimen = "administration site",
      verified = TRUE
    ),
    depot_0101 = list(
      analyte = "PF-05231023 species [0101] (arm 1 C+/N+, arm 2 C+/N+)",
      units = "mg/kg",
      specimen = "administration site",
      verified = TRUE
    ),
    depot_0100 = list(
      analyte = "PF-05231023 species [0100] (arm 1 C+/N+, arm 2 C+/N+)",
      units = "mg/kg",
      specimen = "administration site",
      verified = TRUE
    ),
    depot_0011 = list(
      analyte = "PF-05231023 species [0011] (arm 1 C+/N+, arm 2 C+/N+)",
      units = "mg/kg",
      specimen = "administration site",
      verified = TRUE
    ),
    depot_0010 = list(
      analyte = "PF-05231023 species [0010] (arm 1 C+/N+, arm 2 C+/N+)",
      units = "mg/kg",
      specimen = "administration site",
      verified = TRUE
    ),
    depot_0001 = list(
      analyte = "PF-05231023 species [0001] (arm 1 C+/N+, arm 2 C+/N+)",
      units = "mg/kg",
      specimen = "administration site",
      verified = TRUE
    ),
    depot_0000 = list(
      analyte = "PF-05231023 species [0000] (arm 1 C+/N+, arm 2 C+/N+)",
      units = "mg/kg",
      specimen = "administration site",
      verified = TRUE
    ),
    central_1111 = list(
      analyte = "PF-05231023 species [1111] (arm 1 C+/N+, arm 2 C+/N+)",
      units = "mg/kg",
      specimen = "plasma",
      verified = TRUE
    ),
    central_1110 = list(
      analyte = "PF-05231023 species [1110] (arm 1 C+/N+, arm 2 C+/N+)",
      units = "mg/kg",
      specimen = "plasma",
      verified = TRUE
    ),
    central_1101 = list(
      analyte = "PF-05231023 species [1101] (arm 1 C+/N+, arm 2 C+/N+)",
      units = "mg/kg",
      specimen = "plasma",
      verified = TRUE
    ),
    central_1100 = list(
      analyte = "PF-05231023 species [1100] (arm 1 C+/N+, arm 2 C+/N+)",
      units = "mg/kg",
      specimen = "plasma",
      verified = TRUE
    ),
    central_1011 = list(
      analyte = "PF-05231023 species [1011] (arm 1 C+/N+, arm 2 C+/N+)",
      units = "mg/kg",
      specimen = "plasma",
      verified = TRUE
    ),
    central_1010 = list(
      analyte = "PF-05231023 species [1010] (arm 1 C+/N+, arm 2 C+/N+)",
      units = "mg/kg",
      specimen = "plasma",
      verified = TRUE
    ),
    central_1001 = list(
      analyte = "PF-05231023 species [1001] (arm 1 C+/N+, arm 2 C+/N+)",
      units = "mg/kg",
      specimen = "plasma",
      verified = TRUE
    ),
    central_1000 = list(
      analyte = "PF-05231023 species [1000] (arm 1 C+/N+, arm 2 C+/N+)",
      units = "mg/kg",
      specimen = "plasma",
      verified = TRUE
    ),
    central_0111 = list(
      analyte = "PF-05231023 species [0111] (arm 1 C+/N+, arm 2 C+/N+)",
      units = "mg/kg",
      specimen = "plasma",
      verified = TRUE
    ),
    central_0110 = list(
      analyte = "PF-05231023 species [0110] (arm 1 C+/N+, arm 2 C+/N+)",
      units = "mg/kg",
      specimen = "plasma",
      verified = TRUE
    ),
    central_0101 = list(
      analyte = "PF-05231023 species [0101] (arm 1 C+/N+, arm 2 C+/N+)",
      units = "mg/kg",
      specimen = "plasma",
      verified = TRUE
    ),
    central_0100 = list(
      analyte = "PF-05231023 species [0100] (arm 1 C+/N+, arm 2 C+/N+)",
      units = "mg/kg",
      specimen = "plasma",
      verified = TRUE
    ),
    central_0011 = list(
      analyte = "PF-05231023 species [0011] (arm 1 C+/N+, arm 2 C+/N+)",
      units = "mg/kg",
      specimen = "plasma",
      verified = TRUE
    ),
    central_0010 = list(
      analyte = "PF-05231023 species [0010] (arm 1 C+/N+, arm 2 C+/N+)",
      units = "mg/kg",
      specimen = "plasma",
      verified = TRUE
    ),
    central_0001 = list(
      analyte = "PF-05231023 species [0001] (arm 1 C+/N+, arm 2 C+/N+)",
      units = "mg/kg",
      specimen = "plasma",
      verified = TRUE
    ),
    central_0000 = list(
      analyte = "PF-05231023 species [0000] (arm 1 C+/N+, arm 2 C+/N+)",
      units = "mg/kg",
      specimen = "plasma",
      verified = TRUE
    ),
    ogtt_memory = list(
      analyte = "time-integrated FGF21 signal S(t) convolved with exp(-t/gamma)",
      units = "h",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "mouse (ob/ob, leptin-deficient)",
    n_subjects = NA_integer_,
    n_studies = NA_integer_,
    age_range = "9-10 weeks (male)",
    sex_female_pct = 0,
    disease_state = "Leptin-deficient obese (ob/ob) mice",
    dose_range = "PK: 3 or 10 mg/kg single IV or SC; OGTT: 0.03-10 mg/kg single IV or SC",
    notes = paste(
      "PK in male 9-10 week old ob/ob mice, single IV or SC dose; intact",
      "C-terminus (CT) and intact N-terminus (NT) measured by two",
      "sandwich ELISAs (Methods; Table 1 and Table 2; n = 3-4 per",
      "group). OGTT model fitted to the OGTT glucose-AUC reduction on",
      "days 3, 4, 5 and 6 after single 3 or 10 mg/kg SC doses",
      "(Fig 6a; at least 8 mice per group per day, Fig 3) and checked",
      "against the day-3 dose response to 0.03-10 mg/kg IV or SC",
      "(Fig 2, Fig 6b; 8 mice per group). Group sizes are per time point;",
      "the total number of animals is not reported."
    )
  )

  ini({
    # PK (Methods, 'PK/PD modeling for OGTT'; all fixed from measurements, none fitted)
    lvc <- fixed(log(40)); label("Systemic compartment volume (mL/kg)") # Methods: 'single systemic compartment, representing 4% of the body weight' = 40 mL/kg
    lka <- fixed(log(log(2) / 48)); label("First-order SC absorption rate constant, all species (1/h)") # Methods: 'absorption half-life of 48 hours'
    thalf_ct <- fixed(2.5); label("C-terminus clipping half-life (h)") # Methods: '2.5 and 13 hours for the C and N termini respectively'
    thalf_nt <- fixed(13); label("N-terminus clipping half-life (h)") # Methods: '2.5 and 13 hours for the C and N termini respectively'
    thalf_bb <- fixed(12 * 24); label("CovX-2000 backbone elimination half-life (h)") # Methods: 'observed half-life of the CovX-2000 backbone (~12 days)'; Eq 1 ln(2)/12 days
    mw <- fixed(192000); label("Molecular weight of PF-05231023 (g/mol)") # Not printed; read from Fig 5a: 10 mg/kg IV gives ~1300 nM at t = 0 in 40 mL/kg (0.25 g/L / 1.3e-6 M = 192 kDa)

    # PD: instantaneous signalling (Eq 2) and OGTT convolution (Eq 3)
    ki <- fixed(0.5); label("Half-maximal signalling concentration of the intact pool, KI (nM)") # Methods: 'KI was set to 0.5 nM, reflecting the in vitro pERK response EC50 in hASC cells'
    kn <- 22; label("Half-maximal signalling concentration of the N-only pool, KN (nM)") # Results: 'KN was optimized at 22 nM'
    gamma <- 3 * 24; label("Decay time constant of signal persistence, gamma (h)") # Results: 'gamma was optimized at 3 days'
    ogtt_a <- fixed(812.74); label("OGTT normalisation A ((mg/dL)*min glucose-AUC reduction per h of integrated signal)") # Not printed; Methods/Results: 'A was set by matching the day 3 response to the 10 mg/kg dose' -- solved so the model equals the Fig 6a day-3 10 mg/kg SC value (29900 (mg/dL)*min, digitised)
  })

  model({
    vc <- exp(lvc)
    ka <- exp(lka)
    kct <- log(2) / thalf_ct
    knt <- log(2) / thalf_nt
    kbb <- log(2) / thalf_bb

    # Eq 1 generalised to all 16 species: SC absorption (ka) from the
    # depot, clipping of each intact terminus (kct for C1 and C2, knt for
    # N1 and N2) in both compartments, and backbone clearance (kbb) in the
    # systemic compartment only. IV doses go to central_1111, SC doses to
    # depot_1111.
    d/dt(depot_1111) <- -ka * depot_1111 - kct * depot_1111 - knt * depot_1111 - kct * depot_1111 - knt * depot_1111
    d/dt(depot_1110) <- -ka * depot_1110 - kct * depot_1110 - knt * depot_1110 - kct * depot_1110 + knt * depot_1111
    d/dt(depot_1101) <- -ka * depot_1101 - kct * depot_1101 - knt * depot_1101 - knt * depot_1101 + kct * depot_1111
    d/dt(depot_1100) <- -ka * depot_1100 - kct * depot_1100 - knt * depot_1100 + kct * depot_1110 + knt * depot_1101
    d/dt(depot_1011) <- -ka * depot_1011 - kct * depot_1011 - kct * depot_1011 - knt * depot_1011 + knt * depot_1111
    d/dt(depot_1010) <- -ka * depot_1010 - kct * depot_1010 - kct * depot_1010 + knt * depot_1110 + knt * depot_1011
    d/dt(depot_1001) <- -ka * depot_1001 - kct * depot_1001 - knt * depot_1001 + knt * depot_1101 + kct * depot_1011
    d/dt(depot_1000) <- -ka * depot_1000 - kct * depot_1000 + knt * depot_1100 + kct * depot_1010 + knt * depot_1001
    d/dt(depot_0111) <- -ka * depot_0111 - knt * depot_0111 - kct * depot_0111 - knt * depot_0111 + kct * depot_1111
    d/dt(depot_0110) <- -ka * depot_0110 - knt * depot_0110 - kct * depot_0110 + kct * depot_1110 + knt * depot_0111
    d/dt(depot_0101) <- -ka * depot_0101 - knt * depot_0101 - knt * depot_0101 + kct * depot_1101 + kct * depot_0111
    d/dt(depot_0100) <- -ka * depot_0100 - knt * depot_0100 + kct * depot_1100 + kct * depot_0110 + knt * depot_0101
    d/dt(depot_0011) <- -ka * depot_0011 - kct * depot_0011 - knt * depot_0011 + kct * depot_1011 + knt * depot_0111
    d/dt(depot_0010) <- -ka * depot_0010 - kct * depot_0010 + kct * depot_1010 + knt * depot_0110 + knt * depot_0011
    d/dt(depot_0001) <- -ka * depot_0001 - knt * depot_0001 + kct * depot_1001 + knt * depot_0101 + kct * depot_0011
    d/dt(depot_0000) <- -ka * depot_0000 + kct * depot_1000 + knt * depot_0100 + kct * depot_0010 + knt * depot_0001

    d/dt(central_1111) <- ka * depot_1111 - kbb * central_1111 - kct * central_1111 - knt * central_1111 - kct * central_1111 - knt * central_1111
    d/dt(central_1110) <- ka * depot_1110 - kbb * central_1110 - kct * central_1110 - knt * central_1110 - kct * central_1110 + knt * central_1111
    d/dt(central_1101) <- ka * depot_1101 - kbb * central_1101 - kct * central_1101 - knt * central_1101 - knt * central_1101 + kct * central_1111
    d/dt(central_1100) <- ka * depot_1100 - kbb * central_1100 - kct * central_1100 - knt * central_1100 + kct * central_1110 + knt * central_1101
    d/dt(central_1011) <- ka * depot_1011 - kbb * central_1011 - kct * central_1011 - kct * central_1011 - knt * central_1011 + knt * central_1111
    d/dt(central_1010) <- ka * depot_1010 - kbb * central_1010 - kct * central_1010 - kct * central_1010 + knt * central_1110 + knt * central_1011
    d/dt(central_1001) <- ka * depot_1001 - kbb * central_1001 - kct * central_1001 - knt * central_1001 + knt * central_1101 + kct * central_1011
    d/dt(central_1000) <- ka * depot_1000 - kbb * central_1000 - kct * central_1000 + knt * central_1100 + kct * central_1010 + knt * central_1001
    d/dt(central_0111) <- ka * depot_0111 - kbb * central_0111 - knt * central_0111 - kct * central_0111 - knt * central_0111 + kct * central_1111
    d/dt(central_0110) <- ka * depot_0110 - kbb * central_0110 - knt * central_0110 - kct * central_0110 + kct * central_1110 + knt * central_0111
    d/dt(central_0101) <- ka * depot_0101 - kbb * central_0101 - knt * central_0101 - knt * central_0101 + kct * central_1101 + kct * central_0111
    d/dt(central_0100) <- ka * depot_0100 - kbb * central_0100 - knt * central_0100 + kct * central_1100 + kct * central_0110 + knt * central_0101
    d/dt(central_0011) <- ka * depot_0011 - kbb * central_0011 - kct * central_0011 - knt * central_0011 + kct * central_1011 + knt * central_0111
    d/dt(central_0010) <- ka * depot_0010 - kbb * central_0010 - kct * central_0010 + kct * central_1010 + knt * central_0110 + knt * central_0011
    d/dt(central_0001) <- ka * depot_0001 - kbb * central_0001 - knt * central_0001 + kct * central_1001 + knt * central_0101 + kct * central_0011
    d/dt(central_0000) <- ka * depot_0000 - kbb * central_0000 + kct * central_1000 + knt * central_0100 + kct * central_0010 + knt * central_0001

    amt_intact <- central_1111 + central_1110 + central_1101 + central_1100 + central_1011 + central_0111 + central_0011
    amt_nonly <- central_1001 + central_0110 + central_0101 + central_0100 + central_0001

    # Pools in nM: amount (mg/kg) / volume (mL/kg) = g/L, / mw -> mol/L
    Cintact <- amt_intact / vc * 1e9 / mw
    Cnonly <- amt_nonly / vc * 1e9 / mw
    # What the two ELISAs measure: every molecule with at least one intact
    # C-terminus (CT assay) or at least one intact N-terminus (NT assay)
    amt_ct <- central_1111 + central_1110 + central_1101 + central_1100 + central_1011 + central_1010 + central_1001 + central_1000 + central_0111 + central_0110 + central_0011 + central_0010
    amt_nt <- central_1111 + central_1110 + central_1101 + central_1100 + central_1011 + central_1001 + central_0111 + central_0110 + central_0101 + central_0100 + central_0011 + central_0001
    Cct <- amt_ct / vc * 1e9 / mw
    Cnt <- amt_nt / vc * 1e9 / mw

    # Eq 2: instantaneous signalling by the intact and N-only pools
    sig <- (Cintact / ki + Cnonly / kn) / (1 + Cintact / ki + Cnonly / kn)
    # Eq 3: OGTT(t) = A * integral_0^t S(tau) exp(-(t - tau)/gamma) dtau
    d/dt(ogtt_memory) <- sig - ogtt_memory / gamma
    ogtt <- ogtt_a * ogtt_memory
  })
}
