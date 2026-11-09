Tsirizani_2024_darunavir <- function() {
  description <- paste(
    "Two-compartment population PK model for oral darunavir boosted with",
    "ritonavir (600/100 or 800/100 mg once daily) in African children with",
    "HIV on second-line ART (CHAPAS-4 trial PK sub-study; 59 children aged",
    "3.8-14.7 years and 14.5-47.0 kg from Zambia, Uganda and Zimbabwe).",
    "Absorption is a Savic transit-compartment chain (mean transit time MTT,",
    "non-integer number of transit compartments NN) feeding an absorption",
    "compartment drained first-order at ka. Clearance and both volumes are",
    "allometrically scaled to total body weight with exponents fixed at 0.75",
    "and 1, referenced to 26 kg. Plasma alpha-1-acid glycoprotein (AAG)",
    "lowers apparent clearance exponentially (slope -0.582 per g/L around",
    "0.66 g/L). Between-subject variability on clearance only;",
    "between-occasion variability on bioavailability, MTT and ka. No effect",
    "of the NRTI backbone (tenofovir alafenamide/emtricitabine,",
    "abacavir/lamivudine or zidovudine/lamivudine), age or ritonavir",
    "exposure was retained. All clearances and volumes are apparent (/F)."
  )
  reference <- paste(
    "Tsirizani L, Mohsenian Naghani S, Waalewijn H, Szubert A, Mulenga V,",
    "Chabala C, Bwakura-Dangarembizi M, Chitsamatanga M, Rutebarika DA,",
    "Musiime V, Kasozi M, Lugemwa A, Monkiewicz LN, McIlleron HM, Burger DM,",
    "Gibb DM, Denti P, Wasmann RE, Colbers A; the CHAPAS-4 trial team.",
    "Pharmacokinetics of once-daily darunavir/ritonavir in second-line",
    "treatment in African children with HIV. J Antimicrob Chemother.",
    "2024;79(11):2990-2998. doi:10.1093/jac/dkae319.",
    "Parameter values from Table 3 and its footnote c (the clearance-AAG",
    "equation).",
    sep = " "
  )
  vignette <- "Tsirizani_2024_darunavir"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Total body weight, the body-size descriptor for allometric scaling of clearance and volume",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Reference 26 kg, the cohort median (Table 1; Table 3 footnote b 'the",
        "typical values reported here refer to a child weighing 26 kg').",
        "Exponents fixed at 0.75 on CL/F and Q/F and 1 on Vc/F and Vp/F",
        "(Methods, 'Population pharmacokinetic analysis'). Total body weight",
        "was preferred over fat-free mass (Results, dOFV = 2.62). Range",
        "14.5-47.0 kg.",
        sep = " "
      ),
      source_name = "WT"
    ),
    AAG = list(
      description = "Plasma alpha-1-acid glycoprotein concentration",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Centred at 0.66 g/L in the clearance equation of Table 3 footnote c,",
        "CL = 8.05 * (WT/26)^0.75 * exp(-0.582 * (AAG - 0.66)). The cohort",
        "median printed in Table 1 is 0.7 g/L (range 0.3-1.9); 0.66 is",
        "presumably the unrounded median and is used as printed. Measured by",
        "a validated immunoturbidimetric assay (LLOQ 0.10 g/L; Methods,",
        "'Procedures'). Treated as time-fixed.",
        sep = " "
      ),
      source_name = "AAG"
    ),
    OCC = list(
      description = "Dosing-occasion index for between-occasion variability on the absorption parameters",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Two occasions with a shared BOV variance per parameter, decomposed",
        "inside model() into binary indicators oc1 / oc2. The paper states",
        "only that BOV was included on the absorption parameters (Methods),",
        "with one intensive PK day per child; the occasion definition is taken",
        "from the same group's CHAPAS-4 ritonavir analysis",
        "(Tsirizani_2025_ritonavir.R, whose supplemental control stream sets",
        "OCC = 1 for the unwitnessed dose preceding the sampling window and",
        "OCC = 2 for the witnessed dose on the PK day). Pass OCC = 1 on the",
        "earlier doses and OCC = 2 on the witnessed dose and its observations;",
        "for a single-dose simulation pass OCC = 2.",
        sep = " "
      ),
      source_name = "OCC"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "darunavir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "darunavir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "darunavir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 59,
    n_studies = 1,
    age_range = "3.8-14.7 years",
    age_median = "10.9 years",
    weight_range = "14.5-47.0 kg",
    weight_median = "26.0 kg",
    height_range = "97.8-161 cm",
    height_median = "132 cm",
    sex_female_pct = 55.9,
    race_ethnicity = "not reported by category; all participants enrolled in Zambia, Uganda and Zimbabwe",
    disease_state = "HIV-1 infection failing first-line abacavir- or zidovudine-based ART by WHO virological, CD4 or clinical criteria, switched to second-line darunavir/ritonavir",
    dose_range = "darunavir/ritonavir 600/100 mg once daily (14-24.9 kg) or 800/100 mg once daily (>= 25 kg); 17.0-41.4 mg/kg darunavir",
    regions = "Zambia, Uganda, Zimbabwe",
    co_medication = "NRTI backbone tenofovir alafenamide/emtricitabine (N = 34), abacavir/lamivudine (N = 12) or zidovudine/lamivudine (N = 13); no NRTI effect on darunavir PK",
    aag_range = "0.3-1.9 g/L (median 0.7)",
    notes = paste(
      "Baseline characteristics in Table 1. Darunavir/ritonavir arm of the",
      "CHAPAS-4 intensive PK sub-study: steady-state sampling at week 6",
      "pre-dose and 1, 2, 4, 6, 8, 12 and 24 h after a witnessed dose taken",
      "with breakfast (plus 0.5 h in the TAF arm); 491 darunavir",
      "concentrations, one below the 0.105 mg/L LLOQ.",
      sep = " "
    )
  )

  ini({
    # Disposition. Table 3 typical values refer to a 26 kg child (footnote b)
    # at AAG 0.66 g/L (footnote c). F is fixed to 1, so all clearances and
    # volumes are apparent (/F).
    lcl <- log(8.05); label("Apparent clearance CL/F at 26 kg and AAG 0.66 g/L (L/h)")             # Table 3 'Clearance: CL (L/h) 8.05 (6.85-8.81)'
    lvc <- log(46.3); label("Apparent central volume Vc/F at 26 kg (L)")                            # Table 3 'Central volume of distribution: V (L) 46.3 (37.2-56.0)'
    lq <- log(7.79); label("Apparent intercompartmental clearance Q/F at 26 kg (L/h)")              # Table 3 'Intercompartmental clearance: Q (L/h) 7.79 (6.34-9.05)'
    lvp <- log(82.6); label("Apparent peripheral volume Vp/F at 26 kg (L)")                          # Table 3 'Peripheral volume of distribution: V2 (L) 82.6 (65.1-101.5)'
    lfdepot <- fixed(log(1)); label("Oral bioavailability, reference (fraction)")                  # Table 3 'Bioavailability: F 1 FIXED'

    # Absorption: Savic transit chain into an absorption compartment.
    lmtt <- log(0.763); label("Mean absorption transit time MTT (h)")                                # Table 3 'Mean absorption transit time: MTT (h) 0.763 (0.655-0.910)'
    lnn <- log(11.6); label("Number of absorption transit compartments NN (count)")                 # Table 3 'Number of absorption transit compartments: NN (n) 11.6 (5.03-23.8)'
    lka <- log(0.551); label("First-order absorption rate constant ka (1/h)")                       # Table 3 'First-order absorption rate constant: Ka (1/h) 0.551 (0.437-0.696)'

    # Allometric exponents, fixed (Methods: 'a fixed exponent of 0.75 and 1').
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent of weight on CL/F and Q/F (unitless)")   # Methods 'Population pharmacokinetic analysis'; Table 3 footnote c '(WTi/26)0.75'
    e_wt_vc_vp <- fixed(1); label("Allometric exponent of weight on Vc/F and Vp/F (unitless)")    # Methods 'Population pharmacokinetic analysis'

    # Exponential AAG effect on clearance (Table 3 footnote c).
    e_aag_cl <- -0.582; label("Exponential slope of AAG on CL/F (per g/L)")                         # Table 3 'Slope of AAG concentration effect on CL -0.582 (-0.236 to -0.667)'; footnote c 'e^(-0.582*(AAGi-0.66))'

    # Variability. Table 3 percentages are read as the log-scale omega SD
    # times 100, the convention this group's CHAPAS-4 ritonavir analysis
    # (Tsirizani_2025_ritonavir.R) verified against its control stream, so
    # omega^2 = (pct / 100)^2.
    etalcl ~ 0.104976 # Table 3 'BSV (%) Clearance 32.4 (25.5-40.8)' -> 0.324^2

    # BOV: one variance per parameter shared across two occasions; the
    # occasion-2 eta is fixed to the same value (NONMEM BLOCK(1) SAME).
    etaiov_fdepot_1 ~ 0.197136 # Table 3 'BOV (%) Bioavailability 44.4 (31.7-55.5)' -> 0.444^2
    etaiov_fdepot_2 ~ fixed(0.197136) # same variance, occasion 2
    etaiov_mtt_1 ~ 0.369664 # Table 3 'BOV (%) Mean absorption transit time 60.8 (48.8-79.7)' -> 0.608^2
    etaiov_mtt_2 ~ fixed(0.369664) # same variance, occasion 2
    etaiov_ka_1 ~ 0.188356 # Table 3 'BOV (%) First-order absorption rate constant 43.4 (31.8-58.5)' -> 0.434^2
    etaiov_ka_2 ~ fixed(0.188356) # same variance, occasion 2

    # Combined additive and proportional residual error.
    propSd <- 0.125; label("Proportional residual error (fraction)")                                 # Table 3 'Proportional error (%) 12.5 (10.6-15.0)'
    addSd <- 0.291; label("Additive residual error (mg/L)")                                          # Table 3 'Additive error (mg/L) 0.291 (0.202-0.382)'
  })

  model({
    # Occasion indicators and per-occasion BOV.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2

    # Individual parameters (Table 3 footnote c for CL).
    cl <- exp(lcl + etalcl) * (WT / 26)^e_wt_cl_q * exp(e_aag_cl * (AAG - 0.66))
    vc <- exp(lvc) * (WT / 26)^e_wt_vc_vp
    q <- exp(lq) * (WT / 26)^e_wt_cl_q
    vp <- exp(lvp) * (WT / 26)^e_wt_vc_vp

    fdepot <- exp(lfdepot + iov_fdepot)
    mtt <- exp(lmtt + iov_mtt)
    nn <- exp(lnn)
    ka <- exp(lka + iov_ka)

    # Savic 2007 transit chain: ktr = (NN + 1) / MTT.
    ktr <- (nn + 1) / mtt

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # The dose record targets depot; f(depot) <- 0 suppresses the bolus and
    # the Savic gamma-density input below delivers the bioavailable dose to
    # the absorption compartment. This is the closed form transit(nn, mtt,
    # fdepot) expands to, written out because podo(depot) is not reduced by
    # f(depot) whereas transit()'s internal dose lookup is.
    d/dt(depot) <- exp(
      log(fdepot * podo(depot)) + log(ktr) + nn * log(ktr * tad(depot)) -
        ktr * tad(depot) - lgamma(nn + 1)
    ) -
      ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    f(depot) <- 0

    # Total plasma darunavir; dose mg / volume L = mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
