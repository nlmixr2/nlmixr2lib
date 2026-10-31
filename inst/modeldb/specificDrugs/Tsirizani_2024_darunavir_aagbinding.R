Tsirizani_2024_darunavir_aagbinding <- function() {
  description <- paste(
    "Semi-mechanistic two-compartment population PK model for oral",
    "darunavir boosted with ritonavir in African children with HIV on",
    "second-line ART (CHAPAS-4 PK sub-study; 59 children, 14.5-47.0 kg),",
    "with saturable instantaneous binding of darunavir to alpha-1-acid",
    "glycoprotein (AAG) in the central compartment. Bound darunavir follows",
    "Cbound = Bmax * Cfree / (Kd + Cfree), and the AAG binding capacity Bmax",
    "rises with plasma AAG through a power relationship (exponent 0.405,",
    "reference 0.7 g/L). Only free darunavir is cleared (intrinsic clearance",
    "CLint) and distributes to the peripheral compartment; the observation",
    "is total plasma darunavir. Absorption is a Savic transit chain feeding",
    "an absorption compartment drained at ka; clearances and volumes scale",
    "allometrically with weight (exponents fixed 0.75 and 1, reference",
    "26 kg). This is the alternative model reported in the supplement; the",
    "authors adopted the simpler direct AAG-on-clearance model",
    "(Tsirizani_2024_darunavir) for parsimony and usability."
  )
  reference <- paste(
    "Tsirizani L, Mohsenian Naghani S, Waalewijn H, Szubert A, Mulenga V,",
    "Chabala C, Bwakura-Dangarembizi M, Chitsamatanga M, Rutebarika DA,",
    "Musiime V, Kasozi M, Lugemwa A, Monkiewicz LN, McIlleron HM, Burger DM,",
    "Gibb DM, Denti P, Wasmann RE, Colbers A; the CHAPAS-4 trial team.",
    "Pharmacokinetics of once-daily darunavir/ritonavir in second-line",
    "treatment in African children with HIV. J Antimicrob Chemother.",
    "2024;79(11):2990-2998. doi:10.1093/jac/dkae319.",
    "Parameter values from Supplementary Table 1 and its Bmax footnote;",
    "structure from Supplementary Figure 2.",
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
        "Reference 26 kg (Supplementary Table 1 footnote: 'the typical values",
        "reported here refer to a child weighing 26 kg'). Exponents fixed at",
        "0.75 on CLint and Q and 1 on Vc and Vp (Methods, 'Population",
        "pharmacokinetic analysis'). Range 14.5-47.0 kg.",
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
        "Enters as a power effect on the darunavir binding capacity,",
        "Bmax = 20.7 * (AAG / 0.7)^0.405 (Supplementary Table 1 footnote #),",
        "centred at the Table 1 cohort median 0.7 g/L (range 0.3-1.9).",
        "Treated as time-fixed.",
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
        "inside model() into oc1 / oc2; the same occasion definition as",
        "Tsirizani_2024_darunavir (OCC = 1 earlier doses, OCC = 2 the",
        "witnessed PK-day dose and its observations; OCC = 2 for a single",
        "dose).",
        sep = " "
      ),
      source_name = "OCC"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "darunavir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "darunavir (free + AAG-bound)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "darunavir (free)", units = "mg", specimen = "plasma", verified = TRUE)
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
    co_medication = "NRTI backbone tenofovir alafenamide/emtricitabine (N = 34), abacavir/lamivudine (N = 12) or zidovudine/lamivudine (N = 13)",
    aag_range = "0.3-1.9 g/L (median 0.7)",
    notes = paste(
      "Same cohort and data as Tsirizani_2024_darunavir (Table 1). No unbound",
      "darunavir concentrations were measured; the binding parameters were",
      "informed by priors from Schalkwijk et al. (Discussion).",
      sep = " "
    )
  )

  ini({
    # Disposition of FREE darunavir, at 26 kg (Supplementary Table 1
    # footnote **). F is fixed to 1, so all are apparent (/F). Rows printed
    # with '(-)' have no SIR interval; the Discussion says the binding model
    # used priors from Schalkwijk et al. 'to inform our model and estimate
    # the parameters', so they are estimated, not fixed.
    lclint <- log(97.7); label("Apparent intrinsic clearance of free darunavir CLint/F at 26 kg (L/h)") # Supplementary Table 1 'Intrinsic Clearance - CLint (L/h) 97.7 (-)'
    lvc <- log(40.6); label("Apparent central volume Vc/F at 26 kg (L)")                                # Supplementary Table 1 'Central Volume of distribution - V (L) 40.6 (38.4-60.6)'
    lq <- log(62.7); label("Apparent intercompartmental clearance of free darunavir Q/F at 26 kg (L/h)") # Supplementary Table 1 'Intercompartmental Clearance - Q (L/h) 62.7 (-)'
    lvp <- log(822); label("Apparent peripheral volume Vp/F at 26 kg (L)")                              # Supplementary Table 1 'Peripheral Volume of distribution - V2 (L) 822 (-)'
    lfdepot <- fixed(log(1)); label("Oral bioavailability, reference (fraction)")                      # Supplementary Table 1 'Bioavailability - F 1 FIXED'

    lmtt <- log(0.849); label("Mean absorption transit time MTT (h)")                                   # Supplementary Table 1 'Mean absorption Transit Time - MTT (h) 0.849 (0.642-0.925)'
    lnn <- log(5.31); label("Number of absorption transit compartments NN (count)")                    # Supplementary Table 1 'Number of absorption transit compartments - NN (n) 5.31 (-)'
    lka <- log(0.519); label("First-order absorption rate constant ka (1/h)")                          # Supplementary Table 1 'First-order absorption rate constant - Ka (1/h) 0.519 (0.444-0.767)'

    # Saturable AAG binding in the central compartment.
    lkd <- log(1.31); label("Darunavir-AAG Michaelis-Menten dissociation constant Kd (mg/L)")         # Supplementary Table 1 'AAG Michaelis Menten dissociation constant - Kd (mg/L) 1.31 (-)'
    lbmax <- log(20.7); label("AAG maximum binding capacity for darunavir Bmax at AAG 0.7 g/L (mg/L)") # Supplementary Table 1 'AAG maximum binding capacity - Bmax (mg/L) 20.7 (-)'
    e_aag_bmax <- 0.405; label("Power exponent of AAG/0.7 on Bmax (unitless)")                        # Supplementary Table 1 'Exponent of AAG concentration effect on Bmax +0.405 (+0.219 - +0.731)'; footnote # 'Bmax_i = 20.7 * (AAG_i/0.7)^0.405'

    # Allometric exponents, fixed (Methods: 'a fixed exponent of 0.75 and 1').
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent of weight on CLint/F and Q/F (unitless)")   # Methods 'Population pharmacokinetic analysis'
    e_wt_vc_vp <- fixed(1); label("Allometric exponent of weight on Vc/F and Vp/F (unitless)")       # Methods 'Population pharmacokinetic analysis'

    # Variability percentages read as omega SD x 100 (see
    # Tsirizani_2024_darunavir), so omega^2 = (pct / 100)^2.
    etalclint ~ 0.123201 # Supplementary Table 1 'Between Subject Variability (%) Intrinsic Clearance 35.1 (25.6-41.3)' -> 0.351^2

    etaiov_fdepot_1 ~ 0.2209 # Supplementary Table 1 'Between Occasion Variability (%) Bioavailability 47.0 (31.8-54.9)' -> 0.470^2
    etaiov_fdepot_2 ~ fixed(0.2209) # same variance, occasion 2
    etaiov_mtt_1 ~ 0.389376 # Supplementary Table 1 'Mean absorption Transit Time 62.4 (48.6-73.7)' -> 0.624^2
    etaiov_mtt_2 ~ fixed(0.389376) # same variance, occasion 2
    etaiov_ka_1 ~ 0.149769 # Supplementary Table 1 'First-order absorption rate constant 38.7 (32.1-56.9)' -> 0.387^2
    etaiov_ka_2 ~ fixed(0.149769) # same variance, occasion 2

    propSd <- 0.127; label("Proportional residual error (fraction)")                                    # Supplementary Table 1 'Proportional error (%) 12.7 (10.5-15.8)'
    addSd <- 0.337; label("Additive residual error (mg/L)")                                             # Supplementary Table 1 'Additive error(mg/L) 0.337 (0.179-0.414)'
  })

  model({
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2

    clint <- exp(lclint + etalclint) * (WT / 26)^e_wt_cl_q
    vc <- exp(lvc) * (WT / 26)^e_wt_vc_vp
    q <- exp(lq) * (WT / 26)^e_wt_cl_q
    vp <- exp(lvp) * (WT / 26)^e_wt_vc_vp

    fdepot <- exp(lfdepot + iov_fdepot)
    mtt <- exp(lmtt + iov_mtt)
    nn <- exp(lnn)
    ka <- exp(lka + iov_ka)
    ktr <- (nn + 1) / mtt

    kd <- exp(lkd)
    bmax <- exp(lbmax) * (AAG / 0.7)^e_aag_bmax

    # Instantaneous binding: total = free + bmax * free / (kd + free). The
    # positive root of the resulting quadratic gives the free concentration
    # from the total held in central (Supplementary Figure 2).
    Ctot <- central / vc
    Cu <- 0.5 * ((Ctot - bmax - kd) + sqrt((Ctot - bmax - kd)^2 + 4 * kd * Ctot))

    # Only free drug is eliminated and exchanges with the peripheral
    # compartment, which holds free drug (Supplementary Figure 2 caption:
    # 'only the free drug is available for distribution and elimination').
    d/dt(depot) <- exp(
      log(fdepot * podo(depot)) + log(ktr) + nn * log(ktr * tad(depot)) -
        ktr * tad(depot) - lgamma(nn + 1)
    ) -
      ka * depot
    d/dt(central) <- ka * depot - clint * Cu - q * Cu + q * peripheral1 / vp
    d/dt(peripheral1) <- q * Cu - q * peripheral1 / vp
    f(depot) <- 0

    # Observation: total plasma darunavir (mg/L).
    Cc <- Ctot
    Cc ~ add(addSd) + prop(propSd)
  })
}
