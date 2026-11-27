Rabie_2020_abacavir <- function() {
  description <- paste(
    "Two-compartment population PK model for oral abacavir in HIV-infected",
    "South African children weighing 3-15 kg who were cotreated for",
    "tuberculosis with rifampin and superboosted lopinavir/ritonavir 4:4",
    "(Rabie 2020). Absorption is a Savic analytical transit chain (mean",
    "transit time 3.6 min, 1.03 transit compartments) feeding a first-order",
    "depot; disposition is allometrically scaled on total body weight",
    "(exponents fixed to 0.75 on CL and Q and 1 on Vc and Vp, reference",
    "9.4 kg), and clearance matures with postmenstrual age through a",
    "sigmoidal Hill function. Rifampin plus LPV/r 4:4 cotreatment lowers",
    "bioavailability by 36% relative to standard LPV/r 4:1, and the",
    "unobserved evening dose before each PK visit is absorbed after an",
    "extra 2.62 h delay. Random effects: between-subject variability on CL,",
    "Vp and Q; between-visit variability on CL and on the evening-dose",
    "delay; between-occasion (per-dose) variability on ka, bioavailability",
    "and mean transit time."
  )
  reference <- paste(
    "Rabie H, Tikiso T, Lee J, Fairlie L, Strehlau R, Bobat R, Liberty A,",
    "McIlleron H, Andrieux-Meyer I, Cotton M, Lallemant M, Denti P.",
    "Abacavir Exposure in Children Cotreated for Tuberculosis",
    "with Rifampin and Superboosted Lopinavir-Ritonavir. Antimicrob Agents",
    "Chemother. 2020;64(5):e01923-19. doi:10.1128/AAC.01923-19."
  )
  vignette <- "Rabie_2020_abacavir"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "abacavir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "abacavir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "abacavir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric size descriptor on CL, Q (exponent 0.75) and Vc, Vp",
        "(exponent 1), centred on 9.4 kg, the median weight in the data set",
        "(Rabie 2020 Table 2 footnote c and Methods). Fat-free mass was",
        "tested as an alternative and did not improve the fit."
      ),
      source_name = "TBW"
    ),
    PAGE = list(
      description = "Postmenstrual age",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives the sigmoidal maturation of CL (Rabie 2020 Eq. 1). No",
        "gestational ages were available, so the authors set postmenstrual",
        "age to postnatal age + 9 months for every child; supply",
        "PAGE = postnatal_age_months + 9 unless a real gestational age is",
        "known."
      ),
      source_name = "PMAGE"
    ),
    CONMED_RIF_LPVR4 = list(
      description = "Rifampin-based antituberculosis treatment with superboosted lopinavir/ritonavir 4:4 indicator (1 = rifampin + LPV/r 4:4, 0 = standard LPV/r 4:1 without rifampin)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (standard LPV/r 4:1 without rifampin)",
      notes = paste(
        "1 at PK visits 1 and 2 (intensive and continuation phases of TB",
        "treatment), 0 at PK visit 3 (one month after TB treatment ended,",
        "back on LPV/r 4:1). Rifampin and the extra ritonavir always",
        "coincided in the data, so the effect cannot be attributed to one",
        "drug (Rabie 2020 Results). Multiplicative -36.0% on",
        "bioavailability."
      ),
      source_name = "RIF"
    ),
    DOSETIME_EVENING = list(
      description = "Evening (night) dose indicator (1 = the unobserved evening dose given at home on the evening before a PK visit, 0 = any other dose)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the directly observed morning dose at the clinic)",
      notes = paste(
        "Per dose record. In the source data every PK visit carries one",
        "evening dose given at home the night before (time recorded by the",
        "caregiver) and one observed morning dose; the model adds an",
        "absorption delay of 2.62 h to the evening dose only, which the",
        "authors introduced because the observed morning predose",
        "concentrations exceeded the 10-h concentrations (Rabie 2020",
        "Results and Table 3). The paper attributes the delay to late",
        "dosing, food or diurnal variation without distinguishing them.",
        "Whether the delay should be applied to every evening dose of a",
        "long twice-daily simulation is an extrapolation of the fit."
      ),
      source_name = "night dose"
    ),
    OCC = list(
      description = "Integer dosing-occasion index for the between-occasion and between-visit random effects",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Rabie 2020 Methods: 'Each dose was treated as a separate occasion,",
        "while consecutive evening and morning doses were grouped within the",
        "same visit.' Coded here as six occasions over the three PK visits:",
        "OCC = 2k - 1 is the evening dose and OCC = 2k the observed morning",
        "dose of PK visit k (k = 1, 2, 3). Between-occasion etas on ka, F",
        "and MTT switch with OCC; the between-visit etas on CL and on the",
        "evening-dose delay are shared by occasions 2k - 1 and 2k, so no",
        "separate visit column is needed. Any value outside 1-6 (e.g.",
        "OCC = 0) switches all between-occasion and between-visit",
        "variability off, leaving only between-subject variability."
      ),
      source_name = "OCC"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 87L,
    n_studies = 1L,
    n_observations = 1344L,
    age_range = "Median (IQR) 19 (4 to 64) months at PK visit 1; no child younger than 3 months",
    age_median = "19 months (PK visit 1)",
    weight_range = "Enrolment criterion 3 to 15 kg; median (IQR) 8.7 (3.9 to 14.9) kg at PK visit 1",
    weight_median = "8.7 kg (PK visit 1); 9.4 kg over the whole data set",
    sex_female_pct = 57,
    race_ethnicity = "Black African (South Africa); not reported by category",
    disease_state = "HIV-infected children with clinician-diagnosed tuberculosis on rifampin-based TB treatment and lopinavir/ritonavir-based ART; median weight-for-age Z-score -2.43 at PK visit 1",
    dose_range = "Abacavir oral solution 20 mg/mL, 40-160 mg twice daily by weight band (Rabie 2020 Table 4), median 9.5-9.7 mg/kg per dose",
    regions = "South Africa (5 sites in Western Cape, Gauteng and KwaZulu-Natal)",
    co_medication = "Rifampin, isoniazid, pyrazinamide (with or without ethambutol / ethionamide); lopinavir/ritonavir 4:1 plus extra ritonavir to 4:4 during rifampin and for 2 weeks after; lamivudine",
    notes = paste(
      "Rabie 2020 Table 1 and Results. Intensive sampling (predose and 1, 2,",
      "4, 6, 10 h) at three visits: PK1 (n = 86) in the second month and PK2",
      "(n = 74) in the last month of TB treatment, both on rifampin + LPV/r",
      "4:4, and PK3 (n = 71) 4-6 weeks after LPV/r 4:1 was reinstated.",
      "Percent female is from PK visit 1 (86 analysed, 37 male = 43%).",
      "Lower limit of quantification 0.0243 mg/L; 13-17% of samples BLQ,",
      "handled with the M6 method. Registered as NCT02348177."
    )
  )

  ini({
    # Structural typical values refer to a 9.4 kg child on standard LPV/r
    # 4:1 (Rabie 2020 Table 2 footnote c). All are apparent oral values.
    lcl <- log(9.67)
    label("Apparent clearance CL/F for a 9.4 kg fully mature child (L/h)") # Table 2 Clearance 9.67 (8.27-10.6) liters/h
    lvc <- log(8.76)
    label("Apparent central volume of distribution Vc/F for a 9.4 kg child (L)") # Table 2 Central vol of distribution 8.76 (6.99-9.94) liters
    lvp <- log(3.32)
    label("Apparent peripheral volume of distribution Vp/F for a 9.4 kg child (L)") # Table 2 Peripheral vol of distribution 3.32 (2.82-3.98) liters
    lq <- log(1.35)
    label("Apparent intercompartmental clearance Q/F for a 9.4 kg child (L/h)") # Table 2 Intercompartmental clearance 1.35 (1.09-1.67) liters/h
    lka <- log(2.22)
    label("First-order absorption rate constant from depot to central (1/h)") # Table 2 Absorption rate constant ka 2.22 (1.82-2.80) 1/h
    lmtt <- log(3.60 / 60)
    label("Mean transit time of the absorption transit chain (h)") # Table 2 Absorption mean transit time 3.60 (2.64-6.54) min; /60 to h
    lnn <- log(1.03)
    label("Number of absorption transit compartments (unitless)") # Table 2 No. of transit compartments 1.03 (0.901-1.20)
    lfdepot <- fixed(log(1))
    label("Oral bioavailability on standard LPV/r 4:1 (fraction)") # Table 2 Bioavailability '1 (fixed)'
    ltlag <- log(2.62)
    label("Extra absorption delay applied to the evening dose before a PK visit (h)") # Table 2 Delay in absorption for night dose 2.62 (2.23-2.90) h

    # Allometric exponents on total body weight (Rabie 2020 Methods: 'exponents
    # fixed to 3/4 for clearance parameters and 1 for volumes of distribution')
    e_wt_cl_q <- fixed(0.75)
    label("Allometric exponent of body weight on CL and Q (unitless)") # Methods, Statistical analysis (i)
    e_wt_vc_vp <- fixed(1)
    label("Allometric exponent of body weight on Vc and Vp (unitless)") # Methods, Statistical analysis (i)

    # Maturation of CL on postmenstrual age (Rabie 2020 Eq. 1)
    ltm50_cl <- log(10.7)
    label("Postmenstrual age at 50% maturation of CL (months)") # Table 2 PMAGE50 10.7 (10.5-11.0) months from conception
    lhill_cl <- log(4.35)
    label("Hill coefficient of the CL maturation function (unitless)") # Table 2 gamma 4.35 (3.23-5.50)

    # Cotreatment with rifampin and superboosted LPV/r 4:4
    e_rif_lpvr4_fdepot <- -0.360
    label("Rifampin + LPV/r 4:4 effect on bioavailability relative to LPV/r 4:1 (fraction)") # Table 2 'LPV/r-4:4 + rifampin effect on bioavailability' -36.0 (-38.7 to -32.9) percent

    # Between-subject variability. Table 2 footnote a: random effects
    # 'assumed to be lognormally distributed and are reported as the
    # coefficient of variation', so omega^2 = log(1 + CV^2).
    etalcl ~ log(1 + 0.144^2) # Table 2 Clearance BSV 14.4% CV (10.9-18.5)
    etalvp ~ log(1 + 0.290^2) # Table 2 Peripheral volume BSV 29.0% CV (19.2-33.6)
    etalq ~ log(1 + 0.135^2) # Table 2 Intercompartmental clearance BSV 13.5% CV (12.4-19.9)

    # Between-visit variability (three PK visits), one shared variance each.
    etabvv_cl_1 ~ log(1 + 0.184^2) # Table 2 Clearance BVV 18.4% CV (14.4-21.4)
    etabvv_cl_2 ~ fixed(log(1 + 0.184^2)) # same variance as visit 1
    etabvv_cl_3 ~ fixed(log(1 + 0.184^2)) # same variance as visit 1
    etabvv_tlag_1 ~ log(1 + 0.099^2) # Table 2 Delay in absorption for night dose BVV 9.90% CV (9.52-15.4)
    etabvv_tlag_2 ~ fixed(log(1 + 0.099^2)) # same variance as visit 1
    etabvv_tlag_3 ~ fixed(log(1 + 0.099^2)) # same variance as visit 1

    # Between-occasion variability (each dose an occasion; six occasions).
    etaiov_ka_1 ~ log(1 + 0.555^2) # Table 2 ka BOV 55.5% CV (47.0-72.0)
    etaiov_ka_2 ~ fixed(log(1 + 0.555^2)) # same variance as occasion 1
    etaiov_ka_3 ~ fixed(log(1 + 0.555^2)) # same variance as occasion 1
    etaiov_ka_4 ~ fixed(log(1 + 0.555^2)) # same variance as occasion 1
    etaiov_ka_5 ~ fixed(log(1 + 0.555^2)) # same variance as occasion 1
    etaiov_ka_6 ~ fixed(log(1 + 0.555^2)) # same variance as occasion 1
    etaiov_fdepot_1 ~ log(1 + 0.448^2) # Table 2 Bioavailability BOV 44.8% CV (32.9-45.9)
    etaiov_fdepot_2 ~ fixed(log(1 + 0.448^2)) # same variance as occasion 1
    etaiov_fdepot_3 ~ fixed(log(1 + 0.448^2)) # same variance as occasion 1
    etaiov_fdepot_4 ~ fixed(log(1 + 0.448^2)) # same variance as occasion 1
    etaiov_fdepot_5 ~ fixed(log(1 + 0.448^2)) # same variance as occasion 1
    etaiov_fdepot_6 ~ fixed(log(1 + 0.448^2)) # same variance as occasion 1
    etaiov_mtt_1 ~ log(1 + 1.75^2) # Table 2 Absorption mean transit time BOV 175% CV (141-205)
    etaiov_mtt_2 ~ fixed(log(1 + 1.75^2)) # same variance as occasion 1
    etaiov_mtt_3 ~ fixed(log(1 + 1.75^2)) # same variance as occasion 1
    etaiov_mtt_4 ~ fixed(log(1 + 1.75^2)) # same variance as occasion 1
    etaiov_mtt_5 ~ fixed(log(1 + 1.75^2)) # same variance as occasion 1
    etaiov_mtt_6 ~ fixed(log(1 + 1.75^2)) # same variance as occasion 1

    # Residual error, combined proportional and additive.
    propSd <- 0.230
    label("Proportional residual error (fraction)") # Table 2 Proportional error 23.0 (21.2-25.0) percent
    addSd <- 0.000310
    label("Estimated additive residual error before the LLOQ floor (mg/L)") # Table 2 Additive error 0.310 (0.0661-0.400) ug/liter; /1000 to mg/L
  })

  model({
    # 1. Occasion indicators: OCC 2k - 1 = evening dose, OCC 2k = observed
    # morning dose of PK visit k.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)

    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2 + oc3 * etaiov_ka_3 +
      oc4 * etaiov_ka_4 + oc5 * etaiov_ka_5 + oc6 * etaiov_ka_6
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 +
      oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4 +
      oc5 * etaiov_fdepot_5 + oc6 * etaiov_fdepot_6
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2 + oc3 * etaiov_mtt_3 +
      oc4 * etaiov_mtt_4 + oc5 * etaiov_mtt_5 + oc6 * etaiov_mtt_6

    # Evening and morning doses of the same visit share the visit-level etas.
    bvv_cl <- (oc1 + oc2) * etabvv_cl_1 + (oc3 + oc4) * etabvv_cl_2 +
      (oc5 + oc6) * etabvv_cl_3
    bvv_tlag <- (oc1 + oc2) * etabvv_tlag_1 + (oc3 + oc4) * etabvv_tlag_2 +
      (oc5 + oc6) * etabvv_tlag_3

    # 2. Maturation of CL on postmenstrual age in months (Rabie 2020 Eq. 1),
    # as printed: not normalised to a reference age, so lcl is the fully
    # mature value. Half maturation at 10.7 months PMA (about 1.7 months
    # postnatal); about 0.99 at 2 years postnatal.
    tm50_cl <- exp(ltm50_cl)
    hill_cl <- exp(lhill_cl)
    fmat <- PAGE^hill_cl / (tm50_cl^hill_cl + PAGE^hill_cl)

    # 3. Individual parameters; allometric reference 9.4 kg.
    cl <- exp(lcl + etalcl + bvv_cl) * (WT / 9.4)^e_wt_cl_q * fmat
    vc <- exp(lvc) * (WT / 9.4)^e_wt_vc_vp
    vp <- exp(lvp + etalvp) * (WT / 9.4)^e_wt_vc_vp
    q <- exp(lq + etalq) * (WT / 9.4)^e_wt_cl_q
    ka <- exp(lka + iov_ka)
    mtt <- exp(lmtt + iov_mtt)
    nn <- exp(lnn)
    fdepot <- exp(lfdepot + iov_fdepot) *
      (1 + e_rif_lpvr4_fdepot * CONMED_RIF_LPVR4)
    tlag <- exp(ltlag + bvv_tlag)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODEs. Savic transit absorption: transit() supplies the analytical
    # gamma-shaped input of the most recent dose (podo(), tad()) scaled by
    # fdepot; f(depot) <- 0 stops the dose from also arriving as a bolus.
    d/dt(depot) <- transit(nn, mtt, fdepot) - ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central +
      k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- 0
    # The evening dose before a PK visit is absorbed after an extra delay.
    alag(depot) <- tlag * DOSETIME_EVENING

    # 5. Observation. Dose mg / volume L = mg/L. Rabie 2020 Methods: 'The
    # additive error for all samples was bound to be at least 20% of the
    # lower limit of quantification (LLOQ)'; LLOQ = 0.0243 mg/L, so the
    # effective additive SD is max(addSd, 0.2 * 0.0243 mg/L).
    Cc <- central / vc
    addSd_eff <- max(addSd, 0.2 * 0.0243)
    Cc ~ add(addSd_eff) + prop(propSd)
  })
}
