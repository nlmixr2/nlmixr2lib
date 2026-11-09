Kengo_2024_linezolid <- function() {
  description <- paste(
    "One-compartment population PK model for oral linezolid in South",
    "African adults with drug-resistant tuberculosis (Kengo 2024; n = 21,",
    "PK sampled before and 3 weeks after clofazimine was added to the",
    "regimen). Absorption is a Savic transit chain (MTT 0.928 h, 12.6",
    "transit compartments) feeding a first-order depot (ka 1.14 1/h);",
    "elimination is first order (CL 3.04 L/h, V 38.7 L). Clearance and",
    "volume are scaled allometrically by fat-free mass (reference 42 kg,",
    "exponents 0.75 and 1). Clofazimine co-administration had no effect on",
    "linezolid PK. Random effects are between-subject and between-visit",
    "variability on CL and between-occasion variability on ka, MTT and F,",
    "the latter inflated 1.60-fold for occasions whose dose was not",
    "observed.",
    sep = " "
  )
  reference <- paste(
    "Kengo A, Nabeemeeah F, Denti P, Sabet R, Okyere-Manu G, Abraham P,",
    "Weisner L, Mosala MH, Tshabalala S, Scholefield J, Resendiz-Galvan JE,",
    "Martinson NA, Variava E, Leggett JE.",
    "Assessing potential drug-drug interactions between clofazimine and",
    "other frequently used agents to treat drug-resistant tuberculosis.",
    "Antimicrob Agents Chemother. 2024;68(5):e01583-23.",
    "doi:10.1128/aac.01583-23",
    sep = " "
  )
  vignette <- "Kengo_2024_clofazimine_ddi"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "linezolid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "linezolid", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    FFM = list(
      description = "Fat-free mass",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric size descriptor on CL (exponent 0.75) and V (exponent",
        "1) (Kengo 2024 Methods; Results, dOFV -4.82). The Table 2 typical",
        "values refer to a participant with FFM 42 kg (Table 2 footnote a),",
        "which is therefore the reference. Linezolid/levofloxacin cohort",
        "median 41.5 kg (range 26.0-51.6, Table 1). The paper does not print",
        "its FFM formula; the same group uses the Janmahasatian equation in",
        "its sibling models.",
        sep = " "
      ),
      source_name = "FFM"
    ),
    OCC = list(
      description = "Dosing-occasion index: 1 = visit-1 unobserved dose, 2 = visit-1 observed dose, 3 = visit-2 unobserved dose, 4 = visit-2 observed dose",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Each dose was a separate occasion for the between-occasion",
        "variability (Methods). Each PK visit had a pre-dose sample that",
        "followed the previous day's unobserved dose (odd occasions) and a",
        "2-10 h profile after the observed dose (even occasions). The odd",
        "occasions carry the 1.60-fold inflation of the absorption BOV",
        "(Table 2 footnote h). Visit 1 (day 14, before clofazimine) is",
        "occasions 1 and 2; visit 2 (21 days after clofazimine was started)",
        "is occasions 3 and 4, from which the model derives the visit for",
        "the between-visit variability on CL. For a prospective simulation",
        "in which every dose is observed, use OCC = 2.",
        sep = " "
      ),
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    CONMED_CLOFAZIMINE = list(
      description = "Concomitant clofazimine indicator: 1 = clofazimine co-administered (study visit 2), 0 = not (visit 1)",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Tested as a categorical covariate on clearance and bioavailability",
        "and not retained: a non-significant 2% decrease in clearance at",
        "visit 2 (dOFV -3.7, 1 df, P = 0.054). Simulation-estimation gave",
        "90% power to detect a 30% change in clearance.",
        sep = " "
      )
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Median 52 kg (range 37-74, Table 1). Tested as an allometric size descriptor and replaced by fat-free mass; still needed as an input to the FFM formula."
    ),
    CRCL = list(
      description = "Creatinine clearance (Cockcroft-Gault, not BSA-normalised)",
      units = "mL/min",
      type = "continuous",
      notes = "Median 81.8 mL/min (range 50.5-136, Table 1). Explored in the covariate analysis (Methods) and not retained for linezolid."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 21L,
    n_studies = 1L,
    age_range = "27-68 years (median 39)",
    weight_range = "37-74 kg (median 52)",
    ffm_range = "26.0-51.6 kg (median 41.5)",
    height_range = "1.49-1.80 m (median 1.64)",
    crcl_range = "50.5-136 mL/min Cockcroft-Gault (median 81.8)",
    sex_female_pct = 38.1,
    race_ethnicity = "Not reported; single South African site",
    n_hiv_positive = 16L,
    disease_state = "Newly diagnosed pulmonary drug-resistant (rifampicin-resistant) tuberculosis; 16 of 21 (76%) living with HIV.",
    dose_range = "Linezolid 600 mg orally once daily (median 11.5 mg/kg, range 8.11-16.2), taken with a meal under direct observation.",
    regions = "South Africa (Klerksdorp/Tshepong Hospital, North-West Province)",
    co_medication = "Standard drug-resistant TB regimen; clofazimine 100 mg daily was withheld for the first 2 weeks and added after the first PK visit.",
    notes = paste(
      "Non-interventional prospective cohort. PK sampling pre-dose and 2,",
      "4, 6, 8 and 10 h post-dose on day 14 (visit 1, no clofazimine) and",
      "21 days after starting clofazimine (visit 2). LC-MS/MS assay, LLOQ",
      "0.100 mg/L; BLQ data imputed to LLOQ/2 (Beal M6). NONMEM 7.5.1;",
      "parameter uncertainty by sampling importance resampling.",
      sep = " "
    )
  )

  ini({
    # Typical values from Kengo 2024 Table 2 (linezolid column), for a
    # participant with FFM 42 kg (Table 2 footnote a).
    lcl <- log(3.04)
    label("Clearance at FFM 42 kg (L/h)")                                                 # Table 2 'Clearance (L/h)' 3.04 (2.74-3.87)
    lvc <- log(38.7)
    label("Volume of distribution at FFM 42 kg (L)")                                      # Table 2 'Volume of distribution (central compartment)' 38.7 (37.0-43.1)
    lka <- log(1.14)
    label("First-order absorption rate constant from the depot (1/h)")                    # Table 2 'Absorption rate constant, ka' 1.14 (0.768-1.77); printed unit '/L' is 1/h
    lmtt <- log(0.928)
    label("Mean transit time of the Savic absorption chain (h)")                          # Table 2 'Mean absorption transit time, MTT' 0.928 (0.713-1.36)
    lntr <- log(12.6)
    label("Number of transit compartments in the Savic chain (unitless)")                 # Table 2 'Transit compartments, NN' 12.6 (9.98-37.0)
    lfdepot <- fixed(log(1))
    label("Oral bioavailability (fraction)")                                              # Table 2 'Bioavailability, F' 1 (fixed), footnote b

    e_ffm_cl <- fixed(0.75)
    label("Allometric exponent of fat-free mass on clearance (unitless)")                 # Methods 'fixed exponents of 0.75 and 1'
    e_ffm_vc <- fixed(1)
    label("Allometric exponent of fat-free mass on volume (unitless)")                    # Methods 'fixed exponents of 0.75 and 1'

    # Table 2 footnote g: %CV = sqrt(omega^2) * 100, so each variance below
    # is (%CV / 100)^2.
    etalcl ~ 0.099856 # Table 2 'BSV in clearance' 31.6% (13.2-41.5); 0.316^2

    etabvv_cl_1 ~ 0.110889        # Table 2 'BVV in clearance' 33.3% (25.4-49.2); 0.333^2
    etabvv_cl_2 ~ fixed(0.110889) # visit 2, same variance as visit 1

    # Between-occasion variability on the absorption parameters. Even
    # occasions (observed doses) carry the estimated variance. Odd
    # occasions (the unobserved dose before each pre-dose sample) carry it
    # inflated by the Table 2 'Scaling factor on BOV for unobserved dose'
    # 1.60 (1.07-2.67), footnote h. The factor multiplies the deviate, so
    # the variance is multiplied by 1.60^2 = 2.56. Footnote h applies it to
    # the 'absorption parameters', taken here as ka, MTT and F, the three
    # parameters that carry BOV.
    etaiov_ka_1 ~ fixed(2.768896) # unobserved dose: 2.56 * 1.0816
    etaiov_ka_2 ~ 1.0816          # Table 2 'BOV in ka' 104% (74.5-127); 1.04^2
    etaiov_ka_3 ~ fixed(2.768896) # unobserved dose, visit 2
    etaiov_ka_4 ~ fixed(1.0816)   # observed dose, visit 2; same variance as occasion 2
    etaiov_mtt_1 ~ fixed(1.39429) # unobserved dose: 2.56 * 0.544644
    etaiov_mtt_2 ~ 0.544644       # Table 2 'BOV in MTT' 73.8% (47.4-101); 0.738^2
    etaiov_mtt_3 ~ fixed(1.39429) # unobserved dose, visit 2
    etaiov_mtt_4 ~ fixed(0.544644) # observed dose, visit 2; same variance as occasion 2
    etaiov_fdepot_1 ~ fixed(0.039363) # unobserved dose: 2.56 * 0.015376
    etaiov_fdepot_2 ~ 0.015376        # Table 2 'BOV in F' 12.4% (6.48-18.4); 0.124^2
    etaiov_fdepot_3 ~ fixed(0.039363) # unobserved dose, visit 2
    etaiov_fdepot_4 ~ fixed(0.015376) # observed dose, visit 2; same variance as occasion 2

    propSd <- 0.0731
    label("Proportional residual error (fraction)")                                       # Table 2 'Proportional error (%)' 7.31 (5.36-8.59)
    addSd <- 0.284
    label("Additive residual error (mg/L)")                                               # Table 2 'Additive error (mg/L)' 0.284 (0.142-0.596)
  })

  model({
    # Occasion and visit indicators; visit 1 = occasions 1-2, visit 2 =
    # occasions 3-4.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2 + oc3 * etaiov_ka_3 + oc4 * etaiov_ka_4
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2 + oc3 * etaiov_mtt_3 + oc4 * etaiov_mtt_4
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 + oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4
    bvv_cl <- (oc1 + oc2) * etabvv_cl_1 + (oc3 + oc4) * etabvv_cl_2

    cl <- exp(lcl + etalcl + bvv_cl) * (FFM / 42)^e_ffm_cl
    vc <- exp(lvc) * (FFM / 42)^e_ffm_vc
    ka <- exp(lka + iov_ka)
    mtt <- exp(lmtt + iov_mtt)
    ntr <- exp(lntr)
    fdepot <- exp(lfdepot + iov_fdepot)
    kel <- cl / vc

    # Savic transit chain in closed form, ktr = (ntr + 1) / mtt, driven by
    # the most recent dose amount and the time since it. The bolus into
    # depot is suppressed with f(depot) <- 0 so the chain is the only input.
    tdos <- tad(depot)
    ktr <- (ntr + 1) / mtt
    ktt <- ktr * tdos
    trin <- 0
    if (ktt > 0) {
      trin <- exp(log(fdepot * podo(depot) * ktr) - lgamma(ntr + 1) + ntr * log(ktt) - ktt)
    }

    d/dt(depot) <- trin - ka * depot
    d/dt(central) <- ka * depot - kel * central
    f(depot) <- 0

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
