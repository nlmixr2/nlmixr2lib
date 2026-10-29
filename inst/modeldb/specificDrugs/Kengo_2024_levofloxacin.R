Kengo_2024_levofloxacin <- function() {
  description <- paste(
    "One-compartment population PK model for oral levofloxacin in South",
    "African adults with drug-resistant tuberculosis, with a clofazimine",
    "drug-drug interaction on clearance (Kengo 2024; n = 21, PK sampled",
    "before and 3 weeks after clofazimine was added to the regimen).",
    "Absorption is a Savic transit chain (MTT 1.32 h, 15.7 transit",
    "compartments) feeding a first-order depot (ka 2.05 1/h); elimination",
    "is first order (CL 6.81 L/h, V 98.5 L). Clearance and volume are",
    "scaled allometrically by fat-free mass (reference 42 kg, exponents",
    "0.75 and 1). Clearance rises 4.70% per 10 mL/min of Cockcroft-Gault",
    "creatinine clearance above the 81.8 mL/min median and is 14.8% lower",
    "when clofazimine is co-administered. Random effects are between-subject",
    "and between-visit variability on CL and between-occasion variability",
    "on ka, MTT and F, the latter inflated 1.99-fold for occasions whose",
    "dose was not observed.",
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
    depot = list(analyte = "levofloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "levofloxacin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    FFM = list(
      description = "Fat-free mass",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric size descriptor on CL (exponent 0.75) and V (exponent",
        "1) (Kengo 2024 Methods; Results, dOFV -17.6). The Table 2 typical",
        "values refer to a participant with FFM 42 kg (Table 2 footnote a),",
        "which is therefore the reference. Linezolid/levofloxacin cohort",
        "median 41.5 kg (range 26.0-51.6, Table 1). The paper does not print",
        "its FFM formula; the same group uses the Janmahasatian equation in",
        "its sibling models.",
        sep = " "
      ),
      source_name = "FFM"
    ),
    CRCL = list(
      description = "Creatinine clearance by the Cockcroft-Gault formula, NOT normalised to body surface area",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Raw Cockcroft-Gault creatinine clearance in mL/min (Table 1",
        "footnote a; Methods), supplied in the canonical CRCL column without",
        "BSA normalisation, as for the other raw Cockcroft-Gault models in",
        "the register. Linear effect centred on the cohort median: CL rises",
        "4.70% (95% CI 1.58-7.77) for every 10 mL/min above the median",
        "(Results; Table 2 'Effect of CLCR (+10 mL/min change) on CL'),",
        "dOFV -4.53, 1 df, P = 0.033. The median is taken as the 81.8 mL/min",
        "of the linezolid/levofloxacin column of Table 1 (range 50.5-136).",
        sep = " "
      ),
      source_name = "CLCR"
    ),
    CONMED_CLOFAZIMINE = list(
      description = "Concomitant clofazimine indicator: 1 = clofazimine co-administered, 0 = not",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no clofazimine; study visit 1)",
      notes = paste(
        "Clofazimine 100 mg once daily was withheld for the first 2 weeks",
        "and started after PK visit 1, so CONMED_CLOFAZIMINE = 0 at visit 1",
        "(occasions 1-2) and 1 at visit 2 (occasions 3-4, 21 days after",
        "clofazimine was started). Fractional effect on CL, 1 +",
        "e_conmed_clofazimine_cl * CONMED_CLOFAZIMINE, with",
        "e_conmed_clofazimine_cl = -0.148 (95% CI -0.232 to -0.0588; dOFV",
        "-10.7, 1 df, P = 0.001), raising the median AUC from 140 to 166",
        "mg*h/L. The study was not a cross-over, so the paper cannot exclude",
        "that some of the visit-2 change reflects other factors such as",
        "clinical improvement (Discussion).",
        sep = " "
      ),
      source_name = "CFZ"
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
        "occasions carry the 1.99-fold inflation of the absorption BOV",
        "(Table 2 footnote h). Visit 1 is occasions 1 and 2, visit 2 is",
        "occasions 3 and 4, from which the model derives the visit for the",
        "between-visit variability on CL. For a prospective simulation in",
        "which every dose is observed, use OCC = 2.",
        sep = " "
      ),
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Median 52 kg (range 37-74, Table 1). Tested as an allometric size descriptor and replaced by fat-free mass; still needed as an input to the FFM formula and to the Cockcroft-Gault creatinine clearance."
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
    dose_range = "Levofloxacin 750 or 1000 mg orally once daily (median 19.1 mg/kg, range 13.5-26.3), taken with a meal under direct observation.",
    regions = "South Africa (Klerksdorp/Tshepong Hospital, North-West Province)",
    co_medication = "Standard drug-resistant TB regimen; clofazimine 100 mg daily was withheld for the first 2 weeks and added after the first PK visit.",
    notes = paste(
      "Non-interventional prospective cohort. PK sampling pre-dose and 2,",
      "4, 6, 8 and 10 h post-dose on day 14 (visit 1, no clofazimine) and",
      "21 days after starting clofazimine (visit 2). LC-MS/MS assay, LLOQ",
      "0.0781 mg/L; BLQ data imputed to LLOQ/2 (Beal M6). NONMEM 7.5.1;",
      "parameter uncertainty by sampling importance resampling.",
      sep = " "
    )
  )

  ini({
    # Typical values from Kengo 2024 Table 2 (levofloxacin column), for a
    # participant with FFM 42 kg (Table 2 footnote a), creatinine clearance
    # at the median and no clofazimine.
    lcl <- log(6.81)
    label("Clearance at FFM 42 kg, median creatinine clearance, no clofazimine (L/h)")    # Table 2 'Clearance (L/h)' 6.81 (6.05-7.63)
    lvc <- log(98.5)
    label("Volume of distribution at FFM 42 kg (L)")                                      # Table 2 'Volume of distribution (central compartment)' 98.5 (90.3-107)
    lka <- log(2.05)
    label("First-order absorption rate constant from the depot (1/h)")                    # Table 2 'Absorption rate constant, ka' 2.05 (1.19-3.25); printed unit '/L' is 1/h
    lmtt <- log(1.32)
    label("Mean transit time of the Savic absorption chain (h)")                          # Table 2 'Mean absorption transit time, MTT' 1.32 (1.00-1.61)
    lntr <- log(15.7)
    label("Number of transit compartments in the Savic chain (unitless)")                 # Table 2 'Transit compartments, NN' 15.7 (9.91-27.3)
    lfdepot <- fixed(log(1))
    label("Oral bioavailability (fraction)")                                              # Table 2 'Bioavailability, F' 1 (fixed), footnote b

    e_ffm_cl <- fixed(0.75)
    label("Allometric exponent of fat-free mass on clearance (unitless)")                 # Methods 'fixed exponents of 0.75 and 1'
    e_ffm_vc <- fixed(1)
    label("Allometric exponent of fat-free mass on volume (unitless)")                    # Methods 'fixed exponents of 0.75 and 1'

    # Linear creatinine-clearance effect: +4.70% per +10 mL/min, i.e.
    # 0.00470 per mL/min, centred on the 81.8 mL/min median.
    e_crcl_cl <- 0.00470
    label("Fractional change in CL per mL/min of creatinine clearance from 81.8 mL/min (1/(mL/min))")  # Table 2 'Effect of CLCR (+10 mL/min change) on CL (%)' +4.70 (+1.58 to +7.77)
    e_conmed_clofazimine_cl <- -0.148
    label("Fractional change in CL with clofazimine co-administration (unitless)")        # Table 2 'Effect of clofazimine on clearance (%)' -14.8 (-23.2 to -5.88)

    # Table 2 footnote g: %CV = sqrt(omega^2) * 100, so each variance below
    # is (%CV / 100)^2.
    etalcl ~ 0.033489 # Table 2 'BSV in clearance' 18.3% (13.2-24.0); 0.183^2

    etabvv_cl_1 ~ 0.0289        # Table 2 'BVV in clearance' 17.0% (9.69-22.0); 0.170^2
    etabvv_cl_2 ~ fixed(0.0289) # visit 2, same variance as visit 1

    # Between-occasion variability on the absorption parameters. Even
    # occasions (observed doses) carry the estimated variance. Odd
    # occasions (the unobserved dose before each pre-dose sample) carry it
    # inflated by the Table 2 'Scaling factor on BOV for unobserved dose'
    # 1.99 (1.40-2.92), footnote h. The factor multiplies the deviate, so
    # the variance is multiplied by 1.99^2 = 3.9601. Footnote h applies it
    # to the 'absorption parameters', taken here as ka, MTT and F, the three
    # parameters that carry BOV.
    etaiov_ka_1 ~ fixed(5.514043) # unobserved dose: 3.9601 * 1.3924
    etaiov_ka_2 ~ 1.3924          # Table 2 'BOV in ka' 118% (88.9-156); 1.18^2
    etaiov_ka_3 ~ fixed(5.514043) # unobserved dose, visit 2
    etaiov_ka_4 ~ fixed(1.3924)   # observed dose, visit 2; same variance as occasion 2
    etaiov_mtt_1 ~ fixed(0.780678) # unobserved dose: 3.9601 * 0.197136
    etaiov_mtt_2 ~ 0.197136        # Table 2 'BOV in MTT' 44.4% (31.6-56.0); 0.444^2
    etaiov_mtt_3 ~ fixed(0.780678) # unobserved dose, visit 2
    etaiov_mtt_4 ~ fixed(0.197136) # observed dose, visit 2; same variance as occasion 2
    etaiov_fdepot_1 ~ fixed(0.183056) # unobserved dose: 3.9601 * 0.046225
    etaiov_fdepot_2 ~ 0.046225        # Table 2 'BOV in F' 21.5% (16.5-27.9); 0.215^2
    etaiov_fdepot_3 ~ fixed(0.183056) # unobserved dose, visit 2
    etaiov_fdepot_4 ~ fixed(0.046225) # observed dose, visit 2; same variance as occasion 2

    propSd <- 0.0497
    label("Proportional residual error (fraction)")                                       # Table 2 'Proportional error (%)' 4.97 (3.81-6.01)
    addSd <- 0.238
    label("Additive residual error (mg/L)")                                               # Table 2 'Additive error (mg/L)' 0.238 (0.121-0.375)
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

    crcl_cl <- 1 + e_crcl_cl * (CRCL - 81.8)
    cfz_cl <- 1 + e_conmed_clofazimine_cl * CONMED_CLOFAZIMINE

    cl <- exp(lcl + etalcl + bvv_cl) * (FFM / 42)^e_ffm_cl * crcl_cl * cfz_cl
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
