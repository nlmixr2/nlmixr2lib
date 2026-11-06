Kengo_2024_isoniazid <- function() {
  description <- paste(
    "Two-compartment population PK model for oral high-dose isoniazid in",
    "South African adults with drug-resistant tuberculosis, with first-pass",
    "extraction through a well-stirred liver (Kengo 2024; n = 27, PK sampled",
    "before and 3 weeks after clofazimine was added to the regimen).",
    "Absorption is a Savic transit chain (MTT 0.122 h, 2.32 transit",
    "compartments) feeding a first-order depot (ka 3.28 1/h), all three",
    "fixed to literature values. Intrinsic clearance depends on the NAT2",
    "acetylator phenotype (12.2, 24.9 and 46.7 L/h for slow, intermediate",
    "and rapid acetylators); hepatic plasma flow (72.4 L/h) and the unbound",
    "fraction (0.95) are fixed. All clearance and volume parameters are",
    "scaled allometrically by fat-free mass (reference 42 kg, exponents 0.75",
    "and 1). Clofazimine co-administration had no effect on isoniazid PK.",
    "Random effects are between-subject and between-visit variability on",
    "intrinsic clearance and between-occasion variability on ka, MTT and",
    "pre-hepatic bioavailability.",
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
    depot = list(analyte = "isoniazid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "isoniazid", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "isoniazid", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    FFM = list(
      description = "Fat-free mass",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric size descriptor on every clearance (exponent 0.75: CLint,",
        "Q and hepatic plasma flow) and volume (exponent 1: Vc, Vp) parameter",
        "(Kengo 2024 Methods, 'Pharmacokinetic modeling'; Results, 'Inclusion",
        "of allometric scaling by fat-free mass (FFM) on all disposition",
        "parameters improved the fit'). The Table 2 typical values refer to a",
        "participant with FFM 42 kg (Table 2 footnote a), which is therefore",
        "the reference. Isoniazid cohort median 44.2 kg (range 26.0-56.5,",
        "Table 1). The paper does not print its FFM formula; the same group",
        "uses the Janmahasatian equation in its sibling models.",
        sep = " "
      ),
      source_name = "FFM"
    ),
    NAT2_SLOW = list(
      description = "NAT2 slow-acetylator phenotype indicator: 1 = slow, 0 = intermediate or rapid",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (intermediate or rapid); intermediate is the joint NAT2_SLOW = 0, NAT2_RAPID = 0 state",
      notes = paste(
        "Phenotype from Sanger sequencing of four NAT2 SNPs (rs1801279,",
        "rs1801280, rs1799930, rs1799931): heterozygous for two or more, or",
        "homozygous variant for any, = slow (Methods, 'NAT2 phenotype').",
        "5 of 27 participants (19%) were slow acetylators (Table 1). Selects",
        "the slow-acetylator intrinsic clearance (dOFV -39.0 for the",
        "three-level phenotype, 2 df).",
        sep = " "
      ),
      source_name = "NAT2"
    ),
    NAT2_RAPID = list(
      description = "NAT2 rapid-acetylator phenotype indicator: 1 = rapid, 0 = intermediate or slow",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (intermediate or slow); intermediate is the joint NAT2_SLOW = 0, NAT2_RAPID = 0 state",
      notes = paste(
        "Rapid = homozygous for the common allele of all four typed SNPs;",
        "intermediate = heterozygous for exactly one (Methods). 9 of 27",
        "participants (33%) were rapid and 13 (48%) intermediate (Table 1).",
        sep = " "
      ),
      source_name = "NAT2"
    ),
    OCC = list(
      description = "Dosing-occasion index: 1 = visit-1 unobserved dose, 2 = visit-1 observed dose, 3 = visit-2 unobserved dose, 4 = visit-2 observed dose",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Each dose was a separate occasion for the between-occasion",
        "variability (Methods), and each PK visit had a pre-dose sample that",
        "followed the previous day's dose plus a 2-10 h profile after the",
        "observed dose. Visit 1 (day 14, before clofazimine) is occasions 1",
        "and 2; visit 2 (21 days after clofazimine was started) is",
        "occasions 3 and 4. The model derives the visit from OCC for the",
        "between-visit variability on intrinsic clearance. For a prospective",
        "simulation, OCC = 2 throughout is a single-occasion, single-visit",
        "simulation.",
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
        "The question the study was designed to answer. Tested as a",
        "categorical covariate on clearance and bioavailability and not",
        "retained: the 10% lower clearance at visit 2 was not significant",
        "(dOFV -3.52, 1 df, P = 0.061). Simulation-estimation gave 97% power",
        "to detect a 30% change in clearance.",
        sep = " "
      )
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Median 53 kg (range 37-66, Table 1). Tested as an allometric size descriptor and replaced by fat-free mass; still needed as an input to the FFM formula."
    ),
    CRCL = list(
      description = "Creatinine clearance (Cockcroft-Gault, not BSA-normalised)",
      units = "mL/min",
      type = "continuous",
      notes = "Median 89.3 mL/min (range 35.0-143, Table 1). Explored in the covariate analysis (Methods) and not retained for isoniazid."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 27L,
    n_studies = 1L,
    age_range = "20-62 years (median 39)",
    weight_range = "37-66 kg (median 53)",
    ffm_range = "26.0-56.5 kg (median 44.2)",
    height_range = "1.30-1.81 m (median 1.67)",
    crcl_range = "35.0-143 mL/min Cockcroft-Gault (median 89.3)",
    sex_female_pct = 29.6,
    race_ethnicity = "Not reported; single South African site",
    n_hiv_positive = 23L,
    nat2_phenotype = "5 slow (19%), 13 intermediate (48%), 9 rapid (33%)",
    disease_state = "Newly diagnosed pulmonary drug-resistant (rifampicin-resistant) tuberculosis; 23 of 27 (85%) living with HIV.",
    dose_range = "High-dose isoniazid 450 or 600 mg orally once daily (median 10.0 mg/kg, range 7.35-12.1), taken with a meal under direct observation.",
    regions = "South Africa (Klerksdorp/Tshepong Hospital, North-West Province)",
    co_medication = "Standard drug-resistant TB regimen; clofazimine 100 mg daily was withheld for the first 2 weeks and added after the first PK visit.",
    notes = paste(
      "Non-interventional prospective cohort. PK sampling pre-dose and 2,",
      "4, 6, 8 and 10 h post-dose on day 14 (visit 1, no clofazimine) and",
      "21 days after starting clofazimine (visit 2). LC-MS/MS assay, LLOQ",
      "0.105 mg/L; BLQ data imputed to LLOQ/2 (Beal M6). NONMEM 7.5.1;",
      "parameter uncertainty by sampling importance resampling.",
      sep = " "
    )
  )

  ini({
    # Typical values from Kengo 2024 Table 2 (isoniazid column), for a
    # participant with FFM 42 kg (Table 2 footnote a).
    lclint_slow <- log(12.2)
    label("Intrinsic clearance in NAT2 slow acetylators at FFM 42 kg (L/h)")               # Table 2 'Intrinsic clearance (slow acetylators)' 12.2 (9.53-14.5)
    lclint_int <- log(24.9)
    label("Intrinsic clearance in NAT2 intermediate acetylators at FFM 42 kg (L/h)")       # Table 2 'Intrinsic clearance (intermediate acetylators)' 24.9 (19.0-27.0)
    lclint_rapid <- log(46.7)
    label("Intrinsic clearance in NAT2 rapid acetylators at FFM 42 kg (L/h)")              # Table 2 'Intrinsic clearance (rapid acetylators)' 46.7 (37.9-56.4)
    lvc <- log(57.0)
    label("Central volume of distribution at FFM 42 kg (L)")                               # Table 2 'Volume of distribution (central compartment)' 57.0 (49.0-61.9)
    lq <- log(1.70)
    label("Intercompartmental clearance at FFM 42 kg (L/h)")                               # Table 2 'Inter-compartmental clearance' 1.70 (0.806-2.27)
    lvp <- log(42.4)
    label("Peripheral volume of distribution at FFM 42 kg (L)")                            # Table 2 'Volume of distribution (peripheral compartment)' 42.4 (9.88-52.9)

    # Absorption fixed to literature values (Table 2 footnote e, citing
    # Gausi 2021). Table 2 prints the ka unit as '/L'; it is 1/h.
    lka <- fixed(log(3.28))
    label("First-order absorption rate constant from the depot (1/h)")                     # Table 2 'Absorption rate constant, ka' 3.28 (fixed), footnote e
    lmtt <- fixed(log(0.122))
    label("Mean transit time of the Savic absorption chain (h)")                           # Table 2 'Mean absorption transit time, MTT' 0.122 (fixed), footnote e
    lntr <- fixed(log(2.32))
    label("Number of transit compartments in the Savic chain (unitless)")                  # Table 2 'Transit compartments, NN' 2.32 (fixed), footnote e
    lfdepot <- fixed(log(1))
    label("Pre-hepatic bioavailability (fraction)")                                        # Table 2 'Bioavailability, F' 1 (fixed), footnotes b and d

    # Well-stirred liver constants, fixed (Table 2 footnote f).
    lqh <- fixed(log(72.4))
    label("Hepatic plasma flow at FFM 42 kg (L/h)")                                        # Table 2 'Hepatic blood flow rate, QH' 72.4 (fixed), footnote f
    fu <- fixed(0.95)
    label("Unbound fraction of isoniazid in plasma (fraction)")                            # Table 2 'Unbound fraction, fu (%)' 95 (fixed), footnote f

    e_ffm_cl <- fixed(0.75)
    label("Allometric exponent of fat-free mass on clearances and hepatic flow (unitless)") # Methods 'fixed exponents of 0.75 and 1'
    e_ffm_vc <- fixed(1)
    label("Allometric exponent of fat-free mass on volumes (unitless)")                    # Methods 'fixed exponents of 0.75 and 1'

    # Table 2 footnote g: %CV = sqrt(omega^2) * 100, so each variance below
    # is (%CV / 100)^2.
    etalclint ~ 0.026244 # Table 2 'BSV in clearance' 16.2% (8.10-26.1); 0.162^2

    # Between-visit variability on intrinsic clearance: one variance shared
    # by the two PK visits, the visit-2 slot fixed equal to visit 1.
    etabvv_clint_1 ~ 0.033856        # Table 2 'BVV in clearance' 18.4% (11.5-25.0); 0.184^2
    etabvv_clint_2 ~ fixed(0.033856) # visit 2, same variance as visit 1

    # Between-occasion variability, one variance shared by the four
    # occasions. Unlike linezolid and levofloxacin, isoniazid has no
    # inflation factor for the unobserved doses (Table 2 prints a dash).
    etaiov_ka_1 ~ 0.889249        # Table 2 'BOV in ka' 94.3% (4.27-101); 0.943^2
    etaiov_ka_2 ~ fixed(0.889249) # same variance, occasion 2
    etaiov_ka_3 ~ fixed(0.889249) # same variance, occasion 3
    etaiov_ka_4 ~ fixed(0.889249) # same variance, occasion 4
    etaiov_mtt_1 ~ 1.6641        # Table 2 'BOV in MTT' 129% (87.0-215); 1.29^2
    etaiov_mtt_2 ~ fixed(1.6641) # same variance, occasion 2
    etaiov_mtt_3 ~ fixed(1.6641) # same variance, occasion 3
    etaiov_mtt_4 ~ fixed(1.6641) # same variance, occasion 4
    etaiov_fdepot_1 ~ 0.225625        # Table 2 'BOV in F' 47.5% (34.3-54.6); 0.475^2
    etaiov_fdepot_2 ~ fixed(0.225625) # same variance, occasion 2
    etaiov_fdepot_3 ~ fixed(0.225625) # same variance, occasion 3
    etaiov_fdepot_4 ~ fixed(0.225625) # same variance, occasion 4

    propSd <- 0.244
    label("Proportional residual error (fraction)")                                        # Table 2 'Proportional error (%)' 24.4 (21.1-28.0)
    addSd <- fixed(0.021)
    label("Additive residual error, 20% of the 0.105 mg/L LLOQ (mg/L)")                    # Table 2 'Additive error (mg/L)' 0.021, i.e., 20% LLOQ (fixed), footnote i
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
    bvv_clint <- (oc1 + oc2) * etabvv_clint_1 + (oc3 + oc4) * etabvv_clint_2

    # NAT2 phenotype selects one of three typical intrinsic clearances.
    nat2_int <- (1 - NAT2_SLOW) * (1 - NAT2_RAPID)
    lclint_typ <- NAT2_SLOW * lclint_slow + nat2_int * lclint_int + NAT2_RAPID * lclint_rapid

    clint <- exp(lclint_typ + etalclint + bvv_clint) * (FFM / 42)^e_ffm_cl
    vc <- exp(lvc) * (FFM / 42)^e_ffm_vc
    q <- exp(lq) * (FFM / 42)^e_ffm_cl
    vp <- exp(lvp) * (FFM / 42)^e_ffm_vc
    qh <- exp(lqh) * (FFM / 42)^e_ffm_cl
    ka <- exp(lka + iov_ka)
    mtt <- exp(lmtt + iov_mtt)
    ntr <- exp(lntr)
    fdepot <- exp(lfdepot + iov_fdepot)

    # Well-stirred liver. The absorbed dose passes the liver once (fraction
    # fh survives) and systemic drug is cleared at qh * eh. The liver is
    # taken as instantly equilibrating with no volume of its own, the form
    # the same group's isoniazid control streams use (EH, FH = 1 - EH and
    # K20 = QH * EH / V; see Calderin_2025_isoniazid.R). Kengo 2024 reports
    # no liver volume.
    eh <- (clint * fu) / (clint * fu + qh)
    fh <- 1 - eh
    kel <- qh * eh / vc
    k12 <- q / vc
    k21 <- q / vp

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
    d/dt(central) <- ka * depot * fh - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    f(depot) <- 0

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
