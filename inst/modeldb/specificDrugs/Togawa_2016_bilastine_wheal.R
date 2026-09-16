Togawa_2016_bilastine_wheal <- function() {
  description <- "Sequential pharmacokinetic/pharmacodynamic model for inhibition of the histamine-induced skin WHEAL response by oral bilastine in healthy adult Japanese male volunteers. The two-compartment first-order-absorption disposition model of Togawa 2016 Table 5 is carried over unchanged and held fixed, and plasma bilastine drives a type-I indirect-response (turnover) model in which the zero-order production rate kin of the wheal area is inhibited by an Emax function of plasma concentration with half-maximal inhibitory concentration IC50, while the response is lost first-order at rate kout. The drug-free wheal area is the turnover steady state kin / kout = 0.49 cm^2. Fit by naive pooling of the single-dose (Part I) data from 27 bilastine-treated subjects: the sparse four-point-per-subject wheal sampling schedule did not support any pharmacodynamic random effect, and no pharmacodynamic residual error is reported. Wheal areas were measured by planimetry after a 10 mg/mL histamine prick test."
  reference <- paste(
    "Togawa M, Yamaya H, Rodriguez M, Nagashima H (2016).",
    "Pharmacokinetics, pharmacodynamics and population pharmacokinetic/pharmacodynamic",
    "modelling of bilastine, a second-generation antihistamine, in healthy Japanese subjects.",
    "Clin Drug Investig 36(12):1011-1021.",
    "doi:10.1007/s40261-016-0447-2.",
    sep = " "
  )
  vignette <- "Togawa_2016_bilastine"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL"
    # units$concentration documents the plasma bilastine output Cc. The second
    # output, `wheal`, is an AREA in cm^2 and is not a concentration; the
    # published kin is reported in cm^2/h (Togawa 2016 Table 6), so the state is
    # carried in cm^2. The paper's observed wheal areas are tabulated in mm^2
    # (Table 1), which is 100 x the state value.
    #
    # As in the parent PK model, the central-compartment amount unit that makes
    # `Cc <- central / vc` come out in ng/mL is ug, so event tables dose in ug.
  )

  paper_specific_compartments <- c("wheal")

  compartmentData <- list(
    depot = list(analyte = "bilastine", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "bilastine", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "bilastine", units = "ug", specimen = "plasma", verified = TRUE),
    wheal = list(analyte = "histamine-induced wheal area", units = "cm^2", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age in years at screening (per-arm means 22.9-29.8 years; Togawa 2016 Table 1).",
      units = "years",
      type = "continuous",
      notes = "Togawa 2016 Methods 2.5 applies the same stepwise covariate screen to the pharmacodynamic model as to the pharmacokinetic model; no covariate was retained in either."
    ),
    WT = list(
      description = "Body weight in kg at screening (per-arm means 61.5-64.2 kg; Togawa 2016 Table 1).",
      units = "kg",
      type = "continuous",
      notes = "Screened per Togawa 2016 Methods 2.5; not retained. The pharmacodynamic dataset (27 subjects, four post-dose wheal measurements each) was too sparse to support even a random effect, let alone a covariate effect."
    ),
    HT = list(
      description = "Standing height in cm at screening (per-arm means 170.2-173.9 cm; Togawa 2016 Table 1).",
      units = "cm",
      type = "continuous",
      notes = "Screened per Togawa 2016 Methods 2.5; not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 27L,
    n_studies = 1L,
    age_range = "20-39 years by inclusion criterion; per-arm means 22.9-29.8 years",
    weight_range = ">= 50 kg by inclusion criterion; per-arm means 61.9-62.3 kg in the three Part I bilastine arms",
    sex_female_pct = 0,
    race_ethnicity = c(Japanese = 100),
    disease_state = "Healthy volunteers with normal histamine skin-prick reactivity at screening (subjects with a > 2 mm wheal to saline or a < 3 mm wheal to 1 mg/mL histamine were excluded).",
    dose_range = "Single oral doses of bilastine 10, 20 or 50 mg under fasting conditions (n = 9 per dose).",
    regions = "Japan (single centre, Tokyo)",
    bl_wheal_mm2 = "Per-arm predose means 29.67-34.72 mm^2 (SD 7.78-12.85) across the three Part I bilastine arms; 33.72 mm^2 (SD 11.68) in the Part I placebo arm (Togawa 2016 Table 1)",
    notes = paste(
      "Part I of the study only (Togawa 2016 Table 6 footnote a, n = 27). Wheal",
      "responses were read 15 min after duplicate histamine prick tests",
      "(10 mg/mL histamine) placed at matching sites either side of the spine,",
      "at predose (14 h before dosing) and 1.5, 8, 12 and 24 h postdose; the two",
      "sites were averaged. Japanese subjects were tested with 10 mg/mL",
      "histamine whereas the Caucasian comparison study used 100 mg/mL, which",
      "the paper notes when comparing the two populations' parameters."
    )
  )

  ini({
    # ---------------------------------------------------------------
    # Pharmacokinetic layer -- FIXED, carried unchanged from the
    # population PK model of Togawa 2016 Table 5 (Japanese column).
    # "The structural PK model previously developed was at this stage
    # linked to the PD model" (Results 3.4), i.e. a sequential fit in
    # which the PK parameters are not re-estimated. See the companion
    # model file Togawa_2016_bilastine.R for the PK-only form and for
    # the source trace of each of these values.
    # ---------------------------------------------------------------
    lka <- fixed(log(1.7)); label("First-order absorption rate constant ka (1/h); from the PK model") # Togawa 2016 Table 5: ka = 1.7 1/h
    lcl <- fixed(log(14.4)); label("Apparent oral clearance CL/F (L/h); from the PK model") # Togawa 2016 Table 5: CL = 14.4 L/h
    lvc <- fixed(log(51.2)); label("Apparent central volume Vc/F (L); from the PK model") # Togawa 2016 Table 5: Vc = 51.2 L
    lq <- fixed(log(1.55)); label("Apparent intercompartmental clearance Q/F (L/h); from the PK model") # Togawa 2016 Table 5: Q = 1.55 L/h
    lvp <- fixed(log(20.2)); label("Apparent peripheral volume Vp/F (L); from the PK model") # Togawa 2016 Table 5: Vp = 20.2 L

    # ---------------------------------------------------------------
    # Pharmacodynamic layer (Togawa 2016 Table 6, 'Wheal' block,
    # 'Japanese' column). Kon is the zero-order production rate of the
    # wheal response and Koff its first-order loss rate, so the
    # drug-free steady state is Kon / Koff -- the paper states this
    # explicitly, calling Kon / Koff the "starting baseline of extent"
    # (Results 3.4). IC50 is the plasma bilastine concentration
    # producing 50% inhibition.
    # ---------------------------------------------------------------
    lkin_wheal <- log(0.84); label("Zero-order production rate of the wheal response kin (cm^2/h)") # Togawa 2016 Table 6, Wheal: Kon = 0.84 cm^2/h (%SEE 4)
    lkout_wheal <- log(1.72); label("First-order loss rate of the wheal response kout (1/h)") # Togawa 2016 Table 6, Wheal: Koff = 1.72 1/h (%SEE 7)
    lic50_wheal <- log(1.03); label("Plasma bilastine concentration producing 50% inhibition of wheal production IC50 (ng/mL)") # Togawa 2016 Table 6, Wheal: IC50 = 1.03 ng/mL (%SEE 11)

    # ---------------------------------------------------------------
    # Inter-individual variability.
    #
    # PK etas are carried fixed from Table 5 at their published diagonal
    # variances (omega reported as a percent; variance = (percent/100)^2).
    # The published PK model used a full 4 x 4 block across CL, Vc, Q and
    # Vp, but the off-diagonal covariances are not reported, so these are
    # encoded as independent -- see vignette Errata.
    #
    # PD etas are ZERO: "No random effects could be estimated for wheal,
    # but a naive pooled population model was possible" (Discussion).
    # Table 6 reports no omega column for any pharmacodynamic parameter.
    # ---------------------------------------------------------------
    etalcl ~ fixed(0.0784) # Togawa 2016 Table 5: omega CL = 28% -> variance 0.28^2
    etalvc ~ fixed(0.1156) # Togawa 2016 Table 5: omega Vc = 34% -> variance 0.34^2
    etalq ~ fixed(0.25) # Togawa 2016 Table 5: omega Q = 50% -> variance 0.50^2
    etalvp ~ fixed(0.3721) # Togawa 2016 Table 5: omega Vp = 61% -> variance 0.61^2
    etalka ~ fixed(0.0784) # Togawa 2016 Table 5: omega ka = 28% -> variance 0.28^2

    # ---------------------------------------------------------------
    # Residual error. The plasma channel keeps the published
    # proportional error of the PK model, fixed. The wheal channel has
    # NO published residual error: Table 6 reports point estimates and
    # %SEE only, and the fit was a naive pooling, so the additive SD is
    # encoded as fixed(0) rather than invented (vignette Errata).
    # ---------------------------------------------------------------
    propSd <- fixed(0.21); label("Proportional residual error on plasma bilastine (fraction); from the PK model") # Togawa 2016 Table 5: sigma = 21%
    addSd_wheal <- fixed(0); label("Additive residual SD on the wheal response (cm^2); ZERO - not reported in the source") # Togawa 2016 Table 6 reports no residual error for the naive-pooled wheal fit
  })

  model({
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    kin_wheal <- exp(lkin_wheal)
    kout_wheal <- exp(lkout_wheal)
    ic50_wheal <- exp(lic50_wheal)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc

    # Type-I indirect response: plasma bilastine inhibits the zero-order
    # production of the wheal response. The inhibition factor is 1 in the
    # absence of drug and tends to 0 at high concentration, which is what
    # produces the near-complete wheal suppression the paper reports from
    # 1.5 to 12 h after 20 and 50 mg doses.
    inh_wheal <- 1 - Cc / (ic50_wheal + Cc)

    d/dt(wheal) <- kin_wheal * inh_wheal - kout_wheal * wheal
    # Drug-free turnover steady state, stated in the source as the
    # "starting baseline of extent" Kon / Koff = 0.49 cm^2.
    wheal(0) <- kin_wheal / kout_wheal

    Cc ~ prop(propSd)
    wheal ~ add(addSd_wheal)
  })
}
