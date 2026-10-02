Wang_2020_vitacoxib_cat <- function() {
  description <- "Preclinical/clinical veterinary (cat). Two-compartment population PK model for the COX-2 selective inhibitor vitacoxib in healthy neutered domestic shorthair cats, pooling intravenous and single- and multiple-dose oral data from six studies. Oral absorption is a parallel mixed-order input: a fraction Fr of the bioavailable dose is absorbed first-order through a depot (ka = 0.13 1/h) and the remainder enters the central compartment by a zero-order process of duration Tk0 = 3.76 h; oral bioavailability 57.8%. Clearance is per kg body weight (Table 2 prints 0.11 'L/h'; the paper's Figure 6 simulations and companion NCA identify it as L/h/kg) and body weight is a power covariate on the central volume of distribution (Wang 2020)"
  reference <- paste(
    "Wang J, Schneider BK, Xiao H, Qiu J, Gong X, Seo Y-J, Li J, Mochel JP, Cao X.",
    "Non-Linear Mixed-Effects Pharmacokinetic Modeling of the Novel COX-2 Selective",
    "Inhibitor Vitacoxib in Cats. Front Vet Sci. 2020;7:554033.",
    "doi:10.3389/fvets.2020.554033.",
    sep = " "
  )
  vignette <- "Wang_2020_vitacoxib_cat"
  # Volumes in L and CL in L/h (Table 2), so doses are in mg (mg/kg x body
  # weight) and central/vc is mg/L; the x1000 in Cc reports ng/mL to match the
  # bioanalytical method (UPLC-MS/MS, LLOQ 0.5 ng/mL) and the PD targets.
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "vitacoxib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "vitacoxib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "vitacoxib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Linear (per-kg) scaling of clearance, and a power effect on the central volume of",
        "distribution centred on 2.9 kg, the",
        "cohort mean body weight (Methods 'Animals': 2.9 +/- 0.78 kg). The per-kg reading",
        "of the Table 2 clearance is justified in the lcl comment in ini(). Equation 2 prints",
        "log(V1i) = log(V1pop) + beta_V1_WT0 * WT0i + eta_V1i and the covariate is the",
        "log-normalised weight evaluated in the covariate search (Methods 'Inclusion of",
        "covariate relationships': log(bodyweight / weighted mean bodyweight)); see the",
        "e_wt_vc comment in ini() for why the raw-weight reading is excluded.",
        sep = " "
      ),
      source_name = "WT0"
    ),
    ROUTE_IV = list(
      description = "Indicator for intravenous administration",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (oral administration)",
      notes = paste(
        "Per-dose-record indicator (1 = intravenous bolus, 0 = oral). The same eight",
        "cats (IDs 9-16) received vitacoxib both i.v. (Study 2) and orally (Studies 1, 3",
        "and 6), so the column varies within subject. rxode2 applies one bioavailability",
        "per compartment rather than per administration type, so ROUTE_IV selects the",
        "bioavailability of a dose placed in central: f(central) collapses to 1 for the",
        "i.v. bolus and to F * (1 - Fr) for the zero-order part of an oral dose. An i.v.",
        "record must be a plain bolus into cmt = 'central' (no rate), which also makes",
        "rxode2 ignore the modelled dur(central). An oral administration needs TWO records",
        "at the same time: a bolus into cmt = 'depot' and a record into cmt = 'central'",
        "carrying rate = -2 so the modelled zero-order duration d1 is applied; ROUTE_IV = 0",
        "on both oral records.",
        sep = " "
      ),
      source_name = "Route (Table 1: I.V / P.O)"
    )
  )

  # Sex and feeding status were screened by the Monolix automated covariate
  # search and not retained (Results 'Parameter Estimates'); no estimate is
  # published, so they are documented rather than encoded.
  covariatesDataExcluded <- list(
    FED = list(
      description = "Fed-state indicator (1 = dosed after a meal, 0 = fasted)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (12 h overnight fast)",
      notes = paste(
        "Study 3 dosed 2 h after feeding; all other oral studies after a 12 h overnight",
        "fast (Table 1). Screened by ANOVA in Monolix 2019R2 and not retained: 'neither",
        "age, feeding status, nor sex had a statistically significant effect on the PK of",
        "vitacoxib in cats' (Results 'Parameter Estimates').",
        sep = " "
      ),
      source_name = "Feeding status"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Screened and not retained (Results 'Parameter Estimates'). The sex split of the 16 neutered cats is not reported.",
      source_name = "sex"
    )
  )

  population <- list(
    species = "cat (domestic shorthair; Felis catus)",
    n_subjects = 16L,
    n_studies = 6L,
    age_range = "1-3 years",
    weight_range = "2.9 +/- 0.78 kg (mean +/- SD)",
    weight_mean = "2.9 kg",
    sex_female_pct = NA_real_,
    disease_state = "Healthy neutered domestic shorthair laboratory cats",
    dose_range = "1, 2 and 4 mg/kg single oral dose; 2 mg/kg single i.v. dose; 2 mg/kg oral once daily for 7 days",
    regions = "China (China Agricultural University, Beijing)",
    notes = paste(
      "Methods 'Animals' and Table 1. Two groups of eight cats: IDs 9-16 received",
      "Study 1 (2 mg/kg p.o. fasted), Study 2 (2 mg/kg i.v. fasted), Study 3 (2 mg/kg",
      "p.o. 2 h after feeding) and Study 6 (2 mg/kg p.o. daily for 7 days, fasted); IDs",
      "1-8 received Study 4 (1 mg/kg p.o.) and Study 5 (4 mg/kg p.o.), both fasted, 14",
      "days apart. Washout of 2 weeks between studies. Plasma assayed by UPLC-MS/MS,",
      "LLOQ 0.5 ng/mL; BLQ data handled by the Monolix equivalent of the NONMEM M3",
      "method. Sex split not reported.",
      sep = " "
    )
  )

  ini({
    # ---- Disposition (Table 2) ----------------------------------------------
    lcl <- log(0.11)
    label("Clearance per kg body weight (L/h/kg)")
    # Table 2, row 'Systemic Clearance', CL = 0.11 (RSE 7.9%), printed with the
    # unit 'L/h' (Abstract '110 ml/h'). The value is carried as printed but
    # read as L/h/kg and multiplied by body weight in model(), because the
    # paper's own simulations and data are only consistent with a per-kg
    # clearance: (i) the deterministic (IIV- and error-free) Figure 6
    # time-above-target curves, simulated by the authors from this final model
    # with mlxR, are reproduced to an RMSE of 1.3 h across 16 dose/route/target
    # points with CL = 0.11 x 2.9 kg = 0.319 L/h and missed by an RMSE of 33 h
    # with CL = 0.11 L/h (the best-fitting multiplier is ~2.9, the cohort mean
    # weight; Q, V1 and V2 need no rescaling, and the Figure 6 dose thresholds
    # confirm V1 = 2.88 L absolute); (ii) the companion non-compartmental
    # analysis of the same Study 2 i.v. data (Wang 2019, J Vet Pharmacol Ther
    # 42:294, ref. 16) reports CL = 95.22 +/- 23.53 ml/kg/h and Vd = 1264 +/-
    # 344 ml/kg, i.e. ~0.28 L/h and ~3.7 L for a 2.9 kg cat. The Discussion's
    # derived half-life (~21 h) and extraction ratio (E < 0.01) inherit the
    # unit slip; see the vignette Errata.

    lvc <- log(2.88)
    label("Central volume of distribution for a 2.9 kg cat (L)")
    # Table 2, row 'Central compartment volume of distribution', V1 = 2.88 L
    # (RSE 25.1%).

    lvp <- log(0.54)
    label("Peripheral volume of distribution (L)")
    # Table 2, row 'Peripheral compartment volume of distribution', V2 = 0.54 L
    # (RSE 19.4%). V1 + V2 = 3.42 L, the VSS quoted in the Abstract and Results.

    lq <- log(0.52)
    label("Intercompartmental clearance (L/h)")
    # Table 2, row 'Inter-compartmental clearance', Q = 0.52 L/h (RSE 6.7%).

    # ---- Oral absorption (Table 2) -------------------------------------------
    lka <- log(0.13)
    label("First-order absorption rate constant (1/h)")
    # Table 2, row 'First-order absorption rate constant (P.O)', Ka = 0.13 1/h
    # (RSE 15.4%).

    ld1 <- log(3.76)
    label("Duration of the zero-order absorption input (h)")
    # Table 2, row 'Zero-order absorption rate constant (P.O)', Tk0 = 3.76 h
    # (RSE 6.4%). Labelled a 'rate constant' but carries units of h: it is the
    # Monolix Tk0 duration of the zero-order input (Results 'PK Model
    # Evaluation').

    logitfdepot <- log(0.578 / (1 - 0.578))
    label("Logit of the oral bioavailability F (fraction)")
    # Table 2, row 'Bioavailability (P.O)', F = 57.8% (RSE 7.1%). Methods 'Model
    # Evaluation': absorption parameters were logit-normal. logit(0.578) =
    # 0.3146.

    logitffo <- log(0.20 / (1 - 0.20))
    label("Logit of the fraction of the bioavailable oral dose absorbed first-order (fraction)")
    # Table 2, row 'Fraction absorbed through 1st order', Fr = 0.20 (RSE 16.4%).
    # Results 'PK Model Evaluation': 'The fraction of drug absorbed by the first-
    # and zero-order rate constant was represented by Fr and (1 - Fr),
    # respectively' - a parallel split of the dose (as in the same group's
    # robenacoxib cat model, Pelligand 2016, cited as ref. 17 for this
    # structure). logit(0.20) = -1.3863.

    # ---- Covariate effect ----------------------------------------------------
    e_wt_vc <- 0.41
    label("Power exponent of body weight (normalised to 2.9 kg) on the central volume")
    # Table 2, row 'Bodyweight effect on V1', beta_V1_WT0 = 0.41 (RSE 20.1%);
    # Equation 2: log(V1i) = log(V1pop) + beta_V1_WT0 * WT0i + eta_V1i. WT0 is
    # read as the log-normalised weight log(WT / 2.9) defined in Methods
    # 'Inclusion of covariate relationships'. The raw-weight reading is
    # excluded by the paper itself: exp(0.41 * 2.9) = 3.28 would put the V1 of
    # an average cat at 9.4 L, whereas the Abstract and Discussion quote a total
    # VSS of 3.42 L = V1pop + V2 for the population.

    # ---- Random effects ------------------------------------------------------
    # Table 2 footnote: 'Random effects are expressed in terms of IOV
    # (expressed as CV%) given the nature of the experimental design'; the
    # paper reports no separate IIV magnitudes and states that most variance
    # was within-subject (Results 'Parameter Estimates'). The tabulated
    # between-occasion (study-period) variability is encoded here as a single
    # random effect per parameter, i.e. one draw per animal-occasion; simulate
    # each dosing period of the same animal as a new id to reproduce the
    # occasion structure. Log-normal parameters: omega^2 = log(1 + CV^2).
    # Logit-normal parameters (F, Fr) have no closed-form CV; the tabulated
    # percentage is taken as the logit-scale SD (omega = CV/100).
    etalka ~ fixed(0.01) # Table 2 CV(Ka) = 10.0%; Results: zero-convergent random effects 'were set to 0.1 (CV ~10%)'
    etald1 ~ 0.094605 # Table 2 CV(Tk0) = 31.5%; log(1 + 0.315^2)
    etalogitffo ~ 0.0064 # Table 2 CV(Fr) = 8.0%; logit-scale SD 0.08
    etalvc ~ 0.13092 # Table 2 CV(V1) = 37.4%; log(1 + 0.374^2)
    etalvp ~ 0.78303 # Table 2 CV(V2) = 109%; log(1 + 1.09^2)
    etalq ~ fixed(0.01) # Table 2 CV(Q) = 10.0%; Results: zero-convergent random effects 'were set to 0.1 (CV ~10%)'
    etalcl ~ 0.19194 # Table 2 CV(CL) = 46.0%; log(1 + 0.46^2)
    etalogitfdepot ~ 0.070225 # Table 2 CV(F) = 26.5%; logit-scale SD 0.265

    # ---- Residual error ------------------------------------------------------
    propSd <- 0.30
    label("Proportional residual SD (fraction)")
    # Table 2, row 'Proportional error constant', b = 0.30 (RSE 3.1%); Results
    # 'PK Model Evaluation': 'a proportional error model was used'.
  })

  model({
    # ---- Individual parameters ----------------------------------------------
    # Per-kg clearance scaled by body weight (see the lcl comment in ini()).
    cl <- exp(lcl + etalcl) * WT
    vc <- exp(lvc + e_wt_vc * log(WT / 2.9) + etalvc)
    vp <- exp(lvp + etalvp)
    q <- exp(lq + etalq)
    ka <- exp(lka + etalka)
    d1 <- exp(ld1 + etald1)

    # Oral bioavailability and its split between the two parallel inputs. ffo
    # is the paper's Fr (first-order share); 1 - ffo enters central zero-order.
    fdepot <- expit(logitfdepot + etalogitfdepot)
    ffo <- expit(logitffo + etalogitffo)

    # ---- Micro-constants -----------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- ODE system (Figure 1) -----------------------------------------------
    # Two-compartment mammillary disposition with first-order elimination from
    # central. `depot` carries only the first-order share of an oral dose; the
    # zero-order share is delivered straight into `central` as a
    # modelled-duration input.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot + k21 * peripheral1 - (kel + k12) * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ---- Dose partitioning ---------------------------------------------------
    # An oral administration is two records at the same time: a bolus into
    # `depot` and a rate = -2 record into `central`. An intravenous
    # administration is a single plain bolus into `central`; with no rate on
    # that record rxode2 ignores dur(central) and f(central) = 1.
    f(depot) <- fdepot * ffo
    f(central) <- ROUTE_IV + (1 - ROUTE_IV) * fdepot * (1 - ffo)
    dur(central) <- d1

    # ---- Observation ---------------------------------------------------------
    # central in mg and vc in L give mg/L; x1000 reports ng/mL.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
