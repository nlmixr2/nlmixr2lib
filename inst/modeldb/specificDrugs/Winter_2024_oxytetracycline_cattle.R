Winter_2024_oxytetracycline_cattle <- function() {
  description <- "Preclinical (cattle). Three-compartment population pharmacokinetic model for oxytetracycline in calves and adult cattle, with two parallel first-order absorption depots for intramuscular long-acting formulations (a rapid depot with Ka1 and a slow, lag-timed depot with Ka2 sharing a single ilogit-transformed bioavailability), age (calf vs adult) as a categorical covariate on all three volumes and all three clearances, and full block interindividual variability for the absorption and disposition parameter sets; meta-analysis of 1,730 plasma concentrations from 69 cattle across eight studies, used to derive VetCAST pharmacokinetic-pharmacodynamic cutoffs (Winter 2024)."
  reference   <- "Winter EA, Pelligand L, Toutain P-L, Lees P, Milanova A, Gehring R. Determination of pharmacokinetic-pharmacodynamic cutoff values of oxytetracycline in calves and adult cattle using population pharmacokinetic modeling. Front Microbiol. 2024;15:1498219. doi:10.3389/fmicb.2024.1498219"
  vignette    <- "Winter_2024_oxytetracycline_cattle"

  # The source model is written on a body-weight-normalised scale: volumes in
  # mL/kg, clearances in mL/(kg*h), and doses in micrograms/kg (the raw data
  # set that accompanies the paper codes a 20 mg/kg dose as dose_IM = 20000
  # microg/kg, and the Phoenix control stream in Supplementary Data S1 sets
  # Dose = 20000 for its secondary macro-constant block). With those units,
  # central / vc is directly in microg/mL, which is the unit of the additive
  # residual term stdev0.
  units <- list(time = "h", dosing = "ug/kg", concentration = "ug/mL")

  # buildModelDb()'s dosing heuristic only recognises states literally named
  # "depot" or "central", so the two numbered intramuscular depots must be
  # declared explicitly or the registry reports "central" alone. An
  # intramuscular administration doses depot1 and depot2 simultaneously; an
  # intravenous administration doses central.
  dosing <- c("depot1", "depot2", "central")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Winter 2024 Figure 3 (model schematic)
  # and the Phoenix control stream in Supplementary Data S1.
  compartmentData <- list(
    depot1      = list(analyte = "oxytetracycline", units = "ug/kg", specimen = "administration site", verified = TRUE),
    depot2      = list(analyte = "oxytetracycline", units = "ug/kg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "oxytetracycline", units = "ug/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "oxytetracycline", units = "ug/kg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "oxytetracycline", units = "ug/kg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CHILD = list(
      description        = "Juvenile (calf) versus adult age-cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (adult cattle)",
      notes              = "1 = calf, 0 = adult cattle. Winter 2024 Methods section 2.2.1.2: 'Age was a categorical covariate, as not all ages were known, with calves classified as animals less than 6 months old or when the authors declared the animals to be calves'; Supplementary Table S1 restates the cutoff as < 0.5 years. Coded '0' for adults (control condition) and '1' for calves in the source data sets. Applied as an exponential categorical effect on all three volumes and all three clearances (Winter 2024 Equation 1 and Table 1 footnotes a-f). This is the only covariate retained in the final model.",
      source_name        = "adult_calve"
    )
  )

  # Covariates that Winter 2024 screened but did NOT retain in the final
  # model. Documentation only -- these names are intentionally absent from
  # model(). Breed (dairy/milk vs beef/meat) and 'source' (the combined
  # formulation + analytical-method indicator, seven levels, explored on Ka2
  # only) were also screened; they have no canonical covariate-column entry
  # and are described in population$notes instead.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Cohort was 30 male / 39 female (Winter 2024 section 2.1; per-study split in Supplementary Table S1). Screened as a categorical covariate and rejected: 'Other evaluated covariates, health, sex and breed, were also tested, but did not, in any combination, improve the model' (Winter 2024 section 3.2). No point estimate is published, so no effect is carried."
    ),
    DIS_INFECT_ACTIVE = list(
      description = "Active clinical infection episode indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "14 of 69 animals were infected -- 8 calves in a severe pneumonia model and 6 cows with Trueperella pyogenes metritis (Winter 2024 section 2.1; coded code_health in the source data). Screened and rejected (Winter 2024 section 3.2); the Discussion notes the health-status variability is instead absorbed into the overall between-subject variance used for the PK/PD cutoff. The Phoenix control stream in Supplementary Data S1 carries the corresponding fixef line (dcode_health = -0.0206) commented out, confirming it is not in the final model."
    )
  )

  population <- list(
    species        = "cattle (Bos taurus); the target species of this veterinary indication, not a preclinical surrogate",
    n_subjects     = 69L,
    n_studies      = 8L,
    n_observations = 1730L,
    age_range      = "0.21 - 11 years where reported; 28 calves (< 0.5 years, or declared calves by the study authors) and 41 adult cattle (Winter 2024 section 2.1 and Supplementary Table S1)",
    weight_range   = "70.2 - 500 kg across the eight data sets (Supplementary Table S1; one data set did not report weights)",
    sex_female_pct = 56.5,
    disease_state  = "55 healthy animals; 14 infected (8 calves in a severe experimental pneumonia model, euthanised after 48 h; 6 cows with Trueperella pyogenes metritis)",
    dose_range     = "20 mg/kg intramuscularly (long-acting formulations, eight products) and 20 or 40 mg/kg intravenously",
    regions        = "Three published studies (Bulgaria, United Kingdom, Ireland/United States), two unpublished academic data sets and three unpublished pharmaceutical-company data sets",
    notes          = "Meta-analysis of eight data sets pooled by the authors. Sampling was rich (11-27 samples per animal per administration, at least 2 within the first hour and 7 within 24 h); sampling windows ran 0-48 h to 0-288 h. Twelve of 1,730 samples (0.7%) were below the limit of quantification and were discarded (Beal M1). Breed was 24 dairy / 45 beef (Friesian Holstein, Jersey, Aberdeen Angus cross, Polled Hereford and unspecified); breed, sex and health status were screened and rejected as covariates. A 'source' covariate (formulation plus analytical method, seven levels) was explored on Ka2 only, as a formulation-comparison exercise reported in Supplementary Table S4; it is NOT part of the final model and was deliberately excluded from the authors' Monte Carlo simulations so that the simulated variability reflects the diversity of EU formulations. Routes are carried in the rxode2 cmt column: 'central' for an intravenous dose; a simultaneous pair of dose records to 'depot1' and 'depot2', each carrying the full nominal dose, for an intramuscular dose (bioavailability is applied per depot by f(), exactly as the two Phoenix dosepoint() statements do). All parameters are body-weight-normalised, so amt is a per-kilogram amount in microg/kg (20 mg/kg = 20000)."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters -- Winter 2024 Table 1, "Single run Estimates"
    # column (the final-model point estimates). The bootstrap medians in the
    # same table quantify precision and are not used here. Values in Table 1
    # are the ADULT typical values; calf values follow from the covariate
    # effects below (Table 1 footnotes a-f).
    #
    # NOTE ON THE SUPPLEMENTARY CONTROL STREAM: Supplementary Data S1
    # ("Phoenix code for the final model") is the authoritative source for the
    # model STRUCTURE, but its fixef()/ranef() numbers are the initial
    # estimates carried over from the exploratory 'source'-covariate run
    # (they match Supplementary Table S4, not Table 1). Structure is taken
    # from the control stream; every value is taken from Tables 1 and 2.
    # ------------------------------------------------------------------
    lvc  <- log(126)     ; label("Central volume of distribution, adult (V1, mL/kg)")                              # Winter 2024 Table 1 tvV1 (adult) = 126 mL/kg
    lvp  <- log(914)     ; label("First peripheral volume of distribution, adult (V2, mL/kg)")                     # Winter 2024 Table 1 tvV2 (adult) = 914 mL/kg
    lvp2 <- log(2564)    ; label("Second peripheral volume of distribution, adult (V3, mL/kg)")                    # Winter 2024 Table 1 tvV3 (adult) = 2,564 mL/kg
    lcl  <- log(62.67)   ; label("Plasma clearance, adult (Cl, mL/kg/h)")                                          # Winter 2024 Table 1 tvCl (adult) = 62.67 mL/kg/h
    lq   <- log(485)     ; label("Inter-compartmental clearance central <-> peripheral1, adult (Cld2, mL/kg/h)")   # Winter 2024 Table 1 tvCld2 (adult) = 485 mL/kg/h
    lq2  <- log(19.96)   ; label("Inter-compartmental clearance central <-> peripheral2, adult (Cld3, mL/kg/h)")   # Winter 2024 Table 1 tvCld3 (adult) = 19.96 mL/kg/h

    # Absorption after intramuscular long-acting product. Two parallel
    # first-order sites: a rapid one (Ka1, no lag) and a slow one (Ka2,
    # entered after Tlag). Winter 2024 Figure 3 and Supplementary Data S1
    # dosepoint(Abs1, ...) / dosepoint(Abs2, tlag = (Tlag), ...).
    lka1  <- log(0.214)  ; label("Rapid first-order absorption rate constant (Ka1, 1/h)")                          # Winter 2024 Table 1 tvKa1 = 0.214 1/h (mean absorption time 4.57 h)
    lka2  <- log(0.0441) ; label("Slow first-order absorption rate constant (Ka2, 1/h)")                           # Winter 2024 Table 1 tvKa2 = 0.0441 1/h (mean absorption time 22.1 h)
    ltlag <- log(14.96)  ; label("Lag time before the slow (Ka2) absorption site starts releasing (Tlag, h)")      # Winter 2024 Table 1 tvTlag for ka2 = 14.96 h

    # Bioavailability. Supplementary Data S1 is explicit that only the TOTAL
    # bioavailability is ilogit-transformed --
    #   stparm(F1     = ilogit(tvF1 + nF1))
    #   stparm(Frapid = tvFrapid * exp(nFrapid))
    # -- so tvFrapid is a plain fraction with exponential IIV, whereas tvF1
    # lives on the ilogit scale. Winter 2024 Methods 2.2.1.1 confirms the
    # transformation was introduced for "bioavailability (F)" specifically,
    # "to prevent the bioavailability being estimated greater than 100%".
    # See the vignette Errata for the Table 1 "Units = ilogit" annotation on
    # the tvFrapid row, which the control stream contradicts.
    lffo        <- log(0.757) ; label("Fraction of the bioavailable dose entering the rapid (Ka1) absorption site (Frapid, fraction)")   # Winter 2024 Table 1 tvFrapid = 0.757
    logitfdepot <- 1.281      ; label("Total intramuscular bioavailability on the ilogit scale (F1; expit(1.281) = 0.783)")              # Winter 2024 Table 1 tvF1 = 1.281 (ilogit scale; the table annotates the back-transform as a total F of 78.6%)

    # ------------------------------------------------------------------
    # Covariate effects -- age (calf vs adult), Winter 2024 Equation 1:
    #   P = tvP * exp(dadult_calfP * CALF) * exp(nP)
    # with adults coded 0 (reference) and calves coded 1. All six effects
    # were significantly different from 0 (Winter 2024 section 3.2).
    # ------------------------------------------------------------------
    e_child_vc  <- 0.320 ; label("Effect of calf age cohort on V1, exponential scale (unitless)")     # Winter 2024 Table 1 dadult_calfV1 = 0.320 (footnote d: V1 calves = 126 * exp(0.320) = 174 mL/kg)
    e_child_vp  <- 0.159 ; label("Effect of calf age cohort on V2, exponential scale (unitless)")     # Winter 2024 Table 1 dadult_calfV2 = 0.159 (footnote e: V2 calves = 914 * exp(0.159) = 1,071 mL/kg)
    e_child_vp2 <- 0.358 ; label("Effect of calf age cohort on V3, exponential scale (unitless)")     # Winter 2024 Table 1 dadult_calfV3 = 0.358 (footnote f: V3 calves = 2,564 * exp(0.358) = 3,669 mL/kg)
    e_child_cl  <- 0.548 ; label("Effect of calf age cohort on Cl, exponential scale (unitless)")     # Winter 2024 Table 1 dadult_calfCl = 0.548 (footnote a: Cl calves = 62.67 * exp(0.548) = 108.4 mL/kg/h)
    e_child_q   <- 0.190 ; label("Effect of calf age cohort on Cld2, exponential scale (unitless)")   # Winter 2024 Table 1 dadult_calfCld2 = 0.190 (footnote b: Cld2 calves = 485 * exp(0.190) = 586 mL/kg/h)
    e_child_q2  <- 0.293 ; label("Effect of calf age cohort on Cld3, exponential scale (unitless)")   # Winter 2024 Table 1 dadult_calfCld3 = 0.293 (footnote c: Cld3 calves = 19.96 * exp(0.293) = 26.76 mL/kg/h)

    # ------------------------------------------------------------------
    # Interindividual variability -- Winter 2024 Table 2 "Omega" block.
    # Two full variance/covariance matrices were estimated: one for the
    # absorption and bioavailability parameters, one for the disposition
    # parameters. The tabulated numbers are variances and covariances on the
    # estimation scale (log for every parameter except F1, which is on the
    # ilogit scale). Reproducing the published BSV% column confirms the
    # variance reading: 100 * sqrt(exp(diag) - 1) returns 28.87, 45.76,
    # 62.64, 20.34, 12.68 and 83.34, 17.75, 33.71, 20.14, 22.13, 46.89, i.e.
    # every one of Table 2's eleven BSV% entries. Both blocks are positive
    # definite as printed.
    #
    # Block 1 (Table 2 order nKa1, nKa2, nF1, nTlag, nFrapid). Eta shrinkage
    # 16.7%, 20.0%, 25.6%, 43.9%, 18.2%.
    # ------------------------------------------------------------------
    etalka1 + etalka2 + etalogitfdepot + etaltlag + etalffo ~
      c(0.080095,
        0.103949, 0.19011,
        -0.095537, -0.138489, 0.331033,
        -0.02741, -0.02635, 0.068231, 0.040549,
        0.02069, 0.038294, -0.05192, -0.01144, 0.015948)

    # Block 2 (Table 2 order nV1, nV2, nV3, nCl, nCl2, nCl3). Eta shrinkage
    # 12.8%, 26.0%, 44.6%, 15.7%, 41.1%, 21.4%.
    etalvc + etalvp + etalvp2 + etalcl + etalq + etalq2 ~
      c(0.527444,
        -0.07716, 0.031032,
        -0.07185, 0.02222, 0.107616,
        0.055107, -0.00525, 0.018945, 0.039763,
        0.010349, 0.013185, 0.041094, 0.016754, 0.047794,
        0.080094, 0.009562, 0.055244, 0.045213, 0.029476, 0.19877)

    # ------------------------------------------------------------------
    # Residual error. Phoenix "additive plus multiplicative, mix-ratio"
    # parameterisation (Supplementary Data S1):
    #   error(CEps = stdev0)
    #   observe(CObs = C + CEps * sqrt(1 + C^2 * (CMultStdev / sigma())^2))
    # Because sigma() returns the standard deviation of CEps (= stdev0), the
    # residual variance collapses to stdev0^2 + (C * CMultStdev)^2, which is
    # exactly nlmixr2's add() + prop() combination.
    # ------------------------------------------------------------------
    propSd <- 0.182  ; label("Proportional residual error (fraction)")     # Winter 2024 Table 1 tvCMultStdev = 18.2%
    addSd  <- 0.0069 ; label("Additive residual error SD (ug/mL)")         # Winter 2024 Table 1 stdev0 = 0.0069 ug/mL
  })

  model({
    # Individual parameters. Age enters as an exponential categorical shift
    # on the log scale, which is identical to the source's multiplicative
    # form P = tvP * exp(d * CHILD) * exp(eta) (Winter 2024 Equation 1;
    # Supplementary Data S1 stparm lines).
    vc   <- exp(lvc  + e_child_vc  * CHILD + etalvc)
    vp   <- exp(lvp  + e_child_vp  * CHILD + etalvp)
    vp2  <- exp(lvp2 + e_child_vp2 * CHILD + etalvp2)
    cl   <- exp(lcl  + e_child_cl  * CHILD + etalcl)
    q    <- exp(lq   + e_child_q   * CHILD + etalq)
    q2   <- exp(lq2  + e_child_q2  * CHILD + etalq2)
    ka1  <- exp(lka1  + etalka1)
    ka2  <- exp(lka2  + etalka2)
    tlag <- exp(ltlag + etaltlag)

    # Frapid is a plain fraction with exponential IIV; F1 is on the ilogit
    # scale (see the ini() note and the vignette Errata).
    ffo    <- exp(lffo + etalffo)
    fdepot <- expit(logitfdepot + etalogitfdepot)

    # Three-compartment disposition with two parallel absorption depots.
    # Supplementary Data S1:
    #   deriv(A1 = -(Cl*C) - (Cl2*(C - C2)) - (Cl3*(C - C3)) + Abs1*Ka1 + Abs2*Ka2)
    #   deriv(A2 = Cl2*(C - C2)); deriv(A3 = Cl3*(C - C3))
    #   deriv(Abs1 = -(Abs1*Ka1)); deriv(Abs2 = -(Abs2*Ka2))
    d/dt(depot1)      <- -ka1 * depot1
    d/dt(depot2)      <- -ka2 * depot2
    d/dt(central)     <- ka1 * depot1 + ka2 * depot2 -
      cl * central / vc -
      q  * (central / vc - peripheral1 / vp) -
      q2 * (central / vc - peripheral2 / vp2)
    d/dt(peripheral1) <- q  * (central / vc - peripheral1 / vp)
    d/dt(peripheral2) <- q2 * (central / vc - peripheral2 / vp2)

    # Both intramuscular depots receive the full nominal dose; the split of
    # the bioavailable fraction between them is applied here, matching
    #   dosepoint(Abs1,             bioavail = (Frapid * F1), ...)
    #   dosepoint(Abs2, tlag = Tlag, bioavail = ((1 - Frapid) * F1), ...)
    # An intravenous dose goes to central, which carries no f().
    f(depot1)    <- ffo * fdepot
    f(depot2)    <- (1 - ffo) * fdepot
    alag(depot2) <- tlag

    # Amounts are in ug/kg and vc is in mL/kg, so Cc is in ug/mL (= mg/L)
    # directly, matching Supplementary Data S1 (C = A1 / V1 with Dose = 20000
    # for a 20 mg/kg administration).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
