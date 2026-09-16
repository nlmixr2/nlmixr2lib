vanRongen_2015_midazolam_circadian <- function() {
  description <- paste(
    "Population PK model characterizing 24-hour (circadian) variation in the",
    "pharmacokinetics of oral and intravenous midazolam in 12 healthy adult",
    "male volunteers (van Rongen 2015, CPT Pharmacometrics Syst Pharmacol).",
    "Semi-simultaneous design: a 2 mg oral solution followed 150 minutes later",
    "by a 1 mg intravenous bolus, repeated at six clock times across the",
    "24-hour period (10:00, 14:00, 18:00, 22:00, 02:00, 06:00). Disposition is",
    "three-compartment with the two peripheral volumes constrained equal",
    "(Vperipheral1 = Vperipheral2 = 22.5 L); oral absorption uses a single",
    "transit compartment with the transit and absorption rate constants",
    "constrained equal (Ka = Ktr = 0.053 1/min). Three chronopharmacologic",
    "terms are carried. Oral bioavailability follows a 24-hour cosine about a",
    "mesor of 0.277 with amplitude 0.041 and acrophase 734 min after midnight",
    "(12:14), a 14.7 percent relative amplitude. Clearance follows a 24-hour",
    "cosine about a mesor of 0.379 L/min with amplitude 0.027 L/min and",
    "acrophase 1130 min after midnight (18:50), a 7.2 percent relative",
    "amplitude. The absorption rate constant is multiplied by 1.41 when the",
    "oral dose is given at 14:00. Both cosines are additive on the natural",
    "scale about an individual mesor (paper Eq. 1 and Eq. 2), so the",
    "inter-individual and inter-occasion random effects act on the mesor only",
    "and not on the amplitude. Residual error is proportional and is separate",
    "for the oral-phase (18.0 percent) and intravenous-phase (15.4 percent)",
    "records. The covariate TCLOCK carries the wall-clock hour of day at model",
    "time zero and is what makes the model's dosing time meaningful; see the",
    "covariateData notes."
  )
  reference <- paste(
    "van Rongen A, Kervezee L, Brill MJE, van Meir H, den Hartigh J,",
    "Guchelaar H-J, Meijer JH, Burggraaf J, van Oosterhout F (2015).",
    "Population Pharmacokinetic Model Characterizing 24-Hour Variation in the",
    "Pharmacokinetics of Oral and Intravenous Midazolam in Healthy Volunteers.",
    "CPT Pharmacometrics Syst Pharmacol 4(8):454-464. doi:10.1002/psp4.12007.",
    sep = " "
  )
  vignette <- "vanRongen_2015_midazolam_circadian"
  units <- list(time = "min", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    TCLOCK = list(
      description        = "Wall-clock time of day at model time zero, in decimal hours",
      units              = "h",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Anchors the model's integration axis to the wall clock so the two",
        "24-hour cosine terms can be evaluated. van Rongen 2015 Eq. 2 defines",
        "its TIME as 'the time in minutes starting at midnight of the first",
        "study visit', so the quantity the cosines consume is minutes after",
        "midnight. This model reconstructs that quantity inside model() as",
        "clockTime = TCLOCK * 60 + time, i.e. TCLOCK is the wall-clock hour at",
        "time = 0 and the integration axis supplies the elapsed minutes. It is",
        "therefore supplied time-FIXED (one value per profile) even though the",
        "wall-clock time it represents advances with the solve. Supplying it as",
        "a time-varying column instead would be wrong here: a 0-24 column wraps",
        "at midnight and rxode2's linear covariate interpolation would sweep",
        "backwards through an entire day across that wrap. Set TCLOCK to the",
        "clock hour of the oral dose that starts the profile (10, 14, 18, 22, 2",
        "or 6 for the six administration times of the study); for an",
        "intravenous-only simulation set it to the clock hour of the",
        "intravenous dose. TCLOCK is also what gates the 14:00 absorption-rate",
        "factor, because the oral dose is placed at time = 0."
      ),
      source_name        = "TIME (minutes after midnight of the first study visit; TCLOCK = TIME / 60 at time = 0)"
    ),
    OCC = list(
      description        = "Administration-occasion index for inter-occasion variability on bioavailability",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Integer 1-6 identifying which of the six midazolam administration",
        "occasions a profile belongs to. Each volunteer attended three study",
        "visits and was dosed twice per visit at a 12-hour interval, giving six",
        "occasions spanning the six clock times of administration (van Rongen",
        "2015 Methods, 'Twenty-four hour variation': the IOV represents 'the",
        "variability between the six different times of administration').",
        "Decomposed inside model() into six mutually exclusive indicators that",
        "multiplex the etaiov_fdepot_1..6 slots, per the OCC register entry.",
        "Only bioavailability retained IOV in the final model; the IOV on the",
        "absorption rate constant was removed for 55 percent eta-shrinkage and",
        "the IOV on clearance was removed as being substantially smaller than",
        "the IIV. For a single-occasion simulation set OCC = 1 throughout."
      ),
      source_name        = "OCC"
    ),
    ROUTE_IV = list(
      description        = "Intravenous-phase observation-record indicator selecting the residual error magnitude",
      units              = "(binary)",
      type               = "binary",
      reference_category = paste(
        "0 (oral-phase records). The comparator non-intravenous route here is",
        "ORAL, not subcutaneous -- van Rongen 2015 pools an oral solution arm",
        "with an intravenous bolus arm and no subcutaneous arm exists."
      ),
      notes              = paste(
        "Residual-error-only role, the same role ROUTE_IV plays in",
        "Fanta_2007_ciclosporin.R and Kuroda_2024_quinidine_horse.R: no",
        "structural parameter switches by route in this model. van Rongen 2015",
        "Table 2 reports two proportional residual errors, 'sigma oral' 18.0",
        "percent and 'sigma intravenous' 15.4 percent. The paper does not state",
        "how individual records were assigned to the two strata, and in a",
        "semi-simultaneous design the post-intravenous samples carry",
        "superimposed oral and intravenous contributions, so record assignment",
        "cannot be read off the source; this model takes the only",
        "record-level split the design admits -- ROUTE_IV = 0 on the samples",
        "drawn between the oral dose and the intravenous dose (t = 0 to 148",
        "min) and ROUTE_IV = 1 on the samples drawn after the intravenous dose",
        "(t = 155 min onward). That assumption is recorded in the vignette",
        "Errata. For a simulation of a single route, set ROUTE_IV = 0",
        "throughout for an oral-only profile and 1 throughout for an",
        "intravenous-only profile."
      ),
      source_name        = "not reported as a named column; derived from the sampling schedule (paper Methods)"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "midazolam", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    transit1 = list(
      analyte = "midazolam", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "midazolam", units = "mg",
      specimen = "serum", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "midazolam", units = "mg",
      specimen = "serum", verified = TRUE
    ),
    peripheral2 = list(
      analyte = "midazolam", units = "mg",
      specimen = "serum", verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 12L,
    n_studies      = 1L,
    age_range      = "18-27 years",
    age_median     = "22 years (mean 21.8, SD 3.19)",
    weight_range   = "63.4-92.9 kg",
    weight_median  = "75.4 kg (mean 76.0, SD 8.65)",
    sex_female_pct = 0,
    race_ethnicity = "Caucasian (100 percent by inclusion criterion)",
    disease_state  = paste(
      "Healthy, nonsmoking Caucasian male volunteers with body mass index",
      "18-30 kg/m^2 (observed 18.8-25.8, median 21.9). Excluded for any",
      "clinically significant abnormality in medical history, routine",
      "laboratory tests or 12-lead ECG, for any medication use, for extreme",
      "morning or evening chronotype on the Horne-Ostberg questionnaire, and",
      "for transmeridian flights or shift work in the month before the study."
    ),
    dose_range     = paste(
      "Semi-simultaneous administration of 2 mg oral midazolam solution",
      "followed 150 minutes later by 1 mg intravenous midazolam, given twice",
      "per study visit at a 12-hour interval across three visits, so that oral",
      "administration occurred at 10:00, 14:00, 18:00, 22:00, 02:00 and 06:00.",
      "Washout at least 2 weeks between visits. Serum samples at 0, 15, 30, 45,",
      "58, 65, 70, 75, 80, 90, 120, 148, 155, 165, 180, 210, 240, 270, 330 and",
      "390 min after the oral dose, plus 715 min on the first half of a visit.",
      "Assay LLQ 0.3 ug/L."
    ),
    regions        = "The Netherlands (Centre for Human Drug Research, Leiden)",
    notes          = paste(
      "Demographics from Table 1 of van Rongen 2015 (n = 12). One subject",
      "withdrew consent during the study for personal reasons and was replaced",
      "by another subject dosed on the same randomization order, so 13",
      "individuals were enrolled and 12 complete datasets were analyzed.",
      "Circadian entrainment was controlled and verified: subjects held a",
      "stable sleep-wake schedule for a week before each visit, wore an",
      "actigraph, remained semirecumbent, slept under dimmed lights with an eye",
      "mask from 23:30 to 07:30, and the expected 24-hour rhythms in serum TSH",
      "(29 percent relative amplitude, peak 03:05), heart rate (10 percent) and",
      "diastolic / systolic blood pressure (6.3 / 5.6 percent, peaks near",
      "16:00) were confirmed by cosinor analysis as external validators."
    )
  )

  ini({
    # Final population PK parameters from Table 2 of van Rongen 2015 (page
    # 460), 'Model estimates (RSE%)' column. The 'Bootstrap estimates (95%
    # confidence interval)' column from 250/250 successful resamples brackets
    # every point estimate below and is quoted alongside each value. The paper
    # reports concentrations in ug/L, doses in mg, volumes in L, clearances in
    # L/min and rate constants in 1/min; those units are carried unchanged.

    # ---- Clearance: 24-hour cosine, additive on the natural scale ----
    # Table 2 header row: CL = CLmesor + Amp x cos((2*pi/1440)*(Time -
    # Acrophase)), with Time in minutes after midnight. The mesor is carried
    # on the log scale per library convention; the amplitude is carried
    # linearly because it is an absolute increment in L/min that a cosine
    # drives through zero and negative values, so a log transform is not
    # meaningful for it. Relative amplitude 0.027 / 0.379 = 7.1 percent,
    # matching the 7.2 percent the paper quotes in Results and Figure 4.
    lcl <- log(0.379)
    label("Midazolam clearance mesor CLmesor (L/min)")                          # Table 2 CLmesor = 0.379 L/min (RSE 4.8%); bootstrap 0.380 (0.344-0.417)
    amp_cl <- 0.027
    label("Amplitude of the 24-hour cosine on clearance (L/min)")               # Table 2 Amp = 0.027 L/min (RSE 14.8%); bootstrap 0.028 (0.017-0.039)
    acrophase_cl <- 1130
    label("Acrophase of the 24-hour cosine on clearance (min after midnight)")  # Table 2 Acrophase = 1,130 min (RSE 2.9%); bootstrap 1,130.2 (1,005.3-1,204.7); 1130 min = 18:50, the peak time quoted in Results and Figure 4

    # ---- Disposition ----
    # Three compartments. The two peripheral volumes were constrained equal
    # because the separately estimated values were almost equal and the
    # equalized model gave a similar objective function (Results, 'Population
    # PK model'), so Table 2 prints them as the single row
    # 'Vperipheral1 = Vperipheral2'. vp2 is set from vp inside model().
    lvc <- log(18.2)
    label("Central volume of distribution Vcentral (L)")                        # Table 2 Vcentral = 18.2 L (RSE 5.4%); bootstrap 18.4 (15.3-20.9)
    lvp <- log(22.5)
    label("Peripheral volume of distribution, both peripheral compartments (L)") # Table 2 'Vperipheral1 = Vperipheral2' = 22.5 L (RSE 2.5%); bootstrap 22.4 (20.2-26.2)
    lq <- log(0.27)
    label("Inter-compartmental clearance Q, central to peripheral1 (L/min)")    # Table 2 Q = 0.27 L/min (RSE 6.8%); bootstrap 0.269 (0.209-0.334)
    lq2 <- log(1.31)
    label("Inter-compartmental clearance Q2, central to peripheral2 (L/min)")   # Table 2 Q2 = 1.31 L/min (RSE 8.5%); bootstrap 1.29 (1.08-1.56)

    # ---- Absorption ----
    # One transit compartment, with the transit rate constant Ktr and the
    # absorption rate constant Ka constrained equal (Results, 'Population PK
    # model'), so Table 2 prints them as the single row 'Ka = Ktr'. ktr is set
    # from ka inside model(), which keeps the constraint visible.
    lka <- log(0.053)
    label("Absorption and transit rate constant Ka = Ktr (1/min)")              # Table 2 'Ka = Ktr' = 0.053 1/min (RSE 5.8%); bootstrap 0.053 (0.048-0.061)

    # Multiplicative factor applied to Ka = Ktr when the oral dose is given at
    # 14:00. Named by the e_<cov>_<param> covariate-effect convention against
    # TCLOCK, the wall-clock covariate that gates it. Estimated as a fraction,
    # not a log-scale increment, so it is carried linearly and the effect is
    # written as a factor rather than an exponent in model(). A cosine on Ka
    # did not predict the 14:00 increase adequately, and a half-sine (paper
    # Eqs. 3-4, peak 14:59, amplitude 0.056 1/min, onset 14:12, offset 15:45)
    # was very sensitive to initial estimates and did not significantly improve
    # the objective function, so the multiplication-factor form below is the
    # one the paper selected.
    e_tclock_ka <- 1.41
    label("Factor on Ka = Ktr when the oral dose is given at 14:00 (unitless)") # Table 2 'Fraction Ka at 14:00' = 1.41 (RSE 4.7%); bootstrap 1.41 (1.07-1.78)

    # ---- Oral bioavailability: 24-hour cosine, additive on the natural scale ----
    # Table 2 header row: F = Fmesor + Amp x cos((2*pi/1440)*(Time -
    # Acrophase)). Relative amplitude 0.041 / 0.277 = 14.8 percent, matching
    # the 14.7 percent quoted in Results and Figure 4, and a peak-to-trough
    # relative difference of 2 x 0.041 / 0.277 = 29.6 percent against the 29.4
    # percent quoted in the Discussion.
    lfdepot <- log(0.277)
    label("Oral bioavailability mesor Fmesor (fraction)")                       # Table 2 F = 0.277 (RSE 7.1%); bootstrap 0.275 (0.244-0.313)
    amp_fdepot <- 0.041
    label("Amplitude of the 24-hour cosine on oral bioavailability (fraction)") # Table 2 Amp = 0.041 (RSE 17.3%); bootstrap 0.041 (0.026-0.055)
    acrophase_fdepot <- 734
    label("Acrophase of the 24-hour cosine on oral bioavailability (min after midnight)") # Table 2 Acrophase = 734 min (RSE 5.3%); bootstrap 739.7 (667.0-821.0); 734 min = 12:14, the peak time quoted in Results and Figure 4

    # ---- Inter-individual variability ----
    # Table 2 'Interindividual variability' block reports CV percentages for a
    # log-normally distributed random effect (Methods, 'Structural and
    # statistical model': 'Interindividual variability (IIV) in PK parameters
    # was assumed to be log-normally distributed'), so the internal variance is
    # omega^2 = log(1 + CV^2). At these magnitudes the alternative reading
    # omega = CV differs in the fourth decimal place and changes nothing.
    # IIV was retained on clearance, the absorption rate constant and
    # bioavailability only; Table 2 reports none on the volumes or on Q / Q2.
    etalcl ~ log(1 + 0.162^2)     # Table 2 interindividual variability, row 'CL (%)' = 16.2 (RSE 21); bootstrap 15.2 (9.7-19.6)
    etalka ~ log(1 + 0.191^2)     # Table 2 interindividual variability, row 'Ka (%)' = 19.1 (RSE 21.9); bootstrap 18.7 (10.7-24.2)
    etalfdepot ~ log(1 + 0.233^2) # Table 2 interindividual variability, row 'F (%)' = 23.3 (RSE 22.2); bootstrap 22.7 (15.8-28.8)

    # ---- Inter-occasion variability on bioavailability ----
    # Table 2 'Interoccasion variability' block reports a single magnitude for
    # bioavailability across the six administration occasions. Expanded into
    # six occasion-multiplexed slots because rxode2 cannot simulate the
    # `eta ~ var | OCC` form; occasions 2-6 are held at occasion 1's variance,
    # which is what NONMEM `$OMEGA BLOCK(1) SAME` encodes.
    etaiov_fdepot_1 ~ log(1 + 0.148^2)        # Table 2 interoccasion variability, row 'F (%)' = 14.8 (RSE 10.5); bootstrap 14.5 (11.5-17.9)
    etaiov_fdepot_2 ~ fixed(log(1 + 0.148^2)) # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_fdepot_3 ~ fixed(log(1 + 0.148^2)) # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_fdepot_4 ~ fixed(log(1 + 0.148^2)) # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_fdepot_5 ~ fixed(log(1 + 0.148^2)) # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_fdepot_6 ~ fixed(log(1 + 0.148^2)) # SAME-equivalent: equal to the occasion-1 IOV variance

    # ---- Residual error ----
    # Table 2 'Residual proportional error' block: separate proportional error
    # magnitudes for the oral and intravenous data (Results, 'Population PK
    # model': 'Residual variability was best described by using a proportional
    # error model for both oral and i.v. data'). Selected per record by
    # ROUTE_IV inside model().
    propSdOral <- 0.180
    label("Proportional residual SD on oral-phase records (fraction)")          # Table 2 residual proportional error, row 'sigma oral (%)' = 18.0 (RSE 5.6); bootstrap 17.8 (15.8-19.8)
    propSdIv <- 0.154
    label("Proportional residual SD on intravenous-phase records (fraction)")   # Table 2 residual proportional error, row 'sigma intravenous (%)' = 15.4 (RSE 6.1); bootstrap 15.1 (13.2-17.3)
  })

  model({
    # ---- 1. Wall-clock reconstruction ----
    # van Rongen 2015 Eq. 2 drives both cosines with TIME, 'the time in minutes
    # starting at midnight'. TCLOCK carries the wall-clock hour at time = 0 and
    # the integration axis supplies the elapsed minutes, so clockTime is the
    # paper's TIME. Both cosines have a 1,440 min period, so clockTime running
    # past 1,440 on a multi-day solve is correct and needs no wrapping.
    clockTime <- TCLOCK * 60 + time

    # Additive cosine components, paper Eq. 2. Written as separate `amp_*_t`
    # terms so neither component can be mistaken for the total quantity the
    # ODE consumes, per the parameter-names.md naming rule.
    amp_cl_t <- amp_cl *
      cos(2 * 3.141592653589793 * (clockTime - acrophase_cl) / 1440)
    amp_fdepot_t <- amp_fdepot *
      cos(2 * 3.141592653589793 * (clockTime - acrophase_fdepot) / 1440)

    # ---- 2. Occasion and dosing-time indicators ----
    # Six administration occasions multiplexing the bioavailability IOV.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 +
      oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4 +
      oc5 * etaiov_fdepot_5 + oc6 * etaiov_fdepot_6

    # The 14:00 absorption-rate factor is a step effect estimated at one
    # administration time, not a continuous function of clock time: the paper
    # tested multiplication factors at the other five administration times and
    # none further improved the model. It is therefore gated on the dosing
    # clock hour, which is TCLOCK because the oral dose is placed at time = 0.
    # A half-open half-hour window is used rather than an equality test so a
    # floating-point TCLOCK of 14 matches robustly.
    is1400 <- (TCLOCK >= 13.5) * (TCLOCK < 14.5)

    # ---- 3. Individual parameters ----
    # Paper Eq. 1 puts the random effects on the population mean,
    # theta_ij = theta_mean * exp(eta_i + kappa_ij), and paper Eq. 2 defines
    # the cosine's first term as the INDIVIDUAL mesor (the value 'around which
    # it oscillates'). The random effects therefore act on the mesor only and
    # the amplitude carries no random effect; that is why `amp_cl_t` is added
    # outside the exp() rather than multiplying it.
    cl <- exp(lcl + etalcl) + amp_cl_t
    vc <- exp(lvc)
    vp <- exp(lvp)
    vp2 <- vp
    q <- exp(lq)
    q2 <- exp(lq2)

    # Ka and Ktr are one estimated quantity (Table 2 row 'Ka = Ktr'), so the
    # 14:00 factor scales the shared constant and both absorption steps speed
    # up together. See the vignette Assumptions section.
    ka <- exp(lka + etalka) * (1 + (e_tclock_ka - 1) * is1400)
    ktr <- ka

    # ---- 4. Micro-constants ----
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # ---- 5. ODE system ----
    # Oral doses enter `depot` and pass through one transit compartment before
    # reaching `central`; intravenous doses are placed directly in `central`.
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ka * transit1
    d/dt(central) <- ka * transit1 - kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # ---- 6. Bioavailability ----
    # Applied to the oral depot only; the intravenous dose goes to `central`
    # and is unaffected. rxode2 evaluates f() at the dose event, so with the
    # oral dose at time = 0 the cosine is read at the administration clock
    # time, which is what the paper's model does.
    f(depot) <- exp(lfdepot + etalfdepot + iov_fdepot) + amp_fdepot_t

    # ---- 7. Observation and error ----
    # Amounts are in mg and volumes in L, so central/vc is mg/L; the factor of
    # 1,000 converts to the paper's ug/L serum assay scale.
    Cc <- 1000 * central / vc

    # Route-specific proportional residual error, the same residual-only
    # ROUTE_IV role as Fanta_2007_ciclosporin.R and Kuroda_2024_quinidine_horse.R.
    propSd <- propSdOral * (1 - ROUTE_IV) + propSdIv * ROUTE_IV
    Cc ~ prop(propSd)
  })
}
