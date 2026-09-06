Heida_2024_mycophenolic_acid <- function() {
  description <- paste(
    "Two-compartment population PK model for mycophenolic acid (MPA) in",
    "pediatric kidney transplant recipients receiving oral mycophenolate",
    "mofetil (CellCept) with tacrolimus or everolimus co-treatment (Heida",
    "2024, Radboudumc model-development cohort). Absorption is Erlang-type:",
    "the dose lands in a depot that drains through one transition",
    "compartment into the central compartment, both steps sharing the same",
    "first-order rate constant ktr (1.48 1/h), so the absorption delay is",
    "Erlang-distributed with shape 2. Apparent clearance CL/F is 16.0 L/h,",
    "apparent central volume Vc/F 24.9 L, apparent peripheral volume Vp/F",
    "1590 L and apparent intercompartmental clearance Q/F 36.2 L/h. All",
    "disposition parameters are allometrically scaled on total body weight",
    "to a 70 kg reference with exponents fixed at 0.75 for flows, 1 for",
    "volumes and -0.25 for the rate constant. Serum albumin lowers apparent",
    "clearance through a power term (ALB/34)^-2.49: low albumin raises the",
    "unbound fraction of this highly protein-bound drug, so total-",
    "concentration apparent clearance rises. Random effects are",
    "between-subject variability on CL/F (38.6% CV), Vc/F (320% CV) and Q/F",
    "(63.6% CV) plus inter-occasion variability on relative bioavailability",
    "(46.1% CV); residual error is proportional (47.3%).",
    "IMPORTANT UNIT CONVENTION: doses are in milligrams of mycophenolate",
    "mofetil (the prodrug) while concentrations are milligrams per litre of",
    "mycophenolic acid, so CL/F and the volumes are apparent values",
    "referenced to the MMF dose and already absorb the 0.739 MMF-to-MPA",
    "molecular-weight ratio."
  )
  reference <- paste(
    "Heida A, Jager NGL, Aarnoutse RE, de Winter BCM, de Jong H, Keizer RJ,",
    "Cornelissen EAM, ter Heine R (2024).",
    "Model-informed dose optimization of mycophenolic acid in pediatric",
    "kidney transplant patients.",
    "Eur J Clin Pharmacol 80(11):1761-1771. doi:10.1007/s00228-024-03743-0",
    sep = " "
  )
  vignette <- "Heida_2024_mycophenolic_acid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric size descriptor for every structural parameter, against",
        "a fixed 70 kg reference (Heida 2024 Methods, Model development:",
        "'all volume and flow parameters were allometrically scaled to a",
        "total body weight of 70 kg'; Eqs. 1-5; supplementary control stream",
        "$PK ALLOCL = (WT/70)**0.75, ALLOV = (WT/70), ALLOK =",
        "(WT/70)**(-0.25)). Exponents were fixed, not estimated: 0.75 for",
        "the flow parameters CL/F and Q/F, 1 for the volumes Vc/F and Vp/F,",
        "and -0.25 for the absorption rate constant ktr. Note the 70 kg",
        "reference is a convention, not a cohort statistic -- the cohort",
        "median weight was 38.5 kg (range 12.9-79.9, Table 1), so a typical",
        "patient sits well below the reference and the scaling is an",
        "extrapolation downward."
      ),
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Serum albumin concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only estimated covariate effect in the model: a power term on",
        "apparent clearance, (ALB/34)^-2.49, normalised to 34 g/L. The paper",
        "reports the albumin normalisation constant inconsistently -- Methods",
        "states the covariate 'was normalized by dividing the value by",
        "35 g/L' and Eq. 1 prints (albumin/35)^-2.5, while Table 2 prints",
        "'(Albumin/34)^theta2' with theta2 = -2.49 and the supplementary",
        "NONMEM control stream computes COVALB = (ALB/34)**THETA(6) with",
        "THETA(6) = -2.49. The 34 g/L / -2.49 pairing is used here because",
        "the executable control stream and the results table agree on it and",
        "34 g/L is the cohort median albumin (Table 1); Eq. 1 additionally",
        "rounds the exponent to -2.5, marking it as a restatement rather",
        "than the fitted form. See the vignette Errata.",
        "The paper reports the covariate is already in canonical SI g/L, so",
        "no unit conversion is applied. Direction: albumin below the",
        "reference INCREASES apparent clearance, because MPA is highly",
        "protein bound and a lower albumin raises the unbound fraction, so",
        "the measured total concentration falls (Discussion). Observed",
        "cohort range 24-42 g/L (Table 1); the exponent magnitude of 2.49",
        "makes this term steep, so extrapolating far outside that range is",
        "not supported."
      ),
      source_name        = "ALB"
    ),
    OCC = list(
      description        = "Sampling-occasion index used for the inter-occasion random effect on relative bioavailability",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Ten occasion slots, matching the supplementary control stream's",
        "IF (OCC.EQ.1) ... IF (OCC.EQ.10) multiplexer over ETA(6)-ETA(15),",
        "each drawn from a single shared variance ($OMEGA BLOCK(1) 0.19",
        "followed by nine SAME repeats). The observed cohort contributed a",
        "median of 2 occasions per patient (range 1-6, Table 1), so slots 7",
        "to 10 were allocated but never exercised by the model-development",
        "data. For simulation, set OCC to the index of each dosing occasion;",
        "a single-occasion simulation may use OCC = 1 throughout. Records",
        "with OCC outside 1-10 zero out every indicator and therefore carry",
        "no inter-occasion effect."
      ),
      source_name        = "OCC"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "mycophenolate mofetil", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    transit1 = list(
      analyte = "mycophenolate mofetil", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "mycophenolic acid", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "mycophenolic acid", units = "mg",
      specimen = "plasma", verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 30L,
    n_studies      = 1L,
    age_range      = "4-18 years",
    age_median     = "13 years",
    weight_range   = "12.9-79.9 kg",
    weight_median  = "38.5 kg",
    sex_female_pct = 40,
    race_ethnicity = "Not reported; single-centre Dutch cohort",
    disease_state  = paste(
      "Pediatric kidney transplant recipients on maintenance",
      "immunosuppression. Median serum albumin 34 g/L (range 24-42), median",
      "height 149 cm (range 95-193), median body surface area 1.3 m2 (range",
      "0.58-2.1). Sampling was mostly early after transplantation: median",
      "post-transplant time 9.5 days (range 2-3058). Immunosuppressive",
      "co-medication was tacrolimus in 25 patients (83.3%) and everolimus in",
      "5 (16.7%), with prednisone in 13.3%; neither calcineurin/mTOR partner",
      "is known to affect MMF PK. Patients on ciclosporin were excluded",
      "because it inhibits MPA enterohepatic recirculation, so the model",
      "should not be extrapolated to ciclosporin co-treatment."
    ),
    dose_range     = paste(
      "Oral mycophenolate mofetil (CellCept) 1200 mg/m2 per day for the",
      "first 2 weeks, reduced to 600 mg/m2 per day thereafter, divided over",
      "two daily doses, with subsequent TDM-guided adjustment. Observed",
      "daily dose median 1000 mg (range 500-2000)."
    ),
    regions        = "Netherlands (Amalia Children's Hospital, Radboudumc, Nijmegen)",
    notes          = paste(
      "Retrospective routine-care therapeutic-drug-monitoring data collected",
      "June 2016 to April 2023. 266 MPA plasma concentrations: 20 full PK",
      "curves (approximately 8 samples at 0, 1, 2, 3, 4, 6, 8 and 12 h), 24",
      "limited sampling curves (3 samples at 0, 0.5 and 2 h) and 25 trough",
      "levels; median 9.5 observations per patient (range 3-18). MPA was",
      "assayed by a validated EMIT immunoassay (Roche, Cobas c502), which is",
      "known to over-read MPA relative to HPLC through cross-reactivity with",
      "the acyl-glucuronide metabolite -- so the model is calibrated to",
      "EMIT-scale concentrations (Discussion). Baseline characteristics are",
      "Table 1. The model was externally evaluated in 18 further children",
      "from Erasmus MC (Sophia Children's Hospital, Rotterdam); that external",
      "cohort was not used for estimation and is not described by this",
      "population block."
    )
  )

  ini({
    # --- Structural parameters. Final estimates are Heida 2024 Table 2
    # (95% CI by sampling importance resampling), corroborated value for
    # value by the supplementary NONMEM control stream $THETA block, whose
    # six entries reproduce Table 2 exactly. Values are typical for a 70 kg
    # patient at the 34 g/L reference albumin.
    #
    # The abstract prints wider/garbled confidence intervals for Vc/F
    # ("93.0-6.71E25") and Q/F ("9.63-74.7") than Table 2 does
    # ("6.53-45.8" and "25.8-49.6"); the point estimates agree everywhere
    # and only Table 2's intervals are internally consistent. See the
    # vignette Errata.
    lcl <- log(16.0)
    label("Apparent oral clearance CL/F at 70 kg and 34 g/L albumin (L/h)")            # Table 2 theta1 16.0 L/h (95% CI 10.3-20.4); Eq. 1; control stream $THETA 1 "16.0 ;Cl"
    lvc <- log(24.9)
    label("Apparent central volume of distribution Vc/F at 70 kg (L)")                 # Table 2 Vc/F 24.9 L (95% CI 6.53-45.8); Eq. 2; control stream $THETA 2 "24.9 ;V3"
    lvp <- log(1590)
    label("Apparent peripheral volume of distribution Vp/F at 70 kg (L)")              # Table 2 Vp/F 1590 L (95% CI 651-2994); Eq. 3; control stream $THETA 3 "1590 ;V4"
    lq <- log(36.2)
    label("Apparent intercompartmental clearance Q/F at 70 kg (L/h)")                  # Table 2 Q/F 36.2 L/h (95% CI 25.8-49.6); Eq. 4; control stream $THETA 4 "36.2;Q"
    lktr <- log(1.48)
    label("Erlang absorption transit rate constant ktr at 70 kg (1/h)")                # Table 2 ktr 1.48 1/h (95% CI 1.15-1.84); Eq. 5; control stream $THETA 5 "1.48 ; KA"
    lfdepot <- fixed(log(1))
    label("Relative oral bioavailability F (fraction)")                                # Control stream $PK "F1=1*EXP(IOV)" -- typical value fixed at 1, carrying inter-occasion variability only; no absolute bioavailability data were available (Methods, Model development)

    # --- Allometric exponents, all FIXED rather than estimated (Methods,
    # Model development: "The allometric coefficients were fixed at 0.75 for
    # flow parameters and 1 for volume parameters and, consequently, at
    # -0.25 for rate constants"). The control stream builds three shared
    # multipliers, ALLOCL / ALLOV / ALLOK, which are expanded here to one
    # named exponent per parameter so each scaling is auditable in isolation.
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of body weight on CL/F (unitless)")                     # Methods, Model development; Eq. 1 (weight/70)^0.75; control stream $PK ALLOCL = (WT/70)**0.75
    e_wt_q <- fixed(0.75)
    label("Allometric exponent of body weight on Q/F (unitless)")                      # Methods, Model development; Eq. 4 (weight/70)^0.75; control stream $PK Q = THETA(4)*ALLOCL*...
    e_wt_vc <- fixed(1)
    label("Allometric exponent of body weight on Vc/F (unitless)")                     # Methods, Model development; Eq. 2 (weight/70); control stream $PK ALLOV = (WT/70)
    e_wt_vp <- fixed(1)
    label("Allometric exponent of body weight on Vp/F (unitless)")                     # Methods, Model development; Eq. 3 (weight/70); control stream $PK V3 = THETA(3)*ALLOV*...
    e_wt_ktr <- fixed(-0.25)
    label("Allometric exponent of body weight on the transit rate constant ktr (unitless)")  # Methods, Model development; Eq. 5 (weight/70)^-0.25; control stream $PK ALLOK = (WT/70)**(-0.25)

    # --- Covariate effect. The single estimated covariate relationship:
    # a power term of serum albumin on apparent clearance, normalised to
    # 34 g/L. Adding it dropped the objective function by 15.29 points
    # (Results, Model development). Negative exponent: below-reference
    # albumin raises apparent clearance of total MPA.
    e_alb_cl <- -2.49
    label("Power exponent on (ALB/34) for CL/F (unitless)")                            # Table 2 theta2 -2.49 (95% CI -3.82 to -1.41); control stream $PK COVALB = (ALB/34)**THETA(6), $THETA 6 "-2.49; ALB~CL". Methods and Eq. 1 instead print a 35 g/L normalisation with the exponent rounded to -2.5; see covariateData$ALB$notes and the vignette Errata.

    # --- Random effects. Heida 2024 Table 2 reports every variability term
    # as a percent CV, and Methods gives the transform explicitly:
    # CV(%) = sqrt(exp(omega^2) - 1) * 100. The variances below invert that
    # formula on the Table 2 percentages. The supplementary control stream
    # $OMEGA block lists 0.139, 2.42, 0.337 and 0.19 for the same four
    # terms; the first two match the Table exactly, the last two differ in
    # the third decimal (0.337 -> 63.3% vs the printed 63.6%; 0.19 -> 45.8%
    # vs the printed 46.1%). Table 2 governs here because it is the
    # peer-reviewed result carrying confidence intervals and shrinkage, and
    # because it is what a reader checks a reproduction against. See the
    # vignette Errata for the side-by-side.
    #
    # Only CL/F, Vc/F and Q/F carry between-subject variability. The control
    # stream fixes the Vp/F and ktr etas to zero variance ($OMEGA "0 FIX ;
    # 3 IIV V4" and "0 FIX; 6 IVV KA"), matching Table 2, which lists no
    # inter-individual variability rows for those two parameters. They are
    # therefore omitted entirely rather than declared at zero.
    etalcl ~ 0.13889
    label("Between-subject variability in CL/F (log-scale variance)")                  # Table 2 IIV Cl/F 38.6% CV (95% CI 15.38-73.95, shrinkage 34.6%); log(0.386^2 + 1); control stream $OMEGA 1 "0.139"
    etalvc ~ 2.41955
    label("Between-subject variability in Vc/F (log-scale variance)")                  # Table 2 IIV Vc/F 320% CV (95% CI 109.65-18819.64, shrinkage 31.6%); log(3.20^2 + 1); control stream $OMEGA 2 "2.42"
    etalq ~ 0.33997
    label("Between-subject variability in Q/F (log-scale variance)")                   # Table 2 IIV Q/F 63.6% CV (95% CI 31.75-97.25, shrinkage 28.0%); log(0.636^2 + 1); control stream $OMEGA 4 "0.337"

    # Inter-occasion variability on relative bioavailability, over the ten
    # occasion slots of the control stream's ETA(6)-ETA(15) multiplexer. The
    # source declares one estimated variance followed by nine
    # "$OMEGA BLOCK(1) SAME" repeats; nlmixr2 has no SAME shortcut, so
    # occasions 2-10 are fix()-pinned to the occasion-1 value (the
    # Jonsson_2011_ethambutol / Oosten_2016_fentanyl pattern).
    etaiov_fdepot_1 ~ 0.19267
    label("Inter-occasion variability in relative bioavailability, occasion 1 (log-scale variance)")  # Table 2 IOV F 46.1% CV (95% CI 33.79-60.86, shrinkage 24.5%); log(0.461^2 + 1); control stream $OMEGA BLOCK(1) "0.19  ; OCC 1"
    etaiov_fdepot_2 ~ fix(0.19267)                                                     # $OMEGA BLOCK(1) SAME ; OCC 2
    etaiov_fdepot_3 ~ fix(0.19267)                                                     # $OMEGA BLOCK(1) SAME ; OCC 3
    etaiov_fdepot_4 ~ fix(0.19267)                                                     # $OMEGA BLOCK(1) SAME ; OCC 4
    etaiov_fdepot_5 ~ fix(0.19267)                                                     # $OMEGA BLOCK(1) SAME ; OCC 5
    etaiov_fdepot_6 ~ fix(0.19267)                                                     # $OMEGA BLOCK(1) SAME ; OCC 6
    etaiov_fdepot_7 ~ fix(0.19267)                                                     # $OMEGA BLOCK(1) SAME ; OCC 7
    etaiov_fdepot_8 ~ fix(0.19267)                                                     # $OMEGA BLOCK(1) SAME ; OCC 8
    etaiov_fdepot_9 ~ fix(0.19267)                                                     # $OMEGA BLOCK(1) SAME ; OCC 9
    etaiov_fdepot_10 ~ fix(0.19267)                                                    # $OMEGA BLOCK(1) SAME ; OCC 10

    # --- Residual error. Purely proportional (Results, Model development:
    # "Residual error was best modelled by a proportional error"; control
    # stream $ERROR "Y=IPRED+IPRED*ERR(1)"). Table 2's 47.3% is the residual
    # standard deviation, not a lognormal CV: the control stream $SIGMA of
    # 0.223 gives sqrt(0.223) = 47.2%, whereas the lognormal transform used
    # for the IIV rows would give 50.0%.
    propSd <- 0.473
    label("Proportional residual error (fraction)")                                    # Table 2 residual error 47.3% CV (95% CI 42.47-52.68, shrinkage 11.7%); control stream $SIGMA "0.223;PROP ERR", sqrt(0.223) = 0.472
  })

  model({
    # 1. Occasion indicators. The control stream multiplexes the
    # inter-occasion eta on relative bioavailability with
    # IF (OCC.EQ.k) IOV=ETA(5+k) over ten occasions.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)
    oc7 <- (OCC == 7)
    oc8 <- (OCC == 8)
    oc9 <- (OCC == 9)
    oc10 <- (OCC == 10)

    iov_fdepot <-
      oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 +
      oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4 +
      oc5 * etaiov_fdepot_5 + oc6 * etaiov_fdepot_6 +
      oc7 * etaiov_fdepot_7 + oc8 * etaiov_fdepot_8 +
      oc9 * etaiov_fdepot_9 + oc10 * etaiov_fdepot_10

    # 2. Individual parameters, Heida 2024 Eqs. 1-5. Every structural
    # parameter is allometrically scaled to the 70 kg reference; apparent
    # clearance additionally carries the albumin power term.
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (ALB / 34)^e_alb_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vp <- exp(lvp) * (WT / 70)^e_wt_vp
    q <- exp(lq + etalq) * (WT / 70)^e_wt_q
    ktr <- exp(lktr) * (WT / 70)^e_wt_ktr
    fdepot <- exp(lfdepot + iov_fdepot)

    # 3. Micro-constants. The source is an ADVAN5 general-linear model whose
    # $PK block writes K30 = CL/V3, K34 = Q/V3 and K43 = Q/V4, i.e.
    # elimination and distribution out of the central compartment are scaled
    # by the central volume and the return flow by the peripheral volume.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system. Erlang-type absorption with a single transition
    # compartment: the control stream's $MODEL declares COMP=(DOSE),
    # COMP=(TRAN), COMP=(CENTRAL), COMP=(PERIPHERAL) and sets K12 = K23 =
    # KTR, so the dose passes through two sequential first-order steps of
    # equal rate before reaching plasma. The resulting absorption-time
    # distribution is Erlang with shape 2 and rate ktr; mean absorption time
    # is 2/ktr = 1.35 h at the 70 kg reference.
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(central) <- ktr * transit1 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 5. Bioavailability, applied to the dosing compartment as in the
    # control stream's F1 = 1*EXP(IOV).
    f(depot) <- fdepot

    # 6. Observation and error. Doses enter in milligrams of mycophenolate
    # mofetil; Cc is the mycophenolic acid plasma concentration in mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
