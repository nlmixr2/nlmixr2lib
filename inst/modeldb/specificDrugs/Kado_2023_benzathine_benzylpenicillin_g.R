Kado_2023_benzathine_benzylpenicillin_g <- function() {
  description <- "One-compartment population PK model for penicillin released from a high-dose subcutaneous infusion of benzathine benzylpenicillin G (SCIP) with a fixed elimination rate and a dual-input depot: part of the dose enters the depot as a bolus and the remainder by a zero-order process of estimated duration (DUR, ~44 days, scaled by body mass index), after which absorption proceeds through a single transit compartment (t1/2,tr) and then a slow first-order step (t1/2,abs ~ 11.8 days) into the central compartment. Developed from 400 dried-blood-spot penicillin concentrations collected over 16 weeks in 24 healthy adult volunteers given a single 3.6, 7.2, or 10.8 MIU subcutaneous abdominal infusion (Kado 2023)."
  reference   <- paste(
    "Kado J, Salman S, Hla TK, Enkel S, Henderson R, Hand RM, Hort A,",
    "Page-Sharp M, Batty K, Moore BR, Bennett J, Anderson A, Carapetis J,",
    "Manning L. Subcutaneous infusion of high-dose benzathine penicillin G is",
    "safe, tolerable, and suitable for less-frequent dosing for rheumatic heart",
    "disease secondary prophylaxis: a phase 1 open-label population",
    "pharmacokinetic study. Antimicrob Agents Chemother. 2023;67(12):e00962-23.",
    "doi:10.1128/aac.00962-23.",
    "Structural approach (fixed elimination rate, absorption-limited kinetics)",
    "follows the same group's earlier analysis; see",
    "modellib('Kado_2020_benzathine_benzylpenicillin_g').",
    sep = " "
  )
  vignette    <- "Kado_2023_benzathine_benzylpenicillin_g"
  units       <- list(time = "day", dosing = "mg", concentration = "ng/mL")

  # The final model was fitted to dried-blood-spot (whole blood) concentrations
  # (Results, 'Pharmacokinetic modeling'); plasma and DBS were shown to be
  # equivalent (bootstrap ratio 0.990, 95% CI 0.897-1.08), so the DBS fit is
  # interchangeable with a plasma fit. `depot` and `depot2` are the bolus and
  # zero-order halves of the single subcutaneous depot of Figure 2.
  compartmentData <- list(
    depot    = list(analyte = "benzathine benzylpenicillin g", units = "mg", specimen = "administration site", verified = FALSE),
    depot2   = list(analyte = "benzathine benzylpenicillin g", units = "mg", specimen = "administration site", verified = FALSE),
    transit1 = list(analyte = "benzathine benzylpenicillin g", units = "mg", specimen = "administration site", verified = FALSE),
    central  = list(analyte = "benzathine benzylpenicillin g", units = "mg", specimen = "whole blood", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight, used as the allometric size descriptor.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reference weight 70 kg: Kado 2023 Table 2 reports kel as 'h^-1 . 70 kg^-1'",
        "and V/F as 'L . 70 kg^-1'. Allometric exponents were fixed a priori at 1 for",
        "volume and 3/4 for clearance (Methods, 'Population pharmacokinetic analysis'),",
        "so kel = CL/V scales with exponent 3/4 - 1 = -1/4. WT was selected over",
        "adjusted body weight, fat-free mass, and ideal body weight as the size",
        "descriptor (Results, 'Pharmacokinetic modeling'); the ABW fit was marginally",
        "better (dOFV -0.696) but WT was retained for comparability with the group's",
        "earlier publications. Individual body weights are not tabulated in the paper;",
        "only BMI is reported. Absorption parameters (DUR, t1/2,abs, t1/2,tr, RATIO)",
        "are NOT allometrically scaled.",
        sep = " "
      ),
      source_name        = "WT"
    ),
    BMI = list(
      description        = "Body mass index at baseline; power covariate on the zero-order input duration DUR.",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Applied as the power form DUR_i = DUR_pop * (BMI / 25)^1.13, matching the",
        "general continuous-covariate power model given in the Supplemental",
        "Information ('Population pharmacokinetic analysis'):",
        "theta_i = theta_POP * (x_i / x_tilde)^theta_POWER, where x_tilde is the",
        "population average covariate value. The paper does NOT print the numeric",
        "value of x_tilde. It is reconstructed as 25 kg/m^2: substituting the study",
        "BMI range (21.9-34.0 kg/m^2, Results 'Participants') into the power form with",
        "DUR_pop = 1062 h and exponent 1.13 returns 38.1 and 62.6 days, reproducing the",
        "paper's own statement that BMI shifts DUR over 'a range of 38-63 days for the",
        "lowest and highest BMI in the study' (Results, 'Pharmacokinetic modeling').",
        "25 kg/m^2 is also the study's ideal-vs-higher BMI group boundary (Methods,",
        "'Participants'). Using the reported median BMI of 25.1 instead changes DUR by",
        "under 0.5%. BMI was chosen over waist circumference, hip circumference, and",
        "ultrasonographic subcutaneous fat, which were also significant (dOFV -6.105).",
        "Assumed time-fixed at baseline.",
        sep = " "
      ),
      source_name        = "BMI"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 24L,
    n_studies      = 1L,
    age_range      = "18.0-54.1 years",
    age_median     = "26.9 years",
    weight_range   = "not reported (BMI 21.9-34.0 kg/m^2; enrolment criterion 20.0-34.9 kg/m^2)",
    weight_median  = "not reported (BMI median 25.1 kg/m^2)",
    sex_female_pct = 16.7,
    race_ethnicity = c(European = 33.3, Asian = 33.3, HispanicLatino = 16.7, African = 8.3, Mixed = 8.3),
    disease_state  = "healthy adult volunteers without rheumatic heart disease (no chronic renal impairment or significant liver dysfunction, no penicillin/cephalosporin allergy)",
    dose_range     = "single subcutaneous abdominal infusion of 3.6 MIU (2,700 mg, n = 4), 7.2 MIU (5,400 mg, n = 10), or 10.8 MIU (8,100 mg, n = 10) benzathine penicillin G (Bicillin L-A), infused over up to 30 minutes",
    regions        = "Australia (Perth, Western Australia)",
    notes          = paste(
      "Phase 1 open-label dose-escalation study (ACTRN12621000135819). Enrolment was",
      "balanced by body composition: 12 participants with ideal BMI (20-24.9 kg/m^2)",
      "and 12 with higher BMI (25-34.9 kg/m^2). 400 dried-blood-spot concentrations",
      "were included in the PK analysis; 3.0% fell between the LOD and LOQ and were",
      "used at their measured value, and 7.7% below the LOD were excluded. Plasma and",
      "DBS concentrations were equivalent (bootstrap ratio 0.990, 95% CI 0.897-1.08),",
      "so the final model was fitted to DBS data only. Three participants received",
      "less than the planned dose (two in cohort 1 received 3.1 MIU; one in cohort 3",
      "received 8.2 MIU). Sampling: 2, 6, 12, 24, 48, 72 h and 5, 7, 14, 21, 28, 42,",
      "56, 70, 84, 98, 112 days post-dose. The highest observed benzylpenicillin",
      "concentration was under 300 ng/mL. Baseline demographics: Results,",
      "'Participants'. Final model estimates and bootstrap CIs: Table 2. Model",
      "schematic: Figure 2. pcVPC: Figure 3. Dosing-interval simulations: Figure 4.",
      sep = " "
    )
  )

  ini({
    # ==========================================================================
    # Disposition
    # ==========================================================================
    # kel was FIXED, not estimated (Table 2 bootstrap column reads "Fixed", and
    # Methods 'Population pharmacokinetic analysis' states the elimination rate
    # was fixed following the group's earlier models). Table 2 reports it in
    # h^-1 for a 70 kg individual; converted to 1/day because every absorption
    # parameter in this model is expressed in days.
    lkel <- fixed(log(1.32 * 24)); label("Elimination rate constant, 70 kg reference (1/day)")  # Table 2: kel = 1.32 h^-1 . 70 kg^-1, Fixed -> 31.68 1/day
    lvc  <- log(42.8);             label("Apparent central volume of distribution V/F, 70 kg reference (L)")  # Table 2: V/F = 42.8 L . 70 kg^-1 (bootstrap 42.8, 95% CI 40-46.1)

    # Allometric exponents were applied a priori, not estimated (Methods,
    # 'Population pharmacokinetic analysis': "Allometric scaling for size was
    # employed a priori, with an exponential of 1 for volume (V) and 3/4 for
    # clearance (CL) terms").
    e_wt_vc <- fixed(1);    label("Allometric exponent of WT on V/F (unitless), applied a priori")  # Methods 'Population pharmacokinetic analysis': exponent 1 for volume
    e_wt_cl <- fixed(0.75); label("Allometric exponent of WT on clearance (unitless), applied a priori")  # Methods 'Population pharmacokinetic analysis': exponent 3/4 for clearance

    # ==========================================================================
    # Absorption
    # ==========================================================================
    # Figure 2: Bolus and a zero-order process (DUR) both enter the Depot;
    # Depot -> Transit is governed by t1/2,tr; Transit -> V/F by t1/2,abs.
    # Half-lives are converted to rate constants inside model(): k = log(2)/t1/2.
    lthalf_abs <- log(11.8);        label("Absorption half-life from the transit compartment, t1/2,abs (day)")  # Table 2: t1/2,abs = 11.8 days (bootstrap 11.8, 95% CI 10.8-12.8)
    lthalf_tr  <- log(0.702 / 24);  label("Transit-compartment half-life, t1/2,tr (day)")  # Table 2: t1/2,tr = 0.702 hours (bootstrap 0.69, 95% CI 0.519-0.873) -> 0.02925 days
    ldur       <- log(1062 / 24);   label("Duration of the zero-order input into the depot, DUR (day)")  # Table 2: DUR = 1,062 hours (bootstrap 1,055, 95% CI 832-1,250) -> 44.25 days

    # RATIO is the amount of dose assigned to the bolus path RELATIVE TO the
    # zero-order path (Table 2 footnote b: "RATIO, the relative amount of drug
    # assigned to bolus over zero-order input into the depot compartment").
    # Hence the zero-order fraction is 1/(1 + RATIO) = 1/3.54 = 28.2%, which
    # reproduces the paper's statement that "the population estimate for RATIO
    # corresponds to 28% of the dose passing through the zero-order absorption
    # pathway" (Results, 'Pharmacokinetic modeling').
    lrat_bolus <- log(2.54);        label("Dose-split ratio RATIO, bolus path over zero-order path (unitless)")  # Table 2: RATIO = 2.54 (bootstrap 2.47, 95% CI 1.66-4.65)

    # Power relationship of BMI on DUR. See covariateData$BMI$notes for the
    # reconstruction of the 25 kg/m^2 centering value used in model().
    e_bmi_dur  <- 1.13;             label("Exponent of the power BMI relationship on DUR (unitless)")  # Table 2: exponent of power BMI relationship on DUR = 1.13 (RSE 28%; bootstrap 1.11, 95% CI 0.21-2.22)

    # ==========================================================================
    # Inter-individual variability
    # ==========================================================================
    # Table 2 footnote b: "IIV and RV are presented as 100% x sqrt(variability
    # estimate)". The tabulated percentage is therefore 100 x omega itself, so
    # omega^2 = (percentage / 100)^2 directly -- NOT log(1 + CV^2). IIV was
    # modeled exponentially on every parameter (Methods, 'Population
    # pharmacokinetic analysis': "The IIV was exponentially modeled for all
    # parameters"). No IIV was reported on t1/2,tr.
    etalvc        ~ 0.0225;  # Table 2: IIV in V/F = 15% (shrinkage 9%); omega^2 = 0.15^2
    etaldur       ~ 0.0484;  # Table 2: IIV in DUR = 22% (shrinkage 14%); omega^2 = 0.22^2
    etalthalf_abs ~ 0.0144;  # Table 2: IIV in t1/2,abs = 12% (shrinkage 26%); omega^2 = 0.12^2
    etalrat_bolus ~ 0.7744;  # Table 2: IIV in RAT = 88% (shrinkage 16%); omega^2 = 0.88^2 -- the paper notes this parameter had the highest IIV (88%)

    # ==========================================================================
    # Residual variability
    # ==========================================================================
    # "Residual variability (RV) was estimated as an additive error for the
    # log-transformed data" (Methods, 'Population pharmacokinetic analysis'),
    # which is a proportional error model on the linear concentration scale.
    propSd <- 0.24; label("Proportional residual error (fraction)")  # Table 2: RV = 24% (shrinkage 10%; bootstrap 24, 95% CI 21-26)
  })

  model({
    # ---- Reference values ----------------------------------------------------
    wtRef  <- 70   # kg      -- Table 2 reports kel and V/F per 70 kg
    bmiRef <- 25   # kg/m^2  -- reconstructed population average BMI; see covariateData$BMI$notes

    # ---- Disposition ---------------------------------------------------------
    # Volume carries the estimated allometric exponent of 1 and the only
    # disposition IIV. kel is a fixed rate constant, so it scales with the
    # difference of the clearance and volume exponents (3/4 - 1 = -1/4) and
    # carries no IIV -- an individual with a larger V/F therefore has a
    # proportionally larger CL at the same kel, which is how the paper's
    # fixed-kel parameterisation behaves. The two exponents are applied on
    # separate factors rather than as (e_wt_cl - e_wt_vc) to avoid a
    # theta-minus-theta expression in the exponent.
    wtScale <- WT / wtRef
    vc  <- exp(lvc + etalvc) * wtScale^e_wt_vc
    kel <- exp(lkel) * wtScale^e_wt_cl / wtScale^e_wt_vc

    # ---- Absorption ----------------------------------------------------------
    thalf_abs <- exp(lthalf_abs + etalthalf_abs)
    thalf_tr  <- exp(lthalf_tr)
    rat_bolus <- exp(lrat_bolus + etalrat_bolus)
    # Power BMI relationship on DUR, per the Supplemental Information general
    # power form theta_i = theta_POP * (x_i / x_tilde)^theta_POWER.
    durZero   <- exp(ldur + etaldur) * (BMI / bmiRef)^e_bmi_dur

    ka  <- log(2) / thalf_abs
    ktr <- log(2) / thalf_tr

    # ---- Dose split between the two depot inputs -----------------------------
    # RATIO = bolus / zero-order, so the two fractions sum to 1 (all of the
    # dose enters the depot; the unknown absolute bioavailability is absorbed
    # into the apparent volume V/F).
    fBolus <- rat_bolus / (rat_bolus + 1)
    fZero  <- 1        / (rat_bolus + 1)

    # ---- ODE system ----------------------------------------------------------
    # Figure 2 draws a SINGLE depot fed by both a bolus and a zero-order input.
    # rxode2 applies one bioavailability per compartment, so the two inputs are
    # carried as two states, `depot` (bolus path) and `depot2` (zero-order
    # path), which drain into the transit compartment with the SAME rate
    # constant ktr. Because the depot outflow is first-order and the system is
    # linear, depot + depot2 obeys exactly the single-depot equation of Figure 2
    # -- the split is an implementation device, not a structural change.
    d/dt(depot)    <- -ktr * depot
    d/dt(depot2)   <- -ktr * depot2
    d/dt(transit1) <-  ktr * (depot + depot2) - ka * transit1
    d/dt(central)  <-  ka * transit1 - kel * central

    # ---- Dose partitioning ---------------------------------------------------
    # Dose BOTH depots with the full dose amount; f() partitions it. The
    # `depot2` record must request a modelled duration (rate = -2 in the event
    # table) so that dur(depot2) below delivers fZero * dose at a constant rate
    # over durZero days.
    f(depot)    <- fBolus
    f(depot2)   <- fZero
    dur(depot2) <- durZero

    # ---- Observation ---------------------------------------------------------
    # central in mg and vc in L give mg/L; x1000 reports ng/mL to match Table 2
    # and the 20 / 10 ng/mL PK/PD targets used throughout the paper.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
