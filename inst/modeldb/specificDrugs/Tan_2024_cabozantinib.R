Tan_2024_cabozantinib <- function() {
  description <- "Two-compartment population PK model for oral cabozantinib in adults with metastatic renal cell carcinoma (Tan 2024, n=27 real-world therapeutic-drug-monitoring patients, Leiden University Medical Center). The structure is the FDA cabozantinib registration popPK model reproduced from the registration file: parallel dual lagged first-order absorption, where a fraction F1 of the dose enters a fast depot (rate ka1, lag ALAG1 = 0.459 h) and the remaining (1 - F1) enters a slow depot (rate ka2, lag ALAG2 = 16.8 h), feeding a two-compartment (central + peripheral1) disposition with linear elimination. Absorption rate ka1 scales with the administered dose via a power function (DOSE/60 mg)^-0.5. Covariates on CL/F carried over from the registration model are female sex (21% lower) and Asian race (27% lower). Tan 2024 re-estimated only apparent clearance CL/F (3.11 L/h vs 2.23 L/h in the registration file) and the proportional residual error, holding every other parameter fixed, and reduced the IIV variance on F1 from 0.385 to 0.05 to stop the absorption fraction leaking above 1. A fixed 50% increase in bioavailability for a high-fat meal (FED_HIGHFAT) reproduces the paper's drug-expense-saving simulations."
  reference <- paste(
    "Tan Z, Voller S, Yin A, Rieborn A, Gelderblom AJ, van der Hulle T,",
    "Knibbe CAJ, Moes DJAR.",
    "Population pharmacokinetics of cabozantinib in metastatic renal cell",
    "carcinoma patients: towards drug expenses saving regimens.",
    "Clin Pharmacokinet. 2024;63(7):1015-1025.",
    "doi:10.1007/s40262-024-01379-y.",
    "Structural model and all fixed parameters reproduced from the FDA",
    "cabozantinib registration popPK model (Tan 2024 reference 13);",
    "the NONMEM control stream is reproduced in Tan 2024 Supplementary",
    "Material.",
    sep = " "
  )
  vignette <- "Tan_2024_cabozantinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Both absorption depots must receive a dose record carrying the SAME amount;
  # f(depot) and f(depot2) then split it into the F1 and 1 - F1 fractions, which
  # is how the source control stream's F1 / F2 = 1 - F1 pair works. Declared
  # explicitly because the automatic detection assumes the usual depot/central
  # pair and would otherwise advertise `central` as a dosing compartment, which
  # for this model would silently bypass both absorption arms.
  dosing <- c("depot", "depot2")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against the Supplementary Material control
  # stream $MODEL block: COMP(ABSORP1), COMP(ABSORP2), COMP(CENTRAL),
  # COMP(PERIPH). The control stream's fifth state, COMP(AUC), is a
  # bookkeeping integrator (DADT(5) = A(3)/V2) and is NOT reproduced here --
  # rxode2 exposes the concentration directly, so the vignette integrates the
  # AUC outside the model.
  compartmentData <- list(
    depot       = list(analyte = "cabozantinib", units = "mg", specimen = "administration site", verified = TRUE),
    depot2      = list(analyte = "cabozantinib", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "cabozantinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cabozantinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    DOSE = list(
      description        = "Administered cabozantinib dose level on the current dose record",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-dose-record. Power covariate on the fast-depot absorption rate constant: ka1(DOSE) = ka1_ref * (DOSE / 60 mg)^-0.5. The reference dose of 60 mg is explicit in the Supplementary Material control stream ($PK: KA1 = TVKA1*EXP(ETA(3))*EXP(ALPHA*LOG(DOS/60))), not inferred. The control stream derives DOS from a DOSEFLAG column restricted to the three marketed tablet strengths (20, 40, 60 mg); this model takes the dose in mg directly. Use case (a) of the DOSE canonical.",
      source_name        = "DOSEFLAG"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male, the typical-value reference)",
      notes              = "Time-fixed. Multiplicative effect on CL/F: females have 21% lower CL/F (multiplier 0.79). Tan 2024 Table 2 reports the transformed multiplier 0.79; the Supplementary Material control stream stores the fractional form, THETA(12) = 0.21 with GEND = 1 - THETA(12) when SEX == 0. NOTE THE ORIENTATION: in the Tan 2024 dataset SEX == 0 denotes FEMALE (the control-stream branch that applies the 0.79 multiplier), which is the inverse of the more common 1 = male / 0 = female coding. This model uses the canonical SEXF orientation (1 = female), so SEXF = 1 - SEX relative to the paper's own column. Fixed to the FDA registration value; not re-estimated by Tan 2024.",
      source_name        = "SEX"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian; the typical-value reference)",
      notes              = "Time-fixed. Multiplicative effect on CL/F: Asian subjects have 27% lower CL/F (multiplier 0.73; Tan 2024 Table 2 and section 2.4.1). Fixed to the FDA registration value. Tan 2024's own analysis could not exercise this covariate: section 2.4.2 states that because ethnicity was not recorded in the TDM dataset, all patients were assumed Caucasian, and the control stream correspondingly hardcodes RACE = 1 (i.e. RACE_ASIAN = 0 for every subject). The term is retained here because it is part of the model Tan 2024 Table 2 publishes; set RACE_ASIAN = 0 to reproduce every simulation in the paper.",
      source_name        = "RACE"
    ),
    FED_HIGHFAT = list(
      description        = "High-fat-meal indicator for the current dose record",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted; the reference condition and the state in which the model was fitted)",
      notes              = "Per-dose-record. Multiplies overall bioavailability by 1.5 (a 50% increase). This is a stated simulation ASSUMPTION, not an estimated parameter: Tan 2024 section 2.6 says 'it was assumed that the bioavailability of the drug will increase by 50% when taken with high-fat meals', citing a phase I food-effect study (Tan 2024 reference 18) in which a high-fat high-calorie meal raised cabozantinib Cmax by 41% and AUC by 57%. The term does not appear in either Supplementary Material control stream, both of which simulate the fasted state only. Set FED_HIGHFAT = 0 to recover the fitted model exactly. The food-effect study used the capsule formulation whereas the TDM cohort took tablets; Tan 2024 section 4 flags this as a limitation (the two formulations are similar but not bioequivalent).",
      source_name        = NULL
    )
  )

  # Screened in Tan 2024 but not carried by the model. Section 2.4.3: "due to
  # the limited sample size of the real-world dataset, new covariates
  # exploration was skipped and covariates in the FDA cabozantinib POPPK model
  # were maintained." These demographics and laboratory values are collected
  # and tabulated (Tan 2024 Table 1, Supplementary Table S1) but no effect is
  # estimated for any of them, so there is nothing to encode in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at start of cabozantinib treatment",
      units       = "years",
      type        = "continuous",
      notes       = "Collected and tabulated (Tan 2024 Table 1, mean 65 y, range 39-85 y) but no covariate effect was estimated; new-covariate exploration was skipped for sample-size reasons (section 2.4.3)."
    ),
    WT = list(
      description = "Body weight at start of cabozantinib treatment",
      units       = "kg",
      type        = "continuous",
      notes       = "Collected and tabulated (Tan 2024 Table 1, mean 78 kg, range 49-105 kg) but no covariate effect was estimated. The FDA registration model this one reproduces carries no allometric term either."
    ),
    HT = list(
      description = "Height at start of cabozantinib treatment",
      units       = "cm",
      type        = "continuous",
      notes       = "Present in the control stream $INPUT list and tabulated (Tan 2024 Table 1, mean 178 cm, range 160-196 cm) but unused in $PK."
    ),
    CRCL = list(
      description = "Creatinine clearance by the CKD-EPI equation, normalised to 1.73 m^2 body surface area",
      units       = "mL/min/1.73m^2",
      type        = "continuous",
      notes       = "Collected and tabulated (Tan 2024 Table 1, mean 70 mL/min, range 31-121 mL/min) and present in the control stream $INPUT list, but unused in $PK. Section 2.2 specifies the CKD-EPI equation."
    ),
    ALT = list(
      description = "Alanine aminotransferase at start of cabozantinib treatment",
      units       = "U/L",
      type        = "continuous",
      notes       = "Collected and tabulated (Tan 2024 Table 1, mean 46 U/L, range 14-202 U/L) and present in the control stream $INPUT list as ALAT, but unused in $PK."
    ),
    AST = list(
      description = "Aspartate aminotransferase at start of cabozantinib treatment",
      units       = "U/L",
      type        = "continuous",
      notes       = "Collected and tabulated (Tan 2024 Table 1, mean 51 U/L, range 14-449 U/L) and present in the control stream $INPUT list as ASAT, but unused in $PK."
    ),
    BILI = list(
      description = "Total bilirubin at start of cabozantinib treatment",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Collected and tabulated (Tan 2024 Table 1, mean 8 umol/L, range 3-21 umol/L) and present in the control stream $INPUT list as Bil, but unused in $PK."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 27L,
    n_observations = 75L,
    n_studies      = 1L,
    age_range      = "39-85 years",
    age_median     = "68 years (mean 65 years)",
    weight_range   = "49-105 kg",
    weight_median  = "79 kg (mean 78 kg)",
    height_range   = "160-196 cm",
    bmi_range      = "18-32 kg/m^2 (mean 25 kg/m^2)",
    sex_female_pct = 29.7,
    race_ethnicity = "Not recorded. Section 2.4.2: because ethnicity data were absent from the TDM dataset, all patients were assumed Caucasian on the grounds that roughly 90% of the Leiden University Medical Center population is Caucasian.",
    disease_state  = "Metastatic renal cell carcinoma. IMDC prognosis group favorable 7.4%, intermediate 63.0%, poor 18.5%, unknown 11.1%. WHO performance score 0 in 33.3%, 1 in 40.7%, >1 in 26.0%. Most patients had prior systemic therapy (pazopanib 42.9%, nivolumab +/- ipilimumab 34.3%, sunitinib 11.4%, everolimus 2.9%); 8.5% were treatment-naive.",
    renal_function = "Mean baseline creatinine clearance 70 mL/min/1.73m^2 (range 31-121) by CKD-EPI",
    dose_range     = "20-60 mg once-daily oral cabozantinib tablets (median 40 mg). Starting dose 20 mg in 15%, 40 mg in 37%, 60 mg in 48%; last recorded dose 20 mg in 19%, 40 mg in 60%, 60 mg in 21%. Median treatment duration 75 days (range 11-552).",
    regions        = "Single center, the Netherlands (Leiden University Medical Center)",
    notes          = "Retrospective therapeutic-drug-monitoring cohort treated between August 2018 and December 2021, one routine TDM observation minimum per patient. Baseline demographics from Tan 2024 Table 1; a side-by-side comparison against the 318 mRCC patients of the FDA registration phase III study is in Supplementary Table S1. Median 2 observations per patient (range 1-10); 36 of 75 observations (48%) were troughs. Median time after last dose 25.30 h (range 1.25-267.15 h). Steady state was NOT assumed; all administered doses were carried in the dataset via the NONMEM ADDL and II columns. Observed concentrations median 603 ng/mL (range 135-1471); trough concentrations median 632 ng/mL (range 308-1134). Assay: UPLC-MS/MS validated over 10-4000 ng/mL. IMPORTANT -- this 27-patient cohort was used only to EVALUATE and partially re-estimate the FDA registration model, which was itself developed on 63 healthy participants (phase I) plus 325 mRCC patients (phase III); every parameter other than CL/F and the proportional residual error is inherited from that much larger analysis."
  )

  ini({
    # ---------------------------------------------------------------------
    # Values are the "Estimates of the final POPPK model" column of Tan 2024
    # Table 2, cross-checked line by line against the $THETA / $OMEGA blocks
    # of the second control stream in the Supplementary Material ("NONMEM
    # CONTROL STREAM SIMULATION WITH THE FINAL CABOZANTINIB PK MODEL").
    #
    # Only CL/F and the proportional residual error were estimated by Tan
    # 2024; every other value is starred in Table 2 as "Value fixed to the
    # FDA registration file" and is wrapped in fixed() here. The one
    # exception to that rule is the IIV variance on F1, which Tan 2024 fixed
    # to a value of its own choosing rather than the registration value --
    # see the note on etalffo below.
    #
    # All disposition parameters are APPARENT (oral) quantities: CL/F, Vc/F,
    # Vp/F, Q/F. There is no intravenous reference arm, so absolute
    # bioavailability is not identifiable and lfdepot is a fixed anchor at 1.
    # ---------------------------------------------------------------------

    # ---- Estimated by Tan 2024 ----
    lcl <- log(3.11)          ; label("Apparent oral clearance CL/F (L/h)")                                   # Table 2 final model CL/F = 3.11 L/h (RSE 5%); bootstrap median 3.11 (95% CI 2.73-3.58). Registration-file value was 2.23 L/h.

    # ---- Fixed to the FDA registration file (Table 2 asterisked rows) ----
    lvc <- fixed(log(81.5))   ; label("Apparent central volume of distribution Vc/F (L)")                     # Table 2 Vc/F = 81.5 L; control stream $THETA (81.5) FIX ; V2 [L]
    lvp <- fixed(log(213))    ; label("Apparent peripheral volume of distribution Vp/F (L)")                  # Table 2 Vp/F = 213 L; control stream $THETA (213) FIX ; V3 [L]
    lq  <- fixed(log(14.2))   ; label("Apparent inter-compartmental clearance Q/F (L/h)")                     # Table 2 Q/F = 14.2 L/h; control stream $THETA (14.2) FIX ; Q [L/h]

    # Dual absorption. Depot 1 is the fast arm, depot 2 the slow arm; both are
    # first-order with their own lag time. Tan 2024 Fig. 1 is the schematic.
    lka    <- fixed(log(0.568)); label("Fast-depot first-order absorption rate constant ka1 at the 60 mg reference dose (1/h)") # Table 2 Ka1 = 0.568 /h; control stream $THETA (0.568) FIX ; KA1 [h^-1]
    lka2   <- fixed(log(0.102)); label("Slow-depot first-order absorption rate constant ka2 (1/h)")                            # Table 2 Ka2 = 0.102 /h; control stream $THETA (0.102) FIX ; KA2
    ltlag  <- fixed(log(0.459)); label("Fast-depot absorption lag time ALAG1 (h)")                                             # Table 2 ALAG1 = 0.459 h; control stream $THETA (0.459) FIX ; ALAG1 [h]
    ltlag2 <- fixed(log(16.8)) ; label("Slow-depot absorption lag time ALAG2 (h)")                                             # Table 2 ALAG2 = 16.8 h; control stream $THETA (16.8) FIX ; ALAG2

    # F1 = fraction of the bioavailable dose routed to the fast depot; the
    # remaining (1 - F1) goes to the slow depot, so the two always sum to the
    # whole dose. The control stream parameterises this MULTIPLICATIVELY with
    # exponential IIV -- F1 = TVF1 * EXP(ETA(4)) -- rather than on the logit
    # scale, so the canonical here is lffo (not logitffo). That choice is
    # load-bearing: an exponential eta genuinely admits F1 > 1 in the upper
    # tail, which is exactly the numerical failure Tan 2024 section 3.2
    # reports and works around. Re-encoding on the logit scale would silently
    # change the model. See the etalffo note below.
    lffo <- fixed(log(0.675)) ; label("Fraction of the bioavailable dose routed to the fast depot F1 (unitless)")  # Table 2 F1 = 0.675; control stream $THETA (0.675) FIX ; F1

    # Bioavailability anchor. Neither control stream carries an F term: the
    # dual-absorption fractions sum to 1 and CL/F, Vc/F, Vp/F, Q/F absorb the
    # unknown absolute bioavailability. Fixed at 1 so that the FED_HIGHFAT
    # multiplier below has a defined reference.
    lfdepot <- fixed(log(1))  ; label("Overall oral bioavailability under fasted conditions (reference anchor)")   # Tan 2024 Supplementary Material control streams carry no F term; apparent parameters absorb F

    # ---- Dose-dependent absorption ----
    # Control stream $PK: KA1 = TVKA1*EXP(ETA(3))*EXP(ALPHA*LOG(DOS/60)),
    # i.e. ka1 = ka1_ref * (DOSE / 60)^ALPHA. The 60 mg reference is explicit
    # in the control stream. ALPHA < 0 means a higher dose absorbs more slowly.
    e_dose_ka <- fixed(-0.5)  ; label("Dose power exponent on the fast-depot absorption rate ka1 (unitless)")      # Table 2 "Dose exponent on Ka1" = -0.5; control stream $THETA (-0.5) FIX ; DOSE exponent on Ka

    # ---- Categorical covariate effects on CL/F ----
    # Encoded as the fractional change from the reference category:
    # CL/F = CL/F_ref * (1 + e * IND), so IND = 1 reproduces the multiplier
    # printed in Table 2. Table 2 footnote #: "Transformed estimates
    # correspond to multiplicative change from the typical PK parameter."
    e_sexf_cl       <- fixed(-0.21) ; label("Female (vs male) fractional change on CL/F (unitless)")              # Table 2 "Female on CL/F" = 0.79 multiplier, i.e. 21% lower; control stream $THETA (0.21) FIX ; Gender, applied as GEND = 1 - THETA(12)
    e_race_asian_cl <- fixed(-0.27) ; label("Asian (vs non-Asian) fractional change on CL/F (unitless)")          # Table 2 "Asian on CL/F" = 0.73 multiplier, i.e. 27% lower; section 2.4.1 "Asian race (27% lower CL/F)"

    # ---- Food effect (simulation assumption, not an estimate) ----
    # Section 2.6: "it was assumed that the bioavailability of the drug will
    # increase by 50% when taken with high-fat meals." Absent from both
    # control streams, which simulate fasted dosing only.
    e_fed_highfat_f <- fixed(0.5)   ; label("High-fat meal fractional change on overall bioavailability (unitless)") # Tan 2024 section 2.6 assumed +50% bioavailability with a high-fat meal, based on the phase I food-effect study (reference 18: Cmax +41%, AUC +57%)

    # ---- Inter-individual variability ----
    # All variances are on the log scale (exponential etas throughout the
    # control stream) and all are FIXed. The CL/F-Vc/F pair is a NONMEM
    # $OMEGA BLOCK(2); the off-diagonal 0.44 is a COVARIANCE, not a
    # correlation (Table 2 legend: "Omega2_CL:Vc covariance of CL/F and
    # Vc/F"). It satisfies Cauchy-Schwarz -- sqrt(0.213 * 1.06) = 0.475 >
    # 0.44 -- giving an implied correlation of 0.93, so the block is
    # positive-definite and is reproduced as published.
    etalcl + etalvc ~ fixed(c(0.213,
                              0.44, 1.06))   # Table 2 Omega2_CL = 0.213, Omega2_CL:Vc = 0.44, Omega2_Vc = 1.06; control stream $OMEGA BLOCK(2) FIX
    etalka  ~ fixed(0.437)                   # Table 2 Omega2_Ka = 0.437; control stream $OMEGA 0.437 FIX ; KA1. Applies to ka1 only -- ka2 carries no eta in $PK.

    # Tan 2024's ONE deviation from the registration file's random effects.
    # Section 3.2: the registration value 0.385 "resulted in terminated runs
    # due to a negative F1 fraction, and it was impossible to continue using
    # this value. Therefore, we decrease the variance and set it to a 20%
    # variation of F1". The mechanism is the exponential eta on F1 described
    # above: F1 = 0.675 * exp(eta) exceeds 1 whenever eta > log(1/0.675) =
    # 0.393, which makes the slow-depot fraction (1 - F1) negative. At the
    # registration variance 0.385 (SD 0.620) that is 26% of subjects; at 0.05
    # (SD 0.224) it is still about 4%. The tail is not eliminated, only
    # shrunk, and it is a genuine property of the published model -- see the
    # vignette Errata for how it shows up in simulation.
    etalffo ~ fixed(0.05)                    # Table 2 final-model Omega2_F1 = 0.05 (registration file value 0.385); control stream $OMEGA 0.05 FIX ; F1

    # ---- Residual error ----
    # Control stream $ERROR: W = SQRT(THETA(10)**2 + (THETA(11)*IPRED)**2),
    # Y = IPRED + W*EPS(1), $SIGMA 1 FIX. That is a combined additive +
    # proportional error whose VARIANCES add, which is nlmixr2's default
    # combined2 parameterisation for `prop() + add()`, so the two thetas map
    # across directly as standard deviations.
    propSd <- 0.335           ; label("Proportional residual error (fraction)")                                   # Table 2 final-model "Residual error" = 0.335 (RSE 10%), the estimated value; control stream $THETA (0.335) FIX ; Prop.RE (sd). Registration-file value was 0.254.
    addSd  <- fixed(0.001)    ; label("Additive residual error (ng/mL)")                                          # Control stream $THETA (0.001) FIX ; Add.RE (sd). Not listed in Table 2. At 0.001 ng/mL against concentrations of several hundred ng/mL this is a numerical stabiliser, not a real assay floor (the assay was validated from 10 ng/mL).
  })

  model({
    # ---- Reference values ----
    ref_dose <- 60   # mg; explicit in the control stream, EXP(ALPHA*LOG(DOS/60))

    # ---- Covariate multipliers ----
    # Reference subject: male, non-Asian, fasted. Each reference category
    # gives a multiplier of exactly 1.
    cl_cov <- (1 + e_sexf_cl * SEXF) * (1 + e_race_asian_cl * RACE_ASIAN)
    f_food <- 1 + e_fed_highfat_f * FED_HIGHFAT

    # ---- Individual PK parameters ----
    cl <- exp(lcl + etalcl) * cl_cov
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    q  <- exp(lq)

    ka  <- exp(lka + etalka) * (DOSE / ref_dose)^e_dose_ka
    ka2 <- exp(lka2)

    tlag  <- exp(ltlag)
    tlag2 <- exp(ltlag2)

    # Exponential IIV on the fast-depot fraction, exactly as the control
    # stream writes it (F1 = TVF1 * EXP(ETA(4))). This is deliberately NOT
    # bounded to (0, 1) -- see the lffo and etalffo notes in ini().
    ffo    <- exp(lffo + etalffo)
    fdepot <- exp(lfdepot) * f_food

    # ---- Micro-constants ----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---- ODE system ----
    # Reproduces $DES of the Supplementary Material control stream, with
    # A(1) = depot, A(2) = depot2, A(3) = central, A(4) = peripheral1.
    # A dose event must be recorded on BOTH depot and depot2 with the same
    # amount; the f() multipliers below split it into the F1 and (1 - F1)
    # fractions, matching NONMEM's F1 / F2 = 1 - F1 bioavailability pair.
    d/dt(depot)       <- -ka * depot
    d/dt(depot2)      <- -ka2 * depot2
    d/dt(central)     <-  ka * depot + ka2 * depot2 -
                          kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot)     <- fdepot * ffo
    f(depot2)    <- fdepot * (1 - ffo)
    alag(depot)  <- tlag
    alag(depot2) <- tlag2

    # ---- Plasma concentration ----
    # central is in mg and vc in L, giving mg/L; multiply by 1000 for ng/mL.
    # Matches the control stream's IPRED = A(3)*1000/V2 and its scaling
    # S3 = V2/1000.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd) + add(addSd)
  })
}
