Chan_2023_nirmatrelvir <- function() {
  description <- "Two-compartment population PK model with first-order absorption for nirmatrelvir coadministered with ritonavir 100 mg (PAXLOVID) in adults with and without COVID-19 (Chan 2023; pooled analysis of 8 phase I and phase II/III studies, N = 1237). Allometric baseline body weight on clearances (fixed 0.75) and volumes (fixed 1) referenced to 70 kg; a breakpoint power model for BSA-normalized creatinine clearance on CL (CL scales with nCLCR below the estimated 70.1 mL/min/1.73 m2 breakpoint and is independent of nCLCR at or above it); fractional carbamazepine, itraconazole and COVID-19 effects on CL; a power effect of age on central volume referenced to 45 years; and a relative-bioavailability model combining a dose power function (referenced to 300 mg) with a fractional 150-mg-tablet formulation effect. Combined additive plus proportional residual error with separate proportional magnitudes for the phase I and phase II/III data."
  reference <- paste(
    "Chan PLS, Singh RSP, Cox DS, Shi H, Damle B, Nicholas T.",
    "Dosing recommendation of nirmatrelvir/ritonavir using an integrated",
    "population pharmacokinetic analysis.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(12):1897-1910.",
    "doi:10.1002/psp4.13039",
    sep = " "
  )
  vignette <- "Chan_2023_nirmatrelvir"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    depot       = list(analyte = "nirmatrelvir", units = "mg", specimen = "administration site", verified = FALSE),
    central     = list(analyte = "nirmatrelvir", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "nirmatrelvir", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed baseline weight (BBWT in the source control stream).",
        "Standard allometric scaling referenced to 70 kg with exponents fixed",
        "at 0.75 for CL and Q and 1 for V2 and V3 to reduce collinearity among",
        "covariates (Chan 2023 Results, 'Population PK modeling').",
        "Median 79.4 kg, range 42.0-158 kg (Table 2)."
      ),
      source_name = "BBWT"
    ),
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on central volume V2, normalized to 45 years",
        "(the analysis-set median age). Median 45.0 years, range 18.0-86.0",
        "years (Table 2)."
      ),
      source_name = "AGE"
    ),
    CRCL = list(
      description = paste(
        "Baseline body surface area-normalized creatinine clearance (nCLCR),",
        "i.e. a Cockcroft-Gault-style creatinine clearance rescaled to",
        "1.73 m2 of body surface area"
      ),
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Creatinine-based estimate, BSA-normalized. Enters CL through a",
        "breakpoint power model: CL scales as (nCLCR / breakpoint)^power below",
        "the estimated breakpoint of 70.1 mL/min/1.73 m2 and is independent of",
        "nCLCR at or above it (Chan 2023 Abstract and Results). Median 119,",
        "range 15.8-318 mL/min/1.73 m2 (Table 2). Renal-function categories in",
        "the source: normal >= 90, mild 60 to < 90, moderate 30 to < 60,",
        "severe < 30 mL/min/1.73 m2."
      ),
      source_name = "NCLCR"
    ),
    DOSE_NMV_MG = list(
      description = "Nirmatrelvir dose amount associated with the current record",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The administered dose level enters relative bioavailability as a power",
        "function normalized to 300 mg, F1 = (DOSE / 300)^power, capturing",
        "dose-dependent absorption (Chan 2023 Results and Data S1).",
        "Nirmatrelvir doses contributing to the analysis span 75-750 mg",
        "(Table 1). Set this column to the nirmatrelvir amount in mg on every",
        "record, including observation records. The drug-qualified name is",
        "required rather than the bare DOSE canonical because rxode2's event",
        "translation consumes a data column literally named DOSE, so the model",
        "would fail to solve with 'the following parameter(s) are required for",
        "solving: DOSE'."
      ),
      source_name = "DOSE"
    ),
    CONMED_CBZ = list(
      description = "Concomitant carbamazepine (strong CYP3A inducer) coadministration",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no carbamazepine)",
      notes = paste(
        "Derived from the source DRUG column (DRUG = 3 is",
        "nirmatrelvir/ritonavir + carbamazepine; Data S1 header comments).",
        "Marker for CYP3A induction; increases CL by a fraction."
      ),
      source_name = "DRUG == 3"
    ),
    CONMED_ITRACONAZOLE = list(
      description = "Concomitant itraconazole (strong CYP3A / P-gp inhibitor) coadministration",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no itraconazole)",
      notes = paste(
        "Derived from the source DRUG column (DRUG = 4 is",
        "nirmatrelvir/ritonavir + itraconazole; Data S1 header comments).",
        "Marker for additional CYP3A / P-gp inhibition on top of ritonavir;",
        "decreases CL by a fraction."
      ),
      source_name = "DRUG == 4"
    ),
    DIS_COVID19 = list(
      description = "Adult with COVID-19 (phase II/III EPIC-HR patient) vs phase I participant",
      units = "(binary)",
      type = "binary",
      reference_category = paste(
        "0 (phase I participant: healthy volunteer or a renal- or",
        "hepatic-impairment cohort member)"
      ),
      notes = paste(
        "PTST column in the source control stream (Data S1); 1 for the 1087",
        "nonhospitalized symptomatic adults with COVID-19 in EPIC-HR",
        "(NCT04960202), 0 for the 150 phase I participants. Applied as a",
        "fractional multiplier on CL. The reference group is not uniformly",
        "healthy - it also contains the phase I renal- and hepatic-impairment",
        "cohorts - so this is coded as a COVID-19 patient indicator rather",
        "than as an inverted DIS_HEALTHY column."
      ),
      source_name = "PTST"
    ),
    FORM_NMV_TAB150 = list(
      description = "Nirmatrelvir 150-mg tablet formulation",
      units = "(binary)",
      type = "binary",
      reference_category = paste(
        "0 (oral suspension or 100-mg tablet; both carry the F1 = 1 anchor)"
      ),
      notes = paste(
        "Derived from the source FFORM column (0 = suspension, 1 = 100-mg",
        "tablet, 2 = 150-mg tablet; Data S1 header comments). Only FFORM = 2",
        "carries an estimated effect, so the suspension and the 100-mg tablet",
        "share the reference bioavailability. Chan 2023 Discussion notes the",
        "150-mg tablet was evaluated only in a single single-dose healthy-",
        "participant study, so the formulation effect on F1 and the COVID-19",
        "effect on CL are likely partially confounded."
      ),
      source_name = "FFORM == 2"
    ),
    STUDY_NMV_PHASE23 = list(
      description = "Observation originates from the phase II/III EPIC-HR study",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (phase I study observation)",
      notes = paste(
        "PROT == 1005 in the source control stream (Data S1 $ERROR block).",
        "Selects the phase II/III proportional residual magnitude in place of",
        "the phase I magnitude. Distinct from DIS_COVID19: this column is a",
        "study-design / sampling-scheme property of the observation record",
        "(sparse outpatient sampling), whereas DIS_COVID19 is a subject-level",
        "disease covariate on CL. The two are collinear in the source data set",
        "but enter different parts of the model."
      ),
      source_name = "PROT == 1005"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1237L,
    n_studies = 8L,
    n_observations = "5149 plasma samples (4404 evaluable; 745 BLQ, 14.5%, excluded)",
    age_range = "18.0-86.0 years",
    age_median = "45.0 years",
    weight_range = "42.0-158 kg",
    weight_median = "79.4 kg",
    sex_female_pct = 46.9,
    race_ethnicity = c(
      White = 69.9, Black = 8.5, Asian = 13.1,
      AmericanIndianAlaskaNative = 7.7, Other = 0.3, Unknown = 0.5
    ),
    bmi_median = "27.9 kg/m^2 (range 16.6-58.1)",
    disease_state = paste(
      "1087 (87.9%) nonhospitalized symptomatic adults with COVID-19 at",
      "increased risk of progression to severe illness (phase II/III EPIC-HR);",
      "150 (12.1%) phase I participants comprising healthy volunteers and",
      "renal- or hepatic-impairment cohorts"
    ),
    renal_function = paste(
      "nCLCR median 119 mL/min/1.73 m^2 (range 15.8-318); normal 84.5%,",
      "mild 11.9%, moderate 2.8%, severe 0.8%"
    ),
    hepatic_function = paste(
      "Includes a phase I moderate-hepatic-impairment cohort",
      "(Child-Pugh class B, 7-9 points; NCT05005312)"
    ),
    co_medication = paste(
      "All nirmatrelvir doses coadministered with ritonavir 100 mg as a PK",
      "enhancer. Dedicated DDI studies with carbamazepine, itraconazole,",
      "dabigatran and midazolam contributed data"
    ),
    dose_range = paste(
      "Nirmatrelvir 75-750 mg with ritonavir 100 mg, single and q12h dosing;",
      "oral suspension, 100-mg tablet and 150-mg tablet formulations"
    ),
    regions = "Multinational (EPIC-HR global phase II/III programme; phase I studies in the US)",
    notes = paste(
      "Demographics from Chan 2023 Table 2 (All column); study designs and",
      "dosing regimens from Table 1. Concentrations below the 10 ng/mL LLOQ",
      "were excluded rather than imputed."
    )
  )

  ini({
    # ---- Structural parameters (Table 3, 'Estimate' column) ----
    # Reference subject: 70 kg, 45 years, nCLCR >= 70.1 mL/min/1.73 m^2,
    # no carbamazepine or itraconazole, no COVID-19, 300 mg dose,
    # suspension or 100-mg tablet.
    lcl <- log(9.09); label("Clearance CL (L/h)")                                      # Table 3: CL = 9.09 L/h (%RSE 3.64)
    lvc <- log(56.9); label("Central volume V2 (L)")                                   # Table 3: V2 = 56.9 L (%RSE 4.32)
    lq <- log(1.28); label("Intercompartmental clearance Q (L/h)")                     # Table 3: Q = 1.28 L/h (%RSE 14.2)
    lvp <- log(12.8); label("Peripheral volume V3 (L)")                                # Table 3: V3 = 12.8 L (%RSE 11.1)
    lka <- log(0.873); label("First-order absorption rate constant ka (1/h)")          # Table 3: ka = 0.873 1/h (%RSE 8.94)
    lfdepot <- fixed(log(1)); label("Relative bioavailability F1 at 300 mg (fraction)") # Table 3 Note: "Fixed parameters: F1 = 1, normalized to 300 mg"; Data S1 $THETA 10 = 1 FIX

    # ---- Allometric scaling on baseline body weight, reference 70 kg ----
    # Both exponents fixed, not estimated (Chan 2023 Results: "standard
    # allometric scaling of BBWT with fixed exponents of 0.75 for CLs and 1
    # for volumes to reduce collinearity among covariates").
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent of body weight shared by CL and Q (unitless)")   # Results, Population PK modeling; Data S1 CLWT = (BBWT/70)**0.75
    e_wt_vc_vp <- fixed(1); label("Allometric exponent of body weight shared by V2 and V3 (unitless)")    # Results, Population PK modeling; Data S1 VWT = (BBWT/70)

    # ---- Covariate effects (Table 3) ----
    e_crcl_bp_cl <- 70.1; label("Breakpoint of nCLCR on CL (mL/min/1.73 m^2)")         # Table 3: nCLCR breakpoint = 70.1 mL/min/1.73 m^2 (%RSE 0.03)
    e_crcl_cl <- 1.05; label("Power of nCLCR on CL below the breakpoint (unitless)")   # Table 3: nCLCR power < nCLCR breakpoint = 1.05 (%RSE 8.44)
    e_dose_fdepot <- -0.409; label("Power of dose on relative bioavailability F1, dose normalized to 300 mg (unitless)")  # Table 3: F1 power = -0.409 (%RSE 8.7)
    e_cbz_cl <- 0.740; label("Fractional change in CL with concomitant carbamazepine")     # Table 3: Effect of carbamazepine on CL = 0.740 (%RSE 27.1)
    e_itraconazole_cl <- -0.308; label("Fractional change in CL with concomitant itraconazole")  # Table 3: Effect of itraconazole on CL = -0.308 (%RSE 7.19)
    e_form_nmv_tab150_fdepot <- -0.379; label("Fractional change in F1 for the 150-mg tablet")   # Table 3: Effect of 150-mg tablet on F1 = -0.379 (%RSE 10.1)
    e_age_vc <- -0.425; label("Power of age on V2, age normalized to 45 years (unitless)")       # Table 3: Power of age effect on V2 = -0.425 (%RSE 17.6)
    e_dis_covid19_cl <- -0.341; label("Fractional change in CL in adults with COVID-19")         # Table 3: Effect of COVID-19 on CL = -0.341 (%RSE 10.7)

    # ---- Interindividual variability ----
    # Table 3 footnote: "CV, coefficient of variation (computed as
    # sqrt(omega^2) x 100%)", so omega^2 = (CV/100)^2 directly - NOT the
    # log-normal log(CV^2 + 1) transformation.
    etalcl ~ 0.128881    # Table 3: IIV CL = 35.9 %CV -> 0.359^2 (%RSE 48.8; eta-shrinkage 55.9%)
    etalvc ~ 0.074529    # Table 3: IIV V2 = 27.3 %CV -> 0.273^2 (%RSE 17.6; eta-shrinkage 68.8%)
    etalka ~ 0.368449    # Table 3: IIV ka = 60.7 %CV -> 0.607^2 (%RSE 20.9; eta-shrinkage 63.1%)
    etalvp ~ 0.344569    # Table 3: IIV V3 = 58.7 %CV -> 0.587^2 (%RSE 26.6; eta-shrinkage 79.2%)
    # No IIV on Q: Data S1 sets Q = TVQ*CLWT with no ETA. No OMEGA block was
    # retained because all IIV correlations were < 0.6 (Chan 2023 Results).

    # ---- Residual error ----
    # Data S1 $ERROR fits log-transformed concentrations with
    # W = sqrt(THETA(6)^2 + THETA(7)^2 / IPRED^2) and Y = LPRED + EPS(1)*W,
    # EPS fixed to 1. On the linear scale the residual SD is
    # IPRED * W = sqrt(THETA(6)^2 * IPRED^2 + THETA(7)^2), i.e. exactly a
    # combined proportional + additive error, so it maps to
    # add(addSd) + prop(propSd) in nlmixr2's linear space.
    propSdPhase1 <- 0.324; label("Proportional residual SD for phase I observations (fraction)")        # Table 3: Proportional error phase I = 32.4% (%RSE 5.69)
    propSdPhase23 <- 1.39; label("Proportional residual SD for phase II/III observations (fraction)")   # Table 3: Proportional error phase II/III = 139% (%RSE 3.81)
    addSd <- fixed(10); label("Additive residual SD (ng/mL)")                                           # Table 3 Note: "additive error = 10 ng/ml"; Data S1 $THETA 7 = 10 FIX (fixed to the LLOQ)
  })

  model({
    # 1. Derived covariate terms
    # Standard allometry on baseline body weight, reference 70 kg.
    wtCl <- (WT / 70)^e_wt_cl_q
    wtV <- (WT / 70)^e_wt_vc_vp

    # Fractional (1 + theta) multipliers on CL for the two CYP3A perpetrator
    # comedications and for COVID-19 disease status
    # (Data S1: CLDRUG = 1 + THETA(13) or 1 + THETA(14); CLPTST = 1 + THETA(17)).
    clCov <-
      (1 + e_cbz_cl * CONMED_CBZ) *
      (1 + e_itraconazole_cl * CONMED_ITRACONAZOLE) *
      (1 + e_dis_covid19_cl * DIS_COVID19)

    # Renal-function breakpoint model. CL scales as a power of nCLCR below the
    # breakpoint and is flat at or above it, i.e. the effective nCLCR is capped
    # at the breakpoint (Data S1: CLR = TVCL; IF (NCLCR.LT.THETA(8))
    # CLR = TVCL*(NCLCR/THETA(8))**THETA(9)).
    crclCl <- (min(CRCL, e_crcl_bp_cl) / e_crcl_bp_cl)^e_crcl_cl

    # 2. Individual parameters
    cl <- exp(lcl + etalcl) * clCov * crclCl * wtCl
    vc <- exp(lvc + etalvc) * (AGE / 45)^e_age_vc * wtV
    q <- exp(lq) * wtCl
    vp <- exp(lvp + etalvp) * wtV
    ka <- exp(lka + etalka)

    # Relative bioavailability: dose power function normalized to 300 mg times
    # the fractional 150-mg-tablet effect (Data S1: TVF1 = F1COV*THETA(10);
    # F1 = TVF1*(DOSE/300)**THETA(11)).
    fdepot <-
      exp(lfdepot) *
      (1 + e_form_nmv_tab150_fdepot * FORM_NMV_TAB150) *
      (DOSE_NMV_MG / 300)^e_dose_fdepot

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system (NONMEM ADVAN4 TRANS4)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 5. Bioavailability
    f(depot) <- fdepot

    # 6. Observation and error
    # Amounts are in mg and volumes in L, so central/vc is mg/L; the factor of
    # 1000 converts to ng/mL (Data S1 scaling S2 = V2/1000).
    Cc <- central / vc * 1000

    # Study-dependent proportional residual magnitude (Data S1 $ERROR:
    # IF (PROT.EQ.1005) W uses THETA(12) instead of THETA(6)).
    propSd <- propSdPhase1 * (1 - STUDY_NMV_PHASE23) + propSdPhase23 * STUDY_NMV_PHASE23
    Cc ~ add(addSd) + prop(propSd)
  })
}
