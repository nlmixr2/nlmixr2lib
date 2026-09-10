Parasrampuria_2025_temsavir <- function() {
  description <- paste0(
    "Two-compartment population PK model for temsavir (TMR, GSK2616713, ",
    "formerly BMS-626529), the active moiety of the HIV-1 attachment ",
    "inhibitor prodrug fostemsavir (FTR, GSK3684934, formerly BMS-663068, ",
    "marketed as RUKOBIA extended-release tablets), pooled across seven ",
    "studies (four phase 1 studies in healthy adults plus the phase 2a ",
    "AI438006, phase 2b AI438011 and phase 3 BRIGHTE / AI438047 studies in ",
    "treatment-experienced adults with HIV-1) (Parasrampuria 2025, ",
    "n = 764, 10,236 quantifiable concentrations). Fostemsavir is a ",
    "methyl-phosphate prodrug hydrolysed by gastrointestinal alkaline ",
    "phosphatase to temsavir before absorption, so only temsavir is ",
    "modelled. Absorption is SEQUENTIAL zero-order then first-order: the ",
    "dose enters the depot as a constant-rate input of estimated duration ",
    "DUR = 3.84 h and the depot then empties into the central compartment ",
    "first-order at Ka = 2.33 1/h (the paper's 'dual zero- and first-order ",
    "absorption', modelled as 'an estimated constant infusion of duration ",
    "DUR in absorption depot compartment'); there is no parallel-pathway ",
    "fraction parameter. Retained covariates are allometrically scaled ",
    "body weight on CL/F, Q/F (exponent 0.75 fixed) and V2/F, V3/F ",
    "(exponent 1 fixed) referenced to the cohort-median 72 kg, plus ",
    "multiplicative effects of concomitant moderate CYP3A inducers ",
    "(1.41x CL/F, primarily etravirine) and strong CYP3A inhibitors ",
    "(0.721x CL/F, primarily ritonavir and cobicistat). HIV status, age, ",
    "sex, race, formulation and baseline clinical laboratory parameters ",
    "were screened and rejected. Inter-individual variability is ",
    "exponential on CL/F, V2/F and Ka with a CL/F-V2/F correlation of ",
    "0.609; residual error is additive on the natural-log scale ",
    "(log-normal) with SD 0.613, itself carrying exponential ",
    "inter-individual variability. Three companion exposure-response ",
    "models in the Parasrampuria_2025_temsavir_* family consume this ",
    "model's steady-state Ctau."
  )
  reference <- paste(
    "Parasrampuria R, Thakkar N, Moore K, Ackerman P, Magee M.",
    "Population pharmacokinetics and exposure-response relationship for",
    "temsavir following fostemsavir administration in",
    "treatment-experienced HIV patients.",
    "Pharmacol Res Perspect. 2025;13(3):e70023.",
    "doi:10.1002/prp2.70023.",
    "Studies AI438006 (NCT01009814), AI438011 (NCT01384734) and",
    "AI438047 / BRIGHTE (NCT02362503).",
    sep = " "
  )
  vignette <- "Parasrampuria_2025_temsavir"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometrically scaled against the population PK cohort MEDIAN of",
        "72 kg (Table S2: median 72 kg, range 38-151 kg, N = 764). Both",
        "exponents were FIXED, not estimated: 0.75 on CL/F and Q/F, 1 on",
        "V2/F and V3/F (Table 2, 'Effect of WT on ... 0.75 Fixed (NA)' and",
        "'1 Fixed (NA)', no RSE and no bootstrap CI). BASELINE weight, held",
        "constant per subject -- the analysis pools single-dose phase 1",
        "studies with multi-year phase 3 follow-up but Table 2 gives one",
        "weight term with no time-varying qualifier. Parasrampuria 2025",
        "Results quantify the resulting exposure spread: steady-state Ctau",
        "over the 40-150 kg baseline weight range is 1.4- to 0.7-fold that",
        "of a 72-kg subject, and Table 5 gives the simulated medians",
        "(599 ng/mL at 40 kg, 296 ng/mL at 150 kg, 433 ng/mL for the",
        "unmodified Phase 3 cohort). The authors conclude no dose",
        "adjustment is warranted across this range because the Day 8",
        "virologic response barely moves."
      ),
      source_name        = "WT"
    ),
    CONMED_CYP3A4_IND_MOD = list(
      description        = "Concomitant moderate CYP3A inducer coadministration indicator; 1 = on a moderate CYP3A inducer, 0 = not.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant moderate CYP3A inducer)",
      notes              = paste(
        "The MODERATE stratum specifically. Parasrampuria 2025 Discussion:",
        "'Simulation estimated effect of moderate CYP3A inducers, primarily",
        "etravirine (ETR)' and 'strong CYP3A inducers were not included in",
        "the analysis, as rifampin was found to decrease AUC by 82% in a",
        "separate DDI study in healthy subjects; thus, they were prohibited",
        "in patient studies.' The indicator therefore never carries a",
        "strong inducer, which is why CONMED_CYP3A4_IND_MOD is the correct",
        "canonical rather than the pooled CONMED_CYP3A4_IND. Time-varying",
        "per record in the source dataset (the paper's IND superscript sits",
        "on a per-observation covariate). 37 of 764 population PK subjects",
        "(5%) were exposed (Table S2). Multiplicative on CL/F:",
        "1.41x, i.e. a 41% clearance increase, so exposure falls.",
        "The model-based Ctau reduction of 53% (Table 5, 205 vs 433 ng/mL)",
        "matched the 52% Ctau decrease observed in a dedicated etravirine",
        "DDI study in healthy subjects."
      ),
      source_name        = "IND"
    ),
    CONMED_CYP3A4_INH_STRONG = list(
      description        = "Concomitant strong CYP3A inhibitor coadministration indicator; 1 = on a strong CYP3A inhibitor, 0 = not.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant strong CYP3A inhibitor)",
      notes              = paste(
        "Parasrampuria 2025 Discussion identifies the agents behind this",
        "flag: 'The population PK model estimated impact of strong CYP3A",
        "inhibitors, primarily ritonavir (RTV) and cobicistat (COBI)'.",
        "Both are the pharmacokinetic boosters routinely coadministered",
        "with protease inhibitors in this heavily treatment-experienced",
        "population, which is why exposure to this covariate is high:",
        "293 of 764 population PK subjects (38%) (Table S2).",
        "Time-varying per record. Multiplicative on CL/F: 0.721x, i.e. a",
        "28% clearance decrease, so exposure rises. The model-based Ctau",
        "increase of 79% (Table 5, 775 vs 433 ng/mL) is consistent with",
        "the 44-88% Ctau increases seen in dedicated RTV, atazanavir/RTV",
        "and darunavir/RTV DDI studies."
      ),
      source_name        = "INH"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline.",
      units       = "years",
      type        = "continuous",
      notes       = "Tested on CL/F and V2/F (Table 1); not retained. Population PK cohort median 42 years, range 17-73 (Table S2). Parasrampuria 2025 flags this as a model limitation: only 11 subjects (1.4%) were 65 years or older."
    ),
    SEXF = list(
      description = "Female sex indicator; 1 = female, 0 = male.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested on CL/F, V2/F and Ka (Table 1, printed as 'gender'); not retained. Population PK cohort 216 female (28%), 548 male (72%) (Table S2). The source reports counts by sex but not the coding direction of its index variable; SEXF is the register canonical."
    ),
    RACE_BLACK = list(
      description = "Black / African American race indicator; 1 = Black or African American, 0 = other.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested on CL/F and V2/F (Table 1, as 'race'); not retained. Population PK cohort 177 (23%) (Table S2)."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator; 1 = Asian, 0 = other.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested on CL/F and V2/F (Table 1, as 'race'); not retained. Population PK cohort 5 subjects (1%) (Table S2) -- Parasrampuria 2025 names this explicitly as a model limitation ('few (N = 5) Asian subjects'), so the null result is an absence of information rather than evidence of no effect."
    ),
    HIV_POS = list(
      description = "HIV-1 infection status of the analysis subject; 1 = HIV-1-infected patient, 0 = healthy volunteer.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested on CL/F, V2/F and Ka (Table 1, as 'population (HIV-1 vs. healthy)'); not retained. Population PK cohort 606 HIV-1-infected (79%) and 158 healthy (21%) (Table S2). A published null of practical importance: it is what licenses pooling the healthy-volunteer phase 1 data with the patient data in one model."
    ),
    CRCL = list(
      description = "Baseline creatinine clearance.",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Tested on CL/F and V2/F (Table 1); not retained. Population PK cohort median 118 mL/min, range 5.26-271 (Table S2)."
    ),
    ALT = list(
      description = "Baseline alanine aminotransferase.",
      units       = "U/L",
      type        = "continuous",
      notes       = "Tested on CL/F and V2/F (Table 1); not retained. Population PK cohort median 24 U/L, range 6-240; Table S2 writes the unit as IU/L, used interchangeably with the canonical U/L."
    ),
    AST = list(
      description = "Baseline aspartate aminotransferase.",
      units       = "U/L",
      type        = "continuous",
      notes       = "Tested on CL/F and V2/F (Table 1); not retained. Population PK cohort median 25 U/L, range 10-288 (Table S2, printed as IU/L)."
    ),
    ALP = list(
      description = "Baseline alkaline phosphatase.",
      units       = "U/L",
      type        = "continuous",
      notes       = "Tested on CL/F (Table 1); not retained. Population PK cohort median 80 U/L, range 36-612 (Table S2, printed as IU/L). Mechanistically interesting and still null: intestinal alkaline phosphatase is the enzyme that hydrolyses the fostemsavir prodrug to temsavir, but serum ALP does not track that intestinal activity."
    ),
    DBIL = list(
      description = "Baseline direct (conjugated) bilirubin.",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Tested on CL/F (Table 1, as 'BILI', defined in the Table 1 abbreviations as direct bilirubin); not retained. Population PK cohort median 0.100 mg/dL, range 0.00-30.8 (Table S2). The register canonical DBIL is in umol/L; the source reports mg/dL and no conversion is applied here because the covariate is not in the model."
    ),
    CREAT = list(
      description = "Baseline serum creatinine.",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Tested on CL/F (Table 1, as 'SCr'); not retained. Not summarised for the population PK cohort in Table S2; the exposure-safety cohort median was 77 umol/L, range 31.8-1630."
    ),
    FASTED_STRICT = list(
      description = "Fasted-state indicator at the dose record; 1 = fasted, 0 = fed.",
      units       = "(binary)",
      type        = "binary",
      notes       = "NOT TESTED, for a structural reason rather than a null result: Parasrampuria 2025 Methods state 'The effect of prandial status was not tested as data were not collected in the Phase 3 study.' The food effect was instead brought in post hoc from a dedicated phase 1 food-effect study -- median Ctau under a standard meal was 1.68-fold the fasted value -- by multiplying model-predicted concentrations by 0.595 (= 1/1.68) to mimic fasting in the Table 5 simulations, rather than by fitting a relative-bioavailability parameter. The 0.595 factor is therefore an external post-processing constant and is deliberately NOT encoded in this model; see the vignette."
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "temsavir", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "temsavir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "temsavir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 764L,
    n_studies      = 7L,
    n_observations = "10,236 quantifiable plasma temsavir concentrations; BLQ records (< 5 ng/mL) were neither imputed nor included",
    age_range      = "17-73 years (median 42; 11 subjects (1.4%) aged 65 years or older) (Table S2)",
    weight_range   = "38-151 kg (median 72; baseline BMI median 24.9 kg/m2, range 14.4-52.2) (Table S2)",
    sex_female_pct = 28.3,
    race_ethnicity = c(White = 64.0, `Black or African American` = 23.0, Asian = 1.0, Other = 12.0),
    disease_state  = "pooled healthy adult volunteers (158, 21%) and adults with HIV-1 infection (606, 79%); the phase 2 and phase 3 patients were treatment-experienced, and the phase 3 BRIGHTE cohort was heavily treatment-experienced (HTE) with multidrug-resistant HIV-1 failing their current antiretroviral regimen",
    dose_range     = "fostemsavir extended-release tablets, oral: monotherapy 600-2400 mg twice daily and 1200 mg once daily; in combination with other HIV therapy 400-800 mg BID, 600-1200 mg every 12 h, 1200 mg at bedtime, and 600-1200 mg BID. Up to 10 days in the early-phase studies and 96 weeks or more in the late-phase studies. The approved regimen is 600 mg BID.",
    regions        = "multinational; the phase 3 BRIGHTE randomised cohort was 41% North America, 38% South America, 19% Europe, 3% other (Table S2)",
    co_medication  = "concomitant CYP3A inhibitor 293 subjects (38%), concomitant CYP3A inducer 37 subjects (5%) (Table S2)",
    notes          = paste0(
      "Seven pooled studies: four phase 1 studies in healthy volunteers ",
      "(including a thorough-QT study and DDI studies with darunavir/ritonavir ",
      "and etravirine), the phase 2a proof-of-concept AI438006 (NCT01009814), ",
      "the phase 2b combination-therapy study with a monotherapy sub-study ",
      "AI438011 (NCT01384734), and the phase 3 BRIGHTE study AI438047 ",
      "(NCT02362503); see Table S1. Sparse and/or intensive sampling ",
      "throughout. Bioanalysis by LC-MS/MS validated over 5-5000 ng/mL. ",
      "Estimation with NONMEM 7.2 FOCE-I; 500-replicate bootstrap and VPC. ",
      "All fixed and random effects were estimated with RSE < 15%. Model ",
      "diagnostics showed some overprediction of concentrations in the ",
      "phase 2a study (Figures S1, S2). The model-estimated CL/F of ",
      "51 L/h agreed with the 66.5 L/h CL/F from a dedicated absolute-",
      "bioavailability study."
    )
  )

  ini({
    # ==================================================================
    # Parasrampuria 2025 Table 2, "Estimate (%RSE) [Shrinkage%]" column,
    # with the bootstrap 95% CI (500 replicates) quoted alongside each
    # value. The covariate model is the Table 2 Note, verbatim:
    #
    #   CL/Fi = thetaCL * (WT/72)^0.75 * theta7^INH * theta8^IND * exp(eta1)
    #   V2/Fi = thetaV2 * (WT/72)             * exp(eta2)
    #   Q/Fi  = thetaQ  * (WT/72)^0.75
    #   V3/Fi = thetaV3 * (WT/72)
    #   Kai   = thetaKA * exp(eta3)
    #   DUR   = thetaDUR
    #   ETA (eta): ETA(eta_i) = SQRT(exp(eta_i) - 1)
    #   BSC: Off-diagonal element of the covariance between CL/F and V2/F
    #   Residual error: Y = log(IPRED) + (eps1) * exp(eta4)
    #
    # Note theta7 is the INHIBITOR multiplier and theta8 the INDUCER
    # multiplier, matching the Table 2 rows ("Effect of CYP3A Inhibitor
    # on CL/F" = 0.721, "Effect of CYP3A Inducer on CL/F" = 1.41). Both
    # enter as a base raised to a 0/1 indicator, so for a binary
    # covariate the form is exactly a multiplicative shift.
    #
    # VARIABILITY SCALE. The Table 2 Note supplies the transform for the
    # random-effect rows itself: the printed percentages are
    # 100 * sqrt(exp(omega^2) - 1), the exact log-normal CV. rxode2 wants
    # omega^2, so each row below is back-transformed as
    # omega^2 = log(1 + CV^2). The percentages are NOT omega, and they
    # are NOT variances.
    #
    # BSC IS A CORRELATION, NOT A COVARIANCE, despite the Table 2 Note
    # calling it "the off-diagonal element of the covariance". The
    # paper's own numbers refute the footnote: reading 0.609 as a
    # covariance implies a CL/F-V2/F correlation of
    # 0.609 / sqrt(0.1667645 * 0.2151920) = 3.21, which is impossible.
    # The refutation does not depend on how the diagonals are read --
    # taking the printed 42.6% and 49.0% as omega gives 2.92, and taking
    # them as variances gives 1.33; every reading exceeds 1. Read as a
    # correlation, 0.609 is unremarkable and its 6.54% RSE and tight
    # bootstrap CI (0.576, 0.645) are consistent with a well-estimated
    # correlation. This is the ordinary convention for a table that has
    # already transformed its diagonals to CV%: the off-diagonal is
    # transformed to the matching correlation. The covariance rxode2
    # needs is therefore 0.609 * omega_CL * omega_V2 = 0.1153672.
    # See the vignette Errata.
    # ==================================================================

    # ----- Structural disposition (apparent, i.e. CL/F, V/F) -----
    lcl <- log(51.0)  ; label("Apparent clearance CL/F at 72 kg with no CYP3A inducer or inhibitor (L/h)")           # Parasrampuria 2025 Table 2, CL/F = 51.0 L/h (RSE 2.16%), bootstrap 95% CI 49.1-52.9
    lvc <- log(257)   ; label("Apparent central volume of distribution V2/F at 72 kg (L)")                           # Parasrampuria 2025 Table 2, V2/F = 257 L (RSE 3.18%), bootstrap 95% CI 233-279
    lq  <- log(2.58)  ; label("Apparent inter-compartmental clearance Q/F at 72 kg (L/h)")                           # Parasrampuria 2025 Table 2, Q/F = 2.58 L/h (RSE 7.13%), bootstrap 95% CI 0.973-4.40
    lvp <- log(37.4)  ; label("Apparent peripheral volume of distribution V3/F at 72 kg (L)")                        # Parasrampuria 2025 Table 2, V3/F = 37.4 L (RSE 3.95%), bootstrap 95% CI 26.5-65.6

    # ----- Absorption: zero-order input into the depot, then first-order -----
    # The two processes are SEQUENTIAL, not parallel. Results 3.1: the
    # base model has "dual zero- and first-order absorption (zero order
    # modeled as an estimated constant infusion of duration DUR in
    # absorption depot compartment)". A parallel-pathway model would need
    # a dose-splitting fraction; Table 2 estimates none, and the Figure 1
    # schematic shows one path (dose -> depot -> central).
    lka <- log(2.33)  ; label("First-order absorption rate constant Ka from depot to central (1/h)")                 # Parasrampuria 2025 Table 2, Ka = 2.33 1/h (RSE 13.3%), bootstrap 95% CI 1.84-2.79
    ld1 <- log(3.84)  ; label("Duration of the zero-order input into the depot, DUR (h)")                            # Parasrampuria 2025 Table 2, DUR = 3.84 h (RSE 2.51%), bootstrap 95% CI 3.68-3.99

    # ----- Covariate effects -----
    # Allometry: both exponents carry "Fixed (NA)" in Table 2 with no RSE
    # and no bootstrap CI, so both are wrapped in fixed(). The 0.75 / 1
    # pair is the standard theory-based allometric choice.
    e_wt_cl_q  <- fixed(0.75) ; label("Body-weight allometric exponent shared by CL/F and Q/F, referenced to 72 kg (unitless)")   # Parasrampuria 2025 Table 2, "Effect of WT on CL/F and Q/F" = 0.75 Fixed (NA)
    e_wt_vc_vp <- fixed(1)    ; label("Body-weight allometric exponent shared by V2/F and V3/F, referenced to 72 kg (unitless)")  # Parasrampuria 2025 Table 2, "Effect of WT on V2/F and V3/F" = 1 Fixed (NA)

    # CYP3A perpetrator effects, both multiplicative on CL/F and both
    # estimated (they are theta8 and theta7 in the Table 2 Note).
    e_conmed_cyp3a4_ind_mod_cl    <- 1.41  ; label("Multiplicative effect of a concomitant moderate CYP3A inducer on CL/F (unitless; 1.41 -> 41% higher CL/F)")     # Parasrampuria 2025 Table 2, "Effect of CYP3A Inducer on CL/F" = 1.41 (RSE 2.39%), bootstrap 95% CI 1.25-1.58; Results 3.1 quotes it as a "41% increase"
    e_conmed_cyp3a4_inh_strong_cl <- 0.721 ; label("Multiplicative effect of a concomitant strong CYP3A inhibitor on CL/F (unitless; 0.721 -> 28% lower CL/F)")     # Parasrampuria 2025 Table 2, "Effect of CYP3A Inhibitor on CL/F" = 0.721 (RSE 1.57%), bootstrap 95% CI 0.687-0.764; Results 3.1 quotes it as a "28% decrease"

    # ----- Inter-individual variability -----
    # Diagonals back-transformed from the Table 2 CV% via
    # omega^2 = log(1 + CV^2); off-diagonal from the correlation 0.609
    # (see the BSC note above):
    #   omega^2(CL/F) = log(1 + 0.426^2) = 0.1667645  -> omega 0.4083681
    #   omega^2(V2/F) = log(1 + 0.490^2) = 0.2151920  -> omega 0.4638879
    #   cov           = 0.609 * 0.4083681 * 0.4638879 = 0.1153672
    etalcl + etalvc ~ c(0.1667645,
                        0.1153672, 0.2151920)  # Parasrampuria 2025 Table 2: ETA(CL/F)% = 42.6 (RSE 3.78%, shrinkage 11%), 95% CI 39.8-45.5; ETA(V2/F)% = 49.0 (RSE 5.59%, shrinkage 35%), 95% CI 38.8-56.8; BSC(CL/F, V2/F) = 0.609 (RSE 6.54%), 95% CI 0.576-0.645, read as a correlation

    # omega^2(Ka) = log(1 + 1.27^2) = 0.9604607. The largest random
    # effect in the model by a wide margin, with 43% shrinkage -- most
    # of the pooled data is sparsely sampled and carries little
    # information about the absorption phase.
    etalka ~ 0.9604607  # Parasrampuria 2025 Table 2: ETA(Ka)% = 127 (RSE 8.53%, shrinkage 43%), 95% CI 112-146

    # ----- Residual error -----
    # Y = log(IPRED) + eps1 * exp(eta4): additive on the natural-log
    # scale, i.e. log-normal in linear space, which is nlmixr2's
    # `lnorm(expSd)`. The residual MAGNITUDE itself carries exponential
    # IIV (eta4), so the individual log-scale SD is expSd * exp(etaexpSd)
    # -- see model(). Following the Chandasana 2024 / Tang 2023 /
    # Yamamoto 2023 precedent for IIV on residual error.
    expSd <- 0.613  ; label("Residual error SD, additive on the natural-log scale (unitless)")   # Parasrampuria 2025 Table 2, "Residual Additive SD on Log Scale" = 0.613 (RSE 1.69%), bootstrap 95% CI 0.595-0.628

    # omega^2(residual) = log(1 + 0.332^2) = 0.1045618
    etaexpSd ~ 0.1045618  # Parasrampuria 2025 Table 2: ETA (Residual)% = 33.2 (RSE 3.59%, shrinkage 13%), 95% CI 30.1-35.8
  })

  model({
    # ==================================================================
    # 1. Individual parameters (Table 2 Note, verbatim)
    # ==================================================================
    # The CYP3A terms are written as a base raised to a 0/1 indicator,
    # exactly as printed (theta7^INH * theta8^IND). For a binary
    # covariate this equals a multiplicative shift when the flag is 1
    # and leaves CL/F untouched when it is 0. A subject on both a
    # moderate inducer and a strong inhibitor gets the product,
    # 1.41 * 0.721 = 1.017, i.e. near-cancellation -- which is what
    # Table 5's "600 mg BID + CYP3A inducer + CYP3A inhibitor" row
    # shows (median Ctau 414 vs 433 ng/mL dosed alone).
    cl <- exp(lcl + etalcl) * (WT / 72)^e_wt_cl_q *
      e_conmed_cyp3a4_inh_strong_cl^CONMED_CYP3A4_INH_STRONG *
      e_conmed_cyp3a4_ind_mod_cl^CONMED_CYP3A4_IND_MOD
    vc <- exp(lvc + etalvc) * (WT / 72)^e_wt_vc_vp
    q  <- exp(lq)           * (WT / 72)^e_wt_cl_q
    vp <- exp(lvp)          * (WT / 72)^e_wt_vc_vp

    ka <- exp(lka + etalka)
    d1 <- exp(ld1)

    # ==================================================================
    # 2. Micro-constants
    # ==================================================================
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ==================================================================
    # 3. ODE system
    # ==================================================================
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ==================================================================
    # 4. Zero-order input into the depot
    # ==================================================================
    # Dose records must carry rate = -2 so rxode2 uses this modelled
    # duration; a plain bolus into `depot` would silently drop the
    # zero-order half of the absorption model and shift Tmax earlier.
    dur(depot) <- d1

    # ==================================================================
    # 5. Observation and error
    # ==================================================================
    # Dose is in mg and vc in L, so central/vc is mg/L; the factor of
    # 1000 converts to the ng/mL in which Parasrampuria 2025 reports
    # temsavir concentrations (assay range 5-5000 ng/mL).
    Cc <- 1000 * central / vc

    # Individual residual SD on the log scale: eps1 * exp(eta4) in the
    # Table 2 Note.
    expSdi <- expSd * exp(etaexpSd)
    Cc ~ lnorm(expSdi)
  })
}
