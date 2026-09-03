Gaffney_2026_niraparib <- function() {
  description <- paste(
    "Three-compartment population PK model for oral niraparib with a",
    "three-transit-compartment absorption chain, fitted to pooled data from",
    "six studies in patients with advanced solid tumours or ovarian cancer.",
    "Baseline albumin, alkaline phosphatase and creatinine clearance act on",
    "apparent clearance; baseline albumin and body weight on apparent",
    "central volume; baseline albumin on the first apparent peripheral",
    "volume; and prandial state on relative bioavailability and mean",
    "transit time. This is the niraparib reference population PK model."
  )
  reference <- paste(
    "Gaffney A, Franchetti Y, Desrosiers M, Trame MN, Jewell RC.",
    "Population Pharmacokinetic Modeling of Niraparib to Assess Different",
    "Absorption Models. J Clin Pharmacol. 2026;66(5):e70210.",
    "doi:10.1002/jcph.70210."
  )
  vignette <- "Gaffney_2026_niraparib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")
  # The oral dose enters the head of the transit chain only. Without this
  # field buildModelDb()'s two-name heuristic reports "depot,central", which
  # is truthy but wrong: `central` is never dosed in this model.
  dosing <- "depot"

  covariateData <- list(
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed baseline value (source column ALBBL). Power effects on",
        "three parameters, all normalised to the same reference: CL/F",
        "exponent 0.742, Vc/F exponent 0.363 and Vp1/F exponent 1.02",
        "(Table 3). The paper reports albumin in the US convention g/dL and",
        "states the reference patient's value as 4 g/dL (Figure 5 legend);",
        "the canonical column is SI g/L, so the reference is written here as",
        "40 g/L. The effect enters only as the ratio (ALB / 40), which is",
        "numerically identical to the paper's (ALBBL / 4) in g/dL, so no",
        "value conversion of the coefficient is needed. Reference confirmed",
        "by reproducing the paper's own steady-state AUC ratios: at the 95th",
        "percentile 4.7 g/dL, (4.7/4)^-0.742 = 0.887 vs the reported 0.89;",
        "at the 5th percentile 3.1 g/dL, 1.208 vs 1.21; at the minimum",
        "1.7 g/dL, 1.887 vs 1.89; at the maximum 7.9 g/dL, 0.604 vs 0.60",
        "(all from the Results 'Final Reference Niraparib Population PK",
        "Model' section). Analysis-population range 1.7 to 7.9 g/dL",
        "(17 to 79 g/L)."
      ),
      source_name        = "ALBBL"
    ),
    ALP = list(
      description        = "Baseline serum alkaline phosphatase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed baseline value (source column ALPBL); the paper reports",
        "it in IU/L, which is interchangeable with the canonical U/L. Power",
        "effect on CL/F with exponent -0.074 normalised to 83 U/L (Table 3;",
        "reference value from the Figure 5 legend). ALP was selected as the",
        "single liver-function test carried forward from a set of four",
        "candidates (AST, ALT, ALP, bilirubin) because niraparib is known to",
        "be hepatically eliminated. Reference and exponent confirmed by",
        "reproducing the paper's AUC ratios: (49/83)^0.074 = 0.962 vs the",
        "reported 0.96; (211.8/83)^0.074 = 1.072 vs 1.07; (0.9/83)^0.074 =",
        "0.715 vs 0.72; (1814/83)^0.074 = 1.256 vs 1.25. Analysis-population",
        "range 0.9 to 1814 IU/L."
      ),
      source_name        = "ALPBL"
    ),
    CRCL = list(
      description        = "Baseline creatinine clearance",
      units              = "mL/min (not BSA-normalized)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed baseline value (source column CrCLBL), reported in raw",
        "mL/min rather than the canonical column's BSA-normalized",
        "mL/min/1.73 m^2; the paper does not name the estimating equation.",
        "Power effect on CL/F with exponent 0.287 normalised to",
        "82.37 mL/min (Table 3; reference value from the Figure 5 legend).",
        "Reference and exponent confirmed by reproducing the paper's AUC",
        "ratios: (147.4/82.37)^-0.287 = 0.846 vs the reported 0.85;",
        "(44.3/82.37)^-0.287 = 1.195 vs 1.20; (27.6/82.37)^-0.287 = 1.369 vs",
        "1.37; (298.6/82.37)^-0.287 = 0.691 vs 0.69. Analysis-population",
        "range 27.6 to 298.6 mL/min; 333 patients (19.7%) had moderate renal",
        "impairment (30 to 59 mL/min) and only 4 (0.2%) were below",
        "30 mL/min."
      ),
      source_name        = "CrCLBL"
    ),
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed baseline value (source column WTBL). Power effect on",
        "Vc/F only, exponent 0.577 normalised to 70 kg (Table 3). The 70 kg",
        "reference is stated twice in the Results: 'Typical values of CL/F",
        "and Vc/F for a 70 kg patient were 15.9 L/h and 450 L' and in the",
        "Figure 5 reference-patient legend. Body weight was tested and",
        "rejected on Vp1/F, Vp2/F and F1 (Table 2), so there is no effect of",
        "weight on AUC - the paper makes this point explicitly."
      ),
      source_name        = "WTBL"
    ),
    FED = list(
      description        = "Fed state at the time of the dose",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted or unknown prandial state)",
      notes              = paste(
        "Per-dose-record indicator. The reference level pools the fasted and",
        "the unknown prandial states: the Covariate Analysis section states",
        "that 'the unknown prandial state was grouped with the reference",
        "(fasted) prandial state for mean transit time (MTT) and relative",
        "bioavailability (F1) to avoid negative parameter estimates and",
        "preserve interpretability', and the Table 2 footnote confirms that",
        "after a refinement step the fed state was compared with the",
        "fasted/unknown state as reference. Two multiplicative fractional",
        "effects, both from Table 3: F1 * (1 + 0.236 * FED) and",
        "MTT * (1 + 1.15 * FED). The only deliberately fed data come from",
        "the TABLET study's high-fat-meal food-effect arm, but the paper",
        "applies the flag generically across all six pooled studies (every",
        "record outside that arm is fasted or unknown), so the general FED",
        "canonical applies rather than FED_HIGHFAT."
      ),
      source_name        = "prandial state"
    ),
    STUDY_NOVA = list(
      description        = "NOVA (NCT01847274) study-cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not NOVA; PN001 when all five study indicators are 0)",
      notes              = paste(
        "Selects the NOVA study-specific residual-error magnitude (0.350;",
        "Table 3). Phase 3, recurrent ovarian cancer, 405 patients, 2054",
        "included observations (Table 1). A separate additive residual error",
        "on log-transformed concentrations was estimated for each of the six",
        "pooled studies; PN001 is encoded here as the reference level, so",
        "exactly one of the five study indicators is 1 per record and all",
        "five are 0 for a PN001 record."
      ),
      source_name        = "study"
    ),
    STUDY_QUADRA = list(
      description        = "QUADRA (NCT02354586) study-cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not QUADRA; PN001 when all five study indicators are 0)",
      notes              = paste(
        "Selects the QUADRA study-specific residual-error magnitude (0.381;",
        "Table 3). Phase 2, advanced relapsed ovarian cancer, 455 patients,",
        "1410 included observations (Table 1)."
      ),
      source_name        = "study"
    ),
    STUDY_PRIMA = list(
      description        = "PRIMA (NCT02655016) study-cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not PRIMA; PN001 when all five study indicators are 0)",
      notes              = paste(
        "Selects the PRIMA study-specific residual-error magnitude (0.451;",
        "Table 3). Phase 3, first-line advanced ovarian cancer, 480",
        "patients, 1856 included observations (Table 1)."
      ),
      source_name        = "study"
    ),
    STUDY_TABLET = list(
      description        = "TABLET (NCT03329001) study-cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not TABLET; PN001 when all five study indicators are 0)",
      notes              = paste(
        "Selects the TABLET study-specific residual-error magnitude (0.324;",
        "Table 3). Phase 1 tablet-versus-capsule bioequivalence and food",
        "effect study in advanced solid tumours, 225 patients, 6487 included",
        "observations (46% of the analysis dataset; Table 1). One of the two",
        "intensively sampled studies added in this analysis, and the source",
        "of essentially all the fed records."
      ),
      source_name        = "study"
    ),
    STUDY_HEPATIC = list(
      description        = "HEPATIC (NCT03359850) study-cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not HEPATIC; PN001 when all five study indicators are 0)",
      notes              = paste(
        "Selects the HEPATIC study-specific residual-error magnitude (0.146;",
        "Table 3). Phase 1 hepatic-impairment study in advanced solid",
        "tumours, 17 patients (8 with moderate hepatic impairment, 9 with",
        "normal hepatic function), 203 included observations (Table 1). The",
        "other intensively sampled study added in this analysis. Hepatic",
        "impairment itself was NOT retained as a structural or covariate",
        "effect: the paper states the ability to identify a difference was",
        "limited by the 8-patient moderate-impairment subgroup, so the study",
        "enters the model only through its residual-error magnitude."
      ),
      source_name        = "study"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened on CL/F, Vc/F, Vp1/F and Vp2/F in the stepwise covariate",
        "analysis and rejected on all four (Table 2). Not retained in the",
        "final model and no point estimate is reported."
      )
    ),
    AGE = list(
      description = "Baseline age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened on CL/F and Vc/F and rejected on both (Table 2), even",
        "though age had been identified as a covariate on CL/F in the",
        "previous niraparib population PK model. No point estimate is",
        "reported."
      )
    ),
    PLT = list(
      description = "Baseline platelet count",
      units       = "10^9/L",
      type        = "continuous",
      notes       = paste(
        "Screened on CL/F and rejected (Table 2). Tested because baseline",
        "platelet count is used in determining the niraparib starting dose.",
        "No point estimate is reported."
      )
    ),
    ECOG = list(
      description = "Baseline ECOG performance status",
      units       = "(ordinal 0-4)",
      type        = "categorical",
      notes       = paste(
        "Screened on CL/F and rejected (Table 2). No point estimate is",
        "reported."
      )
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "niraparib", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    transit1 = list(
      analyte = "niraparib", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    transit2 = list(
      analyte = "niraparib", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "niraparib", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "niraparib", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral2 = list(
      analyte = "niraparib", units = "mg",
      specimen = "plasma", verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1686,
    n_studies      = 6,
    n_observations = 14106,
    sex_female_pct = 91.6,
    race_ethnicity = c(White = 86.4),
    disease_state  = paste(
      "Advanced solid tumours or haematologic malignancies (PN001),",
      "platinum-sensitive ovarian cancer (NOVA, PRIMA), relapsed high-grade",
      "serous ovarian cancer (QUADRA), and advanced solid tumours with",
      "normal hepatic function or moderate hepatic impairment (HEPATIC) or",
      "in a tablet-versus-capsule bioequivalence and food-effect design",
      "(TABLET). 83.6% of the pooled analysis population had ovarian cancer",
      "and 86.5% had normal hepatic function."
    ),
    dose_range     = paste(
      "Oral niraparib. HEPATIC and TABLET each administered a single 300 mg",
      "dose; the remaining four studies were the monotherapy dose-escalation",
      "(PN001) and once-daily maintenance (NOVA, QUADRA, PRIMA) programmes",
      "described in the cited primary study reports."
    ),
    hepatic_function = paste(
      "1459 of 1686 patients (86.5%) had normal hepatic function; 8 patients",
      "in the HEPATIC study had moderate hepatic impairment."
    ),
    renal_function = paste(
      "Baseline creatinine clearance 27.6 to 298.6 mL/min; 333 patients",
      "(19.7%) had moderate renal impairment (30 to 59 mL/min) and 4 (0.2%)",
      "were below 30 mL/min."
    ),
    notes          = paste(
      "Six pooled studies: PN001 (NCT00749502), NOVA (NCT01847274), QUADRA",
      "(NCT02354586), PRIMA (NCT02655016), HEPATIC (NCT03359850) and TABLET",
      "(NCT03329001); per-study patient and observation counts in Table 1.",
      "1988 records (12.4%) were excluded before fitting and a further 111",
      "(0.8%) with |CWRES| > 4 were excluded during model building, leaving",
      "13,995 concentrations for the base model. Age, weight and race",
      "distributions are tabulated in Supplementary Table S1, which is not",
      "in the open-access deposit; the percentages recorded here are the",
      "ones quoted in the Results text."
    )
  )

  ini({
    # ---- Structural disposition parameters (Table 3, Estimates column).
    #      All are apparent (/F) values; F1 is fixed to 1 so the model
    #      carries no absolute bioavailability.
    lcl  <- log(15.9); label("Apparent clearance CL/F (L/h)")                                  # Table 3 CL/F 15.9 L/h (RSE 1.3%; bootstrap 15.9, 95% CI 15.5-16.3)
    lvc  <- log(450);  label("Apparent central volume of distribution Vc/F (L)")               # Table 3 Vc/F 450 L (RSE 2.6%; bootstrap 449, 95% CI 422-478)
    lq   <- log(43.8); label("Apparent first intercompartmental clearance Q1/F (L/h)")         # Table 3 Q1/F 43.8 L/h (RSE 5.4%; bootstrap 44, 95% CI 37.9-51.5). The Discussion quotes 43.5 L/h for the same parameter; Table 3 is used - see vignette Errata
    lvp  <- log(395);  label("Apparent first peripheral volume of distribution Vp1/F (L)")     # Table 3 Vp1/F 395 L (RSE 3.2%; bootstrap 394, 95% CI 362-429)
    lq2  <- log(2.19); label("Apparent second intercompartmental clearance Q2/F (L/h)")        # Table 3 Q2/F 2.19 L/h (RSE 6.1%; bootstrap 2.21, 95% CI 1.79-2.78)
    lvp2 <- log(361);  label("Apparent second peripheral volume of distribution Vp2/F (L)")    # Table 3 Vp2/F 361 L (RSE 4.8%; bootstrap 363, 95% CI 327-402)

    # ---- Absorption. The chain is parameterised by its mean transit time;
    #      the single first-order rate constant the paper calls Ka is
    #      derived from it in model().
    lmtt <- log(1.78); label("Mean transit time MTT (h)")                                      # Table 3 MTT 1.78 h (RSE 2.4%; bootstrap 1.77, 95% CI 1.66-1.88)

    # ---- Relative bioavailability, fixed to 1 as a structural anchor but
    #      carrying estimated IIV and a prandial-state effect.
    lfdepot <- fixed(log(1)); label("Relative bioavailability F1 (unitless)")                  # Table 3 "F1 (-) 1 Fixed"; Base Population PK Model section: "Relative bioavailability was fixed to 1"

    # ---- Covariate effects. Continuous covariates enter as power terms
    #      normalised to a reference value, Pki = theta_k * (Xij / M(Xj))^theta_j;
    #      the binary prandial-state covariate enters as a fractional change,
    #      Pki = theta_k * (1 + theta_j)^Xij, i.e. theta_k * (1 + theta_j * Xij)
    #      for a 0/1 indicator (both forms printed on page 4).
    e_alb_cl     <- 0.742;  label("Power exponent on (ALB / 40 g/L) for CL/F (unitless)")      # Table 3 CL/F ~ ALBBL exponent 0.742 (RSE 11.8%; bootstrap 0.734, 95% CI 0.567-0.923)
    e_alp_cl     <- -0.074; label("Power exponent on (ALP / 83 U/L) for CL/F (unitless)")      # Table 3 CL/F ~ ALPBL exponent -0.074 (RSE 18.4%; bootstrap -0.0756, 95% CI -0.118 to -0.0383)
    e_crcl_cl    <- 0.287;  label("Power exponent on (CRCL / 82.37 mL/min) for CL/F (unitless)") # Table 3 CL/F ~ CrCLBL exponent 0.287 (RSE 8.9%; bootstrap 0.286, 95% CI 0.241-0.337)
    e_alb_vc     <- 0.363;  label("Power exponent on (ALB / 40 g/L) for Vc/F (unitless)")      # Table 3 Vc/F ~ ALBBL exponent 0.363 (RSE 28.6%; bootstrap 0.369, 95% CI 0.11-0.609)
    e_wt_vc      <- 0.577;  label("Power exponent on (WT / 70 kg) for Vc/F (unitless)")        # Table 3 Vc/F ~ WTBL exponent 0.577 (RSE 11.7%; bootstrap 0.58, 95% CI 0.455-0.697)
    e_alb_vp     <- 1.02;   label("Power exponent on (ALB / 40 g/L) for Vp1/F (unitless)")     # Table 3 Vp1/F ~ ALBBL exponent 1.02 (RSE 13.7%; bootstrap 1.01, 95% CI 0.643-1.37)
    e_fed_fdepot <- 0.236;  label("Fractional change in F1 when dosed fed (unitless)")         # Table 3 F1 ~ fed fractional change 0.236 (RSE 9.6%; bootstrap 0.231, 95% CI 0.156-0.333)
    e_fed_mtt    <- 1.15;   label("Fractional change in MTT when dosed fed (unitless)")        # Table 3 MTT ~ fed fractional change 1.15 (RSE 2.7%; bootstrap 1.2, 95% CI 0.641-1.82)

    # ---- Interindividual variability. Table 3 reports each IIV as a CV%
    #      with the conversion stated in the table note: CV% =
    #      sqrt(exp(omega^2) - 1) * 100, so omega^2 = log(CV^2 + 1). The
    #      variance-covariance structure investigated during covariate
    #      modelling is not reported as off-diagonal estimates, so the
    #      diagonal form of Table 3 is encoded here.
    etalcl     ~ 0.056002  # Table 3 IIV on CL/F 24.0% CV  (RSE 3.2%; bootstrap 23.9, 95% CI 21.4-26.3; shrinkage 35.1%) -> log(0.240^2 + 1)
    etalvc     ~ 0.055549  # Table 3 IIV on Vc/F 23.9% CV  (RSE 5.8%; bootstrap 23.7, 95% CI 17.8-29.8; shrinkage 59.4%) -> log(0.239^2 + 1)
    etalmtt    ~ 0.273624  # Table 3 IIV on MTT 56.1% CV   (RSE 2.2%; bootstrap 56.1, 95% CI 50.2-63.3; shrinkage 28.0%) -> log(0.561^2 + 1)
    etalvp2    ~ 0.254211  # Table 3 IIV on Vp2/F 53.8% CV (RSE 6.7%; bootstrap 53.6, 95% CI 40.6-71.1; shrinkage 68.4%) -> log(0.538^2 + 1)
    etalfdepot ~ 0.108176  # Table 3 IIV on F1 33.8% CV    (RSE 2.8%; bootstrap 33.7, 95% CI 31.8-36; shrinkage 27.1%) -> log(0.338^2 + 1)

    # ---- Residual error. The Table 3 note states that "an additive
    #      residual error model on log-transformed concentrations was used",
    #      which is nlmixr2's `~ lnorm(expSd)` with the SD on the log scale.
    #      One magnitude was estimated per pooled study; PN001 is the
    #      reference level carried by the bare `expSd` name.
    expSd         <- 0.227; label("Log-scale residual SD, PN001 (unitless)")                   # Table 3 Res variability CV PN001 0.227 (RSE 0.8%; bootstrap 0.225, 95% CI 0.206-0.245; shrinkage 11.3%)
    expSdNova     <- 0.350; label("Log-scale residual SD, NOVA (unitless)")                    # Table 3 Res variability CV NOVA 0.350 (RSE 0.9%; bootstrap 0.35, 95% CI 0.315-0.385)
    expSdQuadra   <- 0.381; label("Log-scale residual SD, QUADRA (unitless)")                  # Table 3 Res variability CV QUADRA 0.381 (RSE 1.5%; bootstrap 0.381, 95% CI 0.34-0.422)
    expSdPrima    <- 0.451; label("Log-scale residual SD, PRIMA (unitless)")                   # Table 3 Res variability CV PRIMA 0.451 (RSE 1.6%; bootstrap 0.451, 95% CI 0.421-0.477)
    expSdTablet   <- 0.324; label("Log-scale residual SD, TABLET (unitless)")                  # Table 3 Res variability CV TABLET 0.324 (RSE 0.3%; bootstrap 0.323, 95% CI 0.305-0.341)
    expSdHepatic  <- 0.146; label("Log-scale residual SD, HEPATIC (unitless)")                 # Table 3 Res variability CV HEPATIC 0.146 (RSE 4.2%; bootstrap 0.146, 95% CI 0.121-0.171)
  })

  model({
    # ---- 1. Individual disposition parameters. Every continuous covariate
    #         is a power term normalised to its reference value (page 4
    #         "Continuous covariates" equation). ALB is held in the
    #         canonical SI unit g/L, so the paper's 4 g/dL reference is
    #         written as 40 g/L; the ratio, and therefore the effect, is
    #         identical either way.
    cl <- exp(lcl + etalcl) *
      (ALB / 40)^e_alb_cl *
      (ALP / 83)^e_alp_cl *
      (CRCL / 82.37)^e_crcl_cl
    vc  <- exp(lvc + etalvc) * (ALB / 40)^e_alb_vc * (WT / 70)^e_wt_vc
    q   <- exp(lq)
    vp  <- exp(lvp) * (ALB / 40)^e_alb_vp
    q2  <- exp(lq2)
    vp2 <- exp(lvp2 + etalvp2)

    # ---- 2. Absorption. The prandial-state effects are the categorical
    #         fractional-change form of page 4: P = theta_k * (1 + theta_j)
    #         when FED = 1 and theta_k when FED = 0.
    mtt <- exp(lmtt + etalmtt) * (1 + e_fed_mtt * FED)

    # The absorption chain (Figure 3) is: oral dose -> Transit -> Transit ->
    # Transit -> Central, with the SAME first-order rate constant (the
    # paper's "Ka, absorption/transfer rate constant") on all three arrows
    # and no separate depot or lag time. There are therefore exactly three
    # first-order steps between the dose record and the central compartment,
    # so the arrival time in `central` is Erlang with shape 3 and rate ktr,
    # whose mean - the mean transit time - is 3 / ktr. Table 3 reports no Ka,
    # which is what pins Ka == ktr rather than a separate terminal rate.
    ktr <- 3 / mtt

    fdepot <- exp(lfdepot + etalfdepot) * (1 + e_fed_fdepot * FED)

    # ---- 3. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # ---- 4. ODE system. The paper's three "Transit" boxes are encoded as
    #         depot + transit1 + transit2 so that the dose lands in the
    #         canonical extravascular dosing compartment; the three ktr
    #         steps are depot -> transit1, transit1 -> transit2 and
    #         transit2 -> central, matching Figure 3 arrow for arrow. The
    #         same depot-plus-(n-1)-transit idiom is used in
    #         Comisar_2025_rimegepant.
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(central)     <-  ktr * transit2 - kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # ---- 5. Bioavailability.
    f(depot) <- fdepot

    # ---- 6. Observation and error. `central` is in mg and `vc` in L, so
    #         central/vc is mg/L = ug/mL; the paper reports plasma niraparib
    #         in ng/mL (Figure 4 y-axis), hence the factor of 1000.
    Cc <- central / vc * 1000

    # A separate residual magnitude per pooled study (Base Population PK
    # Model section: "a separate additive residual error model for each
    # study was estimated"). PN001 is the reference: a record with all five
    # study indicators at 0 gets the PN001 magnitude.
    expSdCc <- expSd * (1 - STUDY_NOVA - STUDY_QUADRA - STUDY_PRIMA -
                          STUDY_TABLET - STUDY_HEPATIC) +
      expSdNova * STUDY_NOVA +
      expSdQuadra * STUDY_QUADRA +
      expSdPrima * STUDY_PRIMA +
      expSdTablet * STUDY_TABLET +
      expSdHepatic * STUDY_HEPATIC

    Cc ~ lnorm(expSdCc)
  })
}
