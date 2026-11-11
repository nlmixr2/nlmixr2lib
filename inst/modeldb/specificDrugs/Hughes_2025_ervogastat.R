Hughes_2025_ervogastat <- function() {
  description <- "Two-compartment population pharmacokinetic model for oral ervogastat (PF-06865571, a diacylglycerol O-acyltransferase 2 inhibitor) in healthy participants, participants with hepatic impairment, and adults with non-alcoholic fatty liver disease (NAFLD) or presumed non-alcoholic steatohepatitis (NASH). Absorption is parallel first-order (up to a capped amount of each dose) and lagged zero-order. Apparent clearance scales allometrically with body weight and with the bioavailability-adjusted total daily dose, and is lower with hepatic impairment and in NAFLD / presumed NASH; the oral suspension raises the first-order amount cap and lowers relative bioavailability. Interindividual variability on CL/F and Vc/F; combined additive and proportional residual error with a larger proportional component in the sparsely sampled phase 2 studies. This is the ervogastat exposure model whose individual average steady-state concentrations drive the Hughes 2025 liver-fat exposure-response model (Hughes_2025_ervogastat_clesacostat_liverfat)."
  reference <- paste(
    "Hughes JH, Amin NB, Wojciechowski J, Vourvahis M.",
    "Exposure-response modeling of liver fat imaging endpoints in",
    "non-alcoholic fatty liver disease populations administered ervogastat",
    "alone and co-administered with clesacostat.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(2):317-330.",
    "doi:10.1002/psp4.13275",
    sep = " "
  )
  vignette <- "Hughes_2025_ervogastat_clesacostat"
  paper_specific_residual_sds <- c("propSd_phase2")
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed baseline body weight. Allometric power function referenced",
        "to 70 kg with exponents fixed to 0.75 on CL/F and Q/F and 1 on Vc/F",
        "and Vp/F (Supplement Table S3; Code S1 ALLMCL / ALLMVD). Code S1",
        "applies the scaling only when BWT is nonzero."
      ),
      source_name = "BWT"
    ),
    DOSE_ERVOGASTAT_MG = list(
      description = "Ervogastat amount per administration",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source DOSE1. Sets the first-order fraction FK1 = AK1 / DOSE1, with",
        "FK1 = 0.999 when AK1 exceeds the dose (Code S1; Supplement Equation",
        "S1 prints the cap as 1, the control stream as 0.999). Must equal the",
        "amt of both dose records of the administration."
      ),
      source_name = "DOSE1"
    ),
    DOSE_ERVOGASTAT_MGD = list(
      description = "Total daily ervogastat dose",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source DOST1 ('actual total-daily ervogastat dose'). CL/F scales as",
        "(DOST1 * FORMFREL / 300)^-0.156342, i.e. on the bioavailability-",
        "adjusted daily dose referenced to 300 mg (Code S1 EFFDOST / DOSTCL;",
        "Table S3 'Effect of total daily dose on CL/F; ref. 300 mg'). The",
        "dose effect is skipped when DOST1 is 0. Equals twice the",
        "per-administration dose for twice-daily and three times it for",
        "every-8-hour regimens."
      ),
      source_name = "DOST1"
    ),
    FORM_ERVOGASTAT_SUSP = list(
      description = "Ervogastat immediate-release oral suspension indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (any tablet: immediate release, slow or fast controlled release, or modified release)",
      notes = paste(
        "1 = immediate-release oral suspension (Code S1 FORM = 1), 0 = any",
        "tablet (FORM = 2-5). Raises the first-order amount cap AK1 by",
        "233.717% and lowers relative bioavailability by 6.55% (Table S3).",
        "The suspension was used in the single- and multiple-ascending-dose",
        "studies (NCT03092232, NCT03230383) and one arm of NCT03372044",
        "(Table S1)."
      ),
      source_name = "FORM"
    ),
    HEPIMP = list(
      description = "Any hepatic impairment indicator (mild, moderate or severe)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (normal hepatic function)",
      notes = paste(
        "1 = mild, moderate or severe hepatic impairment (Code S1 PTST = 1, 2",
        "or 3), pooled into one effect on CL/F (-26.874%, Table S3 'Effect of",
        "any hepatic impairment on CL/F'). Six participants per impairment",
        "category came from the dedicated hepatic-impairment study NCT04091061",
        "(Table 1, Table S1). The paper does not state the classification",
        "system used to grade impairment."
      ),
      source_name = "PTST"
    ),
    DIS_NAFLD = list(
      description = "Non-alcoholic fatty liver disease or presumed non-alcoholic steatohepatitis patient indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy participant or participant in a hepatic-impairment cohort)",
      notes = paste(
        "1 = NAFLD or presumed NASH (Code S1 PTST = 4 or 5), 0 otherwise.",
        "The two categories were combined into one estimate (-18.4998% on",
        "CL/F, Table S3) 'based on the similarity of non-parametric",
        "distributions and physiological plausibility' (Supplement Results).",
        "Code S1 also defines PTST = 6 (confirmed NASH), which is not part",
        "of the effect; no confirmed-NASH participant appears in Table 1."
      ),
      source_name = "PTST"
    ),
    STUDY_PHASE2 = list(
      description = "Phase 2 study indicator (residual-error stratum)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (phase 1 and phase 1b studies, densely sampled)",
      notes = paste(
        "1 = record from one of the phase 2 studies with sparse PK sampling",
        "(Code S1 PROT = 1171002, 3711001 or 3711005, i.e. NCT03248882,",
        "NCT03776175 and NCT04399538; only the last two dosed ervogastat).",
        "Selects the proportional residual SD only: 0.711411 for these",
        "studies versus 0.397557 for the dense designs (Table S3 rows",
        "'Proportional residual error for dense / sparse study designs')."
      ),
      source_name = "PROT"
    )
  )

  covariatesDataExcluded <- list(
    RACE_BLACK = list(
      description = "Black or African American race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Tested on CL/F in the full covariate model and classified as not clinically important; not retained (Supplement Results, Figure S1).",
      source_name = "Race"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Tested on CL/F in the full covariate model and classified as not clinically important; not retained (Supplement Results, Figure S1).",
      source_name = "Race"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL/F in the full covariate model (reference 49 years) and classified as not clinically important; not retained (Supplement Results, Figure S1).",
      source_name = "Age"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "ervogastat", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ervogastat", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ervogastat", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 260L,
    n_studies = 9L,
    n_observations = 3353L,
    age_range = "18-72 years",
    age_median = "49 years",
    weight_range = "58.4-152 kg",
    weight_median = "91 kg",
    sex_female_pct = 27.3,
    race_ethnicity = c(White = 71.9, Black = 19.2, Asian = 7.3, Other = 1.5),
    disease_state = paste(
      "Healthy volunteers (33.1%), mild / moderate / severe hepatic",
      "impairment (2.3% each), NAFLD (38.8%) and presumed NASH (21.2%)"
    ),
    dose_range = paste(
      "5-1500 mg single doses; 30-600 mg every 8 hours; 25-300 mg twice",
      "daily and 300 mg once daily, alone or with clesacostat (Table S1)"
    ),
    regions = "United States; includes healthy Japanese participants (NCT04044053)",
    notes = paste(
      "Demographics from Hughes 2025 Table 1 (ervogastat popPK column).",
      "3648 concentrations measured, 3353 evaluable; the 8% below the 0.25",
      "ng/mL LLOQ were excluded (M1). Seven phase 1 studies, one phase 1b",
      "study in NAFLD (NCT03513588) and two of the three phase 2 studies",
      "(NCT03776175, NCT04399538) contributed ervogastat concentrations.",
      "Estimation by FOCE-I in NONMEM 7.5 (Code S1, run16)."
    )
  )

  ini({
    # Structural parameters (Supplement Code S1 $THETA final estimates;
    # Table S3 prints the same values rounded). Reference: 70 kg, tablet,
    # 300 mg total daily dose, healthy participant.
    lka <- log(0.704108); label("First-order absorption rate constant (1/h)") # Code S1 THETA(1) TVKA = 0.704108; Table S3 'ka' = 0.704
    lr1 <- log(101.155); label("Zero-order input rate into central (mg/h)") # Code S1 THETA(2) TVK0 = 101.155; Table S3 'k0 (mg/hour)' = 101
    lamt_fo <- log(196.508); label("Maximum amount of a dose absorbed by the first-order process (mg)") # Code S1 THETA(3) TVAK1 = 196.508; Table S3 'Ak1 (mg)' = 197
    ltlag2 <- log(1.23515); label("Absorption lag time of the zero-order process (h)") # Code S1 THETA(4) TVLAG2 = 1.23515; Table S3 'LAG2 (hour)' = 1.24
    lcl <- log(33.1973); label("Apparent clearance CL/F (L/h)") # Code S1 THETA(5) TVCL = 33.1973; Table S3 'CL/F' = 33.2
    lvc <- log(76.1031); label("Apparent central volume Vc/F (L)") # Code S1 THETA(6) TVV2 = 76.1031; Table S3 'Vc/F' = 76.1
    lq <- log(0.401027); label("Apparent intercompartmental clearance Q/F (L/h)") # Code S1 THETA(7) TVQ = 0.401027; Table S3 'Q/F' = 0.401
    lvp <- log(19.3238); label("Apparent peripheral volume Vp/F (L)") # Code S1 THETA(8) TVV3 = 19.3238; Table S3 'Vp/F' = 19.3
    lfdepot <- fixed(log(1)); label("Relative bioavailability of the tablet formulations (unitless)") # Code S1 POPF1 = POPF2 = 1*FORMFREL; structural anchor, no THETA

    # Allometric exponents, fixed (Table S3 'Fixed'; Code S1 '0.75 FIX', '1.00 FIX')
    e_wt_cl <- fixed(0.75); label("Allometric exponent of WT/70 on CL/F and Q/F (unitless)") # Code S1 THETA(12) TVALLMCL = 0.75 FIX
    e_wt_vc <- fixed(1.00); label("Allometric exponent of WT/70 on Vc/F and Vp/F (unitless)") # Code S1 THETA(13) TVALLMVD = 1.00 FIX

    # Covariate effects (Code S1: categorical effects enter as 1 + theta)
    e_dose_cl <- -0.156342; label("Power exponent of effective total daily dose / 300 mg on CL/F (unitless)") # Code S1 THETA(14) TVDOSTCL = -0.156342; Table S3 = -0.156
    e_susp_amt_fo <- 2.33717; label("Suspension fractional change in the first-order amount cap (fraction)") # Code S1 THETA(15) TVFORMAK1 = 2.33717; Table S3 = 2.34
    e_susp_fdepot <- -0.0655178; label("Suspension fractional change in relative bioavailability (fraction)") # Code S1 THETA(16) TVFORMFREL = -0.0655178; Table S3 = -0.0655
    e_hepimp_cl <- -0.268740; label("Any hepatic impairment fractional change in CL/F (fraction)") # Code S1 THETA(17) TVPTST123CL = -0.268740; Table S3 = -0.269
    e_nafld_cl <- -0.184998; label("NAFLD or presumed NASH fractional change in CL/F (fraction)") # Code S1 THETA(18) TVPTST45CL = -0.184998; Table S3 = -0.185

    # IIV (Code S1 $OMEGA variances; Table S3 prints omega as SD 0.309 / 0.316)
    etalcl ~ 0.0952774 # Code S1 $OMEGA PPVCL = 0.0952774; Table S3 omega CL/F 0.309
    etalvc ~ 0.0997356 # Code S1 $OMEGA PPVV2 = 0.0997356; Table S3 omega Vc/F 0.316

    # Residual error: W = sqrt((IPRED*PRO)^2 + ADD^2), EPS variance 1 FIX
    propSd <- 0.397557; label("Proportional residual SD, densely sampled studies (fraction)") # Code S1 THETA(9) TVRUVPROPH1 = 0.397557; Table S3 dense = 0.398
    propSd_phase2 <- 0.711411; label("Proportional residual SD, sparsely sampled phase 2 studies (fraction)") # Code S1 THETA(10) TVRUVPROPH2 = 0.711411; Table S3 sparse = 0.711
    addSd <- 0.502185; label("Additive residual SD (ng/mL)") # Code S1 THETA(11) TVRUVADD = 0.502185; Table S3 = 0.502
  })

  model({
    # Allometric scaling referenced to 70 kg (Code S1 ALLMCL / ALLMVD)
    allomCl <- (WT / 70)^e_wt_cl
    allomV <- (WT / 70)^e_wt_vc

    # Formulation effects (Code S1 FORMAK1 / FORMFREL)
    frel <- exp(lfdepot) * (1 + e_susp_fdepot * FORM_ERVOGASTAT_SUSP)

    # Clearance dose effect on the bioavailability-adjusted daily dose
    # (Code S1 EFFDOST = DOST1 * FORMFREL; skipped when DOST1 is 0)
    cl_dose <- 1
    if (DOSE_ERVOGASTAT_MGD > 0) cl_dose <- (DOSE_ERVOGASTAT_MGD * frel / 300)^e_dose_cl

    cl <- exp(lcl + etalcl) * allomCl * cl_dose *
      (1 + e_hepimp_cl * HEPIMP) * (1 + e_nafld_cl * DIS_NAFLD)
    vc <- exp(lvc + etalvc) * allomV
    q <- exp(lq) * allomCl
    vp <- exp(lvp) * allomV
    ka <- exp(lka)
    r1 <- exp(lr1)
    tlag2 <- exp(ltlag2)

    # Split of each dose: up to amt_fo mg first-order from the depot, the
    # rest zero-order into central at r1 mg/h after a lag (Supplement
    # Equation S1; Code S1 caps FK1 at 0.999 rather than 1).
    amt_fo <- exp(lamt_fo) * (1 + e_susp_amt_fo * FORM_ERVOGASTAT_SUSP)
    ffo <- amt_fo / DOSE_ERVOGASTAT_MG
    if (amt_fo > DOSE_ERVOGASTAT_MG) ffo <- 0.999

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Each administration is two dose records with the same amt: one into
    # depot and one into central with rate = -1, so the zero-order input
    # runs at r1 for frel * (1 - ffo) * amt / r1 hours (NONMEM RATE = -1
    # with modelled R2 and F2, Code S1 $INPUT notes).
    f(depot) <- frel * ffo
    f(central) <- frel * (1 - ffo)
    rate(central) <- r1
    alag(central) <- tlag2

    propSdStudy <- propSd * (1 - STUDY_PHASE2) + propSd_phase2 * STUDY_PHASE2

    # Dose in mg, volume in L: x 1000 for ng/mL (Code S1 S2 = V2/1000)
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSdStudy)
  })
}
