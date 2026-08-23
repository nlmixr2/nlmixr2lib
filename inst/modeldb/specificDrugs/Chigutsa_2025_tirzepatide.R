Chigutsa_2025_tirzepatide <- function() {
  description <- paste(
    "Tirzepatide exposure-response body composition model (Chigutsa 2025,",
    "SURMOUNT-1). Two parallel indirect-response (turnover) compartments",
    "carry fat-free mass (FFM) and fat mass as separate dependent",
    "variables; model-predicted total body weight is their sum. Tirzepatide",
    "inhibits the zero-order formation (Kin) of each pool through an Imax",
    "function of plasma concentration, with a separate Imax and IC50 for",
    "FFM and for fat mass -- the estimated three-fold greater effect on fat",
    "mass is what produces the improvement in body composition. An additive",
    "placebo effect on the same fractional-Kin-reduction scale decays",
    "exponentially with a 40.3-week half-life. The upstream two-compartment",
    "tirzepatide PK model (Schneck 2024) is reproduced inline with all",
    "parameters fixed, and is driven by a semimechanistic body-composition",
    "allometry in which the FFM and fat-mass ODE states themselves supply",
    "the size descriptor -- so weight reduction feeds back on clearance and",
    "volume over time. Sex acts on baseline FFM, baseline fat mass, Imax",
    "and IC50; Asian race acts on both baselines."
  )
  reference <- paste(
    "Chigutsa E, Her L, Ma X, Urva S, Schneck K. A pharmacometric method",
    "for quantitative determination of improvement in body composition and",
    "characterization of the exposure-response relationship during",
    "treatment of obesity with tirzepatide.",
    "Clin Pharmacol Ther. 2025;118(6):1489-1497. doi:10.1002/cpt.3750.",
    "PMCID PMC12641085. Model structure and the fixed random-effect",
    "component for placebo waning are taken from the NONMEM control stream",
    "in Supplementary Material S1 (file CPT-118-1489-s001.txt).",
    "PK layer reproduced from the upstream population PK model:",
    "Schneck K, Urva S. Population pharmacokinetics of the GIP/GLP receptor",
    "agonist tirzepatide. CPT Pharmacometrics Syst Pharmacol.",
    "2024;13:494-503. doi:10.1002/psp4.13099. PMCID PMC10962491.",
    "The fat-free mass / fat mass dependent variables were calculated from",
    "total body weight, height and sex using:",
    "Janmahasatian S, Duffull SB, Ash S, Ward LC, Byrne NM, Green B.",
    "Quantification of lean bodyweight. Clin Pharmacokinet.",
    "2005;44(10):1051-1065. doi:10.2165/00003088-200544100-00004."
  )
  vignette <- "Chigutsa_2025_tirzepatide"

  # `fat` and `ffm` are the paper's two body-composition turnover states.
  # They are genuinely paper-mechanistic (body compartments, not drug
  # compartments) and have no canonical equivalent; `fat` follows the
  # precedent set by Bosch_2024_glp1ra_bodyweight.R.
  paper_specific_compartments <- c("fat", "ffm")

  units <- list(
    time          = "week",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  compartmentData <- list(
    depot = list(
      analyte = "tirzepatide", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "tirzepatide", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "tirzepatide", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    fat = list(
      analyte = "fat mass", units = "kg",
      specimen = "not applicable", verified = TRUE
    ),
    ffm = list(
      analyte = "fat-free mass", units = "kg",
      specimen = "not applicable", verified = TRUE
    )
  )

  covariateData <- list(
    SEXF = list(
      description = paste(
        "Biological sex indicator, 1 = female, 0 = male. Acts on baseline",
        "fat-free mass, baseline fat mass, Imax and IC50 for both",
        "endpoints; all six effects are fractional, applied as",
        "(1 + theta * SEXF)."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Supplementary Material S1 sets FEMALE = 1 when the source column",
        "SEX == 1, so the reference category is male and the Table 1",
        "covariate rows read 'in females relative to males'. The",
        "encoding is confirmed arithmetically by the Table 1 footnote b",
        "baseline total body weights: male 73.5 + 45.0 = 118.5 kg and",
        "female 73.5 * (1 - 0.312) + 45.0 * (1 + 0.0539) = 98.0 kg,",
        "reproducing both published values exactly."
      ),
      source_name = "SEX (recoded to FEMALE in Supplementary Material S1)"
    ),
    RACE_ASIAN = list(
      description = paste(
        "Asian race indicator, 1 = Asian, 0 = non-Asian. Acts on baseline",
        "fat-free mass and baseline fat mass only; no Asian effect on drug",
        "response was retained."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian)",
      notes              = paste(
        "Supplementary Material S1 sets ASIAN = 1 when the source column",
        "RACE == 58. Asian participants were estimated to have 11% lower",
        "baseline fat-free mass and 25% lower baseline fat mass",
        "(Chigutsa 2025 Results, 'Model results and diagnostics')."
      ),
      source_name = "RACE (recoded to ASIAN in Supplementary Material S1)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 2539L,
    n_studies      = 1L,
    age_range      = NA_character_,
    weight_median  = "107 kg (mean baseline weight in the Table 2 simulation)",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = paste(
      "Adults with obesity or overweight WITHOUT type 2 diabetes mellitus:",
      "BMI >= 30 kg/m2, or BMI >= 27 kg/m2 with at least one",
      "weight-related complication."
    ),
    dose_range     = paste(
      "Placebo, 5, 10 or 15 mg tirzepatide subcutaneously once weekly for",
      "72 weeks. Dose into the `depot` compartment in mg."
    ),
    regions        = NA_character_,
    notes          = paste(
      "SURMOUNT-1 (NCT04184622), a phase 3 randomised (1:1:1:1) trial with",
      "monthly body weight measurements over 72 weeks. The model was",
      "developed on a random 50% of the data and externally evaluated on",
      "the held-out 50%; it was further evaluated against DXA measurements",
      "in the ~10% sub-population who had them, which were not used in",
      "model development. Fat-free mass and fat mass were not measured but",
      "CALCULATED for every participant from total body weight, height and",
      "sex using the Janmahasatian 2005 lean-bodyweight equations, then",
      "fitted simultaneously as two dependent variables. Baseline",
      "demographics are not tabulated in the paper; the only published",
      "anchors are the Table 1 footnote b typical baseline weights",
      "(male 118.5 kg, female 98.0 kg) and the Table 2 footnote a mean",
      "simulation baseline weight of 107 kg."
    )
  )

  ini({
    # ================================================================
    # PK layer -- tirzepatide two-compartment model with first-order
    # absorption. Every value is FIXED: this analysis used individual
    # post hoc PK parameters from the upstream population PK model and
    # estimated no PK parameter of its own (Chigutsa 2025 Methods,
    # 'Exposure-response modeling'). Rate parameters are converted from
    # /h to /week by x168, exactly as Supplementary Material S1 does
    # ("TIME is in weeks, so convert /h to /week").
    # ================================================================
    lka <- fixed(log(6.2664))
    label("Log of first-order absorption rate constant (1/week)")               # Schneck 2024 Table 3: ka 0.0373 /h * 168
    lcl <- fixed(log(5.5272))
    label("Log of clearance at the reference body size (L/week per 70 kg)")     # Schneck 2024 Table 3: CL 0.0329 L/h/70 kg * 168
    lq <- fixed(log(21.168))
    label("Log of intercompartmental clearance at reference size (L/week per 70 kg)") # Schneck 2024 Table 3: Q 0.126 L/h/70 kg * 168
    lvc <- fixed(log(2.47))
    label("Log of central volume of distribution at reference size (L per 70 kg)")    # Schneck 2024 Table 3
    lvp <- fixed(log(3.98))
    label("Log of peripheral volume of distribution at reference size (L per 70 kg)") # Schneck 2024 Table 3
    lfdepot <- fixed(log(0.8))
    label("Log of subcutaneous bioavailability (fraction)")                     # Schneck 2024 Table 3: bioavailability 0.8 fixed

    e_wt_cl <- fixed(0.8)
    label("Allometric exponent on the body-size descriptor for CL and Q (unitless)")   # Schneck 2024 Table 3 (0.8 fixed); Suppl S1 CLWT = AVLC**0.8
    e_wt_vc <- fixed(1)
    label("Allometric exponent on the body-size descriptor for Vc and Vp (unitless)")  # Schneck 2024 Table 3 (1 fixed); Suppl S1 VDWT = AVLV
    # The fat-mass fractions come from the upstream PK model, NOT from the
    # stale constants hard-coded in Chigutsa 2025 Supplementary Material S1
    # (FATCL = 0.711, FATVD = 0.417). Three independent lines of evidence
    # give 1 and 0.482 instead; see the vignette Errata for the full audit.
    e_fatfrac_cl <- fixed(1)
    label("Fraction of fat mass entering the body-size descriptor for CL and Q (unitless)") # Schneck 2024 Table 3 reports no fat fraction on CL (fixed to 1) and states the estimate "approached a value of 1"; Chigutsa 2025 Methods: "CL best correlated with total body weight (fat-free mass plus fat mass)"
    e_fatfrac_vc <- fixed(0.482)
    label("Fraction of fat mass entering the body-size descriptor for Vc and Vp (unitless)") # Schneck 2024 Table 3: "Fraction of fat mass with effect on volume of distribution 0.482 (0.447, 0.524)"; Chigutsa 2025 Methods rounds it to "48%"

    # ================================================================
    # PD layer -- parallel indirect-response models for FFM and fat mass.
    # Table 1 reports the FINAL estimates; the $THETA block of
    # Supplementary Material S1 holds INITIAL estimates only (e.g. Imax
    # FFM initial 0.2733 vs final 0.144) and is NOT used for values.
    # ================================================================
    lrbase_ffm <- log(73.5)
    label("Log of baseline fat-free mass in a non-Asian male (kg)")             # Chigutsa 2025 Table 1: 73.5 kg (0.592% RSE)
    lrbase_fat <- log(45.0)
    label("Log of baseline fat mass in a non-Asian male (kg)")                  # Chigutsa 2025 Table 1: 45.0 kg (1.54% RSE)
    lkout <- log(0.0314)
    label("Log of the first-order loss rate constant Kout, shared by both pools (1/week)") # Chigutsa 2025 Table 1: 0.0314 /week (ln2/Kout = 22.1 weeks, the published ~22-week weight-reduction half-life)

    limax_ffm <- log(0.144)
    label("Log of maximum fractional inhibition of fat-free mass Kin in males (unitless)") # Chigutsa 2025 Table 1: 0.144 (11.2% RSE)
    limax_fat <- log(0.319)
    label("Log of maximum fractional inhibition of fat mass Kin in males (unitless)")      # Chigutsa 2025 Table 1: 0.319 (9.81% RSE)
    lic50_ffm <- log(1760)
    label("Log of tirzepatide concentration at half-maximal inhibition of fat-free mass Kin in males (ng/mL)") # Chigutsa 2025 Table 1: 1760 ng/mL (9.89% RSE)
    lic50_fat <- log(518)
    label("Log of tirzepatide concentration at half-maximal inhibition of fat mass Kin in males (ng/mL)")      # Chigutsa 2025 Table 1: 518 ng/mL (17.9% RSE)

    # Placebo effects carry an ADDITIVE (not log-normal) random effect in
    # Supplementary Material S1 -- PLAC = MU_6 + ETA(6) -- so they are
    # parameterised on the natural scale rather than log-transformed.
    plac_ffm <- 0.0658
    label("Placebo fractional reduction in fat-free mass Kin at time zero (unitless)")     # Chigutsa 2025 Table 1: 0.0658 (3.63% RSE)
    plac_fat <- 0.213
    label("Placebo fractional reduction in fat mass Kin at time zero (unitless)")          # Chigutsa 2025 Table 1: 0.213 (3.60% RSE)
    lthalfwane <- log(40.3)
    label("Log of the half-life of the waning placebo effect (weeks)")                     # Chigutsa 2025 Table 1: 40.3 weeks (4.32% RSE)

    # Covariate effects -- all fractional, applied as (1 + theta * IND),
    # matching TVE0 = THETA(1)*(1+THETA(8)*FEMALE)*(1+THETA(15)*ASIAN)
    # and the analogous lines in Supplementary Material S1.
    e_sexf_rbase_ffm <- -0.312
    label("Fractional change in baseline fat-free mass in females relative to males")      # Chigutsa 2025 Table 1: -0.312 (1.55% RSE)
    e_sexf_rbase_fat <- 0.0539
    label("Fractional change in baseline fat mass in females relative to males")           # Chigutsa 2025 Table 1: 0.0539 (35.8% RSE)
    e_sexf_imax_ffm <- 1.78
    label("Fractional change in Imax for fat-free mass in females relative to males")      # Chigutsa 2025 Table 1: 1.78 (18.88% RSE)
    e_sexf_imax_fat <- 0.809
    label("Fractional change in Imax for fat mass in females relative to males")           # Chigutsa 2025 Table 1: 0.809 (24.1% RSE)
    e_sexf_ic50_ffm <- 1.14
    label("Fractional change in IC50 for fat-free mass in females relative to males")      # Chigutsa 2025 Table 1: 1.14 (18.9% RSE)
    e_sexf_ic50_fat <- 2.05
    label("Fractional change in IC50 for fat mass in females relative to males")           # Chigutsa 2025 Table 1: 2.05 (24.1% RSE)
    e_asian_rbase_ffm <- -0.110
    label("Fractional change in baseline fat-free mass in Asians relative to non-Asians")  # Chigutsa 2025 Table 1: -0.110 (8.37% RSE)
    e_asian_rbase_fat <- -0.254
    label("Fractional change in baseline fat mass in Asians relative to non-Asians")       # Chigutsa 2025 Table 1: -0.254 (7.91% RSE)

    # ================================================================
    # Interindividual variability. Table 1 reports CV% and correlations;
    # variances below are back-calculated. For the log-normal random
    # effects omega^2 = log(1 + CV^2). For the two placebo terms the
    # random effect is ADDITIVE on the natural-scale parameter, so the
    # reported CV% is SD / typical value and omega^2 = (CV * theta)^2.
    # Each block reproduces the corresponding $OMEGA BLOCK of
    # Supplementary Material S1 to within initial-vs-final movement.
    # ================================================================
    etalrbase_ffm + etalrbase_fat + etalkout ~
      c(0.0120274,
        0.0267842, 0.0823622,
        -0.0148163, -0.0476341, 0.9312188)
    # Chigutsa 2025 Table 1 IIV block: CV 11.0%, 29.3%, 124%; correlations
    # fat~FFM 0.851, FFM~Kout -0.140, fat~Kout -0.172.
    # Suppl S1 $OMEGA BLOCK(3) initials: 0.0121068 / 0.0268708 0.0825637 /
    # -0.0150427 -0.0463987 0.915321.

    etalimax_ffm + etalimax_fat ~ c(0.5388444, 0.4581327, 0.3950214)
    # Chigutsa 2025 Table 1: CV 84.5% and 69.6%, correlation 0.993.
    # Suppl S1 $OMEGA BLOCK(2) initials: 0.488297 / 0.376389 0.301214.

    etaplac_ffm + etaplac_fat ~ c(0.00404023, 0.01252643, 0.04068733)
    # Chigutsa 2025 Table 1: CV 96.6% and 94.7%, correlation 0.977, on
    # ADDITIVE etas: SD = 0.966*0.0658 = 0.06356 and 0.947*0.213 = 0.20171.
    # Suppl S1 $OMEGA BLOCK(2) initials: 0.00386354 / 0.0122032 0.0402101
    # confirm the SD/theta reading of the reported CV%.

    etalic50_ffm + etalic50_fat ~ c(0.0639593, 0.1137870, 0.2026763)
    # Chigutsa 2025 Table 1: CV 25.7% and 47.4%, correlation 0.9994.
    # Suppl S1 $OMEGA BLOCK(2) initials: 0.0798053 / 0.143795 0.259453.

    etalthalfwane ~ fixed(0.0225)
    # Suppl S1 `$OMEGA 0.0225 FIX ; PPV_WANE` (CV 15.1%). Held fixed, which
    # is why Table 1 carries no IIV row for the placebo waning half-life.

    # PK interindividual variability, carried unchanged from the upstream
    # population PK model and held fixed here (this analysis used the
    # individual post hoc PK parameters rather than re-estimating them).
    etalka ~ fixed(0.0493852)
    label("IIV on absorption rate constant")                                   # Schneck 2024 Table 3: ka IIV 22.5% -> log(1 + 0.225^2)
    etalcl ~ fixed(0.0199634)
    label("IIV on clearance")                                                  # Schneck 2024 Table 3: CL IIV 14.2% -> log(1 + 0.142^2)
    etalvc ~ fixed(0.2151920)
    label("IIV on central volume of distribution")                             # Schneck 2024 Table 3: Vc IIV 49.0% -> log(1 + 0.490^2)

    # ================================================================
    # Residual error. Supplementary Material S1 applies ERR(1) to fat
    # mass (CMT 6) and ERR(2) to fat-free mass (CMT 7), both purely
    # proportional. The $SIGMA BLOCK(2) of Supplementary Material S1 is
    # the FINAL estimate: sqrt(0.00112075) = 3.35% and
    # sqrt(9.73119e-5) = 0.985% reproduce Table 1 exactly, and
    # 0.000325097 / (0.03348 * 0.009865) = 98.4% reproduces the reported
    # residual correlation. That correlation is a NONMEM L2 construct and
    # is not representable in nlmixr2; see the vignette Errata.
    # ================================================================
    propSd_fatMass <- 0.0335
    label("Proportional residual error for fat mass (fraction)")               # Chigutsa 2025 Table 1: 3.35%; Suppl S1 $SIGMA sqrt(0.00112075)
    propSd_fatFreeMass <- 0.00985
    label("Proportional residual error for fat-free mass (fraction)")          # Chigutsa 2025 Table 1: 0.985%; Suppl S1 $SIGMA sqrt(9.73119e-5)
  })

  model({
    # ---- 1. Covariate-adjusted typical values -------------------------
    # Supplementary Material S1:
    #   TVE0    = THETA(1)*(1+THETA(8)*FEMALE)*(1+THETA(15)*ASIAN)
    #   TVFE0   = THETA(2)*(1+THETA(9)*FEMALE)*(1+THETA(16)*ASIAN)
    #   TVIMAX  = THETA(4)*(1+THETA(17)*FEMALE)
    #   TVFIMAX = THETA(5)*(1+THETA(18)*FEMALE)
    #   TVIC50  = THETA(10)*(1+THETA(13)*FEMALE)
    #   TVFIC50 = THETA(11)*(1+THETA(14)*FEMALE)
    # with the log-normal random effect applied after the covariates.
    rbase_ffm <- exp(lrbase_ffm + etalrbase_ffm) *
      (1 + e_sexf_rbase_ffm * SEXF) * (1 + e_asian_rbase_ffm * RACE_ASIAN)
    rbase_fat <- exp(lrbase_fat + etalrbase_fat) *
      (1 + e_sexf_rbase_fat * SEXF) * (1 + e_asian_rbase_fat * RACE_ASIAN)

    imax_ffm <- exp(limax_ffm + etalimax_ffm) * (1 + e_sexf_imax_ffm * SEXF)
    imax_fat <- exp(limax_fat + etalimax_fat) * (1 + e_sexf_imax_fat * SEXF)
    ic50_ffm <- exp(lic50_ffm + etalic50_ffm) * (1 + e_sexf_ic50_ffm * SEXF)
    ic50_fat <- exp(lic50_fat + etalic50_fat) * (1 + e_sexf_ic50_fat * SEXF)

    kout <- exp(lkout + etalkout)

    # Placebo effects use additive random effects (Suppl S1 PLAC = MU_6 +
    # ETA(6)), so no exponentiation here.
    placeff_ffm <- plac_ffm + etaplac_ffm
    placeff_fat <- plac_fat + etaplac_fat
    wane <- log(2) / exp(lthalfwane + etalthalfwane)

    # ---- 2. Body-size descriptor and individual PK parameters ---------
    # Suppl S1:
    #   AVLC = (FFMV + FATV*FATCL)/70 ; CLWT = AVLC**0.8
    #   AVLV = (FFMV + FATV*FATVD)/70 ; VDWT = AVLV
    # FFMV / FATV were the time-varying body-composition columns; here the
    # model's own ODE states supply them, so weight reduction feeds back
    # on clearance and volume as treatment proceeds.
    avlc <- (ffm + fat * e_fatfrac_cl) / 70
    avlv <- (ffm + fat * e_fatfrac_vc) / 70
    clwt <- avlc^e_wt_cl
    vdwt <- avlv^e_wt_vc

    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * clwt
    q  <- exp(lq) * clwt
    vc <- exp(lvc + etalvc) * vdwt
    vp <- exp(lvp) * vdwt

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- 3. Drug effect ----------------------------------------------
    # Concentration in ng/mL: amounts are mg and volumes L, so mg/L * 1000.
    cp <- 1000 * central / vc

    # Suppl S1 $DES:
    #   EFF  = PLAC*(EXP(-WANE*T))  + IMAX*CP/(CP+IC50)  ; IF(EFF.GT.1) EFF=1
    #   FEFF = FPLAC*(EXP(-WANE*T)) + FIMAX*CP/(CP+FIC50); IF(FEFF.GT.1) FEFF=1
    # The cap is load-bearing: with 84.5% / 69.6% CV on Imax and the female
    # multipliers, individual Imax draws can exceed 1.
    eff_ffm <- min(1.0, placeff_ffm * exp(-wane * t) +
                     imax_ffm * cp / (cp + ic50_ffm))
    eff_fat <- min(1.0, placeff_fat * exp(-wane * t) +
                     imax_fat * cp / (cp + ic50_fat))

    # ---- 4. ODE system ------------------------------------------------
    # Kin is pinned to the baseline so that each pool starts at steady
    # state: KIN = KOUT*E0 and A_0(7) = E0 (Suppl S1 $PK).
    kin_ffm <- kout * rbase_ffm
    kin_fat <- kout * rbase_fat

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(fat)         <-  kin_fat * (1 - eff_fat) - kout * fat
    d/dt(ffm)         <-  kin_ffm * (1 - eff_ffm) - kout * ffm

    fat(0) <- rbase_fat
    ffm(0) <- rbase_ffm

    f(depot) <- exp(lfdepot)

    # ---- 5. Observations ----------------------------------------------
    # Cc is the exposure driving the PD and is not itself an observed
    # endpoint in this analysis, so it carries no residual error.
    Cc          <- cp
    fatMass     <- fat
    fatFreeMass <- ffm
    bodyWeight  <- fat + ffm

    fatMass     ~ prop(propSd_fatMass)
    fatFreeMass ~ prop(propSd_fatFreeMass)
  })
}
