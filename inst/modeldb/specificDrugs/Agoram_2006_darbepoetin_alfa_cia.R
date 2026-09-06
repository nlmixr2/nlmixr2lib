Agoram_2006_darbepoetin_alfa_cia <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order subcutaneous",
    "absorption for darbepoetin alfa in adults with nonmyeloid malignancies",
    "and chemotherapy-induced anemia (Agoram 2006, AAPS J). Both IV and SC",
    "routes are supported. Body weight modifies clearance and central volume",
    "via a normalized power model (reference 70 kg), and receiving more than",
    "two cycles of platinum-containing chemotherapy during the PK assessment",
    "multiplies clearance by 0.737. Total measured serum concentration is the",
    "sum of the simulated darbepoetin alfa and an individual-specific",
    "endogenous-erythropoietin (eEPO) constant that the ELISA assay",
    "cross-detects. Exponential (log-normal) residual error. This is the",
    "PK layer of the paper's sequentially fitted PK/PD analysis; the",
    "hemoglobin response model is Agoram_2006_darbepoetin_alfa_cia_hemoglobin."
  )
  reference <- paste(
    "Agoram B, Heatherington AC, Gastonguay MR. Development and evaluation of",
    "a population pharmacokinetic-pharmacodynamic model of darbepoetin alfa in",
    "patients with nonmyeloid malignancies undergoing multicycle chemotherapy.",
    "The AAPS Journal. 2006;8(3):E552-E563. doi:10.1208/aapsj080364"
  )
  vignette <- "Agoram_2006_darbepoetin_alfa_cia"
  units <- list(time = "day", dosing = "ug", concentration = "ng/mL")

  compartmentData <- list(
    depot       = list(analyte = "darbepoetin alfa", units = "ug", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "darbepoetin alfa", units = "ug", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "darbepoetin alfa", units = "ug", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Normalized power effect on clearance and central volume with an",
        "explicit reference weight of 70 kg printed in Agoram 2006 Equations",
        "12 and 13 (`(BWT/70)^theta`). PK development cohort mean 70.8 +/-",
        "17.2 kg, range 36-123 kg (Table 2, combined 990146 + 20010162)."
      ),
      source_name        = "BWT"
    ),
    CONMED_PLATIN_GT2 = list(
      description        = "More than two cycles of platinum-containing chemotherapy during the PK assessment window",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "Agoram 2006 defines the underlying covariate PCNT as the count of",
        "platinum-containing chemotherapy cycles administered during the PK",
        "assessment, then dichotomizes it in Equation 12 as `X = 1 if PCNT >",
        "2, else X = 0` -- the split at which the paper observed the greatest",
        "difference in mean clearance. Enters as a multiplicative power term",
        "`(theta1)^X` on CL with theta1 = 0.737, i.e. a 26.3% clearance",
        "reduction in the >2-cycle group. Note that the sibling PD model",
        "Agoram_2006_darbepoetin_alfa_cia_hemoglobin dichotomizes the SAME",
        "underlying count at a DIFFERENT threshold (PCNT > 0, canonical",
        "column CONMED_PLATIN, Equation 14), so a data set supporting both",
        "models must carry both indicator columns. Mutual-consistency",
        "constraint: CONMED_PLATIN_GT2 = 1 implies CONMED_PLATIN = 1."
      ),
      source_name        = "PCNT (dichotomized at > 2)"
    )
  )

  # Screened graphically but not carried into the full covariate model
  # (Agoram 2006 Results, "Pharmacokinetic Model": age, sex, race, total
  # chemotherapy count and baseline hemoglobin were evaluated through a
  # graphical analysis and only BWT and PCNT were taken into the full model).
  # No point estimates are reported for these, so they are documentation only.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline",
      units       = "years",
      type        = "continuous",
      notes       = "Screened graphically on CL / V1 / Ka; not retained in the full PK covariate model (Agoram 2006 Results). No estimate reported."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened graphically; not retained in the full PK covariate model (Agoram 2006 Results). No estimate reported."
    ),
    HGB_BL = list(
      description = "Baseline hemoglobin",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Screened graphically as a PK covariate; not retained in the full PK covariate model (Agoram 2006 Results). No estimate reported. Baseline hemoglobin IS an estimated parameter (Hb0) in the sibling PD model rather than a covariate."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 140L,
    n_studies      = 2L,
    n_observations = 1442L,
    age_range      = "22-87 years (combined PK data set; Agoram 2006 Table 2)",
    age_median     = "Mean 62.3 years, SD 12.9 (median not tabulated)",
    weight_range   = "36-123 kg (combined PK data set; Agoram 2006 Table 2)",
    weight_median  = "Mean 70.8 kg, SD 17.2 (median not tabulated)",
    sex_female_pct = 69,
    race_ethnicity = "Screened as a covariate but not tabulated in Agoram 2006.",
    disease_state  = paste(
      "Adults (>= 18 years) with nonmyeloid malignancies receiving cyclic",
      "chemotherapy, with chemotherapy-induced anemia. Study 990146 required",
      "hemoglobin <= 13.0 g/dL (so some subjects were not anemic by the",
      "<= 11.0 g/dL criterion); study 20010162 required hemoglobin >= 9.0 and",
      "<= 11.0 g/dL. ECOG performance status 0-2 with adequate renal and",
      "hepatic function. Excluded: history of seizures, significant cardiac or",
      "inflammatory disease, primary hematologic disorder causing anemia, prior",
      "rHuEPO, > 2 RBC transfusions within 4 weeks of study drug assignment, or",
      "any RBC transfusion during the current chemotherapy cycle before",
      "randomization. Baseline hemoglobin 10.30 +/- 1.19 g/dL, range 5.5-12.7",
      "(Agoram 2006 Table 2). Tumor types: breast 27%, gynecologic 19%,",
      "gastrointestinal 14%, lymphoma 14%, lung 11%, genitourinary 9%,",
      "other 8%."
    ),
    dose_range     = paste(
      "Study 990146: 2.25 ug/kg QW, IV or SC. Study 20010162: 6.75 ug/kg SC",
      "Q3W, randomized 1:1 to day 1 or day 15 of a 21-day chemotherapy cycle",
      "(Agoram 2006 Table 1)."
    ),
    regions        = "Two Amgen-sponsored clinical studies (990146, 20010162); geographic sites not stated in the publication.",
    notes          = paste(
      "PK model development data set: 140 patients, 1442 serum darbepoetin",
      "alfa concentrations, pooled from Amgen studies 990146 (n = 56) and",
      "20010162 (n = 84). Full PK profiles were collected at weeks 1 and 9",
      "(990146; 5 and 30 min and 6, 24, 48, 72, 96, 120 and 168 h postdose)",
      "and on day 1 (20010162; 6, 24, 48, 72, 96, 120, 168, 174, 192, 216,",
      "312, 336, 360 and 480 h postdose), plus weekly trough samples.",
      "Darbepoetin alfa was quantified with the Quantikine IVD human",
      "erythropoietin ELISA (R&D Systems, Minneapolis MN); standard curve",
      "0.125-5.0 ng/mL, LLOQ 0.14 ng/mL, interassay CV 5-7%. Endogenous",
      "erythropoietin cross-reacts with this assay, so the reported",
      "concentrations include an eEPO contribution which the model carries as",
      "a subject-constant additive term. Estimation used NONMEM V level 1.2",
      "(FOCE) with 1000-replicate nonparametric bootstrap. Demographics from",
      "Table 2, trial summary from Table 1."
    )
  )

  ini({
    # Structural PK parameters (Agoram 2006 Table 3, "Full Pharmacokinetic
    # Model Parameter Estimates"). Reported in mL and mL/day; converted to L
    # and L/day so that (ug in central) / (L) gives ug/L = ng/mL directly.
    lcl <- log(2.010);   label("Typical clearance at WT 70 kg and PCNT <= 2 (L/day)")  # Agoram 2006 Table 3: theta_CL = 2010 mL/day, SE 128 (6.37% RSE)
    lvc <- log(3.390);   label("Typical central volume at WT 70 kg (L)")               # Agoram 2006 Table 3: theta_V1 = 3390 mL, SE 313 (9.23% RSE)
    lvp <- log(0.251);   label("Peripheral volume (L)")                                # Agoram 2006 Table 3: theta_V2 = 251 mL, SE 92.1 (36.7% RSE)
    lka <- log(0.318);   label("First-order SC absorption rate constant (1/day)")      # Agoram 2006 Table 3: theta_Ka = 0.318 1/day, SE 0.0275 (8.65% RSE)

    # Intercompartmental clearance is the only row of Table 3 with a blank SE
    # AND a blank bootstrap column; every estimated parameter in that table has
    # both. Encoded as fixed(); see the vignette Assumptions and deviations.
    lq <- fixed(log(2.900)); label("Intercompartmental clearance (L/day)")             # Agoram 2006 Table 3: theta_Q = 2900 mL/day, no SE and no bootstrap CI reported

    lfdepot <- log(0.443); label("Subcutaneous bioavailability (fraction)")            # Agoram 2006 Table 3: F = 0.443 (18.3% RSE); IIV of F not estimated (IV and SC data came from separate patients)
    leepo   <- log(0.415); label("Endogenous-EPO contribution to measured serum concentration (ng/mL)")  # Agoram 2006 Table 3: theta_C0 = 0.415 ng/mL, SE 0.0391 (9.42% RSE)

    # Covariate effects. Normalized power model on body weight (reference
    # 70 kg) and a multiplicative power term on the platinum-chemotherapy
    # indicator (Agoram 2006 Equations 12 and 13).
    e_conmed_platin_gt2_cl <- 0.737;  label("Multiplicative effect of > 2 concomitant platinum chemotherapy cycles on CL (unitless)")  # Agoram 2006 Table 3: theta_1 (PCNT on CL) = 0.737, asymptotic 95% CI (0.551, 0.923); NULL value 1
    e_wt_cl                <- 0.623;  label("Power exponent of (WT/70 kg) on clearance (unitless)")                                   # Agoram 2006 Table 3: theta_2 (BWT on CL) = 0.623, asymptotic 95% CI (0.172, 1.07); NULL value 0
    e_wt_vc                <- 0.639;  label("Power exponent of (WT/70 kg) on central volume (unitless)")                              # Agoram 2006 Table 3: theta_3 (BWT on V1) = 0.639, asymptotic 95% CI (-0.01, 1.29); NULL value 0

    # IIV - log-normal (Agoram 2006 Equation 1, P_i = TVP * exp(eta)). Table 3
    # reports omega^2 VARIANCES directly. The Omega block carries one
    # off-diagonal, between CL and V1. Q, V2 and F have no IIV: the paper
    # states it was not possible to estimate the IIV of Q and V2, and that IIV
    # of F was not estimated because IV and SC data came from separate patients.
    etalcl + etalvc ~ c(0.181,
                        0.134, 0.225)  # Agoram 2006 Table 3: omega^2_CL = 0.181 (21.9% RSE), omega^2_CL,V1 = 0.134 (60.4% RSE), omega^2_V1 = 0.225 (40.4% RSE)
    etalka   ~ 0.0883                  # Agoram 2006 Table 3: omega^2_Ka = 0.0883 (30.8% RSE)
    etaleepo ~ 0.501                   # Agoram 2006 Table 3: omega^2_C0 = 0.501 (23.4% RSE)

    # Residual error - log-transformed exponential model (Agoram 2006 Equation
    # 2: ln(C_ij) = ln(Chat_ij + Chat_i,eEPO) + eps_ij), i.e. log-normal on the
    # total measured concentration. sigma1^2 = 0.483 -> SD = sqrt(0.483).
    expSd <- sqrt(0.483); label("Log-scale residual SD (~CV 74%)")  # Agoram 2006 Table 3: sigma1^2 = 0.483, SE 0.0348 (7.20% RSE)
  })

  model({
    # Individual PK parameters. Equation 12 applies the platinum indicator and
    # the weight power term to CL; Equation 13 applies the weight power term
    # to V1. Q and V2 carry no covariate effects.
    cl   <- exp(lcl + etalcl) * (e_conmed_platin_gt2_cl^CONMED_PLATIN_GT2) * (WT / 70)^e_wt_cl
    vc   <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vp   <- exp(lvp)
    q    <- exp(lq)
    ka   <- exp(lka + etalka)
    eepo <- exp(leepo + etaleepo)

    # Micro-constants for the two-compartment disposition (Agoram 2006 Figure
    # 1A: K = CL/V1, Q/V1 into and Q/V2 out of the peripheral compartment).
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # IV doses enter central directly (data set cmt = central, f(central) = 1);
    # SC doses enter depot (data set cmt = depot) with bioavailability F.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                               k12 * central - k21 * peripheral1

    f(depot) <- exp(lfdepot)

    # Total measured serum concentration = simulated darbepoetin alfa + the
    # subject's endogenous EPO, which cross-reacts with the ELISA (Agoram 2006
    # Analytical Methods and Equation 2). Dose in ug, central in ug, vc in L
    # -> central/vc in ug/L = ng/mL; eepo also in ng/mL.
    Cc <- central / vc + eepo
    Cc ~ lnorm(expSd)
  })
}
