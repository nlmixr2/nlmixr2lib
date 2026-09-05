Agoram_2006_darbepoetin_alfa_cia_hemoglobin <- function() {
  description <- paste(
    "Sequentially fitted PK/PD model of the hemoglobin response to",
    "darbepoetin alfa in adults with nonmyeloid malignancies and",
    "chemotherapy-induced anemia (Agoram 2006, AAPS J). The two-compartment",
    "PK layer is carried at the fixed-effect estimates of the full covariate",
    "PK model (Agoram_2006_darbepoetin_alfa_cia) with its interindividual and",
    "residual variability deliberately dropped, exactly as the paper fitted",
    "the PD stage. Serum darbepoetin alfa concentration stimulates a",
    "zero-order hemoglobin production rate through an Emax-type term",
    "(incremental Smax = 43.7%, S50 = 3.68 ng/mL), feeding a Friberg-style",
    "catenary chain: one progenitor transit compartment of mean transit time",
    "RBCPT = 4.68 days followed by four equal-transit-time mature red-cell",
    "lifespan compartments summing to RBCLS = 82.2 days. Total hemoglobin is",
    "the sum of the four lifespan compartments. Receiving any concomitant",
    "platinum-containing chemotherapy multiplies S50 by 1.86. Additive",
    "residual error on hemoglobin."
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
    depot       = list(analyte = "darbepoetin alfa", units = "ug",   specimen = "administration site", verified = TRUE),
    central     = list(analyte = "darbepoetin alfa", units = "ug",   specimen = "serum",               verified = TRUE),
    peripheral1 = list(analyte = "darbepoetin alfa", units = "ug",   specimen = "serum",               verified = TRUE),
    precursor1  = list(analyte = "hemoglobin",       units = "g/dL", specimen = "not applicable",      verified = TRUE),
    hb1         = list(analyte = "hemoglobin",       units = "g/dL", specimen = "whole blood",         verified = TRUE),
    hb2         = list(analyte = "hemoglobin",       units = "g/dL", specimen = "whole blood",         verified = TRUE),
    hb3         = list(analyte = "hemoglobin",       units = "g/dL", specimen = "whole blood",         verified = TRUE),
    hb4         = list(analyte = "hemoglobin",       units = "g/dL", specimen = "whole blood",         verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters only through the fixed PK layer, via the normalized power",
        "effects on clearance and central volume with an explicit reference",
        "weight of 70 kg printed in Agoram 2006 Equations 12 and 13",
        "(`(BWT/70)^theta`). Body weight was screened as a direct PD covariate",
        "on RBCLS and S50 and was NOT retained (Agoram 2006 Results, PkPd",
        "Model). PD development cohort mean 69.9 +/- 16.2 kg, range 39-136 kg",
        "(Table 2, combined 20010162 + 980290 + 980291)."
      ),
      source_name        = "BWT"
    ),
    CONMED_PLATIN_GT2 = list(
      description        = "More than two cycles of concomitant platinum-containing chemotherapy during the PK assessment window",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "Enters only through the fixed PK layer (Agoram 2006 Equation 12,",
        "`X = 1 if PCNT > 2, else X = 0`), multiplying clearance by",
        "theta1 = 0.737. This is the SAME underlying platinum-cycle count",
        "PCNT that drives the PD covariate CONMED_PLATIN below, but",
        "dichotomized at a DIFFERENT threshold, so a data set supporting this",
        "model must carry both indicator columns. Mutual-consistency",
        "constraint: CONMED_PLATIN_GT2 = 1 implies CONMED_PLATIN = 1."
      ),
      source_name        = "PCNT (dichotomized at > 2)"
    ),
    CONMED_PLATIN = list(
      description        = "Any concomitant platinum-containing chemotherapy during the assessment window",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "Agoram 2006 Equation 14 dichotomizes the platinum-cycle count PCNT",
        "as `X = 1 if PCNT > 0, else X = 0` and applies it to S50 as a",
        "multiplicative power term `(theta4)^X` with theta4 = 1.86, i.e. an",
        "86% higher half-maximal stimulating concentration (a blunted",
        "hemoglobin response) in patients receiving platinum-containing",
        "chemotherapy. The paper notes that the asymptotic 95% CI of theta4,",
        "(0.982, 2.74), overlaps its NULL value of 1, indicating lack of",
        "information about this effect rather than absence of an effect; it",
        "is retained here because the publication reports it as part of the",
        "full PkPd model (Table 4)."
      ),
      source_name        = "PCNT (dichotomized at > 0)"
    )
  )

  # Screened as PD covariates against the empirical Bayes estimates of the
  # RBCLS and S50 random effects from the base PkPd model and NOT retained
  # (Agoram 2006 Results, "PkPd Model": age, BWT, sex, race, and serum levels
  # of lactate dehydrogenase, ferritin and albumin "did not influence
  # parameter variability"). No point estimates are reported for any of them,
  # so these entries are documentation only. Race is documented in the
  # population notes rather than here because the paper never states which
  # race groups were present.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline",
      units       = "years",
      type        = "continuous",
      notes       = "Screened graphically against the RBCLS and S50 etas; not retained in the full PkPd covariate model (Agoram 2006 Results). No estimate reported."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened graphically against the RBCLS and S50 etas; not retained (Agoram 2006 Results). No estimate reported."
    ),
    LDH = list(
      description = "Serum lactate dehydrogenase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened graphically against the RBCLS and S50 etas; not retained (Agoram 2006 Results). No estimate reported."
    ),
    FERRITIN = list(
      description = "Serum ferritin",
      units       = "ug/L",
      type        = "continuous",
      notes       = "Screened graphically against the RBCLS and S50 etas; not retained (Agoram 2006 Results). No estimate reported."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened graphically against the RBCLS and S50 etas; not retained (Agoram 2006 Results). No estimate reported."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 573L,
    n_studies      = 3L,
    n_observations = 6356L,
    age_range      = "20-91 years (combined PkPd data set; Agoram 2006 Table 2)",
    age_median     = "Mean 60.1 years, SD 12.3 (median not tabulated)",
    weight_range   = "39-136 kg (combined PkPd data set; Agoram 2006 Table 2)",
    weight_median  = "Mean 69.9 kg, SD 16.2 (median not tabulated)",
    sex_female_pct = 70,
    race_ethnicity = "Screened as a PD covariate but not tabulated in Agoram 2006.",
    disease_state  = paste(
      "Adults (>= 18 years) with nonmyeloid malignancies receiving cyclic",
      "chemotherapy, with chemotherapy-induced anemia (hemoglobin >= 9.0 and",
      "<= 11.0 g/dL at entry). ECOG performance status 0-2 with adequate renal",
      "and hepatic function. Excluded: history of seizures, significant cardiac",
      "or inflammatory disease, primary hematologic disorder causing anemia,",
      "prior rHuEPO, > 2 RBC transfusions within 4 weeks of study drug",
      "assignment, or any RBC transfusion during the current chemotherapy cycle",
      "before randomization. Baseline hemoglobin 9.96 +/- 1.00 g/dL, range",
      "5.5-12.3 (Agoram 2006 Table 2). Tumor types: breast 31%, lung 17%,",
      "gastrointestinal 17%, gynecologic 17%, genitourinary 7%, lymphoma 2%,",
      "other 9%. Hemoglobin data from study 990146 were excluded from the PD",
      "analysis because its hemoglobin entry criterion (<= 13.0 g/dL) differed",
      "substantially from the other studies."
    ),
    dose_range     = paste(
      "0.5 to 4.5 ug/kg SC QW, 3.0 to 9.0 ug/kg SC Q2W, 4.5 to 15 ug/kg SC",
      "Q3W and 9.0 to 18.0 ug/kg SC Q4W across studies 20010162, 980290 and",
      "980291; study 20010162 also contributed 6.75 ug/kg SC Q3W (Agoram 2006",
      "Table 1). The PD estimation data set spanned 0.5 ug/kg QW to 15 ug/kg",
      "Q3W."
    ),
    regions        = "Three Amgen-sponsored clinical studies (20010162, 980290, 980291); geographic sites not stated in the publication.",
    notes          = paste(
      "PkPd model development data set: 573 patients and 6356 hemoglobin",
      "observations pooled from Amgen studies 20010162 (n = 84), 980290",
      "(n = 228) and 980291 (n = 261). Hemoglobin was measured weekly, before",
      "darbepoetin alfa administration on dosing days, by standard complete",
      "blood count; assay variability 1-5%. Hemoglobin data for 28 days after",
      "an RBC transfusion were censored. For safety the protocols withheld",
      "dosing when hemoglobin reached 14 g/dL (women) or 15 g/dL (men) and",
      "resumed below 13 g/dL; the paper cautions that this design is why IIV",
      "on Smax could not be estimated and that simulations of titration",
      "algorithms with higher hemoglobin targets should be treated carefully.",
      "Estimation used NONMEM V level 1.2 with the first-order (FO) method,",
      "because FOCE did not converge; no bootstrap was run for the PkPd model",
      "because of its computational cost. The model was evaluated by an",
      "external predictive check against a further 302 patients and 1253",
      "hemoglobin observations from the second parts of studies 980290 and",
      "980291 (Agoram 2006 Figure 4). Demographics from Table 2, trial",
      "summary from Table 1."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Pharmacokinetic layer. Agoram 2006 fitted the PD stage sequentially:
    # "the Pk model consisted of fixed-effects parameters obtained from the
    # full covariate model; random IIV and RRV in the Pk model were ignored"
    # (Results, "PkPd Model"). Every PK entry is therefore fixed() at its
    # Table 3 point estimate and carries no eta. Values converted from
    # mL and mL/day to L and L/day so that (ug in central) / (L) gives
    # ug/L = ng/mL directly. See Agoram_2006_darbepoetin_alfa_cia for the
    # freely estimated PK model with its full random-effects structure.
    # ---------------------------------------------------------------------
    lcl     <- fixed(log(2.010)); label("Typical clearance at WT 70 kg and PCNT <= 2 (L/day)")  # Agoram 2006 Table 3: theta_CL = 2010 mL/day
    lvc     <- fixed(log(3.390)); label("Typical central volume at WT 70 kg (L)")               # Agoram 2006 Table 3: theta_V1 = 3390 mL
    lvp     <- fixed(log(0.251)); label("Peripheral volume (L)")                                # Agoram 2006 Table 3: theta_V2 = 251 mL
    lq      <- fixed(log(2.900)); label("Intercompartmental clearance (L/day)")                 # Agoram 2006 Table 3: theta_Q = 2900 mL/day
    lka     <- fixed(log(0.318)); label("First-order SC absorption rate constant (1/day)")      # Agoram 2006 Table 3: theta_Ka = 0.318 1/day
    lfdepot <- fixed(log(0.443)); label("Subcutaneous bioavailability (fraction)")              # Agoram 2006 Table 3: F = 0.443

    e_conmed_platin_gt2_cl <- fixed(0.737); label("Multiplicative effect of > 2 concomitant platinum chemotherapy cycles on CL (unitless)")  # Agoram 2006 Table 3: theta_1 (PCNT on CL) = 0.737
    e_wt_cl                <- fixed(0.623); label("Power exponent of (WT/70 kg) on clearance (unitless)")                                   # Agoram 2006 Table 3: theta_2 (BWT on CL) = 0.623
    e_wt_vc                <- fixed(0.639); label("Power exponent of (WT/70 kg) on central volume (unitless)")                              # Agoram 2006 Table 3: theta_3 (BWT on V1) = 0.639

    # ---------------------------------------------------------------------
    # Pharmacodynamic layer (Agoram 2006 Table 4, "Full PkPd Model Parameter
    # Estimates"). All five structural PD parameters were estimated.
    # ---------------------------------------------------------------------
    lmtt      <- log(4.68);  label("Progenitor-compartment mean transit time RBCPT (day)")                       # Agoram 2006 Table 4: theta_RBCPT = 4.68 days, SE 1.26 (26.9% RSE)
    lmtt_rbc  <- log(82.2);  label("Mature red-cell lifespan RBCLS, split over four transit compartments (day)")  # Agoram 2006 Table 4: theta_RBCLS = 82.2 days, SE 11.3 (13.7% RSE)
    lemax     <- log(0.437); label("Incremental maximum stimulation of hemoglobin production Smax (fraction)")    # Agoram 2006 Table 4: theta_SMAX = 0.437, SE 0.0699 (16.0% RSE)
    lec50     <- log(3.68);  label("Darbepoetin alfa concentration at half-maximal stimulation S50 (ng/mL)")      # Agoram 2006 Table 4: theta_S50 = 3.68 ng/mL, SE 1.42 (38.6% RSE)
    lrbase_hb <- log(9.92);  label("Baseline total hemoglobin Hb0 (g/dL)")                                        # Agoram 2006 Table 4: theta_Hb0 = 9.92 g/dL, SE 0.0413 (0.416% RSE)

    e_conmed_platin_ec50 <- 1.86; label("Multiplicative effect of any concomitant platinum chemotherapy on S50 (unitless)")  # Agoram 2006 Table 4: theta_4 (PCNT on S50) = 1.86, SE 0.448 (24.1% RSE), asymptotic 95% CI (0.982, 2.74); NULL value 1

    # IIV - log-normal (Agoram 2006 Equation 1, P_i = TVP * exp(eta)); Table 4
    # reports omega^2 VARIANCES directly. The paper was "not able to estimate
    # random IIV in RBCPT and Smax", so those two carry no eta. The Omega
    # block carries one off-diagonal, between S50 and Hb0. The S50 and RBCLS
    # variances are extreme (CV > 100%; omega^2_S50 = 8.05 corresponds to a
    # CV of ~5700%) and the paper attributes this to the high heterogeneity
    # of the CIA population and to S50 absorbing the unestimated Smax
    # variability. Simulated cohorts are consequently very heavy-tailed.
    etalmtt_rbc ~ 1.63                     # Agoram 2006 Table 4: omega^2_RBCLS = 1.63 (30.0% RSE)
    etalec50 + etalrbase_hb ~ c(8.05,
                                0.0539, 0.00720)  # Agoram 2006 Table 4: omega^2_S50 = 8.05 (33.9% RSE), omega^2_S50,Hb0 = 0.0539 (38.8% RSE), omega^2_Hb0 = 0.00720 (7.78% RSE)

    # Residual error on hemoglobin. Table 4 reports sigma2^2 = 0.401 without
    # naming the error model; Agoram 2006 Equation 2 defines a log-scale
    # exponential residual for the PHARMACOKINETIC observations only. An
    # additive reading, SD = sqrt(0.401) = 0.633 g/dL (~6% of the 9.92 g/dL
    # baseline), is the only one consistent with Figure 3A, where the
    # individual-prediction scatter about the line of identity is roughly
    # +/- 1 g/dL; a log-scale reading of the same variance would be a 63% CV
    # and would spread predictions over 3-34 g/dL. See the vignette
    # "Assumptions and deviations" section.
    addSd_hb <- sqrt(0.401); label("Additive residual SD on hemoglobin (g/dL)")  # Agoram 2006 Table 4: sigma2^2 = 0.401, SE 0.0179 (4.46% RSE)
  })

  model({
    # -- PK layer: typical values only, no etas (see ini() note) -----------
    cl <- exp(lcl) * (e_conmed_platin_gt2_cl^CONMED_PLATIN_GT2) * (WT / 70)^e_wt_cl
    vc <- exp(lvc) * (WT / 70)^e_wt_vc
    vp <- exp(lvp)
    q  <- exp(lq)
    ka <- exp(lka)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # -- PD parameters -----------------------------------------------------
    # Equation 14: TVS50 = theta_S50 * (theta_4)^X, X = 1 if PCNT > 0.
    mtt      <- exp(lmtt)
    mtt_rbc  <- exp(lmtt_rbc + etalmtt_rbc)
    emax     <- exp(lemax)
    ec50     <- exp(lec50 + etalec50) * (e_conmed_platin_ec50^CONMED_PLATIN)
    rbase_hb <- exp(lrbase_hb + etalrbase_hb)

    # First-order transit rate constants. Agoram 2006 text after Equation 6:
    # kPT = 1 / RBCPT for the single progenitor compartment, and
    # kLS = 4 / RBCLS for the four equal-transit-time lifespan compartments.
    kpt <- 1 / mtt
    kls <- 4 / mtt_rbc

    # Baseline production rate. Agoram 2006 gives the indirect-response
    # baseline condition in Equations 5 and 6 (dHb/dt = 0 at Hb = Hb0). For
    # the transit-chain realisation the drug-free steady state is
    # precursor1 = Rin / kPT and hb1 = ... = hb4 = Rin / kLS, so
    # Hb0 = 4 * Rin / kLS = Rin * RBCLS, giving Rin = Hb0 / RBCLS.
    rin <- rbase_hb / mtt_rbc

    # -- Two-compartment PK disposition (Agoram 2006 Figure 1A) ------------
    # IV doses enter central directly (data set cmt = central); SC doses
    # enter depot (data set cmt = depot) with bioavailability F.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                               k12 * central - k21 * peripheral1

    f(depot) <- exp(lfdepot)

    # Serum darbepoetin alfa concentration, the PD driver C_DA of Equation 8.
    # Unlike the companion PK model this does NOT add the endogenous-EPO
    # assay contribution: Agoram 2006 states that "the effect of endogenous
    # EPO level fluctuation on the stimulation of hemoglobin production rate
    # was ignored", and Equation 8 is written in terms of the darbepoetin
    # alfa concentration alone.
    Cc <- central / vc

    # Equation 8. Smax is the INCREMENTAL maximum stimulation, so the
    # production multiplier saturates at 1 + Smax = 1.437.
    stim <- 1 + emax * Cc / (ec50 + Cc)

    # -- Erythropoiesis catenary chain (Agoram 2006 Figure 1B, Equations 7,
    # 9, 10 and 11). precursor1 is the paper's progenitor compartment CPT1;
    # hb1-hb4 are the four mature red-cell lifespan compartments LS1-LS4
    # (CPT2-CPT5), each holding the hemoglobin carried by red cells of that
    # age bin.
    d/dt(precursor1) <- rin * stim - kpt * precursor1
    d/dt(hb1)        <- kpt * precursor1 - kls * hb1
    d/dt(hb2)        <- kls * hb1        - kls * hb2
    d/dt(hb3)        <- kls * hb2        - kls * hb3
    d/dt(hb4)        <- kls * hb3        - kls * hb4

    precursor1(0) <- rin * mtt
    hb1(0)        <- rbase_hb / 4
    hb2(0)        <- rbase_hb / 4
    hb3(0)        <- rbase_hb / 4
    hb4(0)        <- rbase_hb / 4

    # Equation 11: total hemoglobin is the sum over the lifespan
    # compartments only; the progenitor compartment does not carry
    # circulating hemoglobin.
    hb <- hb1 + hb2 + hb3 + hb4
    hb ~ add(addSd_hb)
  })
}
