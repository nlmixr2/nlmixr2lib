Farrell_2014_eltrombopag <- function() {
  description <- "Population PK/PD model for eltrombopag in healthy male volunteers (single dose) and adult patients with chronic liver disease (CLD; multiple daily doses) (Farrell 2014). Two-compartment apparent disposition with dual sequential first-order absorption: Ka1 acts on the depot from the end of the absorption lag time (ALAG1) until time MTIME after the dose, and Ka2 acts thereafter. CL/F is reduced in females, in East Asian subjects, and in CLD patients with a linear-in-Child-Pugh-score gradient (HEPIMP_CP_SCORE >= 5). Vc/F is approximately three-fold higher in South/Central Asian subjects. Platelet dynamics use a four-compartment lifespan model (three maturing precursor pools feeding the circulating-platelet pool) with linear stimulation of precursor production by plasma eltrombopag; the slope SLOP is 34% lower in East Asian CLD patients. PD parameters are CLD-specific (median baseline platelet 41 Gi/L)."
  reference <- "Farrell C, Hayes S, Wire M, Zhang J. Population pharmacokinetic/pharmacodynamic modelling of eltrombopag in healthy volunteers and subjects with chronic liver disease. Br J Clin Pharmacol. 2014 May;77(5):717-728. doi:10.1111/bcp.12244."
  vignette <- "Farrell_2014_eltrombopag"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL", platelet = "Gi/L (10^9/L)")

  covariateData <- list(
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Multiplicative fractional change on CL/F: females have 38% lower apparent clearance than males (Farrell 2014 Table 3: 'CL/F ~ Females' multiplier 0.622, 95% CI 0.495-0.749). The discussion attributes the female effect to lower CYP1A2 and UGT activity (Farrell 2014 Discussion paragraph 5).",
      source_name        = "Sex (Female)"
    ),
    HEPIMP_CP_SCORE = list(
      description        = "Child-Pugh composite score for chronic liver disease severity (integer)",
      units              = "(integer score, 5-15 for CLD; 0 by convention for healthy)",
      type               = "count",
      reference_category = "0 (healthy; no CLD)",
      notes              = "Child-Pugh integer score (Class A = 5-6, Class B = 7-9, Class C = 10-15). Healthy subjects without CLD are assigned HEPIMP_CP_SCORE = 0 so the model gates the CP effect off (no CL/F reduction). Linear-in-score effect on CL/F for CLD patients: factor = e_cp_cl * (1 + e_cp_slope * (HEPIMP_CP_SCORE - 5)) = 0.536 * (1 - 0.113 * (CP - 5)) (Farrell 2014 Table 3 footnote and rows 'CL/F ~ CP Score 5' and 'CL/F ~ CP Score > 5'). CP Class A patients (score 5-6) have 46-59% lower CL/F than the healthy reference; Class B patients (score 7-9) have 59-71% lower CL/F (Farrell 2014 Results, paragraph after Table 3).",
      source_name        = "Child-Pugh score"
    ),
    RACE_NEAS = list(
      description        = "North East Asian race indicator (Chinese, Japanese, or Korean heritage)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-North East Asian: White, Black, South/Central Asian, Other)",
      notes              = "Carries TWO multiplicative fractional-change effects in this model: (1) CL/F: East Asians have 52% lower apparent clearance than non-East-Asians (Farrell 2014 Table 3: 'CL/F ~ East Asians' multiplier 0.476, 95% CI 0.395-0.557). (2) SLOP: in CLD patients, East Asians have 34% lower linear-slope drug-induced platelet-production effect (Farrell 2014 Table 6: 'SLOP ~ East Asians' multiplier 0.660, 95% CI 0.529-0.791). In the source data set the East Asian subgroup is exclusively Japanese (Study 2, n=38) plus n=7 East Asian patients from Study 3.",
      source_name        = "East Asian"
    ),
    RACE_ASIAN_SOUTHCENTRAL = list(
      description        = "South / Central Asian race indicator (Indian, Pakistani, Bangladeshi, Sri Lankan, Nepali, or Central Asian heritage)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-South/Central Asian: White, Black, East Asian, Other)",
      notes              = "Multiplicative fractional change on Vc/F: South/Central Asian CLD patients have an apparent central volume of distribution approximately 3-fold higher than all other races (Farrell 2014 Table 3: 'Vc/F ~ South/Central Asians' multiplier 2.98, 95% CI 2.21-3.75). Farrell 2014 Discussion attributes this to lower median serum albumin (29 vs 35 umol/L in non-S/C-Asian CLD patients) and to a higher proportion of serial-sampling subjects in the South/Central Asian subgroup of Study 3. In the source data set, 21 of the 41 PK-substudy patients in Study 3 (51%) were South/Central Asian.",
      source_name        = "South/Central Asian"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 107L,
    n_subjects_cld  = 79L,
    n_studies       = 3L,
    age_range       = "19-81 years (median 50)",
    weight_range    = "40.6-102 kg (median 67.9)",
    sex_female_pct  = 27.0,
    race_ethnicity  = c(White = 35, BlackAfricanAmerican = 2, EastAsian = 42, SouthCentralAsian = 20, Other = 2),
    disease_state   = "Pooled: 28 healthy adult male volunteers (Study 1, single dose) and 79 thrombocytopenic adult male and female patients with chronic liver disease (Studies 2 and 3, repeat once-daily dosing). CLD patients had baseline platelet counts < 50 Gi/L (median 41); Child-Pugh classification A in 37 patients, B in 40, C in 2.",
    dose_range      = "Study 1 healthy volunteers: single oral 30, 50, or 75 mg. Study 2 Japanese CLD: 12.5, 25, or 37.5 mg PO QD x 14 d. Study 3 CLD: 75 mg PO QD x 14 d.",
    regions         = "Multinational; the East Asian (Japanese) CLD cohort is Study 2; Study 3 included a South/Central Asian CLD subgroup (51% of the PK substudy).",
    samples         = "786 plasma eltrombopag concentrations across the 107 subjects; 451 platelet count observations across the 79 CLD patients.",
    notes           = "PD parameters are calibrated to the CLD platelet response: the typical KIN (0.211 Gi/L/h) is much lower than the value previously reported for healthy and ITP populations (1.43); the typical baseline circulating platelet count is 41 Gi/L (Farrell 2014 Table 2 CLD median). For typical-value simulations of healthy subjects use the PK layer only (HEPIMP_CP_SCORE = 0, baseline platelet population not modelled here). Baseline demographics from Farrell 2014 Tables 1 and 2."
  )

  ini({
    # ---- PK structural parameters (Farrell 2014 Table 3, NONMEM point estimates) ----
    # Reference subject: White, male, healthy volunteer (HEPIMP_CP_SCORE = 0, SEXF = 0,
    # RACE_NEAS = 0, RACE_ASIAN_SOUTHCENTRAL = 0).
    lcl    <- log(0.953); label("Apparent clearance CL/F (L/h) for the reference White male healthy volunteer")        # Farrell 2014 Table 3
    lvc    <- log(9.75);  label("Apparent central volume Vc/F (L) for non-South/Central-Asian subjects")               # Farrell 2014 Table 3
    lvp    <- log(9.41);  label("Apparent peripheral volume Vp/F (L)")                                                 # Farrell 2014 Table 3
    lq     <- log(0.633); label("Apparent inter-compartmental clearance Q/F (L/h)")                                    # Farrell 2014 Table 3
    lka1   <- log(0.333); label("First-order absorption rate constant Ka1 (1/h), pre-MTIME slow absorption phase")     # Farrell 2014 Table 3
    lka2   <- log(5.26);  label("First-order absorption rate constant Ka2 (1/h), post-MTIME fast absorption phase")    # Farrell 2014 Table 3
    ltlag    <- log(0.447); label("Absorption lag time ALAG1 (h)")                                                       # Farrell 2014 Table 3
    ltswitch <- log(1.46);  label("Time after dose at which the absorption rate switches from Ka1 to Ka2 (h; paper MTIME)") # Farrell 2014 Table 3

    # ---- PK covariate effects (Farrell 2014 Table 3) ----
    # Fractional-change forms: parameter_i = parameter_pop * (1 + e * X) where X is the binary covariate.
    # Source paper reports the multiplicative factor m directly; here e = m - 1, so the factor at X = 1
    # reproduces the paper's m exactly.
    e_sexf_cl    <- -0.378; label("Fractional change in CL/F for SEXF = 1: factor = 1 + e_sexf_cl * SEXF = 0.622")                                  # Farrell 2014 Table 3 'CL/F ~ Females' multiplier 0.622
    e_neas_cl    <- -0.524; label("Fractional change in CL/F for RACE_NEAS = 1: factor = 1 + e_neas_cl * RACE_NEAS = 0.476")                        # Farrell 2014 Table 3 'CL/F ~ East Asians' multiplier 0.476
    e_cp_cl      <-  0.536; label("CL/F multiplier anchor at HEPIMP_CP_SCORE = 5 (Child-Pugh Class A bottom); CLD-only effect")                     # Farrell 2014 Table 3 'CL/F ~ CP Score 5' multiplier 0.536
    e_cp_slope   <- -0.113; label("Per-unit fractional decrement on CL/F for each unit of HEPIMP_CP_SCORE above 5")                                 # Farrell 2014 Table 3 footnote, 'CL/F ~ CP Score > 5' coefficient
    e_scasian_vc <-  1.98;  label("Fractional change in Vc/F for RACE_ASIAN_SOUTHCENTRAL = 1: factor = 1 + e_scasian_vc * indicator = 2.98")        # Farrell 2014 Table 3 'Vc/F ~ South/Central Asians' multiplier 2.98

    # ---- PK IIV (Farrell 2014 Table 3) ----
    # Correlated block on CL/F and Vc/F: omega^2_CL = 0.161, covar = 0.106, omega^2_Vc = 0.233
    # (CV_CL = 40.1%, CV_Vc = 48.3%, R = 0.547). Single IIV on Ka1 with omega^2 = 1.61 (CV 127%) -
    # the source paper labels this as interoccasion variability (IOV) on Ka1; it is encoded here as
    # IIV because the nlmixr2lib registry does not carry NONMEM-style occasion structure (see
    # vignette 'Assumptions and deviations').
    etalcl + etalvc ~ c(0.161, 0.106, 0.233)   # Farrell 2014 Table 3, correlated block on CL/Vc
    etalka1         ~ 1.61                     # Farrell 2014 Table 3 omega^2_Ka (paper labels IOV on Ka1)

    # ---- PK residual error (Farrell 2014 Table 3) ----
    # Combined proportional + additive on Cc (ug/mL).
    # propSd = sqrt(omega^2_prop) = sqrt(0.0237) = 0.154.
    # addSd  = sqrt(omega^2_add)  = sqrt(1110 ng^2/mL^2) = 33.3 ng/mL = 0.0333 ug/mL.
    # The paper's two supplementary factors sigma_Prop~HV (= 0.453, scaling propSd for healthy-volunteer
    # samples) and sigma_Prop~TAD<4h (= 1.72, scaling propSd for early-time-after-dose samples) are
    # simplified out of the registry model and documented in vignette 'Assumptions and deviations'.
    propSd <- 0.154;  label("Proportional residual error SD on Cc (fraction); CLD steady-state samples")  # Farrell 2014 Table 3 sigma^2_prop = 0.0237
    addSd  <- 0.0333; label("Additive residual error SD on Cc (ug/mL)")                                    # Farrell 2014 Table 3 sigma^2_add = 1110 ng^2/mL^2

    # ---- PD structural parameters (Farrell 2014 Table 6) ----
    # Calibrated to the CLD platelet response.
    lslop  <- log(0.648);  label("Linear slope of eltrombopag effect on platelet-precursor production (mL/ug)")     # Farrell 2014 Table 6 SLOP
    lkin   <- log(0.211);  label("Zero-order platelet-precursor production rate KIN (Gi/L/h)")                       # Farrell 2014 Table 6 KIN
    lkt    <- log(0.0214); label("First-order platelet-precursor maturation/transit rate KT (1/h)")                  # Farrell 2014 Table 6 KT
    lrbase <- log(41);     label("Typical baseline platelet count for CLD patients (Gi/L); kdeg is derived as kin/rbase") # Farrell 2014 Table 2 CLD median baseline platelet 41 Gi/L

    # ---- PD covariate effects (Farrell 2014 Table 6) ----
    e_neas_slop <- -0.340; label("Fractional change in SLOP for RACE_NEAS = 1: factor = 1 + e_neas_slop * RACE_NEAS = 0.660") # Farrell 2014 Table 6 'SLOP ~ East Asians' multiplier 0.660

    # ---- PD IIV (Farrell 2014 Table 6) ----
    # No IIV reported on KIN. Diagonal block on SLOP and KT.
    etalslop ~ 0.287   # Farrell 2014 Table 6 omega^2_SLOP (CV 53.6%)
    etalkt   ~ 0.313   # Farrell 2014 Table 6 omega^2_KT   (CV 55.9%)

    # ---- PD residual error (Farrell 2014 Table 6) ----
    # propSd_PLT = sqrt(omega^2_prop) = sqrt(0.197) = 0.444; matches the reported CV 44.4%.
    propSd_PLT <- 0.444; label("Proportional residual error SD on circulating platelet count (fraction)")          # Farrell 2014 Table 6 sigma^2_prop = 0.197
  })

  model({
    # ---- Derived CL/F covariate multiplier ----
    # CP-score gate: factor = 1 when HEPIMP_CP_SCORE = 0 (healthy / no CLD), else
    # e_cp_cl * (1 + e_cp_slope * (HEPIMP_CP_SCORE - 5)) = 0.536 * (1 - 0.113 * (CP - 5)).
    cp_factor <- 1
    if (HEPIMP_CP_SCORE > 0) cp_factor <- e_cp_cl * (1 + e_cp_slope * (HEPIMP_CP_SCORE - 5))

    # ---- Individual PK parameters ----
    cl  <- exp(lcl + etalcl) *
           (1 + e_sexf_cl * SEXF) *
           (1 + e_neas_cl * RACE_NEAS) *
           cp_factor
    vc  <- exp(lvc + etalvc) * (1 + e_scasian_vc * RACE_ASIAN_SOUTHCENTRAL)
    vp  <- exp(lvp)
    q   <- exp(lq)
    ka1 <- exp(lka1 + etalka1)
    ka2 <- exp(lka2)
    tlag    <- exp(ltlag)
    tswitch <- exp(ltswitch)

    # ---- Two-compartment disposition micro-constants ----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---- Dual sequential first-order absorption ----
    # The depot dose enters at t = T + tlag (alag below). From tad = tlag to tad = tswitch the
    # absorption rate constant is Ka1 (slow); for tad >= tswitch it is Ka2 (fast). tad() resets at
    # each dose event, so the same Ka1 -> Ka2 profile is repeated for every QD administration.
    # (`tswitch` is the paper's MTIME; renamed because `mtime` is reserved in rxode2.)
    ka_now <- ka1
    if (tad() >= tswitch) ka_now <- ka2

    # ---- ODE system: two-compartment PK ----
    d/dt(depot)       <- -ka_now * depot
    d/dt(central)     <-  ka_now * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Absorption lag applied at the dose event on the depot.
    alag(depot) <- tlag

    # Plasma eltrombopag concentration: dose mg / vc L = mg/L = ug/mL.
    Cc <- central / vc

    # ---- Individual PD parameters ----
    slop  <- exp(lslop + etalslop) * (1 + e_neas_slop * RACE_NEAS)
    kin   <- exp(lkin)
    kt    <- exp(lkt + etalkt)
    rbase <- exp(lrbase)
    kdeg  <- kin / rbase   # platelet first-order degradation rate, derived from kin and rbase

    # ---- ODE system: four-compartment platelet lifespan ----
    # Plasma eltrombopag stimulates precursor production linearly (Farrell 2014 Figure 3 schematic):
    #   d/dt(precursor1) = kin * (1 + slop * Cc) - kt * precursor1
    # Three sequential precursor pools (precursor1 -> precursor2 -> precursor3) feed circulating
    # platelets (circ); circulating platelets are removed with first-order rate kdeg. With
    # kdeg = kin/rbase the drug-free steady-state circulating value equals rbase by construction.
    d/dt(precursor1) <- kin * (1 + slop * Cc) - kt * precursor1
    d/dt(precursor2) <- kt * precursor1 - kt * precursor2
    d/dt(precursor3) <- kt * precursor2 - kt * precursor3
    d/dt(circ)       <- kt * precursor3 - kdeg * circ

    # ---- Steady-state initial conditions (drug-free) ----
    precursor1(0) <- kin / kt
    precursor2(0) <- kin / kt
    precursor3(0) <- kin / kt
    circ(0)       <- rbase

    # ---- Observations ----
    # Two algebraic observables: Cc (= central / vc, plasma eltrombopag in ug/mL)
    # and PLT (= circ, circulating platelet count in Gi/L). Event tables address
    # the observables by name (cmt = "Cc" for PK observations, cmt = "PLT" for
    # PD observations); rxode2 auto-injects the dvid mapping. Both columns are
    # returned in the rxSolve output regardless of which observable each
    # observation row addresses.
    PLT <- circ
    Cc  ~ add(addSd) + prop(propSd)
    PLT ~ prop(propSd_PLT)
  })
}
