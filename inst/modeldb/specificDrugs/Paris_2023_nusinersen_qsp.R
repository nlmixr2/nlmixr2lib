Paris_2023_nusinersen_qsp <- function() {
  description <- paste(
    "QSP. Paediatric neurofilament-trafficking model of phosphorylated",
    "neurofilament heavy subunit (pNfH) in spinal muscular atrophy (SMA), with",
    "and without intrathecal nusinersen, from birth to 20 years (Paris 2023).",
    "Seven pNfH concentration states (CNS neurons, PNS motor neurons, brain",
    "ISF, endoneurial fluid, cranial CSF, spinal CSF, blood) carry first-order",
    "synthesis/degradation, leakage and CSF/blood transport, with volume-growth",
    "dilution terms from seven growing compartment volumes driven by",
    "time-varying body weight (WT) and height (HT). CNS synthesis and leakage",
    "follow a synaptic-pruning function of age; SMA adds an age-dependent PNS",
    "leakage r14 = A + B * C^(-age). Nusinersen amounts in the three",
    "spinal-cord segments come from the embedded nine-compartment Biliouris",
    "2018 PK model (typical rate constants) and pass through a four-stage",
    "first-order delay chain whose output inhibits r14 as",
    "1 - ASOeff / (ASOeff + EC50). Model time is postnatal age in days;",
    "simulate with covsInterpolation = 'linear'."
  )
  reference <- paste(
    "Paris A, Bora P, Parolo S, MacCannell D, Monine M, van der Munnik N,",
    "Tong X, Eraly S, Berger Z, Graham D, Ferguson T, Domenici E, Nestorov I,",
    "Marchetti L. A pediatric quantitative systems pharmacology model of",
    "neurofilament trafficking in spinal muscular atrophy treated with the",
    "antisense oligonucleotide nusinersen. CPT Pharmacometrics Syst Pharmacol.",
    "2023;12(2):196-206. doi:10.1002/psp4.12890"
  )
  vignette <- "Paris_2023_nusinersen_qsp"
  units <- list(
    time = "day",
    dosing = "mg",
    concentration = "ng/mL"
  )

  # The pNfH model tracks concentrations (ng/mL) in seven fluid/tissue spaces
  # and the volumes (mL) of those spaces as ODE states; the nusinersen PK
  # sub-model reuses the Biliouris 2018 state names; the four-stage delay
  # chain uses the numbered `effect<n>` family.
  paper_specific_compartments <- c(
    "pnfh_cns",
    "pnfh_pns",
    "pnfh_isf",
    "pnfh_endo",
    "pnfh_csf_cranial",
    "pnfh_csf_spinal",
    "pnfh_blood",
    "vol_cns",
    "vol_pns",
    "vol_isf",
    "vol_endo",
    "vol_csf_cranial",
    "vol_csf_spinal",
    "vol_blood",
    "spinal_cord_cervical",
    "spinal_cord_lumbar",
    "spinal_cord_thoracic",
    "pons"
  )

  compartmentData <- list(
    csf = list(analyte = "nusinersen", units = "ng", specimen = "CSF", verified = TRUE),
    central = list(analyte = "nusinersen", units = "ng", specimen = "plasma", verified = TRUE),
    spinal_cord_cervical = list(analyte = "nusinersen", units = "ng", specimen = "tissue", verified = TRUE),
    brain = list(analyte = "nusinersen", units = "ng", specimen = "tissue", verified = TRUE),
    peripheral1 = list(analyte = "nusinersen", units = "ng", specimen = "tissue", verified = TRUE),
    spinal_cord_lumbar = list(analyte = "nusinersen", units = "ng", specimen = "tissue", verified = TRUE),
    brain_deep = list(analyte = "nusinersen", units = "ng", specimen = "tissue", verified = TRUE),
    spinal_cord_thoracic = list(analyte = "nusinersen", units = "ng", specimen = "tissue", verified = TRUE),
    pons = list(analyte = "nusinersen", units = "ng", specimen = "tissue", verified = TRUE),
    effect1 = list(analyte = "nusinersen", units = "ng", specimen = "not applicable", verified = TRUE),
    effect2 = list(analyte = "nusinersen", units = "ng", specimen = "not applicable", verified = TRUE),
    effect3 = list(analyte = "nusinersen", units = "ng", specimen = "not applicable", verified = TRUE),
    effect4 = list(analyte = "nusinersen", units = "ng", specimen = "not applicable", verified = TRUE),
    pnfh_cns = list(analyte = "pNfH", units = "ng/mL", specimen = "tissue", verified = TRUE),
    pnfh_pns = list(analyte = "pNfH", units = "ng/mL", specimen = "tissue", verified = TRUE),
    pnfh_isf = list(analyte = "pNfH", units = "ng/mL", specimen = "brain ISF", verified = TRUE),
    pnfh_endo = list(analyte = "pNfH", units = "ng/mL", specimen = "tissue", verified = TRUE),
    pnfh_csf_cranial = list(analyte = "pNfH", units = "ng/mL", specimen = "CSF", verified = TRUE),
    pnfh_csf_spinal = list(analyte = "pNfH", units = "ng/mL", specimen = "CSF", verified = TRUE),
    pnfh_blood = list(analyte = "pNfH", units = "ng/mL", specimen = "plasma", verified = TRUE),
    vol_cns = list(analyte = "volume", units = "mL", specimen = "not applicable", verified = TRUE),
    vol_pns = list(analyte = "volume", units = "mL", specimen = "not applicable", verified = TRUE),
    vol_isf = list(analyte = "volume", units = "mL", specimen = "not applicable", verified = TRUE),
    vol_endo = list(analyte = "volume", units = "mL", specimen = "not applicable", verified = TRUE),
    vol_csf_cranial = list(analyte = "volume", units = "mL", specimen = "not applicable", verified = TRUE),
    vol_csf_spinal = list(analyte = "volume", units = "mL", specimen = "not applicable", verified = TRUE),
    vol_blood = list(analyte = "volume", units = "mL", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight at the current age (time-varying)",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives the growth of the CNS, ISF and blood volumes (Appendix S1",
        "Formulas S1, S3, S7). The authors supply the CDC clinical growth",
        "chart 50th-percentile male weight-for-age curve (Appendix S2",
        "Growth_charts/, linearly interpolated). Must be supplied as a",
        "time-varying column and solved with covsInterpolation = 'linear';",
        "the volume growth rates are read off the movement of WT between",
        "records."
      ),
      source_name = "BW"
    ),
    HT = list(
      description = "Body height (length or stature) at the current age (time-varying)",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives the growth of the PNS and endoneurial-fluid volumes",
        "(proportional to height, Formulas S2 and S4) and enters the blood",
        "volume (Formula S7). The authors supply the CDC 50th-percentile male",
        "length/stature-for-age curve. Time-varying; solve with",
        "covsInterpolation = 'linear'."
      ),
      source_name = "BH"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 457,
    n_studies = 5,
    age_range = "birth to 22 years (model simulated from birth to 20 years)",
    disease_state = paste(
      "Spinal muscular atrophy (SMA), untreated and treated with intrathecal",
      "nusinersen; plus healthy paediatric controls"
    ),
    dose_range = "Intrathecal nusinersen 1 to 12.5 mg, individual loading and maintenance schedules",
    administration_routes = "Intrathecal bolus (lumbar puncture) into CSF",
    regions = "Multinational (Biogen SHINE, NURTURE, EMBRACE, ENDEAR, CHERISH trials)",
    notes = paste(
      "98 healthy paediatric subjects (ENDEAR) re-estimated the healthy leakage",
      "and blood-clearance rates; 359 SMA patients (17 untreated controls plus",
      "pre-treatment samples of treated patients) defined a 400-member untreated",
      "virtual population (Table S1); 278 treated patients were fitted",
      "individually for the treatment parameters and 25 held out for validation",
      "in four dose/age groups (Table S2). Measurements are pNfH in CSF and",
      "plasma. Parameters were obtained by normalised least-squares fits of",
      "deterministic trajectories (deSolve/optim), not by nonlinear",
      "mixed-effects estimation."
    )
  )

  ini({
    # ---- Neurofilament synthesis / degradation (r1-r4) ----
    # r2 = r4 = ln(2)/thalf_stat * f_stat + ln(2)/thalf_mov * (1 - f_stat)
    # = 0.0101 1/day (Table S4 r2, r4 '0.01 1/day').
    thalf_nf_stat <- fixed(90); label("Half-life of stationary neurofilament (day)") # Appendix S2 pNfH_SMA_model_parameters.R NF_halflife_stat = 3*30; Table S4 r2/r4 refs [5,6]
    thalf_nf_mov <- fixed(22); label("Half-life of moving neurofilament (day)") # Appendix S2 NF_halflile_moving = 22
    f_nf_stat <- fixed(0.9); label("Fraction of neurofilament that is stationary (unitless)") # Appendix S2 NF_static_fraction = 0.9

    # ---- Leakage from neurons (r5 CNS -> ISF, r6 PNS -> endoneurial fluid) ----
    lk_leak <- log(1.7e-6); label("Neurofilament leakage rate constant, r5 = r6 before the age factor (1/day)") # Table S4 r5 1.7 x 10-6 (healthy data), r6 assumed equal; Table S5; Appendix S2 k_leakage_CNS_motor_H = 1.7e-06

    # ---- Fluid transport (mL/day) ----
    q_isf_csf <- fixed(449.28); label("ISF -> cranial CSF flow, r7 (mL/day)") # Appendix S2 k_ISF_to_CSF = 449.28; Table S4 r8 '449.28 ml/day, assumed equal to r7'
    q_endo_csf <- fixed(449.28); label("Endoneurial fluid -> spinal CSF flow, r8 (mL/day)") # Table S4 r8 449.28 ml/day
    q_cran_spin <- fixed(7.4712); label("Cranial -> spinal CSF flow, r9 (mL/day)") # Table S4 r9 7.4712 ml/day
    q_spin_cran <- fixed(20.148); label("Spinal -> cranial CSF flow, r10 (mL/day)") # Table S4 r10 20.148 ml/day
    q_cran_blood <- fixed(90); label("Cranial CSF -> blood flow, r11 (mL/day)") # Table S4 r11 90 ml/day
    q_spin_blood <- fixed(259.2); label("Spinal CSF -> blood flow, r12 (mL/day)") # Table S4 r12 259.2 ml/day

    # ---- pNfH clearance from blood (r13), SMA untreated population ----
    lkel_pnfh <- log(0.09860226); label("pNfH clearance rate from blood, r13 (1/day)") # Appendix S2 SMA_untreated k_blood_clearance_H = 9.860226e-02; Table 1 / Table S4 0.10 +/- 0.01 (untreated population)

    # ---- SMA-specific PNS leakage r14 = A + B * C^(-age) ----
    r14_a <- 1.6e-7; label("Age-independent SMA PNS leakage term, r14A (1/day)") # Table 1 / Table S4 r14A (1.6 +/- 1.0) x 10-7 (untreated population); Appendix S2 r14_A
    r14_b <- 1.8e-5; label("Age-dependent SMA PNS leakage amplitude, r14B (1/day)") # Table 1 / Table S4 r14B (1.8 +/- 0.6) x 10-5 (untreated population); Appendix S2 r14_B
    lr14_c <- log(3.3); label("Base of the age decay of SMA leakage, r14C (unitless, per year of age)") # Table 1 / Table S4 r14C 3.3 +/- 0.9 (untreated population); Appendix S2 r14_C

    # ---- Synaptic-pruning function f(age) = P1 exp(-P2 age) + P3 exp(-P4 age) + P5 ----
    p1_prune <- -541.4912301; label("Pruning-function coefficient P1 (unitless)") # Appendix S2 P1 = -541.4912301; Table S4 P1 -541.5
    p2_prune <- 0.1597563; label("Pruning-function rate P2 (1/year)") # Appendix S2 P2 = 0.1597563; Table S4 prints 0.161
    p3_prune <- 527.1709528; label("Pruning-function coefficient P3 (unitless)") # Appendix S2 P3 = 527.1709528; Table S4 P3 527.2
    p4_prune <- 0.1341329; label("Pruning-function rate P4 (1/year)") # Appendix S2 P4 = 0.1341329; Table S4 prints 0.131
    p5_prune <- 42.2240297; label("Pruning-function constant P5 (unitless)") # Appendix S2 P5 = 42.2240297; Table S4 P5 42.2

    # ---- Nusinersen effect: four-stage delay chain + inhibition of r14 ----
    lktr <- log(0.09); label("Diffusion rate through the nusinersen delay chain, kdiff (1/day)") # Table 1 / Table S4 kdiff 0.09 +/- 0.08 1/day (individual treatment fits)
    lec50 <- log(550); label("Spinal-cord nusinersen amount giving half-maximal inhibition of r14, EC50 (ng)") # Table 1 / Table S4 EC50 (5.5 +/- 2.1) x 10^2 'ng/ml'; the driver is the spinal-cord AMOUNT (ng), see vignette

    # ---- pNfH concentrations at birth, average SMA untreated scenario (ng/mL) ----
    bl_pnfh_cns <- fixed(1e5); label("pNfH in CNS neurons at birth (ng/mL)") # Table S6 1 x 10^8 pg/ml
    bl_pnfh_pns <- fixed(2e6); label("pNfH in PNS motor neurons at birth (ng/mL)") # Table S6 2 x 10^9 pg/ml
    bl_pnfh_isf <- fixed(0.018); label("pNfH in brain ISF at birth (ng/mL)") # Table S6 18 pg/ml
    bl_pnfh_endo <- fixed(3.4); label("pNfH in endoneurial fluid at birth, SMA (ng/mL)") # Table S6 340 pg/ml (Appendix S2 3.4 ng/ml)
    bl_pnfh_csf_cranial <- fixed(1.15); label("pNfH in cranial CSF at birth, SMA (ng/mL)") # Table S6 1150 pg/ml
    bl_pnfh_csf_spinal <- fixed(5.4); label("pNfH in spinal CSF at birth, SMA (ng/mL)") # Table S6 5400 pg/ml
    bl_pnfh_blood <- fixed(47); label("pNfH in blood at birth, SMA (ng/mL)") # Table S6 47000 pg/ml

    # ---- Compartment volumes at birth (mL) ----
    v0_cns <- fixed(338); label("CNS neuron volume at birth (mL)") # Table S3 338 ml
    v0_pns <- fixed(38.23); label("PNS neuron volume at birth (mL)") # Table S3 38.23 ml
    v0_isf <- fixed(50.7); label("Brain ISF volume at birth (mL)") # Table S3 50.7 ml (Appendix S2 338*0.15)
    v0_endo <- fixed(25.88); label("Endoneurial fluid volume at birth (mL)") # Table S3 25.9 ml (Appendix S2 25.88)
    v0_csf_cranial <- fixed(60); label("Cranial CSF volume at birth (mL)") # Table S3 60 ml
    v0_csf_spinal <- fixed(60); label("Spinal CSF volume at birth (mL)") # Table S3 60 ml
    v0_blood <- fixed(267.329); label("Blood volume at birth (mL)") # Table S3 267 ml (Appendix S2 267.329)

    # ---- Nusinersen PK (Biliouris 2018 typical rate constants, 1/h) ----
    lk_csf_cerv <- fixed(log(1.71e-3)); label("CSF -> cervical spinal cord rate K13 (1/h)") # Appendix S2 PK_ASO_model.R K13 = 1.71e-3
    lk_cerv_csf <- fixed(log(1e-4)); label("Cervical spinal cord -> CSF rate K31 (1/h)") # Appendix S2 K31 = 1e-4
    lk_csf_brain <- fixed(log(6e-3)); label("CSF -> brain rate K14 (1/h)") # Appendix S2 K14 = 6e-3
    lk_brain_csf <- fixed(log(4e-4)); label("Brain -> CSF rate K41 (1/h)") # Appendix S2 K41 = 4e-4
    lk_csf_plasma <- fixed(log(8.91e-2)); label("CSF -> plasma rate K12 (1/h)") # Appendix S2 K12 = 8.91e-2
    lkel <- fixed(log(0.206)); label("Plasma elimination rate K20 (1/h)") # Appendix S2 K20 = 0.206
    lk_plasma_periph <- fixed(log(8.18e-3)); label("Plasma -> peripheral tissue rate K25 (1/h)") # Appendix S2 K25 = 8.18e-3
    lk_periph_plasma <- fixed(log(1e-4)); label("Peripheral tissue -> plasma rate K52 (1/h)") # Appendix S2 K52 = 1e-4
    lk_csf_lumb <- fixed(log(2.86e-3)); label("CSF -> lumbar spinal cord rate K16 (1/h)") # Appendix S2 K16 = 2.86e-3
    lk_lumb_csf <- fixed(log(3e-4)); label("Lumbar spinal cord -> CSF rate K61 (1/h)") # Appendix S2 K61 = 3e-4
    lk_brain_deep <- fixed(log(2.57e-3)); label("Brain -> deep brain tissue rate K47 (1/h)") # Appendix S2 K47 = 2.57e-3
    lk_deep_brain <- fixed(log(1e-4)); label("Deep brain tissue -> brain rate K74 (1/h)") # Appendix S2 K74 = 1e-4
    lk_csf_thor <- fixed(log(2.1e-3)); label("CSF -> thoracic spinal cord rate K18 (1/h)") # Appendix S2 K18 = 2.1e-3
    lk_thor_csf <- fixed(log(4.5e-4)); label("Thoracic spinal cord -> CSF rate K81 (1/h)") # Appendix S2 K81 = 4.5e-4
    lk_csf_pons <- fixed(log(1.57e-3)); label("CSF -> pons rate K19 (1/h)") # Appendix S2 K19 = 1.57e-3
    lk_pons_csf <- fixed(log(2e-4)); label("Pons -> CSF rate K91 (1/h)") # Appendix S2 K91 = 2e-4

    # ---- Between-subject spread (log-normal, moment-matched to Table 1 mean +/- SD) ----
    # omega^2 = log(1 + (SD/mean)^2); r13 and r14 from the untreated
    # virtual-population column, EC50 and kdiff from the treatment-fit column.
    etalkel_pnfh ~ fixed(0.00995) # Table 1 r13 0.10 +/- 0.01 (CV 10%)
    etar14_a ~ fixed(0.3298) # Table 1 r14A (1.6 +/- 1.0) x 10-7 (CV 62.5%)
    etar14_b ~ fixed(0.1054) # Table 1 r14B (1.8 +/- 0.6) x 10-5 (CV 33.3%)
    etalr14_c ~ fixed(0.0717) # Table 1 r14C 3.3 +/- 0.9 (CV 27.3%)
    etalec50 ~ fixed(0.1361) # Table 1 EC50 (5.5 +/- 2.1) x 10^2 (CV 38.2%)
    etalktr ~ fixed(0.5823) # Table 1 kdiff 0.09 +/- 0.08 (CV 88.9%)
  })

  model({
    # ---- Age (years); model time t is postnatal age in days ----
    age_yr <- t / 365
    # The authors evaluate the volume targets one day ahead (finite-difference
    # step delta = 1/365 year) and relax each volume state onto its target at
    # rate 1/delta = 365 per day (Appendix S2 pNfH_SMA_model.R).
    age_lead <- age_yr + 1 / 365
    krelax <- 365

    # ---- Neurofilament turnover ----
    kdeg_nf <- log(2) / thalf_nf_stat * f_nf_stat + log(2) / thalf_nf_mov * (1 - f_nf_stat)

    # Synaptic-pruning function (Formula S16) and its age derivative
    prune <- p1_prune * exp(-p2_prune * age_yr) + p3_prune * exp(-p4_prune * age_yr) + p5_prune
    prune0 <- p1_prune + p3_prune + p5_prune
    prune20 <- p1_prune * exp(-p2_prune * 20) + p3_prune * exp(-p4_prune * 20) + p5_prune
    dprune <- -p1_prune * p2_prune * exp(-p2_prune * age_yr) - p3_prune * p4_prune * exp(-p4_prune * age_yr)

    # Synthesis at equilibrium with degradation for the birth concentrations,
    # scaled by the pruning shape (Formula S17; Appendix S2 applies it to both
    # r1 and r3)
    r1 <- kdeg_nf * bl_pnfh_cns * prune / prune0
    r3 <- kdeg_nf * bl_pnfh_pns * prune / prune0

    # Age factor on leakage (Yilmaz 2017 formula, held at its 20-year value
    # below 20 years; Table S4 f(age) = 0.0975 x 1.031^20)
    fage <- 0.0975 * 1.031^20
    if (age_yr >= 20) {
      fage <- 0.0975 * 1.031^age_yr
    }
    k_leak <- exp(lk_leak)
    r5 <- k_leak * fage * (1 - dprune / prune20) # Formula S18
    r6 <- k_leak * fage

    kel_pnfh <- exp(lkel_pnfh + etalkel_pnfh)

    # ---- SMA leakage r14 (Equation 1) ----
    r14a <- r14_a * exp(etar14_a)
    r14b <- r14_b * exp(etar14_b)
    r14c <- exp(lr14_c + etalr14_c)
    r14 <- r14a + r14b * r14c^(-age_yr)

    # ---- Nusinersen PK (Biliouris 2018), amounts in ng, rate constants 1/h -> 1/day ----
    k_csf_cerv <- 24 * exp(lk_csf_cerv)
    k_cerv_csf <- 24 * exp(lk_cerv_csf)
    k_csf_brain <- 24 * exp(lk_csf_brain)
    k_brain_csf <- 24 * exp(lk_brain_csf)
    k_csf_plasma <- 24 * exp(lk_csf_plasma)
    kel <- 24 * exp(lkel)
    k_plasma_periph <- 24 * exp(lk_plasma_periph)
    k_periph_plasma <- 24 * exp(lk_periph_plasma)
    k_csf_lumb <- 24 * exp(lk_csf_lumb)
    k_lumb_csf <- 24 * exp(lk_lumb_csf)
    k_brain_deep <- 24 * exp(lk_brain_deep)
    k_deep_brain <- 24 * exp(lk_deep_brain)
    k_csf_thor <- 24 * exp(lk_csf_thor)
    k_thor_csf <- 24 * exp(lk_thor_csf)
    k_csf_pons <- 24 * exp(lk_csf_pons)
    k_pons_csf <- 24 * exp(lk_pons_csf)

    d/dt(csf) <- -(k_csf_cerv + k_csf_brain + k_csf_plasma + k_csf_lumb + k_csf_thor + k_csf_pons) * csf +
      k_cerv_csf * spinal_cord_cervical + k_brain_csf * brain + k_lumb_csf * spinal_cord_lumbar +
      k_thor_csf * spinal_cord_thoracic + k_pons_csf * pons
    d/dt(central) <- k_csf_plasma * csf - kel * central - k_plasma_periph * central + k_periph_plasma * peripheral1
    d/dt(spinal_cord_cervical) <- k_csf_cerv * csf - k_cerv_csf * spinal_cord_cervical
    d/dt(brain) <- k_csf_brain * csf - k_brain_csf * brain - k_brain_deep * brain + k_deep_brain * brain_deep
    d/dt(peripheral1) <- k_plasma_periph * central - k_periph_plasma * peripheral1
    d/dt(spinal_cord_lumbar) <- k_csf_lumb * csf - k_lumb_csf * spinal_cord_lumbar
    d/dt(brain_deep) <- k_brain_deep * brain - k_deep_brain * brain_deep
    d/dt(spinal_cord_thoracic) <- k_csf_thor * csf - k_thor_csf * spinal_cord_thoracic
    d/dt(pons) <- k_csf_pons * csf - k_pons_csf * pons
    # Intrathecal dose in mg -> amount in ng
    f(csf) <- 1e6

    # ---- Delay chain (Equations S20-S23) driven by total spinal-cord nusinersen ----
    aso_sc <- spinal_cord_cervical + spinal_cord_lumbar + spinal_cord_thoracic
    ktr <- exp(lktr + etalktr)
    d/dt(effect1) <- ktr * (aso_sc - effect1)
    d/dt(effect2) <- ktr * (effect1 - effect2)
    d/dt(effect3) <- ktr * (effect2 - effect3)
    d/dt(effect4) <- ktr * (effect3 - effect4)
    ec50 <- exp(lec50 + etalec50)
    inh_r14 <- 1 - effect4 / (effect4 + ec50) # Equation S24

    # ---- Volume targets (Formulas S1-S7) and growth rates (mL/day) ----
    vt_cns <- (1.449 - 3.62 / WT) * 1000
    vt_isf <- vt_cns * 0.15
    vt_pns <- HT / 170 * 130
    vt_endo <- HT / 170 * 88
    # CSF cranial = CSF spinal = half of the total paediatric CSF volume,
    # linearly interpolated between 120, 130, 135, 140, 150 mL at
    # 0, 0.375, 0.75, 1.5, 2 years (constant afterwards)
    if (age_lead < 0.375) {
      vt_csf <- (120 + (130 - 120) * age_lead / 0.375) / 2
    } else if (age_lead < 0.75) {
      vt_csf <- (130 + (135 - 130) * (age_lead - 0.375) / 0.375) / 2
    } else if (age_lead < 1.5) {
      vt_csf <- (135 + (140 - 135) * (age_lead - 0.75) / 0.75) / 2
    } else if (age_lead < 2) {
      vt_csf <- (140 + (150 - 140) * (age_lead - 1.5) / 0.5) / 2
    } else {
      vt_csf <- 75
    }
    if (age_lead < 2) {
      vt_blood <- 10^(0.7891 * log10(WT) + 0.004132 * HT + 1.8117)
    } else if (age_lead <= 14) {
      vt_blood <- 10^(0.6459 * log10(WT) + 0.00282 * HT + 2.09)
    } else {
      vt_blood <- (13.1 * HT + 18.05 * WT - 480) / 0.5723
    }

    dv_cns <- krelax * (vt_cns - vol_cns)
    dv_isf <- krelax * (vt_isf - vol_isf)
    dv_pns <- krelax * (vt_pns - vol_pns)
    dv_endo <- krelax * (vt_endo - vol_endo)
    dv_csf_cranial <- krelax * (vt_csf - vol_csf_cranial)
    dv_csf_spinal <- krelax * (vt_csf - vol_csf_spinal)
    dv_blood <- krelax * (vt_blood - vol_blood)
    # Neural and CSF volumes stop growing at 18 years, blood at 20 years
    if (age_yr >= 18) {
      dv_cns <- 0
      dv_isf <- 0
      dv_pns <- 0
      dv_endo <- 0
      dv_csf_cranial <- 0
      dv_csf_spinal <- 0
    }
    if (age_lead >= 20) {
      dv_blood <- 0
    }

    d/dt(vol_cns) <- dv_cns
    d/dt(vol_pns) <- dv_pns
    d/dt(vol_isf) <- dv_isf
    d/dt(vol_endo) <- dv_endo
    d/dt(vol_csf_cranial) <- dv_csf_cranial
    d/dt(vol_csf_spinal) <- dv_csf_spinal
    d/dt(vol_blood) <- dv_blood
    vol_cns(0) <- v0_cns
    vol_pns(0) <- v0_pns
    vol_isf(0) <- v0_isf
    vol_endo(0) <- v0_endo
    vol_csf_cranial(0) <- v0_csf_cranial
    vol_csf_spinal(0) <- v0_csf_spinal
    vol_blood(0) <- v0_blood

    # ---- pNfH concentration ODEs (Equations S9-S15) ----
    d/dt(pnfh_cns) <- r1 - kdeg_nf * pnfh_cns - r5 * pnfh_cns - pnfh_cns * dv_cns / vol_cns
    d/dt(pnfh_pns) <- r3 - kdeg_nf * pnfh_pns - r6 * pnfh_pns - r14 * inh_r14 * pnfh_pns -
      pnfh_pns * dv_pns / vol_pns
    d/dt(pnfh_isf) <- r5 * vol_cns / vol_isf * pnfh_cns - q_isf_csf / vol_isf * pnfh_isf -
      pnfh_isf * dv_isf / vol_isf
    d/dt(pnfh_endo) <- r6 * vol_pns / vol_endo * pnfh_pns - q_endo_csf / vol_endo * pnfh_endo +
      r14 * inh_r14 * pnfh_pns * vol_pns / vol_endo - pnfh_endo * dv_endo / vol_endo
    d/dt(pnfh_csf_cranial) <- q_isf_csf / vol_csf_cranial * pnfh_isf -
      q_cran_spin / vol_csf_cranial * pnfh_csf_cranial + q_spin_cran / vol_csf_cranial * pnfh_csf_spinal -
      q_cran_blood / vol_csf_cranial * pnfh_csf_cranial - pnfh_csf_cranial * dv_csf_cranial / vol_csf_cranial
    d/dt(pnfh_csf_spinal) <- q_cran_spin / vol_csf_spinal * pnfh_csf_cranial -
      q_spin_cran / vol_csf_spinal * pnfh_csf_spinal - q_spin_blood / vol_csf_spinal * pnfh_csf_spinal +
      q_endo_csf / vol_csf_spinal * pnfh_endo - pnfh_csf_spinal * dv_csf_spinal / vol_csf_spinal
    d/dt(pnfh_blood) <- q_cran_blood / vol_blood * pnfh_csf_cranial + q_spin_blood / vol_blood * pnfh_csf_spinal -
      kel_pnfh * pnfh_blood - pnfh_blood * dv_blood / vol_blood
    pnfh_cns(0) <- bl_pnfh_cns
    pnfh_pns(0) <- bl_pnfh_pns
    pnfh_isf(0) <- bl_pnfh_isf
    pnfh_endo(0) <- bl_pnfh_endo
    pnfh_csf_cranial(0) <- bl_pnfh_csf_cranial
    pnfh_csf_spinal(0) <- bl_pnfh_csf_spinal
    pnfh_blood(0) <- bl_pnfh_blood

    # ---- Outputs in the paper's reporting units (pg/mL) ----
    pnfh_blood_pgml <- 1000 * pnfh_blood
    pnfh_csf_pgml <- 1000 * pnfh_csf_spinal
  })
}
