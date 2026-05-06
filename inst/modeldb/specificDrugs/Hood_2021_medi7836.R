Hood_2021_medi7836 <- function() {
  description <- "Population PK-PD binding model for MEDI7836 (anti-IL13 IgG1 lambda-YTE mAb) in healthy adult males (Hood 2021): two-compartment SC PK with first-order absorption, ADA-on-CL covariate, plus IL13 turnover, fixed Kon/Koff binding to MEDI7836:IL13 complex, complex distribution sharing CL/Q/V3 with parent drug, and a serum PD observation modelled as the molar sum of free IL13 and a small fraction of complex."
  reference <- "Hood T, Slob CJ, Slager J, Yates JWT, Bouzom F, Mistry HB. Pharmacokinetic-Pharmacodynamic Modelling of Systemic IL13 Blockade by Monoclonal Antibody Therapy: A Free Assay Disguised as Total. Pharmaceutics. 2021;13(5):613. doi:10.3390/pharmaceutics13050613"
  vignette <- "Hood_2021_medi7836"
  units <- list(time = "day", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    ADA_POS = list(
      description        = "Post-baseline anti-drug antibody (ADA) status",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ADA-negative post-baseline)",
      notes              = "1 = subject became ADA-positive post-baseline at any visit (Hood 2021 Table 1: 21/24 active subjects, 87.5%); 0 = ADA-negative throughout. Linear fractional effect on apparent clearance per Hood 2021 Eq. 5: CL/F = (CL/F)_pop * (1 + e_ada_pos_cl * ADA_POS). The paper used post-baseline ADA status as a fixed (not time-varying) covariate, which the authors discuss as a limitation that may have biased the estimate.",
      source_name        = "ADA"
    )
  )

  population <- list(
    n_subjects     = 32L,
    n_active       = 24L,
    n_placebo      = 8L,
    n_studies      = 1L,
    study          = "NCT02388347 first-in-human phase 1 single-ascending-dose trial",
    age_range      = "18-50 years",
    age_median     = "35 years (Table 1, overall)",
    weight_range   = "60.4-95.2 kg (Table 1, overall IQR)",
    weight_median  = "76.5 kg (Table 1, overall)",
    sex_female_pct = 0,
    sex_note       = "Healthy adult males only (study inclusion criterion)",
    race_ethnicity = "White 19/32 (59.4%), Black/African American 8/32 (25.0%), Asian 5/32 (15.6%) per Table 1.",
    disease_state  = "Healthy adult male volunteers (no asthma, no other indication).",
    dose_range     = "Single SC dose: 30, 105, 300, or 600 mg (six active subjects per cohort) or placebo (eight subjects).",
    regions        = "United Kingdom (single site).",
    ada_status     = "21/24 active subjects became ADA-positive post-baseline (87.5%); 16/24 (67%) classified as persistent ADA-positive.",
    notes          = "Demographics from Hood 2021 Table 1. 13.9 PK samples per individual on average across 24 active subjects; 14.9 PD (IL13) samples per individual across all 32 subjects. Two subjects whose PD profiles appeared swapped (one placebo, one 600 mg) were excluded from the PD analysis dataset. Follow-up extended to Day 281 because of the YTE half-life-extension expectation."
  )

  ini({
    # Structural PK parameters -- apparent (no IV data, F not estimated separately).
    # All values from Hood 2021 Table 2, "Estimate" column; CIs are 5-95% bootstrap.
    lka <- log(0.156);    label("First-order SC absorption rate Ka (1/day)")                              # Hood 2021 Table 2: Ka = 0.156 [0.137-0.183]
    lcl <- log(0.441);    label("Apparent clearance CL/F for ADA-negative subject (L/day)")               # Hood 2021 Table 2: CL/F = 0.441 [0.366-0.587]
    lvc <- log(2.83);     label("Apparent central volume V2/F for free MEDI7836 (L)")                     # Hood 2021 Table 2: V2/F = 2.83 [2.21-4.06]
    lvp <- log(8.03);     label("Apparent peripheral volume V3/F (shared with complex peripheral) (L)")   # Hood 2021 Table 2: V3/F = 8.03 [6.82-9.52]
    lq  <- log(0.825);    label("Apparent inter-compartmental clearance Q/F (L/day, shared with complex)")# Hood 2021 Table 2: Q/F = 0.825 [0.638-1.13]

    # Structural PD parameters (IL13 turnover, complex distribution, assay fraction).
    lkin   <- log(0.0173);  label("IL13 production rate Kin (nM/day; paper labels nmol/d but the steady-state baseline 0.096 pM = Kin/Kout requires nM/day)") # Hood 2021 Table 2: Kin = 0.0173 [0.0136-0.0227]
    lkout  <- log(180);     label("IL13 elimination rate constant Kout (1/day; estimated, not constrained to MEDI7836 kel)") # Hood 2021 Table 2: Kout = 180 [143-227]
    lvc_complex   <- log(13.6);    label("Apparent central volume of IL13:MEDI7836 complex Vcx/F (L)")           # Hood 2021 Table 2: Vcx/F = 13.6 [10.5-16.7]
    lcxfr  <- log(0.0429);  label("Fraction of total IL13:MEDI7836 complex captured by the bioanalytical PD assay (unitless)") # Hood 2021 Table 2: Cx fraction = 4.29% [2.98-5.96%]

    # Fixed in-vitro binding constants from Biacore T100 molecule characterization.
    kon  <- fixed(138.24); label("Binding rate Kon (1/(nM*day)), fixed from Biacore")                     # Hood 2021 Section 2.5: Kon = 138.24 nM^-1 day^-1
    koff <- fixed(0.69);   label("Dissociation rate Koff (1/day), fixed from Biacore")                    # Hood 2021 Section 2.5: Koff = 0.69 day^-1

    # Covariate: ADA-positive post-baseline on CL/F (linear fractional, Hood Eq. 5).
    e_ada_pos_cl <- 0.717; label("Fractional increase in CL/F for ADA-positive post-baseline subjects (unitless)") # Hood 2021 Table 2: ADA CL Increase = 71.7% [17.3-155%]

    # IIV (log-normal); paper reports CV%, so omega^2 = log(1 + CV^2). Diagonal omega
    # used here: Hood 2021 mentions correlations between PK omegas were estimated
    # "if supported by the data," but Table 2 reports only diagonal entries and
    # neither the paper nor the supplement publishes the off-diagonals -- so this
    # implementation uses a diagonal IIV matrix and notes the deviation in the
    # vignette's "Assumptions and deviations" section.
    etalcl   ~ 0.250303   # Hood 2021 Table 2: IIV CL/F = 53.3% -> log(1 + 0.533^2) = 0.250303
    etalvc   ~ 0.424439   # Hood 2021 Table 2: IIV V2/F = 72.7% -> log(1 + 0.727^2) = 0.424439
    etalvp   ~ 0.193389   # Hood 2021 Table 2: IIV V3/F = 46.2% -> log(1 + 0.462^2) = 0.193389
    etalq    ~ 0.395029   # Hood 2021 Table 2: IIV Q/F  = 69.6% -> log(1 + 0.696^2) = 0.395029
    etalkout ~ 0.282271   # Hood 2021 Table 2: IIV Kout = 57.1% -> log(1 + 0.571^2) = 0.282271
    etalcxfr ~ 1.215281   # Hood 2021 Table 2: IIV Cx Fraction = 154% -> log(1 + 1.54^2) = 1.215281

    # Residual error.
    # PK is a combined model (Hood Eq. 2: y = pred * (1 + eps_prop) + eps_add) on
    # linear-space ng/mL concentrations. PD is proportional on the log-transformed
    # observation (Hood Eq. 3: y = log(IL13 + FR*Cx) * (1 + eps_prop)).
    addSd        <- 0.0546; label("Additive PK residual error (ng/mL = ug/L)")                     # Hood 2021 Table 2: e1 = 0.0546 ug/L
    propSd       <- 0.130;  label("Proportional PK residual error (fraction)")                     # Hood 2021 Table 2: e2 = 13.0%
    propSd_pdIL13  <- 0.069;  label("Proportional PD residual error on log(serum IL13 readout) (fraction)") # Hood 2021 Table 2: e3 = 6.9%
  })

  model({
    # Stoichiometric / unit constants for the binding sub-model. MEDI7836 is a
    # ~150 kDa IgG1 lambda-YTE mAb (Hood 2021 Methods); IL13 is the much smaller
    # cytokine target, so binding kinetics are written on the molar (nM) scale
    # while drug compartments are tracked in mass units (mg).
    mwdrug   <- 0.150  # MW of MEDI7836 (mg/nmol = kg/mol/1000); 1 nmol = 0.150 mg = 150 ug
    v_target   <- 1      # Implicit reference volume (L) for the IL13 turnover
                       # compartment so that target (state) is numerically equal to
                       # IL13 concentration in nM. With v_target = 1 L, the paper's
                       # Kin = 0.0173 (labelled "nmol/d") gives the reported
                       # baseline serum free IL13 = Kin/Kout = 0.096 pM exactly.

    # Individual PK parameters.
    ka  <- exp(lka)
    cl  <- exp(lcl + etalcl) * (1 + e_ada_pos_cl * ADA_POS)
    vc  <- exp(lvc + etalvc)
    vp  <- exp(lvp + etalvp)
    q   <- exp(lq  + etalq)

    # Individual PD parameters (no IIV on Kin or Vcx per Hood 2021 Table 2).
    kin  <- exp(lkin)
    kout <- exp(lkout + etalkout)
    vc_complex  <- exp(lvc_complex)
    cxfr <- exp(lcxfr + etalcxfr)

    # Micro-rate constants for the linear distribution / elimination kinetics.
    # The complex shares CL/F and Q/F with the parent drug but has its own central
    # volume Vcx/F; its peripheral compartment shares V3/F (Hood 2021 Section 3.2).
    kel    <- cl / vc
    k12    <- q  / vc
    k21    <- q  / vp
    kel_cx <- cl / vc_complex
    k12_cx <- q  / vc_complex
    k21_cx <- q  / vp

    # Drug mass balance (compartments in mg of MEDI7836).
    # Binding rate, expressed as a drug-mass loss rate, is Kon * central * target
    # (mg/day). Derivation: Cdrug_nM = (central_mg / vc) / mwdrug; binding rate
    # in V2 reference = vc * Kon * Cdrug_nM * Ctarget = Kon * central_mg * target /
    # mwdrug (nmol/day); times mwdrug (mg/nmol) gives Kon * central * target mg/day.
    # Dissociation gives back drug at rate Koff * complex_nmol * mwdrug (mg/day).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1 -
                          kon * central * target + koff * complex * mwdrug
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # IL13 free-target turnover (state in nM via implicit v_target = 1 L).
    # Binding/dissociation amount-rates from the drug-mass balance above are
    # divided by v_target to give the corresponding IL13 concentration rates.
    # Initial condition is the no-drug steady state Kin/Kout (paper: 0.096 pM).
    d/dt(target)        <-  kin - kout * target -
                          (kon * central * target / mwdrug) / v_target +
                          (koff * complex) / v_target
    target(0)           <-  kin / kout

    # IL13:MEDI7836 complex compartments (amounts in nmol).
    d/dt(complex)   <-  kon * central * target / mwdrug - koff * complex -
                          kel_cx * complex - k12_cx * complex + k21_cx * complex_peripheral
    d/dt(complex_peripheral)   <-  k12_cx * complex - k21_cx * complex_peripheral

    # Observation 1: free MEDI7836 in serum, ng/mL = ug/L.
    # central is in mg, vc in L, so central / vc is mg/L = ug/mL. Multiply by
    # 1000 to express in ng/mL = ug/L (the unit reported in Hood 2021 Figure 2 /
    # Table 2 epsilon1).
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)

    # Observation 2: serum IL13 PD readout. Hood 2021 Eq. 3 models the recorded
    # log-IL13 as the molar sum of free IL13 and a small fraction (cxfr) of the
    # central complex concentration; the paper's interpretation is that the
    # supposedly-free IL13 immunoassay also captures ~4% of the bound drug-target
    # complex, biasing the apparent free signal upward when drug is present.
    pdIL13 <- log(target + cxfr * complex / vc_complex)
    pdIL13 ~ prop(propSd_pdIL13)
  })
}
