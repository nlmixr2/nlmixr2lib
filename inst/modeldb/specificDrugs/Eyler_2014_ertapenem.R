Eyler_2014_ertapenem <- function() {
  description <- "Two-compartment population PK model for intravenous ertapenem in critically ill adults with acute kidney injury receiving continuous venovenous hemodialysis (CVVHD) or hemodiafiltration (CVVHDF). PK is parameterised on unbound drug; total serum concentrations are reconstructed via a single-site saturable albumin-binding equation Cb = Bmax * Cu / (KD + Cu). Systemic (body) clearance and a separate dialytic clearance arm are estimated as primary parameters; the dialytic arm is added to body clearance only while the CRRT circuit is running, gated by the time-varying RRT_HEMODIAL_ACTIVE covariate. Eyler 2014, n = 8 subjects, single 1 g IV dose over 30 min."
  reference <- "Eyler RF, Vilay AM, Nader AM, Heung M, Pleva M, Sowinski KM, DePestel DD, Sorgel F, Kinzig M, Mueller BA. Pharmacokinetics of ertapenem in critically ill patients receiving continuous venovenous hemodialysis or hemodiafiltration. Antimicrob Agents Chemother. 2014;58(3):1320-1326. doi:10.1128/AAC.02090-12"
  vignette <- "Eyler_2014_ertapenem"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    RRT_HEMODIAL_ACTIVE = list(
      description        = "CRRT-active indicator (1 while continuous venovenous hemodialysis or hemodiafiltration is running, 0 otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CRRT circuit off / paused)",
      notes              = "Time-varying within subject. Gates the additive dialytic-arm clearance cl_dialysis: the CRRT contribution is added to body systemic clearance only while the circuit is running (Eyler 2014 Methods, Pharmacokinetic analysis: 'An indicator variable, DIAL, with a value of 1 or 0, was used to turn the effluent compartment on and off, respectively, if the CVVHD/F was turned off for any reason'). Three of eight subjects had CRRT paused mid-sampling (135, 142, and 276 minutes; Results paragraph 1) for filter changes or an off-floor procedure; RRT_HEMODIAL_ACTIVE flips to 0 during those windows. For subjects whose circuit ran continuously across sampling, RRT_HEMODIAL_ACTIVE = 1 throughout. The paper does not distinguish CVVHD vs CVVHDF in the model: both modalities are encoded by RRT_HEMODIAL_ACTIVE = 1, with the dialytic-arm clearance (CLdial = 36 mL/min) reflecting the pooled-cohort mean of dialytic solute removal at the cohort mean effluent rate (38 ml/h/kg). This is the time-varying per-session dialysis gate (canonical RRT_HEMODIAL_ACTIVE), distinct from the static RRT_HEMODIAL_STATUS subject-level indicator.",
      source_name        = "DIAL",
      source_alias       = "HEMODIALYSIS (working column name before the 2026-06-19 canonical-register standardization; renamed to the time-varying gate canonical RRT_HEMODIAL_ACTIVE)"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested as a power covariate (centered on the cohort median 80 kg) on CLS, CLdial, and VC during covariate screening (Eyler 2014 Methods, Covariate testing equation 10) but not retained in the final model (Table 2 reports only the structural typical-value estimates; no WT exponents). The cohort weight range was 56.0-119.2 kg (mean 78.9 +/- 19.8 kg, Table 1). Listed in covariatesDataExcluded to preserve the screened-but-not-retained provenance.",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested as a power covariate on CLS and on the CLS interindividual-variability estimate (Eyler 2014 Methods, Covariate testing paragraph 2) but not retained in the final model. The cohort age range was 31-78 years (mean 62 +/- 16 years, Table 1).",
      source_name        = "AGE"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 8L,
    n_studies        = 1L,
    age_range        = "31-78 years",
    age_median       = "70 years (cohort mean 62 +/- 16 years; Table 1)",
    weight_range     = "56.0-119.2 kg",
    weight_median    = "80 kg (cohort mean 78.9 +/- 19.8 kg; Table 1)",
    sex_female_pct   = 62.5,
    race_ethnicity   = "Not reported (single-centre study at the University of Michigan, USA)",
    disease_state    = "Critically ill adult ICU patients with acute kidney injury and suspected or confirmed Gram-negative infection, requiring continuous renal replacement therapy. Severity: APACHE III scores 63-123 (mean 83 +/- 19 for the 7 subjects with reported scores). Mean serum albumin 3.0 +/- 0.5 g/dL. Urine output minimal (< 50 mL / 24 h) in all subjects.",
    dose_range       = "Single 1 g ertapenem (Merck) as a 30-min IV infusion (clinical-care dose for empiric Gram-negative coverage).",
    regions          = "United States (single-centre, University of Michigan)",
    renal_function   = "Acute kidney injury requiring continuous renal replacement therapy. Modalities: CVVHD (n=4, subjects 1-4 with low ultrafiltration) and CVVHDF (n=4, subjects 5-8). Mean blood flow rate 181 +/- 26 mL/min; mean dialysate flow rate 24 +/- 10 mL/h/kg; mean ultrafiltration rate 14 +/- 8 mL/h/kg; mean effluent rate (dialysate + ultrafiltrate) 38 +/- 10 mL/h/kg. Mean sieving coefficient (effluent / pre-filter serum) 0.21 +/- 0.06. Three of eight subjects had CRRT paused (135-276 minutes) during the 24 h sampling interval; the pauses are reflected by RRT_HEMODIAL_ACTIVE = 0 in those windows.",
    notes            = "ClinicalTrials.gov NCT00877370. Cohort enrolled April 2009 - March 2011. Bioanalytical: total ertapenem serum concentrations by HPLC-MS/MS at IBMP Nuernberg-Heroldsberg (Sorgel laboratory); unbound concentrations by equilibrium dialysis; effluent concentrations from CRRT effluent samples. LLOQ 0.5 ug/mL; inter- and intra-day precision < 12% across matrices. Pre-filter serum concentrations were adjusted for citrate dilution (Methods equation 1) before modeling."
  )

  ini({
    # Structural PK parameters (Eyler 2014 Table 2, final-model column).
    # Typical values reported in mL/min for clearances and L for volumes;
    # mL/min -> L/h conversion factor is x60/1000 = x0.06:
    #   CLS_unbound     = 48  mL/min = 2.88 L/h
    #   CLD_unbound     = 115 mL/min = 6.90 L/h
    #   CLdial_unbound  = 36  mL/min = 2.16 L/h
    # All parameters are parameterised on unbound concentration; the
    # 'apparent' VC absorbs the unbound-fraction scaling (the same
    # convention as Fauchet_2015_lopinavir_unbound).
    lcl          <- log(2.88); label("Systemic body unbound clearance CLS (L/h)")               # Eyler 2014 Table 2: CLS = 48 mL/min (RSE 10%)
    lvc          <- log(32);   label("Apparent central volume of distribution VC (L)")          # Eyler 2014 Table 2: VC = 32 L     (RSE 167%)
    lvp          <- log(21);   label("Apparent peripheral volume of distribution VP (L)")       # Eyler 2014 Table 2: VP = 21 L     (RSE 23%)
    lq           <- log(6.9);  label("Inter-compartmental distribution clearance CLD (L/h)")    # Eyler 2014 Table 2: CLD = 115 mL/min (RSE 41%)
    lcl_dialysis <- log(2.16); label("CRRT dialytic clearance CLdial (L/h)")                    # Eyler 2014 Table 2: CLdial = 36 mL/min (RSE 13%)

    # Saturable single-site protein-binding model linking unbound and
    # total serum concentrations (Eyler 2014 Methods, Pharmacokinetic
    # analysis: 'Ertapenem total concentrations were linked to unbound
    # concentrations using a nonlinear maximum binding model').
    # Standard saturable form: Cb = Bmax * Cu / (KD + Cu); Ct = Cu + Cb.
    lbmax <- log(144); label("Maximum protein binding capacity Bmax (mg/L)")                 # Eyler 2014 Table 2: Bmax = 144 ug/mL (RSE 26%)
    lkd   <- log(38);  label("Binding equilibrium dissociation constant KD (mg/L)")          # Eyler 2014 Table 2: KD   = 38  ug/mL (RSE 25%)

    # Interindividual variability (Eyler 2014 Table 2 'Interindividual
    # variability (%)' column; exponential random-effect model, %CV on the
    # natural scale; variance omega^2 = log(CV^2 + 1)).
    #   CLS     : 23% CV -> omega^2 = log(0.23^2 + 1) = 0.05163
    #   VC      : 33% CV -> omega^2 = log(0.33^2 + 1) = 0.10324
    #   VP      : 20% CV -> omega^2 = log(0.20^2 + 1) = 0.03922
    #   CLdial  : 32% CV -> omega^2 = log(0.32^2 + 1) = 0.09752
    #   Bmax    : 17% CV -> omega^2 = log(0.17^2 + 1) = 0.02842
    # CLD and KD: IIV reported as 'NA' (not estimated).
    etalcl          ~ 0.05163  # Eyler 2014 Table 2 (IIV CLS,    23% CV)
    etalvc          ~ 0.10324  # Eyler 2014 Table 2 (IIV VC,     33% CV)
    etalvp          ~ 0.03922  # Eyler 2014 Table 2 (IIV VP,     20% CV)
    etalcl_dialysis ~ 0.09752  # Eyler 2014 Table 2 (IIV CLdial, 32% CV)
    etalbmax        ~ 0.02842  # Eyler 2014 Table 2 (IIV Bmax,   17% CV)

    # Residual error: proportional on both unbound and total outputs
    # (Eyler 2014 Methods, Pharmacokinetic analysis: 'Residual unexplained
    # variability ... was modeled using a proportional error term
    # (equation 9) with separate terms for unbound, total, and effluent
    # concentrations'). Magnitudes are NOT reported in Table 2 or anywhere
    # else in the paper; fixed at a typical 20% as an extraction
    # assumption (see vignette Assumptions and deviations).
    propSd          <- fixed(0.20); label("Proportional residual error on total Cc (fraction; ASSUMED, see vignette Errata)")        # Eyler 2014: magnitude not reported; fixed at 20% as extraction assumption
    propSd_Cunbound <- fixed(0.20); label("Proportional residual error on unbound Cunbound (fraction; ASSUMED, see vignette Errata)") # Eyler 2014: magnitude not reported; fixed at 20% as extraction assumption
  })

  model({
    # Individual PK parameters (apparent unbound values; the structural
    # clearance/volume terms are reported by the paper as parameterised
    # on unbound concentration).
    cl           <- exp(lcl          + etalcl)
    vc           <- exp(lvc          + etalvc)
    vp           <- exp(lvp          + etalvp)
    q            <- exp(lq)                            # IIV reported NA (Table 2)
    cl_dialysis  <- exp(lcl_dialysis + etalcl_dialysis)

    # Saturable binding model parameters.
    bmax <- exp(lbmax + etalbmax)
    kd   <- exp(lkd)                                   # IIV reported NA (Table 2)

    # Total body clearance during the active CRRT session is
    # cl + cl_dialysis (Eyler 2014 Methods, equation 2:
    # dx(1)/dt = R0 - (CLS/VC)*X(1) - (CLD/VC)*X(1) + (CLD/VP)*X(2)
    #                                       - (CLdial/VC)*DIAL*X(1)).
    # The dialytic arm is added to the body baseline only when the
    # RRT_HEMODIAL_ACTIVE covariate is 1 (CRRT circuit running).
    cl_total <- cl + RRT_HEMODIAL_ACTIVE * cl_dialysis

    # Rate-constant form of the two-compartment system.
    kel <- cl_total / vc
    k12 <- q / vc
    k21 <- q / vp

    # ODE system. The central state carries the unbound-equivalent
    # amount of ertapenem (the apparent VC absorbs the unbound-fraction
    # scaling, as in Fauchet 2015 lopinavir unbound). 1 g IV infusion
    # over 30 min enters 'central' via the event-table dose row.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # Unbound serum concentration (mg/L = ug/mL).
    Cunbound <- central / vc

    # Total serum concentration via the single-site saturable albumin-
    # binding equation (Eyler 2014 Methods, Pharmacokinetic analysis;
    # see Fig 3 of the paper for the empirical unbound-fraction vs total
    # relationship). Cb = Bmax * Cu / (KD + Cu); Ct = Cu + Cb.
    Cc <- Cunbound + bmax * Cunbound / (kd + Cunbound)

    # Two-output residual error (proportional on each scale; both
    # magnitudes fixed at 20% as extraction assumptions because the
    # paper does not report the residual-error SDs).
    Cc       ~ prop(propSd)
    Cunbound ~ prop(propSd_Cunbound)
  })
}
