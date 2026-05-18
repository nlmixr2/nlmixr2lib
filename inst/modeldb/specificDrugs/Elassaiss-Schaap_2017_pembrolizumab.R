`Elassaiss-Schaap_2017_pembrolizumab` <- function() {
  description <- "Two-compartment population PK model with parallel linear and Michaelis-Menten clearance plus a direct-response Imax PK/PD model on the ex vivo IL-2 stimulation ratio (PD-1 target engagement) for IV pembrolizumab (anti-PD-1 IgG4 mAb) in adults with advanced solid tumors (Elassaiss-Schaap 2017, KEYNOTE-001 parts A, A1, A2)."
  reference <- paste(
    "Elassaiss-Schaap J, Rossenu S, Lindauer A, Kang SP, de Greef R,",
    "Sachs JR, de Alwis DP. Using Model-Based 'Learn and Confirm' to",
    "Reveal the Pharmacokinetics-Pharmacodynamics Relationship of",
    "Pembrolizumab in the KEYNOTE-001 Trial. CPT Pharmacometrics Syst",
    "Pharmacol. 2017;6(1):21-28. doi:10.1002/psp4.12132"
  )
  vignette <- "Elassaiss-Schaap_2017_pembrolizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 22L,
    n_studies      = 1L,
    age_range      = "adults (>=18 years); specific range not reported",
    age_median     = "not reported",
    weight_range   = "not reported (mg/kg dosing used)",
    weight_median  = "not reported",
    sex_female_pct = NA_real_,
    race_ethnicity = "not reported",
    disease_state  = "Advanced solid tumors (KEYNOTE-001 parts A, A1, A2). PD-1, PD-L1, PD-L2, and CTLA-4 inhibitor-naive; no systemic corticosteroids required.",
    dose_range     = "0.005, 0.02, 0.06, 0.1, 1, 2, 3, or 10 mg/kg IV over 30 min; Q3W (parts A, A1) or escalated within-patient (part A2)",
    regions        = "United States (two sites)",
    notes          = "KEYNOTE-001 (NCT01295827); part A 3 + 3 dose escalation at 1, 3, 10 mg/kg Q2W (n = 9), part A1 expansion at 10 mg/kg Q2W (n = 7 additional), part A2 within-patient escalation from 0.005-0.02 mg/kg up to 2 or 10 mg/kg Q3W (design C; 12 patients). Elassaiss-Schaap 2017 Methods, Results 'Part A2', and Table 2."
  )

  ini({
    # Structural PK parameters (Elassaiss-Schaap 2017 Table 2 final model).
    # The paper notes the final fixed-effect parameters were log-transformed
    # ("now using log-transformation for parameters"), so each lXXX is the
    # log of the reported point estimate.
    lcl   <- log(0.168);    label("Linear clearance CL_lin (L/day)")                                  # Table 2 CL_lin = 0.168 L/d
    lvc   <- log(2.88);     label("Central volume of distribution Vc (L)")                            # Table 2 Vc    = 2.88 L
    lq    <- log(0.384);    label("Inter-compartmental clearance Q (L/day)")                          # Table 2 Q     = 0.384 L/d
    lvp   <- log(2.85);     label("Peripheral volume of distribution Vp (L)")                         # Table 2 Vp    = 2.85 L
    lvmax <- log(0.114);    label("Maximum Michaelis-Menten elimination rate Vmax (mg/day)")          # Table 2 Vmax  = 0.114 mg/d
    lkm   <- log(0.0784);   label("Michaelis-Menten constant Km (ug/mL)")                             # Table 2 Km    = 0.0784 ug/mL
    lfdepot <- fixed(log(1)); label("Bioavailability F (unitless; fixed at 1)")                       # Table 2 F = 1 (fixed)

    # PD parameters (Imax model on IL-2 stimulation ratio; Table 2 footnote
    # d: "Exponent of the estimated parameter" - the listed value is the
    # back-transformed point estimate from the log-transformed theta).
    # Note: paper text confirms IC50 ~ 0.54 mg/L for the final model, matching
    # the listed Table 2 value 0.535 ug/mL directly (i.e., 0.535 is the
    # point estimate, not exp(0.535)).
    lbaseline <- log(2.09);  label("Baseline IL-2 stimulation ratio (unitless)")                       # Table 2 Base  = 2.09
    limax     <- log(0.961); label("Maximal fractional inhibition of margin above 1 (unitless)")       # Table 2 Imax  = 0.961
    lic50     <- log(0.535); label("Pembrolizumab concentration producing 50% Imax (ug/mL)")           # Table 2 IC50  = 0.535 ug/mL

    # IIV (Elassaiss-Schaap 2017 Table 2). omega^2 = log(CV^2 + 1) for
    # log-normal BSV per the paper's stated CV-to-variance formula.
    # Vmax    : BSV 22.7% -> omega^2 = log(0.227^2 + 1) = log(1.05153) = 0.05024
    # F       : IOV 37.7% -> omega^2 = log(0.377^2 + 1) = log(1.14213) = 0.13289
    #   F's variability is reported as interoccasion variability across
    #   dosing cycles. Following the nlmixr2lib convention for paper-reported
    #   IOV when no occasion variable is part of the standalone model file,
    #   it is encoded as conventional IIV on lfdepot (see Liesenfeld 2013
    #   dabigatran for the same pattern).
    # Base    : BSV 12.0% -> omega^2 = log(0.120^2 + 1) = log(1.01440) = 0.01430
    etalvmax     ~ 0.05024                                                                              # Table 2: Vmax BSV 22.7% CV
    etalfdepot   ~ 0.13289                                                                              # Table 2: F   IOV 37.7% (recast as IIV)
    etalbaseline ~ 0.01430                                                                              # Table 2: Base BSV 12.0% CV

    # Residual error (Elassaiss-Schaap 2017 Table 2).
    # RUV_PK is reported as a 29.6% proportional term; RUV_PD is 0.209,
    # interpreted as a proportional residual on the IL-2 stim ratio
    # consistent with the log-transformed parameterization of the PD model.
    propSd       <- 0.296; label("Proportional residual error on pembrolizumab concentration (fraction)") # Table 2 RUV_PK = 29.6%
    propSd_stim_ratio <- 0.209; label("Proportional residual error on IL-2 stimulation ratio (fraction)")     # Table 2 RUV_PD = 0.209
  })
  model({
    # Individual PK parameters
    cl      <- exp(lcl)
    vc      <- exp(lvc)
    q       <- exp(lq)
    vp      <- exp(lvp)
    vmax    <- exp(lvmax + etalvmax)
    km      <- exp(lkm)

    # Individual PD parameters
    baseline <- exp(lbaseline + etalbaseline)
    imax     <- exp(limax)
    ic50     <- exp(lic50)

    # Two-compartment pembrolizumab PK with parallel linear and
    # Michaelis-Menten clearance from central (Elassaiss-Schaap 2017
    # Methods / Results "Part A2" describes the parallel linear + nonlinear
    # CL structure introduced for the final dataset). Dose in mg, volumes
    # in L => central / vc has units mg/L = ug/mL.
    Cc <- central / vc

    d/dt(central)     <- -(cl / vc) * central -
                          vmax * Cc / (km + Cc) -
                          (q / vc) * central +
                          (q / vp) * peripheral1
    d/dt(peripheral1) <-  (q / vc) * central -
                          (q / vp) * peripheral1

    # Bioavailability anchor (F = 1; IOV ~37.7% recast onto lfdepot above).
    f(central) <- exp(lfdepot + etalfdepot)

    Cc ~ prop(propSd)

    # Direct-response Imax PK/PD on the ex vivo IL-2 stimulation ratio.
    # At Cc = 0 the ratio equals baseline (~2); at Cc >> IC50 the ratio
    # asymptotes near 1, consistent with the paper's description that
    # "data from samples with high concentrations of pembrolizumab are
    # centered around a ratio value of ~1". Parameterized as a margin
    # above 1 so the asymptote at saturating drug equals 1 - (1 - Imax)
    # * (baseline - 1), which for Imax = 0.961 and Base = 2.09 gives
    # ~ 1.04, matching the qualitative description.
    inhibition  <- imax * Cc / (ic50 + Cc)
    stim_ratio  <- 1 + (baseline - 1) * (1 - inhibition)

    stim_ratio ~ prop(propSd_stim_ratio)
  })
}
