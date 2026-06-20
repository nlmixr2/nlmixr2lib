Langdon_2010_PF00821385_dog <- function() {
  description <- paste(
    "Preclinical (beagle dog). Translational popPK-PD model for PF-00821385,",
    "a Pfizer HIV-1 gp120 cell-fusion inhibitor candidate (molecular weight",
    "440.49 g/mol) studied in conscious freely-moving Beagle dogs. PK is a",
    "one-compartment disposition model with first-order oral absorption; PD",
    "describes heart rate as the sum of (a) a typical-value baseline HR with",
    "log-normal inter-subject variability, (b) a 24-h cosine circadian rhythm",
    "with typical-value amplitude and a log-normal inter-subject variable",
    "peak time, and (c) a linear drug effect on free plasma concentration",
    "with no IIV. The PD-slope SLOPE = 1.76 bpm per micromolar free drug is",
    "from Langdon 2010 Table 1; plasma unbound fraction fu = 0.64 is FIXED",
    "via back-calculation from the published unbound vs total Cmax ratio at",
    "20 mg/kg oral (paper Introduction); see vignette Errata. PK and PD were",
    "fit sequentially in NONMEM VI using FOCE INTER with individual Bayesian",
    "post hoc PK estimates serving as input to the PD model.",
    sep = " "
  )
  reference <- paste(
    "Langdon G, Davis JD, McFadyen LM, Dewhurst M, Brunton NS, Rawal JK,",
    "Van der Graaf PH, Benson N.",
    "Translational pharmacokinetic-pharmacodynamic modelling; application",
    "to cardiovascular safety data for PF-00821385, a novel HIV agent.",
    "Br J Clin Pharmacol. 2010 Apr;69(4):335-345.",
    "doi:10.1111/j.1365-2125.2009.03594.x.",
    sep = " "
  )
  vignette <- "Langdon_2010_PF00821385"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ng/mL",
    dosing_notes  = paste(
      "Oral gavage. Doses in the source paper are mg/kg; convert to total",
      "mg by multiplying by individual body weight (cohort 16.3 - 16.4 kg).",
      sep = " "
    )
  )

  covariateData <- list()

  population <- list(
    species        = "beagle dog (Pfizer laboratory colony, conscious freely-moving)",
    n_subjects     = 18L,
    n_studies      = 1L,
    weight_range   = "16.3 - 16.4 kg",
    weight_median  = "approximately 16.35 kg",
    sex_female_pct = NA_real_,
    disease_state  = "Healthy laboratory beagle (no induced disease).",
    dose_range     = paste(
      "Oral gavage single doses 0.5, 1.5, 5, 15, 20, 40, 100, and 120 mg/kg.",
      "PD experiments (4 dogs): vehicle plus 1.5, 5, and 15 mg/kg per animal",
      "on separate occasions. Subsequent PK experiment in the same 4 PD",
      "animals at 5 mg/kg oral. Satellite PK experiments in 14 additional",
      "dogs at 0.5 (n = 2), 20 (n = 1), 40 (n = 1), 100 (n = 2), and 120",
      "(n = 2) mg/kg oral; 2 of those animals were also given 0.5 mg/kg IV.",
      "88 oral PD/PK plasma concentrations and 22 IV PK plasma",
      "concentrations contributed to the canine popPK fit.",
      sep = " "
    ),
    regions        = "United Kingdom (Pfizer Sandwich laboratory colony).",
    notes          = paste(
      "Cardiovascular telemetry data (heart rate, blood pressure) acquired",
      "as 1-minute means via implanted Konigsberg pressure transducers and",
      "subcutaneous ECG electrodes; 1-minute means were preprocessed into",
      "15-minute bins yielding 96 HR observations per animal per dose.",
      "Periods of feeding were excluded. Dosing was at 09:00 in the morning.",
      "NONMEM VI with FOCE INTER; PK and PD were fit sequentially (individual",
      "Bayesian post hoc PK served as input to the PD model). Diagnostic",
      "plots in Xpose (S-Plus); see Langdon 2010 Methods 'Canine",
      "pharmacokinetic-pharmacodynamic analysis'.",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # PK structural parameters -- Langdon 2010 Table 1 (canine final NONMEM
    # estimate; SE and 95% CI as reported).
    # =====================================================================
    lfdepot <- log(1.01); label("Oral bioavailability (estimated; near-unity)")               # Table 1: F = 1.01 (SE 0.003; 95% CI 1.00 - 1.02)
    lka     <- log(2.42); label("First-order oral absorption rate constant ka (1/h)")         # Table 1: ka = 2.42 1/h (SE 0.520; 95% CI 1.40 - 3.44)
    lcl     <- log(1.02); label("Apparent clearance CL/F (L/h)")                              # Table 1: CL/F = 1.02 L/h (SE 0.135; 95% CI 0.755 - 1.28)
    lvc     <- log(7.05); label("Apparent central volume V/F (L)")                            # Table 1: V/F = 7.05 L (SE 0.640; 95% CI 5.80 - 8.30)

    # =====================================================================
    # PD structural parameters -- Langdon 2010 Table 1 (canine final NONMEM
    # estimate). Equation 1 (placebo HR cosine) and Equation 2 (linear free-
    # drug effect) on page 338.
    # =====================================================================
    lbase  <- log(68.2); label("Baseline heart rate (bpm)")                                   # Table 1: BASE = 68.2 bpm (SE 1.01; 95% CI 66.2 - 70.2)
    lamp   <- log(12.4); label("Circadian-rhythm amplitude on HR (bpm)")                      # Table 1: AMP = 12.4 bpm (SE 0.706; 95% CI 11.01 - 13.8)
    lpeak  <- log(10.4); label("Time of circadian peak post morning dose (h)")                # Table 1: PEAK = 10.4 h (SE 0.246; 95% CI 9.91 - 10.88)
    lslope <- log(1.76); label("Linear slope of HR on free plasma concentration (bpm per uM free)") # Table 1: SLOPE = 1.76 bpm/uM-free (SE 0.29; 95% CI 1.17 - 2.35)

    # Reference unbound fraction. FIXED via back-calculation from the
    # Introduction's 'unbound Cmax = 56 uM' at 20 mg/kg oral, divided by
    # the typical-value total Cmax = 87.8 uM computed from the Table 1 PK
    # parameters (D = 326 mg in a 16.3 kg dog; ka = 2.42, kel = CL/V =
    # 0.145, F = 1.01; Cmax = D * F * ka / (V * (ka - kel)) *
    # (exp(-kel * tmax) - exp(-ka * tmax)) = 38.67 mg/L = 87.8 uM total).
    # fu = 56 / 87.8 = 0.64. The paper does not state fu directly; see
    # vignette Errata for the full derivation. PF-00821385 was discontinued
    # at first-in-human so author correspondence is unlikely to recover the
    # original NONMEM-supplied value.
    lfu <- fixed(log(0.64)); label("Plasma unbound fraction (FIXED; back-calculated)")        # vignette Errata (not in source); back-calculation from paper Introduction

    # =====================================================================
    # IIV. Per Methods paragraph 3 of Langdon 2010 ('The interindividual
    # variability (IIV) of PK and PD parameters was modelled using
    # multiplicative exponential random effects'), every IIV term below is
    # encoded as a log-normal (additive on log-scale) random effect on the
    # log-transformed structural parameter; the reported variances are
    # omega^2 directly.
    #
    # PK IIV (Methods + Table 1: 'IIV was included on the following
    # structural parameters: absorption rate constant (ka), volume of
    # distribution of the central compartment (V/F), and clearance from
    # the central compartment (CL/F)').
    # =====================================================================
    etalka ~ 0.404   # Table 1: Variance for ka = 0.404 (95% CI 0.128 - 0.680); CV approx 70.7%
    etalcl ~ 0.311   # Table 1: Variance for CL/F = 0.311 (95% CI 0.111 - 0.511); CV approx 60.4%
    etalvc ~ 0.136   # Table 1: Variance for V/F = 0.136 (95% CI 0.003 - 0.269); CV approx 38.8%

    # PD IIV. The main text reads 'IIV was included on the following
    # structural parameters: baseline HR and amplitude of the circadian
    # effect' but Table 1 reports variance for BASE and PEAK (not AMP).
    # Table 1 is followed here; the discrepancy is documented in the
    # vignette Errata.
    etalbase ~ 0.0057  # Table 1: Variance for BASE = 0.0057 (95% CI 0.0024 - 0.0090); CV approx 7.6%, matches text 'approximately 7.5%'
    etalpeak ~ 0.0134  # Table 1: Variance for PEAK = 0.0134 (95% CI 0.0044 - 0.0224); CV approx 11.6%

    # =====================================================================
    # Residual error.
    # PK: combined proportional + additive (Methods + Table 1).
    # PD: additive only on HR (Methods + Table 1).
    # =====================================================================
    propSd   <- 0.181; label("Proportional PK residual error (fraction)")                     # Table 1: PK proportional residual error = 0.181 (SE 0.022; 95% CI 0.138 - 0.224)
    addSd    <- 9.43;  label("Additive PK residual error (ng/mL)")                            # Table 1: PK additive residual error = 9.43 ng/mL (SE 3.53; 95% CI 2.51 - 16.3)
    addSd_HR <- 14.1;  label("Additive PD residual error on HR (bpm)")                        # Table 1: PD additive residual error = 14.1 bpm (SE 4.17; 95% CI 5.93 - 22.3)
  })

  model({
    # Constant: PF-00821385 molecular weight (paper Introduction,
    # 'PF-00821385, molecular weight = 440.49 Da'). Used to convert
    # plasma total ng/mL into micromolar for the free-drug PD argument.
    mw <- 440.49  # g/mol

    # Individual structural parameters.
    fdepot <- exp(lfdepot)
    ka     <- exp(lka + etalka)
    cl     <- exp(lcl + etalcl)
    vc     <- exp(lvc + etalvc)
    fu     <- exp(lfu)
    kel    <- cl / vc

    # Individual PD parameters. BASE_i = BASE * exp(eta_base) etc per
    # multiplicative exponential IIV; encoded via log-scale parents.
    base_i  <- exp(lbase + etalbase)
    amp_i   <- exp(lamp)
    peak_i  <- exp(lpeak + etalpeak)
    slope_i <- exp(lslope)

    # PK ODE -- 1-compartment disposition with first-order absorption.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Bioavailability anchor (typical value approximately 1.01).
    f(depot) <- fdepot

    # Plasma concentration: central is in mg, vc in L so central/vc is in
    # mg/L; multiply by 1000 to express Cc in ng/mL (paper unit).
    Cc <- 1000 * central / vc

    # Free plasma concentration in micromolar (the argument of the SLOPE
    # in Equation 2 of Langdon 2010). Cc_free_uM = fu * Cc_ng_per_mL / MW
    # (g/mol) (1 ng/mL / (g/mol) = 1 nmol/L; the fu factor and the implicit
    # 1000 ng/mL -> mg/L cancel so the conversion is Cc / MW in uM).
    Cc_free_uM <- fu * Cc / mw

    # Equation 1: placebo HR cosine circadian rhythm.
    HR_placebo <- base_i + amp_i * cos(2 * pi / 24 * (t - peak_i))

    # Equation 2: HR = HR_placebo + SLOPE * Cc_free.
    HR <- HR_placebo + slope_i * Cc_free_uM

    # Multi-output residual error.
    Cc ~ prop(propSd) + add(addSd)
    HR ~ add(addSd_HR)
  })
}
