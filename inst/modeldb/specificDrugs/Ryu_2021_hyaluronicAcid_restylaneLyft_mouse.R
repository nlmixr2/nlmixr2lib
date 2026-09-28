Ryu_2021_hyaluronicAcid_restylaneLyft_mouse <- function() {
  description <- paste(
    "Preclinical (mouse). Swelling-degradation kinetic model of the residual",
    "volume of the hyaluronic acid dermal filler Restylane Lyft with Lidocaine",
    "(Galderma; bi-phasic, BDDE-cross-linked, 1.0 x 10^6 Da HA, 20 mg/mL)",
    "after a single 100 uL subcutaneous injection into the dorsal skin of",
    "hairless mice (Ryu 2021). A depot compartment empties by first-order",
    "'swelling' (Kswell) into a subcutaneous observation compartment that is",
    "lost by first-order 'degradation' (Kdeg); the observed filler volume is",
    "the subcutaneous amount. Between-animal variability is carried on Kswell",
    "only, and the Kdeg random effect is the Kswell random effect scaled by an",
    "estimated slope. One of five filler-specific fits from the same paper.",
    sep = " "
  )
  reference <- paste(
    "Ryu H-j, Kwak S-s, Rhee C-h, Yang G-h, Yun H-y, Kang W-h. Model-Based",
    "Prediction to Evaluate Residence Time of Hyaluronic Acid Based Dermal",
    "Fillers. Pharmaceutics. 2021;13(2):133. doi:10.3390/pharmaceutics13020133.",
    "Parameter estimates from Table 1; structure from Equations 1-4 and",
    "Supplementary Text S1 (the NONMEM 7.4 control stream deposited for",
    "Neuramis VOLUME Lidocaine; the same structure was fitted to each filler).",
    sep = " "
  )
  vignette <- "Ryu_2021_hyaluronicAcid_fillers_mouse"
  units <- list(time = "day", dosing = "uL", concentration = "cm3")

  # The fit used AMT = 100 (the injected 100 uL) against observed volumes in
  # cm3 with no unit conversion, and Table 1 prints Kswell in units of
  # 10^-3 / day. Both are recovered from the paper's own Table 2 and Figure 4
  # (see the vignette); the output below is therefore on the cm3 scale when
  # the depot is dosed with amt = 100.
  compartmentData <- list(
    depot = list(
      analyte = "hyaluronic acid dermal filler",
      units = "uL",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "hyaluronic acid dermal filler",
      units = "cm3",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "mouse (SKH1-Hr hr hairless, female)",
    n_subjects = 8,
    n_studies = 1,
    age_range = "6-7 weeks at purchase",
    sex_female_pct = 100,
    disease_state = "healthy (dermal filler residence study per ISO 10993-6)",
    dose_range = "single 100 uL subcutaneous injection into the dorsal skin",
    notes = paste(
      "Methods 2.1: n = 8 mice per filler group. Filler volume was measured",
      "with a PRIMOS Lite 3D imaging system at 0, 1, 4, 7, 21 and 28 days and",
      "monthly from 2 to 18 months (lower detection limit 3 mm^3 = 0.003",
      "cm^3). No covariates were assessed.",
      sep = " "
    )
  )

  ini({
    # Table 1 prints Kswell in 'day^-1', but the deposited Supplementary Text S1
    # stream initialises KSWELL at 0.0036 and only Kswell x 10^-3 reproduces the Table 2 Tcd and
    # the Figure 4 volume peaks (see the vignette). Kswell is the depot ->
    # subcutaneous first-order transfer, i.e. the canonical absorption rate.
    lka <- log(4.74e-3); label("Swelling rate constant Kswell (1/day)") # Table 1 'Restylane Lyft with Lidocaine' 'Kswell' = 4.74 (RSE 13.1%), read as 4.74 x 10^-3 / day
    lkdeg <- log(4.24); label("Degradation rate constant Kdeg (1/day)") # Table 1 'Restylane Lyft with Lidocaine' 'Kdeg' = 4.24 (RSE 11.6%)
    kdeg_eta_scale <- 1.01; label("Slope scaling the Kswell random effect onto Kdeg (unitless)") # Table 1 'Restylane Lyft with Lidocaine' 'Slope' = 1.01 (RSE 18.3%); Equation 2 KDEG = THETA(2)*EXP(THETA(3)*ETA(1))

    etalka ~ 0.080750 # Table 1 'Restylane Lyft with Lidocaine' IIV 29% (RSE 22.3%); omega^2 = log(1 + 0.29^2)

    propSd <- 0.179; label("Proportional residual error (fraction)") # Table 1 'Restylane Lyft with Lidocaine' 'Proportional residual variability CV%' = 17.9% (RSE 11.9%)
  })
  model({
    ka <- exp(lka + etalka)
    kdeg <- exp(lkdeg + kdeg_eta_scale * etalka)

    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - kdeg * central

    # Supplementary Text S1 $ERROR: IPRED = A(2), Y = IPRED + IPRED * EPS(1);
    # the subcutaneous amount is the observed filler volume (cm3).
    central ~ prop(propSd)
  })
}
