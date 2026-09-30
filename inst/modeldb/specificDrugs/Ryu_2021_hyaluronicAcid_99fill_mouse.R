Ryu_2021_hyaluronicAcid_99fill_mouse <- function() {
  description <- paste(
    "Preclinical (mouse). Swelling-degradation kinetic model of the residual",
    "volume of the hyaluronic acid dermal filler 99 fill (Hanmi",
    "Pharmaceutical; mono-phasic, other properties confidential per Table S1)",
    "after a single 100 uL subcutaneous injection into the dorsal skin of",
    "hairless mice (Ryu 2021). A depot compartment empties by first-order",
    "'swelling' (Kswell) into a subcutaneous observation compartment that is",
    "lost by first-order 'degradation' (Kdeg); the observed filler volume is",
    "the subcutaneous amount. Between-animal variability is carried on Kswell",
    "only; the slope scaling that random effect onto Kdeg was fixed to 0 for",
    "this filler. One of five filler-specific fits from the same paper.",
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
    lka <- log(2.2e-3); label("Swelling rate constant Kswell (1/day)") # Table 1 '99 fill' 'Kswell' = 2.2 (RSE 12.6%), read as 2.2 x 10^-3 / day
    lkdeg <- log(1.45); label("Degradation rate constant Kdeg (1/day)") # Table 1 '99 fill' 'Kdeg' = 1.45 (RSE 8.8%)
    kdeg_eta_scale <- fixed(0); label("Slope scaling the Kswell random effect onto Kdeg (unitless)") # Table 1 '99 fill' 'Slope' = 0, footnote: 'The slope for 99 fill was fixed as 0 because of the estimation tendency to zero value'; Equation 2

    etalka ~ 0.031538 # Table 1 '99 fill' IIV 17.9% (RSE 20.3%); omega^2 = log(1 + 0.179^2)

    propSd <- 0.163; label("Proportional residual error (fraction)") # Table 1 '99 fill' 'Proportional residual variability CV%' = 16.3% (RSE 13.1%)
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
