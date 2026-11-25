Figueroa_2018_antiCD33ADC_human <- function() {
  description <- "Human projection (scaled from cynomolgus monkey, not fitted to human data). Two-compartment PK model with linear non-specific clearance and Michaelis-Menten target-mediated (CD33) elimination from the central compartment, plus a receptor-occupancy output RO = Cc / (KM + Cc), for a site-specific PBD anti-CD33 antibody-drug conjugate in patients with acute myeloid leukemia. Linear parameters and KM carried from the monkey fit with a per-kg scaling exponent of 1; Vmax is the monkey value and the source explores 1- to 10-fold higher target capacity (Figueroa 2018)."
  reference <- "Figueroa I, Leipold D, Leong S, Zheng B, Triguero-Carrasco M, Fourie-O'Donohue A, Kozak KR, Xu K, Schutten M, Wang H, Polson AG, Kamath AV. Prediction of non-linear pharmacokinetics in humans of an antibody-drug conjugate (ADC) when evaluation of higher doses in animals is limited by tolerability: Case study with an anti-CD33 ADC. MAbs. 2018;10(5):738-750. doi:10.1080/19420862.2018.1465160"
  vignette <- "Figueroa_2018_antiCD33ADC"
  units <- list(time = "day", dosing = "ug/kg", concentration = "ug/mL")

  covariateData <- list()

  compartmentData <- list(
    central = list(analyte = "anti-CD33 ADC (total antibody)", units = "ug/kg", specimen = "serum", verified = TRUE),
    peripheral1 = list(
      analyte = "anti-CD33 ADC (total antibody)",
      units = "ug/kg",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 0L,
    disease_state = "Projected for patients with acute myeloid leukemia (AML) / CD33-positive myeloid malignancies; no human data were used to build the model.",
    dose_range = "Projections shown for single IV doses of 0.25, 2.4 and 21 ug/kg (about 17 ug to 1.5 mg for a 70 kg patient), chosen for about 10%, 50% and 90% RO at Cmax; clearance-versus-dose curves span 100 to 10000 ug/kg.",
    notes = "Scaled from the cynomolgus monkey fit (Figueroa_2018_antiCD33ADC_cyno). CL, CLD, V1 and V2 are carried per kg (weight-based scaling exponent of 1, from historical in-house monkey-to-human ADC data). KM is carried unchanged because the anti-CD33 mAb has comparable affinity for human and monkey CD33 (Table 1). Vmax, which reflects whole-body target capacity, was explored over 1- to 10-fold the monkey value to reflect a possibly larger and more variable CD33 burden in AML; the packaged value is the 1-fold case. Figueroa 2018 'Model scale-up from monkey to human' and 'Prediction of human PKPD'."
  )

  ini({
    # All values: Figueroa 2018 Table 3, 'Predicted Human' column. They are
    # projections, not estimates, so every value is fixed. Per kg of body
    # weight: a dose in ug/kg gives a concentration in ug/mL.
    lcl <- fixed(log(4.92)); label("Non-specific (linear) clearance (mL/day/kg)") # Table 3 Predicted Human, CL = 4.92 mL/kg day
    lq <- fixed(log(22.7)); label("Distributional clearance (mL/day/kg)") # Table 3 Predicted Human, CLD = 22.7 mL/kg day
    lvc <- fixed(log(36.7)); label("Central volume of distribution (mL/kg)") # Table 3 Predicted Human, V1 = 36.7 mL/kg
    lvp <- fixed(log(18.1)); label("Peripheral volume of distribution (mL/kg)") # Table 3 Predicted Human, V2 = 18.1 mL/kg
    # Table 3 gives a range of 22.4 to 224 ug/day/kg (1- to 10-fold the
    # monkey value); 22.4 is the lower (1-fold) end. Set lvmax to log(67.2)
    # or log(224) for the 3- and 10-fold scenarios of Table 4 and Figure 5.
    lvmax <- fixed(log(22.4)); label("Maximum target-mediated elimination rate (ug/day/kg), 1-fold monkey value") # Table 3 Predicted Human, Vmax = 22.4 to 224 ug/day/kg
    lkm <- fixed(log(10^-1.22)); label("Michaelis-Menten constant of target-mediated elimination (ug/mL)") # Table 3 Predicted Human, alpha = log10(KM) = -1.22, KM = 0.0602 ug/mL

    propSd <- fixed(0); label("Proportional residual error on Cc (fraction); not applicable to a projection")
    addSd_RO <- fixed(0); label("Additive residual error on RO (%); not applicable to a projection")
  })
  model({
    cl <- exp(lcl)
    q <- exp(lq)
    vc <- exp(lvc)
    vp <- exp(lvp)
    vmax <- exp(lvmax)
    km <- exp(lkm)

    Cc <- central / vc

    # Equations 1-2 written in amounts (ug/kg): linear CL, distributional CLD,
    # and target-mediated Vmax * C1 / (KM + C1).
    d/dt(central) <- -cl * Cc - q * Cc + q * peripheral1 / vp - vmax * Cc / (km + Cc)
    d/dt(peripheral1) <- q * Cc - q * peripheral1 / vp

    # Equation 3: receptor occupancy (%); receptor availability RA = 100 - RO.
    RO <- 100 * Cc / (km + Cc)
    RA <- 100 - RO

    Cc ~ prop(propSd)
    RO ~ add(addSd_RO)
  })
}
