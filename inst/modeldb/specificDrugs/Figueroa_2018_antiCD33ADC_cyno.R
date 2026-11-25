Figueroa_2018_antiCD33ADC_cyno <- function() {
  description <- "Preclinical (cynomolgus monkey). Two-compartment PK model with linear non-specific clearance and Michaelis-Menten target-mediated (CD33) elimination from the central compartment, plus a receptor-occupancy output RO = Cc / (KM + Cc), fitted by a naive-pooled approach to total-antibody serum concentrations and CD33 receptor occupancy on circulating monocytes and granulocytes after single IV doses of an unconjugated anti-CD33 IgG1 and used to predict the PK of the corresponding site-specific PBD anti-CD33 antibody-drug conjugate (Figueroa 2018)."
  reference <- "Figueroa I, Leipold D, Leong S, Zheng B, Triguero-Carrasco M, Fourie-O'Donohue A, Kozak KR, Xu K, Schutten M, Wang H, Polson AG, Kamath AV. Prediction of non-linear pharmacokinetics in humans of an antibody-drug conjugate (ADC) when evaluation of higher doses in animals is limited by tolerability: Case study with an anti-CD33 ADC. MAbs. 2018;10(5):738-750. doi:10.1080/19420862.2018.1465160"
  vignette <- "Figueroa_2018_antiCD33ADC"
  units <- list(time = "day", dosing = "ug/kg", concentration = "ug/mL")

  covariateData <- list()

  compartmentData <- list(
    central = list(
      analyte = "anti-CD33 antibody (total antibody)",
      units = "ug/kg",
      specimen = "serum",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "anti-CD33 antibody (total antibody)",
      units = "ug/kg",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "cynomolgus monkey",
    n_subjects = 6L,
    n_studies = 1L,
    disease_state = "Healthy (non-diseased) cynomolgus monkeys in a single-dose PKPD study of the unconjugated anti-CD33 mAb.",
    dose_range = "Single IV bolus of the unconjugated anti-CD33 mAb at 0.5 mg/kg (n = 3) or 15 mg/kg (n = 3) nominal; measured dosing-solution concentrations were 64% and 86% of nominal, so the fit used actual doses of 0.32 and 12.9 mg/kg.",
    notes = "Fitted simultaneously to total-antibody (TAB) serum concentrations and CD33 receptor occupancy on monocytes and granulocytes using a naive-pooled approach (no between-animal variability) in SimBiology (MATLAB 2016a); proportional weighting for PK and constant weighting for RO. ADA were detected in 5 of 6 animals with no apparent effect on TAB. The model was then used, without refitting, to predict the TAB PK of the conjugated anti-CD33 ADC in a separate toxicology study (0.1, 0.2, 0.4 and 1 mg/kg, n = 4, 4, 1, 1). Figueroa 2018 Results and Materials and methods."
  )

  ini({
    # All values: Figueroa 2018 Table 3, 'Fitted Monkey' column. Parameters
    # are per kg of body weight (mL/kg, mL/day/kg, ug/day/kg), so a dose in
    # ug/kg gives a concentration in ug/mL.
    lcl <- log(4.92); label("Non-specific (linear) clearance (mL/day/kg)") # Table 3, CL = 4.92 mL/kg day (SE 0.635)
    lq <- log(22.7); label("Distributional clearance (mL/day/kg)") # Table 3, CLD = 22.7 mL/kg day (SE 8.78)
    lvc <- log(36.7); label("Central volume of distribution (mL/kg)") # Table 3, V1 = 36.7 mL/kg (SE 1.53)
    lvp <- log(18.1); label("Peripheral volume of distribution (mL/kg)") # Table 3, V2 = 18.1 mL/kg (SE 2.75)
    lvmax <- log(22.4); label("Maximum target-mediated elimination rate (ug/day/kg)") # Table 3, Vmax = 22.4 ug/day/kg (SE 9.7)
    # The paper estimated alpha = log10(KM) = -1.22 (SE 0.316); Table 3 prints
    # the derived KM = 0.0602 ug/mL = 10^-1.22, which fixes the log base at 10.
    lkm <- log(10^-1.22); label("Michaelis-Menten constant of target-mediated elimination (ug/mL)") # Table 3, alpha = log(KM) = -1.22, KM = 0.0602 ug/mL

    # Residual-error magnitudes were not reported (proportional weighting for
    # PK and constant weighting for RO, Materials and methods 'PKPD modeling').
    propSd <- fixed(0); label("Proportional residual error on Cc (fraction); not reported")
    addSd_RO <- fixed(0); label("Additive residual error on RO (%); not reported")
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
