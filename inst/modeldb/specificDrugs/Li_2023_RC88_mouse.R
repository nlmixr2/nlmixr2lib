Li_2023_RC88_mouse <- function() {
  description <- paste(
    "Preclinical (mouse, Balb/c nude bearing OVCAR3-MSLN-3# xenografts).",
    "Two-compartment PK model of the mesothelin-targeting antibody-drug",
    "conjugate RC88 (conjugate measured by an anti-MMAE capture ELISA) with",
    "parallel linear and Michaelis-Menten elimination from the central",
    "compartment, after single IV doses of 0.7 and 2.1 mg/kg (Li 2023)."
  )
  reference <- paste(
    "Li Q, Wang L, Zhang J, Zhao G, Liu Z, Ma X, Jiang J.",
    "Translation of the efficacy of antibody-drug conjugates from preclinical",
    "to clinical using a semimechanistic PK/PD model: A case study with RC88.",
    "Clin Transl Sci. 2023;16(7):1232-1242. doi:10.1111/cts.13526"
  )
  vignette <- "Li_2023_RC88"
  units <- list(time = "h", dosing = "mg/kg", concentration = "ug/mL")

  covariateData <- list()

  compartmentData <- list(
    central = list(analyte = "RC88 (antibody-drug conjugate)", units = "mg/kg", specimen = "serum", verified = TRUE),
    peripheral1 = list(
      analyte = "RC88 (antibody-drug conjugate)",
      units = "mg/kg",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "mouse (Balb/c nude, OVCAR3-MSLN-3# xenograft-bearing)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    age_range = "4-5 weeks at purchase",
    sex_female_pct = 100,
    disease_state = "subcutaneous OVCAR3-MSLN-3# (mesothelin-transfected ovarian cancer) xenograft, staged at about 250 mm^3",
    dose_range = "0.7 and 2.1 mg/kg RC88, single IV dose",
    regions = "China (preclinical study, RemeGen)",
    notes = paste(
      "Destructive sampling: three animals per dose group were euthanised at",
      "each of 0.083, 4, 12, 24, 48, 72, 96, 120, 168 and 240 h post-dose",
      "(Methods, 'Mice PK studies'). The total number of animals is not",
      "stated. RC88 does not cross-react with rodent mesothelin, so no target",
      "binding is modelled in the mouse."
    )
  )

  ini({
    # Li 2023 Table 2, 'Mouse' column. Values in parentheses are bootstrap %RSE.
    lcl <- log(0.527); label("Linear clearance CL_ADC (mL/h/kg)") # Table 2: CLADC = 0.527 mL/h/kg (45.9%)
    lq <- log(1.31); label("Intercompartmental clearance CLD_ADC (mL/h/kg)") # Table 2: CLDADC = 1.31 mL/h/kg (29.4%)
    lvc <- log(67.1); label("Central volume of distribution V1 (mL/kg)") # Table 2: V1 = 67.1 mL/kg (4.97%)
    lvp <- log(49.1); label("Peripheral volume of distribution V2 (mL/kg)") # Table 2: V2 = 49.1 mL/kg (42.5%)
    lkm <- log(5.02); label("Michaelis-Menten constant Km (ug/mL)") # Table 2: Km = 5.02 ug/mL (78.5%)
    # Table 2 labels Vmax 'ug/mL/h'. Encoded as an amount rate per kg body
    # weight (ug/h/kg), the only reading under which the printed estimates
    # keep RC88 in serum for the 240 h of Figure 3a; read as a concentration
    # rate (amount rate Vmax*V1) RC88 is eliminated within a day. See the
    # vignette's Assumptions and deviations section.
    lvmax <- log(1.73); label("Maximum Michaelis-Menten elimination rate Vmax (ug/h/kg; Table 2 prints ug/mL/h)") # Table 2: Vmax = 1.73 (131%)

    # Li 2023 reports neither interindividual nor residual variability for
    # this fit (Phoenix NLME, FOCE-ELS). Carried as fixed(0) so the model
    # simulates the published typical-value profile.
    propSd <- fixed(0); label("Proportional residual error (fraction; not reported in the source)")
  })
  model({
    cl <- exp(lcl)
    q <- exp(lq)
    vc <- exp(lvc)
    vp <- exp(lvp)
    km <- exp(lkm)
    vmax <- exp(lvmax)

    # Doses are mg/kg and volumes mL/kg, so central/vc is mg/mL; x1000 gives
    # ug/mL, the unit of Km and of the ELISA concentrations.
    Cc <- 1000 * central / vc

    # Michaelis-Menten elimination is an amount rate in ug/h/kg; /1000 puts
    # it on the mg/kg scale of the states.
    mm_elim <- vmax * Cc / (km + Cc) / 1000

    d/dt(central) <- -(cl / vc) * central - (q / vc) * central + (q / vp) * peripheral1 - mm_elim
    d/dt(peripheral1) <- (q / vc) * central - (q / vp) * peripheral1

    Cc ~ prop(propSd)
  })
}
