Singh_2016_stamulumab_mouse <- function() {
  description <- "Preclinical (mouse, C57/SCID). One-compartment population PK model with parallel linear and Michaelis-Menten elimination for stamulumab (MYO-029), a human IgG1 anti-myostatin (GDF-8) monoclonal antibody, after IV dosing in C57/SCID mice (Singh 2016). Parameters are reported per kg body weight and scale linearly with WT. Typical-value model: no between-animal variability or residual error is reported for the mouse fit."
  reference <- "Singh P, Rong H, Gordi T, Bosley J, Bhattacharya I. Translational Pharmacokinetic/Pharmacodynamic Analysis of MYO-029 Antibody for Muscular Dystrophy. Clin Transl Sci. 2016;9(6):302-310. doi:10.1111/cts.12420. PMID 27700008. PMCID PMC5351001. MYO-029 is the development code of the antibody assigned the INN stamulumab."
  vignette <- "Singh_2016_stamulumab"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Singh 2016 Table 1 reports every mouse parameter per kg body weight",
        "(mL/h/kg, mL/kg, ug/h/kg), so CL, V and Vmax scale linearly with WT",
        "and doses are given as mg/kg * WT. The paper's conversion",
        "'189 ug/h/kg (or 0.11 mg/day)' implies a 0.025 kg mouse."
      ),
      source_name = "per-kg normalisation (no WT column printed)"
    )
  )

  compartmentData <- list(
    central = list(analyte = "stamulumab", units = "mg", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "mouse (C57/SCID)",
    n_subjects = NA,
    n_studies = 1,
    age_range = NA,
    weight_range = "not reported; about 0.025 kg implied by the paper's Vmax unit conversion",
    sex_female_pct = NA,
    race_ethnicity = NA,
    disease_state = "Healthy C57/SCID mice (PK study); the companion efficacy study used female SCID mice",
    dose_range = "1, 5, 20 and 100 mg/kg IV (weekly for 4 weeks)",
    regions = NA,
    notes = paste(
      "Singh 2016 Methods 'Mouse pharmacokinetic studies': 1-2 serum samples",
      "per mouse at 3, 6, 48, 144, 312 and 648 h, so the data were sparse and",
      "no inter-animal variability was estimated (Table 1 mouse column has",
      "no IIV or residual error). Serum MYO-029 by GDF-8-capture ELISA",
      "(quantitation range 10-80 ng/mL in 100% mouse serum)."
    )
  )

  ini({
    # Singh 2016 Table 1, 'Mouse' column (no RSE printed for the mouse fit).
    lcl <- log(0.245); label("Linear clearance per kg body weight (mL/h/kg)") # Table 1, CL 'Clearance from central compartment' mouse = 0.245 mL/h/kg
    lvc <- log(103); label("Central volume per kg body weight (mL/kg)") # Table 1, V1 'Volume of central compartment' mouse = 103 mL/kg
    lvmax <- log(189); label("Maximum rate of the nonlinear elimination route per kg body weight (ug/h/kg)") # Table 1, Vmax mouse = 189 ug/h/kg; text '189 ug/h/kg (or 0.11 mg/day)'
    # Table 1 prints Km = 2.67 'ug/mL'. Taken as printed, Vmax/Km = 70.8
    # mL/h/kg would eliminate a 1 mg/kg dose within hours, contradicting the
    # fitted lines of Figure 1a (about 1e4 ng/mL falling to 1.4e3 ng/mL over
    # 650 h). Read as 2.67 mg/mL (2670 ug/mL), the model reproduces the four
    # Figure 1a fitted lines within 10% on average per dose; see the
    # vignette errata.
    lkm <- log(2.67 * 1000); label("Michaelis-Menten constant of the nonlinear route (mg/L)") # Table 1, Km 'Apparent binding constant' mouse = 2.67, printed 'ug/mL', read as mg/mL (= 2670 mg/L)

    # Residual error not reported for the mouse fit; fixed to zero.
    propSd <- fixed(0); label("Proportional residual error (fraction)") # not reported in Table 1 (mouse column '-')
  })

  model({
    # Per-kg parameters scaled to the individual: mL -> L.
    cl <- exp(lcl) * WT / 1000
    vc <- exp(lvc) * WT / 1000
    # Vmax ug/h/kg -> mg/h for the individual
    vmax <- exp(lvmax) * WT / 1000
    km <- exp(lkm)

    Cc <- central / vc
    d/dt(central) <- -cl / vc * central - vmax * Cc / (km + Cc)

    Cc ~ prop(propSd)
  })
}
