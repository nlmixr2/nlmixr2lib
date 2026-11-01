Michalon_2025_cliramitug_rat <- function() {
  description <- "Preclinical (rat). Two-compartment population PK model with linear elimination for cliramitug (NI006, ALXN2220), a human IgG1 anti-ATTR amyloid depleter antibody, after intravenous dosing in rats (3 and 30 mg/kg single dose, 300 mg/kg weekly x3). These rat estimates were allometrically scaled to the human PK layer of the companion QSP model Michalon_2025_cliramitug_human_qsp."
  reference <- "Michalon A, Renaud L, Machacek M, Cortijo C, Udata C, Mercuri MF, Buller F, Hock C, Nitsch RM, Kahr PC, Grimm J. Prediction of Cardiac ATTR Depletion by NI006 (ALXN2220) Using Mechanistic PK/PD Modeling. Clin Pharmacol Ther. 2025;117(1):261. doi:10.1002/cpt.3455. Rat parameters from Table 1; model structure from Supplementary Data S1 'Rat PK' Mlxtran code."
  vignette <- "Michalon_2025_cliramitug"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list()

  compartmentData <- list(
    central = list(
      analyte = "cliramitug (NI006)",
      units = "mg",
      specimen = "serum",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "cliramitug (NI006)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "rat",
    n_subjects = 36L,
    n_studies = 1L,
    weight_range = "Not reported; 0.3 kg used as the rat reference body weight for allometric scaling",
    sex_female_pct = 50,
    disease_state = "Healthy (cliramitug does not bind native transthyretin; rats carry no ATTR amyloid)",
    dose_range = "3 or 30 mg/kg IV single dose, or 300 mg/kg IV weekly for 3 doses (12 animals per group, 6 male and 6 female)",
    regions = "Not reported (GLP toxicology / PK study)",
    notes = "Serum sampled from 5 min to 28 days post dose with a GLP-qualified assay; in the 300 mg/kg group samples were taken after the first and third doses only (Methods 'Rat PK study'; Figure 1 legend). Fit in Monolix by SAEM. Strain not reported."
  )

  ini({
    # Table 1, 'Rat parameter values' column (between-animal variability as %CV in brackets).
    lcl <- log(0.00245); label("Clearance, Cl (L/day)") # Table 1, Cl = 0.00245 L/day (21%)
    lvc <- log(0.0144); label("Central volume of distribution, V1 (L)") # Table 1, V1 = 0.0144 L (10%)
    lvp <- log(0.0189); label("Peripheral volume of distribution, V2 (L)") # Table 1, V2 = 0.0189 L (15%)
    lq <- log(0.0178); label("Intercompartmental clearance, Q (L/day)") # Table 1, Q = 0.0178 L/day (18%)

    # Monolix log-normal random effects; %CV converted with omega^2 = log(1 + CV^2).
    etalcl ~ 0.04316 # Table 1, Cl between-animal variability 21% CV
    etalvc ~ 0.00995 # Table 1, V1 between-animal variability 10% CV
    etalvp ~ 0.02225 # Table 1, V2 between-animal variability 15% CV
    etalq ~ 0.03189 # Table 1, Q between-animal variability 18% CV

    # The residual-error model and its estimate are not reported anywhere in
    # the paper or supplement; held at zero (see the vignette's Assumptions).
    propSd <- fixed(0); label("Proportional residual error (fraction); not reported") # Not reported in Table 1 or Data S1
  })

  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q <- exp(lq + etalq)

    # Data S1 'Rat PK': k = Cl/V1, k12 = Q/V1, k21 = Q/V2
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # mg/L is numerically ug/mL (Data S1 converts mg/mL to ug/mL for output).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
