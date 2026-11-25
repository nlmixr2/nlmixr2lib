Daryani_2016_fluorouracil_mouse <- function() {
  description <- paste(
    "Preclinical (mouse). Daryani 2016 PK-PD model of intravenous bolus and",
    "continuous-infusion 5-fluorouracil (5-FU) in CD1 nude mice bearing",
    "orthotopic (cerebral) mEP Ephb2 pediatric ependymoma, fitted",
    "sequentially to plasma PK, cerebral-microdialysis tumor extracellular",
    "fluid (tECF) PK and tumor bioluminescence. Plasma: two-compartment model",
    "with Michaelis-Menten elimination. Tumor: the unbound plasma amount",
    "(fu x central) enters a tECF compartment of fixed volume 0.001 L/kg by a",
    "first-order rate constant and returns by another; the tECF exchanges with",
    "an intracellular tumor compartment and has its own Michaelis-Menten",
    "elimination. Tumor growth: Simeoni-type exponential-to-linear growth with",
    "a power switch (exponent psi) of bioluminescence normalised to day 7",
    "after implantation, with a sigmoidal Emax kill rate driven by the tECF",
    "concentration that routes cells through three transit compartments. The",
    "exponential growth rate differs between efficacy experiment one and",
    "experiments two and three (additive covariate effect). Model is per kg",
    "of body weight: dose in mg/kg, volumes in L/kg; concentrations in uM."
  )
  reference <- paste(
    "Daryani VM, Patel YT, Tagen M, Turner DC, Carcaboso AM, Atkinson JM,",
    "Gajjar A, Gilbertson RJ, Wright KD, Stewart CF. Translational",
    "Pharmacokinetic-Pharmacodynamic Modeling and Simulation: Optimizing",
    "5-Fluorouracil Dosing in Children With Pediatric Ependymoma.",
    "CPT Pharmacometrics Syst Pharmacol. 2016;5(4):211-221.",
    "doi:10.1002/psp4.12075.",
    "Structure and the plasma/tECF values carried into the PK-PD fit are from",
    "the authors' NONMEM control streams (Supplementary Material 9)."
  )
  vignette <- "Daryani_2016_fluorouracil"
  units <- list(
    time = "h",
    dosing = "mg/kg",
    concentration = "umol/L (plasma and tECF); tumor bioluminescence is a fold change relative to day 7 after implantation (unitless)"
  )

  compartmentData <- list(
    central = list(analyte = "5-fluorouracil", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "5-fluorouracil", units = "mg/kg", specimen = "tissue", verified = TRUE),
    ecf = list(
      analyte = "5-fluorouracil (unbound, tumor extracellular fluid)",
      units = "mg/kg",
      specimen = "tumor",
      verified = TRUE
    ),
    int_tumor = list(
      analyte = "5-fluorouracil (intracellular tumor)",
      units = "mg/kg",
      specimen = "tissue",
      verified = TRUE
    ),
    cycling_cells = list(
      analyte = "proliferating tumor cells (bioluminescence fold change)",
      units = "unitless",
      specimen = "tumor",
      verified = TRUE
    ),
    damaged_cells1 = list(
      analyte = "drug-damaged tumor cells, transit 1",
      units = "unitless",
      specimen = "tumor",
      verified = TRUE
    ),
    damaged_cells2 = list(
      analyte = "drug-damaged tumor cells, transit 2",
      units = "unitless",
      specimen = "tumor",
      verified = TRUE
    ),
    damaged_cells3 = list(
      analyte = "drug-damaged tumor cells, transit 3",
      units = "unitless",
      specimen = "tumor",
      verified = TRUE
    )
  )

  covariateData <- list(
    STUDY_EXP23 = list(
      description = "1 = animal from efficacy experiment two or three of Daryani 2016; 0 = efficacy experiment one (and the PK studies)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Daryani 2016 Equations 7-8 and control stream variable E23",
        "(IF(EXPER.GT.1) E23 = 1). Adds 0.102 1/h to the exponential tumor",
        "growth rate of experiment one (0.158 1/h)."
      ),
      source_name = "EXPER"
    )
  )

  population <- list(
    species = "mouse (CD1 nude)",
    n_subjects = paste(
      "Plasma PK: 4 (i.v. bolus) + 5 (3-day s.c. infusion) non-tumor-bearing",
      "mice; microdialysis: tumor-bearing mice with plasma and tECF sampling;",
      "efficacy: 55 tumor-bearing mice (experiment one 15, two 29, three 11;",
      "Supplementary Material 1)"
    ),
    n_studies = 3L,
    disease_state = "Orthotopic (cerebral cortex) mEP Ephb2-LUC mouse model of supratentorial pediatric ependymoma",
    dose_range = paste(
      "75 mg/kg i.v. bolus (single, or weekly x 4); 37.5 mg/kg i.v. bolus",
      "weekly x 4; 75 mg/kg s.c. continuous infusion over 3 or 5 days",
      "(single, or every 3 weeks); 13 mg/kg/h s.c. infusion over 24 h"
    ),
    regions = "USA (St. Jude Children's Research Hospital)",
    notes = paste(
      "Fitted sequentially in NONMEM 7.3: plasma PK first, then the tECF",
      "model with the plasma parameters held at their estimates, then the",
      "tumor growth model with all PK parameters held. Fraction unbound in",
      "mouse plasma 0.366 by equilibrium dialysis. Bioluminescence was",
      "normalised to each animal's day-7 value, so model time zero is day 7",
      "after tumor implantation (enrolment) and cycling_cells(0) = 1."
    )
  )

  ini({
    # Plasma PK (Table 1, 'Plasma PK study'; NONMEM THETA(1)-THETA(5)).
    lvmax <- log(2040); label("Maximum plasma elimination rate (umol/h/kg)") # Table 1 Vmax,plasma,m = 2040 +/- 239
    lkm <- log(125); label("Michaelis-Menten constant for plasma elimination (umol/L)") # Table 1 Km,plasma,m = 125 +/- 22
    lvc <- log(0.962); label("Central volume of distribution (L/kg)") # Table 1 V1,m = 0.962 +/- 0.691
    lq <- log(1.67); label("Intercompartmental clearance (L/h/kg)") # Table 1 Qm = 1.67 +/- 1.03
    lvp <- log(0.332); label("Peripheral volume of distribution (L/kg)") # Table 1 V2,m = 0.332 +/- 0.50

    # Tumor ECF PK (Table 1, 'Tumor microdialysis study'). Values are the ones
    # the PK-PD control stream (Supplementary Material 9) holds constant.
    fu <- fixed(0.366); label("Fraction unbound in mouse plasma (fraction)") # Supplementary Material 9 FU = 0.366; Table 1 prints the rounded 0.37
    lvecf <- fixed(log(0.001)); label("Volume of the tumor ECF compartment (L/kg)") # Table 1 V3,m = 0.001 (unit printed as L/hr, a typo for L/kg; Results text)
    lkin_ecf <- log(0.00426); label("Rate constant from plasma (unbound) into tumor ECF (1/h)") # Supplementary Material 9 K13 = 0.00426; Table 1 K13,m = 0.0043 +/- 0.001
    lkout_ecf <- log(3.24); label("Rate constant from tumor ECF back to plasma (1/h)") # Table 1 K31,m = 3.24 +/- 0.58
    lkin_int_tumor <- log(1.82); label("Rate constant from tumor ECF into intracellular tumor (1/h)") # Table 1 K34,m = 1.82 +/- 0.30
    lkout_int_tumor <- log(0.334); label("Rate constant from intracellular tumor back to tumor ECF (1/h)") # Table 1 K43,m = 0.334 +/- 0.077
    lvmax_ecf <- log(0.00633); label("Maximum elimination rate from tumor ECF (umol/h/kg)") # Supplementary Material 9 K30 = 0.00633; Table 1 Vmax,tumor,m = 0.0063 +/- 0.004
    lkm_ecf <- log(0.0107); label("Michaelis-Menten constant for tumor ECF elimination (umol/L)") # Supplementary Material 9 KM2 = 0.0107; Table 1 prints Km,tumor,m = 0.012 +/- 0.023

    # Tumor growth and 5-FU effect (Table 1, 'Efficacy study'; Equations 1-9).
    lkgrow <- log(0.158); label("Exponential tumor growth rate, experiment one (1/h)") # Table 1 theta_Kexp,m,1 = 0.158 (0.079-0.337)
    e_study_exp23_kgrow <- 0.102; label("Additive change in exponential growth rate for experiments two and three (1/h)") # Table 1 theta_Kexp,m,2/3 = 0.102 (0.031-0.287); Equation 8
    lkgrow_lin <- log(214000); label("Linear tumor growth rate (1/h)") # Table 1 Klin,m = 214,000 (66,700-1,480,000)
    psi <- 0.0936; label("Exponent switching tumor growth from exponential to linear (unitless)") # Table 1 Psi_m = 0.0936 (0.085-0.104)
    lkmax <- log(0.658); label("Maximum 5-FU kill rate of proliferating tumor cells (1/h)") # Table 1 Kmax,m = 0.658 (0.316-1.134)
    lec50 <- log(2.12); label("Tumor ECF 5-FU concentration giving half-maximal kill rate (umol/L)") # Table 1 IC50,m = 2.12 (0.28-5.84)
    hill <- 1.22; label("Hill coefficient of the kill-rate Emax model (unitless)") # Table 1 Hm = 1.22 (0.64-1.89)
    lktr <- log(0.946); label("Transit rate constant through the damaged-cell compartments (1/h)") # Table 1 Kdel,m = 0.946 (0.124-1.950)

    # IIV: exponential, Table 1 IIV (% CV) converted as omega^2 = log(1 + CV^2).
    etalvmax ~ 0.01477 # Table 1 IIV Vmax,plasma,m 12.2% CV
    etalvc ~ 0.04769 # Table 1 IIV V1,m 22.1% CV
    etalq ~ 0.09232 # Table 1 IIV Qm 31.1% CV
    etalkin_ecf ~ 0.006066 # Table 1 IIV K13,m 7.8% CV
    etalkout_ecf ~ 0.0659 # Table 1 IIV K31,m 26.1% CV
    etalkgrow ~ 0.1526 # Table 1 IIV Kexp 40.6% CV (one eta on the summed rate, Equation 9)
    etalec50 ~ 0.08839 # Table 1 IIV IC50,m 30.4% CV

    # Residual error: proportional (Table 1). Plasma uses the value from the
    # joint plasma + tECF fit; the plasma-only fit gave 30.7% CV.
    propSd <- 0.359; label("Proportional residual error, plasma 5-FU (fraction)") # Table 1 sigma_plasma,prop 35.9% CV (tumor microdialysis block)
    propSd_Cecf <- 0.279; label("Proportional residual error, tumor ECF 5-FU (fraction)") # Table 1 sigma_tECF,prop 27.9% CV
    propSd_tumor_size <- 0.37; label("Proportional residual error, tumor bioluminescence fold change (fraction)") # Table 1 sigma_efficacy,prop 37.0% CV
  })

  model({
    # 5-FU molecular weight (g/mol): converts the mg/kg amounts to the umol/kg
    # basis of the published Vmax / Km values; Supplementary Table 1 converts
    # 751 mg/h/m^2 to 5773 umol/h/m^2 with the same factor.
    mw <- 130.08

    vmax <- exp(lvmax + etalvmax)
    km <- exp(lkm)
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp)
    vecf <- exp(lvecf)
    kin_ecf <- exp(lkin_ecf + etalkin_ecf)
    kout_ecf <- exp(lkout_ecf + etalkout_ecf)
    kin_int_tumor <- exp(lkin_int_tumor)
    kout_int_tumor <- exp(lkout_int_tumor)
    vmax_ecf <- exp(lvmax_ecf)
    km_ecf <- exp(lkm_ecf)

    # Equations 7-9: additive experiment effect, then one exponential eta.
    kgrow <- (exp(lkgrow) + e_study_exp23_kgrow * STUDY_EXP23) * exp(etalkgrow)
    kgrow_lin <- exp(lkgrow_lin)
    kmax <- exp(lkmax)
    ec50 <- exp(lec50 + etalec50)
    ktr <- exp(lktr)

    k12 <- q / vc
    k21 <- q / vp

    # Concentrations (umol/L) from mg/kg amounts.
    Cc <- central / vc * 1000 / mw
    Cecf <- ecf / vecf * 1000 / mw

    # Kill rate (Equation 6); the control stream floors the tECF
    # concentration at zero before the Hill power.
    Cecf_pd <- Cecf
    if (Cecf_pd < 0) {
      Cecf_pd <- 0
    }
    kkill <- kmax * Cecf_pd^hill / (ec50^hill + Cecf_pd^hill)

    # Exponential-to-linear growth (Equation 1).
    growth_rate <- kgrow * cycling_cells / (1 + (kgrow / kgrow_lin * cycling_cells)^psi)^(1 / psi)

    # Michaelis-Menten eliminations converted from umol/h/kg to mg/h/kg.
    d/dt(central) <- -vmax * mw / 1000 * Cc / (km + Cc) - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(ecf) <- kin_ecf * fu * central - kout_ecf * ecf - kin_int_tumor * ecf + kout_int_tumor * int_tumor -
      vmax_ecf * mw / 1000 * Cecf / (km_ecf + Cecf)
    d/dt(int_tumor) <- kin_int_tumor * ecf - kout_int_tumor * int_tumor

    # Tumor growth and transit (Equations 1-4); W(0) = 1 at day 7.
    d/dt(cycling_cells) <- growth_rate - kkill * cycling_cells
    d/dt(damaged_cells1) <- kkill * cycling_cells - ktr * damaged_cells1
    d/dt(damaged_cells2) <- ktr * damaged_cells1 - ktr * damaged_cells2
    d/dt(damaged_cells3) <- ktr * damaged_cells2 - ktr * damaged_cells3
    cycling_cells(0) <- 1

    # Fold change in bioluminescence relative to day 7 (Equation 5).
    tumor_size <- cycling_cells + damaged_cells1 + damaged_cells2 + damaged_cells3

    Cc ~ prop(propSd)
    Cecf ~ prop(propSd_Cecf)
    tumor_size ~ prop(propSd_tumor_size)
  })
}
