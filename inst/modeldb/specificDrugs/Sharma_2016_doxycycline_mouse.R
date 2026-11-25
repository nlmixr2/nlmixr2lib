Sharma_2016_doxycycline_mouse <- function() {
  description <- "Preclinical (mouse). One-compartment oral PK model with first-order absorption, per kg body weight, for doxycycline in BALB/c SCID mice, linked through a mass-transfer effect compartment (influx k1e from the central amount, efflux ke1) to an Emax model of the percentage depletion of Wolbachia endosymbionts in adult female Brugia malayi filarial worms. EC50 is fixed to the in vitro potency against Wolbachia-infected C6/36 insect cells and Emax to 99.999 %; the effect-site volume is a fixed fraction of V/F. Companion of Sharma_2016_minocycline_mouse (Sharma 2016)"
  reference <- paste(
    "Sharma R, Al Jayoussi G, Tyrer HE, Gamble J, Hayward L, Guimaraes AF,",
    "Davies J, Waterhouse D, Cook DAN, Myhill LJ, Clare RH, Cassidy A,",
    "Steven A, Johnston KL, Ford L, Turner JD, Ward SA, Taylor MJ.",
    "Minocycline as a re-purposed anti-Wolbachia macrofilaricide: superiority",
    "compared with doxycycline regimens in a murine infection model of human",
    "lymphatic filariasis. Sci Rep. 2016;6:23458. doi:10.1038/srep23458.",
    "Corrigendum: Sci Rep. 2018;8:46934 (author name corrected to Aljayyoussi G;",
    "no change to any value).",
    sep = " "
  )
  vignette <- "Sharma_2016_tetracyclines_wolbachia"

  # Every PK parameter is published per kg body weight (CL/F in L/h/kg, V/F in
  # L/kg; Table 2) and every dose in mg/kg (Methods 'Preparation and
  # administration of Drugs/Compounds'). Doses are therefore given in mg/kg and
  # compartment amounts carried in mg/kg, so central / vc is mg/L = ug/mL, the
  # unit of the whole-blood LC-MS assay and of Figure 1.
  units <- list(time = "h", dosing = "mg/kg", concentration = "ug/mL")

  # The Wolbachia depletion output is a percentage reduction of the wsp gene
  # copy number per adult female worm against the vehicle-treated control
  # (Methods 'Pharmacokinetics-Pharmacodynamic Modelling'). It has no
  # registered canonical name; it is declared paper-specific until a second
  # model of this endpoint arrives.
  paper_specific_compartments <- c("wolbachia_reduction")

  compartmentData <- list(
    depot = list(
      analyte = "doxycycline",
      units = "mg/kg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "doxycycline",
      units = "mg/kg",
      specimen = "whole blood",
      verified = TRUE
    ),
    effect = list(
      analyte = "doxycycline",
      units = "mg/kg",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "mouse (male BALB/c SCID)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    weight_range = "24-28 g (PK study)",
    sex_female_pct = 0,
    disease_state = paste(
      "PK: uninfected mice. PD: mice infected intraperitoneally with 100",
      "Brugia malayi L3 larvae, treated from 6 weeks post-infection and",
      "necropsied at 12 weeks post-infection (Fig. 3) or immediately after the",
      "last dose (Fig. 4)."
    ),
    dose_range = "25 mg/kg oral gavage once or twice daily for 7, 14, 21, 28 or 42 days",
    regions = "United Kingdom (Liverpool School of Tropical Medicine / University of Liverpool)",
    notes = paste(
      "PK: rich tail-vein microsampling up to 8 h after a single 25 mg/kg dose and",
      "on day 7 of 25 mg/kg once daily; Table 2 reports the day-7 (steady-state)",
      "fit, and day-1 parameters did not differ (Mann-Whitney, n = 5). The number",
      "of PK mice is not stated beyond that n = 5. PD: 4-12 mice per arm, with",
      "10-40 adult female worms per arm assayed for Wolbachia wsp copy number by",
      "qPCR (Fig. 3b legend; Table 3). Modelling in Pmetrics (non-parametric",
      "adaptive grid), so Table 2 reports means and SDs of the parameter",
      "distributions rather than parametric typical values and omegas."
    )
  )

  ini({
    # ---- PK (Methods eq. 1a/1b; Table 2, Doxycycline column) ----------------
    lka <- log(2.53)
    label("First-order absorption rate constant ka (1/h)")
    # NOT reported in the paper. Digitised from Figure 1 (left panel): the median
    # day-7 profile peaks at 3.02 ug/mL; with the Table 2 CL/F and V/F under
    # 25 mg/kg once daily for 7 days, ka = 2.53 /h reproduces that Cmax.

    lcl <- log(1.2)
    label("Apparent clearance CL/F (L/h/kg)")
    # Table 2, 'CL/F (L/hr/Kg)' Doxycycline = 1.2 (SD 0.4)

    lvc <- log(6.9)
    label("Apparent volume of distribution V/F (L/kg)")
    # Table 2, 'V/F (L/Kg)' Doxycycline = 6.9 (SD 3.9)

    # ---- PD (Methods eq. 2a/2b/3; Table 4, Doxycycline row) -----------------
    lk1e <- log(4.36e-7)
    label("Transfer rate constant from the central amount to the effect compartment K12 (1/h)")
    # Table 4, 'K12/K21' Doxycycline = 4.36e-7 / 8.06e-4; time unit not printed, read as 1/h (the PK time unit)

    lke1 <- log(8.06e-4)
    label("Transfer rate constant from the effect compartment back to central K21 (1/h)")
    # Table 4, 'K12/K21' Doxycycline = 4.36e-7 / 8.06e-4

    lec50 <- fixed(log(15.32))
    label("Effect-site EC50 for Wolbachia depletion (nM)")
    # Table 4, 'Effect site EC50 (nM)' Doxycycline = 15.32; fixed to the in vitro C6/36 (wAlbB) EC50 (Table 1, 15.3 nM; Methods)

    lemax <- fixed(log(99.999))
    label("Maximum Wolbachia depletion Emax (%)")
    # Methods 'Pharmacokinetics-Pharmacodynamic Modelling': 'Emax 99.999%'

    lveff_ratio <- fixed(log(0.023))
    label("Effect-site volume Vf as a fraction of V/F (unitless)")
    # NOT reported. Methods: 'Vf is assumed to be proportional to the volume of
    # distribution of the central compartment'. The constant is back-solved from
    # Table 4 'Ratio of AUCPD curve (0->28d)' = 1.7 (minocycline / doxycycline,
    # 25 mg/kg bid for 28 days), one value shared by both drugs: the ratio is
    # 1.70 at 0.0230 and spans 1.65-1.75 for 0.021-0.025. See the vignette.

    # ---- IIV ----------------------------------------------------------------
    # Pmetrics reports the mean and SD of each non-parametric marginal
    # distribution (Table 2). They are carried as log-normal variances with the
    # same CV: omega^2 = log(1 + (SD / mean)^2). No dispersion is reported for
    # ka, K12, K21, so those carry no IIV.
    etalcl ~ 0.10536 # Table 2 CL/F SD 0.4 / mean 1.2 = CV 33.3 %; log(1 + 0.3333^2)
    etalvc ~ 0.27723 # Table 2 V/F SD 3.9 / mean 6.9 = CV 56.5 %; log(1 + 0.5652^2)

    # ---- Residual error -----------------------------------------------------
    propSd <- fixed(0)
    label("Proportional residual error on the blood concentration (fraction); not reported")
    # Pmetrics assay-error polynomial and gamma/lambda are not reported
    addSd_wolbachia_reduction <- fixed(0)
    label("Additive residual error on Wolbachia depletion (%); not reported")
    # not reported
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    kel <- cl / vc
    k1e <- exp(lk1e)
    ke1 <- exp(lke1)
    ec50 <- exp(lec50)
    emax <- exp(lemax)
    # Methods: Vf proportional to the central volume of distribution
    veff <- exp(lveff_ratio) * vc

    # Methods eq. 1a, 2a, 2b: the effect compartment exchanges drug AMOUNT with
    # central, so K12 also removes drug from central (negligible here,
    # K12 / kel = 2.5e-6).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k1e * central + ke1 * effect
    d/dt(effect) <- k1e * central - ke1 * effect

    Cc <- central / vc

    # Methods eq. 3 with the effect-site concentration converted from mg/L to nM
    # (doxycycline free base, C22H24N2O8, 444.43 g/mol) to match the EC50 units.
    Ce_nM <- effect / veff / 444.43 * 1e6
    wolbachia_reduction <- emax * Ce_nM / (ec50 + Ce_nM)

    Cc ~ prop(propSd)
    wolbachia_reduction ~ add(addSd_wolbachia_reduction)
  })
}
