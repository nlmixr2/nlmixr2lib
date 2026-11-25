Srivastava_2017_linezolid <- function() {
  description <- paste(
    "One-compartment oral population pharmacokinetic model of linezolid",
    "used by Srivastava 2017 for Monte Carlo dose selection in 10,000 adult",
    "tuberculosis patients (Table 2), with log-normal between-subject",
    "variability on CL, V and ka and no covariates. The typical values and",
    "IIV are those of Abe et al. 2009 (435 patients) as reprinted in Table 2.",
    "Adds the paper's translation to the site of action: the free",
    "epithelial-lining-fluid concentration Celf = Cc x 3.29 (ELF:serum",
    "AUC0-24 penetration ratio) x 0.7 (fraction unbound, 30% protein",
    "binding), whose AUC0-24 divided by the MIC is compared with the",
    "hollow-fiber EC80 target of 119 (Srivastava_2017_linezolid_hfs_aucmic),",
    "and the unbound serum concentration Cu = 0.7 x Cc, whose AUC0-24 is",
    "compared with the 96 mg*h/L mitochondrial-toxicity threshold. No residual",
    "error was used in the simulations, so propSd is fixed at 0."
  )
  reference <- paste(
    "Srivastava S, Magombedze G, Koeuth T, Sherman C, Pasipanodya JG,",
    "Raj P, Wakeland E, Deshpande D, Gumbo T. Linezolid dose that maximizes",
    "sterilizing effect while minimizing toxicity and resistance emergence for",
    "tuberculosis. Antimicrob Agents Chemother. 2017;61(8):e00751-17.",
    "doi:10.1128/AAC.00751-17. PMCID: PMC5527615. Parameter values from",
    "Table 2, reprinted from Abe S, Chiba K, Cirincione B, Grasela TH, Ito K,",
    "Suwa T. Population pharmacokinetic analysis of linezolid in patients with",
    "infectious disease: application to lower body weight and elderly",
    "patients. J Clin Pharmacol. 2009;49(9):1071-1078.",
    "doi:10.1177/0091270009337947."
  )
  vignette <- "Srivastava_2017_linezolid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "linezolid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "linezolid", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = 435L,
    n_studies = NA_integer_,
    disease_state = "adults with infectious disease (Abe 2009 source population); simulated as adult tuberculosis patients",
    dose_range = "oral linezolid 150, 300, 600 and 1200 mg once daily (simulated); 1800 mg/day also shown in Figure 6",
    regions = NA_character_,
    notes = paste(
      "Srivastava 2017 reports only the means and IIV (as %CV) of CL, V and",
      "ka from Abe et al. 2009 and fed them to the ADAPT 5 PRIOR subroutine",
      "for a 10,000-patient Monte Carlo simulation; Abe's covariate model",
      "(body weight, age) is not reported and not encoded. The paper's",
      "simulated mean AUC0-24 at 600 mg/day was 107.5 +/- 30.16 mg*h/L",
      "(Results)."
    )
  )

  ini({
    lka <- log(0.58)
    label("Absorption rate constant (1/h)") # Table 2, 'Absorption constant (h-1)' mean 0.58 (435 patients)
    lcl <- log(6.0)
    label("Apparent clearance CL (L/h)") # Table 2, 'Clearance (liters/h)' mean 6.0 (435 patients)
    lvc <- log(54.53)
    label("Apparent volume of distribution V (L)") # Table 2, 'Vol (liters)' mean 54.53 (435 patients)

    # Site-of-action translation used for target attainment (Methods, 'Dose selection').
    lrelf <- fixed(log(3.29))
    label("ELF:serum AUC0-24 penetration ratio (unitless)") # Methods, 'The linezolid ELF-to-serum AUC0-24 penetration ratio of 3.29 was also taken into account' (sensitivity analysis: 1.52)
    fu <- fixed(0.7)
    label("Fraction unbound in serum (unitless)") # Methods, 'as was the 30% protein binding'

    # IIV printed as %CV; omega^2 = log(1 + CV^2).
    etalcl ~ 0.19645 # Table 2, 'IIV (%)' for clearance 46.6 -> log(1 + 0.466^2)
    etalvc ~ 0.064930 # Table 2, 'IIV (%)' for volume 25.9 -> log(1 + 0.259^2)
    etalka ~ 0.14842 # Table 2, 'IIV (%)' for absorption constant 40 -> log(1 + 0.40^2)

    propSd <- fixed(0)
    label("Proportional residual error (not used in the simulations)") # No residual error reported or used in the Monte Carlo simulation
  })

  model({
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    relf <- exp(lrelf)
    kel <- cl / vc

    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - kel * central

    Cc <- central / vc
    # Unbound serum concentration: driver of the 96 mg*h/L toxicity threshold.
    Cu <- fu * Cc
    # Free epithelial-lining-fluid concentration, the driver of the EC80
    # target. The paper applies both the ELF:serum ratio and the 30% protein
    # binding, so this is relf x the UNBOUND serum concentration.
    Celf <- relf * Cu

    Cc ~ prop(propSd)
  })
}
