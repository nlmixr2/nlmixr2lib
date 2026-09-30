Kim_2020_voriconazole <- function() {
  description <- paste(
    "One-compartment population PK model for oral and intravenous voriconazole in hospitalised",
    "patients treated for invasive aspergillosis (Kim 2020), fitted jointly to paired plasma and",
    "saliva concentrations. Oral doses enter a gut depot with bioavailability 0.849 and first-order",
    "absorption (ka fixed at 0.858 1/h); intravenous doses enter the central (plasma) compartment",
    "directly. Saliva is not a separate kinetic compartment: the salivary voriconazole concentration",
    "is the plasma concentration multiplied by an estimated saliva:plasma scale factor of 0.501, a",
    "structure the authors selected over a separate saliva compartment (dOFV = -102.7). IIV on",
    "clearance only; no covariates; separate proportional residual errors for plasma and saliva.",
    sep = " "
  )
  reference <- paste(
    "Kim HY, Martson A-G, Dreesen E, Spriet I, Wicha SG, McLachlan AJ, Alffenaar J-W.",
    "Saliva for Precision Dosing of Antifungal Drugs: Saliva Population PK Model for Voriconazole",
    "Based on a Systematic Review. Front Pharmacol. 2020;11:894. doi:10.3389/fphar.2020.00894",
    sep = " "
  )
  vignette <- "Kim_2020_voriconazole"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    # Oral doses only (Figure 2: 'Oral dose -> Bioavailability (F) -> Gut').
    depot = list(analyte = "voriconazole", units = "mg", specimen = "administration site", verified = TRUE),
    # Plasma. Intravenous doses enter here directly (Figure 2: 'IV dose -> Plasma').
    # Saliva is NOT a compartment: Csaliva is an algebraic rescaling of the plasma
    # concentration by fsaliva (Figure 2 dashed 'Saliva' box; Results).
    central = list(analyte = "voriconazole", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 11,
    n_studies = 1,
    age_range = "30-66 years in the 10 adults (median 55), plus one 9-year-old child",
    age_median = "55 years (adults)",
    weight_range = "not reported; adult mean 65.9 kg (SD 20.1)",
    sex_female_pct = 36.4,
    disease_state = "patients on adult oncology/haematology (8) or respiratory (3) wards or a paediatric ward treated with voriconazole for invasive aspergillosis",
    dose_range = "3.7 +/- 0.4 mg/kg voriconazole orally (tablet) or intravenously every 12 h for at least 4 days (Table 2, Vanstraelen 2015 row)",
    sampling = "steady-state paired plasma and saliva at pre-dose, 0.5, 1, 1.5, 2, 6 and 12 h post-dose; 69 plasma and 68 saliva concentrations",
    notes = paste(
      "Kim 2020 Results 'Data Retrieval From Authors' Studies': individual data re-analysed from",
      "Vanstraelen et al. 2015 (Ther Drug Monit 37:766-771). 7 male / 4 female. The paediatric",
      "patient, excluded from the original study, was included after removing one outlier saliva",
      "concentration. Saliva collected with Salivette; concentrations by LC-MS/MS."
    )
  )

  ini({
    lka <- fixed(log(0.858)); label("Absorption rate constant (1/h)")                 # Table 4 theta_3 = 0.858 'fixed to model estimate'; Results: fixed because RSE was 129% when estimated
    lcl <- log(4.56); label("Clearance (L/h)")                                         # Table 4 theta_1 = 4.56 L/h (RSE 16%; bootstrap 4.39 [3.23-5.98])
    lvc <- log(60.7); label("Volume of distribution (L)")                              # Table 4 theta_2 = 60.7 L (RSE 12%; bootstrap 57.9 [41.4-72.3])
    lfdepot <- log(0.849); label("Oral bioavailability (fraction)")                    # Table 4 theta_4 = 0.849 (RSE 14%; bootstrap 0.819 [0.577-0.983])

    # Saliva:plasma concentration scale factor applied to the plasma compartment
    # (Figure 2; Results: scale factor preferred over a separate saliva compartment, dOFV -102.7).
    lfsaliva <- log(0.501); label("Saliva:plasma concentration scale factor (unitless)") # Table 4 theta_5 = 0.501 (RSE 4%; bootstrap 0.499 [0.458-0.541])

    # IIV on CL only; IIV on V was not estimable (99% shrinkage, Results).
    etalcl ~ 0.136  # Table 4 CL omega^2 = 0.136 (36.9% CV; RSE 42%; bootstrap 0.115)

    # Residual error: Table 4 reports sigma^2 (variances); SD = sqrt(sigma^2).
    propSd <- sqrt(0.057); label("Proportional residual SD for plasma Cc (fraction)")             # Table 4 sigma^2 proportional, plasma = 0.057 (RSE 25%) -> SD 0.239
    propSd_Csaliva <- sqrt(0.078); label("Proportional residual SD for saliva Csaliva (fraction)") # Table 4 sigma^2 proportional, saliva = 0.078 (RSE 26%) -> SD 0.279
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    fsaliva <- exp(lfsaliva)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    f(depot) <- exp(lfdepot)

    # Plasma is the central concentration; saliva is that concentration scaled
    # by fsaliva (Figure 2).
    Cc <- central / vc
    Csaliva <- fsaliva * Cc

    Cc ~ prop(propSd)
    Csaliva ~ prop(propSd_Csaliva)
  })
}
