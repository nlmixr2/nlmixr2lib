Dimelow_2018_ceftazidime <- function() {
  description <- paste(
    "Three-compartment intravenous population PK model for ceftazidime in the plasma of healthy male",
    "volunteers, jointly fitted with a direct-response (instantaneous-equilibrium) saturable plasma-to-",
    "epithelial-lining-fluid (ELF) link. Plasma disposition is parameterized by CL, V1, V2, V3, Q12 and",
    "Q13 with log-normal between-subject variability on every structural parameter. The ELF observable",
    "Celf is algebraic, not a compartment: an effect-site equilibration model was fitted and rejected",
    "(ELF half-life 13 min, no OFV improvement), so ELF concentration is an instantaneous",
    "Michaelis-Menten function of the total plasma concentration, Celf = emax_elf * Cc / (km_elf + Cc).",
    "The ELF:plasma penetration ratio therefore falls as plasma concentration rises: emax_elf / km_elf =",
    "63.3% as Cc approaches zero, 52% at the efficacy-relevant Cc of 15.3 mg/L, and 32.1% at a typical",
    "Cmax of 70 mg/L. Ceftazidime plasma proportional residual noise itself carries between-subject",
    "variability (median 0.101, 54% CV); the ELF residual could not be identified from one ELF sample",
    "per subject and was FIXED to the plasma median. Penetration ratios are relative to TOTAL plasma and",
    "TOTAL ELF concentrations.",
    sep = " "
  )
  reference <- paste(
    "Dimelow R, Wright JG, MacPherson M, Newell P, Das S.",
    "Population pharmacokinetic modelling of ceftazidime and avibactam in the plasma and epithelial",
    "lining fluid of healthy volunteers. Drugs R D. 2018;18(3):221-230.",
    "doi:10.1007/s40268-018-0241-0 (Sects. 2.2, 2.2.1, 3.1, 3.2; Tables 1 and 2).",
    "Subject data are from the phase I open-label ELF study NCT01395420.",
    sep = " "
  )
  vignette <- "Dimelow_2018_ceftazidime_avibactam_elf"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "ceftazidime", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ceftazidime", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "ceftazidime", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 42,
    n_studies = 1,
    sex_female_pct = 0,
    disease_state = "healthy volunteers",
    dose_range = paste(
      "ceftazidime 2000 mg + avibactam 500 mg (cohort A, n = 22) or ceftazidime 3000 mg +",
      "avibactam 1000 mg (cohort B, n = 20), each given as a 2 h intravenous infusion every 8 h for",
      "3 days (nine doses per subject)",
      sep = " "
    ),
    notes = paste(
      "Sect. 2.1: a previously reported phase I open-label study (NCT01395420) enrolled 43 healthy male",
      "volunteers; PK data were available for 42 (22 cohort A, 20 cohort B). One bronchoalveolar-lavage",
      "ELF sample was taken per subject, which is why the ELF residual error is not identifiable.",
      "The paper reports no age, weight or race distribution and screened no covariates, so none are",
      "encoded.",
      sep = " "
    )
  )

  ini({
    # Plasma disposition -- Table 1, ceftazidime rows. Values are population
    # medians; the trailing +/- figure in the source table is the relative
    # standard error of the estimate, not a standard deviation.
    lcl <- log(6.55); label("Clearance (L/h)") # Table 1 ceftazidime CL = 6.55 L/h (RSE 2%)
    lvc <- log(10.32); label("Central volume of distribution, V1 (L)") # Table 1 ceftazidime V1 = 10.32 L (RSE 8%)
    lvp <- log(5.82); label("First peripheral volume of distribution, V2 (L)") # Table 1 ceftazidime V2 = 5.82 L (RSE 9%)
    lvp2 <- log(0.64); label("Second peripheral volume of distribution, V3 (L)") # Table 1 ceftazidime V3 = 0.64 L (RSE 8%)
    lq <- log(6.87); label("Intercompartmental clearance to peripheral1, Q12 (L/h)") # Table 1 ceftazidime Q12 = 6.87 L/h (RSE 24%)
    lq2 <- log(0.040); label("Intercompartmental clearance to peripheral2, Q13 (L/h)") # Table 1 ceftazidime Q13 = 0.040 L/h (RSE 33%)

    # Plasma-to-ELF link -- Table 2, ceftazidime rows. Direct response with a
    # saturable (Michaelis-Menten) link function, selected over the proportional
    # (35 OFV units worse) and power (8 OFV units worse) link functions in
    # Sect. 3.2.
    lemax_elf <- log(45.4); label("Maximum attainable ELF ceftazidime concentration, EMAX (mg/L)") # Table 2 ceftazidime EMAX = 45.4 mg/L (RSE 12%)
    lkm_elf <- log(71.7); label("Plasma concentration giving half the maximum ELF concentration, KM (mg/L)") # Table 2 ceftazidime KM = 71.7 mg/L (RSE 22%)

    # Between-subject variability. Sect. 2.2: 'Between-subject variability in
    # parameters was modelled as being log-normally distributed'. The source
    # reports it as a %CV (Table 1 / Table 2, final column), so the variance is
    # recovered as log(1 + CV^2); the arithmetic is left inline so the published
    # CV stays visible in the source trace. Sect. 3.1 states plasma BSV used a
    # FULL covariance matrix, but no off-diagonal element is published, so the
    # block is carried as diagonal (see vignette Errata).
    etalcl ~ log(1 + 0.10^2) # Table 1 ceftazidime CL, BSV CV = 10% (5th-95th percentile 5.59-7.68 L/h)
    etalvc ~ log(1 + 0.18^2) # Table 1 ceftazidime V1, BSV CV = 18% (7.73-13.78 L)
    etalvp ~ log(1 + 0.18^2) # Table 1 ceftazidime V2, BSV CV = 18% (4.37-7.75 L)
    etalvp2 ~ log(1 + 0.25^2) # Table 1 ceftazidime V3, BSV CV = 25% (0.43-0.95 L)
    etalq ~ log(1 + 0.66^2) # Table 1 ceftazidime Q12, BSV CV = 66% (2.54-18.55 L/h)
    etalq2 ~ log(1 + 0.61^2) # Table 1 ceftazidime Q13, BSV CV = 61% (0.016-0.101 L/h)
    etalemax_elf ~ log(1 + 0.24^2) # Table 2 ceftazidime EMAX, BSV CV = 24% (30.9-66.5 mg/L)
    etalkm_elf ~ log(1 + 0.97^2) # Table 2 ceftazidime KM, BSV CV = 97% (18.8-273.3 mg/L)

    # Residual error. Sect. 2.2: plasma residual noise is proportional to the
    # predicted plasma concentration. Table 1 reports RESMp as a population
    # MEDIAN with its own between-subject variability, so the proportional SD is
    # itself log-normally distributed across subjects.
    propSd <- 0.101; label("Median plasma proportional residual SD (fraction)") # Table 1 ceftazidime RESMp = 0.101 (RSE 9%)
    etapropSd ~ log(1 + 0.54^2) # Table 1 ceftazidime RESMp, BSV CV = 54% (0.044-0.233)

    # Sect. 2.2 and Sect. 4: with one ELF observation per individual the ELF
    # residual noise could not be estimated and was fixed to the plasma median.
    # Table 2 confirms it carries no between-subject variability (CV 0%).
    propSd_Celf <- fixed(0.101); label("ELF proportional residual SD, carried from the plasma median (fraction)") # Table 2 ceftazidime RESM_ELF = 0.101, CV 0%; Sect. 2.2 states it was fixed to the plasma RESMp median
  })

  model({
    # Individual plasma disposition parameters
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    vp2 <- exp(lvp2 + etalvp2)
    q <- exp(lq + etalq)
    q2 <- exp(lq2 + etalq2)

    # Individual plasma-to-ELF link parameters
    emax_elf <- exp(lemax_elf + etalemax_elf)
    km_elf <- exp(lkm_elf + etalkm_elf)

    # Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # Sect. 2.2: 'coded as a series of differential equations'
    d/dt(central) <- -(kel + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # Observations. Cc is the TOTAL plasma concentration; Celf is the TOTAL ELF
    # concentration. Sect. 2.2.1.1 saturable link: Celf = EMAX * Cp / (KM + Cp).
    # Celf is algebraic rather than a compartment because the effect-site
    # equilibration model was rejected (Sect. 3.2).
    Cc <- central / vc
    Celf <- emax_elf * Cc / (km_elf + Cc)

    # The plasma proportional SD is itself a per-subject quantity (Table 1
    # reports a median RESMp with 54% CV), so the residual magnitude is scaled
    # by its own eta.
    propSdi <- propSd * exp(etapropSd)

    Cc ~ prop(propSdi)
    Celf ~ prop(propSd_Celf)
  })
}
