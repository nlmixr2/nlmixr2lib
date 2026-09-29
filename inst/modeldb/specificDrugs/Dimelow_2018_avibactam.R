Dimelow_2018_avibactam <- function() {
  description <- paste(
    "Three-compartment intravenous population PK model for avibactam in the plasma of healthy male",
    "volunteers, jointly fitted with a direct-response (instantaneous-equilibrium) power-function",
    "plasma-to-epithelial-lining-fluid (ELF) link. Plasma disposition is parameterized by CL, V1, V2,",
    "V3, Q12 and Q13 with log-normal between-subject variability on every structural parameter. The ELF",
    "observable Celf is algebraic, not a compartment: an effect-site equilibration model was fitted and",
    "rejected (ELF half-life 8 min, no OFV improvement), so ELF concentration is an instantaneous power",
    "function of the total plasma concentration, Celf = relf * Cc^pow_elf, where relf is the ELF:plasma",
    "penetration ratio at a plasma concentration of 1 mg/L. Because pow_elf is below 1, penetration",
    "falls as plasma concentration rises: 47.2% at Cc = 1 mg/L, 42% at the efficacy-relevant Cc of",
    "2.4 mg/L, and 33.2% at a typical Cmax of 12 mg/L. Avibactam plasma proportional residual noise",
    "itself carries between-subject variability (median 0.117, 45% CV); the ELF residual could not be",
    "identified from one ELF sample per subject and was FIXED to the plasma median. Penetration ratios",
    "are relative to TOTAL plasma and TOTAL ELF concentrations.",
    sep = " "
  )
  reference <- paste(
    "Dimelow R, Wright JG, MacPherson M, Newell P, Das S.",
    "Population pharmacokinetic modelling of ceftazidime and avibactam in the plasma and epithelial",
    "lining fluid of healthy volunteers. Drugs R D. 2018;18(3):221-230.",
    "doi:10.1007/s40268-018-0241-0 (Sects. 2.2, 2.2.1, 3.1, 3.3; Tables 1 and 2).",
    "Subject data are from the phase I open-label ELF study NCT01395420.",
    sep = " "
  )
  vignette <- "Dimelow_2018_ceftazidime_avibactam_elf"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "avibactam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "avibactam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "avibactam", units = "mg", specimen = "plasma", verified = TRUE)
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
    # Plasma disposition -- Table 1, avibactam rows. Values are population
    # medians; the trailing +/- figure in the source table is the relative
    # standard error of the estimate, not a standard deviation.
    lcl <- log(12.5); label("Clearance (L/h)") # Table 1 avibactam CL = 12.5 L/h (RSE 1%)
    lvc <- log(15.10); label("Central volume of distribution, V1 (L)") # Table 1 avibactam V1 = 15.10 L (RSE 10%)
    lvp <- log(6.52); label("First peripheral volume of distribution, V2 (L)") # Table 1 avibactam V2 = 6.52 L (RSE 15%)
    lvp2 <- log(1.58); label("Second peripheral volume of distribution, V3 (L)") # Table 1 avibactam V3 = 1.58 L (RSE 5%)
    lq <- log(5.43); label("Intercompartmental clearance to peripheral1, Q12 (L/h)") # Table 1 avibactam Q12 = 5.43 L/h (RSE 42%)
    lq2 <- log(0.14); label("Intercompartmental clearance to peripheral2, Q13 (L/h)") # Table 1 avibactam Q13 = 0.14 L/h (RSE 41%)

    # Plasma-to-ELF link -- Table 2, avibactam rows. Direct response with a
    # power link function, selected over the proportional (22 OFV units worse)
    # and saturable (6 OFV units worse) link functions in Sect. 3.3.
    lrelf <- log(0.472); label("ELF / total-plasma penetration ratio at a plasma concentration of 1 mg/L (unitless)") # Table 2 avibactam 'EPR (1 mg/l)' = 0.472 (RSE 11%)
    lpow_elf <- log(0.860); label("Power exponent of total plasma concentration in the plasma-to-ELF link (unitless)") # Table 2 avibactam POW = 0.860 (RSE 5%)

    # Between-subject variability. Sect. 2.2: 'Between-subject variability in
    # parameters was modelled as being log-normally distributed'. The source
    # reports it as a %CV (Table 1 / Table 2, final column), so the variance is
    # recovered as log(1 + CV^2); the arithmetic is left inline so the published
    # CV stays visible in the source trace. Sect. 3.1 states plasma BSV used a
    # FULL covariance matrix, but no off-diagonal element is published, so the
    # block is carried as diagonal (see vignette Errata).
    etalcl ~ log(1 + 0.07^2) # Table 1 avibactam CL, BSV CV = 7% (5th-95th percentile 11.1-14.17 L/h)
    etalvc ~ log(1 + 0.09^2) # Table 1 avibactam V1, BSV CV = 9% (13.08-17.44 L)
    etalvp ~ log(1 + 0.25^2) # Table 1 avibactam V2, BSV CV = 25% (4.32-9.84 L)
    etalvp2 ~ log(1 + 0.20^2) # Table 1 avibactam V3, BSV CV = 20% (1.15-2.17 L)
    etalq ~ log(1 + 0.60^2) # Table 1 avibactam Q12, BSV CV = 60% (2.19-13.45 L/h)
    etalq2 ~ log(1 + 0.75^2) # Table 1 avibactam Q13, BSV CV = 75% (0.05-0.43 L/h)
    etalrelf ~ log(1 + 0.68^2) # Table 2 avibactam EPR (1 mg/l), BSV CV = 68% (0.172-1.298)
    etalpow_elf ~ log(1 + 0.23^2) # Table 2 avibactam POW, BSV CV = 23% (0.588-1.259)

    # Residual error. Sect. 2.2: plasma residual noise is proportional to the
    # predicted plasma concentration. Table 1 reports RESMp as a population
    # MEDIAN with its own between-subject variability, so the proportional SD is
    # itself log-normally distributed across subjects.
    propSd <- 0.117; label("Median plasma proportional residual SD (fraction)") # Table 1 avibactam RESMp = 0.117 (RSE 7%)
    etapropSd ~ log(1 + 0.45^2) # Table 1 avibactam RESMp, BSV CV = 45% (0.058-0.239)

    # Sect. 2.2 and Sect. 4: with one ELF observation per individual the ELF
    # residual noise could not be estimated and was fixed to the plasma median.
    # Table 2 confirms it carries no between-subject variability (CV 0%).
    propSd_Celf <- fixed(0.117); label("ELF proportional residual SD, carried from the plasma median (fraction)") # Table 2 avibactam RESM_ELF = 0.117, CV 0%; Sect. 2.2 states it was fixed to the plasma RESMp median
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
    relf <- exp(lrelf + etalrelf)
    pow_elf <- exp(lpow_elf + etalpow_elf)

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
    # concentration. Sect. 2.2.1.1 power link: Celf = EPR(1 mg/L) * Cp^POW.
    # Celf is algebraic rather than a compartment because the effect-site
    # equilibration model was rejected (Sect. 3.3).
    Cc <- central / vc
    Celf <- relf * Cc^pow_elf

    # The plasma proportional SD is itself a per-subject quantity (Table 1
    # reports a median RESMp with 45% CV), so the residual magnitude is scaled
    # by its own eta.
    propSdi <- propSd * exp(etapropSd)

    Cc ~ prop(propSdi)
    Celf ~ prop(propSd_Celf)
  })
}
