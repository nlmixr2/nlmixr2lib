Hammer_2020_acetaminophen <- function() {
  description <- "Two-compartment population PK model for intravenous acetaminophen (paracetamol) in neonates, infants, children and adolescents (Hammer 2020), updating the Zuppa/Palmer pediatric model with neonate and infant data from a randomized placebo-controlled postoperative-pain study. All clearances scale allometrically with body weight (fixed exponent 0.75) and both volumes linearly with weight (reference 70 kg); systemic clearance additionally carries an exponential postmenstrual-age maturation function (1 - 0.611 * exp(-(PMA - 40) * ln 2 / 32.6 weeks)) and a 0.524-fold multiplier in the placebo arm, whose subjects had low residual acetaminophen concentrations. Log-normal residual error; correlated IIV on CL and Vc."
  reference <- paste(
    "Hammer GB, Maxwell LG, Taicher BM, Visoiu M, Cooper DS, Szmuk P,",
    "Pheng LH, Gosselin NH, Lu J, Devarakonda K (2020).",
    "Randomized population pharmacokinetic analysis and safety of",
    "intravenous acetaminophen for acute postoperative pain in neonates",
    "and infants. J Clin Pharmacol 60(1):16-27.",
    "doi:10.1002/jcph.1508. PMID 31448420; PMCID PMC6973014.",
    sep = " "
  )
  vignette <- "Hammer_2020_acetaminophen"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "acetaminophen", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "acetaminophen", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling with reference 70 kg: exponent 0.75 on CL and on the inter-compartmental clearance CLp, exponent 1 on Vc and Vp (Hammer 2020 Table 3, current-model column). The screening weight (Table 2) is the only weight reported.",
      source_name = "WT"
    ),
    PAGE = list(
      description = "Postmenstrual age (gestational age + postnatal age)",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "Drives the exponential maturation function on CL, 1 - 0.611 * exp(-(PMA - 40) * ln(2) / 32.6), whose 40-week pivot and 32.6-week half-life are on the WEEKS scale (Hammer 2020 Table 3 and Figure 4 axis). Supply PAGE in weeks.",
      source_name = "PMA"
    ),
    PLACEBO = list(
      description = "Randomized placebo-arm membership indicator (1 = placebo control arm, 0 = intravenous acetaminophen arm)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (active intravenous acetaminophen arm)",
      notes = "Hammer 2020 Table 3 'Placebo on CL' = x0.524: systemic CL in the placebo control groups (C + D) is 0.524 times the active-arm value, encoded as e_placebo_cl^PLACEBO. The placebo subjects received saline study drug yet had low residual acetaminophen concentrations (Results 'Population PK', Figure 1); the paper does not describe the dosing history that produced them. For simulating the approved IV regimen set PLACEBO = 0.",
      source_name = "TREATMENT (Placebo vs Active)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 283L,
    n_studies = 3L,
    age_range = "Current study: neonates and infants under 24 months (postnatal age 1-725 days). Pooled with the Zuppa (neonates to adolescents) and Palmer (neonates and infants) data, PMA about 33.7 to 906 weeks (about 0 to 17 years; Figure 4).",
    age_median = "Current study efficacy population: postnatal age 195 days (Table 2)",
    weight_range = "Current study efficacy population: 0.8-14.3 kg (Table 2)",
    weight_median = "Current study efficacy population: 7.0 kg (Table 2)",
    sex_female_pct = 35.5,
    race_ethnicity = "Current study efficacy population: White 68.5%, Black or African American 15.2%, Asian 6.6%, American Indian or Alaska Native 0.5%, other 6.6%, missing 2.5% (Table 2)",
    disease_state = "Current study: surgical neonates and infants with acute postoperative pain expected to need at least 24 h of pain management (none had traumatic injury), randomized to IV acetaminophen low dose, high dose, or saline placebo, all with standard-of-care opioids. Pooled prior data: pediatric patients from Zuppa et al. and Palmer et al.",
    dose_range = "IV acetaminophen as a 15-minute infusion every 6 h for 4 doses. Low dose (group A): 7.5 mg/kg (extreme preterm neonates) to 12.5 mg/kg (infants); high dose (group B): 10 mg/kg (extreme preterm neonates) to 15 mg/kg (infants).",
    regions = "United States (multicentre, NCT01635101)",
    notes = "The final model was fitted to 581 samples from 158 subjects of the current study combined with the 1260 samples from 125 pediatric subjects (neonates through adolescents) of the previously developed model (Zuppa and Palmer studies), 1841 samples in total. Demographics are from Table 2 (efficacy population, n = 197); the 158-subject PK population was reported not to deviate significantly from it."
  )

  ini({
    # Hammer 2020 Table 3, 'Current Model' column (page 20). Typical values
    # are at the 70 kg reference weight.
    lcl <- log(18.9); label("Systemic clearance at 70 kg, fully matured, active arm (L/h)") # Table 3 current model: CL = 18.9 x (WT/70)^0.75
    lvc <- log(23.0); label("Central volume of distribution at 70 kg (L)") # Table 3 current model: Vc = 23.0 x (WT/70)
    lq <- log(47.7); label("Inter-compartmental clearance at 70 kg (L/h)") # Table 3 current model: CLp = 47.7 x (WT/70)^0.75
    lvp <- log(45.5); label("Peripheral volume of distribution at 70 kg (L)") # Table 3 current model: Vp = 45.5 x (WT/70)

    # Allometric exponents printed without uncertainty in Table 3 (a standard
    # fixed theory-based allometry).
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent on CL and CLp (unitless)") # Table 3: (WT/70)^0.75 under CL and CLp
    e_wt_vc_vp <- fixed(1); label("Allometric exponent on Vc and Vp (unitless)") # Table 3: (WT/70) under Vc and Vp

    # Placebo-arm multiplier on CL.
    e_placebo_cl <- 0.524; label("Multiplier on CL in the placebo control arm, applied as base^PLACEBO (unitless)") # Table 3 current model: 'Placebo on CL' x0.524

    # Maturation of CL: 1 - e_page_cl * exp(-(PMA - 40) * ln(2) / thalf_cl).
    # The 40-week pivot is a hardcoded constant of the printed equation.
    e_page_cl <- 0.611; label("Fractional reduction of CL at 40 weeks PMA due to immaturity (unitless)") # Table 3 current model maturation function: 1 - 0.611 x exp(...)
    thalf_cl <- fixed(32.6); label("Maturational half-life of the CL deficit (weeks PMA)") # Table 3 current model: ln(2)/32.6 carries the footnote marker, 'Fixed at previous value derived in.'

    # IIV: Table 3 reports omega^2 (log-scale variances) on CL and Vc only,
    # with the note 'the correlation between CL and Vc BSV was 51.4%'.
    # Covariance = 0.514 * sqrt(0.127 * 0.993) = 0.18253.
    etalcl + etalvc ~ c(
      0.127,
      0.18253, 0.993
    ) # Table 3: omega2 CL 0.127, omega2 Vc 0.993, correlation 51.4 percent

    # Residual error: a single 'Log residual error' row = 0.221. Encoded as
    # the log-scale SD (not a variance); the Figure 2 VPC tails discriminate
    # the reading (see vignette 'Assumptions and deviations').
    expSd <- 0.221; label("Log-scale residual error SD (unitless)") # Table 3 current model: Log residual error 0.221
  })

  model({
    # Size and maturation (Hammer 2020 Table 3).
    allom_cl <- (WT / 70)^e_wt_cl_q
    allom_v <- (WT / 70)^e_wt_vc_vp
    fmat <- 1 - e_page_cl * exp(-(PAGE - 40) * log(2) / thalf_cl)

    cl <- exp(lcl + etalcl) * allom_cl * fmat * e_placebo_cl^PLACEBO
    vc <- exp(lvc + etalvc) * allom_v
    q <- exp(lq) * allom_cl
    vp <- exp(lvp) * allom_v

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, volume in L -> mg/L (= ug/mL).
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
