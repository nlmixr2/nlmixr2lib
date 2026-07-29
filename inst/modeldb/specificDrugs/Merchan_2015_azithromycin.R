Merchan_2015_azithromycin <- function() {
  description <- "Population PK model for intravenous azithromycin in preterm neonates at risk for Ureaplasma respiratory tract colonization (Merchan 2015). Pooled re-analysis of three studies (single 10 mg/kg, single 20 mg/kg, and 3 daily doses of 20 mg/kg). Two-compartment linear model with all PK parameters allometrically scaled on body weight: fixed exponent 0.75 on CL and Q, fixed exponent 1.0 on V1 and V2, reference body weight 1 kg."
  reference <- "Merchan LM, Hassan HE, Terrin ML, Waites KB, Kaufman DA, Ambalavanan N, Donohue P, Dulkerian SJ, Schelonka R, Magder LS, Shukla S, Eddington ND, Viscardi RM. Pharmacokinetics, microbial response, and pulmonary outcomes of multidose intravenous azithromycin in preterm infants at risk for Ureaplasma respiratory colonization. Antimicrob Agents Chemother. 2015;59(1):570-578. doi:10.1128/AAC.03951-14"
  vignette <- "Merchan_2015_azithromycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric exponents fixed a priori at 0.75 on CL and Q and at 1.0 on V1 and V2 (Merchan 2015 Methods 'Pharmacokinetic data analysis'). Reference weight 1 kg, indicated by the parameter units CL = L/h/kg^0.75 and V1 = L/kg reported in Table 3 and confirmed by the worked half-life calculation 'a typical neonate weighing 1 kg' in Results.",
      source_name        = "WT"
    )
  )

  population <- list(
    species         = "human (preterm neonates)",
    n_subjects      = 40,
    n_studies       = 3,
    age_range       = "postnatal age <72 h at first dose; gestational age at birth 24-28 weeks (multidose cohort)",
    weight_range    = "preterm neonate birth weight, multidose cohort 856 +/- 202 g (Ureaplasma positive) and 929 +/- 285 g (Ureaplasma negative)",
    sex_female_pct  = NA_real_,
    race_ethnicity  = "Multidose cohort (n=15): White ~67%, Black ~33%; ~13% Hispanic ethnicity (Table 1).",
    disease_state   = "Mechanically ventilated preterm infants at risk for Ureaplasma respiratory tract colonization and bronchopulmonary dysplasia.",
    dose_range      = "Pooled across three studies: single 10 mg/kg IV (n=12), single 20 mg/kg IV (n=13), and 3 doses of 20 mg/kg IV q24h (n=15). All infusions over 60 min.",
    n_observations  = "239 plasma azithromycin concentrations from 40 subjects.",
    regions         = "United States (six neonatal intensive care sites, December 2011 - June 2012 enrolment for the multidose cohort; the two earlier single-dose studies are cited as references 14 and 15).",
    fda_ind         = "FDA IND 78990 (multidose study); related trial NCT01778634.",
    notes           = "Demographic counts and ranges reproduced from Merchan 2015 Methods and Table 1. Covariates explored but not retained in the final model: gestational age, sex, height, body surface area (Methods 'Pharmacokinetic data analysis'); interoccasion variability across the three studies was also tested and dropped. The multidose cohort detail in Table 1 is the only stratified-demographic table; baseline characteristics from the two earlier single-dose studies are summarised in their primary references (14 and 15) and are not reproduced here."
  )

  ini({
    # Structural parameters at the paper's reference weight 1 kg (Table 3).
    lcl <- log(0.15);  label("Clearance, CL, at 1 kg (L/h)")                              # Merchan 2015 Table 3: 0.15 L/h/kg^0.75 (10% RSE)
    lvc <- log(1.88);  label("Central volume of distribution, V1, at 1 kg (L)")            # Merchan 2015 Table 3: 1.88 L/kg (11% RSE)
    lq  <- log(1.79);  label("Inter-compartmental clearance, Q, at 1 kg (L/h)")            # Merchan 2015 Table 3: 1.79 L/h/kg^0.75 (10% RSE)
    lvp <- log(13.00); label("Peripheral volume of distribution, V2, at 1 kg (L)")         # Merchan 2015 Table 3: 13.00 L/kg (12% RSE)

    # Allometric exponents held fixed at the paper's a priori values
    # (Merchan 2015 Methods 'Pharmacokinetic data analysis').
    allo_cl <- fixed(0.75); label("Allometric exponent on CL and Q (unitless)")            # Merchan 2015 Methods 'Pharmacokinetic data analysis'
    allo_v  <- fixed(1.0);  label("Allometric exponent on V1 and V2 (unitless)")           # Merchan 2015 Methods 'Pharmacokinetic data analysis'

    # IIV - Merchan 2015 reports ISV% as sqrt(omega^2) x 100 (Table 3
    # footnote a). Variances stored here are (ISV%/100)^2 directly; this is
    # NOT the same convention as %CV for a log-normal exponential model.
    etalcl ~ 0.3376  # Merchan 2015 Table 3: ISV 58.1% -> 0.581^2 = 0.3376
    etalvc ~ 0.6115  # Merchan 2015 Table 3: ISV 78.2% -> 0.782^2 = 0.6115
    etalq  ~ 0.4135  # Merchan 2015 Table 3: ISV 64.3% -> 0.643^2 = 0.4135
    etalvp ~ 0.6100  # Merchan 2015 Table 3: ISV 78.1% -> 0.781^2 = 0.6100

    # Proportional residual error - Merchan 2015 Methods 'Pharmacokinetic
    # data analysis' (Yij_obs = Yij_pred * (1 + eps)); Table 3 reports
    # residual error 28% (24% RSE).
    propSd <- 0.28; label("Proportional residual error (fraction)")                        # Merchan 2015 Table 3: 28% (24% RSE)
  })

  model({
    # Reference weight (Merchan 2015 Results: "a typical neonate weighing 1 kg").
    ref_wt <- 1

    # Individual structural parameters with allometric scaling on body weight.
    cl <- exp(lcl + etalcl) * (WT / ref_wt)^allo_cl
    vc <- exp(lvc + etalvc) * (WT / ref_wt)^allo_v
    q  <- exp(lq  + etalq)  * (WT / ref_wt)^allo_cl
    vp <- exp(lvp + etalvp) * (WT / ref_wt)^allo_v

    # Micro-constants for the explicit two-compartment ODE system.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system - intravenous dosing routes directly into the central
    # compartment (no depot); the 60-min infusion duration is supplied by
    # the user via the rate column of the event table.
    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
