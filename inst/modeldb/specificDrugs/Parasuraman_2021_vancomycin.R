Parasuraman_2021_vancomycin <- function() {
  description <- "Parallel one-compartment plasma and one-compartment CSF population PK model for intravenous and intraventricular vancomycin in extremely preterm infants (<28 weeks gestation) treated for ventriculitis (Parasuraman 2021). Plasma CL and V are allometrically scaled to 70 kg (fixed exponents 0.632 on CL, 1 on V) with a fixed Rhodin postmenstrual-age maturation sigmoid on CL; the CSF compartment receives intraventricular doses and has its own clearance and volume, with no plasma-CSF transfer."
  reference <- paste(
    "Parasuraman JM, Kloprogge F, Standing JF, Albur M, Heep A.",
    "Population pharmacokinetics of intraventricular vancomycin in neonatal ventriculitis, a preterm pilot study.",
    "Eur J Pharm Sci. 2021;158:105643. doi:10.1016/j.ejps.2020.105643.",
    "Size and maturation form from Germovsek E, Barker CIS, Sharland M, Standing JF.",
    "Scaling clearance in paediatric pharmacokinetics: all models are wrong, which are useful?",
    "Br J Clin Pharmacol. 2017;83(4):777-790. doi:10.1111/bcp.13160 (Table 1, model 9b).",
    sep = " "
  )
  vignette <- "Parasuraman_2021_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Current body weight at the time of dosing/sampling",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying. Allometric scaling of plasma CL (fixed exponent 0.632) and V (exponent 1) to a 70 kg reference (Methods 'NLME modelling'; Table 3 units L/h/70 kg and L/70 kg). The CSF parameters are NOT weight-scaled (Table 3 units L/h and L). Current weight was not tabulated; Table 2 reports birth weight only (median 0.78 kg, range 0.517-1.13 kg).",
      source_name = "WT"
    ),
    PAGE = list(
      description = "Postmenstrual age (gestational age + postnatal age)",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying. Drives the fixed postmenstrual-age maturation sigmoid on plasma CL, PAGE^3.33 / (55.4^3.33 + PAGE^3.33), which is the form Germovsek 2017 (Table 1 model 9b) pairs with the fixed 0.632 weight exponent. Weeks, per the neonatal maturation convention noted in the PAGE register entry. Table 2: median 34.4 weeks (range 30.2-48.1).",
      source_name = "PMA"
    )
  )

  covariatesDataExcluded <- list(
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on plasma CL; not retained (Results: 'None of the covariates provided a significant reduction in the OFV'). Table 2: median 25 umol/L (range 16-47).",
      source_name = "creatinine"
    ),
    CSF_TPRO = list(
      description = "Cerebrospinal-fluid total protein concentration",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CSF clearance, centred on the cohort median 3.97 g/L; not retained. Table 2: median 3.97 g/L (range 0.62-24.34).",
      source_name = "CSF protein"
    ),
    VI = list(
      description = "Cranial-ultrasound ventricular index (falx to lateral wall of the anterior horn, coronal plane)",
      units = "mm",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CSF clearance and CSF volume, centred on the cohort median 17.3 mm; not retained. Table 2: median 17.3 mm (range 14.1-34.6). Not a register covariate; documented here only because it was screened.",
      source_name = "VI"
    )
  )

  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE),
    csf = list(analyte = "vancomycin", units = "mg", specimen = "CSF", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 8L,
    n_studies = 1L,
    age_range = "Postnatal age 3.9-23.1 weeks; postmenstrual age 30.2-48.1 weeks",
    age_median = "Postnatal age 8.7 weeks; postmenstrual age 34.4 weeks",
    gestational_age_range = "23.9-27.7 weeks (median 25.3)",
    weight_range = "Birth weight 0.517-1.13 kg (current weight at sampling not reported)",
    weight_median = "Birth weight 0.78 kg",
    sex_female_pct = 62.5,
    race_ethnicity = "Not reported",
    disease_state = "Extremely preterm infants (<28 weeks gestation) with ventriculitis (CSF white cell count >20/mm^3 or positive CSF culture, mainly coagulase-negative staphylococci) and a ventricular access device (Ommaya reservoir) for CSF drainage.",
    dose_range = "Intraventricular vancomycin starting doses of 3, 5, 10 or 15 mg (10 mg/mL slow bolus over 2 min via the reservoir), repeated when the CSF level fell below 10 mg/L; 6/8 infants also received IV vancomycin 15 mg/kg every 24 h (<29 weeks PMA) or 12 h (29-35 weeks PMA).",
    regions = "United Kingdom (Southmead Hospital NICU, Bristol), 2009-2016",
    samples_plasma = "37 plasma concentrations from 5 infants (median 5 per infant, range 2-18)",
    samples_csf = "67 CSF concentrations from 8 infants (median 7.5 per infant, range 4-16)",
    notes = "Retrospective single-centre case review; demographics from Table 2. CSF drainage was fairly constant at 5-10 mL/kg/day (Discussion). One infant's second treatment course was excluded."
  )

  ini({
    # Plasma parameters (Table 3, final model), standardised to 70 kg.
    lcl <- log(9.29); label("Plasma clearance standardised to 70 kg, fully mature (L/h)") # Table 3: CL = 9.29 L/h/70 kg (RSE 5.2%)
    lvc <- log(54.4); label("Plasma volume of distribution standardised to 70 kg (L)") # Table 3: V = 54.4 L/70 kg (RSE 11.5%)

    # CSF parameters (Table 3, final model); not weight-scaled.
    lcl_csf <- log(0.002); label("Clearance from the CSF compartment (L/h)") # Table 3: CLCSF = 0.002 L/h (RSE 12.9%)
    lvcsf <- log(0.109); label("CSF compartment volume (L)") # Table 3: VCSF = 0.109 L (RSE 23.8%)

    # A priori size scaling (Methods 'NLME modelling': weight standardised to
    # 70 kg, allometric exponent 0.632); V linear in weight (Table 3 units L/70 kg).
    e_wt_cl <- fixed(0.632); label("Allometric exponent on plasma CL (unitless)") # Methods 'NLME modelling': exponent 0.632
    e_wt_vc <- fixed(1); label("Allometric exponent on plasma V (unitless)") # Table 3: V reported per 70 kg (linear)

    # A priori postmenstrual-age maturation of CL (Methods cite Germovsek 2016,
    # 2017). The paper does not print the constants; Germovsek 2017 Table 1
    # model 9b pairs the fixed 0.632 exponent with TM50 = 55.4 weeks and
    # Hill = 3.33 (Rhodin 2009 GFR fit with the estimated 0.632 exponent), the
    # same pairing as the Germovsek 2016 gentamicin model the paper cites.
    tmat50 <- fixed(55.4); label("Postmenstrual age at 50% CL maturation (weeks)") # Germovsek 2017 Table 1 model 9b: 55.4 (fixed)
    hill <- fixed(3.33); label("Hill coefficient of the CL maturation sigmoid (unitless)") # Germovsek 2017 Table 1 model 9b: 3.33 (fixed)

    # IIV. Table 3 footnote: 'CV is coefficient of variation (100xOmega)', so
    # omega = CV/100 and the variance is (CV/100)^2. IIV on V was 0 FIX and
    # is omitted.
    etalcl ~ 0.0676 # Table 3: IIV CL = 26% CV -> 0.26^2
    etalcl_csf ~ 0.0841 # Table 3: IIV CLCSF = 29% CV -> 0.29^2
    etalvcsf ~ 0.36 # Table 3: IIV VCSF = 60% CV -> 0.60^2

    # Residual error: proportional on each matrix (Results: 'proportional
    # error on plasma and CSF data').
    propSd <- 0.3131; label("Plasma proportional residual error (fraction)") # Table 3: sigmaprop plasma = 31.31% CV
    propSd_Ccsf <- 0.458; label("CSF proportional residual error (fraction)") # Table 3: sigmaprop CSF = 45.8% CV
  })

  model({
    # 1. Postmenstrual-age maturation of plasma CL (Germovsek 2017 model 9b form)
    fmat <- PAGE^hill / (tmat50^hill + PAGE^hill)

    # 2. Individual parameters
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * fmat
    vc <- exp(lvc) * (WT / 70)^e_wt_vc
    cl_csf <- exp(lcl_csf + etalcl_csf)
    vcsf <- exp(lvcsf + etalvcsf)

    kel <- cl / vc
    kel_csf <- cl_csf / vcsf

    # 3. Parallel one-compartment models with no plasma-CSF transfer
    # (Methods equations dA(p)/dt and dA(csf)/dt; Results: transfer in either
    # direction did not improve the fit). IV doses go to central,
    # intraventricular doses to csf.
    d/dt(central) <- -kel * central
    d/dt(csf) <- -kel_csf * csf

    # 4. Observations (mg / L = mg/L)
    Cc <- central / vc
    Ccsf <- csf / vcsf

    Cc ~ prop(propSd)
    Ccsf ~ prop(propSd_Ccsf)
  })
}
