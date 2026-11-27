Hassan_2020_meropenem <- function() {
  description <- "Two-compartment IV-infusion population PK model for meropenem in infants and children from birth to 17 years (Hassan 2020), fitted to pooled observed neonatal data (Smith 2011) and simulated older-child data (replicating Du 2006). All four disposition parameters scale with body weight normalised to 70 kg (one estimated exponent shared by CL and Q, another by V1 and V2); CL and Q share a sigmoid maturation function of postnatal age; CL also carries power effects of serum creatinine and of gestational age normalised to 40 weeks. Exponential IIV on CL, V1, Q and V2; proportional residual error."
  reference <- "Hassan HE, Ivaturi V, Gobburu J, Green TP. Dosage Regimens for Meropenem in Children with Pseudomonas Infections Do Not Meet Serum Concentration Targets. Clin Transl Sci. 2020;13(2):301-308. doi:10.1111/cts.12710"
  vignette <- "Hassan_2020_meropenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects on CL and Q (shared exponent dCLdWT) and on V1 and V2 (shared exponent dVdWT), all normalised to 70 kg (final-model equations, Results). Pooled data set: median 2.36 kg, range 0.39-88.3 kg (Table 2).",
      source_name = "WT"
    ),
    PNA = list(
      description = "Postnatal age",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = "Hassan 2020 carries postnatal age in YEARS (Table 2 'PNA, years'; Table 3 Age50 unit 'y'). model() converts the canonical months to years with PNA / 12 before evaluating the maturation term PNA^gamma / (PNA^gamma + Age50) on CL and Q. Pooled data set: median 0.11 years, range 0.002-16.3 years (Table 2).",
      source_name = "PNA"
    ),
    GA = list(
      description = "Gestational age at birth",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL, normalised to 40 weeks: (GA/40)^dCLdGA. Missing gestational age for infants and children older than 120 days was imputed as 40 weeks (Methods 'Data'), so for older children supply GA = 40 when unknown. Pooled data set: median 33 weeks, range 22.5-40 (Table 2).",
      source_name = "GA"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Un-normalised power effect on CL: SCR^dCLdSCR with SCR in mg/dL (final-model CL equation; there is no reference value in the printed equation, so the typical CL is the value at SCR = 1 mg/dL). Pooled data set: median 0.50 mg/dL, range 0.1-1.9 (Table 2).",
      source_name = "SCR"
    )
  )

  compartmentData <- list(
    central = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = FALSE)
  )

  population <- list(
    species = "human",
    n_subjects = 288L,
    n_studies = 4L,
    n_concentrations = 1482L,
    age_range = "Postnatal age 0.002-16.3 years (median 0.11 years)",
    age_median = "0.11 years",
    ga_range = "22.5-40 weeks (median 33)",
    weight_range = "0.39-88.3 kg (median 2.36 kg)",
    weight_median = "2.36 kg",
    renal_function = "Serum creatinine 0.1-1.9 mg/dL (median 0.50); subjects selected to have normal renal function (Discussion)",
    disease_state = "Infants with suspected or complicated intra-abdominal infection (Smith 2011) and infants and children with bacterial infection (Blumer 1995, Parker 1995, Du 2006)",
    dose_range = "Intravenous meropenem; Smith 2011 infants received weight- and age-based dosing; the simulated Du 2006 cohort received 20 mg/kg every 8 hours infused over 30 minutes",
    regions = "United States (Pediatric Trials Network data) plus the published Du 2006 / Blumer 1995 / Parker 1995 cohorts",
    notes = "Pooled data set of 188 observed infants (23-40 weeks GA, 1-92 days PNA; Smith 2011 via the NICHD Pediatric Trials Network repository) and 100 virtual subjects (2-200 months of age) simulated by the authors from the published Du 2006 model to stand in for the 99 children of Du 2006 / Blumer 1995 / Parker 1995 whose data were not available (Methods 'Data'; Table 2). Phoenix NLME 7.0, FOCE-ELS."
  )

  ini({
    # Structural parameters (Hassan 2020 Table 3, 'Parametric Estimate (SE)').
    # Typical values are for a 70-kg subject with the maturation term,
    # SCR = 1 mg/dL and GA = 40 weeks as printed in the final-model equations.
    lcl <- log(37.1)
    label("Clearance for a 70-kg subject at SCR 1 mg/dL, GA 40 weeks, before maturation (L/h)") # Table 3: TVCL = 37.1 L/hour/70 kg (SE 3.8)
    lvc <- log(29.1)
    label("Central volume for a 70-kg subject (L)") # Table 3: TVV = 29.1 L/70 kg (SE 2.9)
    lq <- log(6.10)
    label("Intercompartmental clearance for a 70-kg subject, before maturation (L/h)") # Table 3: TVCL2 = 6.10 L/hour/70 kg (SE 1.1)
    lvp <- log(31.9)
    label("Peripheral volume for a 70-kg subject (L)") # Table 3: TVV2 = 31.9 L/70 kg (SE 4.3)

    # Body-weight exponents (Table 3). One exponent is shared by CL and Q,
    # another by V1 and V2 (Results 'PopPK model development and qualification').
    e_wt_cl_q <- 1.20
    label("Power exponent of WT/70 on CL and Q (unitless)") # Table 3: dCLdWT = 1.20 (SE 0.05)
    e_wt_vc_vp <- 0.97
    label("Power exponent of WT/70 on V1 and V2 (unitless)") # Table 3: dVdWT = 0.97 (SE 0.03)

    # Postnatal-age maturation on CL and Q: PNA^gamma / (PNA^gamma + Age50),
    # PNA in years (final-model CL and CL2 equations).
    lpna50 <- log(0.39)
    label("Age50 of the postnatal-age maturation function on CL and Q (years)") # Table 3: Age 50 = 0.39 y (SE 0.05)
    lhill <- log(0.27)
    label("Hill coefficient of the postnatal-age maturation function (unitless)") # Table 3: gamma = 0.27 (SE 0.04)

    # Renal-function and gestational-age effects on CL (Table 3).
    e_creat_cl <- -0.25
    label("Power exponent of serum creatinine (mg/dL) on CL (unitless)") # Table 3: dCLdSCR = -0.25 (SE 0.03)
    e_ga_cl <- 0.84
    label("Power exponent of GA/40 on CL (unitless)") # Table 3: dCLdGA = 0.84 (SE 0.16)

    # IIV (Table 3 'Intersubject variability'; exponential model, Appendix S1
    # 'Variability models'). Values are read as omega^2 variances, the
    # quantity Phoenix NLME reports.
    etalcl ~ 0.315 # Table 3: CL = 0.315 (SE 0.220), shrinkage 10%
    etalvc ~ 0.186 # Table 3: V = 0.186 (SE 0.095), shrinkage 58%
    etalq ~ 1.14 # Table 3: CL2 = 1.14 (SE 0.56), shrinkage 39%
    etalvp ~ 0.319 # Table 3: V2 = 0.319 (SE 0.10), shrinkage 76%

    # Residual error: CObs = C * (1 + eps) (final-model equations).
    propSd <- 0.38
    label("Proportional residual error (fraction)") # Table 3: Proportional = 0.38 (SE 0.01)
  })
  model({
    # Postnatal age in years (canonical PNA is in months).
    pna_yr <- PNA / 12
    pna50 <- exp(lpna50)
    hill <- exp(lhill)
    pnah <- pna_yr^hill
    # Maturation as printed in the final-model equations and Appendix S1:
    # Age50 enters the denominator un-exponentiated.
    f_mat <- pnah / (pnah + pna50)

    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q * f_mat * CREAT^e_creat_cl * (GA / 40)^e_ga_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q <- exp(lq + etalq) * (WT / 70)^e_wt_cl_q * f_mat
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
