Leroux_2016_cefotaxime <- function() {
  description <- "Two-compartment IV population PK model for cefotaxime in neonates and young infants (Leroux 2016). Clearance, central volume, peripheral volume, and inter-compartmental clearance are allometrically scaled to current body weight (fixed exponents 0.75 on CL and Q, 1.0 on V1 and V2; reference weight 1.665 kg). Clearance carries a power-form maturation function on gestational age (reference 30 weeks) and postnatal age (reference 12 days). Only CL has inter-individual variability; residual error is proportional."
  reference <- "Leroux S, Roue J-M, Gouyon J-B, Biran V, Zheng H, Zhao W, Jacqz-Aigrain E. A Population and Developmental Pharmacokinetic Analysis To Evaluate and Optimize Cefotaxime Dosing Regimen in Neonates and Young Infants. Antimicrob Agents Chemother. 2016;60(11):6626-6634. doi:10.1128/AAC.01045-16"
  vignette <- "Leroux_2016_cefotaxime"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Current body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Used in allometric scaling (WT / 1.665)^exponent on CL, V1, V2 and Q. Reference weight 1.665 kg from Leroux 2016 Table 4 (close to the cohort median current weight of 1.6475 kg).",
      source_name        = "CW"
    ),
    GA = list(
      description        = "Gestational age at birth",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Enters the cefotaxime CL maturation factor FGA = (GA / 30)^2.27 (Leroux 2016 Table 4). Cohort range 23.0 to 42.0 weeks (Table 2).",
      source_name        = "GA"
    ),
    PNA = list(
      description        = "Postnatal age (chronological since birth)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Leroux 2016 Table 4 expresses PNA in days with reference 12 days; the canonical PNA carries months, so the in-model term is reparameterised as FPNA = (PNA / 0.3943)^0.28 with reference 0.3943 months = 12 / 30.4375 days. Cohort range 0 to 69 days (Table 2).",
      source_name        = "PNA"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 100L,
    n_studies       = 1L,
    age_range       = "GA 23.0-42.0 weeks at birth; PNA 0-69 days; PMA 25-44 weeks",
    age_median      = "GA 31.5 weeks; PNA 9 days; PMA 33 weeks",
    weight_range    = "current weight 0.530-4.200 kg; birth weight 0.512-3.990 kg",
    weight_median   = "current weight 1.6475 kg; birth weight 1.415 kg",
    sex_female_pct  = 60.0,
    race_ethnicity  = "Not reported (three French neonatal intensive care units)",
    disease_state   = "Neonates and young infants (postmenstrual age <= 44 weeks) receiving intravenous cefotaxime as part of routine clinical care for suspected or proven neonatal sepsis. Twenty-five subjects had positive blood cultures (2 early-onset: Streptococcus agalactiae, Escherichia coli; 23 late-onset, predominantly coagulase-negative staphylococci). 53 of 100 subjects required invasive ventilation; 21 received vasopressors; 74 received an aminoglycoside.",
    dose_range      = "Cefotaxime 50 mg/kg per dose (mean 47.7 mg/kg, SD 8.2) given as a direct IV injection or 15-30 min IV infusion two, three, or four times daily; one patient received a single 100 mg/kg dose.",
    regions         = "France (multicentre: Robert Debre Paris, Brest, Saint Pierre de la Reunion University Hospitals).",
    n_observations  = 185L,
    samples_per_subject = "median 1.0 (range 1-6); mean 1.8 samples per patient",
    notes           = "Open-label opportunistic-sampling popPK study; cefotaxime concentrations measured by HPLC-MS/MS with LLOQ 0.05 mg/L. Demographics from Table 2; pharmacokinetic estimates from Table 4. Forward and backward covariate selection retained CW (allometric, exponents fixed), GA, and PNA on CL; serum creatinine did not survive backward elimination."
  )

  ini({
    # Structural parameters at the reference covariate values (CW = 1.665 kg,
    # GA = 30 weeks, PNA = 12 days; FGA = FPNA = 1 at these references).
    # Values from Leroux 2016 Table 4, final-model column.
    lcl <- log(0.21); label("Clearance at reference covariates (L/h)")                              # Leroux 2016 Table 4: theta_1 = 0.21 L/h
    lvc <- log(0.71); label("Central volume of distribution at reference weight (L)")               # Leroux 2016 Table 4: theta_2 = 0.71 L
    lvp <- log(0.35); label("Peripheral volume of distribution at reference weight (L)")            # Leroux 2016 Table 4: theta_3 = 0.35 L
    lq  <- log(0.27); label("Inter-compartmental clearance at reference weight (L/h)")              # Leroux 2016 Table 4: theta_4 = 0.27 L/h

    # Allometric exponents (Leroux 2016 Results 'Model building': "Allometric
    # exponents of 0.75 and 1 were fixed for CL and V, respectively (19);
    # their estimation by use of the model ... did not significantly improve
    # the fit of the data."). Applied to both CL and Q (clearance terms) and
    # to both V1 and V2 (volume terms).
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL (unitless)")  # Leroux 2016 Results: fixed at 0.75
    e_wt_q  <- fixed(0.75); label("Allometric exponent on Q (unitless)")   # Leroux 2016 Results: fixed at 0.75
    e_wt_vc <- fixed(1.00); label("Allometric exponent on V1 (unitless)")  # Leroux 2016 Results: fixed at 1
    e_wt_vp <- fixed(1.00); label("Allometric exponent on V2 (unitless)")  # Leroux 2016 Results: fixed at 1

    # Maturation effects on CL (Leroux 2016 Table 4 final model). Both
    # exponents estimated (not fixed). FGA = (GA / 30)^theta_5,
    # FPNA = (PNA_days / 12)^theta_6. The PNA term is reparameterised
    # in model() to use the canonical PNA in months (1.665 kg reference
    # weight stays as-is for WT in kg).
    e_ga_cl  <- 2.27; label("Gestational-age exponent on CL (unitless); FGA = (GA / 30)^e_ga_cl")  # Leroux 2016 Table 4: theta_5 = 2.27
    e_pna_cl <- 0.28; label("Postnatal-age exponent on CL (unitless); FPNA = (PNA_months / 0.3943)^e_pna_cl") # Leroux 2016 Table 4: theta_6 = 0.28

    # Inter-individual variability. Leroux 2016 Table 4 reports "IIV (%) for
    # CL = 21.0" (CV%); under the exponential eta model documented in
    # Methods 'Step 1: model building', the corresponding variance on the
    # log-normal eta scale is log(1 + 0.21^2) = 0.04316. Only CL carried IIV
    # in the final model ("Interindividual variability ... could be
    # estimated only for CL.", Results 'Model building').
    etalcl ~ 0.04316  # Leroux 2016 Table 4: IIV CL = 21.0 % CV -> log(1 + 0.210^2) = 0.04316

    # Residual error. Leroux 2016 Results 'Model building': "A proportional
    # model best described the residual unexplained variability." Table 4
    # reports residual variability 35.1 % which translates directly to a
    # proportional SD of 0.351 in nlmixr2's prop() form.
    propSd <- 0.351; label("Proportional residual error (fraction)")  # Leroux 2016 Table 4: residual variability = 35.1 %
  })

  model({
    # Maturation factors on CL (Leroux 2016 Table 4).
    # GA: paper uses weeks with reference 30 weeks (canonical GA is weeks).
    # PNA: paper uses days with reference 12 days; canonical PNA is months,
    #   so the reference is converted to 12 / 30.4375 = 0.3943 months and
    #   the input column is expected to carry PNA in months.
    fga  <- (GA / 30)^e_ga_cl
    fpna <- (PNA / 0.3943)^e_pna_cl

    # Individual PK parameters. WT (current weight) in kg with reference
    # 1.665 kg = 1665 g per Leroux 2016 Table 4.
    cl <- exp(lcl + etalcl) * (WT / 1.665)^e_wt_cl * fga * fpna
    vc <- exp(lvc)          * (WT / 1.665)^e_wt_vc
    vp <- exp(lvp)          * (WT / 1.665)^e_wt_vp
    q  <- exp(lq)           * (WT / 1.665)^e_wt_q

    # Two-compartment IV; cefotaxime is given as a direct IV injection or
    # short (15-30 min) IV infusion into the central compartment.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Observation: cefotaxime plasma concentration in mg/L (dose in mg, V in L).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
