Ku_2018_diazepam <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous diazepam in children",
    "aged 3 months to 18 years treated for status epilepticus. Clearance,",
    "central volume, inter-compartmental clearance, and peripheral volume",
    "scale allometrically with total body weight referenced to a 70 kg adult",
    "(fixed exponents 0.75 on CL and Q; 1 on V1 and V2). IIV is estimated on",
    "CL and V1 only; IIV on Q and V2 was held fixed at 0 in the final model",
    "to avoid >50% shrinkage. Proportional residual error."
  )
  reference <- paste(
    "Ku LC, Hornik CP, Beechinor RJ, Chamberlain JM, Guptill JT, Harper B,",
    "Capparelli EV, Martz K, Anand R, Cohen-Wolkowiez M, Gonzalez D, on",
    "behalf of the Best Pharmaceuticals for Children Act - Pediatric Trials",
    "Network Steering Committee.",
    "Population Pharmacokinetics and Exploratory Exposure-Response",
    "Relationships of Diazepam in Children Treated for Status Epilepticus.",
    "CPT Pharmacometrics Syst Pharmacol. 2018;7(11):718-727.",
    "doi:10.1002/psp4.12349.",
    sep = " "
  )
  vignette <- "Ku_2018_diazepam"
  units    <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric power scaling with reference weight 70 kg and",
        "fixed theory-based exponents (Ku 2018 Methods / Results 'Population",
        "PK model development'): 0.75 on CL and Q, 1.0 on V1 and V2.",
        "Estimating the exponents resulted in high CL-IIV shrinkage (28.1% ->",
        "99.7%) and was rejected. Cohort median 15 kg, range 5-89 kg",
        "(Table 1)."
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 87L,
    n_studies       = 1L,
    n_observations  = 162L,
    age_range       = "0.4-17.8 years",
    age_median      = "3.9 years",
    age_groups      = "3 months to <3 years 45%; 3 to <13 years 44%; 13-18 years 11%",
    weight_range    = "5-89 kg",
    weight_median   = "15 kg",
    sex_female_pct  = 48,
    race_ethnicity  = c(White = 57, Hispanic_ethnicity = 38),
    disease_state   = paste(
      "Generalized convulsive status epilepticus presenting to the emergency",
      "department. Patients excluded for pregnancy, hypotensive shock,",
      "significant cardiac dysrhythmia, need for emergent surgery or general",
      "anesthesia, known contraindication to benzodiazepines, or",
      "benzodiazepine use within 7 days of presentation."
    ),
    dose_range      = paste(
      "Initial IV diazepam 0.2 mg/kg (maximum 8 mg) by slow push over 1",
      "minute. Patients with continued convulsions received a second dose of",
      "0.1 mg/kg (maximum 4 mg). 28% of subjects received a second dose."
    ),
    regions         = "United States (11 large academic pediatric hospitals)",
    notes           = paste(
      "Multicenter double-blind randomized clinical trial NCT00621478 /",
      "IND #79,010 comparing IV diazepam vs IV lorazepam for pediatric",
      "status epilepticus (Chamberlain 2014 JAMA). Patients enrolled under",
      "21 CFR 50.24 (Exception from Informed Consent for Emergency",
      "Research). Sparse PK sampling: up to 3 samples within 48 h of the",
      "first diazepam dose; median 2 samples per patient (range 1-4)."
    )
  )

  ini({
    # Final-model fixed-effect estimates from Ku 2018 Table 2 ("Final model
    # parameter estimates and bootstrap results"). Reference subject: 70 kg
    # adult body weight. Allometric exponents are theory-based and were
    # held fixed during estimation (Results section: estimating them
    # produced 99.7% CL-IIV shrinkage and was rejected; fixed-exponent
    # model "in line with previous literature" -- Anderson & Holford 2008).
    lcl <- log(2.36) ; label("Clearance at WT = 70 kg (L/h)")                                # Ku 2018 Table 2 (CL_70KG = 2.36 L/h, final)
    lvc <- log(42)   ; label("Central volume of distribution at WT = 70 kg (L)")             # Ku 2018 Table 2 (V1_70KG = 42 L, final)
    lq  <- log(22.6) ; label("Inter-compartmental clearance at WT = 70 kg (L/h)")            # Ku 2018 Table 2 (Q_70KG = 22.6 L/h, final)
    lvp <- log(56.5) ; label("Peripheral volume of distribution at WT = 70 kg (L)")          # Ku 2018 Table 2 (V2_70KG = 56.5 L, final)

    # Theory-based fixed allometric exponents (Ku 2018 Methods + Results).
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent of (WT/70) on CL (unitless; fixed)") # Ku 2018 Methods / Results
    e_wt_q  <- fixed(0.75) ; label("Allometric exponent of (WT/70) on Q (unitless; fixed)")  # Ku 2018 Methods / Results
    e_wt_vc <- fixed(1)    ; label("Allometric exponent of (WT/70) on V1 (unitless; fixed)") # Ku 2018 Methods / Results
    e_wt_vp <- fixed(1)    ; label("Allometric exponent of (WT/70) on V2 (unitless; fixed)") # Ku 2018 Methods / Results

    # IIV. Ku 2018 Table 2 footnote c states "IIV terms are shown as
    # variance"; the per-row CV% column in the source text (49.9% on CL,
    # 115% on V1) matches sqrt(variance) within rounding, confirming the
    # values are NONMEM omega^2 on the log-eta scale and may be passed
    # directly as eta variances in nlmixr2. IIV on Q and V2 was held fixed
    # at 0 in the final model after producing high shrinkage when
    # estimated (Results "Population PK model development").
    etalcl ~ 0.249           # Ku 2018 Table 2 (IIV CL variance = 0.249)
    etalvc ~ 1.31            # Ku 2018 Table 2 (IIV V1 variance = 1.31)
    etalq  ~ fixed(0)        # Ku 2018 Table 2 (IIV Q = 0 FIX)
    etalvp ~ fixed(0)        # Ku 2018 Table 2 (IIV V2 = 0 FIX)

    # Proportional residual error. Table 2 reports the proportional error
    # as 0.132 in NONMEM SIGMA-variance units (the additive component did
    # not improve fit and was dropped). Converted to nlmixr2 propSd
    # (linear-scale standard deviation): sqrt(0.132) = 0.3633.
    propSd <- 0.3633 ; label("Proportional residual SD (fraction)")  # Ku 2018 Table 2 (proportional variance 0.132 -> SD = sqrt(0.132))
  })
  model({
    # Body-weight scaling reference: 70 kg adult (Ku 2018 Results).
    wt70 <- WT / 70

    # Individual PK parameters. Allometric scaling on CL, Q (exponent 0.75)
    # and V1, V2 (exponent 1.0); IIV on CL and V1 only.
    cl <- exp(lcl + etalcl) * wt70^e_wt_cl
    vc <- exp(lvc + etalvc) * wt70^e_wt_vc
    q  <- exp(lq  + etalq)  * wt70^e_wt_q
    vp <- exp(lvp + etalvp) * wt70^e_wt_vp

    # Two-compartment IV disposition (no depot; doses load central via
    # cmt = central in the event table).
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Observed plasma diazepam concentration. Dose in mg, volume in L
    # yields mg/L = ug/mL; multiply by 1000 to convert to ng/mL.
    Cc <- central / vc * 1000

    Cc ~ prop(propSd)
  })
}
