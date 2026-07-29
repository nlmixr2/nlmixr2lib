Tremoulet_2014_ampicillin <- function() {
  description <- "One-compartment IV population PK model for ampicillin in preterm and term neonates (Tremoulet 2014; opportunistic POPS / PTN study). Clearance is allometrically scaled linearly to body weight and modulated by a serum-creatinine power factor (0.6/SCR)^0.428 and a postmenstrual-age power factor (PMA/37)^1.34. Central volume scales linearly with body weight (0.399 L/kg). Inter-individual variability is supported on CL only; residual variability is proportional."
  reference <- "Tremoulet A, Le J, Poindexter B, Sullivan JE, Laughon M, Delmore P, Salgado A, Chong SI, Melloni C, Gao J, Benjamin DK Jr, Capparelli EV, Cohen-Wolkowiez M; Administrative Core Committee of the Best Pharmaceuticals for Children Act-Pediatric Trials Network. Characterization of the population pharmacokinetics of ampicillin in neonates using an opportunistic study design. Antimicrob Agents Chemother. 2014;58(6):3013-3020. doi:10.1128/AAC.02374-13"
  vignette <- "Tremoulet_2014_ampicillin"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight (current; time-varying weight as measured at PK-sample visits)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Used for linear weight scaling on CL and Vc with no exponent (per-kg parameterisation in Tremoulet 2014: V = theta1 * WTKG; CL = theta2 * WTKG * ...). Source column WTKG.",
      source_name        = "WTKG"
    ),
    PAGE = list(
      description        = "Postmenstrual age (PMA = gestational age at birth + postnatal age)",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Drives the PMA power factor (PAGE/37)^e_page_cl on CL. The canonical PAGE convention in covariate-columns.md is months; this model declares weeks to match the Tremoulet 2014 published equation directly (reference 37 weeks). Same units precedent as Germovsek_2018_meropenem. Source column PMA.",
      source_name        = "PMA"
    ),
    CREAT = list(
      description        = "Serum creatinine concentration",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Used in the renal-function power factor (0.6/CREAT)^e_creat_cl on CL with reference 0.6 mg/dL. Source column SCR.",
      source_name        = "SCR"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 73L,
    n_studies      = 1L,
    age_range      = "Postnatal age 0-25 days (median 5); gestational age 24-41 weeks (median 36)",
    age_median     = "PNA 5 days; GA 36 weeks",
    weight_range   = "Body weight not tabulated in Tremoulet 2014 Table 1; cohort spans extremely premature (GA 24 weeks) to term (GA 41 weeks) neonates",
    weight_median  = "Not tabulated",
    sex_female_pct = 47.9,
    race_ethnicity = "Race: 77% White, 16% Black, 4% Other, 1% not reported. Ethnicity: 18% Hispanic/Latino, 77% non-Hispanic, 6% not reported (Tremoulet 2014 Table 1).",
    disease_state  = "Hospitalised neonates receiving ampicillin per standard of care in the neonatal intensive care unit; indications included presumed/confirmed infection, sepsis, necrotizing enterocolitis, abdominal procedure, meconium ileus with peritonitis, and pneumonia.",
    dose_range     = "Median 200 mg/kg/day (range 100-350); intervals 6, 8, or 12 h. Routes IV (one excluded subject received intramuscular). 11% q6h, 34% q8h, 55% q12h overall.",
    gestational_age_range = "24-41 weeks (median 36)",
    postnatal_age_range = "0-25 days (median 5)",
    postmenstrual_age_range = "Derived (PMA = GA + PNA/7); range spans approximately 24-45 weeks given the GA + PNA ranges above",
    samples_plasma = "142 plasma samples from 73 neonates (median 2.1 per subject; 17 (23%) of subjects contributed more than 2 samples)",
    regions        = "United States (9 NICU centres; Pediatric Trials Network sites under POPS protocol NICHD-2011-POP01)",
    notes          = "Demographics from Tremoulet 2014 Table 1. Original enrolment was 75 neonates; 2 excluded (one with a single below-quantitative-limit sample and PNA > 28 days; one with a suspected dose-history error) leaving 73 in the analysis. 14 of the 156 samples (9%) from the 73 retained subjects were excluded for assay or timing reasons. Median observed concentration 123 ug/mL (range 0.85-464). LLOQ 0.05 ug/mL; assay range 0.05-50 ug/mL by HPLC-MS/MS."
  )

  ini({
    # Structural parameters: per-kg form from Tremoulet 2014 Table 4 (final
    # model). Typical CL and Vc are reported on a per-kg basis so the
    # weight scaling enters as a simple multiplicative WT term in model().
    lcl <- log(0.078); label("Typical clearance per kg body weight at SCR = 0.6 mg/dL and PMA = 37 weeks (CL, L/h/kg)")  # Tremoulet 2014 Table 4: theta_CL = 0.078
    lvc <- log(0.399); label("Typical central volume of distribution per kg body weight (Vc, L/kg)")                       # Tremoulet 2014 Table 4: theta_V  = 0.399

    # Covariate effects on CL (Tremoulet 2014 Table 4).
    e_creat_cl <- 0.428; label("Exponent on the (0.6 / CREAT) ratio for CL (unitless)")  # Tremoulet 2014 Table 4: theta_CL,SCR = 0.428
    e_page_cl  <- 1.34;  label("Exponent on the (PAGE / 37) PMA ratio for CL (unitless)") # Tremoulet 2014 Table 4: theta_CL,PMA = 1.34

    # Inter-individual variability. Paper reports CV% = 22.8 for CL only
    # (Table 4). The log-normal mapping omega^2 = log(1 + CV^2) gives the
    # variance on the ETA scale.
    etalcl ~ 0.05068  # Tremoulet 2014 Table 4: omega^2(CL) reported as 22.8% CV -> log(1 + 0.228^2)

    # Residual error. Paper reports sigma^2 = 33.9 (CV%) (Table 4) and
    # explored proportional + additive models in Methods but retained
    # proportional only in the final model.
    propSd <- 0.339; label("Proportional residual error (fraction)")  # Tremoulet 2014 Table 4: sigma^2 = 33.9% CV
  })

  model({
    # Per-kg PK with renal-function and PMA-maturation covariates.
    cl <- exp(lcl + etalcl) * WT * (0.6 / CREAT)^e_creat_cl * (PAGE / 37)^e_page_cl
    vc <- exp(lvc) * WT

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Concentration: dose in mg, vc in L -> mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
