Colombo_2006_atazanavir <- function() {
  description <- "One-compartment first-order-absorption population PK model with absorption lag-time for orally administered atazanavir in HIV-1 infected adults; binary low-dose ritonavir (RTV) coadministration reduces apparent oral clearance by 46% (Colombo 2006)."
  reference <- "Colombo S, Buclin T, Cavassini M, Decosterd LA, Telenti A, Biollaz J, Csajka C. Population pharmacokinetics of atazanavir in patients with human immunodeficiency virus infection. Antimicrob Agents Chemother. 2006;50(11):3801-3808. doi:10.1128/aac.00098-06"
  vignette <- "Colombo_2006_atazanavir"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CONMED_RTV = list(
      description        = "Concomitant low-dose ritonavir (RTV) indicator (pharmacokinetic booster)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (unboosted atazanavir regimen)",
      notes              = "1 = subject receives low-dose ritonavir (typically 100 mg q.d.) as a pharmacokinetic booster of atazanavir; 0 = unboosted ATV. Drives a 46% reduction in apparent oral clearance via the linear-deviation form cl = exp(lcl) * (1 + e_rtv_cl * CONMED_RTV) per Colombo 2006 Table 2 / Table 3 (theta_b = -0.46, RSE 18%). 167 of 214 model-building subjects (78%) were on RTV at 100 mg q.d.",
      source_name        = "RTV"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Median 69 kg (range 43-117). Screened on CL in Table 2 with theta_b = 0.2 per relative body-weight deviation; delta-OF = -2.2 (below the 3.84 significance threshold for one additional parameter). Not retained in the final model."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Median 42 years (range 19-78). Screened on CL in Table 2 with theta_b = -0.3 per relative age deviation; delta-OF = -2.6. Not retained in the final model."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = "60 of 214 (28%) female. Screened in Table 2 as male = 1 with theta_b = 0.038; delta-OF = -0.1. Not retained. Canonical SEXF (1 = female) inverts the source paper's male = 1 encoding.",
      source_name = "SEX"
    ),
    RACE_WHITE = list(
      description = "Caucasian (White) race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "183 of 214 (86%) Caucasian. Screened in Table 2 as a Caucasian-vs-non-Caucasian contrast with theta_b = 6.78; delta-OF = -3.0. Not retained.",
      source_name = "CAUCASIAN"
    ),
    RACE_BLACK = list(
      description = "African (Black) race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "20 of 214 (9%) African. Screened in Table 2 as an African-vs-non-African contrast with theta_b = 8.4; delta-OF = -3.1. Not retained.",
      source_name = "AFRICAN"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "7 of 214 (3%) Asian. Screened in Table 2 with theta_b = 4.8; delta-OF = -1.5. Not retained.",
      source_name = "ASIAN"
    ),
    RACE_HISPANIC = list(
      description = "Hispanic ethnicity indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "4 of 214 (2%) Hispanic. Screened in Table 2 with theta_b = 12.9; delta-OF = -1.9. Not retained.",
      source_name = "HISPANIC"
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Derived from baseline serum creatinine (median 83 mmol/L, range 22-165). Screened on CL in Table 2 with theta_b = -0.2 per relative deviation; delta-OF = -1.7. Not retained.",
      source_name = "CLCR"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 214L,
    n_studies      = 2L,
    n_observations = 574L,
    age_range      = "19-78 years",
    age_median     = "42 years",
    weight_range   = "43-117 kg",
    weight_median  = "69 kg",
    sex_female_pct = 28,
    race_ethnicity = c(Caucasian = 86, African = 9, Asian = 3, Hispanic = 2),
    disease_state  = "HIV-1 infected adults under routine antiretroviral therapy including atazanavir, either alone or boosted with low-dose ritonavir.",
    dose_range     = "Oral atazanavir 300 mg q.d. boosted with ritonavir 100 mg q.d. (n=167; 78%) or unboosted atazanavir 400 mg q.d. (n=47; 22%); always combined with other antiretroviral agents.",
    regions        = "Switzerland (Swiss HIV Cohort Study, Lausanne center)",
    notes          = "Model built from 574 plasma samples: 346 sparse routine-TDM samples from 201 patients (median 1 sample per subject; range 1-6) plus 228 intensive-PK samples from a 13-patient drug-interaction substudy at steady state (12-24 timepoints per subject, 12 timepoints from pre-dose to 24 h post-dose). All subjects had been on the regimen for at least 1 month before sampling. Observed concentrations 50-6680 ng/mL. Model fit with NONMEM V (FOCE). External validation set: 112 sparse observations from 78 additional patients (model-validation cohort). Typical median elimination half-life is 4.6 h without RTV and 8.8 h with RTV per the paper Discussion."
  )

  ini({
    # Structural parameters (Table 3 final-model column).
    # CL/F is the typical apparent oral clearance in the absence of ritonavir
    # (the CONMED_RTV = 0 reference stratum); the CONMED_RTV covariate effect
    # reduces this by 46% when CONMED_RTV = 1.
    lcl <- log(12.9); label("Apparent oral clearance at CONMED_RTV = 0 (CL/F, L/h)")  # Table 3: CL/F = 12.9 L/h (RSE 17%)
    lvc <- log(88.3); label("Apparent volume of distribution (V/F, L)")               # Table 3: V/F = 88.3 L (RSE 9.5%)

    # ka and lag time were fixed to the rich-data substudy estimates because the
    # sparse cohort had too few early-phase observations to estimate them.
    # Results page 3804: "the population mean estimate and variance of ka and
    # the mean estimate of the lag time were fixed to the final values obtained
    # from the rich data set."
    lka   <- fixed(log(0.405)); label("First-order absorption rate constant (ka, 1/h)")  # Table 3: ka = 0.405 1/h, fixed (rich-data substudy)
    ltlag <- fixed(log(0.876)); label("Absorption lag-time (tlag, h)")                   # Results p3804: lag = 0.876 h, fixed mean (Table 3 reports rounded 0.88 h, RSE 10.3%)

    # F: relative bioavailability in the sparse routine-TDM cohort, accounting
    # for undercompliance vs the rich-data substudy (F_rich was set to 1
    # because no IV reference was available, Table 3 footnote f). The sparse
    # F is the one used by the authors in their dosage-adaptation simulations
    # (Discussion p3805) and is the appropriate value for forward simulation.
    lfdepot <- log(0.81); label("Relative bioavailability (F)")  # Table 3: F_sparse = 0.81 (RSE 75%)

    # Covariate effect of CONMED_RTV on CL/F (linear-deviation form per Table 2)
    e_rtv_cl <- -0.46; label("CONMED_RTV effect on CL/F (fractional change when RTV coadministered)")  # Table 3: theta_ritonavir = -0.46 (RSE 18.0%)

    # Inter-individual variability (omega^2 = log(CV^2 + 1))
    etalcl     ~ 0.065413        # Table 3: CV(CL/F) = 26% (RSE 56%); log(1 + 0.26^2) = 0.065413
    etalvc     ~ 0.080750        # Table 3: CV(V/F)  = 29% (RSE 80%); log(1 + 0.29^2) = 0.080750
    etalka     ~ fixed(0.911640) # Table 3: CV(ka)   = 122%, fixed (rich-data substudy); log(1 + 1.22^2) = 0.911640
    etalfdepot ~ 0.184403        # Table 3: CV(F)    = 45% (RSE 49%); log(1 + 0.45^2) = 0.184403

    # Residual error (sparse routine-TDM cohort, Table 3): combined proportional
    # + additive. The rich-data substudy uses a separate residual-error pair
    # (CV 19%, SD 0.370 mg/L) documented in the vignette Assumptions section.
    propSd <- 0.30;  label("Proportional residual error (CV, fraction)")  # Table 3 sparse: 30% (RSE 35%)
    addSd  <- 0.542; label("Additive residual error (mg/L)")              # Table 3 sparse: SD = +/-542 ng/mL = 0.542 mg/L
  })

  model({
    # Individual PK parameters; CONMED_RTV enters linearly on CL (Table 2 form)
    ka     <- exp(lka + etalka)
    cl     <- exp(lcl + etalcl) * (1 + e_rtv_cl * CONMED_RTV)
    vc     <- exp(lvc + etalvc)
    tlag   <- exp(ltlag)
    fdepot <- exp(lfdepot + etalfdepot)

    # Micro-constant
    kel <- cl / vc

    # ODE system: one-compartment with first-order absorption from depot
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Bioavailability and absorption lag-time on the depot compartment
    f(depot)    <- fdepot
    alag(depot) <- tlag

    # Observation: dose in mg, V in L -> central / vc gives mg/L (= ug/mL).
    # Combined proportional + additive residual error.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
