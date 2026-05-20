Nemoto_2017_ethanol <- function() {
  description <- "Bayesian population PK model for orally ingested ethanol (alcohol) in 34 healthy Japanese adults (Nemoto 2017). One-compartment model with first-order absorption and Michaelis-Menten elimination; covariates: sex, age, body weight, ALDH2 and ADH1B genotypes. Final model fit by a fully conditional MCMC Bayesian analysis with informative priors derived from Seng et al. 2014 (Chinese + Indian cohort)."
  reference   <- "Nemoto A, Masaaki M, Yamaoka K. A Bayesian Approach for Population Pharmacokinetic Modeling of Alcohol in Japanese Individuals. Curr Ther Res Clin Exp. 2017;85:1-7. doi:10.1016/j.curtheres.2017.04.001"
  vignette    <- "Nemoto_2017_ethanol"
  units       <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form scaling with reference 61.3 kg per Nemoto 2017 Table II (structural model row for the final model).",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form scaling with reference 29.4 y (the median age of the analysis cohort) per Nemoto 2017 Results page 4 and Table II.",
      source_name        = "age"
    ),
    SEXF = list(
      description        = "Female sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Additive shift on ka (-1.3 1/h) and on Vd/F (-12.2 L) when SEXF = 1, per Nemoto 2017 Table II final-model structural equations.",
      source_name        = "FEMALE"
    ),
    ALDH2_S2_CARRIER = list(
      description        = "Carrier of at least one ALDH2*2 inactive variant allele (1 = ALDH2*1/*2 or *2/*2; 0 = ALDH2*1/*1 wild-type)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ALDH2*1/*1 wild-type)",
      notes              = "Additive shift on Vd/F (-20.4 L) when ALDH2_S2_CARRIER = 1 per Nemoto 2017 Table II final model. The analysis cohort comprised 21/34 (62%) ALDH2*1/*1 wild-type and 13/34 (38%) ALDH2*1/*2 heterozygous subjects; no ALDH2*2/*2 homozygotes were enrolled, so only the heterozygote-vs-wild-type contrast is informed by data. Applying the indicator to ALDH2*2/*2 homozygotes is an untested extrapolation. Time-fixed per subject (germline genotype).",
      source_name        = "ALDH2"
    ),
    ADH1B_S2_HOM = list(
      description        = "Homozygous ADH1B*2/*2 indicator (1 = *2/*2; 0 = *2/*1 heterozygous)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ADH1B*2/*1 heterozygous)",
      notes              = "Switch between two Vmax values per Nemoto 2017 Table II structural model: Vmax = 7790 mg/h for *2/*1 carriers (reference) vs 7966 mg/h for *2/*2 homozygotes. The Nemoto 2017 Japanese cohort distribution across ADH1B genotypes is not tabulated in the paper; the structural model is silent on ADH1B*1/*1 wild-type subjects, who would conventionally be assigned ADH1B_S2_HOM = 0 (the heterozygous Vmax value) when applying the model. Time-fixed per subject (germline genotype).",
      source_name        = "ADH1B"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 34L,
    n_observations  = 157L,
    n_studies       = 1L,
    age_range       = "20-62 years",
    age_median      = "29.4 years (centring value used in the model)",
    weight_range    = "not reported; reference WT 61.3 kg used as centring value",
    weight_median   = "61.3 kg (centring value used in the model)",
    sex_female_pct  = 38.2,
    race_ethnicity  = "Japanese (single-cohort)",
    disease_state   = "Healthy adult volunteers",
    dose_range      = "14 g ethanol (350 mL beer) ingested over 10 minutes in the fasted state",
    regions         = "Japan (Osaka University cohort)",
    notes           = "ALDH2 distribution: 21/34 (62%) *1/*1, 13/34 (38%) *1/*2, 0/34 (0%) *2/*2 (Nemoto 2017 Methods 'Dataset', page 2). ADH1B genotype distribution within the Nemoto 2017 cohort is not tabulated; the Vmax structural model uses two values, one per *2/*1 vs *2/*2 ADH1B genotype, with no value provided for ADH1B*1/*1. Blood alcohol concentration observations sampled at 5, 10, 20, 30, and 60 minutes post-dose (early-absorption and around the peak); no elimination-phase data available (Nemoto 2017 Table I and Methods). Posterior estimates from a fully conditional MCMC Bayesian analysis with priors derived from Seng et al. 2014 popPK of ethanol in Chinese + Indian subjects (Nemoto 2017 'Prior information' and Table II)."
  )

  ini({
    # Structural parameters - Nemoto 2017 Table II final-model column.
    # Reference covariates: 29.4-year-old male, 61.3 kg, ALDH2*1/*1 wild-type, ADH1B*2/*1.
    lka   <- log(3.0);   label("Absorption rate ka for a 29.4-year-old male (1/h)")          # Table II final-model row: ka = 3.0 (95% CI 2.4 to 3.9)
    lvc   <- log(49.3);  label("Apparent Vd/F for a 29.4 y, 61.3 kg, ALDH2*1/*1 male (L)")    # Table II final-model row: Vd/F = 49.3 (47.4 to 51.2)
    lvmax <- log(7790);  label("Vmax for ADH1B*2/*1 carriers at 61.3 kg reference (mg/h)")    # Table II final-model row: Vmax_ADH1B*2/*1 = 7790 (7403 to 8264)
    lkm   <- log(0.074); label("Michaelis-Menten constant Km (mg/L)")                          # Table II final-model row: Km = 0.074 (0.001 to 0.391)

    # Covariate effects - additive on linear scale (Nemoto 2017 Table II structural model).
    e_sexf_ka             <- -1.3;  label("Additive shift on ka in females (1/h)")               # Table II final-model: theta_SEX(ka) = -1.3 (-2.1 to -0.56)
    e_sexf_vc             <- -12.2; label("Additive shift on Vd/F in females (L)")               # Table II final-model: theta_SEX(Vd/F) = -12.2 (-15.0 to -9.41)
    e_aldh2_s2_carrier_vc <- -20.4; label("Additive shift on Vd/F for ALDH2*2 carriers (L)")     # Table II final-model: theta_ALDH2(Vd/F) = -20.4 (-27.7 to -10.9)
    e_adh1b_s2_hom_vmax   <- 176;   label("Additive shift on Vmax for ADH1B*2/*2 homozygotes (mg/h)") # Table II final-model: 7966 - 7790 = 176 (difference between the two reported Vmax thetas)

    # Power-form covariate exponents (Nemoto 2017 Table II structural model).
    e_age_ka  <- 2.7;   label("Power exponent of (AGE/29.4) on ka (unitless)")     # Table II final-model: theta_AGE(ka) = 2.7 (2.1 to 3.4)
    e_age_vc  <- 0.52;  label("Power exponent of (AGE/29.4) on Vd/F (unitless)")   # Table II final-model: theta_AGE(Vd/F) = 0.52 (0.19 to 0.83)
    e_wt_vc   <- 0.78;  label("Power exponent of (WT/61.3) on Vd/F (unitless)")    # Table II final-model: theta_WT(Vd/F) = 0.78 (0.60 to 0.95)
    e_wt_vmax <- 0.78;  label("Power exponent of (WT/61.3) on Vmax (unitless)")    # Table II final-model: theta_WT(Vmax) = 0.78 (0.66 to 0.90)

    # IIV - log-normal between-subject variance from Nemoto 2017 Table II final-model.
    etalka   ~ 0.37   # omega^2_ka    = 0.37  (0.24 to 0.55)
    etalvc   ~ 0.029  # omega^2_Vd/F  = 0.029 (0.018 to 0.048)
    etalvmax ~ 0.027  # omega^2_Vmax  = 0.027 (0.017 to 0.042)
    etalkm   ~ 1.13   # omega^2_Km    = 1.13  (0.69 to 1.83)

    # Residual error - proportional; sigma^2 = 0.028 -> SD = sqrt(0.028) ~ 0.167.
    propSd <- sqrt(0.028); label("Proportional residual SD on Cc (fraction)") # Table II final-model: sigma^2 = 0.028 (0.020 to 0.038)
  })

  model({
    # Individual PK parameters - additive sex / ALDH2 / ADH1B effects on the
    # linear-scale typical value, power-form WT / AGE effects, multiplicative
    # log-normal IIV applied at the end. Mirrors Nemoto 2017 Table II structural
    # equations: ka_i = (ka + theta_SEX(ka)*FEMALE) * (age/29.4)^theta_AGE(ka) * exp(eta_ka),
    # Vd/F_i = (Vd/F + theta_SEX(Vd/F)*FEMALE + theta_ALDH2(Vd/F)*ALDH2) *
    #          (WT/61.3)^theta_WT(Vd/F) * (age/29.4)^theta_AGE(Vd/F) * exp(eta_Vd/F),
    # Vmax_i = Vmax_ADH1B[genotype] * (WT/61.3)^theta_WT(Vmax) * exp(eta_Vmax),
    # Km_i = Km * exp(eta_Km).
    ka   <- (exp(lka) + e_sexf_ka * SEXF) * (AGE / 29.4)^e_age_ka * exp(etalka)
    vc   <- (exp(lvc) + e_sexf_vc * SEXF + e_aldh2_s2_carrier_vc * ALDH2_S2_CARRIER) *
            (WT / 61.3)^e_wt_vc * (AGE / 29.4)^e_age_vc * exp(etalvc)
    vmax <- (exp(lvmax) + e_adh1b_s2_hom_vmax * ADH1B_S2_HOM) *
            (WT / 61.3)^e_wt_vmax * exp(etalvmax)
    km   <- exp(lkm + etalkm)

    # 1-compartment first-order absorption with Michaelis-Menten elimination.
    # Mass-unit consistency: dose mg, vc L, Cc mg/L, vmax mg/h, km mg/L.
    Cc <- central / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - vmax * Cc / (km + Cc)

    Cc ~ prop(propSd)
  })
}
