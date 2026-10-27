Dari_2023b_ad26cov2s <- function() {
  description <- "QSP. Seven-compartment mechanistic humoral-immune-response model of spike-protein-binding antibody (S-ELISA, EU/mL) and wild-type-virus neutralizing antibody (wtVNA, IC50 titer) after a single dose or a homologous booster (2, 3 or 6 months after the prime) of the COVID-19 vaccine Ad26.COV2.S (5 x 10^10 viral particles per dose) in SARS-CoV-2-seronegative adults (Dari 2023, CPT Pharmacometrics Syst Pharmacol). A K-PD virtual antigen compartment drives two B-cell pathways: logistic memory-B-cell expansion that differentiates into short-lived plasma cells, and germinal-centre formation that differentiates into long-lived plasma cells. Plasma cells secrete antibody distributing between serum and a peripheral compartment; the neutralizing titer is linked to the binding antibody through a decaying IgM contribution and a log-linear slope. Age lowers both antibody production rates (power model, reference 55 years); female sex raises the long-lived-plasma-cell production rate. Time zero must be the prime dose. The antigen dose is entered in units of 10^8 viral particles (amt = 500 per 5 x 10^10 vp dose); the paper does not state the dosing scale, and the maintainers selected this reading as the only power-of-ten scaling of the nominal dose that reproduces the Figure 3 VPC and the Table 2 predicted persistence."
  reference <- paste(
    "Dari A, Boulton M, Neyens M, Le Gars M, Valenzuela B, Shukarev G,",
    "Cardenas V, Ruiz-Guinazu J, Sadoff J, Hoetelmans RMW, Perez Ruixo JJ.",
    "Mechanistic modeling projections of antibody persistence after homologous",
    "booster regimens of COVID-19 vaccine Ad26.COV2.S in humans.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(10):1485-1498.",
    "doi:10.1002/psp4.13025.",
    sep = " "
  )
  vignette <- "Dari_2023b_ad26cov2s"
  units <- list(
    time = "day",
    dosing = "1e8 viral particles",
    concentration = "EU/mL"
  )

  # Mechanistic immune-response states of the Dari 2023 Figure 1 schematic;
  # none has a canonical compartment name.
  paper_specific_compartments <- c(
    "antigen",
    "bcell_memory",
    "germinal_center",
    "plasma_short",
    "plasma_long",
    "antibody",
    "antibody_periph"
  )

  compartmentData <- list(
    antigen = list(
      analyte = "SARS-CoV-2 spike protein (virtual K-PD antigen)",
      units = "1e8 viral particles",
      specimen = "not applicable",
      verified = TRUE
    ),
    bcell_memory = list(
      analyte = "spike-specific memory B cells (M1)",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    germinal_center = list(
      analyte = "spike-specific germinal-centre B cells (M2)",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    plasma_short = list(
      analyte = "spike-specific short-lived plasma cells",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    plasma_long = list(
      analyte = "spike-specific long-lived plasma cells",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    antibody = list(
      analyte = "anti-spike binding IgG (S-ELISA)",
      units = "EU/mL",
      specimen = "serum",
      verified = TRUE
    ),
    antibody_periph = list(
      analyte = "anti-spike binding IgG (S-ELISA)",
      units = "EU/mL",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    AGE = list(
      description = "Age at vaccination",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects on both antibody production rates, reference 55 years (the median age of the S-ELISA dataset): kAbS * (AGE/55)^e_age_kab_s and kAbL * (AGE/55)^e_age_kab_l (Table 1; Appendix S1 'Covariate analysis').",
      source_name = "Age"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Multiplicative fractional effect on the long-lived-plasma-cell antibody production rate kAbL: kAbL * (1 + e_sexf_kab_l * SEXF) (Table 1 'kAbL_SEX'; Appendix S1 'Covariate analysis').",
      source_name = "Sex"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on kAbS, kAbL and the intercept but not retained in the final model (Methods 'Covariate analysis'; Appendix S1). Race and country were screened the same way and not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 978L,
    n_studies = 5L,
    n_observations = 7490L,
    age_range = "18-84 years",
    age_median = "55 years (S-ELISA); 47 years (wtVNA subset)",
    weight_range = "40.8-160 kg",
    weight_median = "73.8 kg",
    sex_female_pct = 48.6,
    race_ethnicity = c(
      `White non-Hispanic or Latino` = 72.8,
      Asian = 12.3,
      Other = 14.9
    ),
    disease_state = "Healthy SARS-CoV-2-seronegative adults; participants who became infected were withdrawn from that point",
    dose_range = "Ad26.COV2.S 5 x 10^10 viral particles, single dose (day 1) or prime plus homologous booster at 2 months (days 25-70), 3 months (days 70-130) or 6 months (day > 130)",
    regions = "Multinational (COV1001, COV1002, COV2001, COV3001, COV3009)",
    notes = "Pooled phase I-III trials COV1001 (NCT04436276), COV1002 (NCT04509947), COV2001 (NCT04535453), COV3001 (NCT04505722) and COV3009 (NCT04614948), data up to 15 months. 5608 S-ELISA samples from 978 participants and 1882 wtVNA samples from a subset of 375 participants (median age 47, range 18-80, 45.9% women). Demographics from Table S1. LLOQ 50.3 EU/mL (S-ELISA) and 58 IC50 (wtVNA)."
  )

  ini({
    # Antigen (virtual K-PD compartment)
    k_ag <- fixed(0.065)
    label("Apparent decay rate of the vaccine-generated S-protein antigen (1/day)") # Table 1 'Ad26 antigen decay, kAg' = 0.065, fixed

    # Memory-B-cell (M1) and germinal-centre (M2) pathways
    kp_ag <- fixed(1)
    label("Antigen-dependent production rate of memory B cells and germinal centres (cells per antigen unit per day)") # Table 1 'kpAg' = 1, fixed
    lkp_m <- log(0.0167)
    label("Antigen-independent memory-B-cell growth rate kpM (1/day)") # Table 1 'kpM' = 0.0167 (95% CI 0.0121-0.0213)
    ke_m <- fixed(0)
    label("Memory-B-cell elimination rate keM1 (1/day)") # Table 1 'keM1 and keM2' = 0, fixed due to short follow-up
    ke_gc <- fixed(0)
    label("Germinal-centre elimination rate keM2 (1/day)") # Table 1 'keM1 and keM2' = 0, fixed due to short follow-up
    n_mem <- fixed(10.7e4)
    label("Maximum memory-B-cell count N1 (cells)") # Table 1 'N1' = 10.7 x 10^4, fixed to the mean of multiple runs
    ln_gc <- log(0.129e4)
    label("Maximum germinal-centre B-cell count N2 (cells)") # Table 1 'N2' = 0.129 x 10^4 (95% CI 0.033-0.225 x 10^4)

    # Plasma-cell differentiation and loss
    kp_s <- fixed(0.5e-5)
    label("Short-lived plasma-cell differentiation rate kpS (1/day per antigen unit)") # Table 1 'kpS' = 0.5 x 10^-5, fixed (sensitivity analysis)
    kp_l <- fixed(0.5e-5)
    label("Long-lived plasma-cell differentiation rate kpL (1/day per antigen unit)") # Table 1 'kpL' = 0.5 x 10^-5, fixed (sensitivity analysis)
    ke_s <- fixed(0.578)
    label("Short-lived plasma-cell elimination rate keS (1/day)") # Table 1 'keS' = 0.578 /day, fixed
    lke_l <- log(19.1e-4)
    label("Long-lived plasma-cell elimination rate keL (1/day)") # Table 1 'keL' = 19.1 x 10^-4 (95% CI 15.2-23.0 x 10^-4); unit printed as /year, read as /day (see vignette)

    # Antibody production and disposition
    lkab_s <- log(3.01)
    label("Antibody production rate per short-lived plasma cell kAbS at age 55 years ((EU/mL)/cell/day)") # Table 1 'kAbS' = 3.01 (95% CI 2.60-3.42)
    lkab_l <- log(0.394)
    label("Antibody production rate per long-lived plasma cell kAbL, male aged 55 years ((EU/mL)/cell/day)") # Table 1 'kAbL' = 0.394 (95% CI 0.143-0.645)
    kab_p <- fixed(0.2075)
    label("Serum-to-peripheral antibody distribution rate kAbP (1/day)") # Table 1 'kAbP' = 0.2075, fixed
    kp_ab <- fixed(0.2716)
    label("Peripheral-to-serum antibody return rate kPAb (1/day)") # Table 1 'kPAb' = 0.2716, fixed
    ke_ab <- fixed(0.0556)
    label("Serum antibody elimination rate keAb (1/day)") # Table 1 'keAb' = 0.0556, fixed

    # Neutralizing-antibody link (Eq. 8-9)
    lalpha_wtvna <- log(0.732)
    label("Slope of log10 wtVNA titer on log10 neutralizing-equivalent binding antibody (unitless)") # Table 1 'alpha' = 0.732 (95% CI 0.691-0.773)
    lbeta_wtvna <- log(3.85)
    label("Constant ratio of neutralizing IgG to binding antibody (IC50 titer per EU/mL)") # Table 1 'beta' = 3.85 (95% CI 2.27-5.43)
    ligm_wtvna <- log(6.10)
    label("Baseline IgM contribution to the neutralizing-to-binding ratio (IC50 titer per EU/mL)") # Table 1 'IgM' = 6.10 (95% CI 2.98-9.22)
    lkigm <- log(0.0368)
    label("First-order decay rate of the IgM contribution kIgM (1/day)") # Table 1 'kIgM' = 0.0368 (95% CI 0.00858-0.065)

    # Covariate effects (Appendix S1: power model for continuous, theta * (1 + estimate * indicator) for categorical)
    e_age_kab_s <- -0.221
    label("Power exponent of (AGE/55) on kAbS (unitless)") # Table 1 'KAbS_AGE' = -0.221 (95% CI -0.421 to -0.0211)
    e_age_kab_l <- -0.867
    label("Power exponent of (AGE/55) on kAbL (unitless)") # Table 1 'KAbL_AGE' = -0.867 (95% CI -1.06 to -0.677)
    e_sexf_kab_l <- 0.816
    label("Fractional change in kAbL for females (unitless)") # Table 1 'kAbL_SEX' = 0.816 (95% CI 0.536-1.10)

    # IIV; omega^2 = log(1 + CV^2) from the Table 1 CV% (footnote a)
    etalkab_s ~ 0.9605 # Table 1 IIV on kAbS CV 127%
    etalkab_l ~ 0.9507 # Table 1 IIV on kAbL CV 126%

    # Residual error: additive on log10 titers (Appendix S1 'Additive error
    # model'; Table 1 footnote c). The Table 1 'Proportional error' entries are
    # read as log10-scale SDs x 100, as for the single-dose model of the same
    # authors, so the natural-log SD is 0.21 * log(10) = 0.4835 (S-ELISA) and
    # 0.238 * log(10) = 0.5480 (wtVNA). See the vignette.
    expSd <- 0.4835
    label("Log-scale residual SD for S-ELISA titers (natural-log units)") # Table 1 'Proportional error' S-ELISA = 21 (log10 SD 0.21)
    expSd_wtvna <- 0.5480
    label("Log-scale residual SD for wtVNA titers (natural-log units)") # Table 1 'Proportional error' wtVNA = 23.8 (log10 SD 0.238)
  })

  model({
    # Individual parameters
    kp_m <- exp(lkp_m)
    n_gc <- exp(ln_gc)
    ke_l <- exp(lke_l)
    kab_s <- exp(lkab_s + etalkab_s) * (AGE / 55)^e_age_kab_s
    kab_l <- exp(lkab_l + etalkab_l) * (AGE / 55)^e_age_kab_l * (1 + e_sexf_kab_l * SEXF)
    alpha_wtvna <- exp(lalpha_wtvna)
    beta_wtvna <- exp(lbeta_wtvna)
    igm_wtvna <- exp(ligm_wtvna)
    kigm <- exp(lkigm)

    # Eq. 1-7. Seronegative participants: every state starts at zero and each
    # vaccine dose is a bolus into antigen.
    d/dt(antigen) <- -k_ag * antigen
    d/dt(bcell_memory) <- (kp_ag * antigen + kp_m * bcell_memory) * (1 - bcell_memory / n_mem) -
      (kp_s * antigen + ke_m) * bcell_memory
    d/dt(germinal_center) <- kp_ag * antigen * (1 - germinal_center / n_gc) -
      (kp_l * antigen + ke_gc) * germinal_center
    d/dt(plasma_short) <- kp_s * antigen * bcell_memory - ke_s * plasma_short
    d/dt(plasma_long) <- kp_l * antigen * germinal_center - ke_l * plasma_long
    d/dt(antibody) <- kab_s * plasma_short + kab_l * plasma_long - ke_ab * antibody -
      kab_p * antibody + kp_ab * antibody_periph
    d/dt(antibody_periph) <- kab_p * antibody - kp_ab * antibody_periph

    # Observed serum binding-antibody concentration (S-ELISA, EU/mL). Named Cc,
    # as for the anti-RBD IgG output of Dasti_2025_mrna1273_qsp, although it is
    # an antibody response rather than a concentration of the administered
    # vaccine.
    Cc <- antibody

    # Eq. 8 scales the binding antibody by the neutralizing-to-binding ratio,
    # whose IgM part decays with time since the prime; Eq. 9 is the log10 slope
    # applied to that neutralizing-equivalent level (IC50 titer), i.e.
    # log10(wtvna) = alpha * log10(ratio * antibody).
    wtvna <- ((beta_wtvna + igm_wtvna * exp(-kigm * t)) * antibody)^alpha_wtvna

    Cc ~ lnorm(expSd)
    wtvna ~ lnorm(expSd_wtvna)
  })
}
