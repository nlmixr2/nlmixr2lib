Zhu_2021_ropeginterferon <- function() {
  description <- paste0(
    "One-compartment quasi-equilibrium target-mediated drug disposition ",
    "(QE-TMDD) population pharmacokinetic model for subcutaneous ",
    "ropeginterferon alfa-2b (a mono-PEGylated interferon alfa-2b) in 57 ",
    "healthy adult volunteers pooled from two single-dose phase I studies: ",
    "30 Caucasian men (A09-102, 24-270 ug) and 27 Chinese men and women ",
    "(A17-101, 90-270 ug). First-order absorption with a lag time into a ",
    "single serum compartment carrying linear clearance plus saturable ",
    "binding (KD) to a turnover receptor pool (R0, kdeg), with the ",
    "drug-receptor complex internalised at kint. Body weight acts on the ",
    "linear clearance as a power function referenced to 70 kg; ethnicity ",
    "had no significant effect."
  )
  reference <- paste(
    "Zhu M, Wang M-X, Li Z-R, Wang W, Su X, Jiao Z.",
    "Population Pharmacokinetics of Ropeginterferon Alfa-2b: A Comparison",
    "Between Healthy Caucasian and Chinese Subjects.",
    "Front Pharmacol. 2021;12:673492.",
    "doi:10.3389/fphar.2021.673492.",
    sep = " "
  )
  vignette <- "Zhu_2021_ropeginterferon"
  units <- list(
    time = "hour (Zhu 2021 Table 2 reports CL/F in L/day and ka in 1/day but tlag, kint and kdeg in hours; the day-unit values are divided by 24 inside ini() so that every rate in the model is per hour -- see the vignette Assumptions for why the mixed printed units are genuine)",
    dosing = "ug (micrograms of ropeginterferon alfa-2b, subcutaneous)",
    concentration = "ug/L total serum ropeginterferon alfa-2b (Cc), numerically identical to the ng/mL Zhu 2021 reports; doses in ug over a volume in L give ug/L directly"
  )

  compartmentData <- list(
    depot = list(analyte = "ropeginterferon alfa-2b", units = "ug", specimen = "administration site", verified = FALSE),
    central = list(
      analyte = "ropeginterferon alfa-2b (total, free plus receptor-bound)",
      units = "ug",
      specimen = "serum",
      verified = FALSE
    ),
    total_target = list(
      analyte = "interferon alfa receptor (total, free plus drug-bound)",
      units = "ug/L (= ng/mL, drug-equivalent binding capacity)",
      specimen = "serum",
      verified = FALSE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The only covariate retained in the final model (Zhu 2021 Results,",
        "Covariate Model; Supplementary Table 3 forward step dOFV -8.309,",
        "backward step +8.309). Enters linear clearance as Equation 15,",
        "CL/F = 0.778 * (WT/70)^0.927. The abstract and Discussion both",
        "quote the typical CL/F 'in 70-kg subjects', so the reference",
        "weight is 70 kg as printed in Equation 15. Pooled mean 74 +/- 11.2",
        "kg (Table 1): Caucasian 79.4 +/- 9.41, Chinese 68.1 +/- 10."
      ),
      source_name = "WT"
    )
  )

  covariatesDataExcluded <- list(
    RACE_CHINESE = list(
      description = "Chinese-heritage indicator; 1 = Chinese (study A17-101), 0 = Caucasian (study A09-102).",
      units = "(binary)",
      type = "binary",
      notes = "Screened on CL/F and ka (Supplementary Table 3). Significant on ka in forward step 1 (dOFV -4.326) but not after weight entered CL/F (step 2 dOFV -1.534, p = 0.216) and not retained. The paper's headline conclusion is that there is no ethnic difference after adjusting for body weight."
    ),
    SEXF = list(
      description = "Female-sex indicator; 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      notes = "Screened on CL/F and ka as the Equation 14 proportional categorical effect (Supplementary Table 3); not retained. All 30 Caucasian subjects were male; 12 of the 27 Chinese subjects were female (Table 1)."
    ),
    AGE = list(
      description = "Age at baseline.",
      units = "years",
      type = "continuous",
      notes = "Screened by visual inspection only; not carried into the stepwise search. Pooled mean 32.3 +/- 6.52 years (Table 1)."
    ),
    BMI = list(
      description = "Baseline body mass index.",
      units = "kg/m^2",
      type = "continuous",
      notes = "Correlated with CL/F on visual inspection, but body weight was chosen among the collinear weight / BMI / BSA trio (Results, Covariate Model). Pooled mean 25.2 +/- 2.68 kg/m^2 (Table 1)."
    ),
    BSA = list(
      description = "Baseline body surface area.",
      units = "m^2",
      type = "continuous",
      notes = "Correlated with CL/F on visual inspection, but body weight was chosen among the collinear weight / BMI / BSA trio. Pooled mean 1.89 +/- 0.182 m^2 (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 57L,
    n_studies = 2L,
    n_observations = 894L,
    age_range = "18-45 years by protocol; mean 32.3 +/- 6.52 (Zhu 2021 Table 1)",
    weight_range = "mean 74 +/- 11.2 kg pooled; Caucasian 79.4 +/- 9.41, Chinese 68.1 +/- 10 (Zhu 2021 Table 1)",
    sex_female_pct = 21.05,
    race_ethnicity = "Caucasian 30 (52.6%, all male, A09-102, Canada); Chinese 27 (47.4%, 15 male / 12 female, A17-101, Beijing)",
    disease_state = "Healthy volunteers",
    dose_range = "Single subcutaneous dose: 24, 48, 90, 180, 225 or 270 ug (A09-102, six per cohort before drop-out) and 90, 180 or 270 ug (A17-101, ten per cohort before drop-out)",
    regions = "Canada (Montreal) and China (Beijing)",
    notes = paste(
      "66 subjects received ropeginterferon alfa-2b; 9 (6 Caucasian, 3",
      "Chinese) were excluded for drop-out, leaving 57 subjects with 894",
      "serum concentrations (456 Caucasian, 438 Chinese). Sampling pre-dose",
      "to 672 h post-dose. Double-antibody sandwich assay, LLOQ 50 pg/mL,",
      "cross-validated between the two laboratories. NONMEM 7.3, FOCE-I."
    )
  )

  ini({
    # ==================================================================
    # Zhu 2021 Table 2, "Parameter estimates of the final model and
    # bootstrap evaluation". NONMEM 7.3, FOCE with eta-epsilon
    # interaction.
    #
    # UNITS. Table 2 prints CL/F in L/day and ka in 1/day, but tlag in
    # h and kint / kdeg in 1/h; Supplementary Table 2 repeats the same
    # mix (kint 'h^-1', CL/F 'L/day', ka 'day^-1'). The model runs in
    # HOURS, so the two per-day values are divided by 24 here. The mix
    # is genuine, not a typesetting slip: with kint / kdeg per hour the
    # typical-value AUC0-inf reproduces the Supplementary Table 1 NCA
    # at every dose from 24 to 270 ug, including its greater-than-
    # proportional rise; read per day, the target arm becomes
    # negligible and the 24 ug AUC is overpredicted about 2-fold (see
    # the vignette Assumptions).
    # ==================================================================

    # ----- Absorption -----
    lka   <- log(0.14 / 24)  ; label("First-order absorption rate constant (1/h)")              # Zhu 2021 Table 2: ka 0.14 1/day (RSE 14%); /24 to 1/h
    ltlag <- log(0.426)      ; label("Absorption lag time (h)")                                 # Zhu 2021 Table 2: tlag 0.426 h (RSE 9%)

    # ----- Linear disposition -----
    lcl   <- log(0.778 / 24) ; label("Linear apparent clearance CL/F for a 70 kg subject (L/h)") # Zhu 2021 Table 2: CL/F 0.778 L/day (RSE 12%); /24 to L/h
    lvc   <- log(2.32)       ; label("Apparent volume of distribution V/F (L)")                 # Zhu 2021 Table 2: V/F 2.32 L (RSE 14%)

    # ----- Target-mediated disposition (quasi-equilibrium) -----
    lrbase <- log(0.111)       ; label("Baseline total receptor concentration R0 (ng/mL)")                  # Zhu 2021 Table 2: R0 0.111 ng/mL (RSE 31%)
    lkint  <- fixed(log(0.0788)); label("Internalisation rate constant of the drug-receptor complex (1/h)") # Zhu 2021 Table 2: kint 0.0788 1/h (Fixed); Results: 'kint was fixed at 0.0788 h-1'
    lkdeg  <- log(0.544)       ; label("Free receptor degradation rate constant (1/h)")                     # Zhu 2021 Table 2: kdeg 0.544 1/h (RSE 44%)
    lkd    <- fixed(log(0.142)) ; label("Equilibrium dissociation constant KD (ng/mL)")                     # Zhu 2021 Table 2: KD 0.142 ng/mL (Fixed); Results: 'The KD was fixed at 0.142 ng/ml'

    # ----- Covariate effect -----
    e_wt_cl <- 0.927 ; label("Power exponent of (WT/70) on CL/F (unitless)")  # Zhu 2021 Table 2: 'Impact of body weight' 0.927 (RSE 43%); Equation 15

    # ----- Between-subject variability -----
    # Table 2 reports BSV as CV% for exponential etas (Equation 9);
    # omega^2 = log(CV^2 + 1). No covariances are reported.
    etalka ~ 0.334133  # Zhu 2021 Table 2: BSV ka 63.5% CV (shrinkage 12%); log(0.635^2 + 1)
    etalcl ~ 0.120128  # Zhu 2021 Table 2: BSV CL/F 35.7% CV (shrinkage 15%); log(0.357^2 + 1)
    etalvc ~ 0.601330  # Zhu 2021 Table 2: BSV V/F 90.8% CV (shrinkage 5%); log(0.908^2 + 1)

    # ----- Residual unexplained variability -----
    # Equation 12, Y = IPRED * (1 + eps_prop) + eps_add, with separate
    # (diagonal) epsilons -- the nlmixr2 combined2 default.
    propSd <- 0.187 ; label("Proportional residual SD (fraction)")  # Zhu 2021 Table 2: 'Proptional (%)' 18.7 (RSE 2%)
    addSd  <- 0.342 ; label("Additive residual SD (ng/mL)")         # Zhu 2021 Table 2: 'Additive (ng/ml)' 0.342 (RSE 3%)
  })

  model({
    # ----- Individual parameters (Equation 9, exponential BSV) -----
    ka    <- exp(lka + etalka)
    tlag  <- exp(ltlag)
    # Equation 15: CL/F = 0.778 * (WT/70)^0.927
    cl    <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc    <- exp(lvc + etalvc)
    rbase <- exp(lrbase)
    kint  <- exp(lkint)
    kdeg  <- exp(lkdeg)
    kd    <- exp(lkd)
    # Equation 5: Rtotal(0) = R0 = ksyn / kdeg
    ksyn  <- rbase * kdeg

    # ----- Quasi-equilibrium free drug (Equation 6) -----
    ctot  <- central / vc
    disc  <- ctot - total_target - kd
    cfree <- 0.5 * (disc + sqrt(disc * disc + 4 * kd * ctot))

    total_target(0) <- rbase

    # ----- ODEs (Equations 2-4) -----
    # central holds the total drug amount Atotal; Afree = cfree * vc and
    # Atotal - Afree = (ctot - cfree) * vc is the bound amount.
    d/dt(depot)        <- -ka * depot
    d/dt(central)      <-  ka * depot - cl * cfree - kint * (ctot - cfree) * vc
    d/dt(total_target) <-  ksyn - kdeg * total_target - (kint - kdeg) * (ctot - cfree)
    alag(depot)        <- tlag

    # ----- Observation -----
    # The sandwich immunoassay measures total serum drug; the model is
    # fit to Ctotal.
    Cc <- ctot
    Cc ~ add(addSd) + prop(propSd)
  })
}
