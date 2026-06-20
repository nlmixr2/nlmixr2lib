Mensing_2017_ritonavir <- function() {
  description <- "One-compartment population PK model for oral ritonavir (co-dosed with paritaprevir as a CYP3A4 pharmacokinetic enhancer) in HCV genotype-1 infected adults receiving the 3D regimen (Mensing 2017). First-order absorption, linear elimination, combined proportional + additive residual error, IIV on CL/F only. The author's final model retained gender, creatinine clearance, and HCV genotype (1a vs 1b) as significant covariates on CL/F, but the paper does not publish point estimates for these covariate coefficients (only graphical exposure-ratio forest plots in Figure 2); the implemented model is the structural typical-value model with covariate coefficients omitted (documented in covariatesDataExcluded)."
  reference <- paste(
    "Mensing S, Eckert D, Sharma S, Polepally AR, Khatri A, Podsadecki TJ,",
    "Awni WM, Menon RM, Dutta S. (2017).",
    "Population pharmacokinetics of paritaprevir, ombitasvir, dasabuvir,",
    "ritonavir and ribavirin in hepatitis C virus genotype 1 infection:",
    "analysis of six phase III trials.",
    "Br J Clin Pharmacol 83(3):527-539.",
    "doi:10.1111/bcp.13138",
    sep = " "
  )
  vignette <- "Mensing_2017_3D_HCV_regimen"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list()

  covariatesDataExcluded <- list(
    SEXF = list(
      description        = "Sex (1 = female, 0 = male) indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Retained in the author's final ritonavir model as a significant covariate on CL/F (Table 3, Mensing 2017). Figure 2 reports a Cmax,ss ratio of 1.13 (0.95, 1.33) and AUC24,ss ratio of 1.15 (0.94, 1.39) for females vs males (<=15% higher exposures); paper does not publish covariate coefficient point estimates so the effect is not encoded in model().",
      source_name        = "SEX"
    ),
    CRCL = list(
      description        = "Baseline creatinine clearance (Cockcroft-Gault)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Retained in the author's final ritonavir model as a significant covariate on CL/F (Table 3, Mensing 2017). Figure 2 reports exposure ratios at CrCL 75 mL/min and 105 mL/min vs the DAA-pharmacokinetic-dataset median of 104 mL/min (Cmax,ss / AUC24,ss ratios of 1.13/1.12 at 75 mL/min and 1.00/1.00 at 105 mL/min); paper does not publish covariate coefficient point estimates so the effect is not encoded in model(). Cohort range 37.0-281.4 mL/min, median 104.0 mL/min.",
      source_name        = "CRCL"
    ),
    HCV_GT1B = list(
      description        = "HCV genotype-1 subtype indicator (1 = GT1B, 0 = GT1A)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (GT1a)",
      notes              = "Retained in the author's final ritonavir model as a significant covariate on CL/F (Table 3, Mensing 2017). Figure 2 reports a Cmax,ss ratio of 1.28 (1.28, 1.28) and AUC24,ss ratio of 1.33 (1.33, 1.33) for GT1a vs GT1b (28-33% higher exposures); per Mensing 2017 Discussion this is considered a chance finding 'based on random variation and not clinically meaningful' (no physiological basis for HCV subtype to influence ritonavir PK). Paper does not publish covariate coefficient point estimates so the effect is not encoded in model(). Note: the canonical HCV_GT1B = 1 for GT1B; the Figure 2 ratio is reported with GT1a in the numerator (i.e., GT1A vs GT1B), so the canonical encoding would flip the direction.",
      source_name        = "GT1A"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 2348L,
    n_studies      = 7L,
    age_range      = "18-71 years",
    age_median     = "54 years",
    weight_range   = "42-129 kg",
    weight_median  = "76 kg",
    sex_female_pct = 42,
    race_ethnicity = c(Asian = 2, Black = 7, Nonblack_NonAsian = 91),
    ethnicity_hispanic_latino_pct = 6,
    disease_state  = "Adults with chronic hepatitis C virus (HCV) genotype 1 infection (HCV RNA > 10,000 IU/mL). 16% had compensated cirrhosis (Child-Pugh A); none had moderate or severe hepatic impairment. 34% were peg-IFN/RBV treatment-experienced.",
    dose_range     = "Ritonavir 100 mg orally once daily coformulated with paritaprevir 150 mg and ombitasvir 25 mg, in combination with dasabuvir 250 mg twice daily and optional weight-based ribavirin; 12-week or 24-week treatment courses. Ritonavir functions here as a CYP3A4 pharmacokinetic enhancer for paritaprevir; the 100 mg ritonavir dose is sub-therapeutic for HIV protease-inhibition but sufficient to boost paritaprevir exposure ~30-fold.",
    regions        = "Multinational phase II (NCT01911845) and phase III (PEARL-II/III/IV, SAPPHIRE-I/II, TURQUOISE-II) studies.",
    notes          = "Demographic and clinical baseline characteristics from Mensing 2017 Table 2 (DAA pharmacokinetic data column, n = 2348). HCV genotype distribution: 53% GT1a, 47% GT1b."
  )

  ini({
    # Structural parameters -- Mensing 2017 Table 3 final-model estimates for
    # ritonavir (one-compartment model with first-order absorption and
    # elimination, combined proportional + additive residual error).
    lka <- log(2.32)  ; label("Absorption rate constant Ka (1/day)")                       # Mensing 2017 Table 3: ka = 2.32 (1.47, 2.77) day^-1
    lcl <- log(439)   ; label("Apparent clearance CL/F (L/day)")                            # Mensing 2017 Table 3: CL/F = 439 (369, 554) L/day = 18.3 L/h
    lvc <- log(21.5)  ; label("Apparent central volume of distribution Vc/F (L)")           # Mensing 2017 Table 3: Vc/F = 21.5 (6.85, 43.9) L

    # Inter-individual variability -- Mensing 2017 Table 3 reports IIV on CL/F
    # only (variance on the log scale = 0.810). No IIV reported on Vc/F.
    etalcl ~ 0.810                                                                          # Mensing 2017 Table 3: IIV CL/F = 0.810 (0.679, 1.11) variance on log scale -> ~107% CV

    # Residual unexplained variability -- Mensing 2017 Table 3 footnotes c
    # (proportional variance) and d (additive variance). Convert NONMEM
    # variance estimates to nlmixr2 SD parameters:
    #   propSd = sqrt(0.533) = 0.730   (fraction)
    #   addSd  = sqrt(4e-6)  = 0.002   (mg/L = ug/mL); = 2.0 ng/mL
    propSd <- 0.730   ; label("Proportional residual error (fraction)")                     # Mensing 2017 Table 3 footnote c: prop variance = 0.533 (0.465, 0.558); SD = sqrt(0.533) = 0.730
    addSd  <- 0.00200 ; label("Additive residual error (ug/mL)")                            # Mensing 2017 Table 3 footnote d: add variance = 4e-6 (2e-7, 8e-6); SD = sqrt(4e-6) = 0.00200 ug/mL = 2.00 ng/mL
  })

  model({
    # Individual PK parameters -- ka and Vc/F have no IIV per Mensing 2017
    # Table 3; only CL/F has between-subject variability.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)

    # Micro-constants
    kel <- cl / vc

    # One-compartment ODE with first-order absorption from a depot
    # compartment into the central compartment.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Plasma ritonavir concentration in ug/mL (= mg/L).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
