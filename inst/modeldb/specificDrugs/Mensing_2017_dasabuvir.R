Mensing_2017_dasabuvir <- function() {
  description <- "Two-compartment population PK model for oral dasabuvir in HCV genotype-1 infected adults receiving the 3D regimen (Mensing 2017). First-order absorption, linear elimination, combined proportional + additive residual error, IIV on CL/F only. The author's final model retained cirrhosis, gender, creatinine clearance, and body weight as significant covariates on CL/F (and age, body weight on Vc/F and Vp/F), but the paper does not publish point estimates for these covariate coefficients (only graphical exposure-ratio forest plots in Figure 2); the implemented model is the structural typical-value model with covariate coefficients omitted (documented in covariatesDataExcluded)."
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
    HEPIMP_MILD = list(
      description        = "Compensated cirrhosis (Child-Pugh stage A) indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no cirrhosis)",
      notes              = "Retained in the author's final dasabuvir model as a significant covariate on CL/F (Table 3, Mensing 2017). Figure 2 reports a Cmax,ss ratio of 1.29 (1.22, 1.36) and AUC24,ss ratio of 1.39 (1.30, 1.49) for cirrhotic vs noncirrhotic patients (29-39% higher exposures); paper does not publish covariate coefficient point estimates so the effect is not encoded in model().",
      source_name        = "CIRR"
    ),
    SEXF = list(
      description        = "Sex (1 = female, 0 = male) indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Retained in the author's final dasabuvir model as a significant covariate on CL/F (Table 3, Mensing 2017). Figure 2 reports a Cmax,ss ratio of 1.16 (1.12, 1.21) and AUC24,ss ratio of 1.21 (1.15, 1.28) for females vs males (16-21% higher exposures); paper does not publish covariate coefficient point estimates so the effect is not encoded in model().",
      source_name        = "SEX"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Retained in the author's final dasabuvir model as a significant covariate on Vc/F and Vp/F (Table 3, Mensing 2017). Figure 2 reports exposure ratios at Age 44 years and Age 64 years vs the cohort median of 54 years (Cmax,ss / AUC24,ss ratios of 0.97/1.00 at age 44 and 1.02/1.00 at age 64); paper does not publish covariate coefficient point estimates so the effect is not encoded in model().",
      source_name        = "AGE"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Retained in the author's final dasabuvir model as a significant covariate on BOTH CL/F and Vc/F and Vp/F (Table 3, Mensing 2017). Figure 2 reports exposure ratios at body weight 66 kg and 86 kg vs the cohort median of 76 kg (Cmax,ss / AUC24,ss ratios of 1.05/1.04 at 66 kg and 0.96/0.97 at 86 kg); paper does not publish covariate coefficient point estimates so the effect is not encoded in model().",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Baseline creatinine clearance (Cockcroft-Gault)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Retained in the author's final dasabuvir model as a significant covariate on CL/F (Table 3, Mensing 2017). Figure 2 reports exposure ratios at CrCL 75 mL/min and 105 mL/min vs the DAA-pharmacokinetic-dataset median of 104 mL/min (Cmax,ss / AUC24,ss ratios of 1.06/1.10 at 75 mL/min and 1.00/1.00 at 105 mL/min); paper does not publish covariate coefficient point estimates so the effect is not encoded in model(). Cohort range 37.0-281.4 mL/min, median 104.0 mL/min.",
      source_name        = "CRCL"
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
    dose_range     = "Dasabuvir 250 mg orally twice daily, in combination with paritaprevir/ritonavir/ombitasvir (150/100/25 mg once daily) and optional weight-based ribavirin; 12-week or 24-week treatment courses.",
    regions        = "Multinational phase II (NCT01911845) and phase III (PEARL-II/III/IV, SAPPHIRE-I/II, TURQUOISE-II) studies.",
    notes          = "Demographic and clinical baseline characteristics from Mensing 2017 Table 2 (DAA pharmacokinetic data column, n = 2348)."
  )

  ini({
    # Structural parameters -- Mensing 2017 Table 3 final-model estimates for
    # dasabuvir (two-compartment model with first-order absorption and
    # elimination, combined proportional + additive residual error).
    lka <- log(4.61) ; label("Absorption rate constant Ka (1/day)")                        # Mensing 2017 Table 3: ka = 4.61 (3.99, 5.45) day^-1
    lcl <- log(1150) ; label("Apparent clearance CL/F (L/day)")                             # Mensing 2017 Table 3: CL/F = 1150 (1100, 1200) L/day = 47.9 L/h
    lvc <- log(110)  ; label("Apparent central volume of distribution Vc/F (L)")            # Mensing 2017 Table 3: Vc/F = 110 (93.3, 133) L
    lq  <- log(182)  ; label("Apparent inter-compartmental clearance Q/F (L/day)")          # Mensing 2017 Table 3: Q/F = 182 (111, 295) L/day
    lvp <- log(286)  ; label("Apparent peripheral volume of distribution Vp/F (L)")         # Mensing 2017 Table 3: Vp/F = 286 (190, 408) L

    # Inter-individual variability -- Mensing 2017 Table 3 reports IIV on CL/F
    # only (variance on the log scale = 0.263). No IIV reported on ka, Vc/F,
    # Vp/F, or Q/F.
    etalcl ~ 0.263                                                                          # Mensing 2017 Table 3: IIV CL/F = 0.263 (0.213, 0.299) variance on log scale -> ~56% CV

    # Residual unexplained variability -- Mensing 2017 Table 3 footnotes c
    # (proportional variance) and d (additive variance). Convert NONMEM
    # variance estimates to nlmixr2 SD parameters:
    #   propSd = sqrt(0.260)  = 0.510   (fraction)
    #   addSd  = sqrt(0.004)  = 0.0632  (mg/L = ug/mL); = 63.2 ng/mL
    propSd <- 0.510  ; label("Proportional residual error (fraction)")                      # Mensing 2017 Table 3 footnote c: prop variance = 0.260 (0.242, 0.283); SD = sqrt(0.260) = 0.510
    addSd  <- 0.0632 ; label("Additive residual error (ug/mL)")                             # Mensing 2017 Table 3 footnote d: add variance = 4.00e-3 (1.00e-3, 7.00e-3); SD = sqrt(4.00e-3) = 0.0632 ug/mL = 63.2 ng/mL
  })

  model({
    # Individual PK parameters -- only CL/F has between-subject variability
    # per Mensing 2017 Table 3.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    # Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment ODE with first-order absorption from a depot
    # compartment into the central compartment.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # Plasma dasabuvir concentration in ug/mL (= mg/L).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
