Mensing_2017_paritaprevir <- function() {
  description <- "One-compartment population PK model for oral paritaprevir (co-dosed with ritonavir) in HCV genotype-1 infected adults receiving the 3D regimen (Mensing 2017). First-order absorption with a fixed absorption lag time, linear elimination, additive residual error on log-transformed concentrations (encoded as log-normal Cc ~ lnorm(expSd)), IIV on CL/F only. The author's final model retained cirrhosis, gender, age, opioid use, and antidiabetic-agent use as significant covariates on CL/F (and age, body weight on Vc/F), but the paper does not publish point estimates for these covariate coefficients (only graphical exposure-ratio forest plots in Figure 2); the implemented model is the structural typical-value model with covariate coefficients omitted (documented in covariatesDataExcluded)."
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
      notes              = "Retained in the author's final paritaprevir model as a significant covariate on CL/F (Table 3, Mensing 2017). Figure 2 reports a Cmax,ss ratio of 2.22 (1.97, 2.54) and AUC24,ss ratio of 2.40 (2.11, 2.79) for cirrhotic vs noncirrhotic patients (122-140% higher exposures); the paper does not publish the covariate coefficient point estimate so the effect is not encoded in model().",
      source_name        = "CIRR"
    ),
    SEXF = list(
      description        = "Sex (1 = female, 0 = male) indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Retained in the author's final paritaprevir model as a significant covariate on CL/F (Table 3, Mensing 2017). Figure 2 reports a Cmax,ss ratio of 1.92 (1.74, 2.11) and AUC24,ss ratio of 1.96 (1.76, 2.15) for females vs males (92-96% higher exposures); the paper does not publish the covariate coefficient point estimate so the effect is not encoded in model().",
      source_name        = "SEX"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Retained in the author's final paritaprevir model as a significant covariate on BOTH CL/F and Vc/F (Table 3, Mensing 2017). Figure 2 reports exposure ratios at Age 44 years and Age 64 years vs the cohort median of 54 years (Cmax,ss / AUC24,ss ratios of 0.81/0.83 at age 44 and 1.18/1.17 at age 64); paper does not publish covariate coefficient point estimates so the effect is not encoded in model().",
      source_name        = "AGE"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Retained in the author's final paritaprevir model as a significant covariate on Vc/F (Table 3, Mensing 2017). Figure 2 reports near-unity exposure ratios at body weight 66 kg and 86 kg vs the cohort median of 76 kg (Cmax,ss / AUC24,ss ratios of 1.01/1.00 at 66 kg and 0.99/1.00 at 86 kg); paper does not publish covariate coefficient point estimates so the effect is not encoded in model().",
      source_name        = "WT"
    ),
    CONMED_OPIOID = list(
      description        = "Concomitant opioid use indicator (methadone, buprenorphine, with or without naloxone)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant opioid)",
      notes              = "Retained in the author's final paritaprevir model as a significant covariate on CL/F (Table 3, Mensing 2017). Figure 2 reports a Cmax,ss ratio of 1.48 (1.27, 1.74) and AUC24,ss ratio of 1.56 (1.34, 1.86) for opioid users vs non-users (48-56% higher exposures); the paper does not publish the covariate coefficient point estimate so the effect is not encoded in model(). The phase II study (NCT01911845) enrolled subjects on stable opioid replacement therapy; methadone/buprenorphine use was 2% in the overall DAA pharmacokinetic data.",
      source_name        = "OPIOID"
    ),
    CONMED_ANTIDIAB = list(
      description        = "Concomitant antidiabetic-agent use indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant antidiabetic agent)",
      notes              = "Retained in the author's final paritaprevir model as a significant covariate on CL/F (Table 3, Mensing 2017). Figure 2 reports a Cmax,ss ratio of 1.40 (1.18, 1.75) and AUC24,ss ratio of 1.46 (1.20, 1.86) for antidiabetic users vs non-users (40-46% higher exposures); the paper does not publish the covariate coefficient point estimate so the effect is not encoded in model().",
      source_name        = "ANTIDIAB"
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
    dose_range     = "Paritaprevir 150 mg orally once daily coformulated with ritonavir 100 mg and ombitasvir 25 mg, in combination with dasabuvir 250 mg twice daily and optional weight-based ribavirin; 12-week or 24-week treatment courses. Ritonavir is co-dosed as a CYP3A4 pharmacokinetic enhancer for paritaprevir.",
    regions        = "Multinational phase II (NCT01911845) and phase III (PEARL-II/III/IV, SAPPHIRE-I/II, TURQUOISE-II) studies.",
    notes          = "Demographic and clinical baseline characteristics from Mensing 2017 Table 2 (DAA pharmacokinetic data column, n = 2348). The 2% of subjects on methadone/buprenorphine and the 4% on antidiabetic agents drove the opioid-use and antidiabetic-agent covariate effects on paritaprevir CL/F (Mensing 2017 Comedication analysis)."
  )

  ini({
    # Structural parameters -- Mensing 2017 Table 3 final-model estimates for
    # paritaprevir (one-compartment model with first-order absorption,
    # absorption lag time, and linear elimination). Ka and ALAG were fixed
    # at literature values during NONMEM estimation per Table 3 (see "(fixed)"
    # annotations in the ka and ALAG rows of the Paritaprevir column).
    lka   <- fixed(log(1.74))   ; label("Absorption rate constant Ka (1/day)")              # Mensing 2017 Table 3: ka = 1.74 (fixed) day^-1
    ltlag <- fixed(log(0.0400)) ; label("Absorption lag time tlag (day)")                   # Mensing 2017 Table 3: ALAG = 0.0400 (fixed) day = 0.96 h
    lcl   <- log(1580)          ; label("Apparent clearance CL/F (L/day)")                  # Mensing 2017 Table 3: CL/F = 1580 (1450, 1710) L/day = 65.9 L/h
    lvc   <- log(16.7)          ; label("Apparent central volume of distribution Vc/F (L)") # Mensing 2017 Table 3: Vc/F = 16.7 (11.8, 22.6) L

    # Inter-individual variability -- Mensing 2017 Table 3 reports IIV on CL/F
    # only (variance on the log scale = 1.18). No IIV reported on Vc/F or on
    # the (fixed) ka or ALAG terms. The very large IIV reflects the sparse
    # outpatient sampling and MEMS-cap reconstructed dosing histories.
    etalcl ~ 1.18                                                                            # Mensing 2017 Table 3: IIV CL/F = 1.18 (1.10, 1.26) variance on log scale -> ~150% CV (sqrt(exp(1.18) - 1))

    # Residual unexplained variability -- Mensing 2017 Methods state that
    # paritaprevir RUV was modelled using "an additive error model on log-
    # transformed data": log(C_obs) = log(C_pred) + eps, eps ~ N(0, sigma^2).
    # In nlmixr2 this maps to Cc ~ lnorm(expSd) with expSd as the SD on the
    # log scale. NONMEM reports the variance estimate (sigma^2) = 1.14 in
    # Table 3 footnote b ("Additive error on log-transformed data RUV"),
    # so expSd = sqrt(1.14) = 1.068.
    expSd <- 1.068    ; label("Log-scale (additive on log-transformed) residual error SD")  # Mensing 2017 Table 3 footnote b: variance = 1.14 (1.09, 1.20); SD = sqrt(1.14) = 1.068 on log scale
  })

  model({
    # Individual PK parameters -- ka, ALAG, and Vc/F have no IIV per Mensing
    # 2017 Table 3; only CL/F has between-subject variability.
    ka   <- exp(lka)
    tlag <- exp(ltlag)
    cl   <- exp(lcl + etalcl)
    vc   <- exp(lvc)

    # Micro-constants
    kel <- cl / vc

    # One-compartment ODE with first-order absorption from a depot
    # compartment into the central compartment, with absorption lag time
    # applied to the depot compartment (NONMEM ALAG1 equivalent).
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central
    alag(depot)   <- tlag

    # Plasma paritaprevir concentration in ug/mL (= mg/L). The paper plots
    # concentrations in ng/mL; multiply by 1000 in downstream simulation
    # code to match the published axis.
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
