Mensing_2017_ombitasvir <- function() {
  description <- "One-compartment population PK model for oral ombitasvir in HCV genotype-1 infected adults receiving the 3D (paritaprevir/ritonavir + ombitasvir + dasabuvir) +/- ribavirin regimen (Mensing 2017). First-order absorption, linear elimination, combined proportional + additive residual error, IIV on CL/F only. The author's final model retained cirrhosis, gender, age, and body weight as significant covariates on CL/F (and age, body weight on Vc/F), but the paper does not publish point estimates for these covariate coefficients (only graphical exposure-ratio forest plots in Figure 2); the implemented model is the structural typical-value model with covariate coefficients omitted (documented in covariatesDataExcluded)."
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
      notes              = "Retained in the author's final ombitasvir model as a significant covariate on CL/F (Table 3, Mensing 2017). The Figure 2 forest plot reports a Cmax,ss ratio of 0.92 (0.88, 0.95) and AUC24,ss ratio of 0.90 (0.86, 0.95) for cirrhotic vs noncirrhotic patients; the paper does not publish the covariate coefficient point estimate so the effect is not encoded in model(). All study subjects had no cirrhosis or compensated cirrhosis (Child-Pugh A) only; moderate/severe hepatic impairment was excluded.",
      source_name        = "CIRR"
    ),
    SEXF = list(
      description        = "Sex (1 = female, 0 = male) indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Retained in the author's final ombitasvir model as a significant covariate on CL/F (Table 3, Mensing 2017). The Figure 2 forest plot reports a Cmax,ss ratio of 1.46 (1.41, 1.51) and AUC24,ss ratio of 1.54 (1.48, 1.60) for females vs males; the paper does not publish the covariate coefficient point estimate so the effect is not encoded in model(). Gender was the only covariate with a notable effect on ombitasvir exposures (46-54% higher in females).",
      source_name        = "SEX"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Retained in the author's final ombitasvir model as a significant covariate on BOTH CL/F and Vc/F (Table 3, Mensing 2017). Figure 2 reports exposure ratios at Age 44 years and Age 64 years vs the cohort median of 54 years (Cmax,ss / AUC24,ss ratios of 0.89/0.91 at age 44 and 1.10/1.08 at age 64); paper does not publish covariate coefficient point estimates so the effect is not encoded in model(). Cohort range 18-71 years, median 54 years.",
      source_name        = "AGE"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Retained in the author's final ombitasvir model as a significant covariate on BOTH CL/F and Vc/F (Table 3, Mensing 2017). Figure 2 reports exposure ratios at body weight 66 kg and 86 kg vs the cohort median of 76 kg (Cmax,ss / AUC24,ss ratios of 1.08/1.10 at 66 kg and 0.93/0.93 at 86 kg); paper does not publish covariate coefficient point estimates so the effect is not encoded in model(). Cohort range 42-129 kg, median 76 kg.",
      source_name        = "WT"
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
    dose_range     = "Ombitasvir 25 mg orally once daily coformulated with paritaprevir 150 mg and ritonavir 100 mg, in combination with dasabuvir 250 mg twice daily and optional weight-based ribavirin (1000 mg/day for body weight < 75 kg, 1200 mg/day for body weight >= 75 kg, divided into twice-daily doses); 12-week or 24-week treatment courses.",
    regions        = "Multinational phase II (NCT01911845) and phase III (PEARL-II/III/IV, SAPPHIRE-I/II, TURQUOISE-II) studies.",
    notes          = "Demographic and clinical baseline characteristics from Mensing 2017 Table 2 (DAA pharmacokinetic data column, n = 2348). HCV genotype distribution: 53% GT1a, 47% GT1b. Hispanic/Latino ethnicity 6%. Methadone/buprenorphine use 2%. Median creatinine clearance 104.0 mL/min (range 37.0-281.4) at baseline."
  )

  ini({
    # Structural parameters -- Mensing 2017 Table 3 final-model estimates for
    # ombitasvir (one-compartment model with first-order absorption and
    # elimination, combined proportional + additive residual error).
    # 95% confidence intervals reported in Table 3 are from a 500-replicate
    # bootstrap; reproduced here in comments.
    lka <- log(1.08)  ; label("Absorption rate constant Ka (1/day)")                       # Mensing 2017 Table 3: ka = 1.08 (1.01, 1.14) day^-1
    lcl <- log(453)   ; label("Apparent clearance CL/F (L/day)")                            # Mensing 2017 Table 3: CL/F = 453 (441, 467) L/day = 18.9 L/h
    lvc <- log(50.1)  ; label("Apparent central volume of distribution Vc/F (L)")           # Mensing 2017 Table 3: Vc/F = 50.1 (44.9, 55.8) L

    # Inter-individual variability -- Mensing 2017 Table 3 reports IIV on CL/F
    # only (variance on the log scale = 0.143 from the NONMEM $OMEGA block).
    # No IIV reported on Vc/F.
    etalcl ~ 0.143                                                                          # Mensing 2017 Table 3: IIV CL/F = 0.143 (0.131, 0.155) variance on log scale -> ~39% CV

    # Residual unexplained variability -- Mensing 2017 Table 3 footnotes c
    # (proportional, variance on linear scale) and d (additive, variance in
    # (concentration unit)^2). Mensing 2017's combined error model in linear
    # space:  C_obs = C_pred * (1 + eps_prop) + eps_add. Convert NONMEM
    # variance estimates to nlmixr2 SD parameters:
    #   propSd = sqrt(0.107)   = 0.327     (fraction)
    #   addSd  = sqrt(2.4e-5)  = 0.00490   (mg/L = ug/mL); = ~4.9 ng/mL
    propSd <- 0.327   ; label("Proportional residual error (fraction)")                     # Mensing 2017 Table 3 footnote c: prop variance = 0.107 (0.096, 0.122); SD = sqrt(0.107) = 0.327
    addSd  <- 0.00490 ; label("Additive residual error (ug/mL)")                            # Mensing 2017 Table 3 footnote d: add variance = 2.4e-5 (6e-6, 4.1e-5); SD = sqrt(2.4e-5) = 0.00490 ug/mL = 4.90 ng/mL
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

    # Plasma ombitasvir concentration. Dose in mg, vc in L gives Cc in
    # mg/L = ug/mL. The paper plots concentrations in ng/mL; multiply by
    # 1000 in downstream simulation code to match the published axis.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
