Mensing_2017_ribavirin <- function() {
  description <- "Two-compartment population PK model for oral ribavirin in HCV genotype-1 infected adults receiving the 3D + ribavirin regimen (Mensing 2017). First-order absorption, linear elimination, combined proportional + additive residual error, IIV on CL/F and a shared IIV on Vc/F + Vp/F. The author's final model retained cirrhosis, gender, and creatinine clearance as significant covariates on CL/F (and gender on Vc/F and Vp/F), but the paper does not publish point estimates for these covariate coefficients (only graphical exposure-ratio forest plots in Figure 2); the implemented model is the structural typical-value model with covariate coefficients omitted (documented in covariatesDataExcluded). Mensing 2017 reports correlated IIV on CL/F and Vc/Vp; the correlation coefficient is not given in Table 3, so this implementation encodes the random effects as independent (a documented deviation; see vignette Errata)."
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
      notes              = "Retained in the author's final ribavirin model as a significant covariate on CL/F (Table 3, Mensing 2017). Figure 2 reports a Cmax,ss ratio of 0.95 (0.93, 0.96) and AUC24,ss ratio of 0.94 (0.92, 0.96) for cirrhotic vs noncirrhotic patients; paper does not publish covariate coefficient point estimates so the effect is not encoded in model(). Mensing 2017 Results: ribavirin exposures were 'similar in subjects with or without compensated cirrhosis'.",
      source_name        = "CIRR"
    ),
    SEXF = list(
      description        = "Sex (1 = female, 0 = male) indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Retained in the author's final ribavirin model as a significant covariate on CL/F AND on Vc/F and Vp/F (Table 3, Mensing 2017). Figure 2 reports a Cmax,ss ratio of 1.31 (1.26, 1.35) and AUC24,ss ratio of 1.29 (1.24, 1.33) for females vs males (29-31% higher exposures); paper does not publish covariate coefficient point estimates so the effect is not encoded in model().",
      source_name        = "SEX"
    ),
    CRCL = list(
      description        = "Baseline creatinine clearance (Cockcroft-Gault)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Retained in the author's final ribavirin model as a significant covariate on CL/F (Table 3, Mensing 2017). Figure 2 reports exposure ratios at CrCL 75 mL/min and 105 mL/min vs the ribavirin-pharmacokinetic-dataset median of 105 mL/min (Cmax,ss / AUC24,ss ratios of 1.08/1.08 at 75 mL/min and 1.00/1.00 at 105 mL/min); paper does not publish covariate coefficient point estimates so the effect is not encoded in model(). Cohort range 37.0-241.0 mL/min, median 111.2 mL/min for the ribavirin pharmacokinetic dataset (Table 2).",
      source_name        = "CRCL"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1841L,
    n_studies      = 6L,
    age_range      = "18-71 years",
    age_median     = "54 years",
    weight_range   = "42-129 kg",
    weight_median  = "77 kg",
    sex_female_pct = 41,
    race_ethnicity = c(Asian = 2, Black = 6, Nonblack_NonAsian = 92),
    ethnicity_hispanic_latino_pct = 6,
    disease_state  = "Adults with chronic hepatitis C virus (HCV) genotype 1 infection (HCV RNA > 10,000 IU/mL) receiving 3D + ribavirin therapy. 21% had compensated cirrhosis (Child-Pugh A); none had moderate or severe hepatic impairment. 38% were peg-IFN/RBV treatment-experienced.",
    dose_range     = "Ribavirin orally twice daily, weight-based dosing (1000 mg/day total for body weight < 75 kg, 1200 mg/day total for body weight >= 75 kg), in combination with paritaprevir/ritonavir/ombitasvir 150/100/25 mg once daily and dasabuvir 250 mg twice daily; 12-week or 24-week treatment courses.",
    regions        = "Multinational phase II (NCT01911845) and phase III (PEARL-II/III/IV, SAPPHIRE-I/II, TURQUOISE-II) studies; excludes the no-ribavirin arms.",
    notes          = "Demographic and clinical baseline characteristics from Mensing 2017 Table 2 (ribavirin pharmacokinetic data column, n = 1841). HCV genotype distribution: 57% GT1a, 43% GT1b. Median creatinine clearance 111.2 mL/min (range 37.0-241.0)."
  )

  ini({
    # Structural parameters -- Mensing 2017 Table 3 final-model estimates for
    # ribavirin (two-compartment model with first-order absorption and
    # elimination, combined proportional + additive residual error, IIV on
    # CL/F and a shared IIV on Vc/F + Vp/F).
    lka <- log(21.3) ; label("Absorption rate constant Ka (1/day)")                        # Mensing 2017 Table 3: ka = 21.3 (18.7, 24.1) day^-1
    lcl <- log(427)  ; label("Apparent clearance CL/F (L/day)")                             # Mensing 2017 Table 3: CL/F = 427 (419, 436) L/day = 17.8 L/h
    lvc <- log(1100) ; label("Apparent central volume of distribution Vc/F (L)")            # Mensing 2017 Table 3: Vc/F = 1100 (983, 1230) L
    lq  <- log(877)  ; label("Apparent inter-compartmental clearance Q/F (L/day)")          # Mensing 2017 Table 3: Q/F = 877 (791, 977) L/day
    lvp <- log(3230) ; label("Apparent peripheral volume of distribution Vp/F (L)")         # Mensing 2017 Table 3: Vp/F = 3230 (3070, 3380) L

    # Inter-individual variability -- Mensing 2017 Table 3 reports IIV on CL/F
    # (variance 0.062) and a shared IIV on Vc/F + Vp/F (variance 0.197).
    # Mensing 2017 Results state "correlated IIV on CL/F and Vc/Vp"; the
    # correlation value is not given in Table 3 or the supplement (S1-S5
    # contain only P-values for covariates), so the random effects are
    # encoded as independent here. This is a documented deviation from the
    # published model; see vignette Errata.
    etalcl ~ 0.062                                                                          # Mensing 2017 Table 3: IIV CL/F = 0.062 (0.057, 0.067) variance on log scale -> ~25% CV
    etalvc ~ 0.197                                                                          # Mensing 2017 Table 3: IIV Vc/F + Vp/F (shared) = 0.197 (0.171, 0.222) variance on log scale -> ~47% CV. Applied to both Vc and Vp via the same eta in model().

    # Residual unexplained variability -- Mensing 2017 Table 3 footnotes c
    # (proportional variance) and d (additive variance). Convert NONMEM
    # variance estimates to nlmixr2 SD parameters:
    #   propSd = sqrt(0.0170) = 0.130   (fraction)
    #   addSd  = sqrt(0.0390) = 0.1975  (mg/L = ug/mL); = 197.5 ng/mL
    propSd <- 0.130   ; label("Proportional residual error (fraction)")                     # Mensing 2017 Table 3 footnote c: prop variance = 0.0170 (0.0140, 0.0190); SD = sqrt(0.0170) = 0.130
    addSd  <- 0.1975  ; label("Additive residual error (ug/mL)")                            # Mensing 2017 Table 3 footnote d: add variance = 0.0390 (0.0300, 0.0480); SD = sqrt(0.0390) = 0.1975 ug/mL = 197.5 ng/mL
  })

  model({
    # Individual PK parameters -- only CL/F and the shared Vc/F+Vp/F eta have
    # between-subject variability per Mensing 2017 Table 3. ka and Q/F have
    # no IIV. The shared volume eta (etalvc) is applied to BOTH vc and vp so
    # that a subject with a larger central volume also has a proportionally
    # larger peripheral volume (Mensing 2017's "IIV on Vc/Vp"); the
    # correlation between this volume eta and etalcl is not encoded here
    # (see Errata).
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq)
    vp <- exp(lvp + etalvc)

    # Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment ODE with first-order absorption from a depot
    # compartment into the central compartment.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # Plasma ribavirin concentration in ug/mL (= mg/L).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
