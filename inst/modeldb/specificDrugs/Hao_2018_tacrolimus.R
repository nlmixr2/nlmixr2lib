Hao_2018_tacrolimus <- function() {
  description <- "One-compartment population PK model with first-order absorption (no lag) and first-order elimination for twice-daily oral immediate-release tacrolimus (Prograf) in paediatric nephrotic-syndrome patients aged 2.7-17.3 years (Hao 2018). Apparent oral clearance CL/F scales allometrically with body weight at a fixed exponent of 0.75 referenced to a 70 kg adult; apparent volume of distribution V/F scales linearly with body weight at a fixed exponent of 1.0 referenced to 70 kg; ka has no body-weight scaling. CL/F additionally varies with CYP3A5 expresser status (multiplicative factor 1.60 for *1/*1 or *1/*3 carriers vs the *3/*3 nonexpresser reference). Inter-individual variability is diagonal on ka, V/F, and CL/F (exponential / log-normal model). Residual unexplained variability is proportional (paper text: 'The proportional model best described residual variability'; Table 2 reports it under the 'Residual variability (exponential)' label, which is the standard NONMEM additive-on-log-scale parameterisation equivalent to proportional in linear space)."
  reference <- "Hao GX, Huang X, Zhang DF, Zheng Y, Shi HY, Li Y, Jacqz-Aigrain E, Zhao W. Population pharmacokinetics of tacrolimus in children with nephrotic syndrome. Br J Clin Pharmacol. 2018;84(8):1748-1756. doi:10.1111/bcp.13605"
  vignette <- "Hao_2018_tacrolimus"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline per Hao 2018. Allometric power scaling on CL/F with reference 70 kg and the theory-based exponent 0.75 fixed (Hao 2018 Methods 'Covariate analysis': 'the allometric coefficients fixed at 0.75 for CL and 1 for V'); linear scaling on V/F with the same reference 70 kg and exponent 1.0 fixed. Study population mean (SD) 36.5 (17.4) kg, range 12.9-81.0 kg (Table 1).",
      source_name        = "WT"
    ),
    CYP3A5_EXPR = list(
      description        = "CYP3A5 expresser indicator: 1 if the patient carries at least one functional CYP3A5*1 allele (genotype *1/*1 or *1/*3), 0 if homozygous *3/*3.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A5 *3/*3 nonexpresser)",
      notes              = "Time-fixed germline genotype derived from rs776746 (CYP3A5 A6986G); the *1 (A) allele encodes functional CYP3A5 protein, the *3 (G) allele creates a cryptic splice site and yields nonfunctional protein. In the Hao 2018 cohort (n = 28) the genotype distribution was *3/*3 = 21 (75%), *1/*3 = 6 (21.4%), *1/*1 = 1 (3.6%); CYP3A5_EXPR = 1 for the 7 *1 carriers, 0 for the 21 nonexpressers. Multiplicative effect on CL/F as `e_cyp3a5_expr_cl ^ CYP3A5_EXPR` with `e_cyp3a5_expr_cl = 1.60` (60% higher CL/F in expressers). Hao 2018 Table 2 final model equation: F_CYP3A5 = theta_3 = 1.60 if *1/*1 or *1/*3; F_CYP3A5 = 1 if *3/*3.",
      source_name        = "FLAG1"
    )
  )

  population <- list(
    species               = "human",
    n_subjects            = 28L,
    n_studies             = 1L,
    n_observations        = 148L,
    age_range             = "2.7-17.3 years",
    age_median            = "9.4 years",
    age_mean_sd           = "9.5 (4.4) years",
    weight_range          = "12.9-81.0 kg",
    weight_median         = "30.0 kg",
    weight_mean_sd        = "36.5 (17.4) kg",
    sex_female_pct        = 32.1,
    race_ethnicity        = "Not reported in source paper (single-centre Chinese cohort, Children's Hospital of Hebei Province, Shijiazhuang).",
    disease_state         = "Paediatric patients with nephrotic syndrome treated with twice-daily oral tacrolimus (Prograf, Astellas, Japan) as initial immunosuppressant. Children under 18 years of age were enrolled prospectively; participants with concomitant medical conditions that posed unacceptable additional risk were excluded. Study registered at ClinicalTrials.gov NCT03347357.",
    dose_range            = "Starting dose 0.05 mg/kg twice daily (Methods: 'Tacrolimus (Prograf, Astellas, Japan), was administered orally at a dose of 0.05 mg kg-1 dose-1 twice daily'); actual administered tacrolimus per-dose range 1.0-8.0 mg twice daily (Table 1; median 4.0 mg, mean 3.8 mg SD 2.2 mg); weight-normalised range 0.0222-0.3876 mg/kg twice daily (Table 1; median 0.0909, mean 0.1199 SD 0.0860).",
    regions               = "China (Children's Hospital of Hebei Province, Shijiazhuang).",
    cyp3a5_distribution   = "*3/*3 = 21 (75%); *1/*3 = 6 (21.4%); *1/*1 = 1 (3.6%). Total n = 28. Determined from rs776746 (CYP3A5 A6986G) by TaqMan allelic discrimination (Hao 2018 Methods 'Analytical method of tacrolimus and genotyping').",
    sampling_window       = "Steady-state full concentration-time profiles were obtained during hospitalisation after a steady-state condition was achieved. Per Methods, samples were drawn predose and at 1, 2, 3, 6, 9, and 12 h after a tacrolimus dose. 148 tacrolimus concentrations in total (1-7 samples per patient).",
    assay                 = "Whole-blood tacrolimus measured by HPLC-MS/MS over the range 2.0-100 ng/mL (lower limit of quantification 2.0 ng/mL; intraday CV 4.4%, interday CV 7.2%). Concentrations below the LOQ (n = 11) were replaced with half-LOQ (1.0 ng/mL) in PK modelling per Methods 'Model building'.",
    target_trough_window  = "Recommended target predose concentration (C0) of 5-10 ng/mL for paediatric nephrotic syndrome (Hao 2018 Discussion and references [8, 54]); the simulation-based dosing recommendation in this paper is 0.10 mg/kg twice daily for *3/*3 nonexpressers and 0.25 mg/kg twice daily for *1 carriers.",
    notes                 = "Single-centre, prospective, open-label trial 2015-2017. The model is intended for paediatric nephrotic-syndrome patients receiving the twice-daily oral immediate-release formulation; it is NOT validated for adults, kidney transplant recipients, or the once-daily extended-release formulation. The Discussion (Limitations paragraph) explicitly notes that only internal validation was performed."
  )

  ini({
    # Final-model fixed-effect estimates from Hao 2018 Table 2.
    # Reference subject for the typical-value structural parameters:
    # WT = 70 kg, CYP3A5 *3/*3 nonexpresser (F_CYP3A5 = 1, the equation
    # reference). All apparent clearances in L/h; apparent volumes in L;
    # ka in 1/h.
    lka <- log(5.21)  ; label("Absorption rate constant ka (1/h)")                                # Hao 2018 Table 2 final ka = 5.21 h^-1 (SE 17.1%)
    lcl <- log(30.9)  ; label("Apparent oral clearance CL/F at WT = 70 kg, CYP3A5 *3/*3 (L/h)")    # Hao 2018 Table 2 final CL/F theta_2 = 30.9 L/h (SE 9.2%)
    lvc <- log(411)   ; label("Apparent volume of distribution V/F at WT = 70 kg (L)")             # Hao 2018 Table 2 final V/F theta_1 = 411 L (SE 20.9%)

    # Allometric exponents -- Hao 2018 Methods 'Covariate analysis':
    # "The allometric size approach was used by implementing the body weight
    # into the basic model (the allometric coefficients fixed at 0.75 for CL
    # and 1 for V)." Theory-based and held fixed during estimation.
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent of (WT/70) on CL/F (unitless; fixed at theory value)")   # Hao 2018 Methods 'Covariate analysis'
    e_wt_vc <- fixed(1)    ; label("Allometric exponent of (WT/70) on V/F (unitless; fixed at theory value)")    # Hao 2018 Methods 'Covariate analysis'

    # Covariate effect on CL/F -- Hao 2018 Table 2 final model equation:
    #   CL/F = theta_2 * (WT/70)^0.75 * F_CYP3A5
    # with F_CYP3A5 = theta_3 = 1.60 if *1/*1 or *1/*3; F_CYP3A5 = 1 if *3/*3.
    # Equivalent encoding: F_CYP3A5 = e_cyp3a5_expr_cl ^ CYP3A5_EXPR.
    e_cyp3a5_expr_cl <- 1.60 ; label("CYP3A5*1-carrier multiplicative factor on CL/F (expressers have 60% higher CL/F)")  # Hao 2018 Table 2 theta_3 = 1.60 (SE 23.8%)

    # Diagonal inter-individual variability on ka, V/F, CL/F. Hao 2018 Table 2
    # reports IIV as %CV from an exponential (log-normal) model. Variances on
    # the internal log-scale are computed as omega^2 = log(1 + CV^2):
    #   ka   CV 79.1% -> log(1 + 0.791^2) = log(1.62568) = 0.48581
    #   V/F  CV 99.4% -> log(1 + 0.994^2) = log(1.98804) = 0.68714
    #   CL/F CV 43.8% -> log(1 + 0.438^2) = log(1.19184) = 0.17555
    # The paper does not report any inter-eta correlations; the IIV is
    # encoded as a diagonal block.
    etalka ~ 0.48581   # Hao 2018 Table 2 IIV ka 79.1% CV (SE 52.8%)
    etalvc ~ 0.68714   # Hao 2018 Table 2 IIV V/F 99.4% CV (SE 34.3%)
    etalcl ~ 0.17555   # Hao 2018 Table 2 IIV CL/F 43.8% CV (SE 33.0%)

    # Residual unexplained variability. Hao 2018 Methods 'Model building':
    # "The proportional model best described residual variability." Table 2
    # labels the value 'Residual variability (exponential) (%) 25.9' -- the
    # 'exponential' label is the standard NONMEM additive-on-log-scale
    # parameterisation, which is equivalent to a proportional model on the
    # linear concentration scale for small variances.
    propSd <- 0.259 ; label("Proportional residual error (fraction)")   # Hao 2018 Table 2 Residual variability 25.9% (SE 22.2%)
  })

  model({
    # Body-weight scaling reference: 70 kg adult.
    wt70 <- WT / 70

    # CYP3A5 multiplier on CL/F: 1.60 if CYP3A5_EXPR = 1 (carrier); 1.0 if 0.
    f_cyp3a5 <- e_cyp3a5_expr_cl ^ CYP3A5_EXPR

    # Individual PK parameters with Hao 2018 covariate equations.
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * wt70 ^ e_wt_cl * f_cyp3a5
    vc <- exp(lvc + etalvc) * wt70 ^ e_wt_vc

    # One-compartment oral disposition with first-order absorption and
    # first-order elimination. Dose lands in `depot`; bioavailability is
    # implicit in the apparent CL/F and V/F parameterisation.
    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Tacrolimus whole-blood concentrations reported in ng/mL. Dose in mg,
    # vc in L, so central/vc is in mg/L = ug/mL; multiply by 1000 to convert
    # to ng/mL.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
