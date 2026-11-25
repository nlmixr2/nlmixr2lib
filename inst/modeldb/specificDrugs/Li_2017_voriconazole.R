Li_2017_voriconazole <- function() {
  description <- "One-compartment population pharmacokinetic model with first-order absorption (ka fixed at 1.1 1/h) and linear elimination for voriconazole in Chinese adult renal transplant recipients, fitted to therapeutic-drug-monitoring trough concentrations (Li 2017); aspartate aminotransferase enters as a power covariate on clearance and CYP2C19 phenotype (rapid / extensive / intermediate / poor metabolizer) as a multiplicative covariate on volume of distribution."
  reference <- "Li ZW, Peng FH, Yan M, Liang W, Liu XL, Wu YQ, Lin XB, Tan SL, Wang F, Xu P, Fang PF, Liu YP, Xiang DX, Zhang BK. Impact of CYP2C19 Genotype and Liver Function on Voriconazole Pharmacokinetics in Renal Transplant Recipients. Ther Drug Monit. 2017;39(4):422-428. doi:10.1097/FTD.0000000000000425"
  vignette <- "Li_2017_voriconazole"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot = list(analyte = "voriconazole", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "voriconazole", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    AST = list(
      description = "Serum aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL normalised to 33 U/L, the cohort mean AST of Li 2017 Table 1 (33.0 +/- 48.1 U/L, range 7.9-311.5): CL = 4.76 * (AST/33)^-0.23. Whether AST was entered as a per-sample (time-varying) or per-subject value is not stated.",
      source_name = "AST"
    ),
    CYP2C19_RM = list(
      description = "CYP2C19 rapid-metabolizer phenotype indicator (*1/*17)",
      units = "(binary)",
      type = "binary",
      reference_category = "0. Li 2017's own reference category for V is the rapid metabolizer: the typical V of 22.47 L applies when none of EM, IM or PM is 1.",
      notes = "The rapid metabolizer (*1/*17, n = 2) is the source-paper reference group for V (Li 2017 Table 5 footnote). The extensive-metabolizer indicator the paper's equation needs is derived inside model() as 1 - CYP2C19_RM - CYP2C19_IM - CYP2C19_PM - CYP2C19_UM, so exactly one of the four phenotype rows must be coded per subject.",
      source_name = "RM"
    ),
    CYP2C19_IM = list(
      description = "CYP2C19 intermediate-metabolizer phenotype indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (rapid metabolizer is the reference for V in Li 2017)",
      notes = "Li 2017 pools *1/*2, *1/*3 and *2/*17 into IM (Methods 'DNA Sequencing and CYP2C19 Genetic Polymorphism'; Table 3). V multiplier 1 + 4.67 = 5.67 versus the rapid-metabolizer reference.",
      source_name = "IM"
    ),
    CYP2C19_PM = list(
      description = "CYP2C19 poor-metabolizer phenotype indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (rapid metabolizer is the reference for V in Li 2017)",
      notes = "Li 2017 defines PM as *2/*2, *2/*3 or *3/*3; all five PM subjects were *2/*2 (Table 3). V multiplier 1 + 3.30 = 4.30 versus the rapid-metabolizer reference.",
      source_name = "PM"
    ),
    CYP2C19_UM = list(
      description = "CYP2C19 ultrarapid-metabolizer phenotype indicator (*17/*17)",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Li 2017 defines a UM (*17/*17) category but enrolled no UM subjects, so the paper gives UM no V coefficient. The model assigns UM to the rapid-metabolizer reference group (V = 22.47 L), the other *17 gain-of-function stratum, rather than to the *1/*1 extensive-metabolizer group; this is an extrapolation outside the fitted cohort.",
      source_name = "UM"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 56L,
    n_studies = 1L,
    n_observations = 125L,
    age_range = "18-60 years",
    age_mean = "40 +/- 8 years",
    weight_range = "33-78 kg",
    weight_mean = "55 +/- 10 kg",
    sex_female_pct = 30.4,
    race_ethnicity = c(Chinese = 100),
    cyp2c19_phenotype = c(RM_pct = 3.6, EM_pct = 42.9, IM_pct = 44.6, PM_pct = 8.9),
    hepatic_function = "AST 33.0 +/- 48.1 U/L (7.9-311.5); ALT 45.7 +/- 108.7 U/L; TBIL 8.9 +/- 4.2 umol/L",
    renal_function = "CLcr 51.8 +/- 24.7 mL/min (11.6-111.6)",
    co_medication = "Proton pump inhibitor 46.4%; glucocorticoid 69.6%",
    disease_state = "Adult renal transplant recipients receiving voriconazole for prevention or treatment of invasive fungal infection; 78.6% within the first year after transplantation.",
    dose_range = "Usual starting dose 200 mg twice daily, adjusted after the first trough; intravenous and oral dosage forms both used. Trough concentrations (Cmin) drawn 30 minutes before the next dose.",
    regions = "Single centre: Second Xiangya Hospital, Central South University, Changsha, China.",
    notes = "Retrospective TDM study, September 2015 - June 2016 (ChiCTR-RRC-16008314). 125 trough concentrations (median 2.18, range 0.16-9.59 ug/mL) from 56 patients; fitted with Phoenix NLME 1.4 (FOCE ELS). Demographics from Li 2017 Table 1; final estimates from Tables 4 and 5."
  )

  ini({
    # Absorption: ka fixed (Methods 'Structural Model': "the absorption rate
    # was fixed to a value of 1.1 h", citing refs 18-20; the unit is 1/h).
    lka <- fixed(log(1.1)); label("Absorption rate constant (1/h)") # Methods 'Structural Model': ka fixed 1.1

    lcl <- log(4.76); label("Clearance at AST = 33 U/L (L/h)") # Table 5 final model: theta_CL = 4.76 L/h
    lvc <- log(22.47); label("Volume of distribution, CYP2C19 rapid-metabolizer reference (L)") # Table 5 final model: theta_V = 22.47 L

    # V = theta_V * [1 + theta_1*(EM=1)] * [1 + theta_2*(IM=1)] * [1 + theta_3*(PM=1)]
    e_cyp2c19_em_vc <- 2.21; label("Fractional increase in V for CYP2C19 EM vs RM (unitless)") # Table 5 final model: theta_1 = 2.21
    e_cyp2c19_im_vc <- 4.67; label("Fractional increase in V for CYP2C19 IM vs RM (unitless)") # Table 5 final model: theta_2 = 4.67
    e_cyp2c19_pm_vc <- 3.30; label("Fractional increase in V for CYP2C19 PM vs RM (unitless)") # Table 5 final model: theta_3 = 3.30

    # CL = theta_CL * (AST/33)^theta_4
    e_ast_cl <- -0.23; label("Power exponent of AST/33 on CL (unitless)") # Table 5 final model: theta_4 = -0.23

    # Table 4/5 print omega as an SD; the Abstract's final-model equation
    # prints the variances exp(0.96) and exp(0.14), which equal 0.98^2 and
    # 0.37^2.
    etalvc ~ 0.9604 # Table 5: omega_V = 0.98 (SD); variance 0.98^2 = 0.96 as in the Abstract
    etalcl ~ 0.1369 # Table 5: omega_CL = 0.37 (SD); variance 0.37^2 = 0.14 as in the Abstract

    # Proportional error Cobs = Cpred * (1 + eps). Table 5 sigma = 0.39 is the
    # SD; Results 'Population Pharmacokinetic Analysis' prints the variance
    # (1 + 0.15), 0.39^2 = 0.152.
    propSd <- 0.39; label("Proportional residual error (fraction)") # Table 5 final model: sigma = 0.39
  })

  model({
    # Extensive metabolizer (*1/*1) is the phenotype not flagged by any other
    # indicator; the rapid metabolizer is Li 2017's reference for V.
    cyp2c19_em <- 1 - CYP2C19_RM - CYP2C19_IM - CYP2C19_PM - CYP2C19_UM

    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (AST / 33)^e_ast_cl
    vc <- exp(lvc + etalvc) *
      (1 + e_cyp2c19_em_vc * cyp2c19_em) *
      (1 + e_cyp2c19_im_vc * CYP2C19_IM) *
      (1 + e_cyp2c19_pm_vc * CYP2C19_PM)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # mg / L = ug/mL, the paper's concentration unit.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
