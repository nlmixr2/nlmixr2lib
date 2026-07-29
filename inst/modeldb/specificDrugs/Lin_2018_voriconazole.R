Lin_2018_voriconazole <- function() {
  description <- "One-compartment population pharmacokinetic model with first-order absorption for intravenous and oral voriconazole in Chinese adult renal transplant recipients receiving therapeutic drug monitoring (Lin 2018); CYP2C19 phenotype enters as a covariate on clearance, postoperative time as a covariate on oral bioavailability, and body weight as a power-form covariate on volume of distribution."
  reference <- "Lin XB, Li ZW, Yan M, Zhang BK, Liang W, Wang F, Xu P, Xiang DX, Xie XB, Yu SJ, Lan GB, Peng FH. Population pharmacokinetics of voriconazole and CYP2C19 polymorphisms for optimizing dosing regimens in renal transplant recipients. Br J Clin Pharmacol. 2018;84(7):1587-1597. doi:10.1111/bcp.13595"
  vignette <- "Lin_2018_voriconazole"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on V centered at the Lin 2018 cohort-median 56.1 kg (Table 1) with estimated exponent 1.30: V = TVV * (WT/56.1)^1.30. The reference is the dataset median rather than a canonical 70 kg adult value; the exponent is data-estimated, not the allometric-theory 1.0.",
      source_name        = "WT"
    ),
    CYP2C19_IM = list(
      description        = "CYP2C19 intermediate-metabolizer phenotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (EM/UM phenotype; the implicit reference when both CYP2C19_IM and CYP2C19_PM are 0)",
      notes              = "Lin 2018 Table 3 uses PM as the source-paper reference category and reports IM and EM as multiplicative log-scale shifts on CL relative to PM (theta_2 = 0.45 for IM vs PM, theta_3 = 0.80 for EM vs PM). The model file is reparameterized to the canonical convention where EM/UM (both CYP2C19_IM and CYP2C19_PM = 0) is the implicit reference; the equivalent EM/UM-referenced shifts are e_im_cl = theta_2 - theta_3 = -0.35 and e_pm_cl = -theta_3 = -0.80. See the in-file source-trace comments for the derivation. Lin 2018 enrolled 49 IM subjects (46.7%); IM genotypes pooled by Lin 2018 were *1/*2, *1/*3, *2/*17 (per Methods 'DNA purification and CYP2C19 genotyping').",
      source_name        = "CYP2C19 IM phenotype"
    ),
    CYP2C19_PM = list(
      description        = "CYP2C19 poor-metabolizer phenotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (EM/UM phenotype; the implicit reference when both CYP2C19_IM and CYP2C19_PM are 0)",
      notes              = "Companion to CYP2C19_IM. See CYP2C19_IM notes for the reparameterization rationale (Lin 2018 reports PM as the source-paper reference; the model file uses EM/UM as the canonical implicit reference). Lin 2018 enrolled 12 PM subjects (11.4%); PM genotypes pooled by Lin 2018 were *2/*2, *2/*3, *3/*3.",
      source_name        = "CYP2C19 PM phenotype"
    ),
    POD = list(
      description        = "Post-operative day (days since renal transplantation)",
      units              = "days",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Lin 2018 parameterizes the postoperative-time effect on F as four mutually-exclusive categorical bins (POT1, POT2, POT3, POT4) defined by month boundaries: <= 1 month (reference), 1-6 months, 6-12 months, > 1 year. The exact day-cutpoint values are not stated in Lin 2018. The model file uses the conventional 30 / 180 / 365 day boundaries to derive binary bin indicators from POD inside model(). The POT effect multiplies F = exp(lfdepot) and therefore applies only when the dose enters via the depot compartment (oral route).",
      source_name        = "POT (categorical 1/2/3/4 in Lin 2018; reconstructed from POD via 30/180/365 day cutpoints)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 105L,
    n_studies      = 1L,
    n_observations = 342L,
    age_range      = "18-58 years",
    age_mean       = "36 +/- 9 years (median 36)",
    weight_range   = "38.9-87.5 kg",
    weight_median  = "56.1 kg",
    sex_female_pct = 20.0,
    race_ethnicity = c(Chinese = 100),
    cyp2c19_phenotype = c(EM_pct = 41.5, IM_pct = 46.2, PM_pct = 11.3, RM_pct = 0.9),
    postoperative_time_distribution = c(
      WithinOneMonth_pct      = 31.4,
      OneToSixMonths_pct      = 33.3,
      SixToTwelveMonths_pct   = 21.0,
      OverOneYear_pct         = 14.3
    ),
    disease_state  = "Adult renal transplant recipients receiving voriconazole (intravenous or oral) for prevention or treatment of invasive fungal infections after kidney transplantation. All patients received tacrolimus or cyclosporine as primary immunosuppression. CYP2C19 genotyping and routine therapeutic drug monitoring were performed.",
    dose_range     = "Intravenous and oral voriconazole administered twice daily after a loading dose, with maintenance dose adjusted by surgeons per clinical response and TDM. Initial dose per voriconazole manufacturer package insert. 28 (26.7%) of patients received oral dosing only and 77 (73.3%) switched from intravenous to oral after stabilization. Trough samples (Cmin) collected 30 min before the next dose at steady state (day 5 or later, or day 2 with loading doses).",
    regions        = "Single center: Second Xiangya Hospital, Central South University, Changsha, Hunan, China.",
    notes          = "Prospective single-center clinical study, March 2016 - January 2017, Chinese Clinical Trial Registry ChiCTR-IPR-16008277. 129 patients screened, 106 included in the dataset (105 retained in the PPK analysis; one rapid metabolizer was excluded due to insufficient sample size). 342 voriconazole plasma concentrations measured by automated two-dimensional HPLC. CYP2C19 alleles tested: *2, *3, *17 (allele frequencies 29.2%, 5.2%, 0.5%). Baseline demographics per Lin 2018 Table 1; final-model parameter estimates per Lin 2018 Table 3."
  )

  ini({
    # Structural parameters. Reference subject for the typical-value
    # equation is a CYP2C19 EM/UM phenotype (the canonical implicit
    # reference when both CYP2C19_IM and CYP2C19_PM = 0) at WT = 56.1
    # kg (cohort median) and POD <= 30 days (POT1 bin).

    # Absorption: ka fixed at 1.1/h per Lin 2018 Methods 'Structural
    # model', citing the literature reference [21] (Hyland 2003).
    lka <- fixed(log(1.1)); label("Absorption rate constant (1/h), fixed")  # Lin 2018 Methods: "The absorption rate constant was fixed at 1.1 h-1 based on the literature report [21]"

    # CL: Lin 2018 Table 3 reports theta_CL = 2.88 L/h with PM as the
    # paper's reference category, with exp(0.80) for EM relative to PM
    # (Table 3 final-model equation: CL = theta_CL * exp(PM=0) *
    # exp[theta_2 * (IM=1)] * exp[theta_3 * (EM=1)] * exp(eta_CL)).
    # The model file is reparameterized to the canonical implicit
    # reference EM/UM (both CYP2C19 indicators = 0); the typical-value
    # CL for the EM/UM reference is 2.88 * exp(0.80) = 6.41 L/h.
    lcl <- log(6.41); label("Clearance for the CYP2C19 EM/UM reference (L/h)")  # Derived: Lin 2018 Table 3 theta_CL * exp(theta_3) = 2.88 * exp(0.80) = 6.4096

    # V: Lin 2018 Table 3 theta_V = 169.27 L (the typical value at the
    # WT = 56.1 kg reference).
    lvc <- log(169.27); label("Volume of distribution at WT = 56.1 kg (L)")  # Lin 2018 Table 3 final model: theta_V = 169.27 L

    # F: Lin 2018 Table 3 theta_F = 0.58 (58%) at the POT1 (POD <= 30
    # days) reference bin.
    lfdepot <- log(0.58); label("Oral bioavailability at POD <= 30 days (fraction)")  # Lin 2018 Table 3 final model: theta_F = 0.58

    # Power exponent for WT on V (Lin 2018 Table 3 theta_1). Note that
    # the estimated value 1.30 is not the allometric-theory 1.0; it is
    # the value estimated from data in the renal-transplant cohort.
    e_wt_vc <- 1.30; label("Power exponent for WT on V (unitless)")  # Lin 2018 Table 3 final model: theta_1 = 1.30

    # CYP2C19 phenotype effects on CL, reparameterized to the canonical
    # EM/UM-implicit-reference convention (see CYP2C19_IM notes in
    # covariateData). Lin 2018 reports PM-referenced exponents
    # theta_2 = 0.45 (IM vs PM) and theta_3 = 0.80 (EM vs PM); the
    # equivalent EM/UM-referenced exponents are:
    #   e_im_cl = theta_2 - theta_3 = 0.45 - 0.80 = -0.35
    #   e_pm_cl = -theta_3 = -0.80
    # These reproduce Lin 2018's per-phenotype CL exactly:
    #   EM/UM (CYP2C19_IM=0, CYP2C19_PM=0): exp(log(6.41)) = 6.41 L/h
    #   IM    (CYP2C19_IM=1, CYP2C19_PM=0): exp(log(6.41) - 0.35) = 4.52 L/h, matches 2.88 * exp(0.45) = 4.52
    #   PM    (CYP2C19_IM=0, CYP2C19_PM=1): exp(log(6.41) - 0.80) = 2.88 L/h, matches 2.88
    e_im_cl <- -0.35; label("Log-scale CL shift for CYP2C19_IM vs EM/UM (unitless)")  # Derived: Lin 2018 Table 3 theta_2 - theta_3 = 0.45 - 0.80
    e_pm_cl <- -0.80; label("Log-scale CL shift for CYP2C19_PM vs EM/UM (unitless)")  # Derived: Lin 2018 Table 3 -theta_3 = -0.80

    # Postoperative-time effects on F. Lin 2018 Table 3 reports the
    # exponents theta_4, theta_5, theta_6 for POT2, POT3, POT4
    # relative to the POT1 (POD <= 30 days) reference. POT3 and POT4
    # happen to share the same point estimate (theta_5 = theta_6 = 0.57).
    e_pot2_fdepot <- 0.43; label("Log-scale F shift for POD 30-180 days vs <= 30 (unitless)")  # Lin 2018 Table 3 final model: theta_4 = 0.43
    e_pot3_fdepot <- 0.57; label("Log-scale F shift for POD 180-365 days vs <= 30 (unitless)")  # Lin 2018 Table 3 final model: theta_5 = 0.57
    e_pot4_fdepot <- 0.57; label("Log-scale F shift for POD > 365 days vs <= 30 (unitless)")    # Lin 2018 Table 3 final model: theta_6 = 0.57

    # IIV. Lin 2018 reports exponential IIV (P_i = P_pop * exp(eta_i))
    # with reported approximate-CV% values 39% on V, 42% on CL, and
    # 22% on F (Lin 2018 Table 3 and Discussion: "the interindividual
    # variability of V, CL and F was as high as 39%, 42% and 22%").
    # Using the standard NONMEM convention CV% ~= sqrt(omega^2) * 100,
    # the internal variances are 0.39^2 = 0.1521 for V, 0.42^2 = 0.1764
    # for CL, and 0.22^2 = 0.0484 for F. The exact log-normal variance
    # log(1 + CV^2) gives 0.1419, 0.1631, 0.0473 respectively; the
    # approximate form is used here as the more common popPK reporting
    # convention.
    etalvc      ~ 0.1521  # Lin 2018 Table 3: omega_V  = 0.39 (39% CV); var = 0.39^2
    etalcl      ~ 0.1764  # Lin 2018 Table 3: omega_CL = 0.42 (42% CV); var = 0.42^2
    etalfdepot  ~ 0.0484  # Lin 2018 Table 3: omega_F  = 0.22 (22% CV); var = 0.22^2

    # Residual error. Lin 2018 Methods 'Statistical model' explicitly
    # selects the additive error model (Cobs = Cpred + epsilon) over
    # the proportional / combined / exponential alternatives. Sigma is
    # interpreted as the SD on the linear concentration scale in
    # ug/mL, the units reported throughout the paper's Cmin tables.
    addSd <- 0.57; label("Additive residual error (ug/mL)")  # Lin 2018 Table 3 final model: sigma = 0.57
  })

  model({
    # Derive POT bin indicators from POD (post-operative day, days).
    # Lin 2018 bins POT by month boundaries (<= 1 mo / 1-6 mo / 6-12
    # mo / > 1 yr). The model file uses 30 / 180 / 365 day cutpoints
    # to translate the conceptual month bins into a continuous-POD
    # decomposition; POT1 (POD <= 30) is the implicit reference when
    # pot_bin2 = pot_bin3 = pot_bin4 = 0.
    pot_bin2 <- (POD > 30) * (POD <= 180)
    pot_bin3 <- (POD > 180) * (POD <= 365)
    pot_bin4 <- (POD > 365)

    # Individual PK parameters
    ka     <- exp(lka)
    cl     <- exp(lcl + e_im_cl * CYP2C19_IM + e_pm_cl * CYP2C19_PM + etalcl)
    vc     <- exp(lvc + e_wt_vc * log(WT / 56.1) + etalvc)
    fdepot <- exp(lfdepot + e_pot2_fdepot * pot_bin2 + e_pot3_fdepot * pot_bin3 + e_pot4_fdepot * pot_bin4 + etalfdepot)

    # One-compartment with first-order oral absorption; intravenous
    # doses bypass the depot compartment by dosing directly into
    # central.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - (cl / vc) * central

    # Oral bioavailability applies only when the dose enters via the
    # depot compartment.
    f(depot) <- fdepot

    # Observation. Concentration units mg / L = ug/mL match the
    # paper's plasma-Cmin reporting scale.
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
