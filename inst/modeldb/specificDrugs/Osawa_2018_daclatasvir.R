Osawa_2018_daclatasvir <- function() {
  description <- "One-compartment population PK model with first-order absorption and linear apparent clearance for oral daclatasvir (DCV, hepatitis C virus NS5A replication-complex inhibitor) in Japanese adults with chronic genotype-1 HCV infection (Osawa 2018). Fit by NONMEM 7.2 FOCE to 3801 plasma concentrations from 336 subjects across four trials (AI444021, AI444022, AI447017, AI447026) receiving daclatasvir 10 or 60 mg once daily, either with asunaprevir (the all-oral DUAL regimen) or with peginterferon-alfa plus ribavirin. Typical apparent clearance is 5.29 L/h and apparent central volume 64.2 L. Correlated inter-individual variability (correlation 0.94) is carried on CL/F and V/F, with independent IIV on Ka; the residual error is additive on the natural-log scale (log-transform-both-sides) and itself carries inter-individual variability on its magnitude. Four covariates survived backward elimination: baseline creatinine clearance (power 0.235) plus female sex and the peginterferon/ribavirin regimen (exponential) on CL/F, and baseline body weight (power 0.605) on V/F. All of these effects lie within or overlap the 80-125 percent boundaries, so the authors judged none of them clinically relevant."

  reference <- paste(
    "Osawa M, Ueno T, Ishikawa H, Imai Y, Garimella T. (2018).",
    "Population Pharmacokinetic Analysis for Daclatasvir and Asunaprevir",
    "in Japanese Subjects With Chronic Hepatitis C Virus Infection.",
    "The Journal of Clinical Pharmacology 58(11):1468-1478.",
    "doi:10.1002/jcph.1274.",
    sep = " "
  )

  vignette <- "Osawa_2018_daclatasvir_asunaprevir"

  # CL/F is reported in L/h and V/F in L (Table 2), so amounts in mg give
  # concentrations in mg/L = ug/mL.
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    CRCL_BASE = list(
      description = "Baseline creatinine clearance, time-fixed per subject. Enters CL/F as the median-normalised power function (CRCL_BASE / 86.48)^0.235.",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference 86.48 mL/min, stated in the text immediately below the daclatasvir covariate equations. Table 1 reports the cohort median as 86.5 mL/min (range 39.56-185.96), consistent with that reference to the printed precision. NOTE: the Figure 3 caption prints the reference as '88.48 mL/min'; that is a typographical error in the caption -- the equation text, the Table 1 median, and the printed 5th/95th-percentile fold-changes all agree on 86.48. The paper's forest plot quotes CL/F changes of approximately -10 percent and +10 percent at the 5th and 95th percentiles (51.36 and 144.42 mL/min); the encoded power function returns -11.5 percent and +12.8 percent. The covariate was retained in preference to age and baseline body weight on CL/F, which were moderately-to-highly correlated with it and tested separately.",
      source_name = "BCRCL (baseline creatinine clearance)"
    ),
    WT_BASE = list(
      description = "Baseline body weight, time-fixed per subject. Enters V/F as the median-normalised power function (WT_BASE / 56)^0.605.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference 56 kg, stated in the text below the covariate equations and equal to the Table 1 cohort median (range 36-93 kg). Registered as WT_BASE rather than WT because the source paper names the covariate BBWT -- baseline body weight -- and deliberately distinguishes baseline from time-varying laboratory covariates elsewhere in the same analysis (time-varying AST is used in the companion asunaprevir model). This is the Verrest 2023 / Yang 2025 shape of the canonical: baseline weight is preferred as the size descriptor, without Wahlby's within-subject delta form. The encoded exponent reproduces the paper's stated forest-plot values exactly: at the 5th and 95th percentiles of 42.8 and 78 kg, V/F is 15.0 percent lower and 22.2 percent higher than at the 56 kg reference, versus the paper's 'approximately 15 percent lower' and 'approximately 22 percent higher'.",
      source_name = "BBWT (baseline body weight)"
    ),
    SEXF = list(
      description = "Female sex indicator. 1 = female; 0 = male. Multiplies CL/F by exp(-0.110), i.e. daclatasvir apparent clearance is about 10 percent lower in women than in men.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male; SEXref is stated as male in the text below the covariate equations).",
      notes = "Time-fixed per subject. 63.7 percent of the daclatasvir cohort was female (214 of 336, Table 1). The source column orientation matches the canonical SEXF directly -- the paper's reference category is male, so the tabulated theta6 = -0.110 is the female-versus-male effect and needs no sign inversion. The paper groups this with the treatment effect: 'the daclatasvir CL/F was reduced by approximately 10 percent for subjects receiving pegIFN/RBV ... and for women relative to men'; exp(-0.110) = 0.896.",
      source_name = "SEX"
    ),
    CONMED_PEGIFN_RBV = list(
      description = "Concomitant peginterferon-alfa plus ribavirin backbone indicator. 1 = daclatasvir given with pegIFN/RBV; 0 = daclatasvir given with asunaprevir (the all-oral DUAL regimen). Multiplies CL/F by exp(-0.122), i.e. daclatasvir apparent clearance is about 11 percent lower on the pegIFN/RBV backbone.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (daclatasvir + asunaprevir DUAL regimen; TXref is stated as 'asunaprevir' in the text below the covariate equations).",
      notes = "Time-fixed per subject in this analysis -- each subject received one backbone regimen for the duration of the study. 21.1 percent of the daclatasvir cohort received daclatasvir + pegIFN/RBV (71 of 336) and 78.9 percent the DUAL regimen (265 of 336, Table 1). The paper labels the covariate 'TX / treatment description'; it is registered under the CONMED_ family because what it encodes is which drugs were co-administered alongside daclatasvir, and the canonical orientation is chosen to keep the DUAL regimen as the reference so the tabulated theta8 = -0.122 can be quoted unchanged. exp(-0.122) = 0.885.",
      source_name = "TX (treatment description)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "daclatasvir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "daclatasvir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 336L,
    n_studies = 4L,
    n_observations = 3801L,
    age_range = "21-75 years",
    age_median = "61 years",
    weight_range = "36-93 kg",
    weight_median = "56 kg",
    sex_female_pct = 63.7,
    race_ethnicity = "Japanese (all subjects; eligibility required Japanese men or women aged 20 years or older)",
    disease_state = "Chronic hepatitis C virus genotype-1 infection (99.4 percent genotype 1b, 0.6 percent genotype 1a). 6.5 percent had compensated cirrhosis. By prior-treatment history: 42.6 percent nonresponder / null responder / partial responder, 46.7 percent pegIFN/RBV-ineligible naive or intolerant, 10.7 percent treatment naive.",
    hepatic_function = "Baseline AST median 51.0 U/L (range 13-377); baseline ALT median 48.0 U/L (range 13-595). Neither cirrhosis status nor AST/ALT survived backward elimination for daclatasvir, consistent with the dedicated hepatic-impairment study.",
    renal_function = "Baseline creatinine clearance median 86.5 mL/min (range 39.56-185.96)",
    dose_range = "Daclatasvir 10 or 60 mg once daily, orally",
    regimens = "Daclatasvir + asunaprevir (DUAL, all-oral) in 78.9 percent; daclatasvir + peginterferon-alfa/ribavirin in 21.1 percent",
    regions = "Japan",
    notes = "Baseline demographics from Osawa 2018 Table 1 (daclatasvir panel, n = 336). Data pooled from four trials -- AI444021, AI444022, AI447017 and AI447026 -- listed in Supplementary Table 1, which is not on disk; the per-trial breakdown is therefore not reproduced here. Of 3808 collected pharmacokinetic records, 7 were excluded as below the lower limit of quantification, leaving 3801 for model development. Samples below LLOQ were discarded rather than imputed, and missing actual sampling times were imputed from nominal times."
  )

  ini({
    # =========================================================================
    # Structural fixed effects -- Osawa 2018 Table 2, daclatasvir panel.
    # The typical values are those of the reference subject: male, on the
    # daclatasvir + asunaprevir DUAL regimen, baseline creatinine clearance
    # 86.48 mL/min, baseline body weight 56 kg.
    # =========================================================================
    lcl <- log(5.29); label("Apparent clearance CL/F at the reference covariate values (L/h)") # Table 2 daclatasvir theta1 = 5.29 L/h (RSE 3.04%, bootstrap 95% CI 4.98-5.59)
    lvc <- log(64.2); label("Apparent central volume of distribution V/F at the reference covariate values (L)") # Table 2 daclatasvir theta2 = 64.2 L (RSE 3.18%, bootstrap 95% CI 60.1-68.2)
    lka <- log(0.865); label("First-order absorption rate constant Ka (1/h)") # Table 2 daclatasvir theta3 = 0.865 1/h (RSE 5.97%, bootstrap 95% CI 0.753-0.974)

    # =========================================================================
    # Covariate effects -- Osawa 2018 Table 2 and the two covariate equations
    # printed in the Results 'Population PK Model Development' subsection:
    #
    #   CL/F_TV = CL/F_TV,ref * (BCRCLb / BCRCLref)^CL/F_BCRCL
    #                         * exp(SEX * CL/F_SEX + TX * CL/F_TX)
    #   V/F_TV  = V/F_TV,ref  * (BBWTb / BBWTref)^V/F_BBWT
    #
    # with BBWTref 56 kg, BCRCLref 86.48 mL/min, SEXref male, TXref
    # asunaprevir (the DUAL regimen).
    # =========================================================================
    e_crcl_base_cl <- 0.235; label("Power exponent of baseline creatinine clearance on CL/F, normalised to 86.48 mL/min (unitless)") # Table 2 daclatasvir theta10 = 0.235 (RSE 19.7%, bootstrap 95% CI 0.142-0.333)
    e_sexf_cl <- -0.110; label("Exponential effect of female sex on CL/F, male reference (unitless)") # Table 2 daclatasvir theta6 = -0.110 (RSE 24.9%, bootstrap 95% CI -0.197 to -0.053)
    e_conmed_pegifn_rbv_cl <- -0.122; label("Exponential effect of the peginterferon/ribavirin backbone on CL/F, DUAL-regimen reference (unitless)") # Table 2 daclatasvir theta8 = -0.122 (RSE 27.5%, bootstrap 95% CI -0.189 to 0.132 as printed)
    e_wt_base_vc <- 0.605; label("Power exponent of baseline body weight on V/F, normalised to 56 kg (unitless)") # Table 2 daclatasvir theta14 = 0.605 (RSE 16.3%, bootstrap 95% CI 0.405-0.974)

    # =========================================================================
    # Inter-individual variability -- Osawa 2018 Table 2 'Random effects'.
    # Table 2 footnote b: diagonal elements are printed as variance (standard
    # deviation) and off-diagonal elements as covariance (correlation), so the
    # first number in each cell is the value nlmixr2 wants. The CL/F-V/F block
    # is c(var_CL, cov, var_V); check: 0.141 / (0.394 * 0.381) = 0.939, which
    # reproduces the printed correlation of 0.941 to rounding.
    # =========================================================================
    etalcl + etalvc ~ c(0.155, 0.141, 0.145) # Table 2 daclatasvir omega1,1 = 0.155 (SD 0.394); omega1,2 = 0.141 (correlation 0.941); omega2,2 = 0.145 (SD 0.381)
    etalka ~ 0.756 # Table 2 daclatasvir omega3,3 = 0.756 (SD 0.869)

    # =========================================================================
    # Residual error -- Osawa 2018 Table 2 'Residual error' plus the omega4,4
    # row that Table 2 places in the 'Random effects' block under the symbol
    # sigma. The Methods residual model is log-transform-both-sides:
    #
    #   ln(y_ij) = ln(yhat_ij) + theta_ADD * epsilon_ij
    #
    # so the residual is additive on the natural-log scale, which is exactly
    # nlmixr2's lnorm() error structure and the canonical `expSd` name.
    # theta4 is that log-scale SD; omega4,4 is inter-individual variability on
    # its magnitude, i.e. the NONMEM W = THETA(4) * EXP(ETA(4)) construct, so
    # the per-subject residual SD is expSd * exp(etaexpSd). The two rows carry
    # distinct standard errors (RSE 2.72% for theta4, 18.0% for omega4,4),
    # confirming they are separately estimated parameters rather than one value
    # printed twice.
    # =========================================================================
    expSd <- 0.375; label("Residual error SD, additive on the natural-log scale (log-scale SD)") # Table 2 daclatasvir theta4 = 0.375 (RSE 2.72%, bootstrap 95% CI 0.358-0.408)
    etaexpSd ~ 0.107 # Table 2 daclatasvir omega4,4 = 0.107 (SD 0.327), listed under 'Random effects' with the symbol sigma; IIV on the residual-error magnitude
  })

  model({
    # -----------------------------------------------------------------------
    # 1. Individual PK parameters.
    #
    # The two covariate equations are transcribed verbatim from the Results
    # section; the reference constants 86.48 mL/min and 56 kg are the values
    # stated in the sentence immediately following them.
    # -----------------------------------------------------------------------
    cl <- exp(lcl + etalcl) *
      (CRCL_BASE / 86.48)^e_crcl_base_cl *
      exp(e_sexf_cl * SEXF + e_conmed_pegifn_rbv_cl * CONMED_PEGIFN_RBV)
    vc <- exp(lvc + etalvc) * (WT_BASE / 56)^e_wt_base_vc
    ka <- exp(lka + etalka)

    # -----------------------------------------------------------------------
    # 2. Micro-constant and ODE system. One compartment with first-order
    # absorption from an oral depot and linear elimination. Bioavailability is
    # not identifiable from oral-only data and is absorbed into CL/F and V/F,
    # so no f(depot) term appears.
    # -----------------------------------------------------------------------
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # -----------------------------------------------------------------------
    # 3. Observation and residual error. The per-subject residual SD carries
    # its own eta (see the ini() note on omega4,4).
    # -----------------------------------------------------------------------
    Cc <- central / vc
    expSdInd <- expSd * exp(etaexpSd)
    Cc ~ lnorm(expSdInd)
  })
}
