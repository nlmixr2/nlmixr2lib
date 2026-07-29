# Population pharmacokinetic model of oral piperaquine from an
# individual-participant-data meta-analysis of 11 clinical studies
# (8,776 plasma piperaquine samples from 728 individuals; Hoglund 2017,
# PLoS Med 14(1):e1002212; doi:10.1371/journal.pmed.1002212).

Hoglund_2017_piperaquine <- function() {
  description <- paste(
    "Population PK model for oral piperaquine in adults, children, and",
    "healthy volunteers across 11 pooled clinical studies (Hoglund 2017;",
    "individual-participant-data meta-analysis, n = 728). Two-transit-",
    "compartment absorption with kA = kTR feeding a three-compartment",
    "disposition model. Allometric body weight scaling on all clearances",
    "(fixed exponent 0.75) and volumes (fixed exponent 1.0) with reference",
    "weight 54 kg. Enzyme maturation function on elimination clearance",
    "(Hill-type sigmoid with MF50 = 0.575 y, Hill = 5.51). Dose-occasion",
    "effect adds 23.7% to relative bioavailability per consecutive",
    "dose. Bioavailability anchored at 1 with IIV. Predictions are venous",
    "plasma piperaquine base concentrations (ng/mL); a separately",
    "estimated capillary-to-venous scale of 106% is documented but not",
    "applied because only venous output is simulated.",
    sep = " "
  )
  reference <- paste(
    "Hoglund RM, Workman L, Edstein MD, Thanh NX, Quang NN, Zongo I, et al.",
    "(2017). Population Pharmacokinetic Properties of Piperaquine in",
    "Falciparum Malaria: An Individual Participant Data Meta-Analysis.",
    "PLoS Medicine 14(1):e1002212.",
    "doi:10.1371/journal.pmed.1002212.",
    sep = " "
  )
  vignette <- "Hoglund_2017_piperaquine"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling with fixed exponents 0.75 on apparent",
        "clearances (CL/F, Q1/F, Q2/F) and 1.0 on apparent volumes",
        "(Vc/F, Vp1/F, Vp2/F), centred on the published reference",
        "weight of 54 kg (Hoglund 2017 Table 3 footnote: 'Population",
        "estimates are given for a typical adult patient weighting 54 kg',",
        "and Methods page 5: 'Body weight was evaluated by adding it as an",
        "allometric function to all clearance (power of 0.75) and volume",
        "of distribution (power of 1) parameters'). Treated as time-fixed",
        "for forward simulation.",
        sep = " "
      ),
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Used in the enzyme maturation function on elimination clearance",
        "(Hoglund 2017 Eq 3 and Methods page 5):",
        "CL_i = theta_CL * AGE^Hill / (AGE^Hill + MF50^Hill).",
        "The maturation factor is ~0 for very young infants, 0.5 at",
        "AGE = MF50 = 0.575 y, and approaches 1 for older children and",
        "adults. Adults in the pooled cohort (16-55 y) are at full",
        "maturation; the maturation correction is active only below",
        "~2 y of age.",
        sep = " "
      ),
      source_name        = "AGE"
    ),
    OCC = list(
      description        = "Dose-occasion counter (1 = first dose, 2 = second, 3 = third)",
      units              = "(count)",
      type               = "count",
      reference_category = NULL,
      notes              = paste(
        "Hoglund 2017 retained a categorical dose-occasion covariate on",
        "relative bioavailability (Methods page 5, Results page 7,",
        "Table 3: 'Dose F (percent) = 23.7'). The published model",
        "interprets this as an increase in F between consecutive doses",
        "within a single 3-day course of treatment. This file encodes the",
        "increment additively: F_OCC = 1 + 0.237 * (OCC - 1), so",
        "OCC = 1 -> F = 1.000, OCC = 2 -> F = 1.237, OCC = 3 -> F = 1.474.",
        "The paper text is not fully explicit about additive-vs-",
        "compounded form; see the vignette Assumptions section.",
        "OCC is carried on each dose event row; observation rows inherit",
        "the most recent dose's OCC value but the value is only used by",
        "f(depot) at dose-administration time.",
        sep = " "
      ),
      source_name        = "OCC"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 728L,
    n_studies       = 11L,
    n_pregnant      = 36L,
    n_healthy       = 50L,
    age_range       = "0.56-55 years (pooled across 11 studies, Table 1)",
    age_min         = "0.56 y (infant cohort, Sambol 2009 and Tarning 2014)",
    weight_range    = "5.1-81 kg (pooled across 11 studies, Table 1)",
    sex_female_pct  = 43.3,
    pregnant_pct    = 4.9,
    disease_state   = paste(
      "Uncomplicated Plasmodium falciparum malaria infection in adults,",
      "children, and pregnant women (n = 678), plus healthy volunteer",
      "cohorts from Viet Nam (n = 50). Most pediatric data come from",
      "African sites (Burkina Faso, Kenya, Uganda); adult / pregnant",
      "cohorts span Thailand, Sudan, Viet Nam."
    ),
    dose_range      = paste(
      "Dihydroartemisinin-piperaquine fixed-dose combination administered",
      "as 1-3 daily doses (one course over 1, 2, or 3 days depending on",
      "the contributing study). Piperaquine administered as piperaquine",
      "phosphate, converted to piperaquine base by a 57.7% scale factor",
      "(Methods page 4). Manufacturer-recommended dose regimens",
      "(Sigma-Tau and Beijing Holley-Cotec) summarised in Table 2."
    ),
    regions         = "Thailand, Sudan, Viet Nam, Burkina Faso, Kenya, Uganda",
    venous_pct      = "Capillary samples used in pediatric African cohorts (Kenya, Burkina Faso, Uganda); venous in adult and pregnancy cohorts.",
    notes           = paste(
      "Demographics from Hoglund 2017 Table 1. Total observations 8,776",
      "(141 below the lower limit of quantification omitted). 11 clinical",
      "studies pooled from the WWARN repository. Inclusion criteria:",
      "prospective dihydroartemisinin-piperaquine study in patients with",
      "uncomplicated P. falciparum infection or in healthy volunteers,",
      "validated capillary and/or venous plasma piperaquine measurements",
      "(Methods page 3)."
    )
  )

  ini({
    # Structural parameters from Hoglund 2017 Table 3 ("Population
    # Estimate" column). Values are apparent (relative to F = 1) and given
    # on the linear scale; log() applied here for the nlmixr2 internal
    # log-scale. Typical values correspond to the paper's allometric
    # reference weight of 54 kg (Table 3 footnote: "Population estimates
    # are given for a 'typical' adult patient weighting 54 kg with acute
    # falciparum malaria").
    lcl  <- log(55.4)
    label("Apparent piperaquine elimination clearance CL/F at WT = 54 kg, fully mature (L/h)")
    # Hoglund 2017 Table 3: CL/F = 55.4 L/h (RSE 4.22%; 95% CI 51.2-60.6)

    lvc  <- log(2910)
    label("Apparent central volume of distribution Vc/F at WT = 54 kg (L)")
    # Hoglund 2017 Table 3: Vc/F = 2910 L (RSE 6.98%; 95% CI 2540-3340)

    lq   <- log(310)
    label("Apparent inter-compartmental clearance to first peripheral compartment Q1/F at WT = 54 kg (L/h)")
    # Hoglund 2017 Table 3: Q1/F = 310 L/h (RSE 8.03%; 95% CI 266-366)

    lvp  <- log(4910)
    label("Apparent first peripheral volume of distribution Vp1/F at WT = 54 kg (L)")
    # Hoglund 2017 Table 3: Vp1/F = 4910 L (RSE 5.85%; 95% CI 4390-5510)

    lq2  <- log(105)
    label("Apparent inter-compartmental clearance to second peripheral compartment Q2/F at WT = 54 kg (L/h)")
    # Hoglund 2017 Table 3: Q2/F = 105 L/h (RSE 4.98%; 95% CI 95.1-115)

    lvp2 <- log(30900)
    label("Apparent second peripheral volume of distribution Vp2/F at WT = 54 kg (L)")
    # Hoglund 2017 Table 3: Vp2/F = 30900 L (RSE 4.79%; 95% CI 28300-34200)

    lmtt <- log(2.11)
    label("Mean transit time of the 2-transit-compartment absorption chain MTT (h)")
    # Hoglund 2017 Table 3: MTT = 2.11 h (RSE 4.54%; 95% CI 1.94-2.30).
    # Number of transit compartments = 2 (fixed; Table 3). With kA = kTR
    # (Results page 6), the absorption chain depot -> transit1 -> transit2
    # -> central has 3 equal transitions; ktr = (NN + 1)/MTT = 3/MTT
    # (Savic & Karlsson 2007 convention, same as Hoglund_2012_piperaquine).

    # Relative bioavailability anchor fixed at 1 (Hoglund 2017 Methods
    # page 5: "The bioavailability was fixed to unity for the population";
    # Table 3: F (percent) = 100 fix).
    lfdepot <- fixed(log(1))
    label("Relative bioavailability F at OCC = 1 (unitless, fixed at 1)")
    # Hoglund 2017 Table 3: F = 100% fix

    # Allometric exponents fixed by the source paper. Hoglund 2017 Methods
    # page 5: "Body weight was evaluated by adding it as an allometric
    # function to all clearance (power of 0.75) and volume of distribution
    # (power of 1) parameters". Estimating the elimination-clearance
    # exponent did not significantly drop OFV (Results page 6), so the
    # fixed values are retained.
    e_wt_cl <- fixed(0.75)
    label("Allometric WT exponent on all apparent clearance parameters (CL/F, Q1/F, Q2/F)")
    # Hoglund 2017 Methods page 5: clearance power = 0.75 (fixed)

    e_wt_vc <- fixed(1.00)
    label("Allometric WT exponent on all apparent volume parameters (Vc/F, Vp1/F, Vp2/F)")
    # Hoglund 2017 Methods page 5: volume power = 1 (fixed)

    # Maturation function on CL (Hoglund 2017 Eq 3, Methods page 5,
    # Table 3 Covariate relationships): CL_i = theta_CL *
    # AGE^Hill / (AGE^Hill + MF50^Hill).
    mat_mf50 <- 0.575
    label("Maturation half-time on CL: age at which CL reaches 50% of fully mature value (years)")
    # Hoglund 2017 Table 3: MF50 = 0.575 y (RSE 13.6%; 95% CI 0.413-0.711)

    mat_hill <- 5.51
    label("Hill coefficient governing the slope of the CL maturation function (unitless)")
    # Hoglund 2017 Table 3: HillMF = 5.51 (RSE 29.6%; 95% CI 3.22-9.95; upper limit 10)

    # Dose-occasion effect on relative bioavailability. Hoglund 2017
    # Table 3 Covariate relationships: Dose F = 23.7%. Methods page 5 and
    # Results page 7 describe this as an additive proportional increase
    # in F at each subsequent dose. Encoded here additively in OCC:
    # F_OCC = 1 + 0.237 * (OCC - 1). See vignette Assumptions and
    # deviations for the additive-vs-compounded reading.
    e_doseocc_f <- 0.237
    label("Increment in relative bioavailability per additional dose occasion (fraction)")
    # Hoglund 2017 Table 3: Dose F = 23.7% (RSE 17.8%; 95% CI 15.8-32.5)

    # IIV. Hoglund 2017 Table 3 footnote: "Coefficients of variation for
    # inter-individual variability (IIV) and inter-occasion variability
    # (IOV) were calculated as 100 * (e^variance - 1)^(1/2)", so the
    # internal log-scale variance is recovered as omega^2 = log(CV^2 + 1).
    #
    #   CL   IIV 27.9%   -> omega^2 = log(0.279^2 + 1) = 0.074960
    #   Vc   IIV 67.0%   -> omega^2 = log(0.670^2 + 1) = 0.370805
    #   Vp1  IIV 24.0%   -> omega^2 = log(0.240^2 + 1) = 0.056002
    #   Q2   IIV 23.6%   -> omega^2 = log(0.236^2 + 1) = 0.054200
    #   Vp2  IIV 34.7%   -> omega^2 = log(0.347^2 + 1) = 0.113694
    #   MTT  IIV 38.0%   -> omega^2 = log(0.380^2 + 1) = 0.134880
    #   F    IIV 41.4%   -> omega^2 = log(0.414^2 + 1) = 0.158196
    #
    # No IIV on Q1 (Results page 6: "Inter-individual variability was
    # retained on all parameters (except inter-compartment clearance)";
    # the only inter-compartment clearance without IIV in Table 3 is Q1).
    # MTT IOV (46.4%) and F IOV (53.5%) are reported in Table 3 but not
    # encoded as separate eta slots here, matching the Hoglund_2012_
    # piperaquine precedent; see vignette Assumptions and deviations.
    etalcl     ~ 0.074960
    # Hoglund 2017 Table 3: IIV on CL/F = 27.9% CV (RSE 7.43%; 95% CI 23.7-31.7)

    etalvc     ~ 0.370805
    # Hoglund 2017 Table 3: IIV on Vc/F = 67.0% CV (RSE 15.5%; 95% CI 49.7-89.4)

    etalvp     ~ 0.056002
    # Hoglund 2017 Table 3: IIV on Vp1/F = 24.0% CV (RSE 44.2%; 95% CI 1.35-40.4)

    etalq2     ~ 0.054200
    # Hoglund 2017 Table 3: IIV on Q2/F = 23.6% CV (RSE 15.6%; 95% CI 16.2-30.2)

    etalvp2    ~ 0.113694
    # Hoglund 2017 Table 3: IIV on Vp2/F = 34.7% CV (RSE 7.21%; 95% CI 29.6-39.4)

    etalmtt    ~ 0.134880
    # Hoglund 2017 Table 3: IIV on MTT = 38.0% CV (RSE 15.8%; 95% CI 26.7-50.3)

    etalfdepot ~ 0.158196
    # Hoglund 2017 Table 3: IIV on F = 41.4% CV (RSE 8.65%; 95% CI 34.0-48.2)

    # Residual error. Hoglund 2017 Methods page 5: "The unknown
    # variability in concentration was described by an additive error on
    # the individually predicted logarithmic concentrations (i.e.,
    # equivalent to an exponential error on non-transformed
    # concentrations)." This NONMEM additive-on-log-scale residual maps
    # to a proportional residual in linear concentration space (see
    # references/parameter-names.md Residual error and the matching
    # comment in Hoglund_2012_piperaquine.R). Table 3 reports the
    # log-scale variance RUV = 0.115; the SD is sqrt(0.115) = 0.339116,
    # which equals the proportional CV to first order.
    propSd <- sqrt(0.115)
    label("Proportional residual SD for piperaquine plasma concentration (SD on log scale)")
    # Hoglund 2017 Table 3: RUV = 0.115 (variance, RSE 3.43%; 95% CI 0.108-0.123; epsilon shrinkage 14.6%)
  })

  model({
    # Maturation factor on elimination clearance (Hoglund 2017 Eq 3).
    # AGE = 0 -> maturation = 0; AGE = MF50 -> maturation = 0.5; large
    # AGE -> maturation -> 1.
    maturation_cl <- AGE^mat_hill / (AGE^mat_hill + mat_mf50^mat_hill)

    # Individual PK parameters. Allometric scaling on clearances
    # (exponent 0.75) and volumes (exponent 1) centred on the published
    # 54 kg reference. Maturation factor applied multiplicatively to CL.
    # Q1 carries no IIV (Hoglund 2017 Table 3).
    cl  <- exp(lcl  + etalcl)  * (WT / 54)^e_wt_cl * maturation_cl
    vc  <- exp(lvc  + etalvc)  * (WT / 54)^e_wt_vc
    q   <- exp(lq)             * (WT / 54)^e_wt_cl
    vp  <- exp(lvp  + etalvp)  * (WT / 54)^e_wt_vc
    q2  <- exp(lq2  + etalq2)  * (WT / 54)^e_wt_cl
    vp2 <- exp(lvp2 + etalvp2) * (WT / 54)^e_wt_vc

    # Mean transit time and chain rate constant. With NN = 2 transit
    # compartments and kA = kTR, the absorption chain
    # depot -> transit1 -> transit2 -> central has (NN + 1) = 3
    # equal-rate transitions; ktr = (NN + 1)/MTT = 3/MTT
    # (Savic & Karlsson 2007 convention).
    mtt <- exp(lmtt + etalmtt)
    ktr <- 3 / mtt

    # Three-compartment disposition micro-constants.
    kel <- cl  / vc
    k12 <- q   / vc
    k21 <- q   / vp
    k13 <- q2  / vc
    k31 <- q2  / vp2

    # ODE system: 2-transit-compartment absorption feeding into a
    # 3-compartment disposition model (Hoglund 2017 Fig 2). The same ktr
    # propagates the dose through the depot and two transit compartments
    # into central.
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(central)     <-  ktr * transit2 - kel * central -
                          k12 * central + k21 * peripheral1 -
                          k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central  - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central  - k31 * peripheral2

    # Relative bioavailability. Anchor lfdepot = log(1) is fixed (paper
    # convention) and dose-occasion-dependent multiplier
    # F_OCC = 1 + e_doseocc_f * (OCC - 1) increments F by 23.7% per
    # successive dose (Hoglund 2017 Table 3 Dose F coefficient). IIV
    # captured by etalfdepot.
    f_occ <- 1 + e_doseocc_f * (OCC - 1)
    f(depot) <- f_occ * exp(lfdepot + etalfdepot)

    # Piperaquine plasma concentration. Dose in mg of piperaquine base
    # (Methods page 4: piperaquine phosphate converted to piperaquine
    # base by a 57.7% scale factor). Vc in L gives central / vc in mg/L;
    # the 1000 factor converts to ng/mL for direct comparison against
    # Hoglund 2017 Table 3 (Cmax, day-7 concentrations all in ng/mL).
    # Predictions correspond to venous plasma piperaquine concentrations
    # (the matrix in which the typical-adult parameters in Table 3 were
    # reported). Capillary plasma values are 106% higher per Hoglund 2017
    # Methods page 5 and Table 3 ("Scale (percent) = 106"); the model
    # does not apply this scaling because only venous output is
    # simulated. Users who wish to predict capillary concentrations
    # should multiply Cc by 2.06 post-hoc.
    Cc <- 1000 * central / vc

    # Proportional residual error on the linear-concentration scale.
    Cc ~ prop(propSd)
  })
}
