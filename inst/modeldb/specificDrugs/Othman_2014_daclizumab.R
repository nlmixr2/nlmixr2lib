Othman_2014_daclizumab <- function() {
  description <- "Two-compartment population PK model with first-order subcutaneous absorption and lag time for daclizumab high-yield process (HYP) in healthy volunteers (Othman 2014)"
  reference <- "Othman AA, Tran JQ, Tang MT, Dutta S. Population Pharmacokinetics of Daclizumab High-Yield Process in Healthy Volunteers: Integrated Analysis of Intravenous and Subcutaneous, Single- and Multiple-Dose Administration. Clin Pharmacokinet. 2014;53(10):907-918. doi:10.1007/s40262-014-0159-9"
  vignette <- "Othman_2014_daclizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used for allometric scaling of CL, Q, Vc, and Vp with reference 70 kg; exponents estimated (not fixed to 0.75/1.0): 0.54 for CL and Q, 0.64 for Vc and Vp (Othman 2014 Table 2).",
      source_name        = "WT"
    ),
    DOSE_50MG = list(
      description        = "Record-level indicator for the 50 mg SC dose (1 = 50 mg SC, 0 = any other SC dose or any IV dose)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (100, 150, 200, or 300 mg SC dose, or any IV dose)",
      notes              = "Othman 2014 estimated two separate absolute bioavailabilities for SC administration because dose-normalized exposure at 50 mg SC was lower than at higher SC doses. F = 0.84 for the clinical 100-300 mg SC range and F = 0.57 for the 50 mg SC cohort. Encoded as a dose-record-level covariate so that `e_dose_50mg_f * DOSE_50MG` reduces bioavailability only on 50 mg SC dose events. For Phase III clinical simulations (150 mg SC every 4 weeks) leave DOSE_50MG = 0. Derived from AMT; not a column in the source dataset.",
      source_name        = "(derived from AMT)"
    )
  )

  population <- list(
    n_subjects     = 70,
    n_studies      = 3,
    study_names    = c("Study 1 (SC single dose 50/150/300 mg)",
                       "Study 2 (SC multiple dose 100 or 200 mg biweekly with 200 mg load)",
                       "Study 3 (IV single dose 200 or 400 mg)"),
    age_range      = "18-66 years",
    age_mean       = "35.9 years (SD 15.4)",
    weight_range   = "55.7-127 kg",
    weight_mean    = "77.7 kg (SD 16.1)",
    bmi_range      = "18.1-44.2 kg/m^2",
    bmi_mean       = "26.7 kg/m^2 (SD 5.3)",
    sex_female_pct = 50.7,
    race_ethnicity = c(`Caucasian/Hispanic` = 88.7, Asian = 9.9, Other = 1.4),
    disease_state  = "Healthy volunteers (Phase I safety / tolerability / PK)",
    dose_range     = "50, 150, or 300 mg SC single dose; 100 or 200 mg SC every 2 weeks (200 mg loading dose); 200 or 400 mg IV single dose",
    regions        = "Australia (all three studies conducted at CMAX, Adelaide)",
    n_observations = 925,
    notes          = "Healthy-volunteer integrated analysis pooling three Phase I studies (N = 71 dosed, 70 analyzable after exclusion of one subject with a likely 150 mg SC dosing error). Baseline demographics per Othman 2014 Table 1."
  )

  ini({
    # Structural PK parameters — reference values for a 70 kg adult.
    # Paper reports values in L/h and /h; converted to /day by x 24 so this
    # model matches the nlmixr2lib convention (`units$time = "day"`).
    lka      <- log(0.009 * 24); label("Absorption rate constant (Ka, 1/day; 0.009 /h SC)")         # Table 2 (ka = 0.009 /h)
    lcl      <- log(0.010 * 24); label("Clearance for a 70 kg adult (CL, L/day; 10 mL/h)")          # Table 2 (CL = 0.010 L/h)
    lvc      <- log(3.89);       label("Central volume of distribution for a 70 kg adult (Vc, L)")  # Table 2 (Vc = 3.89 L)
    lvp      <- log(2.52);       label("Peripheral volume of distribution for a 70 kg adult (Vp, L)")  # Table 2 (Vp = 2.52 L)
    lq       <- log(0.044 * 24); label("Inter-compartmental clearance for a 70 kg adult (Q, L/day; 44 mL/h)")  # Table 2 (Q = 0.044 L/h)
    lfdepot  <- log(0.84);       label("Subcutaneous bioavailability for 100-300 mg doses (F, fraction)")  # Table 2 (F_100-300mg = 84%)
    lalag    <- log(2 / 24);     label("Absorption lag time for SC doses (Tlag, day; 2 h)")         # Table 2 (Lag time = 2.0 h)

    # Allometric exponents — estimated, not fixed (Othman 2014 reported that
    # fixing to the classical 0.75 / 1.0 was statistically less favorable than
    # the estimated values).
    allo_cl <- 0.54; label("Allometric exponent on CL and Q (SFCL, unitless)")    # Table 2 (SFCL = 0.54)
    allo_v  <- 0.64; label("Allometric exponent on Vc and Vp (SFV, unitless)")    # Table 2 (SFV = 0.64)

    # Relative bioavailability effect for the 50 mg SC cohort (multiplicative
    # fractional change on F). Derived as F_50mg / F_100-300mg - 1 = 0.57/0.84 - 1.
    e_dose_50mg_f <- -0.32143; label("Relative change in bioavailability for 50 mg SC vs 100-300 mg SC (fraction)")  # Table 2 (F_50mg = 57%, F_100-300mg = 84%)

    # Inter-individual variability (omega^2 = log(CV^2 + 1)). Values are the SC
    # cohort IIV from Table 2 because the clinical route of administration for
    # daclizumab HYP is SC (Phase III regimen is 150 mg SC every 4 weeks); IV
    # cohort IIV was smaller (15% CV on CL, 16% CV on Vc) and is noted in the
    # vignette's Assumptions and deviations section. The 50 mg SC cohort had an
    # inflated CL ISV (scaling factor 1.82 over SC baseline, Table 2); that
    # inflation is specific to the low dose and is not carried into the
    # simulation-oriented typical variance here.
    # ka (58% CV)  <-> omega^2 = log(1 + 0.58^2) = 0.29003
    # cl SC (27%)  <-> omega^2 = log(1 + 0.27^2) = 0.07038
    # corr(ka,cl SC) = -0.72 -> cov = -0.72 * sqrt(0.29003 * 0.07038) = -0.10290
    etalka + etalcl ~ c(0.29003,
                        -0.10290, 0.07038)
    # Vc SC (31%) <-> omega^2 = log(1 + 0.31^2) = 0.09175
    etalvc ~ 0.09175

    # Residual error. Paper reports a combined proportional + additive error
    # with the proportional component as exp(eps1) (log-normal multiplicative)
    # and Table 2 giving "Proportional error (%) = 22 (rprop x 100)" and
    # "Additive error SD (ug/mL) = 0.33". propSd = 0.22 maps to rprop on
    # nlmixr2's proportional scale; for |rprop| << 1 the log-normal and linear
    # proportional forms are numerically indistinguishable.
    propSd <- 0.22; label("Proportional residual error (fraction)")    # Table 2 (rprop = 0.22)
    addSd  <- 0.33; label("Additive residual error (ug/mL)")           # Table 2 (radd = 0.33 ug/mL)
  })
  model({
    # Individual PK parameters with allometric weight scaling (reference 70 kg)
    ka  <- exp(lka + etalka)
    cl  <- exp(lcl + etalcl) * (WT / 70)^allo_cl
    vc  <- exp(lvc + etalvc) * (WT / 70)^allo_v
    vp  <- exp(lvp)          * (WT / 70)^allo_v
    q   <- exp(lq)           * (WT / 70)^allo_cl

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Subcutaneous absorption lag time; bioavailability for SC (reduced on
    # 50 mg SC records via DOSE_50MG). For IV administration, users route the
    # dose directly to `central` (cmt = 2) so f(depot) and alag(depot) have no
    # effect.
    alag(depot) <- exp(lalag)
    f(depot)    <- exp(lfdepot) * (1 + e_dose_50mg_f * DOSE_50MG)

    # Concentration: dose in mg, volumes in L -> mg/L = ug/mL
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
