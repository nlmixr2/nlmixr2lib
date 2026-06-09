Jorga_2000_tolcapone_nonfluctuators <- function() {
  description <- "Two-compartment population PK model with first-order absorption and absorption lag for tolcapone in parkinsonian patients with stable (non-fluctuating) levodopa response, with effects of creatinine clearance on clearance, serum protein on central volume, and concomitant food on bioavailability (Jorga 2000, nonfluctuator dataset, n=60)"
  reference <- paste(
    "Jorga K, Fotteler B, Banken L, Snell P, Steimer JL.",
    "Population pharmacokinetics of tolcapone in parkinsonian patients in dose finding studies.",
    "Br J Clin Pharmacol. 2000;49(1):39-48.",
    "doi:10.1046/j.1365-2125.2000.00113.x"
  )
  vignette <- "Jorga_2000_tolcapone"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance (raw Cockcroft-Gault, NOT BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form covariate on CL with reference 68 mL/min (population median, Table 1). Computed in the paper as CL_Cr = factor * (140 - age) * BW / serum_creatinine with factor = 1.23 for males / 1.04 for females, BW in kg, serum creatinine in umol/L (Lott & Hayton 1978; Methods). Raw Cockcroft-Gault without BSA normalization, same convention as Delattre 2010 amikacin. Source column 'CL_Cr' in the paper.",
      source_name        = "CL_Cr"
    ),
    TPRO = list(
      description        = "Total serum protein",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form covariate on Vc with reference 72 g/L (population median, Table 1). Source column 'Protein' in the paper.",
      source_name        = "Protein"
    ),
    FED = list(
      description        = "Fed-vs-fasted dose-record indicator (1 = dose taken with concomitant food, 0 = fasted)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = "Paper indicator I_Food. Multiplicative effect (1 + e_food_f * FED) applied on F1; e_food_f = -0.17 corresponds to a ~17% reduction in relative bioavailability in the fed state for the nonfluctuator cohort (Jorga 2000 Table 3 theta_Food = 0.83; Discussion: 15-20% reduction in nonfluctuators). Reference fasted F1 fixed at 0.6.",
      source_name        = "Food"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 60L,
    n_studies      = 1L,
    study_names    = c("Nonfluctuator 200 mg t.i.d. (n = 33)",
                       "Nonfluctuator 400 mg t.i.d. (n = 29)"),
    n_observations = 433L,
    age_range      = "47-83 years",
    age_median     = "67 years",
    weight_range   = "44-110 kg",
    weight_median  = "72 kg",
    lbm_range      = "34-75 kg",
    lbm_median     = "55 kg",
    crcl_range     = "41-141 mL/min",
    crcl_median    = "70 mL/min",
    sex_female_pct = NA_real_,
    race_ethnicity = c(Caucasian = 98.0, Black = 0.24, Asian = 0.73, Other = 1.46),
    disease_state  = "Parkinson's disease with stable (non-fluctuating) motor response to levodopa/AADC inhibitor therapy ('nonfluctuators')",
    dose_range     = "200 or 400 mg tolcapone three times daily for 6 weeks, with the levodopa dose reduced 33-43% on the first day of test treatment; ongoing levodopa-carbidopa (Sinemet) or levodopa-benserazide (Madopar) therapy",
    regions        = "One multicentre Phase II dose-finding study (Phase II nonfluctuator arm)",
    sampling       = "Sparse: 5-8 plasma samples per patient on 2-5 occasions; samples taken pre-dose, near Cmax, and during the decline phase. Only four nonfluctuator samples were obtained 0.5-1 h after dosing across the entire cohort (Discussion).",
    notes          = "Demographics from Table 1, 'Non-fluctuators' column (n=97 enrolled in the nonfluctuator arm; n=60 with usable PK)."
  )

  ini({
    # Structural PK parameters -- typical values for the 200 mg t.i.d.
    # nonfluctuator at the population-median creatinine clearance (68 mL/min)
    # and serum protein (72 g/L), under fasted conditions. From Jorga 2000
    # Table 3 (Non-fluctuator model, Final estimate column).
    lka      <- log(0.7);  label("Absorption rate constant (ka, 1/h)")                        # Table 3 (ka = 0.7 /h)
    lcl      <- log(4.5);  label("Apparent clearance for the reference covariates (CL/F1, L/h)")  # Table 3 (CL = 4.5 L/h)
    lvc      <- log(3.5);  label("Central volume of distribution for the reference covariates (Vc, L)")  # Table 3 (Vc = 3.5 L)
    lvp      <- log(24);   label("Peripheral volume of distribution (Vp, L)")                 # Table 3 (Vp = 24 L)
    lq       <- log(7.7);  label("Inter-compartmental clearance (Q, L/h)")                    # Table 3 (Q = 7.7 L/h)
    ltlag    <- log(0.4);  label("Absorption lag time (Tlag, h)")                             # Table 3 (Tlag = 0.4 h)
    # Fasted absolute bioavailability fixed at 0.6 from upstream IV/PO single-
    # dose study (Jorga et al. 1998, Eur J Clin Pharmacol 54:443-447; ref [22]).
    lfdepot  <- fixed(log(0.6)); label("Fasted absolute bioavailability (F1, fraction)")     # Methods: F1 fixed to 0.6 [ref 22]

    # Covariate effects on CL: power-form (CRCL/68)^e_crcl_cl
    e_crcl_cl   <-  1.19; label("Creatinine clearance power exponent on CL")    # Table 3 (theta_CLCr(CL) = 1.19)

    # Covariate effects on Vc: power-form (TPRO/72)^e_tpro_vc
    e_tpro_vc   <- -7.34; label("Total serum protein power exponent on Vc")     # Table 3 (theta_Protein(Vc) = -7.34)

    # Food effect on F1: (1 + e_food_f * FED) multiplier on fasted F1 = 0.6.
    # Paper theta_Food = 0.83 -> e_food_f = -0.17.
    e_food_f    <- -0.17; label("Relative change in F1 for fed vs fasted dosing (fraction)")  # Table 3 (theta_Food(F) = 0.83 -> 0.83 - 1)

    # Inter-individual variability (paper reports omega^2 on the log-normal
    # internal scale; sqrt(exp(omega^2) - 1) reproduces the parenthetical CV%).
    # omega^2 = 0.06 -> CV ~= 25% (Table 3)
    # omega^2 = 1.34 -> CV ~= 168% (Table 3); poorly identified, see vignette
    # Assumptions and deviations section.
    etalcl ~ 0.06
    etalvc ~ 1.34

    # Residual error -- combined log-normal multiplicative + additive (paper
    # equation DV = CP * exp(eps_mult) + eps_add; Table 3 reports both
    # variances for the nonfluctuator model). omega^2(eps_mult) = 0.18 -> SD
    # = sqrt(0.18) ~= 0.424; omega^2(eps_add) = 0.52 -> SD = sqrt(0.52)
    # ~= 0.721 ug/mL.
    propSd <- 0.424; label("Proportional residual error (fraction)")     # Table 3 (eps_mult variance = 0.18 -> SD = sqrt(0.18))
    addSd  <- 0.721; label("Additive residual error (ug/mL)")            # Table 3 (eps_add variance = 0.52 -> SD = sqrt(0.52))
  })

  model({
    # Individual PK parameters with covariate effects.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) *
          (CRCL / 68)^e_crcl_cl
    vc <- exp(lvc + etalvc) *
          (TPRO / 72)^e_tpro_vc
    vp <- exp(lvp)
    q  <- exp(lq)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 2-compartment first-order absorption with absorption lag, per Jorga 2000
    # Results "Final estimates" paragraph; the lag was retained in the
    # nonfluctuator model (Table 2: excluding tlag worsened OF by 12 units).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Fasted bioavailability with multiplicative food effect; applied to the
    # depot compartment (oral dosing). Absorption lag time.
    f(depot)    <- exp(lfdepot) * (1 + e_food_f * FED)
    alag(depot) <- exp(ltlag)

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
