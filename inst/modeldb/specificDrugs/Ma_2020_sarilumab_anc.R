Ma_2020_sarilumab_anc <- function() {
  description <- "Indirect-response PopPK/PD model for absolute neutrophil count (ANC) following subcutaneous sarilumab in adults with rheumatoid arthritis (Ma 2020). Sarilumab concentrations drive stimulation of ANC elimination (margination); PK backbone is Xu 2019."
  reference   <- paste(
    "Ma L, Xu C, Paccaly A, Kanamaluru V. Population Pharmacokinetic-Pharmacodynamic Relationships of Sarilumab Using Disease Activity Score 28-Joint C-Reactive Protein and Absolute Neutrophil Counts in Patients with Rheumatoid Arthritis. Clin Pharmacokinet. 2020;59(11):1451-1466. doi:10.1007/s40262-020-00899-7 (PMID 32451909). PK backbone:",
    "Xu C, Su Y, Paccaly A, Kanamaluru V. Population Pharmacokinetics of Sarilumab in Patients with Rheumatoid Arthritis. Clin Pharmacokinet. 2019;58(11):1455-1467. doi:10.1007/s40262-019-00765-1 (PMID 31055792)."
  )
  vignette <- "Ma_2020_sarilumab_anc"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L", ANC = "10^9/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline weight; used on Kout with a power function centred at the median 71 kg (Ma 2020 Table 4 footnote b).",
      source_name        = "WT"
    ),
    SMOKE = list(
      description        = "Current smoker at baseline",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-smoker)",
      notes              = "Power-form covariate on baseline ANC: BASE = BASE_typ * 1.15^SMOKE (Ma 2020 Table 4).",
      source_name        = "Smoking"
    ),
    PRICORT = list(
      description        = "Prior corticosteroid treatment at baseline",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no prior corticosteroid)",
      notes              = "Power-form covariate on Emax: Emax = Emax_typ * 0.819^PRICORT (Ma 2020 Table 4).",
      source_name        = "PRICORT"
    )
  )

  population <- list(
    n_subjects     = 1672,
    n_studies      = 5,
    age_mean_sd    = "51.7 (12.1) years",
    weight_mean_sd = "74.1 (18.7) kg",
    weight_median  = "71 kg (reference for Kout)",
    sex_female_pct = 82.2,
    race_ethnicity = c(Caucasian = 84.8),
    disease_state  = "Moderate-to-severe rheumatoid arthritis (MTX-IR or TNF-IR)",
    dose_range     = "Sarilumab 50-150 mg SC qw, or 100-200 mg SC q2w (labelled dose is 200 mg q2w with step-down to 150 mg q2w for neutropenia)",
    regions        = "Pooled phase I-III studies (NCT01011959, NCT01061736, NCT01709578, NCT01768572)",
    baseline_anc   = "Mean 5.38 x 10^9/L (CV 32.1%)",
    smokers_pct    = 14.2,
    mtx_pct        = 97.7,
    pri_cort_pct   = 63.8,
    notes          = "Demographics from Ma 2020 Table 2 (ANC final dataset, n = 1672). Baseline ANC summary from Table 4."
  )

  ini({
    # ---------------------------------------------------------------------
    # PK backbone (typical-value parameters from Xu 2019 Table 3; reference
    # patient is a 71 kg female, albumin 38 g/L, CrCl 100 mL/min, baseline
    # CRP 14.2 mg/L, ADA-negative, drug product DP1 or DP3). PK covariates
    # (weight, albumin, CrCl, BLCRP, ADA, DP2, sex) are NOT included here;
    # they live in the companion Xu_2019_sarilumab.R file. This file
    # simulates ANC for the typical PK patient.
    # ---------------------------------------------------------------------
    lvmax <- log(8.06);   label("Sarilumab Michaelis-Menten maximum elimination rate (Vm, mg/day)")      # Xu 2019 Table 3
    lkm   <- log(0.939);  label("Sarilumab Michaelis-Menten constant (Km, mg/L)")                        # Xu 2019 Table 3
    lvc   <- log(2.08);   label("Apparent central volume of distribution (Vc/F, L)")                     # Xu 2019 Table 3
    lcl   <- log(0.260);  label("Apparent linear clearance (CLO/F, L/day)")                              # Xu 2019 Table 3
    lka   <- log(0.136);  label("Absorption rate constant (Ka, 1/day)")                                  # Xu 2019 Table 3
    lq    <- log(0.156);  label("Apparent intercompartmental clearance (Q/F, L/day)")                    # Xu 2019 Table 3
    lvp   <- log(5.23);   label("Apparent peripheral volume of distribution (Vp/F, L)")                  # Xu 2019 Table 3

    # ---------------------------------------------------------------------
    # PD parameters (Ma 2020 Table 4 ANC PopPK/PD)
    # ---------------------------------------------------------------------
    lbase  <- log(5.38); label("Typical baseline ANC (10^9/L)")                            # Ma 2020 Table 4
    lemax  <- log(1.50); label("Maximum drug-induced stimulation of ANC elimination (Emax, unitless)")  # Ma 2020 Table 4
    lec50  <- log(10.3); label("Sarilumab concentration at 50% of Emax (EC50, mg/L)")      # Ma 2020 Table 4

    # NOTE: Ma 2020 Table 4 prints this parameter as
    #       "211 (1.67-2.88)"
    # which is almost certainly a decimal-point typo: the bootstrap
    # 95% CI brackets 2.11 but not 211 (a 100x discrepancy). A neutrophil
    # turnover rate constant of ~2/day (half-life ~8 h) is physiologically
    # consistent; 211/day is not. We implement Kout=2.11 as the point
    # estimate; the CI from Table 4 is preserved as published.
    # See vignette Assumptions and deviations for full discussion.
    lkout <- log(2.11);   label("First-order ANC elimination rate constant (Kout, 1/day)")  # Ma 2020 Table 4 (corrected decimal typo; published bootstrap median "211" -> 2.11)

    lgamma <- log(0.862); label("Hill coefficient for sigmoidicity of ANC effect (unitless)")  # Ma 2020 Table 4

    # Covariate effect parameters (power-form exponents)
    allo_kout       <- 0.875; label("Weight exponent on Kout (ref 71 kg, unitless)")                           # Ma 2020 Table 4
    e_smoke_base    <- 1.15;  label("Smoking multiplier on baseline ANC (power-form: BASE * 1.15^SMOKE)")      # Ma 2020 Table 4
    e_pricort_emax  <- 0.819; label("Prior corticosteroid multiplier on Emax (power-form: Emax * 0.819^PRICORT)")  # Ma 2020 Table 4

    # ---------------------------------------------------------------------
    # Inter-individual variability (omega^2 = log(CV^2 + 1))
    # ---------------------------------------------------------------------
    # PK IIV (Xu 2019 Table 3)
    # Vm 32.4% CV, CLO 55.3% CV with correlation r = -0.566 (block); Vc 37.3%
    # CV; Ka 32.1% CV.
    # Block covariance = r * sqrt(omega2_vm * omega2_cl)
    #                  = -0.566 * sqrt(0.09982 * 0.26691) = -0.09238
    etalvmax + etalcl ~ c(0.09982,
                         -0.09238, 0.26691)                                                    # Xu 2019 Table 3
    etalvc  ~ 0.13023  # 37.3% CV                                                              # Xu 2019 Table 3
    etalka  ~ 0.09807  # 32.1% CV                                                              # Xu 2019 Table 3

    # PD IIV (Ma 2020 Table 4)
    # BASE 32.1%, Emax 61.9%, EC50 36.9%, Kout 227% (very high; preserved as
    # reported), gamma 80.4%.
    etalbase  ~ 0.09807                                                                        # Ma 2020 Table 4 (32.1% CV)
    etalemax  ~ 0.32423                                                                        # Ma 2020 Table 4 (61.9% CV)
    etalec50  ~ 0.12762                                                                        # Ma 2020 Table 4 (36.9% CV)
    etalkout  ~ 1.81707                                                                        # Ma 2020 Table 4 (227% CV)
    etalgamma ~ 0.49851                                                                        # Ma 2020 Table 4 (80.4% CV)

    # ---------------------------------------------------------------------
    # Residual error
    # ---------------------------------------------------------------------
    # PK: Xu 2019 fitted log-transformed concentrations with an additive error
    # on the log-scale of variance sigma^2 = 0.395; this maps to a
    # proportional error of sqrt(0.395) = 0.6285 in linear space.
    CcpropSd  <- 0.6285; label("Proportional residual error on sarilumab concentration (fraction)")   # Xu 2019 Table 3
    # PD: Ma 2020 Table 4 reports a proportional residual error of 28.2%.
    ANCpropSd <- 0.282;  label("Proportional residual error on ANC (fraction)")                      # Ma 2020 Table 4
  })

  model({
    # -------------------------------------------------------------------
    # PK: 2-compartment with first-order absorption and parallel linear
    # + Michaelis-Menten elimination from the central compartment.
    # -------------------------------------------------------------------
    vmax <- exp(lvmax + etalvmax)
    km   <- exp(lkm)
    vc   <- exp(lvc + etalvc)
    cl   <- exp(lcl + etalcl)
    ka   <- exp(lka + etalka)
    q    <- exp(lq)
    vp   <- exp(lvp)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    Cc  <- central / vc

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1 - vmax * Cc / (km + Cc)
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # -------------------------------------------------------------------
    # PD: indirect-response model, sarilumab stimulates ANC elimination
    # (margination of functional neutrophils out of circulation). The
    # `effect` compartment holds ANC; its initial value is baseline.
    #   d/dt(effect) = Kin - Kout * (1 + Eff) * effect
    #   Eff          = Emax * Cc^gamma / (EC50^gamma + Cc^gamma)
    # At baseline (Cc = 0): Kin = Kout * BASE, so effect(0) = BASE.
    # -------------------------------------------------------------------
    base  <- exp(lbase  + etalbase)  * (e_smoke_base)^SMOKE
    emax  <- exp(lemax  + etalemax)  * (e_pricort_emax)^PRICORT
    ec50  <- exp(lec50  + etalec50)
    kout  <- exp(lkout  + etalkout)  * (WT / 71)^allo_kout
    gamma <- exp(lgamma + etalgamma)

    kin <- kout * base

    eff <- emax * Cc^gamma / (ec50^gamma + Cc^gamma)

    d/dt(effect) <- kin - kout * (1 + eff) * effect
    effect(0)    <- base

    ANC <- effect

    Cc  ~ prop(CcpropSd)
    ANC ~ prop(ANCpropSd)
  })
}
