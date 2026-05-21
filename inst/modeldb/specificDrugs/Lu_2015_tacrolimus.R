Lu_2015_tacrolimus <- function() {
  description <- "Two-compartment population PK model with first-order absorption and lag time for oral tacrolimus in pooled Chinese healthy volunteers and adult orthotopic liver-transplant recipients (Lu 2015). Apparent peripheral volume V3/F is fixed at the healthy-volunteer-only estimate (916 L). Apparent clearance CL/F is reduced multiplicatively in liver-transplant recipients and further modulated by an exponential serum ALT effect that applies only to the transplant cohort."
  reference   <- "Lu YX, Su QH, Wu KH, Ren YP, Li L, Zhou TY, Lu W. A population pharmacokinetic study of tacrolimus in healthy Chinese volunteers and liver transplant patients. Acta Pharmacol Sin. 2015;36(2):281-288. doi:10.1038/aps.2014.110"
  vignette    <- "Lu_2015_tacrolimus"
  units       <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    DIS_HEALTHY = list(
      description        = "Healthy-participant cohort indicator: 1 = healthy Chinese volunteer (single 2 mg oral dose under bioequivalence study), 0 = adult orthotopic liver-transplant recipient on chronic oral tacrolimus immunosuppression.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (liver-transplant recipient)",
      notes              = "Time-fixed per subject. Lu 2015 uses the source name `SubPop` with the opposite orientation (SubPop = 0 for healthy volunteer, SubPop = 1 for liver-transplant patient); the canonical `DIS_HEALTHY` flip is applied via the identity SubPop = 1 - DIS_HEALTHY when reproducing the paper's Eq. 10. The reference complement here is the liver-transplant patient cohort (n = 112), not the union of all non-healthy indications referenced elsewhere in the register.",
      source_name        = "SubPop"
    ),
    ALT = list(
      description        = "Serum alanine aminotransferase activity at the time of the observation. Time-varying per subject in the transplant cohort (daily clinical-chemistry panel); near-normal and not load-bearing in the healthy-volunteer cohort.",
      units              = "IU/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Lu 2015 Eq. 10 enters ALT via an exponential effect Exp(ALT / 40 * theta_ALT) applied only to liver-transplant patients (gated on SubPop = 1, i.e., DIS_HEALTHY = 0). The normalisation factor 40 IU/L is interpreted as the clinical upper limit of normal for serum ALT used by the central laboratory; the paper does not state it explicitly but uses ALT / 40 verbatim in the printed equation. Patient ALT mean 146.4 +/- 290.0 IU/L (range 5 - 6300, Table 1); healthy-volunteer ALT mean 33.0 +/- 20.9 IU/L (range 8.5 - 125.3).",
      source_name        = "ALT"
    )
  )

  population <- list(
    species              = "human",
    n_subjects           = 152L,
    n_studies            = 1L,
    age_range            = "24 - 78 years",
    age_median           = "Liver-transplant patients 58.4 +/- 11.6 years; healthy volunteers 28.7 +/- 3.47 years (Table 1 means +/- SD)",
    weight_range         = "44 - 97 kg",
    weight_median        = "Liver-transplant patients 69.0 +/- 11.8 kg; healthy volunteers 62.5 +/- 6.46 kg (Table 1 means +/- SD)",
    sex_female_pct       = 17.1,
    sex_distribution     = "Liver-transplant patients: 86 male / 26 female (76.8% / 23.2%). Healthy volunteers: 40 male / 0 female. Pooled cohort 126 male / 26 female (82.9% / 17.1%).",
    race_ethnicity       = "Chinese (single-country: Beijing, China; healthy volunteers from PLA Second Artillery General Hospital, transplant recipients from General Hospital of Armed Police Forces).",
    disease_state        = "Pooled cohort of (i) 112 adult Chinese orthotopic liver-transplant recipients in the early postoperative period (POD 2 - 137, mean 19.4 days) receiving the triple immunosuppressive regimen of tacrolimus + mycophenolate mofetil + corticosteroids; (ii) 40 healthy adult Chinese male volunteers in a bioequivalence study receiving a single 2 mg oral dose. Underlying liver disease in the transplant arm: primary hepatic carcinoma (50.9%), liver cirrhosis or chronic severe hepatitis (30.6%), hepatic cancer recurrence (6.5%), alcoholic cirrhosis (4.6%), other (< 5%).",
    dose_range           = "Liver-transplant patients: initial 0.05 mg/kg/day in two divided oral doses, titrated by TDM (concentration target window not stated). Healthy volunteers: single 2 mg oral capsule.",
    formulations         = "Tacrolimus (Prograf, FK506; Astellas Pharma China) oral capsules 0.5 mg or 1 mg; identical formulation across cohorts.",
    n_observations       = 1951L,
    n_observations_breakdown = "1100 trough samples from 112 liver-transplant patients (sparse, microparticle enzyme immunoassay, 1.5 - 30 ug/L linear range); 851 dense post-dose samples from 40 healthy volunteers (HPLC-MS, 0.1 - 25 ug/L linear range).",
    regions              = "China (single-country)",
    notes                = "Cohort and demographic details from Lu 2015 Table 1. Of the 112 transplant patients, 36 were >= 65 years (elderly subgroup); age showed no significant CL/F effect in either the forward-inclusion or backward-elimination steps (Table 3) and was not retained in the final model. Co-medications in the transplant arm included methylprednisolone (perioperative bolus 500 - 1000 mg, oral taper to 4 mg/day) and MMF dispersible tablets 750 mg bid; the analysis does not test these as covariates."
  )

  ini({
    # Structural PK estimates from Lu 2015 Table 2 final-model column.
    # Time in hours, apparent clearances (CL/F, Q/F) in L/h, apparent volumes
    # (V2/F = central, V3/F = peripheral) in L, ka in 1/h, ALAG1 in h.
    #
    # Re-parameterisation note (see vignette Assumptions and deviations):
    # Lu 2015 Eq. 10 writes CL/F = theta1 * theta7^SubPop * Exp(ALT/40 * theta8)^SubPop * eta
    # with theta1 = 32.8 (healthy reference), theta7 = 0.562, theta8 = -0.0237,
    # and SubPop = 1 for liver-transplant patient, 0 for healthy volunteer.
    # The canonical covariate-columns register orients DIS_HEALTHY = 1 for
    # healthy and 0 for patient, with DIS_HEALTHY = 0 as the reference category.
    # The encoding below holds the patient cohort as the structural typical:
    #   lcl              = log(theta1 * theta7) = log(18.4); patient CL/F at ALT = 0.
    #   e_dis_healthy_cl = 1 / theta7         = 1.78; healthy/patient CL ratio.
    #   e_alt_cl         = theta8             = -0.0237; gated on (1 - DIS_HEALTHY).
    # This is mathematically equivalent to the paper's Eq. 10.
    lcl   <- log(32.8 * 0.562)  ; label("Apparent CL/F (L/h); typical liver-transplant patient at ALT = 0")  # Lu 2015 Table 2 final: CL/F_healthy = 32.8 L/h x SubPop multiplier theta7 = 0.562 -> 18.4 L/h
    lvc   <- log(22.7)          ; label("Apparent central volume V2/F (L)")                                  # Lu 2015 Table 2 final: V2/F = 22.7 L
    lq    <- log(76.3)          ; label("Apparent inter-compartmental clearance Q/F (L/h)")                  # Lu 2015 Table 2 final: Q/F = 76.3 L/h
    lvp   <- fixed(log(916))    ; label("Apparent peripheral volume V3/F (L); fixed from healthy-only fit")  # Lu 2015 Table 2 final: V3/F = 916 L (fixed per Methods "Population pharmacokinetic model development")
    lka   <- log(0.419)         ; label("Absorption rate constant ka (1/h)")                                 # Lu 2015 Table 2 final: ka = 0.419 1/h
    ltlag <- log(0.404)         ; label("Absorption lag time ALAG1 (h)")                                     # Lu 2015 Table 2 final: ALAG1 = 0.404 h

    # Covariate effects on CL/F (Lu 2015 Eq. 10).
    e_dis_healthy_cl <- 1 / 0.562  ; label("Multiplicative factor on CL/F for DIS_HEALTHY = 1 (healthy vs patient at ALT = 0)")  # Lu 2015 Table 2 final: CL_SubPop = theta7 = 0.562; re-expressed as 1/theta7 for the DIS_HEALTHY = 1 -> healthy direction
    e_alt_cl         <- -0.0237    ; label("Exponent in exp(e_alt_cl * ALT / 40) on CL/F for DIS_HEALTHY = 0 (liver-transplant patient)")  # Lu 2015 Table 2 final: CL_ALT = theta8 = -0.0237 (gated on SubPop = 1 = 1 - DIS_HEALTHY per Eq. 10)

    # Inter-individual variability -- Lu 2015 reports IIV as %CV with an
    # exponential / log-normal eta model (Eqs. 4 - 9); omega^2 = log(CV^2 + 1).
    etalcl ~ 0.1965  # Lu 2015 Table 2 final: CL-IIV = 46.6% CV  -> log(0.466^2 + 1) = 0.1965
    etalvc ~ 0.2841  # Lu 2015 Table 2 final: V2-IIV = 57.3% CV  -> log(0.573^2 + 1) = 0.2841
    etalq  ~ 0.1919  # Lu 2015 Table 2 final: Q-IIV  = 46.0% CV  -> log(0.460^2 + 1) = 0.1919
    etalvp ~ 0.6280  # Lu 2015 Table 2 final: V3-IIV = 93.5% CV  -> log(0.935^2 + 1) = 0.6280
    # No eta on ka (Lu 2015 Table 2 final: ka-IIV = 0%*, fixed at zero) and no eta on ALAG1 (Eq. 15 has no random term).

    # Residual unexplained variability (Lu 2015 final residual model: pure proportional, no additive component).
    propSd <- 0.398  ; label("Proportional residual error (fraction)")  # Lu 2015 Table 2 final: Proportional = 39.8%
  })

  model({
    # Individual PK parameters. CL/F follows Lu 2015 Eq. 10 with the canonical
    # DIS_HEALTHY orientation applied: paper's SubPop = 1 - DIS_HEALTHY.
    #   patient (DIS_HEALTHY = 0): cl = exp(lcl + etalcl) * exp(e_alt_cl * ALT/40)
    #   healthy (DIS_HEALTHY = 1): cl = exp(lcl + etalcl) * e_dis_healthy_cl
    cl   <- exp(lcl + etalcl) *
            e_dis_healthy_cl ^ DIS_HEALTHY *
            exp(e_alt_cl * (ALT / 40) * (1 - DIS_HEALTHY))
    vc   <- exp(lvc + etalvc)
    q    <- exp(lq + etalq)
    vp   <- exp(lvp + etalvp)
    ka   <- exp(lka)
    tlag <- exp(ltlag)

    # Micro-constants for the two-compartment system (mirrors NONMEM ADVAN4).
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment oral PK with first-order absorption and absorption lag.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    # Tacrolimus whole-blood concentration (ng/mL = ug/L). Doses are in mg, vc
    # in L, so central / vc has units mg/L; multiply by 1000 to ng/mL.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
