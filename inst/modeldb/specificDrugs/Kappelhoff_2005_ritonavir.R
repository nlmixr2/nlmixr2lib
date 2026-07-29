Kappelhoff_2005_ritonavir <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption,",
    "an absorption lag time, and first-order elimination for oral",
    "ritonavir in HIV-1-infected adults (186 patients, 1228 plasma",
    "concentrations; Kappelhoff 2005). Concomitant lopinavir is the",
    "only retained covariate and multiplies apparent oral clearance",
    "by 2.72-fold (power form: CL/F = exp(lcl) * 2.72^CONMED_LPV).",
    "Inter-individual variability on apparent CL/F, V/F, and ka, with",
    "correlated etas for V and ka (rho = 0.868). Residual error has",
    "a single 15.4% proportional component and a mixture-model",
    "additive component (subpopulation P1, 64.8% of subjects: 0.0600",
    "mg/L; subpopulation P2, 35.2%: 0.199 mg/L), gated by the binary",
    "covariate MIX_LARGE_RUV. Interoccasion variability on apparent",
    "bioavailability (59.1% in the source) is not propagated -- see",
    "the validation vignette Assumptions and deviations section."
  )
  reference <- paste(
    "Kappelhoff BS, Huitema ADR, Crommentuyn KML, Mulder JW, Meenhorst PL,",
    "van Gorp ECM, Mairuhu ATA, Beijnen JH. Development and validation of",
    "a population pharmacokinetic model for ritonavir used as a booster",
    "or as an antiviral agent in HIV-1-infected patients.",
    "Br J Clin Pharmacol. 2005;59(2):174-182.",
    "doi:10.1111/j.1365-2125.2004.02241.x."
  )
  vignette <- "Kappelhoff_2005_ritonavir"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CONMED_LPV = list(
      description        = paste(
        "1 = subject is receiving concomitant lopinavir as part of an",
        "antiretroviral regimen at the time of the observation;",
        "0 = subject is not on lopinavir."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant lopinavir)",
      notes              = paste(
        "Effect on apparent oral clearance is a multiplicative power",
        "form, CL/F = exp(lcl) * e_lpv_cl^CONMED_LPV, with",
        "e_lpv_cl = 2.72 (Kappelhoff 2005 Results equation following",
        "Table 2). 36 of the 186 subjects were on lopinavir/ritonavir",
        "in the source cohort (Table 1 'Lopinavir/Ritonavir (booster",
        "(n,%))' = 36 (19.4%)). Indinavir and saquinavir co-",
        "administration were tested and not retained, so CONMED_LPV",
        "is the only co-medication covariate in the final model.",
        "Time-fixed within an evaluated regimen in the source cohort;",
        "if a regimen switch occurs in simulation use it can be",
        "supplied as a time-varying covariate column."
      ),
      source_name        = "LPV"
    ),
    MIX_LARGE_RUV = list(
      description        = paste(
        "Latent mixture-model class indicator for the additive residual",
        "error magnitude: 1 = subject classified to the minority larger",
        "additive-RUV subpopulation P2 (addSd = 0.199 mg/L); 0 = subject",
        "classified to the majority smaller additive-RUV subpopulation",
        "P1 (addSd = 0.0600 mg/L). Both subpopulations share the same",
        "15.4% proportional RUV component."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (smaller-additive-RUV subpopulation P1)",
      notes              = paste(
        "Not a measured patient covariate -- this is the per-subject",
        "posterior class-membership index from Kappelhoff 2005's NONMEM",
        "$MIX block (paper Methods, 'Basic pharmacokinetic model':",
        "'Subpopulations were estimated using the $MIX function in the",
        "control stream'). Population probability of MIX_LARGE_RUV = 0",
        "is 0.648 (Table 2 'Fraction in P1 (%)' = 64.8); probability of",
        "MIX_LARGE_RUV = 1 is 0.352. For typical-value simulation set",
        "MIX_LARGE_RUV = 0 (dominant subpopulation). For population",
        "simulation, draw MIX_LARGE_RUV ~ Bernoulli(0.352) per subject.",
        "The paper (Discussion) notes the mechanistic explanation could",
        "not be identified."
      ),
      source_name        = "$MIX class assignment"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 186L,
    n_studies      = 1L,
    n_observations = 1228L,
    age_range      = "median 39.4 years; IQR 35.0-46.0 (Table 1)",
    weight_range   = "median 71.5 kg; IQR 63.0-79.8 (Table 1)",
    sex_female_pct = 12.4,
    race_ethnicity = c(
      White         = 64.5,
      Black         = 10.2,
      Asian         = 4.3,
      Hispanic      = 5.4,
      Missing       = 15.6
    ),
    disease_state  = paste(
      "Ambulatory HIV-1-infected adults from the outpatient clinic at",
      "Slotervaart Hospital, Amsterdam, sampled between January 1999",
      "and June 2003. Baseline median CD4 240 cells/mm^3 (IQR 110-380),",
      "baseline median log10 HIV-1 RNA 4.86 copies/mL (IQR 3.48-5.47).",
      "Hepatitis B co-infection 3.8% (7/186); hepatitis C co-infection",
      "9.1% (17/186). Gender, race, ALAT, ASAT, GGT, AP, TBR, CD4, CD8,",
      "viral load, age, weight, serum creatinine, HBV / HCV status,",
      "and concomitant saquinavir, indinavir or lopinavir use were all",
      "tested in covariate screening; only concomitant lopinavir was",
      "retained in the final model."
    ),
    dose_range     = paste(
      "Oral ritonavir 100-750 mg once or twice daily. 115 subjects",
      "(62%) received ritonavir 100 mg QD, 100 mg BID, 133 mg BID, or",
      "200 mg BID as a booster of another protease inhibitor (49 with",
      "indinavir, 78 with saquinavir, 36 with lopinavir); 71 subjects",
      "(38%) received ritonavir as an antiviral agent at 300, 400,",
      "500, 600, or 750 mg twice daily."
    ),
    regions        = "Netherlands (Slotervaart Hospital, Amsterdam)",
    notes          = paste(
      "Sparse per-visit single-time-point samples (505) supplemented",
      "with 55 full PK profiles (12-15 sampling time points each);",
      "average 3-4 samples per patient over 7-12 months follow-up",
      "(up to 28 months in some patients). All samples collected at",
      "steady state, at least 2 weeks after initiation of a ritonavir-",
      "containing regimen. Plasma ritonavir measured by isocratic",
      "reversed-phase ion-pair HPLC with UV detection (assay range",
      "0.05-25 mg/L). Plasma concentrations < 0.01 mg/L were excluded",
      "as likely non-adherence."
    )
  )

  ini({
    # ============================================================
    # Structural PK -- one-compartment, first-order absorption with
    # a lag time, first-order elimination. All point estimates from
    # Kappelhoff 2005 Table 2 (final model, 'Estimate' column).
    # ============================================================
    lcl    <- log(10.5);   label("Apparent oral clearance without LPV, CL/F (L/h)")  # Table 2 row "CL/F (L h-1)" = 10.5, RSE 5.55%
    lvc    <- log(96.6);   label("Apparent volume of distribution, V/F (L)")          # Table 2 row "V/F (L)" = 96.6, RSE 10.7%
    lka    <- log(0.871);  label("First-order absorption rate, ka (1/h)")             # Table 2 row "ka (h-1)" = 0.871, RSE 23.1%
    ltlag  <- log(0.778);  label("Absorption lag time, Tlag (h)")                     # Table 2 row "Lag-time (h)" = 0.778, RSE 4.91%

    # Concomitant-lopinavir effect on apparent oral clearance.
    # Kappelhoff 2005 Results equation (after Table 2):
    #   CL/F = 10.5 * 2.72^LPV
    # where LPV = 1 if subject is on lopinavir, 0 otherwise. Encoded
    # as the multiplicative power form e_lpv_cl^CONMED_LPV inside
    # model() so that exp(lcl) is preserved as the LPV = 0 typical
    # value.
    e_lpv_cl <- 2.72;      label("Multiplicative factor on CL/F when CONMED_LPV = 1 (unitless)")  # Table 2 row "theta_LPV" = 2.72, RSE 11.7%

    # ============================================================
    # Inter-individual variability. Source reports CV%; convert to
    # log-normal variance via omega^2 = log(1 + CV^2). The V/F and
    # ka random effects share a correlated block; the off-diagonal
    # covariance is rho * sqrt(var_vc * var_ka) with rho = 0.868.
    # ============================================================
    etalcl ~ log(1 + 0.383^2)
    # Table 2 row "IIV CL/F (%)" = 38.3; omega^2 = log(1 + 0.383^2) approx 0.137

    etalvc + etalka ~ c(
      log(1 + 0.800^2),
      0.868 * sqrt(log(1 + 0.800^2) * log(1 + 1.690^2)),
      log(1 + 1.690^2)
    )
    # Table 2 rows "IIV V/F (%)" = 80.0, "IIV ka (%)" = 169,
    # "Correlation eta_V - eta_ka" = 0.868.

    # ============================================================
    # Residual unexplained variability. The source uses a mixture
    # model on the ADDITIVE component: subpopulation P1 (64.8% of
    # subjects) has additive SD = 0.0600 mg/L; subpopulation P2
    # (35.2%) has additive SD = 0.199 mg/L. Both subpopulations
    # share the same 15.4% proportional component. The per-subject
    # mixture class is supplied via the binary covariate
    # MIX_LARGE_RUV (0 = P1, 1 = P2) and selects the active addSd
    # inside model().
    # ============================================================
    addSd_p1 <- 0.0600;    label("Additive residual SD, subpopulation P1 (mg/L)")  # Table 2 row "Additive error P1 (mg l-1)" = 0.0600, RSE 13.5%
    addSd_p2 <- 0.199;     label("Additive residual SD, subpopulation P2 (mg/L)")  # Table 2 row "Additive error P2 (mg l-1)" = 0.199, RSE 15.2%
    propSd   <- 0.154;     label("Proportional residual SD, both subpopulations (fraction)")  # Table 2 row "Proportional error (%)" = 15.4, RSE 23.8%
  })

  model({
    # Individual PK parameters. Lopinavir-driven increase in apparent
    # CL/F follows the multiplicative power form printed in the
    # Results section: CL/F_i = exp(lcl + etalcl) * 2.72^CONMED_LPV.
    cl   <- exp(lcl + etalcl) * e_lpv_cl^CONMED_LPV
    vc   <- exp(lvc + etalvc)
    ka   <- exp(lka + etalka)
    tlag <- exp(ltlag)

    # Micro-constant (first-order elimination)
    kel  <- cl / vc

    # ODE system: oral depot + central compartment with first-order
    # absorption (with lag time) and first-order elimination.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central
    alag(depot)   <- tlag

    # Plasma concentration of ritonavir (mg/L)
    Cc <- central / vc

    # Mixture-gated additive residual error magnitude. The proportional
    # component is shared across both subpopulations; only the additive
    # SD switches between P1 (0.0600 mg/L) and P2 (0.199 mg/L) based on
    # the per-subject indicator MIX_LARGE_RUV. The pattern matches the
    # Allegaert 2015 paracetamol OCC = 5 piecewise-RUV implementation.
    add_sd_eff <- addSd_p1 * (1 - MIX_LARGE_RUV) + addSd_p2 * MIX_LARGE_RUV
    Cc ~ prop(propSd) + add(add_sd_eff)
  })
}
