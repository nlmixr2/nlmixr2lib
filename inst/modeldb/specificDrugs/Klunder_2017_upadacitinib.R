Klunder_2017_upadacitinib <- function() {
  description <- "Two-compartment population PK model with first-order absorption and an absorption lag time for oral upadacitinib (ABT-494), a selective JAK1 inhibitor, in healthy adults and adults with rheumatoid arthritis (Klunder 2017, pooled phase I + phase IIb analysis). Statistically significant covariates retained in the final model: population (RA vs healthy) on CL/F, sex on CL/F and Vc/F, baseline creatinine clearance on CL/F (raw Cockcroft-Gault, not BSA-normalized), and total body weight on Vc/F. ISV is reported separately for healthy subjects and RA patients on CL/F and Vc/F, and is encoded here as paired healthy / RA structural means with cohort-specific log-normal random effects gated by DIS_HEALTHY."
  reference   <- "Klunder B, Mohamed M-EF, Othman AA. Population pharmacokinetics of upadacitinib in healthy subjects and subjects with rheumatoid arthritis: analyses of phase I and II clinical trials. Clin Pharmacokinet. 2018;57(8):977-988. doi:10.1007/s40262-017-0605-6"
  vignette    <- "Klunder_2017_upadacitinib"
  units       <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight at baseline.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on Vc/F (exponent 0.50). Reference body weight 74 kg, the male-population median in the analysis dataset (Klunder 2017 Table 3 footnote d). Body weight was tested but not retained as a covariate on CL/F.",
      source_name        = "WT"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male). Klunder 2017 Table 3 footnotes c and d use males as the reference category for both the CL/F (theta2 = 1 for males) and Vc/F (theta4 = 1 for males) sex contrasts.",
      notes              = "Multiplicative effects on CL/F (theta2 = 0.86; females have 14% lower CL/F than males) and on Vc/F (theta4 = 0.75; females have 25% lower Vc/F than males) per Klunder 2017 Table 3. Encoded canonically (SEXF = 1 for females; same orientation as the paper).",
      source_name        = "SEX"
    ),
    CRCL = list(
      description        = "Baseline creatinine clearance estimated by the Cockcroft-Gault formula (raw mL/min, NOT BSA-normalized).",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on CL/F (exponent 0.32). Reference CrCL 107 mL/min, the healthy-male population median in the analysis dataset (Klunder 2017 Table 3 footnote c). Klunder 2017 reports raw Cockcroft-Gault CrCL (mL/min), distinct from the BSA-normalized eGFR variants of CRCL used in some other registered models; the analysis range was 41-241 mL/min. Serum creatinine itself was tested but not retained -- only the Cockcroft-Gault-derived CrCL was statistically significant on CL/F.",
      source_name        = "CrCL"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant cohort indicator (1 = healthy adult, 0 = adult with rheumatoid arthritis).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (rheumatoid arthritis patient) under the canonical orientation. Klunder 2017 Table 3 footnote c reports the contrast in the paper's PATIENT-coded direction (theta1 = 0.76 = CL_RA / CL_healthy, applied as theta1^RA with theta1 = 1 for healthy volunteers); the canonical DIS_HEALTHY register uses DIS_HEALTHY = 1 for healthy participants. Because every eta<X> must pair with a structural lX in ini(), the disease-state contrast is encoded via paired healthy / RA structural means (lcl_h / lcl_ra and lvc_h / lvc_ra; the RA means are reconstructed as log(39.7 * 0.76) and log(146) to preserve Klunder 2017 Table 3 verbatim). DIS_HEALTHY then gates which mean and which log-normal random effect applies for each subject.",
      notes              = "Time-fixed per subject. RA patients had 24% lower upadacitinib CL/F than matched healthy subjects (corresponding to ~32% higher steady-state AUC), attributed by the authors to a composite of older age, possibly elevated IL-6 / inflammation suppressing CYP3A activity, and lower metabolic capacity in the RA cohort (Klunder 2017 Discussion). The cohort split also drives the heteroscedastic ISV: healthy %ISV(CL/F) = 16% vs RA %ISV(CL/F) = 26%; healthy %ISV(Vc/F) = 14% vs RA %ISV(Vc/F) = 27%.",
      source_name        = "POP (1 = RA in the paper; DIS_HEALTHY = 1 - POP)"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 573,
    n_studies       = 5,
    age_range       = "19-85 years (pooled); healthy 19-56, RA 19-85",
    age_median      = "52 years (pooled cohort, Klunder 2017 Table 2)",
    weight_range    = "42-134 kg (pooled cohort)",
    weight_median   = "76 kg (pooled cohort, Klunder 2017 Table 2; reference for Vc/F covariate is 74 kg, the male-population median per Table 3 footnote d)",
    sex_female_pct  = 65,
    race_ethnicity  = c(White = 83, Black = 7, Asian = 6, Other = 3),
    disease_state   = "Pooled cohort: 107 healthy subjects and 466 adult patients with rheumatoid arthritis (RA). RA cohort recruited under MTX-IR (study 4) or anti-TNF-IR (study 5) inadequate-response definitions; all RA subjects had >=3 months of active disease.",
    dose_range      = "Upadacitinib immediate-release capsules. Phase I: single doses 1-48 mg and twice-daily multiple doses 3-24 mg for up to 14 days (study 2) or 26 days (study 2 RA cohort). Phase IIb: 3, 6, 12, 18 mg BID and 24 mg QD for 12 weeks (study 4 MTX-IR); 3, 6, 12, 18 mg BID for 12 weeks (study 5 anti-TNF-IR).",
    regions         = "Multinational (USA-based phase I plus phase IIb at international RA sites). Study 3 enrolled Japanese and Chinese subjects.",
    renal_function  = "Mild and moderate renal impairment (CrCL 41-241 mL/min) included; subjects with eGFR < 40 mL/min/1.73 m^2 excluded by protocol.",
    n_observations  = "6,399 quantifiable plasma concentrations. ~4% of post-baseline samples were below the LLOQ; intra-dose BLQ values and the first post-last-dose BLQ were imputed at LLOQ/2 and subsequent BLQ values were censored. ~0.9% of records were excluded as outliers consistent with inaccurate dosing records.",
    co_medication   = "Concomitant CYP3A inhibitors (moderate/weak), CYP3A inducers, CYP2D6 inhibitors, and pH-modifying medications were tested as covariates but none were retained. Strong CYP3A inhibitors and inducers were prohibited in studies 4 and 5.",
    notes           = "Clinical-trial registry identifiers: NCT01741493, NCT02066389, NCT01960855. NONMEM 7.3 with FOCEI for the structural / covariate model and SAEM + importance sampling for the final model incorporating box-cox-transformed IOV on ALAG1 in phase II visits."
  )

  ini({
    # ---------------- Structural parameters ----------------------------------
    # Reference covariate values for the typical healthy male (Klunder 2017
    # Table 3 footnotes c and d): SEXF = 0 (male), DIS_HEALTHY = 1 (healthy),
    # WT = 74 kg, CrCL = 107 mL/min. The disease-state contrast is encoded via
    # paired healthy / RA structural means (lcl_h / lcl_ra and lvc_h / lvc_ra)
    # so that every eta<X> in this ini() block has a corresponding fixed
    # effect <X> for convention compliance. The RA means are reconstructed
    # from the paper's Table 3 healthy reference and theta1 ratio so the
    # values trace cleanly to the source: lcl_ra = log(39.7 * 0.76) =
    # log(30.172), and lvc_ra = log(146) because Klunder 2017 reports no
    # disease-state effect on the typical Vc/F (only on the Vc/F ISV variance).
    lka     <- log(12.3);         label("First-order absorption rate constant Ka (1/h)")                      # Klunder 2017 Table 3: Ka = 12.3 1/h, modeled as exp(theta) with theta = 2.51 (RSE 5.3%)
    lcl_h   <- log(39.7);         label("Apparent oral clearance CL/F in healthy males at CrCL=107 (L/h)")    # Klunder 2017 Table 3: CL/F = 39.7 L/h, RSE 2% (healthy / male / CrCL=107 reference)
    lcl_ra  <- log(39.7 * 0.76);  label("Apparent oral clearance CL/F in RA males at CrCL=107 (L/h)")         # Klunder 2017 Table 3: CL/F * theta1 = 39.7 * 0.76 = 30.172 L/h (theta1 = 0.76, RSE 2.9%)
    lvc_h   <- log(146);          label("Apparent central volume Vc/F in healthy males at WT=74 kg (L)")      # Klunder 2017 Table 3: Vc/F = 146 L, RSE 2.2%
    lvc_ra  <- log(146);          label("Apparent central volume Vc/F in RA males at WT=74 kg (L)")           # Klunder 2017 Table 3: no disease-state effect on typical Vc/F; the cohort split affects only the Vc/F ISV variance (etalvc_h vs etalvc_ra below)
    lvp     <- log(64.3);         label("Apparent peripheral volume Vp/F (L)")                                # Klunder 2017 Table 3: Vp/F = 64.3 L, RSE 1.1%
    lq      <- log(3.23);         label("Apparent intercompartmental clearance Q/F (L/h)")                    # Klunder 2017 Table 3: Q/F = 3.23 L/h, RSE 2.5%
    ltlag   <- log(0.48);         label("Absorption lag time ALAG1 (h)")                                      # Klunder 2017 Table 3: ALAG1 = 0.48 h, RSE 0.017%

    # ---------------- Covariate effects --------------------------------------
    # Klunder 2017 Table 3 footnotes c and d:
    #   CL/F_typ = CL/F * theta1^RA * theta2^female * (CrCL/107)^theta3
    #   Vc/F_typ = Vc/F * theta4^female * (WT/74)^theta5
    # The theta1 (RA vs healthy) contrast is absorbed into lcl_ra above; the
    # remaining effects act identically on the healthy and RA cohorts.
    e_sexf_cl  <- log(0.86); label("Log ratio of CL/F (female vs male)")                                      # Klunder 2017 Table 3: theta2 = 0.86, RSE 2.8% (CL/F female vs male)
    e_crcl_cl  <- 0.32;      label("Power exponent: CrCL on CL/F (unitless)")                                 # Klunder 2017 Table 3: theta3 = 0.32, RSE 11% (CrCL on CL/F)
    e_sexf_vc  <- log(0.75); label("Log ratio of Vc/F (female vs male)")                                      # Klunder 2017 Table 3: theta4 = 0.75, RSE 1.6% (Vc/F female vs male)
    e_wt_vc    <- 0.50;      label("Power exponent: WT on Vc/F (unitless)")                                   # Klunder 2017 Table 3: theta5 = 0.50, RSE 15.8% (body weight on Vc/F)

    # ---------------- IIV (variances on log scale) ---------------------------
    # Klunder 2017 reports %ISV in Table 3 with the footnote
    # "%ISV was calculated as SQRT(omega^2) * 100", i.e. variance = (%ISV/100)^2
    # on the log scale (the standard NONMEM convention reported as %CV for
    # log-normal random effects). ISV on CL/F and Vc/F is split by cohort:
    # the healthy etas pair with lcl_h / lvc_h and are gated by DIS_HEALTHY,
    # and the RA etas pair with lcl_ra / lvc_ra and are gated by
    # (1 - DIS_HEALTHY) in model().
    etalka    ~ 2.25                                                                                          # Klunder 2017 Table 3: %ISV Ka = 150%, RSE 39%; 1.50^2
    etalcl_h  ~ 0.0256                                                                                        # Klunder 2017 Table 3: %ISV CL/F healthy = 16%, RSE 8.5%; 0.16^2
    etalcl_ra ~ 0.0676                                                                                        # Klunder 2017 Table 3: %ISV CL/F RA = 26%, RSE 2.9%; 0.26^2
    etalvc_h  ~ 0.0196                                                                                        # Klunder 2017 Table 3: %ISV Vc/F healthy = 14%, RSE 14%; 0.14^2
    etalvc_ra ~ 0.0729                                                                                        # Klunder 2017 Table 3: %ISV Vc/F RA = 27%, RSE 5.5%; 0.27^2

    # ---------------- Residual error (combined additive + proportional) ------
    # Klunder 2017 Methods Eq. 2-3 specify a combined additive + proportional
    # residual error model: C_obs = C_pred * (1 + eps_prop) + eps_add.
    addSd  <- 0.18; label("Additive residual SD (ng/mL)")                                                     # Klunder 2017 Table 3: additive residual error SD = 0.18 ng/mL, RSE 13%
    propSd <- 0.31; label("Proportional residual SD (fraction)")                                              # Klunder 2017 Table 3: proportional residual error SD = 31%, RSE 13%
  })

  model({
    # Reference covariate values for the typical healthy male (Klunder 2017
    # Table 3 footnotes c and d).
    ref_wt   <- 74    # kg, male-population median (Table 3 footnote d)
    ref_crcl <- 107   # mL/min, healthy-male population median (Table 3 footnote c)

    # ---------------- Individual parameters ----------------------------------
    # Disease-state-specific structural means and IIV are selected by the
    # DIS_HEALTHY indicator. Each subject contributes exactly one of the two
    # log-normal random effects per parameter (the other is multiplied by 0),
    # so the per-subject log-normal variance on CL/F is omega_h^2 = 0.0256
    # for healthy subjects and omega_ra^2 = 0.0676 for RA patients (mirror
    # for Vc/F: 0.0196 vs 0.0729).
    lcl_eff    <- lcl_h * DIS_HEALTHY + lcl_ra * (1 - DIS_HEALTHY)
    etalcl_eff <- etalcl_h * DIS_HEALTHY + etalcl_ra * (1 - DIS_HEALTHY)
    lvc_eff    <- lvc_h * DIS_HEALTHY + lvc_ra * (1 - DIS_HEALTHY)
    etalvc_eff <- etalvc_h * DIS_HEALTHY + etalvc_ra * (1 - DIS_HEALTHY)

    ka <- exp(lka + etalka)
    cl <- exp(lcl_eff + etalcl_eff
              + e_sexf_cl * SEXF
              + e_crcl_cl * log(CRCL / ref_crcl))
    vc <- exp(lvc_eff + etalvc_eff
              + e_sexf_vc * SEXF
              + e_wt_vc * log(WT / ref_wt))
    vp   <- exp(lvp)
    q    <- exp(lq)
    tlag <- exp(ltlag)

    # Micro-constants for the explicit two-compartment ODE form.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---------------- ODE system ---------------------------------------------
    # Two-compartment disposition with first-order absorption from a depot
    # compartment and an absorption lag time applied at dose entry (Klunder
    # 2017 Methods 'Pharmacokinetic Model Development': two-compartment model
    # with ALAG1).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    # ---------------- Observation -------------------------------------------
    # Dose is in mg and vc in L, so central / vc yields mg/L = ug/mL. The
    # source paper reports plasma concentrations in ng/mL, so the model output
    # is rescaled by 1000 to match the source tables and figures.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
