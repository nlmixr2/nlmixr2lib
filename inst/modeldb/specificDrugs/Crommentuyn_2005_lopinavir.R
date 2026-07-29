Crommentuyn_2005_lopinavir <- function() {
  description <- paste(
    "One-compartment first-order-absorption population PK model for oral",
    "lopinavir co-administered with ritonavir in 122 HIV-1-infected adults",
    "on BID lopinavir/ritonavir 400-666/100-166 mg. Apparent oral clearance",
    "CL/F follows an inverse-saturable function of per-subject ritonavir",
    "AUC over the 12 h dosing interval (CONMED_RTV_AUC_12h, mg*h/L, computed",
    "from the upstream Kappelhoff 2005 ritonavir popPK model) plus a pooled",
    "+39% NNRTI co-medication factor (efavirenz or nevirapine, encoded as",
    "the CONMED_NNRTI class indicator). IIV is estimated on ka, CL/F, and",
    "V/F as a full 3x3 correlated block; residual error is combined",
    "additive plus proportional. The reported IOV on relative",
    "bioavailability F (17.5% CV) is NOT encoded structurally (Brooks 2021",
    "precedent); downstream users who want IOV can add an OCC covariate and",
    "a per-occasion eta in rxode2 (Crommentuyn 2005)."
  )
  reference <- paste(
    "Crommentuyn KML, Kappelhoff BS, Mulder JW, Mairuhu ATA, van Gorp ECM,",
    "Meenhorst PL, Huitema ADR, Beijnen JH. Population pharmacokinetics of",
    "lopinavir in combination with ritonavir in HIV-1-infected patients.",
    "Br J Clin Pharmacol. 2005 Oct;60(4):378-389.",
    "doi:10.1111/j.1365-2125.2005.02455.x.",
    "Per-subject ritonavir AUC over the 12 h dosing interval is computed",
    "from the upstream Kappelhoff et al. 2005 ritonavir popPK model",
    "(Br J Clin Pharmacol 2005;59:174-82) and supplied as the time-fixed",
    "CONMED_RTV_AUC_12h covariate; the upstream ritonavir model is not",
    "structurally re-instantiated here (consistent with the Dickinson 2009",
    "atazanavir precedent for an AUC-of-ritonavir-as-covariate encoding)."
  )
  vignette <- "Crommentuyn_2005_lopinavir"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CONMED_RTV_AUC_12h = list(
      description        = "Per-subject ritonavir AUC over the 12 h dosing interval (BID ritonavir)",
      units              = "mg*h/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-subject (time-fixed within an evaluated regimen) ritonavir",
        "AUC over the 12 h dosing interval. In Crommentuyn 2005 the value",
        "is computed per-subject as DOSE_RTV / CL_RTV using individual",
        "Bayesian CL_RTV estimates from the upstream Kappelhoff 2005",
        "ritonavir popPK model (paper reference 19). Enters lopinavir",
        "CL/F via the inverse-saturable form",
        "cl = exp(lcl) * (auc50 / (auc50 + CONMED_RTV_AUC_12h)) * IND",
        "with auc50 = 2.26 mg*h/L (Crommentuyn 2005 Methods Equation 2,",
        "Results page 6, Table 2 row AUC50). Crommentuyn 2005 cohort",
        "median 3.58 mg*h/L (range 0.85-18.77); at that median, lopinavir",
        "CL/F = 14.8 * 2.26/(2.26+3.58) = 5.73 L/h (matches the value",
        "reported on Results page 6). For simulation users without",
        "observed ritonavir AUC, the cohort median 3.58 reproduces",
        "typical-value behaviour."
      ),
      source_name        = "AUC12h"
    ),
    CONMED_NNRTI = list(
      description        = "Concomitant CYP3A4-inducing NNRTI (efavirenz or nevirapine) indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant NNRTI; typically NRTI-only nucleoside backbone)",
      notes              = paste(
        "Time-varying per occasion in Crommentuyn 2005 (8% of patients",
        "on efavirenz across 59 samples, 13% on nevirapine across 73",
        "samples; per-occasion flagging). Pooled class indicator -- the",
        "paper tested separate factors for efavirenz (1.41) and nevirapine",
        "(1.52) and found no improvement in fit over a single pooled NNRTI",
        "factor (1.39); see Results page 7. CONMED_NNRTI = 1 = subject on",
        "EFV or NVP at the observation occasion; 0 = on the NRTI-only",
        "nucleoside backbone without an NNRTI. Effect on apparent lopinavir",
        "CL/F is multiplicative: cl *= (1 + 0.39 * CONMED_NNRTI)",
        "(Crommentuyn 2005 Table 2 row IND = 1.39 and Results page 7)."
      ),
      source_name        = "IND"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 122L,
    n_studies      = 1L,
    n_observations = 748L,
    age_range      = "median 42 years (IQR 36-46)",
    weight_range   = "median 72 kg (IQR 63-80)",
    sex_female_pct = 15,
    race_ethnicity = c(Caucasian = 75, Black = 12, Asian = 9, Latino = 4),
    disease_state  = "HIV-1 infection on a lopinavir/ritonavir-containing antiretroviral regimen",
    dose_range     = paste(
      "Oral co-formulated lopinavir/ritonavir capsules (133/33 mg each),",
      "400/100 mg to 666/166 mg twice daily (3-5 capsules BID)."
    ),
    regions        = "Netherlands (Slotervaart Hospital outpatient clinic, Amsterdam)",
    iov_structure  = paste(
      "Crommentuyn 2005 Table 2 reports IOV on relative bioavailability F",
      "of 17.5% CV. This model file does NOT encode IOV structurally --",
      "the model-library convention is to omit IOV when no occasion column",
      "is defined for downstream simulation (Andrews 2017 / Brooks 2021",
      "precedent); see vignette Assumptions and deviations."
    ),
    notes          = paste(
      "Retrospective TDM analysis (February 2001 - March 2004). 748",
      "lopinavir + 748 ritonavir plasma concentrations: 14 full",
      "pharmacokinetic profiles plus 568 single-time-point random samples;",
      "median 4 samples per patient (range 1-14). Sampling at steady state",
      "(minimum 2 weeks on treatment). 23% of patients (28 of 122) on",
      "concomitant tenofovir; 60% pre-treated with protease inhibitors;",
      "11% chronic HCV; 7% chronic HBV. Lopinavir / ritonavir quantified",
      "by HPLC-MS/MS (LoQ lopinavir 0.1 mg/L, ritonavir 0.05 mg/L).",
      "Concentrations below 0.7 mg/L lopinavir flagged for compliance",
      "review; 11 patients with confirmed non-compliance or unknown",
      "post-dose time were excluded before model fitting."
    )
  )

  # Documented-but-not-retained covariates (Crommentuyn 2005 covariate analysis,
  # Results page 7 / Discussion): screened in the stepwise build, did NOT meet
  # the dOFV/clinical-relevance retention criteria, and are NOT referenced in
  # model(). Recorded here for provenance per the extract-literature-model
  # skill's covariatesDataExcluded convention so they do not trigger the
  # "declared but unused" convention warning.
  covariatesDataExcluded <- list(
    WT       = list(description = "Body weight", units = "kg", type = "continuous",
                    notes = "Median 72 kg (IQR 63-80) per Table 1; screened but not retained on CL/F or V/F."),
    AGE      = list(description = "Subject age", units = "years", type = "continuous",
                    notes = "Median 42 years (IQR 36-46) per Table 1; screened but not retained."),
    SEXF     = list(description = "Sex (1 = female, 0 = male)", units = "(binary)", type = "binary",
                    notes = "15% female (18 of 122) per Table 1; screened but not retained."),
    RACE_BLACK    = list(description = "Black-race indicator", units = "(binary)", type = "binary",
                         notes = "12% (14 of 122) per Table 1; screened but not retained."),
    RACE_ASIAN    = list(description = "Asian-race indicator", units = "(binary)", type = "binary",
                         notes = "9% (11 of 122) per Table 1; screened but not retained."),
    RACE_LATINO   = list(description = "Latino-race indicator", units = "(binary)", type = "binary",
                         notes = "4% (5 of 122) per Table 1; screened but not retained."),
    ALAT     = list(description = "Alanine aminotransferase", units = "U/L", type = "continuous",
                    notes = "Median 35 U/L (IQR 22-49) per Table 1; screened but not retained."),
    ASAT     = list(description = "Aspartate aminotransferase", units = "U/L", type = "continuous",
                    notes = "Median 30 U/L (IQR 22-43) per Table 1; screened but not retained."),
    AP       = list(description = "Alkaline phosphatase", units = "U/L", type = "continuous",
                    notes = "Median 80 U/L (IQR 63-94) per Table 1; screened but not retained."),
    TBR      = list(description = "Serum total bilirubin", units = "umol/L", type = "continuous",
                    notes = "Median 11 umol/L (IQR 8-14) per Table 1; screened but not retained."),
    HCV      = list(description = "Chronic hepatitis C indicator", units = "(binary)", type = "binary",
                    notes = "11% (13 of 122) per Table 1; screened but not retained."),
    HBV      = list(description = "Chronic hepatitis B indicator", units = "(binary)", type = "binary",
                    notes = "7% (8 of 122) per Table 1; screened but not retained."),
    CONMED_TENOFOVIR = list(description = "Concomitant tenofovir indicator", units = "(binary)", type = "binary",
                            notes = "23% (28 of 122) per Table 1; tested on CL/F and V/F per Results page 7 and not retained.")
  )

  ini({
    # ============================================================
    # Structural PK parameters -- Crommentuyn 2005 Table 2
    # (final-model "Est" column).
    # The CL/F value (14.8 L/h) is the typical value of apparent
    # oral clearance in the ABSENCE of ritonavir, per the paper's
    # final-model equation (Results page 7); the inverse-saturable
    # term in model() reduces this to the observed CL/F at non-zero
    # ritonavir AUC (5.73 L/h at the cohort-median AUC of 3.58
    # mg*h/L, matching Results page 6).
    # ============================================================
    lka     <- log(0.564); label("First-order absorption rate ka (1/h)")                                    # Table 2 row 1: ka = 0.564 /h
    lcl     <- log(14.8);  label("Apparent clearance CL/F in the absence of ritonavir (L/h)")               # Table 2 row 2: CL/F = 14.8 L/h
    lvc     <- log(61.6);  label("Apparent volume of distribution V/F (L)")                                 # Table 2 row 5: V/F = 61.6 L
    lauc50  <- log(2.26);  label("Ritonavir AUC12h producing 50% reduction of CL/F (mg*h/L)")               # Table 2 row 3: AUC50 = 2.26 mg*h/L

    # ============================================================
    # Covariate effects.
    # NNRTI co-medication: pooled efavirenz / nevirapine class
    # factor; IND = 1.39 in the paper means CL/F is +39% during
    # treatment with EFV or NVP. Encoded as the additive fraction
    # so the same parameter name pattern as other CONMED_* effects
    # carries over: IND = 1 + e_nnrti_cl * CONMED_NNRTI.
    # ============================================================
    e_nnrti_cl <- 0.39;   label("Fractional increase in CL/F per NNRTI co-medication (unitless)")           # Table 2 row 4: IND = 1.39 (= 1 + 0.39) and Results page 7

    # ============================================================
    # Inter-individual variability -- correlated full 3x3 block on
    # ka, CL/F, V/F. The paper reports IIV as CV%; convert to the
    # log-normal variance scale via omega^2 = log(1 + CV^2). The
    # off-diagonal covariances use cov = rho * sqrt(var_X * var_Y)
    # with the paper-reported correlations from Table 2.
    #
    #   Table 2:  IIV ka = 97.8% CV          --> var_ka = log(1 + 0.978^2)
    #             IIV CL = 17.2% CV          --> var_cl = log(1 + 0.172^2)
    #             IIV V  = 63.8% CV          --> var_vc = log(1 + 0.638^2)
    #             rho(ka, CL) = -0.267
    #             rho(ka, V)  =  0.822
    #             rho(CL, V)  =  0.242
    # ============================================================
    etalka + etalcl + etalvc ~ c(
      log(1 + 0.978^2),
      -0.267 * sqrt(log(1 + 0.978^2) * log(1 + 0.172^2)), log(1 + 0.172^2),
       0.822 * sqrt(log(1 + 0.978^2) * log(1 + 0.638^2)),  0.242 * sqrt(log(1 + 0.172^2) * log(1 + 0.638^2)), log(1 + 0.638^2)
    )
    # Note: IOV on F (17.5% CV per Table 2 row 9) is NOT encoded
    # structurally -- see vignette Assumptions and deviations.

    # ============================================================
    # Residual error -- combined additive + proportional. Paper
    # text (Results page 6 and Discussion page 9) reports the
    # additive 1.15 mg/L component was deliberately retained to
    # down-weight low / suspected-non-compliance concentrations.
    # ============================================================
    addSd   <- 1.15;       label("Additive residual error (mg/L)")                                          # Table 2 row 14: Additive error = 1.15 mg/L
    propSd  <- 0.0755;     label("Proportional residual error (fraction)")                                  # Table 2 row 15: Proportional error = 7.55%
  })

  model({
    # ------------------------------------------------------------
    # NNRTI co-medication multiplicative factor on CL/F.
    # CONMED_NNRTI = 1 when subject is on EFV or NVP; = 0
    # otherwise. With e_nnrti_cl = 0.39 the on-NNRTI factor is 1.39
    # (a +39% multiplicative increase) per Crommentuyn 2005 Results
    # page 7 / Table 2 row IND.
    # ------------------------------------------------------------
    ind <- 1 + e_nnrti_cl * CONMED_NNRTI

    # ------------------------------------------------------------
    # Inverse-saturable ritonavir-AUC dependence on lopinavir CL/F:
    #   CL/F = exp(lcl + etalcl) * AUC50 / (AUC50 + AUC_RTV) * IND
    # (Crommentuyn 2005 Methods Equation 2 / Results page 7). At
    # AUC_RTV = 0 the saturation term evaluates to 1 and CL/F
    # reduces to the no-ritonavir typical value exp(lcl) = 14.8
    # L/h; at the cohort-median AUC_RTV = 3.58 mg*h/L the typical
    # CL/F is 14.8 * 2.26 / (2.26 + 3.58) = 5.73 L/h, matching
    # the value the paper reports.
    # ------------------------------------------------------------
    auc50 <- exp(lauc50)
    cl    <- exp(lcl + etalcl) * (auc50 / (auc50 + CONMED_RTV_AUC_12h)) * ind
    vc    <- exp(lvc + etalvc)
    ka    <- exp(lka + etalka)

    # Micro-constant
    kel   <- cl / vc

    # ------------------------------------------------------------
    # ODE system: 1-compartment with first-order absorption from a
    # depot compartment. Crommentuyn 2005 Results page 6 reports
    # that lag-time / zero-order absorption forms did not improve
    # fit, so absorption is plain first-order without lag.
    # ------------------------------------------------------------
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ------------------------------------------------------------
    # Observation and combined additive + proportional residual
    # error.
    # ------------------------------------------------------------
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
