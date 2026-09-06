# Two-compartment population PK of meltdose (LCP) tacrolimus (Envarsus) in
# elderly de novo kidney transplant recipients, with a three-compartment
# first-order absorption transit chain, a two-class absorption-lag mixture,
# inter-occasion variability on apparent clearance, and sample-matrix-specific
# proportional residual error for whole blood vs dried blood spot
# (Kamp 2023, Pharmaceutics 16(1):17; doi:10.3390/pharmaceutics16010017).

Kamp_2023_tacrolimus <- function() {
  description <- paste(
    "Two-compartment population PK model for meltdose (LCP, prolonged-release)",
    "tacrolimus in elderly de novo kidney transplant recipients. Oral",
    "absorption is a chain of three first-order compartments (the dosing depot",
    "plus two absorption transit compartments) sharing a single transit rate",
    "constant Ktr, feeding a two-compartment disposition with linear",
    "elimination. Apparent oral bioavailability is fixed at 1, so CL, Q, Vc and",
    "Vp are apparent (/F) values. The peripheral volume is fixed at 500 L per",
    "70 kg because it was not identifiable. Absorption lag time follows a",
    "two-class latent mixture gated by MIX_LAGGED_ABS: the lagged class (74% of",
    "the population) carries a 2.29 h lag and the no-lag class (26%) carries",
    "none. Allometric scaling on body",
    "weight (exponent 0.75 on CL and Q, 1 on Vc and Vp, reference 70 kg) was",
    "included a priori; no clinical covariate was retained in the final model.",
    "Inter-individual variability is log-normal on Ktr and on a correlated",
    "CL/Vc block; inter-occasion variability is carried on CL across five",
    "occasions. The proportional residual error switches per observation",
    "between whole blood (venous) and dried blood spot (capillary) sampling.",
    sep = " "
  )

  reference <- paste(
    "Kamp J, Zwart TC, Meziyerh S, van der Boog PJM, Nijgh EE, van Duin K,",
    "de Vries APJ, Moes DJAR. Meltdose Tacrolimus Population Pharmacokinetics",
    "and Limited Sampling Strategy Evaluation in Elderly Kidney Transplant",
    "Recipients. Pharmaceutics. 2024;16(1):17. doi:10.3390/pharmaceutics16010017.",
    "Published online 21 December 2023; PMCID PMC10819724. Structural model,",
    "final parameter estimates and the NONMEM control stream are taken from",
    "Table 2 and Supplementary Data S1.",
    sep = " "
  )

  vignette <- "Kamp_2023_tacrolimus"

  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling was included a priori on biological plausibility",
        "(Kamp 2023 Sect. 2.3.2): CL/F and Q/F scale as (WT/70)^0.75 and Vc/F",
        "and Vp/F as (WT/70)^1. Reference weight 70 kg. Cohort median 77.0 kg",
        "(IQR 68.1-84.5; Table 1). The a priori inclusion of allometric scaling",
        "improved the model by 9.7 OFV points (Sect. 3.3.2).",
        sep = " "
      ),
      source_name        = "WT"
    ),
    MIX_LAGGED_ABS = list(
      description        = paste(
        "Per-subject latent mixture-model class indicator for oral absorption",
        "lag (1 = lagged-absorption class, 0 = no-lag class)."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no-lag class)",
      notes              = paste(
        "Kamp 2023 fitted a NONMEM $MIX two-subpopulation model on the",
        "absorption lag time (Sect. 3.3.1; Supplementary Data S1). The lagged",
        "class carries the estimated 2.29 h lag and is the 74% MAJORITY; the",
        "no-lag class (lag structurally fixed at 0) is the 26% minority. So",
        "MIX_LAGGED_ABS = 1 with probability 0.74. The 'Pop parameter' of 26%",
        "in Table 2 is the fraction WITHOUT a lag, per the Table 2 and Table S2",
        "footnotes and per the $THETA comments of the control stream ('0 FIX ;",
        "Absorption lag time for group 1' with P(1) = THETA(11) = 0.26).",
        "CAUTION: the $PK IF-block as printed in Supplementary Data S1 reads",
        "'IF(MIXNUM.EQ.2) ALAG1 = THETA(4) ELSE ALAG1 = THETA(5)', which maps",
        "the 26% probability onto the 2.29 h lag and therefore contradicts both",
        "the footnotes and its own $THETA comments. The footnote reading is",
        "adopted because the paper's own Figure 4C -- the median of 1000",
        "simulated 7 mg q.d. profiles in a 70 kg subject from this final model",
        "-- has its median on the baseline through 3 h and peaks at about 6.1 h,",
        "which the 74%-lagged reading reproduces (simulated Tmax 6.25 h) and",
        "the 26%-lagged reading does not (simulated Tmax 4.25 h, with a median",
        "already at 5.4 ug/L by 2 h). See the vignette's Assumptions and",
        "deviations section for the full argument. Not a measured clinical",
        "covariate: the paper tested diabetes mellitus as a mechanistic",
        "correlate of the lag (diabetic gastroparesis hypothesis) and found no",
        "significant effect (Sect. 4). For typical-value simulation set",
        "MIX_LAGGED_ABS = 1 to reproduce the dominant 74% lagged phenotype; for",
        "population simulation draw MIX_LAGGED_ABS ~ Bernoulli(0.74) per",
        "subject.",
        sep = " "
      ),
      source_name        = "MIXNUM (NONMEM $MIX class index; MIX_LAGGED_ABS = as.integer(MIXNUM == 1))"
    ),
    OCC = list(
      description        = "Integer occasion index for the inter-occasion variability on apparent clearance",
      units              = "(count)",
      type               = "categorical",
      reference_category = "n/a - decomposed into mutually exclusive indicators inside model()",
      notes              = paste(
        "One occasion per AUC profile. Kamp 2023 Supplementary Data S1 declares",
        "five IOV slots via $ABBREVIATED REPLACE ETA(OCC_CL)=ETA(4,5,6,7,8) and",
        "$OMEGA BLOCK(1) 0.224 followed by four BLOCK(1) SAME blocks; Table 2",
        "reports five per-occasion shrinkage values (33; 13; 77; 100; 100),",
        "confirming five occasions. Each patient contributed a median of 2 AUCs",
        "(range 1-5; Table 1). Values outside 1-5 switch all IOV indicators off.",
        sep = " "
      ),
      source_name        = "OCC"
    ),
    SAMPLE_CAPILLARY = list(
      description        = "Per-observation blood sampling-matrix indicator (1 = dried blood spot / capillary finger-prick sample, 0 = venous whole blood sample)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (venous whole blood sample)",
      notes              = paste(
        "Record-level indicator that switches the proportional residual-error",
        "magnitude per observation: whole blood 20.8% and dried blood spot",
        "30.7% (Kamp 2023 Table 2). Supplementary Data S1 $ERROR selects",
        "PROPERR from THETA(1) or THETA(2) on the DBS0 data flag. Introducing",
        "the differential residual errors for the sample matrix improved the",
        "model by 8.2 OFV points (Sect. 3.3.1). Both matrices were quantified",
        "with the same LC-MS/MS method validated for whole blood and DBS",
        "(Sect. 2.2), so the contrast is the collection matrix rather than the",
        "assay. DBS samples were 42.5% of the AUCs in the dataset. May vary",
        "within a subject.",
        sep = " "
      ),
      source_name        = "DBS0"
    )
  )

  covariatesDataExcluded <- list(
    HCT = list(
      description = "Hematocrit",
      units       = "L/L",
      type        = "continuous",
      notes       = paste(
        "Screened on CL. Univariate analysis found a significant (p < 0.05)",
        "association (increasing hematocrit associated with decreased CL) and",
        "the stepwise covariate search retained it, but it reduced only the CL",
        "inter-occasion variability (8.5%) and not the inter-individual",
        "variability, so it was judged to lack clinical significance and was",
        "excluded from the final model (Kamp 2023 Sect. 3.3.2). The alternative",
        "covariate model reported in Table S2 has a hematocrit effect on CL of",
        "-4.94; that model is not the final model and is not encoded here.",
        "Cohort median 0.334 L/L (IQR 0.297-0.362; Table 1).",
        sep = " "
      )
    ),
    CYP3A5_EXPR = list(
      description = "CYP3A5 expresser status (1 = CYP3A5*1 carrier, 0 = CYP3A5*3/*3 non-expresser)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened on CL. CYP3A5 expression was associated with a 62.5%",
        "increase in CL and was retained by the stepwise covariate search, but",
        "it reduced the CL inter-individual variability by only 4.4% and was",
        "therefore excluded from the final model (Kamp 2023 Sect. 3.3.2). The",
        "alternative covariate model in Table S2 carries a CL - CYP3A5",
        "expressor effect of 0.625; that model is not the final model and is",
        "not encoded here. Cohort: 26 CYP3A5*3/*3 non-expressers (76.5%), 6",
        "CYP3A5*1/*3 (17.6%) and 2 CYP3A5*1/*1 (5.9%) expressers (Table 1).",
        sep = " "
      )
    ),
    PRED_DOSE = list(
      description = "Concomitant oral prednisolone daily dose at the time of AUC sampling",
      units       = "mg/day",
      type        = "continuous",
      notes       = paste(
        "Screened on CL. Univariate analysis found CL increased with",
        "increasing corticosteroid dose (p < 0.05), but the effect failed to",
        "explain the CL inter-individual variability and was not retained by",
        "the stepwise covariate search (Kamp 2023 Sect. 3.3.2). Median 10 mg",
        "(range 5-50; Table 1).",
        sep = " "
      )
    ),
    AGE = list(
      description = "Age at the time of AUC sampling",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened on CL and not retained. The paper attributes the null result",
        "to the narrow age range of the cohort (all recipients aged >= 65 y at",
        "transplantation) rather than to an absence of an age effect, and names",
        "this its one important limitation (Kamp 2023 Sect. 4).",
        sep = " "
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on CL and not retained (Kamp 2023 Sect. 2.3.2 and 3.3.2). 12 of 34 recipients (35.3%) were female (Table 1)."
    ),
    CONMED_CCB = list(
      description = "Concomitant calcium-channel blocker coadministration indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened on CL and not retained (Kamp 2023 Sect. 2.3.2 and 3.3.2).",
        "The recorded agents were nifedipine, barnidipine, lercanidipine,",
        "diltiazem and verapamil (Sect. 2.1); 20 of 34 recipients (58.8%) used",
        "one (Table 1).",
        sep = " "
      )
    ),
    DIS_DIAB = list(
      description = "Diabetes mellitus indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Collected as a proxy for delayed gastric passage and tested on the",
        "absorption lag time in an early two-compartment pilot model; no",
        "significant effect on the lag parameter was found, which is why the",
        "lag heterogeneity is carried by the latent MIX_LAGGED_ABS mixture",
        "rather than by a measured covariate (Kamp 2023 Sect. 4). 14 of 34",
        "recipients (41.2%) had diabetes (Table 1).",
        sep = " "
      )
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "tacrolimus", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    transit1 = list(
      analyte = "tacrolimus", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    transit2 = list(
      analyte = "tacrolimus", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "tacrolimus", units = "mg",
      specimen = "whole blood", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "tacrolimus", units = "mg",
      specimen = "whole blood", verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 34L,
    n_studies      = 1L,
    age_range      = ">= 65 years at transplantation",
    age_median     = "71.5 years (IQR 68.8-73.2)",
    weight_median  = "77.0 kg (IQR 68.1-84.5)",
    sex_female_pct = 35.3,
    disease_state  = "de novo kidney transplant recipients aged 65 years or older, on a tacrolimus / mycophenolate / prednisolone or tacrolimus / everolimus / prednisolone regimen after basiliximab induction",
    dose_range     = "Meltdose tacrolimus (Envarsus) 7 mg once daily from the day of transplantation, subsequently individualised by AUC-guided therapeutic drug monitoring; median dose at AUC sampling 6.0 mg (IQR 3.25-7.0)",
    regions        = "Netherlands (Leiden University Medical Center)",
    genotype       = "CYP3A5*3/*3 non-expressers 26 (76.5%), CYP3A5*1/*3 6 (17.6%), CYP3A5*1/*1 2 (5.9%)",
    renal_function = "Serum creatinine median 132 umol/L (IQR 99-177) at AUC sampling",
    co_medication  = "Mycophenolate mofetil 500 mg b.i.d.; prednisolone 50 mg b.i.d. tapered to 25 mg b.i.d. at day 4, 10 mg q.d. after day 4 and 5 mg q.d. after three months; calcium-channel blockers in 20 (58.8%). Patients with severe systemic infection were excluded, so no CYP inhibitors such as fluconazole were coadministered.",
    n_observations = "546 tacrolimus concentrations over 87 AUC profiles (median 2 per subject, range 1-5); 37 of the AUCs (42.5%) were dried blood spot profiles obtained in 20 (58.8%) of the patients",
    notes          = paste(
      "Add-on pharmacokinetic study to the OPTIMIZE trial (NCT03497196) at the",
      "Leiden University Medical Center site. 36 patients were enrolled and 2",
      "excluded after their transplantation was postponed. A full whole blood",
      "AUC0-24h (T = 0, 1, 2, 3, 4, 6, 8, 12, 24 h) was planned in the second",
      "week after tacrolimus initiation, and an abbreviated dried blood spot",
      "AUC (T = 0, 1, 2, 3, 6 h) at 6 weeks after transplantation for patients",
      "able to sample at home; additional post-discharge routine-care AUCs were",
      "also included. Baseline characteristics are in Kamp 2023 Table 1. The",
      "first AUC was obtained a median of 7.0 days post transplantation (IQR",
      "5.0-9.0), and the paper's own steady-state simulation concluded that at",
      "least 5-7 days are needed to reach 90% of steady state, so not all",
      "profiles in the fitting dataset were at steady state.",
      sep = " "
    )
  )

  ini({
    # =======================================================================
    # Final population PK model parameter estimates. Point estimates are the
    # "Final Model Estimate" column of Kamp 2023 Table 2; each is reproduced
    # verbatim as the corresponding $THETA / $OMEGA initial value in the
    # NONMEM control stream of Supplementary Data S1, which is the final
    # model with its estimates substituted back in.
    #
    # Omega values below are on the NONMEM variance scale and are taken
    # directly from the control stream's $OMEGA blocks. The relationship to
    # the CV% column of Table 2 is CV = sqrt(exp(omega^2) - 1); each is
    # checked in the trailing comment.
    # =======================================================================

    # ---- Absorption -------------------------------------------------------
    # A single first-order rate constant governs the whole absorption chain:
    # Supplementary Data S1 sets KA = TVKA * exp(ETA(1)) and then
    # KTR = KA, K23 = KTR, K34 = KTR, so the depot and both transit
    # compartments empty at the same rate. Table 2 names this parameter
    # "Ktr, oral transit constant"; the Results narrative calls the same
    # number "Ka".
    lktr <- log(0.752); label("Oral absorption transit rate constant (Ktr, 1/h)")   # Table 2 Ktr = 0.752 (RSE 21%); S1 $THETA(6)

    # Absorption lag time of the lagged mixture class (74% of the population).
    # The no-lag class has its lag structurally fixed at 0 (S1 $THETA(4) =
    # "0 FIX ; Absorption lag time for group 1") and is gated in model() by
    # MIX_LAGGED_ABS rather than by a second parameter. The mixture proportion
    # itself is metadata, not an ini() entry -- see
    # covariateData[[MIX_LAGGED_ABS]]$notes.
    ltlag <- log(2.29); label("Absorption lag time of the lagged mixture class (h)") # Table 2 Lag = 2.29 h (RSE 2%); S1 $THETA(5)

    # ---- Disposition ------------------------------------------------------
    # All disposition parameters are apparent (/F) values at the 70 kg
    # allometric reference; oral bioavailability is fixed at 1 (see lfdepot).
    lcl <- log(19.6); label("Apparent elimination clearance at 70 kg (CL/F, L/h)")   # Table 2 CL/F = 19.6 (RSE 8%); S1 $THETA(7)
    lvc <- log(123);  label("Apparent central volume of distribution at 70 kg (Vc/F, L)") # Table 2 Vc/F = 123 (RSE 38%); S1 $THETA(8)
    lq  <- log(74.9); label("Apparent intercompartmental clearance at 70 kg (Q/F, L/h)") # Table 2 Q/F = 74.9 (RSE 9%); S1 $THETA(9)

    # The peripheral volume was not identifiable and was fixed to the
    # literature value of 500 L per 70 kg (Sect. 3.3.1: "identifiability
    # issues occurred for the peripheral volume of the distribution parameter
    # estimate. We therefore fixed the peripheral volume of the distribution
    # parameter to 500 L/70 kg, based on the literature's values [10,20]").
    # Kamp 2023 reference [10] is Martial 2021 (meltdose tacrolimus in stable
    # adult liver transplant recipients, Br J Clin Pharmacol 87:4262-4272) and
    # reference [20] is Benkali 2010 (once-daily tacrolimus in renal
    # transplant recipients, Clin Pharmacokinet 49:683-692; see
    # modellib('Benkali_2010_tacrolimus')).
    lvp <- fixed(log(500)); label("Apparent peripheral volume of distribution at 70 kg (Vp/F, L); literature value taken from Martial 2021 and Benkali 2010") # Table 2 Vp/F = 500 (Fixed); S1 $THETA(10) 500 FIX

    # Oral bioavailability was fixed to 1, which is what makes every
    # clearance and volume above an apparent (/F) quantity and what the
    # paper's reference AUC, AUCref = (F * D * 1000) / CL, assumes.
    lfdepot <- fixed(log(1)); label("Oral bioavailability of the depot (F)")         # Table 2 F = 1 (Fixed); Sect. 2.3.1

    # ---- Allometric scaling (a priori, not estimated) ---------------------
    # Sect. 2.3.2: "The flow parameters CL/F and Q/F were allometrically
    # scaled as theta_i = theta_pop * (WT/70)^0.75 ... The distribution
    # parameters Vc and Vp were allometrically scaled as
    # theta_i = theta_pop * (WT/70)^1." Reported without uncertainty and
    # matching S1 (ALLO_CL = (WT/70)**0.75, ALLO_V = WT/70), so both
    # exponents are structural constants rather than estimates.
    e_wt_cl_q  <- fixed(0.75); label("Allometric exponent on body weight for CL/F and Q/F (unitless)") # Sect. 2.3.2; S1 $PK ALLO_CL
    e_wt_vc_vp <- fixed(1);    label("Allometric exponent on body weight for Vc/F and Vp/F (unitless)") # Sect. 2.3.2; S1 $PK ALLO_V

    # ---- Inter-individual variability -------------------------------------
    # S1 $OMEGA 0.312 -> CV = sqrt(exp(0.312) - 1) = 60.5%, matching Table 2.
    etalktr ~ 0.312                                                                  # Table 2 IIV Ka = 60.5% CV (RSE 19%, shrinkage 0%); S1 $OMEGA

    # S1 $OMEGA BLOCK(2): 0.0979 / 0.23 / 0.609. Diagonals reproduce the
    # Table 2 CV%: sqrt(exp(0.0979) - 1) = 32.1% for CL/F and
    # sqrt(exp(0.609) - 1) = 91.6% for Vc/F. The off-diagonal is a
    # covariance, and 0.23 / sqrt(0.0979 * 0.609) = 0.942 reproduces the
    # "94% covariance between CL/F and Vc/F" quoted in Sect. 3.3.1, which is
    # therefore the correlation coefficient.
    etalcl + etalvc ~ c(0.0979,
                        0.23, 0.609)                                                 # Table 2 IIV CL 32.1% / Vc 91.6% CV, correlation 0.94; S1 $OMEGA BLOCK(2)

    # ---- Inter-occasion variability on CL/F -------------------------------
    # S1 declares five IOV slots (ETA(4,5,6,7,8)) as $OMEGA BLOCK(1) 0.224
    # followed by four BLOCK(1) SAME, i.e. one shared variance across all
    # five occasions. sqrt(exp(0.224) - 1) = 50.1% CV, matching Table 2.
    # Occasions 2-5 are fixed to the occasion-1 value to encode SAME.
    etaiov_cl_1 ~ 0.224                                                              # Table 2 IOV CL = 50.1% CV (RSE 12%); S1 $OMEGA BLOCK(1)
    etaiov_cl_2 ~ fixed(0.224)                                                       # S1 $OMEGA BLOCK(1) SAME
    etaiov_cl_3 ~ fixed(0.224)                                                       # S1 $OMEGA BLOCK(1) SAME
    etaiov_cl_4 ~ fixed(0.224)                                                       # S1 $OMEGA BLOCK(1) SAME
    etaiov_cl_5 ~ fixed(0.224)                                                       # S1 $OMEGA BLOCK(1) SAME

    # ---- Residual error ---------------------------------------------------
    # S1 $ERROR: W = sqrt(PROPERR**2 * IPRED**2 + THETA(3)**2) with PROPERR
    # selected from THETA(1) (whole blood) or THETA(2) (DBS) on the DBS0
    # flag, and $SIGMA 1 FIX. The proportional terms are therefore SDs on
    # the fraction scale and switch per observation.
    propSdVenous    <- 0.208; label("Proportional residual error SD for venous whole blood samples (fraction)") # Table 2 Proportional Error Whole blood = 20.8% (RSE 8%); S1 $THETA(1)
    propSdCapillary <- 0.307; label("Proportional residual error SD for dried blood spot samples (fraction)")   # Table 2 DBS = 30.7% (RSE 12%); S1 $THETA(2)

    # The additive term was fixed at 1e-4 ug/L in S1 ($THETA(3) "0.0001 FIX")
    # purely to keep W away from zero; it is not reported in Table 2 and is
    # numerically negligible against a ug/L observation scale.
    addSd <- fixed(1e-4); label("Additive residual error SD (ug/L); numerical stabiliser") # S1 $THETA(3) 0.0001 FIX
  })

  model({
    # ---- 1. Mixture and occasion indicators -------------------------------
    # Absorption-lag mixture (S1 $MIX, Table 2 "Pop parameter" 26%): the
    # no-lag class is the 26% minority and the lagged class the 74%
    # majority. MIX_LAGGED_ABS = 1 selects the lagged class and applies the
    # 2.29 h lag; MIX_LAGGED_ABS = 0 gives no lag.
    tlag <- exp(ltlag) * MIX_LAGGED_ABS

    # Five-occasion IOV on apparent clearance (S1 $ABBREVIATED REPLACE
    # ETA(OCC_CL) = ETA(4,5,6,7,8)).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 +
      oc4 * etaiov_cl_4 + oc5 * etaiov_cl_5

    # ---- 2. Individual parameters -----------------------------------------
    # S1 $PK: CL = TVCL * ALLO_CL * EXP(ETA(2) + ETA(OCC_CL));
    #         V2 = TVV2 * ALLO_V  * EXP(ETA(3));
    #         Q  = TVQ  * ALLO_CL; V3 = TVV3 * ALLO_V.
    # Q and Vp carry no eta.
    ktr <- exp(lktr + etalktr)
    cl  <- exp(lcl + etalcl + iov_cl) * (WT / 70)^e_wt_cl_q
    vc  <- exp(lvc + etalvc)          * (WT / 70)^e_wt_vc_vp
    q   <- exp(lq)                    * (WT / 70)^e_wt_cl_q
    vp  <- exp(lvp)                   * (WT / 70)^e_wt_vc_vp

    # ---- 3. Micro-constants -----------------------------------------------
    # S1 $PK: K40 = CL/V2, K45 = Q/V2, K54 = Q/V3.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- 4. ODE system ----------------------------------------------------
    # S1 $MODEL / $DES. The dosing compartment (TRANS1, DEFDOSE) and the two
    # absorption transit compartments (TRANS2, TRANS3) all empty at ktr, so
    # the input into central is a three-stage Erlang chain. Figure 2 names
    # the first of the three the depot and the other two "Absorption transit
    # 1" and "Absorption transit 2".
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(central)     <-  ktr * transit2 - kel * central -
      k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central  - k21 * peripheral1

    # ---- 5. Bioavailability and lag ---------------------------------------
    f(depot)    <- exp(lfdepot)
    alag(depot) <- tlag

    # ---- 6. Observation and error -----------------------------------------
    # S1 $PK sets S4 = V2/1000 with doses in mg and Vc in L, so the scaled
    # prediction is in ug/L (= ng/mL), the unit used throughout Kamp 2023.
    Cc <- 1000 * central / vc

    # The proportional SD switches per observation between venous whole blood
    # and dried blood spot sampling (S1 $ERROR on the DBS0 flag).
    propSd <- propSdCapillary * SAMPLE_CAPILLARY +
      propSdVenous * (1 - SAMPLE_CAPILLARY)
    Cc ~ add(addSd) + prop(propSd)
  })
}
