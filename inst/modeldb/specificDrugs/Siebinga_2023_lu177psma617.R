Siebinga_2023_lu177psma617 <- function() {
  description <- paste(
    "Six-compartment population PK dosimetry model for the PSMA-targeted",
    "radioligand [177Lu]Lu-PSMA-617 in men with low-volume metastatic prostate",
    "cancer, built from SPECT/CT-derived organ and tumor activity plus blood",
    "samples over two treatment cycles (~3 and ~6 GBq). States are blood (1),",
    "salivary gland (2), kidney (3), liver (4), tumor (5) and a lumped",
    "other-tissue compartment (6); every state carries a decay-corrected",
    "radioactivity amount (MBq) and only the blood compartment has an estimated",
    "volume (V1 = 10.3 L), so blood is observed as a concentration (MBq/L) and",
    "the tissues as amounts. Salivary-gland uptake is saturable in the bound",
    "amount (Bmax 40.4 MBq, IIV 59.3% CV); all other tissue exchange is first",
    "order. The tumor uptake rate constant carries a power effect of total tumor",
    "volume (exponent 0.705), interindividual variability (58.8% CV) and",
    "interoccasion variability across the two cycles (43.5% CV), and the blood",
    "excretion rate constant is scaled a priori by Cockcroft-Gault creatinine",
    "clearance under an assumption of complete renal elimination. The blood",
    "compartment has two observation channels with separate residual error:",
    "blood samples, and SPECT-derived blood activity linearly recalibrated to",
    "the blood-sample scale.",
    sep = " "
  )
  reference <- paste(
    "Siebinga H, Prive BM, Peters SMB, Nagarajah J, Dorlo TPC, Huitema ADR,",
    "de Wit-van der Veen BJ, Hendrikx JJMA. Population pharmacokinetic",
    "dosimetry model using imaging data to assess variability in",
    "pharmacokinetics of 177Lu-PSMA-617 in prostate cancer patients.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(8):1060-1071.",
    "doi:10.1002/psp4.12914",
    sep = " "
  )
  vignette <- "Siebinga_2023_lu177psma617"
  units <- list(time = "h", dosing = "MBq", concentration = "MBq/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    blood          = list(analyte = "[177Lu]Lu-PSMA-617", units = NA_character_, specimen = "blood cell", verified = FALSE),
    salivary_gland = list(analyte = "[177Lu]Lu-PSMA-617", units = NA_character_, specimen = "administration site", verified = FALSE),
    kidney         = list(analyte = "[177Lu]Lu-PSMA-617", units = NA_character_, specimen = "tissue", verified = FALSE),
    liver          = list(analyte = "[177Lu]Lu-PSMA-617", units = NA_character_, specimen = "tissue", verified = FALSE),
    tumor          = list(analyte = "[177Lu]Lu-PSMA-617", units = NA_character_, specimen = "tumor", verified = FALSE),
    other          = list(analyte = "[177Lu]Lu-PSMA-617", units = NA_character_, specimen = "tissue", verified = FALSE)
  )

  covariateData <- list(
    CRCL = list(
      description = paste(
        "Creatinine clearance calculated with the Cockcroft-Gault equation and",
        "NOT normalised to body surface area (raw mL/min), as reported in",
        "Table 1 (median 87.9, range 50.2-110 mL/min).",
        sep = " "
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Added a priori on the blood excretion rate constant k10 with a linear",
        "function under the assumption of complete renal elimination (Methods,",
        "'Structural effects'). The paper reports no coefficient for this",
        "effect in Table 2, which is consistent with a slope fixed to unity:",
        "kel scales in direct proportion to CRCL. Encoded as",
        "(CRCL / 87.9)^e_crcl_kel with e_crcl_kel fixed at 1 and the reference",
        "set to the Table 1 population median, because the paper does not state",
        "a centering value and the Table 2 k10 estimate of 0.288 1/h is the",
        "typical value for the cohort. Baseline value, held constant per",
        "subject.",
        sep = " "
      ),
      source_name = "creatinine clearance"
    ),
    TUM_VOL = list(
      description = paste(
        "Total tumor volume summed over all PSMA-positive lesions, delineated",
        "on the CT component of the SPECT/CT acquisition. Reported in mL in",
        "Table 1; carried here in the canonical mm^3 (1 mL = 1000 mm^3).",
        sep = " "
      ),
      units = "mm^3",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying between cycles: determined before the start of each",
        "treatment cycle (Table 1 median 2.14 mL in cycle 1 and 0.92 mL in",
        "cycle 2; overall range 0.27-78.0 mL). Enters kin_tumor as a power",
        "function (TUM_VOL / 1730)^0.705. The reference of 1730 mm^3 is the",
        "1.73 mL 'median tumor volume' the authors used for the absorbed-dose",
        "simulations; the covariate section of the paper states only 'compared",
        "to the median'. The exponent is confirmed by the paper's own check",
        "that a two-fold volume increase gives a 1.63-fold k15 increase",
        "(2^0.705 = 1.630). Clinical volumetric measurement, not the",
        "caliper-derived preclinical variant that founded the TUM_VOL entry.",
        sep = " "
      ),
      source_name = "total tumor volume"
    ),
    OCC = list(
      description = "Treatment-cycle occasion index for interoccasion variability (1 = first cycle, 2 = second cycle).",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Each of the two [177Lu]Lu-PSMA-617 cycles (8 weeks apart) is a",
        "separate occasion (Methods, 'Dosimetry model'). Decomposed inside",
        "model() into binary indicators multiplexing the IOV etas on the tumor",
        "uptake rate constant kin_tumor.",
        sep = " "
      ),
      source_name = "treatment cycle"
    )
  )

  # Screened as structural effects but NOT retained in the final model, so they
  # are documented here rather than in covariateData (they are never referenced
  # in model()).
  covariatesDataExcluded <- list(
    PSA = list(
      description = "Baseline prostate-specific antigen (Table 1 median 1.75, range 0.43-20 ng/mL).",
      units = "ng/mL",
      type = "continuous",
      notes = paste(
        "Tested as a structural effect on the tumor uptake rate constant k15",
        "and gave a significant OFV drop (dOFV -8.38), but the effect was less",
        "pronounced than total tumor volume, goodness-of-fit plots were worse,",
        "and PSA is correlated with tumor volume; the authors retained only",
        "tumor volume (Results, 'Structural effects'). No coefficient is",
        "reported, so the effect cannot be encoded.",
        sep = " "
      ),
      source_name = "PSA"
    ),
    CYCLE = list(
      description = "Treatment-cycle counter used as a dichotomous second-cycle indicator.",
      units = "(count)",
      type = "count",
      notes = paste(
        "The second treatment cycle was tested as a dichotomous structural",
        "effect on salivary-gland uptake (k12) and on Bmax, on the hypothesis",
        "that radiation-induced salivary-gland cell death reduces uptake in",
        "later cycles; model fits did not improve and no effect was retained",
        "(Results, 'Structural effects'). The cycle index is still used, as",
        "OCC, to carry interoccasion variability on kin_tumor.",
        sep = " "
      ),
      source_name = "treatment cycle"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 10,
    n_studies = 1,
    age_range = "61-77 years",
    age_median = "67 years",
    weight_range = "59.4-97.0 kg",
    weight_median = "84.0 kg",
    height_median = "178 cm (range 174-182)",
    sex_female_pct = 0,
    disease_state = paste(
      "Low-volume early-stage metastatic prostate cancer with PSMA-positive",
      "lesions on diagnostic [68Ga]Ga-PSMA PET/CT; median Gleason score 8.5",
      "(range 7-9); median baseline PSA 1.75 ng/mL (0.43-20).",
      sep = " "
    ),
    renal_function = "Cockcroft-Gault creatinine clearance median 87.9 mL/min (range 50.2-110)",
    dose_range = paste(
      "Two intravenous cycles of [177Lu]Lu-PSMA-617 8 weeks apart: median",
      "injected activity 3064 MBq (3025-3155) in cycle 1 and 6039 MBq",
      "(4972-6073) in cycle 2.",
      sep = " "
    ),
    tumor_volume = "Total tumor volume median 2.14 mL (0.27-76.6) before cycle 1 and 0.92 mL (0.27-78.0) before cycle 2",
    regions = "Netherlands (Radboud University Medical Center; NCT03828838)",
    notes = paste(
      "Baseline demographics from Table 1. Observations comprised 180 blood",
      "samples (5, 30, 60, 120 and 180 min and 24, 48, 72 and 168 h",
      "post-injection) and 491 SPECT/CT-derived observations (scans at 1, 24,",
      "48, 72 and 168 h). All activity records were decay-corrected back to the",
      "time of injection. Kidney observations from the ~1 h scan were excluded",
      "as predominantly urinary activity, and liver blood-volume activity was",
      "subtracted from measured liver activity.",
      sep = " "
    )
  )

  ini({
    # --------------------------------------------------------------------------
    # Structural rate constants and volume -- all from Table 2 ("Final model
    # parameter estimates ... for the six-compartment model"). The paper's
    # numeric compartment indices are 1 = blood, 2 = salivary gland,
    # 3 = kidney, 4 = liver, 5 = tumor, 6 = other tissue; the canonical
    # kin_<tissue> / kout_<tissue> names below carry the same meaning without
    # relying on the index.
    # --------------------------------------------------------------------------
    lkel <- log(0.288); label("Excretion (elimination) rate constant from blood, k10 (1/h)")                    # Table 2: k10 = 0.288 (RSE 7.6%)
    lvc <- log(10.3); label("Blood (central) volume of distribution, V1 (L)")                                   # Table 2: V1 = 10.3 (RSE 4.5%)

    lkin_salivary_gland <- log(0.0238); label("Blood-to-salivary-gland uptake rate constant, k12 (1/h)")        # Table 2: k12 = 0.0238 (RSE 12.4%)
    lkout_salivary_gland <- log(0.0307); label("Salivary-gland-to-blood rate constant, k21 (1/h)")              # Table 2: k21 = 0.0307 (RSE 5.8%)
    lbmax <- log(40.4); label("Maximum salivary-gland binding capacity, Bmax (MBq)")                            # Table 2: Bmax compartment 2 = 40.4 (RSE 12.3%)

    lkin_kidney <- log(0.00867); label("Blood-to-kidney uptake rate constant, k13 (1/h)")                       # Table 2: k13 = 0.00867 (RSE 8.6%)
    lkout_kidney <- log(0.0141); label("Kidney-to-blood rate constant, k31 (1/h)")                              # Table 2: k31 = 0.0141 (RSE 4.7%)

    lkin_liver <- log(0.0238); label("Blood-to-liver uptake rate constant, k14 (1/h)")                          # Table 2: k14 = 0.0238 (RSE 7.9%)
    lkout_liver <- log(0.0283); label("Liver-to-blood rate constant, k41 (1/h)")                                # Table 2: k41 = 0.0283 (RSE 4.6%)

    lkin_tumor <- log(0.000248); label("Blood-to-tumor uptake rate constant, k15 (1/h)")                        # Table 2: k15 = 0.000248 (RSE 14.0%)
    lkout_tumor <- log(0.00902); label("Tumor-to-blood rate constant, k51 (1/h)")                               # Table 2: k51 = 0.00902 (RSE 11.4%)

    lkin_other <- log(1.05); label("Blood-to-other-tissue rate constant, k16 (1/h)")                            # Table 2: k16 = 1.05 (RSE 15.3%)
    lkout_other <- log(0.744); label("Other-tissue-to-blood rate constant, k61 (1/h)")                          # Table 2: k61 = 0.744 (RSE 7.9%)

    # --------------------------------------------------------------------------
    # Structural (covariate) effects.
    # --------------------------------------------------------------------------
    e_tum_vol_kin_tumor <- 0.705; label("Power exponent of total tumor volume on kin_tumor (unitless)")          # Table 2 "Structural effects Tumor volume on k15" = 0.705 (RSE 12.3%); 2^0.705 = 1.630 reproduces the paper's "two-fold volume -> 1.63-fold k15"
    e_crcl_kel <- fixed(1); label("Power exponent of creatinine clearance on kel (unitless)")                   # Methods "Structural effects": CRCL added a priori on k10 with a linear function assuming complete renal elimination; no coefficient in Table 2, so the slope is unity (see covariateData$CRCL$notes)

    # --------------------------------------------------------------------------
    # Blood measurement model. Alpha and beta linearly recalibrate the
    # SPECT-derived blood prediction to the blood-sample scale (Equation 1);
    # gamma is the structural measurement effect added to predicted blood
    # concentrations (Equation 2). All three were estimated in the blood-only
    # model and then held fixed for the six-compartment model.
    # --------------------------------------------------------------------------
    cal_slope_spect <- fixed(0.828); label("SPECT-to-blood-sample recalibration slope, beta (unitless)")        # Results "Blood compartment": beta = 0.828, then "These parameter values were fixed in further model development"
    cal_int_spect <- fixed(6.27); label("SPECT-to-blood-sample recalibration intercept, alpha (MBq/L)")         # Results "Blood compartment": alpha = 6.27 MBq/L, subsequently fixed
    cal_bias_blood <- fixed(0.273); label("Structural blood measurement effect, gamma (MBq/L)")                 # Results "Blood compartment": gamma = 0.273 MBq/L, subsequently fixed

    # --------------------------------------------------------------------------
    # Interindividual variability. Table 2 reports IIV as CV%; the internal
    # log-normal variances below are omega^2 = log(CV^2 + 1) (Equation 4 is the
    # exponential model P_i = P_pop * exp(eta_i)). IIV was estimated only on
    # k10, k13, k41, k15 and Bmax (Results, "Dosimetry model").
    # --------------------------------------------------------------------------
    etalkel ~ 0.029155          # Table 2 IIV k10 = 17.2% CV -> log(0.172^2 + 1)
    etalkin_kidney ~ 0.025591   # Table 2 IIV k13 = 16.1% CV -> log(0.161^2 + 1)
    etalkout_liver ~ 0.008985   # Table 2 IIV k41 =  9.5% CV -> log(0.095^2 + 1)
    etalkin_tumor ~ 0.296947    # Table 2 IIV k15 = 58.8% CV -> log(0.588^2 + 1)
    etalbmax ~ 0.301325         # Table 2 IIV Bmax compartment 2 = 59.3% CV -> log(0.593^2 + 1)

    # Interoccasion variability on the tumor uptake rate constant, one occasion
    # per treatment cycle. nlmixr2 has no NONMEM $OMEGA BLOCK(1) SAME shortcut,
    # so occasion 2 gets its own eta with the variance fixed equal to the
    # occasion-1 estimate (Jonsson_2011_ethambutol.R / Chen_2023_nemonoxacin.R
    # precedent).
    etaiov_kin_tumor_1 ~ 0.173302        # Table 2 IOV k15 = 43.5% CV -> log(0.435^2 + 1)
    etaiov_kin_tumor_2 ~ fix(0.173302)   # SAME-equivalent: equal to the occasion-1 IOV variance

    # --------------------------------------------------------------------------
    # Residual unexplained variability. Table 2 reports the proportional terms
    # directly as CV%, which for the Equation 2 / Equation 5 mixed model is the
    # proportional SD on the linear scale. The single additive term applies to
    # both blood channels (Table 2 footnote a) and is reported without an RSE or
    # confidence interval, unlike every other entry in the table, so it is
    # encoded as fixed.
    # --------------------------------------------------------------------------
    propSd <- 0.193; label("Proportional residual error, blood samples (fraction)")                            # Table 2 RUV compartment 1 (blood samples) = 19.3% CV (RSE 11.5%)
    addSd <- fixed(0.25); label("Additive residual error, blood (MBq/L)")                                      # Table 2 RUV (additive) compartment 1 = 0.25 MBq/L, no RSE / CI reported -> fixed
    propSd_Cc_spect <- 0.560; label("Proportional residual error, SPECT-derived blood activity (fraction)")    # Table 2 RUV compartment 1 (SPECT data) = 56.0% CV (RSE 12.4%)
    addSd_Cc_spect <- fixed(0.25); label("Additive residual error, SPECT-derived blood activity (MBq/L)")      # Table 2 RUV (additive) footnote a: the 0.25 MBq/L additive term applies to blood samples and SPECT data alike
    propSd_Asalivary_gland <- 0.246; label("Proportional residual error, salivary-gland activity (fraction)")  # Table 2 RUV compartment 2 = 24.6% CV (RSE 11.7%)
    propSd_Akidney <- 0.304; label("Proportional residual error, kidney activity (fraction)")                  # Table 2 RUV compartment 3 = 30.4% CV (RSE 11.2%)
    propSd_Aliver <- 0.850; label("Proportional residual error, liver activity (fraction)")                    # Table 2 RUV compartment 4 = 85.0% CV (RSE 12.9%)
    propSd_Atumor <- 0.458; label("Proportional residual error, tumor activity (fraction)")                    # Table 2 RUV compartment 5 = 45.8% CV (RSE 12.9%)
  })

  model({
    # Occasion indicators -- one occasion per treatment cycle. IOV enters the
    # tumor uptake rate constant only (Results, "Structural effects").
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_kin_tumor <- oc1 * etaiov_kin_tumor_1 + oc2 * etaiov_kin_tumor_2

    # Individual parameters. Equation 4 is the exponential IIV / IOV model.
    # CRCL scales the excretion rate constant in direct proportion (complete
    # renal elimination); total tumor volume scales the tumor uptake rate
    # constant as a power function referenced to the 1.73 mL cohort median.
    kel <- exp(lkel + etalkel) * (CRCL / 87.9)^e_crcl_kel
    vc <- exp(lvc)

    kin_salivary_gland <- exp(lkin_salivary_gland)
    kout_salivary_gland <- exp(lkout_salivary_gland)
    bmax <- exp(lbmax + etalbmax)

    kin_kidney <- exp(lkin_kidney + etalkin_kidney)
    kout_kidney <- exp(lkout_kidney)

    kin_liver <- exp(lkin_liver)
    kout_liver <- exp(lkout_liver + etalkout_liver)

    kin_tumor <- exp(lkin_tumor + etalkin_tumor + iov_kin_tumor) * (TUM_VOL / 1730)^e_tum_vol_kin_tumor
    kout_tumor <- exp(lkout_tumor)

    kin_other <- exp(lkin_other)
    kout_other <- exp(lkout_other)

    # Blood measurement-model constants. These are aliased to locals because
    # rxode2 rejects an expression that combines more than one population
    # parameter inside a mu-referenced model ("2+ single population parameters
    # in a single mu-referenced expression"); the reported parameter names in
    # that error are misleading, the offending line is the Cc_spect channel
    # below. Aliasing changes nothing numerically.
    slope_spect <- cal_slope_spect
    int_spect <- cal_int_spect
    gbias_blood <- cal_bias_blood

    # Saturable salivary-gland binding, Equation 3: uptake into the target is
    # scaled by the unoccupied fraction of the maximum binding capacity, while
    # release back to blood stays first order. Only the salivary gland was
    # retained as saturable -- Bmax could not be estimated for liver and did not
    # improve the kidney or tumor fits.
    flux_salivary_gland <- kin_salivary_gland * blood * (1 - salivary_gland / bmax)

    # Six-compartment system (Figure 1). Every state is a decay-corrected
    # radioactivity amount in MBq; only blood is converted to a concentration.
    d/dt(blood) <-
      -kel * blood -
      flux_salivary_gland + kout_salivary_gland * salivary_gland -
      kin_kidney * blood + kout_kidney * kidney -
      kin_liver * blood + kout_liver * liver -
      kin_tumor * blood + kout_tumor * tumor -
      kin_other * blood + kout_other * other
    d/dt(salivary_gland) <- flux_salivary_gland - kout_salivary_gland * salivary_gland
    d/dt(kidney) <- kin_kidney * blood - kout_kidney * kidney
    d/dt(liver) <- kin_liver * blood - kout_liver * liver
    d/dt(tumor) <- kin_tumor * blood - kout_tumor * tumor
    d/dt(other) <- kin_other * blood - kout_other * other

    # Observation channels. Blood is measured two ways: directly in blood
    # samples, and via the aortic region-of-interest on the SPECT/CT scan. The
    # SPECT channel is passed through the Equation 1 linear recalibration; both
    # channels then carry the Equation 2 structural measurement effect gamma and
    # their own proportional error over a shared additive term. Tissue channels
    # are the Equation 5 proportional-error amounts in MBq.
    Cblood <- blood / vc
    Cc <- Cblood + gbias_blood
    Cc_spect <- slope_spect * Cblood + int_spect + gbias_blood
    Asalivary_gland <- salivary_gland
    Akidney <- kidney
    Aliver <- liver
    Atumor <- tumor

    Cc ~ add(addSd) + prop(propSd)
    Cc_spect ~ add(addSd_Cc_spect) + prop(propSd_Cc_spect)
    Asalivary_gland ~ prop(propSd_Asalivary_gland)
    Akidney ~ prop(propSd_Akidney)
    Aliver ~ prop(propSd_Aliver)
    Atumor ~ prop(propSd_Atumor)
  })
}
