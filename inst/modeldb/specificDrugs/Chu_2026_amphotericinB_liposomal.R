Chu_2026_amphotericinB_liposomal <- function() {
  description <- paste(
    "Joint population PK-PD model of liposomal amphotericin B (LAmB),",
    "oral miltefosine and skin parasite burden in 60 patients with",
    "post-kala-azar dermal leishmaniasis (PKDL) on the Indian subcontinent",
    "(Chu 2026), treated with LAmB 4 mg/kg on days 1, 4, 8, 11 and 15",
    "alone or combined with 21 days of oral miltefosine. LAmB: two",
    "compartments (plasma and the mononuclear phagocyte system, MPS) with",
    "linear clearance from plasma and SATURABLE uptake into the MPS,",
    "kin * A_central * (1 - A_MPS / Bmax(t)), with first-order release",
    "kout. The MPS capacity Bmax(t) is a body-weight-scaled physiological",
    "baseline (10 mg at 45 kg, fixed) plus a disease component of 16 mg",
    "per log10 unit of the current parasite score, so the capacity falls",
    "as treatment clears the parasites and systemic LAmB exposure rises",
    "over the course of treatment. Miltefosine: two-compartment oral model",
    "allometrically scaled on fat-free mass (reference 40 kg), ka fixed,",
    "with relative bioavailability a power function (exponent fixed at",
    "-0.13) of the cumulative miltefosine dose in mg/kg. Parasite score",
    "(skin qPCR parasites/ug DNA times PKDL lesion score) is modelled on",
    "the log10 scale: the log10 score (baseline 4.7) decays first-order",
    "toward zero at a rate proportional to the LAmB amount in the MPS plus",
    "a rate proportional to the miltefosine plasma concentration, with no",
    "growth term, and is observed with additive error."
  )
  reference <- paste(
    "Chu WY, Singh OP, Sundar S, Mondal D, Pandey K, Das P, Roseboom IC,",
    "Raja S, Torres A, Carrillo E, Huitema ADR, Alves F, Dorlo TPC.",
    "Bidirectional Interaction Between Liposomal Amphotericin B",
    "Pharmacokinetics and Parasite Dynamics in Patients With",
    "Post-Kala-Azar Dermal Leishmaniasis: Potential Implications for",
    "Optimal Dosing. Clin Pharmacol Ther. 2026;119(2):437-446.",
    "doi:10.1002/cpt.70124. Trial registration CTRI/2017/04/008421."
  )
  vignette <- "Chu_2026_amphotericinB_liposomal"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  dosing <- c("central", "depot_miltefosine")

  # mps: the paper's own structural element (Eq. 1, Figure 2), the
  # macrophage pool that takes up liposomes; no canonical compartment
  # carries that role. log10_parasites: the parasite-score state is held
  # on the log10 scale (see model()), so the canonical linear-scale
  # `parasites` pool does not describe it.
  paper_specific_compartments <- c("mps", "log10_parasites")

  compartmentData <- list(
    central = list(analyte = "amphotericin B (total, liposomal)", units = "mg", specimen = "plasma", verified = TRUE),
    mps = list(analyte = "amphotericin B (total, liposomal)", units = "mg", specimen = "tissue", verified = TRUE),
    depot_miltefosine = list(analyte = "miltefosine", units = "mg", specimen = "administration site", verified = TRUE),
    central_miltefosine = list(analyte = "miltefosine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_miltefosine = list(analyte = "miltefosine", units = "mg", specimen = "plasma", verified = TRUE),
    log10_parasites = list(
      analyte = "Leishmania donovani parasite score (skin qPCR parasites/ug DNA x PKDL lesion score)",
      units = "log10((parasites/ug DNA) x (squares on the body map))",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric power scaling of LAmB clearance (exponent fixed 0.75),",
        "LAmB central volume (exponent fixed 1) and the physiological",
        "baseline MPS capacity Bmax,base (exponent fixed 1), normalized to",
        "the cohort median of 45 kg (Table 2 footnote a; Results).",
        "Estimating the exponents did not improve the fit."
      ),
      source_name = "WT"
    ),
    FFM = list(
      description = "Fat-free mass",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric power scaling of miltefosine CL/F (exponent fixed",
        "0.75) and of Vc/F and Vp/F (exponent fixed 1), normalized to the",
        "cohort median of 40 kg (Table 2 footnote c; Supplementary Methods",
        "S3). Q/F is not scaled: footnote c marks only CL/F, Vc/F and",
        "Vp/F. The paper does not name the FFM formula; Table 1 reports a",
        "median of 39.2 kg (IQR 34.1-45)."
      ),
      source_name = "FFM"
    ),
    DOSE_MF_CUM_MGKG = list(
      description = "Cumulative administered miltefosine dose per kg body weight",
      units = "mg/kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "TIME-VARYING: the sum of the weight-normalized daily miltefosine",
        "doses given up to and including the current day (Supplementary",
        "Methods S3, 'CumDose represents cumulative miltefosine dose",
        "(mg/kg), derived from the sum of weight normalized daily dose').",
        "Drives relative bioavailability F = 1 * CumDose^-0.13 with no",
        "normalizing value (Equation S1), so the value carried on a dose",
        "record must already include that day's dose and must be",
        "positive -- at CumDose = 0 the power function is infinite. Only",
        "evaluated on miltefosine dose records; LAmB-monotherapy subjects",
        "can carry any positive placeholder."
      ),
      source_name = "CumDose"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on all PK parameters; no significant effect (Results).",
      source_name = "AGE"
    ),
    BMI = list(
      description = "Body mass index (nutritional status)",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on all PK parameters as a nutritional-status marker; no significant effect (Results).",
      source_name = "BMI"
    ),
    REGION_BANGLADESH = list(
      description = "Geographical region (Bangladesh vs India)",
      units = "(binary)",
      type = "binary",
      reference_category = "India",
      notes = paste(
        "Screened on all PK parameters (no effect) and on the baseline",
        "parasite score, where a lower Bangladesh baseline improved the",
        "OFV significantly but was dropped for model instability and low",
        "precision; no effect on parasite clearance (Results). No",
        "estimate is reported."
      ),
      source_name = "region"
    ),
    CONMED_MILTEFOSINE = list(
      description = "Treatment arm (LAmB plus miltefosine vs LAmB monotherapy)",
      units = "(binary)",
      type = "binary",
      reference_category = "LAmB monotherapy",
      notes = "Screened on CL, Vc and kin of LAmB; no significant difference between arms (Results).",
      source_name = "arm"
    ),
    LESION_TYPE = list(
      description = "PKDL lesion type (macular, papular/nodular, mixed)",
      units = "(categorical)",
      type = "categorical",
      reference_category = NULL,
      notes = "Screened on the baseline parasite score and parasite clearance; no significant effect (Results).",
      source_name = "lesion type"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 60L,
    n_studies = 1L,
    age_range = "median 22.5 years (IQR 15.8-36.3); 18 (30%) aged 18 years or younger",
    weight_range = "median 45.5 kg (IQR 40.3-55)",
    ffm_range = "median 39.2 kg (IQR 34.1-45)",
    height_range = "median 155 cm (IQR 149-162)",
    bmi_range = "median 19.0 kg/m^2 (IQR 17.5-20.6)",
    sex_female_pct = 53.3,
    race_ethnicity = "South Asian (India 50, Bangladesh 10); no further breakdown reported",
    disease_state = paste(
      "Post-kala-azar dermal leishmaniasis (Leishmania donovani). Baseline",
      "skin qPCR median 1050 parasites/ug DNA (IQR 156-5920); baseline",
      "PKDL lesion score median 175 squares on the body map (IQR",
      "75.5-262)."
    ),
    dose_range = paste(
      "LAmB (AmBisome) 4 mg/kg IV on days 1, 4, 8, 11 and 15 (total 20",
      "mg/kg) in both arms; in the combination arm plus oral miltefosine",
      "twice daily for 21 days -- allometric dosing for 30 kg or less,",
      "100 mg/day for >30 to 44 kg, 150 mg/day for 45 kg and above",
      "(about 2.5 mg/kg/day, maximum 150 mg/day)."
    ),
    regions = "India (KAMRC Muzaffarpur; RMRIMS Patna) and Bangladesh (icddr,b)",
    co_medication = "Randomized 30/30 to LAmB monotherapy or LAmB plus miltefosine.",
    samples = paste(
      "646 total amphotericin B plasma samples (66 below the 0.5 mg/L LLOQ,",
      "set to LLOQ/2), 143 miltefosine plasma samples and 177 parasite",
      "scores (baseline, day 30 and 3-month follow-up)."
    ),
    notes = paste(
      "Phase II randomized trial CTRI/2017/04/008421. NONMEM 7.5, FOCE-I;",
      "parameter precision by sampling importance resampling. The",
      "miltefosine sub-model follows the earlier Sudanese PKDL model of",
      "the same group (Chu 2026 reference 29), from which ka and the",
      "cumulative-dose exponent are fixed. Skin amphotericin B",
      "concentrations were measured but excluded from this analysis."
    )
  )

  ini({
    # ============================================================
    # LAmB disposition -- Table 2 'Pharmacokinetics of liposomal
    # amphotericin B (LAmB)'. CL, Vc and Bmax,base are given for the
    # median body weight of 45 kg (footnote a). Amounts are mg and
    # volumes L, so concentrations are mg/L.
    # ============================================================
    lcl <- log(0.3)
    label("LAmB clearance from the central compartment CL at 45 kg (L/h)") # Table 2 'CL LAmB' 0.3 L/h (SIR 95% CI 0.26-0.33); Results 0.30 L/h
    lvc <- log(3.1)
    label("LAmB central volume Vc at 45 kg (L)") # Table 2 'Vc LAmB' 3.1 L (2.8-3.4)
    lkin_mps <- log(0.53)
    label("First-order rate constant of LAmB uptake into the MPS kin (1/h)") # Table 2 'k in' 0.53 1/h (0.44-0.65)
    lkout_mps <- log(0.018)
    label("First-order rate constant of LAmB release from the MPS kout (1/h)") # Table 2 'k out' 0.018 1/h (0.014-0.024)
    lbmax_mps <- fixed(log(10))
    label("Physiological baseline MPS uptake capacity Bmax,base at 45 kg (mg)") # Table 2 'B max,base' 10 mg (fixed); Results 'first estimated and later fixed to 10 mg'
    lbmax_mps_disease <- log(16)
    label("Disease-related MPS capacity Bmax,disease (mg per log10 unit of parasite score)") # Table 2 'B max,disease' 16 mg/log10(parasite score) (15-18)

    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of body weight on LAmB CL (unitless)") # Table 2 footnote a; Results 'fixed exponent of 0.75, 1 and 1'
    e_wt_vc <- fixed(1)
    label("Allometric exponent of body weight on LAmB Vc (unitless)") # Table 2 footnote a
    e_wt_bmax_mps <- fixed(1)
    label("Allometric exponent of body weight on Bmax,base (unitless)") # Table 2 footnote a

    # ============================================================
    # Miltefosine -- Table 2 'Pharmacokinetics of miltefosine (MF)'.
    # CL/F, Vc/F and Vp/F at the median fat-free mass of 40 kg
    # (footnote c); Q/F is unscaled. ka and the cumulative-dose
    # exponent are fixed to the earlier Sudanese PKDL model (footnote
    # d). TVF fixed to 1 (Supplementary Methods S3).
    # ============================================================
    lcl_miltefosine <- log(0.086)
    label("Miltefosine apparent clearance CL/F at FFM 40 kg (L/h)") # Table 2 'CL/F MF' 0.086 L/h (0.079-0.095)
    lvc_miltefosine <- log(20)
    label("Miltefosine apparent central volume Vc/F at FFM 40 kg (L)") # Table 2 'Vc/F MF' 20 L (18-22)
    lq_miltefosine <- log(0.0012)
    label("Miltefosine apparent intercompartmental clearance Q/F (L/h)") # Table 2 'Q/F MF' 0.0012 L/h (0.0009-0.0017)
    lvp_miltefosine <- log(1.9)
    label("Miltefosine apparent peripheral volume Vp/F at FFM 40 kg (L)") # Table 2 'Vp/F MF' 1.9 L (1.3-3.3)
    lka_miltefosine <- fixed(log(0.225))
    label("Miltefosine first-order absorption rate constant ka (1/h)") # Table 2 'k a' 0.225 1/h (fixed, footnote d); Supplementary Methods S3
    lfdepot_miltefosine <- fixed(log(1))
    label("Miltefosine typical relative bioavailability TVF (unitless)") # Supplementary Methods S3 'TVF ... was fixed to 1'
    e_dose_mf_cum_mgkg_fdepot <- fixed(-0.13)
    label("Power exponent of cumulative miltefosine dose (mg/kg) on relative bioavailability (unitless)") # Table 2 'COV F,CumDose' -0.13 (fixed, footnote d); Equation S1

    e_ffm_cl_miltefosine <- fixed(0.75)
    label("Allometric exponent of fat-free mass on miltefosine CL/F (unitless)") # Table 2 footnote c; Supplementary Methods S3
    e_ffm_vc_vp_miltefosine <- fixed(1)
    label("Allometric exponent of fat-free mass on miltefosine Vc/F and Vp/F (unitless)") # Table 2 footnote c; Supplementary Methods S3

    # ============================================================
    # Parasite dynamics -- Table 2 'Pharmacodynamics of LAmB and MF
    # on parasite score' and Eq. 4. Table 2 gives lambda_LAmB in
    # 1/(mg h), which matches Eq. 4 (it multiplies the MPS AMOUNT);
    # the Results prose writes 'L/mg h' for it, a units slip.
    # The parasite-score state is the LOG10 score (see model()), so
    # the baseline 4.7 log10 units is its initial value and lrbase is
    # the natural log of 4.7.
    # ============================================================
    lslope_amphotericinb <- log(1.6e-4)
    label("Parasite clearance per mg of LAmB in the MPS lambda_LAmB (1/(mg*h))") # Table 2 'lambda LAmB' 1.6e-4 1/mg h (1.3e-4 to 2.1e-4)
    lslope_miltefosine <- log(4.9e-5)
    label("Parasite clearance per mg/L of plasma miltefosine lambda_MF (L/(mg*h))") # Table 2 'lambda MF' 4.9e-5 L/mg h (1.3e-5 to 8.2e-5)
    lrbase <- log(4.7)
    label("Baseline log10 parasite score (log10 units)") # Table 2 'Baseline parasite score' 4.7 log10(parasite score) (4.4-5.1)

    # ============================================================
    # Between-subject variability -- Table 2 'BSV in CV%' column.
    # Methods: 'exponential variance model', so omega^2 =
    # log(1 + CV^2) throughout. CL and Vc are correlated at 92%
    # (footnote b).
    # ============================================================
    etalcl + etalvc ~ c(
      0.16966,
      0.11832, 0.09749
    )
    # Table 2 BSV CL LAmB 43% -> log(1 + 0.43^2) = 0.16966; Vc LAmB 32% -> log(1 + 0.32^2) = 0.09749; footnote b correlation 0.92 -> cov = 0.92 * sqrt(0.16966 * 0.09749) = 0.11832
    etalkin_mps ~ 0.29870
    # Table 2 BSV k in 59% -> log(1 + 0.59^2) = 0.29870
    etalbmax_mps_disease ~ 0.14842
    # Table 2 BSV B max,disease 40% -> log(1 + 0.40^2) = 0.14842
    etalcl_miltefosine ~ 0.03546
    # Table 2 BSV CL/F MF 19% -> log(1 + 0.19^2) = 0.03546
    etalslope_amphotericinb ~ 0.21519
    # Table 2 BSV lambda LAmB 49% -> log(1 + 0.49^2) = 0.21519
    etalrbase ~ 2.58400
    # Table 2 BSV baseline parasite score 350% (SIR 95% CI printed as 185-316) -> log(1 + 3.5^2) = 2.58400. The point estimate lies above its own printed CI; encoded as printed, see the vignette

    # ============================================================
    # Residual error -- Table 2 'Residual variability'.
    # ============================================================
    propSd <- 0.22
    label("Proportional residual error on LAmB plasma concentration (fraction)") # Table 2 'Proportional error in LAmB' 22% (21-24)
    addSd <- fixed(0.25)
    label("Additive residual error on LAmB plasma concentration (mg/L)") # Table 2 'Additive error in LAmB' 0.25 mg/L (fixed) = LLOQ/2
    propSd_miltefosine <- 0.27
    label("Proportional residual error on miltefosine plasma concentration (fraction)") # Table 2 'Proportional error in MF' 27% (23-31)
    addSd_log10_parasites <- 1.1
    label("Additive residual error on log10 parasite score (log10 units)") # Table 2 'Additive error in parasite score' 1.1 log10 (0.89-1.3)
  })

  model({
    # ------------------------------------------------------------
    # Individual parameters.
    # ------------------------------------------------------------
    cl <- exp(lcl + etalcl) * (WT / 45)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 45)^e_wt_vc
    kin_mps <- exp(lkin_mps + etalkin_mps)
    kout_mps <- exp(lkout_mps)
    bmax_mps <- exp(lbmax_mps) * (WT / 45)^e_wt_bmax_mps
    bmax_mps_disease <- exp(lbmax_mps_disease + etalbmax_mps_disease)

    cl_miltefosine <- exp(lcl_miltefosine + etalcl_miltefosine) * (FFM / 40)^e_ffm_cl_miltefosine
    vc_miltefosine <- exp(lvc_miltefosine) * (FFM / 40)^e_ffm_vc_vp_miltefosine
    q_miltefosine <- exp(lq_miltefosine)
    vp_miltefosine <- exp(lvp_miltefosine) * (FFM / 40)^e_ffm_vc_vp_miltefosine
    ka_miltefosine <- exp(lka_miltefosine)

    slope_amphotericinb <- exp(lslope_amphotericinb + etalslope_amphotericinb)
    slope_miltefosine <- exp(lslope_miltefosine)
    rbase <- exp(lrbase + etalrbase)

    # ------------------------------------------------------------
    # MPS capacity, Eq. 3: Bmax(t) = Bmax,base + Bmax,disease *
    # A_parasite(t). Bmax,disease carries units of mg per
    # log10(parasite score) (Table 2), so A_parasite is the log10
    # score, which is also the scale of the estimated baseline (4.7)
    # and of the residual error.
    # ------------------------------------------------------------
    bmax_t <- bmax_mps + bmax_mps_disease * log10_parasites

    # ------------------------------------------------------------
    # LAmB: linear elimination from central, saturable uptake into
    # the MPS (Eq. 1) and first-order release back to central.
    # ------------------------------------------------------------
    kel <- cl / vc
    uptake_mps <- kin_mps * central * (1 - mps / bmax_t)
    d/dt(central) <- -kel * central - uptake_mps + kout_mps * mps
    d/dt(mps) <- uptake_mps - kout_mps * mps

    # ------------------------------------------------------------
    # Miltefosine: two-compartment oral model. Relative
    # bioavailability F = TVF * CumDose^COV (Equation S1, days > 7
    # form; the first-week malabsorption term was not supported in
    # this cohort -- Supplementary Results S1).
    # ------------------------------------------------------------
    kel_miltefosine <- cl_miltefosine / vc_miltefosine
    k12_miltefosine <- q_miltefosine / vc_miltefosine
    k21_miltefosine <- q_miltefosine / vp_miltefosine
    d/dt(depot_miltefosine) <- -ka_miltefosine * depot_miltefosine
    d/dt(central_miltefosine) <- ka_miltefosine * depot_miltefosine -
      (kel_miltefosine + k12_miltefosine) * central_miltefosine +
      k21_miltefosine * peripheral1_miltefosine
    d/dt(peripheral1_miltefosine) <- k12_miltefosine * central_miltefosine -
      k21_miltefosine * peripheral1_miltefosine
    f(depot_miltefosine) <- exp(lfdepot_miltefosine) * DOSE_MF_CUM_MGKG^e_dose_mf_cum_mgkg_fdepot

    Cc <- central / vc
    Cc_miltefosine <- central_miltefosine / vc_miltefosine

    # ------------------------------------------------------------
    # Parasite score, Eq. 4: drug-induced first-order clearance driven
    # by the LAmB amount in the MPS and the plasma miltefosine
    # concentration (mg/L); no growth term. A_parasite is the log10
    # score, as in Eq. 3, so the score decays exponentially toward
    # log10 = 0 and levels off once the MPS drains -- the shape of
    # Figure 3. Integrating Eq. 4 on the linear score instead leaves
    # the typical patient near 2 log10 at 3 months against the ~0.1-0.3
    # Figure 3 shows.
    # ------------------------------------------------------------
    d/dt(log10_parasites) <- -(slope_amphotericinb * mps + slope_miltefosine * Cc_miltefosine) * log10_parasites
    log10_parasites(0) <- rbase

    Cc ~ add(addSd) + prop(propSd)
    Cc_miltefosine ~ prop(propSd_miltefosine)
    log10_parasites ~ add(addSd_log10_parasites)
  })
}
