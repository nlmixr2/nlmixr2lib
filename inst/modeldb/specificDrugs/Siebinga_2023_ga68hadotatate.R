Siebinga_2023_ga68hadotatate <- function() {
  description <- paste(
    "Semi-physiological (lumped-PBPK) population PK model for the diagnostic",
    "somatostatin-receptor radioligand [68Ga]Ga-HA-DOTATATE in patients with",
    "neuroendocrine tumors. Six compartments carry the decay-corrected peptide",
    "amount (ug): blood (1), spleen (2), kidney (3), tumor lesions (4), a lumped",
    "SSTR-expressing organ compartment (5; lungs, pancreas, stomach, thyroid and",
    "liver) and a lumped rest compartment (6). All transfers are unidirectional",
    "blood-to-tissue, with elimination by renal excretion from blood (kel) and by",
    "degradation from every tissue (kdeg). Uptake into the four SSTR-expressing",
    "compartments is receptor-capacity limited: the blood-to-tissue flux is",
    "scaled by the unoccupied fraction (1 - A / bmax) of the maximum binding",
    "capacity, and only the unbound fraction (fu = 0.69) is available for uptake.",
    "Because bmax is a concentration, a larger tumor carries a proportionally",
    "larger binding capacity. Spleen uptake carries a tumor-sink effect of total",
    "tumor volume and tumor uptake a power effect of the segmented target-lesion",
    "volume, which also sets the tumor compartment volume. Every parameter in",
    "this diagnostic model is fixed: it was built bottom-up from a previously",
    "published whole-body [68Ga]Ga-HA-DOTATATE PBPK model and patient imaging",
    "data were used only for evaluation. The companion therapeutic model is",
    "Siebinga_2023_lu177hadotatate.",
    sep = " "
  )
  reference <- paste(
    "Siebinga H, de Wit-van der Veen BJ, Beijnen JH, Stokkel MPM, Dorlo TPC,",
    "Huitema ADR, Hendrikx JJMA. Predicting [177Lu]Lu-HA-DOTATATE kidney and",
    "tumor accumulation based on [68Ga]Ga-HA-DOTATATE diagnostic imaging using",
    "semi-physiological population pharmacokinetic modeling.",
    "EJNMMI Phys. 2023;10(1):48. doi:10.1186/s40658-023-00565-4",
    sep = " "
  )
  vignette <- "Siebinga_2023_hadotatate"
  units <- list(time = "h", dosing = "ug", concentration = "ug/L")

  covariateData <- list(
    TUM_VOL_TARGET = list(
      description = paste(
        "Volume summed over the segmented target tumor lesions only, i.e. the",
        "subset of lesions carried as the model's tumor compartment. Lesions",
        "were eligible if the diameter exceeded 2 cm, to a maximum of five",
        "segmented lesions and two per organ system per patient. Delineated on",
        "[68Ga]Ga-HA-DOTATATE PET/CT with a semi-automatic 50% SUVmax threshold",
        "segmentation (Methods, 'Patient population and imaging data').",
        "Table 1 median 80.0 mL (range 7.81-212).",
        sep = " "
      ),
      units = "mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Load-bearing twice over. (1) It sets the tumor compartment volume",
        "directly: 'The tumor compartment volume was based on individual",
        "measured tumor volumes and consisted of the sum of volumes of all",
        "segmented tumors', so v_tumor = TUM_VOL_TARGET / 1000 L. Because bmax",
        "is a concentration, this also makes the tumor binding capacity",
        "proportional to tumor volume ('for compartment four, larger tumors had",
        "a higher maximal SSTR binding capacity'). (2) It drives the Equation 5",
        "power effect on tumor uptake, (TUM_VOL_TARGET / 80.0)^eff, referenced",
        "to the Table 1 cohort median of 80.0 mL, which the paper names",
        "explicitly as 'V tumor cmt, median'. For this diagnostic model eff is",
        "fixed to 1, so tumor uptake is directly proportional to target-lesion",
        "volume; the therapeutic sibling estimates eff = 0.67. Baseline value,",
        "held constant per subject.",
        sep = " "
      ),
      source_name = "tumor volume of target tumors (representing the tumor compartment)"
    ),
    TUM_VOL_TOTAL = list(
      description = paste(
        "Total tumor volume summed over ALL lesions, not only the segmented",
        "target lesions, delineated on [68Ga]Ga-HA-DOTATATE PET/CT with the same",
        "semi-automatic 50% SUVmax threshold segmentation. Table 1 median 283 mL",
        "(range 22.4-644).",
        sep = " "
      ),
      units = "mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives the tumor-sink effect of Equation 4 on spleen uptake only:",
        "kin_spleen = kin_spleen_pop * exp(0.4 * (-TUM_VOL_TOTAL_in_L)). A high",
        "tumor burden sequesters radioligand and reduces uptake into healthy",
        "tissue. The coefficient 0.4 is fixed, not estimated -- 'the extent of",
        "this effect was based on previous PBPK simulations'. The exponent is",
        "dimensional, so the volume unit is load-bearing: in litres the",
        "Table 1 range 22.4-644 mL gives a 1-23% reduction in spleen uptake,",
        "whereas millilitres would give exp(-113) and annihilate spleen uptake",
        "entirely. Litres is therefore the only self-consistent reading, and the",
        "model file converts from the canonical mL accordingly. Distinct from",
        "TUM_VOL_TARGET, which is a strict subset of these lesions and drives a",
        "different effect with a different coefficient. Baseline value, held",
        "constant per subject.",
        sep = " "
      ),
      source_name = "total tumor volume"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 9,
    n_studies = 1,
    age_range = "44-76 years",
    age_median = "70 years",
    weight_range = "55.0-108 kg",
    weight_median = "75.0 kg",
    height_median = "174 cm (range 160-189)",
    sex_female_pct = 56,
    disease_state = paste(
      "Neuroendocrine tumors (NET) with somatostatin-receptor-positive lesions",
      "on diagnostic [68Ga]Ga-HA-DOTATATE PET/CT, referred for peptide receptor",
      "radionuclide therapy. Tumor volume of target tumors median 80.0 mL",
      "(7.81-212); total tumor volume median 283 mL (22.4-644).",
      sep = " "
    ),
    renal_function = "Creatinine clearance median 68.8 mL/min (range 53.6-111)",
    dose_range = paste(
      "Single intravenous [68Ga]Ga-HA-DOTATATE administration of median 96.0 MBq",
      "(75.6-102), corresponding to a median peptide amount of 5.23 ug",
      "(3.01-9.64), given within six months before the first",
      "[177Lu]Lu-HA-DOTATATE therapy cycle.",
      sep = " "
    ),
    notes = paste(
      "Baseline demographics from Table 1. Ten patients were enrolled and nine",
      "were evaluable: patient one was excluded because all tumor lesions were",
      "under 2 cm in diameter, and for patient two the kidney data were excluded",
      "because only one kidney could be quantified. Diagnostic imaging was a",
      "single PET/CT at approximately 45 min post-injection, so there is exactly",
      "one observation per compartment per patient and no blood samples were",
      "drawn. Observations are decay-corrected peptide concentrations (ug/L)",
      "computed from the measured radioactivity per volume and the administered",
      "specific activity. Model development used a bottom-up approach: all",
      "parameters were fixed from the upstream PBPK model and the patient",
      "imaging data were used only for visual model evaluation. NONMEM 7.5.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # Structural rate constants -- Table 2, column "[68Ga]Ga-HA-DOTATATE". The
    # paper's numeric micro-constants map to the canonical role-based names as
    # 1 = blood, 2 = spleen, 3 = kidney, 4 = tumor, 5 = SSTR-expressing group,
    # 6 = rest: k10 -> kel, k12 -> kin_spleen, k13 -> kin_kidney,
    # k14 -> kin_tumor, k15 -> kin_sstr, k16 -> kin_other, and k20/k30/k40/k50/
    # k60 -> the matching kdeg_<tissue>. Every value in this diagnostic model is
    # fixed: "All model parameters were fixed, since our single-time-point data
    # were limited to estimate PK parameters, and most parameters were based on
    # the previously developed PBPK model [26]."
    # ------------------------------------------------------------------------
    lkel <- fixed(log(0.25)); label("Renal excretion rate constant from blood, k10 (1/h)")            # Table 2: k10 = 0.25; Methods "renal excretion ... fixed to ... 0.25 h-1"
    lkin_spleen <- fixed(log(0.21)); label("Blood-to-spleen uptake rate constant, k12 (1/h)")         # Table 2: k12 = 0.21
    lkin_kidney <- fixed(log(0.22)); label("Blood-to-kidney uptake rate constant, k13 (1/h)")         # Table 2: k13 = 0.22
    lkin_tumor <- fixed(log(0.11)); label("Blood-to-tumor uptake rate constant, k14 (1/h)")           # Table 2: k14 = 0.11
    lkin_sstr <- fixed(log(2.5)); label("Blood-to-SSTR-organ-group uptake rate constant, k15 (1/h)")  # Table 2: k15 = 2.5
    lkin_other <- fixed(log(1)); label("Blood-to-rest-compartment rate constant, k16 (1/h)")          # Table 2: k16 = 1

    lkdeg_spleen <- fixed(log(0.01)); label("Spleen degradation rate constant, k20 (1/h)")            # Table 2: k20 = 0.01; Methods: degradation "fixed to 0.01 h-1"
    lkdeg_kidney <- fixed(log(0.01)); label("Kidney degradation rate constant, k30 (1/h)")            # Table 2: k30 = 0.01
    lkdeg_tumor <- fixed(log(0.01)); label("Tumor degradation rate constant, k40 (1/h)")              # Table 2: k40 = 0.01
    lkdeg_sstr <- fixed(log(0.01)); label("SSTR-organ-group degradation rate constant, k50 (1/h)")    # Table 2: k50 = 0.01
    lkdeg_other <- fixed(log(0.01)); label("Rest-compartment degradation rate constant, k60 (1/h)")   # Table 2: k60 = 0.01

    lvc <- fixed(log(4)); label("Blood (central) volume, V1 (L)")                                     # Table 2: V1 = 4; Methods: organ volumes from ICRP Publication 89

    # ------------------------------------------------------------------------
    # Receptor binding capacity -- Table 2 reports bmax as a CONCENTRATION in
    # nmol/L for the four SSTR-expressing compartments. It is converted to an
    # amount inside model() by multiplying by the compartment volume, per
    # "B MAX amounts were calculated based on the compartment volume". The rest
    # compartment (6) has no bmax in Table 2 and its uptake is not gated.
    # ------------------------------------------------------------------------
    lbmax_spleen <- fixed(log(16.7)); label("Maximum spleen binding capacity (nmol/L)")               # Table 2: BMAX compartment 2 = 16.7
    lbmax_kidney <- fixed(log(6.7)); label("Maximum kidney binding capacity (nmol/L)")                # Table 2: BMAX compartment 3 = 6.7
    lbmax_tumor <- fixed(log(30)); label("Maximum tumor binding capacity (nmol/L)")                   # Table 2: BMAX compartment 4 = 30
    lbmax_sstr <- fixed(log(2.4)); label("Maximum SSTR-organ-group binding capacity (nmol/L)")        # Table 2: BMAX compartment 5 = 2.4

    fu <- fixed(0.69); label("Fraction unbound in plasma (unitless)")                                 # Table 2: fraction unbound in plasma = 0.69; Methods: "fixed to 0.69"
    mw <- fixed(1628.5); label("HA-DOTATATE molar mass (g/mol)")                                      # NOT in this paper: inherited from the on-disk upstream framework paper Siebinga 2023 EJNMMI Res 13:8 (doi:10.1186/s13550-023-00958-7), its Table 3 "Molecular weight" for [68Ga]Ga-HA-DOTATATE. Needed to reconcile bmax (nmol/L) with the ug/L observation scale.

    # ------------------------------------------------------------------------
    # Structural (covariate) effects, Equations 4 and 5. Both are fixed here.
    # ------------------------------------------------------------------------
    e_tum_vol_total_kin_spleen <- fixed(0.4); label("Tumor-sink coefficient of total tumor volume (in L) on kin_spleen (1/L)")  # Equation 4: k12,cov = k12,pop * e^(0.4 * (-Vtumortotal,i)); Methods: "the extent of this effect was based on previous PBPK simulations [26]" -> fixed
    e_tum_vol_target_kin_tumor <- fixed(1); label("Power exponent of target-lesion tumor volume on kin_tumor (unitless)")       # Equation 5 exponent; Table 2 "Tumor volume on k14" = 1; Methods: "structural effect of tumor volume on k14 ... fixed to ... 1"

    # ------------------------------------------------------------------------
    # Interindividual variability. Table 2 reports IIV as CV%; the log-normal
    # variances below are omega^2 = log(CV^2 + 1) for the Equation 2 exponential
    # model P_i = P_pop * exp(eta_i). For [68Ga]Ga-HA-DOTATATE the IIV sits on
    # the uptake rate constants (Table 2 footnote **), because at the very low
    # administered peptide amount "receptors are not close to full occupancy".
    # All are fixed: "IIV was fixed based on assumed population variability".
    # ------------------------------------------------------------------------
    etalkel ~ fix(0.0951793)         # Table 2 IIV k10 = 31.6% CV -> log(0.316^2 + 1)
    etalkin_spleen ~ fix(0.2231436)  # Table 2 IIV k12 = 50% CV -> log(0.50^2 + 1)
    etalkin_kidney ~ fix(0.2231436)  # Table 2 IIV k13 = 50% CV -> log(0.50^2 + 1)
    etalkin_tumor ~ fix(0.2231436)   # Table 2 IIV k14 = 50% CV -> log(0.50^2 + 1)
    etalkin_sstr ~ fix(0.0951793)    # Table 2 IIV k15 = 31.6% CV -> log(0.316^2 + 1)

    # ------------------------------------------------------------------------
    # Residual unexplained variability -- Equation 3 is a pure proportional
    # error model. Table 2 reports a SINGLE 31.6% CV proportional error shared
    # by every observed compartment, with no RSE (the header states RSE is
    # given for estimated values only), so it is fixed. nlmixr2 requires a
    # distinct endpoint parameter per observation channel, so the one published
    # value is replicated across the three observed compartments; this is a
    # transcription of one estimate, not three independent estimates.
    # ------------------------------------------------------------------------
    propSd_Cspleen <- fixed(0.316); label("Proportional residual error, spleen (fraction)")   # Table 2 RUV proportional error = 31.6% CV (single shared value)
    propSd_Ckidney <- fixed(0.316); label("Proportional residual error, kidney (fraction)")   # Table 2 RUV proportional error = 31.6% CV (single shared value)
    propSd_Ctumor <- fixed(0.316); label("Proportional residual error, tumor (fraction)")     # Table 2 RUV proportional error = 31.6% CV (single shared value)
  })

  model({
    # Alias the fixed structural constants to locals. rxode2 rejects an
    # expression combining more than one population parameter inside a
    # mu-referenced model ("2+ single population parameters in a single
    # mu-referenced expression"); aliasing breaks the pattern and changes
    # nothing numerically (Siebinga_2023_lu177psma617.R precedent).
    fup <- fu
    mwt <- mw
    esink <- e_tum_vol_total_kin_spleen
    epow <- e_tum_vol_target_kin_tumor

    # Compartment volumes (L). Fixed for all patients from ICRP Publication 89
    # (Table 2), except the tumor compartment, whose volume is the individual
    # segmented target-lesion volume converted from the canonical mL.
    v_blood <- exp(lvc)
    v_spleen <- 0.21    # Table 2: V2 = 0.21 L
    v_kidney <- 0.3     # Table 2: V3 = 0.3 L
    v_sstr <- 4         # Table 2: V5 = 4 L
    v_other <- 50       # Table 2: V6 = 50 L
    v_tumor <- TUM_VOL_TARGET / 1000

    # Individual parameters, Equation 2 (exponential IIV). Equation 4 applies
    # the tumor-sink effect to spleen uptake exactly as printed -- on k12 alone,
    # with the total tumor volume in litres. Equation 5 applies the power effect
    # of target-lesion volume to tumor uptake, referenced to the 80.0 mL cohort
    # median that the paper names as V tumor cmt, median.
    kel <- exp(lkel + etalkel)
    kin_spleen <- exp(lkin_spleen + etalkin_spleen) * exp(esink * (-TUM_VOL_TOTAL / 1000))
    kin_kidney <- exp(lkin_kidney + etalkin_kidney)
    kin_tumor <- exp(lkin_tumor + etalkin_tumor) * (TUM_VOL_TARGET / 80.0)^epow
    kin_sstr <- exp(lkin_sstr + etalkin_sstr)
    kin_other <- exp(lkin_other)

    kdeg_spleen <- exp(lkdeg_spleen)
    kdeg_kidney <- exp(lkdeg_kidney)
    kdeg_tumor <- exp(lkdeg_tumor)
    kdeg_sstr <- exp(lkdeg_sstr)
    kdeg_other <- exp(lkdeg_other)

    # Maximum binding capacity as an AMOUNT (ug): the Table 2 concentration in
    # nmol/L is converted to ug/L with the molar mass and then multiplied by the
    # compartment volume. For the tumor this makes the capacity scale with the
    # individual tumor volume, which is the mechanism behind "larger tumors had
    # a higher maximal SSTR binding capacity".
    bmax_spleen <- exp(lbmax_spleen)
    bmax_kidney <- exp(lbmax_kidney)
    bmax_tumor <- exp(lbmax_tumor)
    bmax_sstr <- exp(lbmax_sstr)

    bmax_amt_spleen <- bmax_spleen * mwt / 1000 * v_spleen
    bmax_amt_kidney <- bmax_kidney * mwt / 1000 * v_kidney
    bmax_amt_tumor <- bmax_tumor * mwt / 1000 * v_tumor
    bmax_amt_sstr <- bmax_sstr * mwt / 1000 * v_sstr

    # Equation 1 uptake fluxes: first-order in the blood amount, available only
    # to the unbound fraction, and gated by the unoccupied fraction of the
    # receptor pool. The rest compartment has no binding capacity in Table 2, so
    # its uptake is un-gated; the unbound-fraction restriction is still applied
    # there because "only unbound molecules were available for uptake into
    # compartments" is stated for uptake generally (see vignette Errata).
    flux_spleen <- kin_spleen * fup * blood * (1 - spleen / bmax_amt_spleen)
    flux_kidney <- kin_kidney * fup * blood * (1 - kidney / bmax_amt_kidney)
    flux_tumor <- kin_tumor * fup * blood * (1 - tumor / bmax_amt_tumor)
    flux_sstr <- kin_sstr * fup * blood * (1 - sstr / bmax_amt_sstr)
    flux_other <- kin_other * fup * blood

    # Six-compartment system (Figure 2). Transfers are unidirectional: Table 2
    # lists no return rate constants (k21, k31, k41, k51, k61) and every arrow
    # in Figure 2 other than k10 leaves the system from a tissue. Blood loses
    # exactly what each tissue gains, so the system is mass conserving apart
    # from renal excretion and tissue degradation. States are peptide amounts in
    # ug; the dose is the administered peptide amount, not the radioactivity.
    d/dt(blood) <- -kel * blood - flux_spleen - flux_kidney - flux_tumor - flux_sstr - flux_other
    d/dt(spleen) <- flux_spleen - kdeg_spleen * spleen
    d/dt(kidney) <- flux_kidney - kdeg_kidney * kidney
    d/dt(tumor) <- flux_tumor - kdeg_tumor * tumor
    d/dt(sstr) <- flux_sstr - kdeg_sstr * sstr
    d/dt(other) <- flux_other - kdeg_other * other

    # Observations are decay-corrected peptide concentrations in ug/L. Spleen,
    # kidney and tumor were quantified on the PET/CT image; blood was not
    # sampled, so Cc is provided as a prediction only and carries no residual
    # error model.
    Cc <- blood / v_blood
    Cspleen <- spleen / v_spleen
    Ckidney <- kidney / v_kidney
    Ctumor <- tumor / v_tumor

    Cspleen ~ prop(propSd_Cspleen)
    Ckidney ~ prop(propSd_Ckidney)
    Ctumor ~ prop(propSd_Ctumor)
  })
}
