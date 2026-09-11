Rich_2026_momelotinib <- function() {
  description <- paste(
    "Joint parent + metabolite population pharmacokinetic model for oral",
    "momelotinib (a JAK1 / JAK2 / ACVR1 inhibitor approved for myelofibrosis",
    "with anemia) and its major active metabolite M21, in 661 participants",
    "from four studies in patients with myelofibrosis (the phase II",
    "translational-biology study GS-US-352-1672 and the phase III SIMPLIFY-1,",
    "SIMPLIFY-2 and MOMENTUM trials) plus three phase I clinical-pharmacology",
    "studies in healthy participants and in renal or hepatic impairment",
    "(Rich 2026). Momelotinib is described by a two-compartment model with six",
    "transit absorption compartments and first-order elimination; the chain is",
    "depot -> transit1 ... transit6 at the common transit rate ktr, then a",
    "slower first-order ka step into the central compartment. M21 is formed",
    "from the eliminated momelotinib in proportion to the fraction metabolised",
    "fm (held at 0.640 from a human mass-balance study) and is itself described",
    "by a two-compartment model with first-order elimination. Retained",
    "momelotinib covariates: NCI-ODWG hepatic impairment and concomitant",
    "moderate or strong CYP3A4 inducers on apparent clearance, and concomitant",
    "OATP1B1/1B3 inhibitors on relative bioavailability. Retained M21",
    "covariates: a power effect of the individual momelotinib apparent",
    "clearance and of baseline creatinine clearance on apparent M21 clearance,",
    "and NCI-ODWG hepatic impairment on fm on the logit scale. Residual error",
    "is proportional in the phase III studies and proportional-plus-additive in",
    "the phase I/II studies, for both analytes. The model also returns the",
    "total active moiety tam, the potency-weighted sum of the two",
    "concentrations that the companion exposure-response analyses use. The",
    "paper fitted momelotinib and M21 SEQUENTIALLY (the M21 run consumed the",
    "parent model's post hoc estimates); this file couples them into one",
    "rxode2 model through fm, which reproduces the published M21 exposure.",
    "See vignette Assumptions and deviations."
  )
  reference <- paste(
    "Rich B, Srinivasan M, Ho YL, Visser SAG, Ferron-Brady G, Vlasakakis G.",
    "Population pharmacokinetics and exposure-response analyses of",
    "momelotinib, its active metabolite (M21), and total active moiety in",
    "myelofibrosis.",
    "Clin Pharmacol Ther. 2026;119(3):629-641. doi:10.1002/cpt.70076.",
    "Structural and covariate parameter estimates are from Table 1 and its",
    "PK-parameter-equations footnote; the residual-error stratification, the",
    "absorption-chain topology and the simulated exposure metrics used to",
    "validate this implementation are from Supplemental Tables S3, S5 and S6",
    "(supplement file CPT-119-629-s001).",
    "The fraction metabolised fm = 0.640 originates in the human",
    "mass-balance study Zheng J et al. Drug Metab Dispos. 2018;46:237-247,",
    "doi:10.1124/dmd.117.078030, which also reports the M21 relative potency",
    "of approximately 0.4 used to form the total active moiety.",
    "The nine exposure-response regressions the same paper reports are NOT",
    "extracted: Tables S9 and S10 print no intercept for any of the seven",
    "logistic models, and neither table nor the main text reports the",
    "centering constant of the log2-transformed exposure metric, so the",
    "absolute level of every one of those regressions is unidentified from",
    "the published record. See the vignette Errata."
  )
  vignette <- "Rich_2026_momelotinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    HEPIMP_MILD = list(
      description        = "Mild hepatic impairment per National Cancer Institute Organ Dysfunction Working Group (NCI-ODWG) criteria; 1 = mild, 0 = otherwise",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal hepatic function, when HEPIMP_MOD and HEPIMP_SEV are also 0)",
      notes              = paste(
        "NCI-ODWG classification (Rich 2026 Results, Momelotinib; reference 29",
        "Patel 2004). The three indicators HEPIMP_MILD / HEPIMP_MOD /",
        "HEPIMP_SEV are mutually exclusive; all three 0 selects the",
        "normal-hepatic-function reference. Acts multiplicatively on apparent",
        "momelotinib clearance AND additively on the logit of the fraction",
        "metabolised to M21, so a single indicator moves both analytes in",
        "opposite directions on exposure. Baseline distribution in the",
        "population PK analysis set (Table S7): normal 494 (74.7%), mild 127",
        "(19.2%), moderate 29 (4.4%), severe 10 (1.5%), missing 1 (0.2%).",
        "Time-fixed at baseline in this analysis."
      ),
      source_name        = "NCI-ODWG hepatic impairment"
    ),
    HEPIMP_MOD = list(
      description        = "Moderate hepatic impairment per NCI-ODWG criteria; 1 = moderate, 0 = otherwise",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal hepatic function, when HEPIMP_MILD and HEPIMP_SEV are also 0)",
      notes              = paste(
        "NCI-ODWG group 3 (total bilirubin > 1.5-3 x ULN with any AST).",
        "See HEPIMP_MILD notes for the shared reference category, the dual",
        "action on momelotinib clearance and on logit(fm), and the baseline",
        "distribution."
      ),
      source_name        = "NCI-ODWG hepatic impairment"
    ),
    HEPIMP_SEV = list(
      description        = "Severe hepatic impairment per NCI-ODWG criteria; 1 = severe, 0 = otherwise",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal hepatic function, when HEPIMP_MILD and HEPIMP_MOD are also 0)",
      notes              = paste(
        "NCI-ODWG group 4 (total bilirubin > 3 x ULN with any AST).",
        "Only 10 of 661 participants (1.5%) were severely impaired, 8 of them",
        "from the dedicated hepatic-impairment study GS-US-352-1153",
        "(Table S7). Severe impairment raises momelotinib Cavg,ss by 110% and",
        "lowers M21 Cavg,ss by 51%, for a net 42% rise in total-active-moiety",
        "Cavg,ss (Rich 2026 Results, Simulated covariate effects).",
        "See HEPIMP_MILD notes for the shared reference category."
      ),
      source_name        = "NCI-ODWG hepatic impairment"
    ),
    CONMED_CYP3A4_IND_MOD = list(
      description        = "Concomitant moderate CYP3A4 inducer at the observation; 1 = yes, 0 = no",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant moderate CYP3A4 inducer)",
      notes              = paste(
        "Rich 2026 does not enumerate which agents were classified as moderate",
        "versus strong CYP3A4 inducers, nor the number of participants in each",
        "stratum; that reporting gap is recorded in the vignette Errata.",
        "Mutually exclusive with CONMED_CYP3A4_IND_STRONG in this model.",
        "Time-varying in principle (induction requires enzyme turnover), but",
        "the source enters it as a per-observation indicator without an",
        "induction-onset lag."
      ),
      source_name        = "Concomitant CYP3A4 inducer (moderate)"
    ),
    CONMED_CYP3A4_IND_STRONG = list(
      description        = "Concomitant strong CYP3A4 inducer at the observation; 1 = yes, 0 = no",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant strong CYP3A4 inducer)",
      notes              = paste(
        "Strong CYP3A4 induction doubles apparent momelotinib clearance",
        "(factor 2.01, Table 1), i.e. a 50% fall in Cavg,ss, consistent with",
        "the dedicated rifampin drug-drug-interaction study which showed a 46%",
        "fall in momelotinib AUCinf (Rich 2026 Discussion; reference 24 Ho",
        "2024). Rich 2026 attributes the effect to induction of CYP3A / 2C8 /",
        "2C19 jointly rather than CYP3A4 alone. See",
        "CONMED_CYP3A4_IND_MOD notes for the unreported agent list."
      ),
      source_name        = "Concomitant CYP3A4 inducer (strong)"
    ),
    CONMED_OATP1B_INH = list(
      description        = "Concomitant OATP1B1 / OATP1B3 inhibitor at the observation; 1 = yes, 0 = no",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant OATP1B1/1B3 inhibitor)",
      notes              = paste(
        "Acts on relative bioavailability rather than on clearance in the",
        "source model (Table 1: 'OATP1B1/1B3 inhibitor on relative",
        "bioavailability' = 1.64), which raises momelotinib, M21 and",
        "total-active-moiety Cavg,ss by the same 64% because the factor",
        "scales the whole absorbed dose. Consistent with the dedicated",
        "drug-drug-interaction study showing a 57% rise in momelotinib AUCinf",
        "with single-dose rifampin, evidence that momelotinib is an OATP1B1/1B3",
        "substrate (Rich 2026 Discussion; reference 24 Ho 2024). Rich 2026",
        "does not enumerate the specific inhibitors pooled into the 1",
        "category; recorded in the vignette Errata."
      ),
      source_name        = "Concomitant OATP1B1/1B3 inhibitor"
    ),
    CRCL_BASE = list(
      description        = "Baseline creatinine clearance by the Cockcroft-Gault equation",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Enters apparent M21 clearance as the power",
        "term (CRCL_BASE / 72)^0.418, where 72.0 mL/min is the median baseline",
        "creatinine clearance of the four myelofibrosis studies and is the",
        "reference used for the published forest plots (Supplemental Methods,",
        "Simulations and forest plots). Baseline distribution in the population",
        "PK analysis set (Table S7): median 75.9 mL/min, mean 82.6 (SD 32.0),",
        "range 21.6-229. Momelotinib clearance itself is NOT a function of",
        "renal function -- M21 clearance is primarily renal while momelotinib",
        "is cleared hepatically (Rich 2026 Discussion)."
      ),
      source_name        = "CrCL at baseline"
    ),
    STUDY_PHASE3 = list(
      description        = "Phase III study indicator; 1 = SIMPLIFY-1, SIMPLIFY-2 or MOMENTUM, 0 = the phase I or phase II studies",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (phase I clinical-pharmacology studies GS-US-352-1151 / -1152 / -1153 and the phase II study GS-US-352-1672)",
      notes              = paste(
        "Selects the residual-error model only; it has no effect on any",
        "structural or covariate parameter. Rich 2026 estimated a",
        "proportional-only residual error for the phase III studies (which",
        "contributed mostly trough samples) and a combined",
        "proportional-plus-additive error for the phase I/II studies (which",
        "contributed the rich profiles). Table S3 shows the split was accepted",
        "on a 242.5-point objective-function drop for one parameter. Simulating",
        "a rich profile should use STUDY_PHASE3 = 0; reproducing phase III",
        "trough scatter should use STUDY_PHASE3 = 1."
      ),
      source_name        = "Study phase"
    )
  )

  compartmentData <- list(
    depot            = list(analyte = "momelotinib", units = "mg", specimen = "administration site", verified = TRUE),
    transit1         = list(analyte = "momelotinib", units = "mg", specimen = "administration site", verified = TRUE),
    transit2         = list(analyte = "momelotinib", units = "mg", specimen = "administration site", verified = TRUE),
    transit3         = list(analyte = "momelotinib", units = "mg", specimen = "administration site", verified = TRUE),
    transit4         = list(analyte = "momelotinib", units = "mg", specimen = "administration site", verified = TRUE),
    transit5         = list(analyte = "momelotinib", units = "mg", specimen = "administration site", verified = TRUE),
    transit6         = list(analyte = "momelotinib", units = "mg", specimen = "administration site", verified = TRUE),
    central          = list(analyte = "momelotinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1      = list(analyte = "momelotinib", units = "mg", specimen = "plasma", verified = TRUE),
    central_m21      = list(analyte = "M21 (momelotinib metabolite)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_m21  = list(analyte = "M21 (momelotinib metabolite)", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 661L,
    n_studies      = 7L,
    n_observations = "4508 measurable momelotinib and 4516 measurable M21 plasma concentrations (Rich 2026 Results, Population PK model development)",
    age_range      = "18-92 years",
    age_median     = "66 years",
    weight_range   = "34.2-136 kg",
    weight_median  = "74.0 kg",
    sex_female_pct = 38.6,
    race_ethnicity = c(White = 82.9, Black = 3.9, Asian = 5.1, Other = 1.8, Missing = 6.2),
    disease_state  = paste(
      "547 patients with intermediate- or high-risk myelofibrosis (primary,",
      "post-essential-thrombocythemia or post-polycythemia-vera) from the",
      "phase II study GS-US-352-1672 and the phase III studies SIMPLIFY-1",
      "(JAK-inhibitor-naive), SIMPLIFY-2 and MOMENTUM (both",
      "JAK-inhibitor-experienced; MOMENTUM also symptomatic and anemic), plus",
      "114 healthy participants or participants with renal or hepatic",
      "impairment from the phase I studies GS-US-352-1151 / -1152 / -1153."
    ),
    renal_function = "Baseline creatinine clearance (Cockcroft-Gault) median 75.9 mL/min, mean 82.6 (SD 32.0), range 21.6-229 mL/min (Table S7)",
    hepatic_function = "NCI-ODWG: normal 494 (74.7%), mild 127 (19.2%), moderate 29 (4.4%), severe 10 (1.5%), missing 1 (0.2%) (Table S7)",
    dose_range     = "100-200 mg once daily oral commercial tablet, plus single 200 mg doses in the phase I drug-drug-interaction and organ-impairment studies; 200 mg once daily in all three phase III trials (Table S1)",
    regions        = "Global (international multicentre phase III programme)",
    notes          = paste(
      "Baseline demographics are from Supplemental Table S7 (population PK",
      "analysis set, N = 661). Only participants who received the commercial",
      "tablet formulation and had at least one measurable postdose",
      "concentration were included. Race percentages are of the full 661 and",
      "include a 6.2% missing stratum, so they do not sum to 100 over the",
      "reported categories. The companion exposure-response analyses used a",
      "subset of 417 momelotinib-randomized phase III patients (Table S8)."
    )
  )

  # Rich 2026 estimated a different residual-error model for the phase III
  # studies than for the phase I/II studies, for each analyte. nlmixr2 takes a
  # single residual-SD symbol per output, so the stratum-specific SDs are
  # separate ini() parameters combined into one symbol inside model() using the
  # STUDY_PHASE3 indicator. This is the Friberg_2012_voriconazole.R pattern
  # (expSdStdy1 / expSdStdy2 / expSdStdy34), declared here so
  # checkModelConventions() does not read the stratum suffixes as deviant
  # residual-error names.
  paper_specific_residual_sds <- c(
    "propSdPh3", "propSdPh12", "addSdPh12",
    "propSdPh3_m21", "propSdPh12_m21", "addSdPh12_m21"
  )

  ini({
    # ==================================================================
    # MOMELOTINIB (parent) STRUCTURAL PARAMETERS
    # Rich 2026 Table 1, 'Momelotinib / Typical values'. All parameters
    # are apparent (per unit bioavailability) because every study dosed
    # the oral tablet.
    # ==================================================================

    lcl <- log(64.7)
    label("Momelotinib apparent clearance CL/F (L/h)")                      # Table 1 CL/F = 64.7 (RSE 3.48%, 95% CI 60.3 to 69.2)

    lq <- log(36.1)
    label("Momelotinib apparent intercompartmental clearance Q/F (L/h)")    # Table 1 Q/F = 36.1 (RSE 7.57%, 95% CI 30.7 to 41.4)

    # Rich 2026 parameterised disposition as a TOTAL apparent volume with an
    # estimated central fraction (Table S3 run mmb-2cmt-erlang6-ka-altv,
    # 'Total V parametrization'), and put ONE random effect on the total.
    # The Table 1 footnote gives:
    #   Vc_i = (0.289) x (383) x exp(eta_Vtot,i)
    #   Vp_i = (1 - 0.289) x (383) x exp(eta_Vtot,i)
    # so the two canonical volumes below carry the printed total volume and
    # central fraction literally, and share the single etalvc random effect
    # inside model().
    lvc <- log(0.289 * 383)
    label("Momelotinib apparent central volume Vc/F (L)")                   # Table 1 Fraction V central = 0.289 (RSE 9.60%) x Total V/F = 383 L (RSE 5.87%); footnote Vc_i

    lvp <- log((1 - 0.289) * 383)
    label("Momelotinib apparent peripheral volume Vp/F (L)")                # Table 1 (1 - 0.289) x Total V/F = 383 L; footnote Vp_i

    lktr <- log(8.63)
    label("Momelotinib absorption transit rate constant ktr (1/h)")         # Table 1 ktr = 8.63 (RSE 6.88%, 95% CI 7.46 to 9.79)

    lka <- log(0.303)
    label("Momelotinib absorption rate constant ka out of the last transit compartment (1/h)")  # Table 1 ka = 0.303 (RSE 8.53%, 95% CI 0.253 to 0.354)

    # Relative bioavailability anchor. Table 1 footnote: 'F_rel,i,j = 1.64 if
    # concomitant use of OATP1B1/1B3 inhibitors (otherwise 1)', so the
    # reference F_rel is structurally 1 and only the OATP1B1/1B3 effect below
    # moves it. Absolute bioavailability is not identifiable from oral data.
    lfdepot <- fixed(log(1))
    label("Momelotinib relative bioavailability F_rel reference value (unitless)")  # Table 1 footnote: F_rel = 1 in the absence of an OATP1B1/1B3 inhibitor

    # ==================================================================
    # MOMELOTINIB COVARIATE EFFECTS
    # Table 1 prints these as multiplicative factors on CL/F (and on
    # F_rel). The Supplemental Methods categorical-covariate form is
    #   P(X_i) = Pbar x exp(sum_j delta_j x I(X_i = c_j)),
    # so each delta is the natural log of the printed factor and is
    # written below as log(<printed factor>) to keep the paper's number
    # literally in the file.
    # ==================================================================

    e_hepimp_mild_cl <- log(0.914)
    label("Log factor on momelotinib CL/F for mild NCI-ODWG hepatic impairment (unitless)")      # Table 1 'Mild hepatic dysfunction (NCI) on CL/F' = 0.914 (RSE 6.93%)

    e_hepimp_mod_cl <- log(0.779)
    label("Log factor on momelotinib CL/F for moderate NCI-ODWG hepatic impairment (unitless)")  # Table 1 'Moderate hepatic dysfunction (NCI) on CL/F' = 0.779 (RSE 13.3%)

    e_hepimp_sev_cl <- log(0.477)
    label("Log factor on momelotinib CL/F for severe NCI-ODWG hepatic impairment (unitless)")    # Table 1 'Severe hepatic dysfunction (NCI) on CL/F' = 0.477 (RSE 21.8%)

    e_cyp3a4_ind_mod_cl <- log(1.39)
    label("Log factor on momelotinib CL/F for a concomitant moderate CYP3A4 inducer (unitless)") # Table 1 'Moderate CYP3A4 inducer on CL/F' = 1.39 (RSE 12.6%)

    # NOTE the Table 1 footnote mis-prints this effect as 'x1.64 if concomitant
    # use of strong CYP3A inducers', reusing the OATP1B1/1B3 bioavailability
    # factor. The Table 1 body row is authoritative at 2.01, and the Results
    # narrative confirms it: a factor of 2.01 on CL/F gives 1/2.01 = 0.498,
    # i.e. the '50% lower Cavg,ss' with strong inducers that Rich 2026 reports
    # (the mis-printed 1.64 would give 39%). See vignette Errata.
    e_cyp3a4_ind_strong_cl <- log(2.01)
    label("Log factor on momelotinib CL/F for a concomitant strong CYP3A4 inducer (unitless)")   # Table 1 'Strong CYP3A4 inducer on CL/F' = 2.01 (RSE 3.90%, 95% CI 1.86 to 2.17)

    e_oatp1b_inh_fdepot <- log(1.64)
    label("Log factor on momelotinib relative bioavailability for a concomitant OATP1B1/1B3 inhibitor (unitless)")  # Table 1 'OATP1B1/1B3 inhibitor on relative bioavailability' = 1.64 (RSE 4.07%)

    # ==================================================================
    # MOMELOTINIB INTERINDIVIDUAL VARIABILITY
    # The Table 1 'Interindividual variability' column holds omega on the
    # STANDARD-DEVIATION scale, not the variance scale. Two independent
    # checks against Supplemental Table S5 force this reading:
    #   CL/F: an SD of 0.650 implies a Cavg,ss geometric CV of
    #     sqrt(exp(0.650^2) - 1) = 74.5%, exactly the 74.5% printed for
    #     momelotinib Cavg,ss; reading 0.650 as a variance implies 95.7%,
    #     which exceeds the TOTAL observed spread and so is falsified.
    #   CLm/F: SDs of 0.371 (M21) and 0.650 (parent), propagated through
    #     the 0.481 power term plus the observed baseline-CrCL spread,
    #     imply an M21 Cavg,ss CV of about 54%, against 49.4% observed;
    #     the variance reading implies about 85%.
    # nlmixr2 eta blocks take VARIANCES, so each printed SD is squared
    # below and the paper's number stays visible as '<sd>^2'.
    # ==================================================================

    etalcl  ~ 0.650^2   # Table 1 'On CL/F' = 0.650 omega SD (RSE 3.48%); shrinkage 9.02%
    etalvc  ~ 0.476^2   # Table 1 'On total V/F' = 0.476 omega SD (RSE 9.80%); shrinkage 48.1%; shared by Vc and Vp in model()
    etalktr ~ 0.910^2   # Table 1 'On ktr' = 0.910 omega SD (RSE 6.04%); shrinkage 52.5%
    etalka  ~ 0.532^2   # Table 1 'On ka' = 0.532 omega SD (RSE 8.27%); shrinkage 49.6%

    # ==================================================================
    # MOMELOTINIB RESIDUAL ERROR
    # Table 1 'Residual error' block. Phase III contributed mostly trough
    # samples and carries a larger proportional error with no additive
    # term; phase I/II carried the rich profiles and a combined error.
    # ==================================================================

    propSdPh3 <- 0.584
    label("Momelotinib proportional residual SD in the phase III studies (fraction)")   # Table 1 'Proportional error phase III studies' = 58.4% (RSE 2.14%)

    propSdPh12 <- 0.343
    label("Momelotinib proportional residual SD in the phase I/II studies (fraction)")  # Table 1 'Proportional error phase I/II studies' = 34.3% (RSE 1.93%)

    addSdPh12 <- 1.47
    label("Momelotinib additive residual SD in the phase I/II studies (ng/mL)")         # Table 1 'Additive error phase I/II studies' = 1.47 ng/mL (RSE 11.0%)

    # ==================================================================
    # M21 (metabolite) STRUCTURAL PARAMETERS
    # Rich 2026 Table 1, 'M21 / Typical values'.
    # ==================================================================

    lcl_m21 <- log(24.9)
    label("M21 apparent clearance CLm/F (L/h)")                             # Table 1 'CLm/F of M21' = 24.9 (RSE 2.59%, 95% CI 23.7 to 26.2)

    lvc_m21 <- log(2.77)
    label("M21 apparent central volume Vcm/F (L)")                          # Table 1 'Vcm/F of M21' = 2.77 (RSE 18.6%, 95% CI 1.76 to 3.78)

    lvp_m21 <- log(45.7)
    label("M21 apparent peripheral volume Vpm/F (L)")                       # Table 1 'Vpm/F of M21' = 45.7 (RSE 7.83%, 95% CI 38.7 to 52.8)

    lq_m21 <- log(8.84)
    label("M21 apparent intercompartmental clearance Qm/F (L/h)")           # Table 1 'Qm/F of M21' = 8.84 (RSE 9.90%, 95% CI 7.12 to 10.6)

    # Fraction of momelotinib metabolised to M21. Rich 2026 Methods,
    # Population PK model development: 'The value of fraction metabolized
    # (Fmet) in plasma was fixed at 0.64 based on a previous human
    # mass-balance study that showed M21 was the most abundant momelotinib
    # metabolite accounting for 64.2% of the AUC of total radioactivity in
    # plasma.' Held on the logit scale because the hepatic-impairment
    # effects below are additive logit shifts.
    logitfm <- fixed(log(0.640 / (1 - 0.640)))
    label("Logit of the fraction of momelotinib metabolised to M21 at normal hepatic function (unitless)")  # Table 1 footnote F_met = 0.640 if normal hepatic function; value from the Zheng 2018 mass-balance study

    # ==================================================================
    # M21 COVARIATE EFFECTS
    # ==================================================================

    e_clmmb_cl_m21 <- 0.481
    label("Power exponent of (momelotinib CL/F / 64.7) on M21 CLm/F (unitless)")   # Table 1 'CL/F of momelotinib on CLm/F' = 0.481 (RSE 7.66%, 95% CI 0.409 to 0.554)

    e_crcl_base_cl_m21 <- 0.418
    label("Power exponent of (CRCL_BASE / 72) on M21 CLm/F (unitless)")            # Table 1 'CrCL at baseline on CLm/F' = 0.418 (RSE 14.1%, 95% CI 0.302 to 0.533)

    e_hepimp_mild_fm <- -0.359
    label("Logit-additive shift on the fraction metabolised for mild NCI-ODWG hepatic impairment (unitless)")      # Table 1 'Mild hepatic dysfunction (NCI) on fraction metabolized (logit scale)' = -0.359 (RSE 38.0%)

    e_hepimp_mod_fm <- -0.679
    label("Logit-additive shift on the fraction metabolised for moderate NCI-ODWG hepatic impairment (unitless)")  # Table 1 'Moderate hepatic dysfunction (NCI) on fraction metabolized (logit scale)' = -0.679 (RSE 49.2%)

    e_hepimp_sev_fm <- -1.84
    label("Logit-additive shift on the fraction metabolised for severe NCI-ODWG hepatic impairment (unitless)")    # Table 1 'Severe hepatic dysfunction (NCI) on fraction metabolized (logit scale)' = -1.84 (RSE 10.2%)

    # ==================================================================
    # M21 INTERINDIVIDUAL VARIABILITY
    # Same SD-scale reading as the parent (see the momelotinib IIV block).
    # Vcm/F is poorly identified: 52.8% shrinkage, and because Qm/F is
    # large relative to Vcm/F the central and peripheral M21 compartments
    # equilibrate fast, so M21 exposure is nearly insensitive to Vcm/F.
    # That is why a very large omega coexists with the tight 27.9%
    # observed M21 Cmax,ss CV of Table S5.
    # ==================================================================

    etalcl_m21 ~ 0.371^2   # Table 1 'On CLm/F of M21' = 0.371 omega SD (RSE 3.71%); shrinkage 13.7%
    etalvc_m21 ~ 2.39^2    # Table 1 'On Vcm/F of M21' = 2.39 omega SD (RSE 7.18%); shrinkage 52.8%

    # ==================================================================
    # M21 RESIDUAL ERROR
    # ==================================================================

    propSdPh3_m21 <- 0.517
    label("M21 proportional residual SD in the phase III studies (fraction)")   # Table 1 M21 'Proportional error phase III studies' = 51.7% (RSE 3.02%)

    propSdPh12_m21 <- 0.337
    label("M21 proportional residual SD in the phase I/II studies (fraction)")  # Table 1 M21 'Proportional error phase I/II studies' = 33.7% (RSE 3.93%)

    addSdPh12_m21 <- 1.86
    label("M21 additive residual SD in the phase I/II studies (ng/mL)")         # Table 1 (Continued) M21 'Additive error phase I/II studies' = 1.86 ng/mL (RSE 40.2%)
  })

  model({
    # ------------------------------------------------------------------
    # Constants and reference values
    # ------------------------------------------------------------------
    # Reference momelotinib CL/F for the power term on M21 clearance, equal
    # to the typical value itself, so the term is 1 for a typical individual
    # (Table 1 footnote: (CL_i,j / 64.7)^0.481).
    ref_cl <- 64.7      # L/h
    # Reference baseline creatinine clearance: the median of the four
    # myelofibrosis studies, and the reference used for the published forest
    # plots (Supplemental Methods, Simulations and forest plots).
    ref_crcl <- 72      # mL/min
    # Relative potency of M21 versus momelotinib used to form the total
    # active moiety. Not estimated by Rich 2026 -- carried from the Zheng
    # 2018 human mass-balance / activity-index work (doi:10.1124/dmd.117.078030)
    # and stated in Rich 2026 Simulations as 'R_p ... estimated to be ~0.4'.
    rp_m21 <- 0.4       # unitless

    # ------------------------------------------------------------------
    # 1. Momelotinib individual parameters
    # ------------------------------------------------------------------
    # Table 1 footnote:
    #   CL_i,j = 64.7 x exp(eta_CL,i) x 0.914^mild x 0.779^moderate
    #            x 0.477^severe x 1.39^modInducer x 2.01^strongInducer
    # written here in the Supplemental Methods exp(sum of log factors) form.
    cl <- exp(lcl + etalcl) *
      exp(e_hepimp_mild_cl       * HEPIMP_MILD +
          e_hepimp_mod_cl        * HEPIMP_MOD +
          e_hepimp_sev_cl        * HEPIMP_SEV +
          e_cyp3a4_ind_mod_cl    * CONMED_CYP3A4_IND_MOD +
          e_cyp3a4_ind_strong_cl * CONMED_CYP3A4_IND_STRONG)

    # A SINGLE random effect acts on the TOTAL apparent volume, so the same
    # etalvc enters both volumes (Table 1 footnote Vc_i / Vp_i).
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvc)
    q  <- exp(lq)

    ktr <- exp(lktr + etalktr)
    ka  <- exp(lka + etalka)

    # Relative bioavailability: 1 at reference, 1.64 under a concomitant
    # OATP1B1/1B3 inhibitor (Table 1 footnote F_rel,i,j).
    fdepot <- exp(lfdepot) * exp(e_oatp1b_inh_fdepot * CONMED_OATP1B_INH)

    # ------------------------------------------------------------------
    # 2. M21 individual parameters
    # ------------------------------------------------------------------
    # Table 1 footnote:
    #   CLm_i,j = 24.9 x exp(eta_CLm,i) x (CL_i,j / 64.7)^0.481
    #             x (CrCL_i / 72)^0.418
    # The dependence on the individual momelotinib CL/F is how Rich 2026
    # carried the positive parent-metabolite clearance correlation of
    # Figure S7 into the sequentially fitted M21 model, in place of an
    # estimated off-diagonal omega.
    cl_m21 <- exp(lcl_m21 + etalcl_m21) *
      (cl / ref_cl)^e_clmmb_cl_m21 *
      (CRCL_BASE / ref_crcl)^e_crcl_base_cl_m21

    vc_m21 <- exp(lvc_m21 + etalvc_m21)
    vp_m21 <- exp(lvp_m21)
    q_m21  <- exp(lq_m21)

    # Fraction metabolised, with additive hepatic-impairment shifts on the
    # logit scale. Table 1 footnote b verifies the three impaired values:
    #   expit(logit(0.64) - 0.359) = 0.554
    #   expit(logit(0.64) - 0.679) = 0.474
    #   expit(logit(0.64) - 1.84)  = 0.220
    fm <- expit(logitfm +
                  e_hepimp_mild_fm * HEPIMP_MILD +
                  e_hepimp_mod_fm  * HEPIMP_MOD +
                  e_hepimp_sev_fm  * HEPIMP_SEV)

    # ------------------------------------------------------------------
    # 3. Micro-constants
    # ------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    kel_m21 <- cl_m21 / vc_m21
    k12_m21 <- q_m21  / vc_m21
    k21_m21 <- q_m21  / vp_m21

    # ------------------------------------------------------------------
    # 4. ODE system
    # ------------------------------------------------------------------
    # Absorption chain. Table S3 selected run mmb-2cmt-erlang6-ka-altv,
    # described as '6 transit comp., 6 ktr transitions followed by ka
    # transition': the oral dose lands in depot, six successive ktr
    # transitions carry it through transit1 ... transit6, and a final,
    # much slower ka step (0.303 vs 8.63 1/h) delivers it to central.
    # ktr is therefore NOT rate-limiting -- ka is.
    d/dt(depot)    <- -ktr * depot
    d/dt(transit1) <-  ktr * depot    - ktr * transit1
    d/dt(transit2) <-  ktr * transit1 - ktr * transit2
    d/dt(transit3) <-  ktr * transit2 - ktr * transit3
    d/dt(transit4) <-  ktr * transit3 - ktr * transit4
    d/dt(transit5) <-  ktr * transit4 - ktr * transit5
    d/dt(transit6) <-  ktr * transit5 - ka  * transit6

    d/dt(central)     <-  ka * transit6 - kel * central -
      k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # M21 formation. Rich 2026 fitted M21 sequentially on the parent's post
    # hoc predictions with the fraction metabolised fm as the link, so the
    # molar formation rate is fm times the momelotinib elimination rate.
    # The construction is exact at steady state: total M21 formed per
    # interval is fm x F_rel x Dose, hence
    #   AUCtau,ss(M21) = fm x F_rel x Dose / CLm,
    # which for a typical individual gives 0.640 x 200 mg / 24.9 L/h =
    # 5141 ng.h/mL against the 5080 ng.h/mL median of Table S5. No
    # molecular-weight correction is applied, matching the source (fm was
    # derived from an AUC-of-total-radioactivity share).
    d/dt(central_m21)     <-  fm * kel * central - kel_m21 * central_m21 -
      k12_m21 * central_m21 + k21_m21 * peripheral1_m21
    d/dt(peripheral1_m21) <-  k12_m21 * central_m21 - k21_m21 * peripheral1_m21

    # ------------------------------------------------------------------
    # 5. Bioavailability
    # ------------------------------------------------------------------
    f(depot) <- fdepot

    # ------------------------------------------------------------------
    # 6. Observations, total active moiety, and residual error
    # ------------------------------------------------------------------
    # Amounts are in mg and volumes in L, so dividing by (V / 1000)
    # converts mg/L to the ng/mL that Rich 2026 reports.
    Cc     <- central     / (vc     / 1000)
    Cc_m21 <- central_m21 / (vc_m21 / 1000)

    # Total active moiety (Rich 2026 Simulations):
    #   TAM(t) = C_MMB(t) + R_p x C_M21(t)
    # Reported as a derived variable rather than a fitted observation --
    # TAM was never measured, and it is the exposure metric the companion
    # exposure-response models consume as their CAV covariate.
    tam <- Cc + rp_m21 * Cc_m21

    # Stratum-specific residual error selected by STUDY_PHASE3 (see the
    # covariateData entry). Phase III is proportional-only, so the additive
    # term is switched off there.
    propSdCc     <- propSdPh3     * STUDY_PHASE3 + propSdPh12     * (1 - STUDY_PHASE3)
    addSdCc      <- addSdPh12     * (1 - STUDY_PHASE3)
    propSdCc_m21 <- propSdPh3_m21 * STUDY_PHASE3 + propSdPh12_m21 * (1 - STUDY_PHASE3)
    addSdCc_m21  <- addSdPh12_m21 * (1 - STUDY_PHASE3)

    Cc     ~ add(addSdCc)     + prop(propSdCc)
    Cc_m21 ~ add(addSdCc_m21) + prop(propSdCc_m21)
  })
}
