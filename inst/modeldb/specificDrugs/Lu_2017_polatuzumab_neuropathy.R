Lu_2017_polatuzumab_neuropathy <- function() {
  description <- "Time-to-event hazard model for the onset of grade >= 2 peripheral neuropathy (PN) during polatuzumab vedotin treatment in adults with relapsed/refractory B-cell non-Hodgkin lymphoma (Lu 2017). The PN hazard is driven by a hypothetical effect compartment receiving plasma antibody-conjugated MMAE (acMMAE) with first-order distribution k1e in and ke0 = k1e out, modulated by a Weibull time function on the drug-effect potency (alpha drug-effect, beta shape) and by twelve baseline-covariate proportional-hazard terms (age, body weight, sex, active grade 1 PN at baseline, prior radiotherapy, prior vinca alkaloid, prior platinum-based chemotherapy, rituximab combination, tumor histology DLBCL vs other-non-FL, baseline tumor sum of products of perpendicular diameters, baseline serum albumin). The acMMAE plasma driver is inlined from the published Lu 2019 integrated two-analyte popPK (acMMAE side only; see Lu_2019_polatuzumab.R) per the standing policy of reusing a published same-drug PK when the originally-used PK source (Lu 2015 ASCPT poster, unpublished) is not on disk. Both the instantaneous hazard and the cumulative hazard / survival outputs are exposed for direct VPC simulation of the Kaplan-Meier curve."
  reference <- paste(
    "Lu D, Gillespie WR, Girish S, Agarwal P, Li C, Hirata J,",
    "Chu Y-W, Kagedal M, Leon L, Maiya V, Jin JY.",
    "Time-to-event analysis of polatuzumab vedotin-induced",
    "peripheral neuropathy to assist in the comparison of",
    "clinical dosing regimens.",
    "CPT Pharmacometrics Syst Pharmacol. 2017;6(6):401-408.",
    "doi:10.1002/psp4.12192. PMID 28294568.",
    "Upstream PK driver (acMMAE side only) from:",
    "Lu D, Lu T, Gibiansky L, Li X, Li C, Agarwal P, Shemesh CS,",
    "Shi R, Dere RC, Hirata J, Miles D, Chanu P, Girish S, Jin JY.",
    "Integrated Two-Analyte Population Pharmacokinetic Model of",
    "Polatuzumab Vedotin in Patients With Non-Hodgkin Lymphoma.",
    "CPT Pharmacometrics Syst Pharmacol. 2020;9(1):48-59.",
    "doi:10.1002/psp4.12482. PMID 31749251.",
    sep = " "
  )
  vignette <- "Lu_2017_polatuzumab_neuropathy"
  units <- list(
    time          = "hour",
    dosing        = "ug",
    concentration = "ng/mL"
  )

  covariateData <- list(
    # ----- Shared between PK and TTE PD layers -----
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Two distinct reference values: (1) Lu 2019 PK side uses 75 kg as the allometric reference (power effects on acMMAE CL_SS exponent 0.73 and on Vc/Vp/Q shared exponent 0.50); (2) Lu 2017 TTE side uses 80 kg as the centering reference (linear effect on log-hazard via `e_wt_haz * (WT - 80)`). Same physical column, two independent normalizations.",
      source_name        = "BWT"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-fixed. Used on both PK side (Lu 2019 effects on Vc and CL_SS, with value-coding flip from the source NONMEM `SEX = SEXN - 1` male-indicator convention -- see Lu_2019_polatuzumab.R covariateData[[SEXF]]$notes) and TTE PD side (Lu 2017 paper's `female` indicator equals SEXF directly, log-hazard effect `e_sexf_haz * SEXF`).",
      source_name        = "SEXN"
    ),
    ALB = list(
      description        = "Baseline serum albumin concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value (SI units). Two distinct reference values: (1) Lu 2019 PK side uses 35 g/L (power effects on CL_SS exponent -0.247); (2) Lu 2017 TTE side uses 39 g/L as the centering reference (linear effect on log-hazard via `e_alb_haz * (ALB - 39)`). Same physical column, two independent normalizations.",
      source_name        = "BALBUM"
    ),
    TUMSZ = list(
      description        = "Baseline tumor sum of the products of perpendicular dimensions (SPD)",
      units              = "mm^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Two distinct reference values: (1) Lu 2019 PK side uses 5000 mm^2 for the linear effect on CL_SS and as the half-maximal Michaelis-Menten anchor for the effect on CL_TIME (`TUMSZ / (1150 + TUMSZ)`); (2) Lu 2017 TTE side uses 3000 mm^2 as the log-normalization reference (linear effect on log-hazard via `e_tumsz_haz * log(TUMSZ / 3000)`). Same physical column, two independent normalizations.",
      source_name        = "BTMBD"
    ),
    CONMED_RITUX = list(
      description        = "Concomitant rituximab combination indicator, 1 = pola co-administered with rituximab, 0 = single-agent pola",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-fixed in the Lu 2017 cohort. Used on both PK side (interpreted as the COMBO_RG anti-CD20-combination indicator; this is exact for the Lu 2017 sub-cohort because no Lu 2017 patient received obinutuzumab combination -- the dataset predates the Pola+G GO29365 / GO29044 studies that motivated the Lu 2019 obinutuzumab branch) and TTE PD side (Lu 2017 paper's `rituximab` indicator equals CONMED_RITUX directly, log-hazard effect `e_conmed_ritux_haz * CONMED_RITUX`).",
      source_name        = "COMBO"
    ),
    # ----- PK-side-only covariates (Lu 2019) -----
    LINE_1L = list(
      description        = "First-line-therapy indicator: 1 = treatment-naive, 0 = relapsed/refractory (>= second line)",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-fixed. PK-side covariate carried over from the Lu 2019 popPK; not used by the Lu 2017 TTE. Set to 0 for the all-R/R Lu 2017 cohort when simulating.",
      source_name        = "RRFN"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator, 1 = Asian, 0 = non-Asian",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-fixed. PK-side covariate carried over from the Lu 2019 popPK; not used by the Lu 2017 TTE.",
      source_name        = "RACEN"
    ),
    BLBCELL = list(
      description        = "Baseline CD19+ B cell count",
      units              = "cells/uL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. PK-side covariate carried over from the Lu 2019 popPK (two effects: power on CL_SS via `max(1, BLBCELL)^0.0212`, power on CL_TIME via `max(1, BLBCELL/121)^0.578`); not used by the Lu 2017 TTE.",
      source_name        = "BBCC"
    ),
    # ----- TTE-side-only covariates (Lu 2017) -----
    AGE = list(
      description        = "Subject age at study entry",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Centered at 65 years for the Lu 2017 TTE log-hazard effect `e_age_haz * (AGE - 65)`. Not used by the Lu 2019 PK side (Lu 2019 did not retain age as a covariate).",
      source_name        = "AGE"
    ),
    BL_PN_GR1 = list(
      description        = "Active grade 1 peripheral neuropathy at study entry indicator, 1 = active grade 1 PN at baseline, 0 = no active PN",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-fixed baseline AE status. The Lu 2017 protocol allowed enrollment of patients with grade 1 PN at study entry; this indicator captures that subset. Log-hazard effect `e_bl_pn_gr1_haz * BL_PN_GR1`. Distinct from a Markov-state AE-grade covariate (see `PREV_AE_SCORE` in the canonical register) -- this is time-fixed at baseline and does not update during the analysis window. The Lu 2017 paper reports a sensitivity analysis in which this indicator was replaced by a broader 'history of prior PN' indicator with similar (inconclusive) results.",
      source_name        = "BLPN"
    ),
    PRIOR_RADIATION = list(
      description        = "Prior radiotherapy exposure indicator, 1 = received prior radiotherapy, 0 = radiotherapy-naive",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-fixed at baseline. Log-hazard effect `e_prior_radiation_haz * PRIOR_RADIATION`.",
      source_name        = "PRRADIOTX"
    ),
    PRIOR_VINCA = list(
      description        = "Prior vinca-alkaloid chemotherapy indicator, 1 = received any prior vinca alkaloid (vincristine, vinblastine, vinorelbine, vindesine), 0 = vinca-naive",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-fixed at baseline. The Lu 2017 hypothesis is that prior antimicrotubule exposure could sensitize patients to subsequent vc-MMAE-induced PN. Log-hazard effect `e_prior_vinca_haz * PRIOR_VINCA`.",
      source_name        = "PRVINCA"
    ),
    PRIOR_PLATIN = list(
      description        = "Prior platinum-based chemotherapy indicator, 1 = received any prior platinum agent (cisplatin, carboplatin, oxaliplatin, etc.), 0 = platinum-naive",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-fixed at baseline. The Lu 2017 hypothesis is that platinum-induced sensory PN history could compound subsequent antimicrotubule-induced PN. Log-hazard effect `e_prior_platin_haz * PRIOR_PLATIN`.",
      source_name        = "PRPLATIN"
    ),
    TUMTP_DLBCL = list(
      description        = "Diffuse large B-cell lymphoma indicator, 1 = DLBCL, 0 = other tumor histology",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-fixed. Part of the Lu 2017 three-level tumor-histology decomposition (FL reference, DLBCL, otherNonFL). FL is the implicit reference when both TUMTP_DLBCL = 0 and TUMTP_OTHER_NHL = 0. Log-hazard effect `e_tumtp_dlbcl_haz * TUMTP_DLBCL`.",
      source_name        = "DLBCL"
    ),
    TUMTP_OTHER_NHL = list(
      description        = "Non-FL non-DLBCL NHL histology indicator, 1 = any NHL histology other than FL or DLBCL, 0 = FL or DLBCL",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-fixed. Captures the 'others' bucket in Lu 2017's three-level tumor-histology decomposition. Log-hazard effect `e_tumtp_other_nhl_haz * TUMTP_OTHER_NHL`. The Lu 2017 paper does not enumerate the specific histologies in this residual bucket; in the source phase I/II datasets these include marginal-zone lymphoma, small lymphocytic lymphoma, mantle cell lymphoma, and chronic lymphocytic leukemia (CLL).",
      source_name        = "OTHERNONFL"
    )
  )

  # Follicular lymphoma is the IMPLICIT reference category for tumor histology
  # in the Lu 2017 covariate model -- both TUMTP_DLBCL = 0 and TUMTP_OTHER_NHL
  # = 0 corresponds to FL. The TUMTP_FL canonical column is therefore declared
  # here as documented-but-not-referenced rather than in covariateData (which
  # would trigger a "declared but not referenced" convention warning).
  covariatesDataExcluded <- list(
    TUMTP_FL = list(
      description = "Follicular lymphoma indicator, 1 = FL, 0 = other tumor histology",
      units       = "(binary)",
      type        = "binary",
      notes       = "Implicit reference category for the Lu 2017 tumor-histology covariate decomposition (FL reference, DLBCL, otherNonFL). Recorded here so the FL reference is visible in the model file metadata without triggering a convention warning."
    )
  )

  population <- list(
    n_subjects     = 155L,
    n_studies      = 2L,
    n_events       = "First grade >= 2 PN event or right-censoring across 155 patients with R/R B-cell NHL pooled from phase I (DCS4968g, NCT01290549) and phase II (GO27834 / ROMULUS, NCT01691898) studies.",
    age_range      = "Adults with R/R B-cell NHL; Lu 2017 does not tabulate age quantiles in the main paper but the Figure 3 sensitivity panels reference 50 vs 65 years and 80 years as covariate-effect anchor points.",
    weight_range   = "Lu 2017 Figure 3 references 60 / 80 / 100 kg as covariate-effect anchor points (5th-95th percentile inferred). 80 kg used as the TTE-covariate centering reference.",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = "Relapsed/refractory B-cell non-Hodgkin lymphoma (NHL): follicular lymphoma (FL), diffuse large B-cell lymphoma (DLBCL), and other non-FL non-DLBCL histologies (mantle cell lymphoma, marginal-zone lymphoma, small lymphocytic lymphoma, transformed FL, CLL).",
    dose_range     = "Polatuzumab vedotin 0.1-2.4 mg/kg IV every 3 weeks (Q3W); the analyzed dose levels for the model-application section are 1.8 and 2.4 mg/kg Q3W. Phase II patients were treated until progression or for a maximum of 17 cycles; phase I patients were treated until progression or unacceptable toxicity. Patients at 2.4 mg/kg were permitted a dose reduction to 1.8 mg/kg for AEs.",
    regions        = "Multi-regional phase I/II studies.",
    studies        = "DCS4968g (NCT01290549, phase I single-agent and Pola+R dose-escalation; N=76 of 77 analyzed), GO27834 / ROMULUS (NCT01691898, phase II Pola+R Q3W; N=19+60 = 79 analyzed across FL 1.8 mg/kg arm, FL 2.4 mg/kg arm, and DLBCL 2.4 mg/kg arm).",
    notes          = "Data-cut July 2014. The TTE PD parameters in this file were fit by Lu 2017 NONMEM v7.3 using empirical-Bayes individual PK parameters from a previously developed (unpublished) population PK model presented at ASCPT 2015 (Lu D et al, 116th Annual Meeting, 2015; reference 19 of Lu 2017). That upstream Lu 2015 popPK is not on disk for this extraction; per the standing nlmixr2lib policy of reusing a published same-drug PK as the upstream layer, this file inlines the acMMAE side of the published Lu 2019 integrated two-analyte popPK model. The Lu 2019 PK is more sophisticated than the Lu 2015 poster popPK (Lu 2019 adds Hill-shaped time-decay on CL_NS, separate exponential decay on CL_t, Michaelis-Menten elimination, and several additional covariates) and therefore produces somewhat different acMMAE concentration profiles than were used to fit the Lu 2017 TTE PD parameters; this difference is documented in the validation vignette's Assumptions and deviations section."
  )

  ini({
    # ============================================================
    # ===== PK side: Lu 2019 integrated 2-analyte popPK (acMMAE only) =====
    # ============================================================
    # See Lu_2019_polatuzumab.R for the full integrated 2-analyte model and
    # the per-parameter source-trace comments. Reference subject: 75 kg, ALB
    # 35 g/L, TUMSZ 5000 mm^2, B-cell 1 cell/uL, male, R/R, non-Asian,
    # single-agent.

    # ----- acMMAE structural parameters (Lu 2019 Table 1, theta1-theta11) -----
    lkdes      <- log(0.0046);  label("Rate constant of CL_TIME exponential decay (kdes, 1/hour)")     # Lu 2019 Table 1, theta1
    lcl_time   <- log(0.00623); label("Initial CL_TIME at time 0 for the reference subject (CL_TIME, L/hour)") # Lu 2019 Table 1, theta2
    lcl        <- log(0.0344);  label("acMMAE nonspecific linear clearance after repeated dosing (CL_SS, L/hour)") # Lu 2019 Table 1, theta3
    lvc        <- log(3.15);    label("acMMAE central volume (Vc, L)")                                 # Lu 2019 Table 1, theta4
    lvp        <- log(3.98);    label("acMMAE peripheral volume (Vp, L)")                              # Lu 2019 Table 1, theta5
    lq         <- log(0.0145);  label("acMMAE intercompartmental clearance (Q, L/hour)")               # Lu 2019 Table 1, theta6
    lvmax      <- log(0.0203);  label("acMMAE Michaelis-Menten maximum elimination rate (Vmax, ng/mL/hour)") # Lu 2019 Table 1, theta7
    lkm_ac     <- log(0.604);   label("acMMAE Michaelis-Menten constant (KM, ng/mL)")                  # Lu 2019 Table 1, theta8
    clss_emax  <- 0.223;        label("Maximum fractional effect of cycle on CL_NS (CLSSEMAX, unitless)") # Lu 2019 Table 1, theta9
    lt50_mo    <- log(3.53);    label("Time of half-maximal cycle effect on CL_NS (T50, months)")      # Lu 2019 Table 1, theta10 (converted to hours in model() via T50_hr = T50_mo * 24 * 30)
    gamma_ns   <- 2.27;         label("Sigmoidicity of the CL_NS(t) Hill function (gamma, unitless)")  # Lu 2019 Table 1, theta11

    # ----- Covariate effects on acMMAE PK parameters (Lu 2019 Table 2, theta22-theta37) -----
    e_wt_cl              <-  0.73;     label("Power exponent of WT on CL_SS (unitless)")                                                     # Lu 2019 Table 2, theta22
    e_wt_vc              <-  0.50;     label("Shared power exponent of WT on Vc, Vp, Q (unitless)")                                          # Lu 2019 Table 2, theta23
    e_sexf_vc            <-  1 / 1.20; label("Multiplicative effect of female sex on Vc (ratio female:male, unitless)")                      # Lu 2019 Table 2, theta24 inverted
    e_asian_vc           <-  0.929;    label("Multiplicative effect of Asian race on Vc (unitless; 7.1% lower V1 in Asian patients)")        # Lu 2019 Table 2, theta25
    e_line1l_vc          <-  1.20;     label("Multiplicative effect of treatment-naive (first-line) status on Vc (unitless)")                # Lu 2019 Table 2, theta26
    e_sexf_cl            <-  1 / 1.10; label("Multiplicative effect of female sex on CL_SS (ratio female:male, unitless)")                   # Lu 2019 Table 2, theta27 inverted
    e_alb_cl             <- -0.247;    label("Power exponent of ALB on CL_SS (unitless; reference 35 g/L for the PK side)")                  # Lu 2019 Table 2, theta28
    e_conmed_ritux_cl    <-  0.844;    label("Multiplicative effect of rituximab combination (interpreted as the COMBO_RG indicator for the polatuzumab-only Lu 2017 cohort) on CL_SS (unitless)") # Lu 2019 Table 2, theta29
    e_blbcell_cl         <-  0.0212;   label("Power exponent of max(1, BLBCELL) on CL_SS (unitless; B-cell count in cells/uL floored at 1)") # Lu 2019 Table 2, theta30
    e_tumsz_cl           <-  0.0521;   label("Linear coefficient of (TUMSZ/5000 - 1) on CL_SS (unitless; PK-side reference 5000 mm^2 SPD)")  # Lu 2019 Table 2, theta31
    e_line1l_kdes        <-  3.38;     label("Multiplicative effect of treatment-naive status on kdes (unitless)")                           # Lu 2019 Table 2, theta32
    e_conmed_ritux_kdes  <-  0.932;    label("Multiplicative effect of rituximab combination on kdes (unitless)")                            # Lu 2019 Table 2, theta33
    e_line1l_cl_time     <-  3.53;     label("Multiplicative effect of treatment-naive status on CL_TIME (unitless)")                        # Lu 2019 Table 2, theta34
    tmbd50_cl_time       <-  1150;     label("Half-maximal-effect TUMSZ on CL_TIME (mm^2 SPD; effect = TUMSZ / (tmbd50_cl_time + TUMSZ))")   # Lu 2019 Table 2, theta35
    bcell_thr_cl_time    <-  121;      label("B-cell threshold below which BLBCELL has no effect on CL_TIME (cells/uL)")                     # Lu 2019 Table 2, theta36
    e_blbcell_cl_time    <-  0.578;    label("Power exponent of max(1, BLBCELL/threshold) on CL_TIME (unitless)")                            # Lu 2019 Table 2, theta37

    # ----- Inter-individual variability on acMMAE PK parameters (Lu 2019 Table S3) -----
    etalcl_time ~ 1.89                       # Lu 2019 Table S3, Omega11 (CV 138%)
    etalcl      ~ 0.0376                     # Lu 2019 Table S3, Omega22 (CV 19.5%)
    etalvc      ~ 0.0151                     # Lu 2019 Table S3, Omega33 (CV 12.3%)
    etalvp      ~ 0.107                      # Lu 2019 Table S3, Omega44 (CV 32.7%)
    etalq       ~ 0.0538                     # Lu 2019 Table S3, Omega55 (CV 23.2%)
    etalvmax    ~ 0.462                      # Lu 2019 Table S3, Omega66 (CV 67.9%)

    # ----- acMMAE residual error (Lu 2019 Table S3, Sigma11) -----
    propSd <- sqrt(0.0254); label("Proportional residual error on acMMAE Cc (fraction)") # Lu 2019 Table S3, Sigma11 = 0.0254

    # ============================================================
    # ===== TTE PD side: Lu 2017 final Model 4 =====
    # ============================================================
    # Eq. 3 of Lu 2017 (page 404):
    #   h(t)        = beta * Edrug(t)^beta * t^(beta - 1)
    #   Edrug(t)    = alpha * Ce(t)
    #   dCe(t)/dt   = k1e * C(t) - ke0 * Ce(t)
    #   ke0         = k1e   (a single rate constant; the per-cycle parametric
    #                        constraint reported by Lu 2017 Table 1 footnote a)
    # alpha and beta and k1e are estimated by NONMEM in the log domain
    # (Table 1 RSE% column); covariate effects are estimated in the normal
    # domain (Table 1 SE column).

    # ----- Lu 2017 TTE PD structural parameters (Table 1) -----
    # No IIV was estimated for the TTE PD parameters in Lu 2017 (the paper
    # reports a single point estimate per parameter and a sequential
    # PK-then-PD fitting strategy with no omega block on alpha / beta / k1e;
    # see Lu 2017 Discussion and Conclusion). The model file records this by
    # not adding `eta*` terms on these parameters (standing nlmixr2lib policy:
    # unreported IIV with structural values present -> typical-value only when
    # no fixed(0) variance is required).
    lalpha_haz <- log(2.26e-6); label("Drug effect parameter alpha (1/(hour*ng/mL); log domain)")                # Lu 2017 Table 1: alpha = 2.26e-6 (RSE 49.2%)
    lbeta_haz  <- log(1.37);    label("Weibull function shape parameter beta (unitless; log domain)")            # Lu 2017 Table 1: beta = 1.37 (RSE 15.1%)
    lk1e_haz   <- log(3.60e-4); label("Effect-compartment distribution rate constant k1e (1/hour; log domain; ke0 = k1e)") # Lu 2017 Table 1: k1e = 3.60e-4 (RSE 73.8%); ke0 = k1e

    # ----- Lu 2017 TTE PD covariate effects (Table 1, THETA(4)-THETA(15); normal domain) -----
    # All covariate effects on the log-hazard are entered as the paper's
    # printed SE-with-sign point estimates, not log-transformed; they
    # multiply the centered / log-normalized covariate inside model() via
    # `exp(sum_i e_i * X_i)`.
    e_age_haz             <- -2.55e-3;  label("Effect of (AGE - 65 years) on log-hazard (1/year)")                                 # Lu 2017 Table 1: THETA(4) = -2.55e-3 (SE 0.0120)
    e_wt_haz              <-  0.0219;   label("Effect of (WT - 80 kg) on log-hazard (1/kg)")                                       # Lu 2017 Table 1: THETA(5) = 0.0219 (SE 0.0111)
    e_sexf_haz            <-  0.296;    label("Effect of SEXF (female sex indicator) on log-hazard (unitless)")                    # Lu 2017 Table 1: THETA(6) = 0.296 (SE 0.373)
    e_bl_pn_gr1_haz       <- -0.222;    label("Effect of BL_PN_GR1 (active grade 1 PN at baseline) on log-hazard (unitless)")      # Lu 2017 Table 1: THETA(7) = -0.222 (SE 0.324)
    e_prior_radiation_haz <- -7.94e-3;  label("Effect of PRIOR_RADIATION on log-hazard (unitless)")                                # Lu 2017 Table 1: THETA(8) = -7.94e-3 (SE 0.319)
    e_prior_vinca_haz     <- -0.102;    label("Effect of PRIOR_VINCA (prior vinca alkaloid) on log-hazard (unitless)")             # Lu 2017 Table 1: THETA(9) = -0.102 (SE 0.469)
    e_prior_platin_haz    <-  0.159;    label("Effect of PRIOR_PLATIN (prior platinum-based treatment) on log-hazard (unitless)")  # Lu 2017 Table 1: THETA(10) = 0.159 (SE 0.345)
    e_conmed_ritux_haz    <- -0.577;    label("Effect of CONMED_RITUX (rituximab combination) on log-hazard (unitless)")           # Lu 2017 Table 1: THETA(11) = -0.577 (SE 0.325)
    e_tumtp_dlbcl_haz     <- -0.0697;   label("Effect of TUMTP_DLBCL on log-hazard relative to FL reference (unitless)")           # Lu 2017 Table 1: THETA(12) = -0.0697 (SE 0.365)
    e_tumtp_other_nhl_haz <-  0.688;    label("Effect of TUMTP_OTHER_NHL on log-hazard relative to FL reference (unitless)")       # Lu 2017 Table 1: THETA(13) = 0.688 (SE 0.758)
    e_tumsz_haz           <-  0.169;    label("Effect of log(TUMSZ / 3000) on log-hazard (unitless; TTE-side reference 3000 mm^2 SPD)") # Lu 2017 Table 1: THETA(14) = 0.169 (SE 0.178)
    e_alb_haz             <-  0.0582;   label("Effect of (ALB - 39 g/L) on log-hazard (1/(g/L); TTE-side reference 39 g/L)")       # Lu 2017 Table 1: THETA(15) = 0.0582 (SE 0.0362)

    # No residual error on the TTE side. The Lu 2017 NONMEM run uses
    # F_FLAG / LIKE estimation in which the likelihood IS the survival /
    # event-density function itself; no observation-error model is
    # specified for PN events. The model exposes the instantaneous hazard
    # `hazard`, the cumulative hazard `cumhaz`, and the survival probability
    # `surv` as derived outputs for forward simulation.
  })

  model({
    # ===== 1. Derived covariate terms (PK side; Lu 2019 supplement Notations) =====
    bcel_cl_time <- max(1, BLBCELL / bcell_thr_cl_time)  # BCEL  for CL_TIME (Lu 2019)
    bcel_cl      <- max(1, BLBCELL)                      # BCEL1 for CL_SS  (Lu 2019)

    # Aggregate covariate multipliers on acMMAE PK parameters (Lu 2019 supplement).
    cov_vc <- (WT / 75)^e_wt_vc *
              e_sexf_vc^SEXF *
              e_asian_vc^RACE_ASIAN *
              e_line1l_vc^LINE_1L

    cov_cl <- (WT / 75)^e_wt_cl *
              e_sexf_cl^SEXF *
              (ALB / 35)^e_alb_cl *
              e_conmed_ritux_cl^CONMED_RITUX *
              bcel_cl^e_blbcell_cl *
              (1 + e_tumsz_cl * (TUMSZ / 5000 - 1))

    cov_kdes <- e_line1l_kdes^LINE_1L *
                e_conmed_ritux_kdes^CONMED_RITUX

    cov_cl_time <- e_line1l_cl_time^LINE_1L *
                   (TUMSZ / (tmbd50_cl_time + TUMSZ)) *
                   bcel_cl_time^e_blbcell_cl_time

    # ===== 2. Individual PK parameters (Lu 2019 acMMAE side) =====
    kdes         <- exp(lkdes) * cov_kdes                      # 1/hour
    cl_time_init <- exp(lcl_time + etalcl_time) * cov_cl_time  # CL_TIME initial value at t = 0 (L/hour)
    cl_inf       <- exp(lcl + etalcl) * cov_cl                 # CL_SS, L/hour
    vc           <- exp(lvc + etalvc) * cov_vc                 # L
    vp           <- exp(lvp + etalvp) * (WT / 75)^e_wt_vc      # L
    q            <- exp(lq + etalq) * (WT / 75)^e_wt_vc        # L/hour
    vmax         <- exp(lvmax + etalvmax)                      # ng/mL/hour
    km_ac        <- exp(lkm_ac)                                # ng/mL
    t50_hr       <- exp(lt50_mo) * 24 * 30                     # T50 in hours (paper reports months)
    t50_hr_gam   <- t50_hr^gamma_ns

    # ===== 3. Time-dependent PK quantities (Lu 2019 $DES; acMMAE only) =====
    # Hill function for CL_NS(t): at t=0, CL_NS = cl_inf * (1 + clss_emax);
    # as t -> infinity, CL_NS -> cl_inf.
    tgam  <- time^gamma_ns
    cl_ns <- cl_inf * (1 + clss_emax * t50_hr_gam / (t50_hr_gam + tgam))

    # Exponential decay of CL_TIME: CLT = CL_TIME_init * exp(-kdes * t).
    cl_t <- cl_time_init * exp(-kdes * time)

    # ===== 4. acMMAE micro-constants =====
    k12 <- q / vc
    k21 <- q / vp
    k10 <- (cl_t + cl_ns) / vc  # Linear (non-MM) elimination of acMMAE

    # ===== 5. acMMAE ODE system (Lu 2019; sign-corrected vs the published
    #         supplement Equations panel per the NONMEM DADT(2) form) =====
    d/dt(central)     <- -k10 * central -
                          k12 * central +
                          k21 * peripheral1 -
                          vmax * central / (km_ac + central / vc)
    d/dt(peripheral1) <-  k12 * central -
                          k21 * peripheral1

    # acMMAE plasma concentration (Cc, ng/mL = ug/L)
    Cc <- central / vc

    # ===== 6. Lu 2017 effect compartment (Eq. 3) =====
    # dCe/dt = k1e * Cc - ke0 * Ce, with ke0 = k1e (single rate constant).
    # 'effect' is the rxode2 state; the bare-name 'ce' below is the
    # algebraic synonym used in the hazard equation.
    alpha_pn <- exp(lalpha_haz)
    beta_pn  <- exp(lbeta_haz)
    k1e      <- exp(lk1e_haz)

    d/dt(effect) <- k1e * Cc - k1e * effect

    # ===== 7. Lu 2017 base Weibull hazard (Eq. 3) =====
    # Edrug(t) = alpha * Ce(t); h_base(t) = beta * Edrug^beta * t^(beta - 1).
    # The DEL = 1e-6 small-time offset on `time` matches the standard TTE
    # idiom (see Zecchin_2016_survival.R and the NONMEM DEL convention) and
    # keeps t^(beta - 1) finite at t = 0 without affecting the integrated
    # cumulative hazard. At t = 0, effect(0) = 0 so Edrug(0) = 0 and the
    # base hazard is identically zero regardless of DEL.
    del      <- 1e-6
    edrug    <- alpha_pn * effect
    haz_base <- beta_pn * edrug^beta_pn * (time + del)^(beta_pn - 1)

    # ===== 8. Lu 2017 covariate hazard (Eq. 2) =====
    # hcov = exp(sum_i theta_i * x_i) with covariates centered or
    # log-normalized at the Lu 2017 references (age 65 y, WT 80 kg,
    # ALB 39 g/L, TUMSZ 3000 mm^2; binary indicators enter unscaled).
    haz_cov <- exp(e_age_haz             * (AGE - 65) +
                   e_wt_haz              * (WT - 80) +
                   e_sexf_haz            * SEXF +
                   e_bl_pn_gr1_haz       * BL_PN_GR1 +
                   e_prior_radiation_haz * PRIOR_RADIATION +
                   e_prior_vinca_haz     * PRIOR_VINCA +
                   e_prior_platin_haz    * PRIOR_PLATIN +
                   e_conmed_ritux_haz    * CONMED_RITUX +
                   e_tumtp_dlbcl_haz     * TUMTP_DLBCL +
                   e_tumtp_other_nhl_haz * TUMTP_OTHER_NHL +
                   e_tumsz_haz           * log(TUMSZ / 3000) +
                   e_alb_haz             * (ALB - 39))

    hazard <- haz_base * haz_cov

    # ===== 9. Cumulative hazard and survival probability =====
    d/dt(cumhaz) <- hazard
    cumhaz(0)    <- 0
    surv         <- exp(-cumhaz)

    # ===== 10. Observation and residual error (PK side only) =====
    Cc ~ prop(propSd)
  })
}
