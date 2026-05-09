Valitalo_2017_ketorolac <- function() {
  description <- "Three-compartment population PK model for IV ketorolac in adults, jointly fit to R-ketorolac and S-ketorolac plasma concentrations after racemic IV dosing in women at delivery, postpartum women, nonpregnant women, and men (Valitalo 2017 BJCP). Body-weight allometric scaling on clearance and volumes (reference 71 kg) plus proportional pregnancy-at-delivery and male-sex effects on clearance (and pregnancy-at-delivery on volumes), shared between enantiomers."
  reference <- paste(
    "Valitalo PA, Kemppainen H, Kulo A, Smits A, van Calsteren K, Olkkola KT,",
    "de Hoon J, Knibbe CAJ, Allegaert K. (2017).",
    "Body weight, gender and pregnancy affect enantiomer-specific ketorolac",
    "pharmacokinetics. Br J Clin Pharmacol 83(9):1966-1975.",
    "doi:10.1111/bcp.13311.",
    sep = " "
  )
  vignette <- "Valitalo_2017_ketorolac"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used for allometric scaling of CL (exponent 0.536) and all three volumes (V1, V2, V3; shared exponent 0.807) for both enantiomers, with reference weight 71 kg (Valitalo 2017 Table 2 footnotes a-d). The allometric exponents were estimated, not fixed.",
      source_name        = "WT"
    ),
    PREG = list(
      description        = "Indicator for woman currently in the at-delivery state (1 = at delivery / immediately post-Caesarean section, 0 = otherwise).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (nonpregnant; covers postpartum women, nonpregnant females, and males)",
      notes              = "Source NONMEM column is `WAD` (women at delivery) which carries the pregnancy-physiology effect at the time of Caesarean delivery; renamed to the canonical PREG per inst/references/covariate-columns.md. Postpartum women 4-5 months after delivery had reverted to PREG = 0. Multiplicative coefficient 0.554 on CL (+55%) and 0.273 on each of V1/V2/V3 (+27%) for both enantiomers (Valitalo 2017 Table 2 final model).",
      source_name        = "WAD"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Source NONMEM column is `MS` (male subject) with values inverted: SEXF = 1 - MS. The Valitalo 2017 typical-value parameters in Table 2 are reported for a 71-kg nonpregnant female (MS = 0, WAD = 0); to preserve those published structural values we apply the male-sex effect via (1 + e_sexf_cl * (1 - SEXF)) with e_sexf_cl = 0.363 (+36% CL in men vs women of the same body weight; Valitalo 2017 Table 2 final model). Effect applied identically to R- and S-enantiomer clearances.",
      source_name        = "MS"
    )
  )

  population <- list(
    n_subjects     = 67L,
    n_studies      = 2L,
    age_range      = "19-44 years",
    weight_range   = "40-106 kg",
    sex_female_pct = 82,
    race_ethnicity = NULL,
    disease_state  = paste(
      "Postoperative analgesia. Pooled analysis of two single-IV-bolus studies:",
      "(i) Leuven (Belgium): 41 women immediately post-Caesarean (at delivery),",
      "8 of those women restudied 4-5 months postpartum (paired session, same",
      "subjects), and 8 nonpregnant healthy female volunteers; (ii) Helsinki",
      "(Finland): 12 men and 6 nonpregnant women, ASA physical status I-II,",
      "undergoing minor eye surgery."
    ),
    dose_range     = paste(
      "Single IV bolus of racemic ketorolac tromethamine: 30 mg (Leuven; equivalent",
      "to 20.345 mg pure ketorolac per dose) or 0.5 mg/kg over 30 s (Helsinki).",
      "Each dose was modelled as half going to the R-enantiomer compartment and",
      "half to the S-enantiomer compartment (no interconversion assumed)."
    ),
    regions        = "Leuven, Belgium and Helsinki, Finland",
    notes          = paste(
      "Demographics from Valitalo 2017 Table 1: women at delivery n = 41, age",
      "33 (25-44) y, weight 73.9 (40-106) kg; postpartum n = 8 (subset of",
      "the 41 women restudied), age 31.5 (25-35) y, weight 60.8 (48.8-87.2) kg;",
      "nonpregnant females n = 14, age 31 (19-43) y, weight 61.2 (49-74.5) kg;",
      "males n = 12, age 32 (23-40) y, weight 74 (64-99) kg. The 8 postpartum",
      "sessions are repeats on 8 of the 41 women at delivery, so 41 + 14 + 12 =",
      "67 unique individuals (sex_female_pct = (41 + 14) / 67 = 82%).",
      "Race / ethnicity composition not reported. The model file encodes the",
      "Leuven-cohort residual-error magnitudes; see vignette Assumptions and",
      "deviations for the deviation from the per-study residual-error structure",
      "in the source NONMEM model (Valitalo 2017 Methods, equations 2-4)."
    )
  )

  ini({
    # R-ketorolac structural parameters (Valitalo 2017 Table 2 final model column,
    # typical values for a 71-kg nonpregnant female: WAD = 0, MS = 0).
    lcl_r  <- log(1.12)  ; label("R-ketorolac clearance for a 71-kg nonpregnant female (CL_R, L/h)")        # Valitalo 2017 Table 2 final, CL_R = 1.12 (RSE 5.95%)
    lvc_r  <- log(3.4)   ; label("R-ketorolac central volume for a 71-kg nonpregnant female (V1_R, L)")     # Valitalo 2017 Table 2 final, V1_R = 3.4 (RSE 5.99%)
    lq_r   <- log(3.86)  ; label("R-ketorolac shallow inter-compartmental clearance (Q1_R, L/h)")            # Valitalo 2017 Table 2 final, Q1_R = 3.86 (RSE 9.94%)
    lvp_r  <- log(2.29)  ; label("R-ketorolac shallow peripheral volume for a 71-kg nonpregnant female (V2_R, L)")  # Valitalo 2017 Table 2 final, V2_R = 2.29 (RSE 9.09%)
    lq2_r  <- log(0.58)  ; label("R-ketorolac deep inter-compartmental clearance (Q2_R, L/h)")               # Valitalo 2017 Table 2 final, Q2_R = 0.58 (RSE 11.8%)
    lvp2_r <- log(4.88)  ; label("R-ketorolac deep peripheral volume for a 71-kg nonpregnant female (V3_R, L)")  # Valitalo 2017 Table 2 final, V3_R = 4.88 (RSE 20.1%)

    # S-ketorolac structural parameters (Valitalo 2017 Table 2 final model column,
    # typical values for a 71-kg nonpregnant female).
    lcl_s  <- log(3.96)  ; label("S-ketorolac clearance for a 71-kg nonpregnant female (CL_S, L/h)")        # Valitalo 2017 Table 2 final, CL_S = 3.96 (RSE 5.7%)
    lvc_s  <- log(3.74)  ; label("S-ketorolac central volume for a 71-kg nonpregnant female (V1_S, L)")     # Valitalo 2017 Table 2 final, V1_S = 3.74 (RSE 6.95%)
    lq_s   <- log(20)    ; label("S-ketorolac shallow inter-compartmental clearance (Q1_S, L/h)")            # Valitalo 2017 Table 2 final, Q1_S = 20 (RSE 10.7%)
    lvp_s  <- log(3.86)  ; label("S-ketorolac shallow peripheral volume for a 71-kg nonpregnant female (V2_S, L)")  # Valitalo 2017 Table 2 final, V2_S = 3.86 (RSE 8.62%)
    lq2_s  <- log(1.3)   ; label("S-ketorolac deep inter-compartmental clearance (Q2_S, L/h)")               # Valitalo 2017 Table 2 final, Q2_S = 1.3 (RSE 20.5%)
    lvp2_s <- log(3.3)   ; label("S-ketorolac deep peripheral volume for a 71-kg nonpregnant female (V3_S, L)")  # Valitalo 2017 Table 2 final, V3_S = 3.3 (RSE 9.71%)

    # Allometric weight exponents (estimated; identical for both enantiomers).
    e_wt_cl    <- 0.536  ; label("Body-weight allometric exponent on CL (shared between R- and S-enantiomers)")  # Valitalo 2017 Table 2 final, WT_CL = 0.536 (RSE 36.5%)
    e_wt_vc_vp <- 0.807  ; label("Body-weight allometric exponent shared across V1, V2, V3 of both enantiomers")  # Valitalo 2017 Table 2 final, WT_V = 0.807 (RSE 20.1%); applied to vc, vp, vp2 of both enantiomers in model() per Table 2 footnotes b and d.

    # Categorical covariate effects (proportional, applied identically to both enantiomers).
    e_preg_cl    <- 0.554 ; label("Proportional CL increase for woman at delivery (vs nonpregnant); applied to both enantiomers")    # Valitalo 2017 Table 2 final, WAD_CL = 0.554 (RSE 20.4%)
    e_sexf_cl    <- 0.363 ; label("Proportional CL increase for males vs nonpregnant females of the same body weight; applied via (1 - SEXF) and shared between enantiomers")  # Valitalo 2017 Table 2 final, MS_CL = 0.363 (RSE 30.5%); source MS = 1 - SEXF.
    e_preg_vc_vp <- 0.273 ; label("Proportional volume increase (V1, V2, V3) for woman at delivery (vs nonpregnant); applied to both enantiomers")  # Valitalo 2017 Table 2 final, WAD_V = 0.273 (RSE 35.6%); applied uniformly to vc, vp, vp2 of both enantiomers in model().

    # IIV variances (omega^2; the paper reports omega = SD in Valitalo 2017 Table 2).
    # Independent etas per parameter; the V eta is shared across V1, V2, V3 of the
    # corresponding enantiomer per Methods: "Random effects were included for
    # clearances, and jointly for central and peripheral volumes of distribution
    # (one random effect affecting multiple parameters). Separate random effects
    # were estimated for each enantiomer."
    etalcl_r ~ 0.083521  # 0.289^2; Valitalo 2017 Table 2 final omega(CL_R) = 0.289
    etalvc_r ~ 0.051529  # 0.227^2; Valitalo 2017 Table 2 final omega(V_R) = 0.227 (shared across V1_R, V2_R, V3_R)
    etalcl_s ~ 0.053824  # 0.232^2; Valitalo 2017 Table 2 final omega(CL_S) = 0.232
    etalvc_s ~ 0.056644  # 0.238^2; Valitalo 2017 Table 2 final omega(V_S) = 0.238 (shared across V1_S, V2_S, V3_S)

    # Residual error: combined proportional + additive model derived from the
    # Valitalo 2017 paper's joint residual structure (equations 2-4):
    #   eps ~ N(0, sigma^2 * f^2 * f_prop + sigma^2 * (1 - f_prop)),
    # where f_prop = exp(theta_scale) / (1 + exp(theta_scale)) with theta_scale =
    # 5.1 from Table 2, giving f_prop ~= 0.9939 (~99.4% proportional, ~0.6%
    # additive). The paper estimated FOUR study-specific sigma values; we encode
    # the Leuven residual error here (the larger study, with all delivery /
    # postpartum data) and document the Helsinki values in the validation
    # vignette.
    # propSd = sigma * sqrt(f_prop), addSd = sigma * sqrt(1 - f_prop).
    propSd_r <- 0.0791  ; label("R-ketorolac proportional residual SD (fraction)")      # 0.0794 * sqrt(0.9939); Valitalo 2017 Table 2 sigma(Leuven, R) = 0.0794 (RSE 8.5%)
    addSd_r  <- 0.00618 ; label("R-ketorolac additive residual SD (mg/L)")              # 0.0794 * sqrt(0.00606)
    propSd_s <- 0.1097  ; label("S-ketorolac proportional residual SD (fraction)")      # 0.110 * sqrt(0.9939); Valitalo 2017 Table 2 sigma(Leuven, S) = 0.110 (RSE 20.5%)
    addSd_s  <- 0.00857 ; label("S-ketorolac additive residual SD (mg/L)")              # 0.110 * sqrt(0.00606)
  })

  model({
    # Reference body weight (median of the pooled population, Valitalo 2017
    # Table 2 footnotes a-d: "weight/71").
    wt_ref <- 71

    # Shared covariate factors. WT is allometric (power form); PREG and SEXF
    # apply linearly (1 + coef * indicator). SEXF has reference category 0
    # (male), so the male-vs-nonpregnant-female effect is applied via
    # (1 + e_sexf_cl * (1 - SEXF)), preserving the paper's published structural
    # values for a 71-kg nonpregnant female.
    cov_wt_cl    <- (WT / wt_ref)^e_wt_cl
    cov_wt_vc_vp <- (WT / wt_ref)^e_wt_vc_vp
    cov_cat_cl   <- (1 + e_preg_cl    * PREG) * (1 + e_sexf_cl * (1 - SEXF))
    cov_cat_v    <- (1 + e_preg_vc_vp * PREG)

    # R-ketorolac individual parameters. The single eta on volumes (etalvc_r)
    # is reused across vc_r, vp_r, vp2_r per Methods (joint volume IIV).
    cl_r  <- exp(lcl_r  + etalcl_r) * cov_wt_cl    * cov_cat_cl
    vc_r  <- exp(lvc_r  + etalvc_r) * cov_wt_vc_vp * cov_cat_v
    q_r   <- exp(lq_r)
    vp_r  <- exp(lvp_r  + etalvc_r) * cov_wt_vc_vp * cov_cat_v
    q2_r  <- exp(lq2_r)
    vp2_r <- exp(lvp2_r + etalvc_r) * cov_wt_vc_vp * cov_cat_v

    # S-ketorolac individual parameters (same covariate structure, separate etas).
    cl_s  <- exp(lcl_s  + etalcl_s) * cov_wt_cl    * cov_cat_cl
    vc_s  <- exp(lvc_s  + etalvc_s) * cov_wt_vc_vp * cov_cat_v
    q_s   <- exp(lq_s)
    vp_s  <- exp(lvp_s  + etalvc_s) * cov_wt_vc_vp * cov_cat_v
    q2_s  <- exp(lq2_s)
    vp2_s <- exp(lvp2_s + etalvc_s) * cov_wt_vc_vp * cov_cat_v

    # R-ketorolac micro-constants and three-compartment IV ODEs.
    kel_r <- cl_r / vc_r
    k12_r <- q_r  / vc_r
    k21_r <- q_r  / vp_r
    k13_r <- q2_r / vc_r
    k31_r <- q2_r / vp2_r

    d/dt(central_r)     <- -kel_r * central_r -
                            k12_r * central_r + k21_r * peripheral1_r -
                            k13_r * central_r + k31_r * peripheral2_r
    d/dt(peripheral1_r) <-  k12_r * central_r - k21_r * peripheral1_r
    d/dt(peripheral2_r) <-  k13_r * central_r - k31_r * peripheral2_r

    # S-ketorolac micro-constants and three-compartment IV ODEs.
    kel_s <- cl_s / vc_s
    k12_s <- q_s  / vc_s
    k21_s <- q_s  / vp_s
    k13_s <- q2_s / vc_s
    k31_s <- q2_s / vp2_s

    d/dt(central_s)     <- -kel_s * central_s -
                            k12_s * central_s + k21_s * peripheral1_s -
                            k13_s * central_s + k31_s * peripheral2_s
    d/dt(peripheral1_s) <-  k12_s * central_s - k21_s * peripheral1_s
    d/dt(peripheral2_s) <-  k13_s * central_s - k31_s * peripheral2_s

    # Plasma concentrations. Each enantiomer is dosed by half of the racemic
    # IV bolus into its own central compartment (no interconversion modelled);
    # the user supplies two dose records per administration with cmt =
    # central_r / central_s and amt = 0.5 * total racemic-ketorolac dose.
    Cc_r <- central_r / vc_r
    Cc_s <- central_s / vc_s

    Cc_r ~ prop(propSd_r) + add(addSd_r)
    Cc_s ~ prop(propSd_s) + add(addSd_s)
  })
}
