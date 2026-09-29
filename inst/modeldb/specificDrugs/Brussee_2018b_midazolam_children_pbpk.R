Brussee_2018b_midazolam_children_pbpk <- function() {
  description <- paste(
    "PBPK (semi-physiological; well-stirred liver + Qgut gut wall) population",
    "PK model for midazolam and its primary metabolite 1-OH-midazolam in 264",
    "post-operative children 1-18 years of age after a single oral dose.",
    "Physiological gut-wall, portal-vein and liver compartments carry the",
    "first-pass and systemic CYP3A metabolism for both analytes, and feed",
    "empirical central plus two peripheral compartments for midazolam and a",
    "single central compartment for 1-OH-midazolam. Tissue volumes, organ",
    "blood flows, plasma albumin, hematocrit, intestinal surface area and the",
    "blood:plasma ratio are all derived inside the model from body weight,",
    "age and sex (height and body surface area are derived first); the only",
    "estimated quantities are the whole-organ intrinsic clearances of the gut",
    "wall and liver, the inter-compartmental clearances, and their body-weight",
    "power exponents. Distribution volumes are fixed from an adult analysis",
    "and scaled linearly with body weight. The model runs in molar units, so",
    "the fraction of midazolam metabolised to 1-OH-midazolam is 1 with no",
    "molecular-weight conversion."
  )
  reference <- paste(
    "Brussee JM, Yu H, Krekels EHJ, Palic S, Brill MJE, Barrett JS,",
    "Rostami-Hodjegan A, de Wildt SN, Knibbe CAJ (2018). Characterization of",
    "Intestinal and Hepatic CYP3A-Mediated Metabolism of Midazolam in Children",
    "Using a Physiological Population Pharmacokinetic Modelling Approach.",
    "Pharm Res 35(9):182. doi:10.1007/s11095-018-2458-6.",
    sep = " "
  )
  vignette <- "Brussee_2018b_midazolam_children_pbpk"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  compartmentData <- list(
    depot = list(analyte = "midazolam", units = "nmol", specimen = "administration site", verified = TRUE),
    gut = list(analyte = "midazolam", units = "nmol", specimen = "tissue", verified = TRUE),
    portal = list(analyte = "midazolam", units = "nmol", specimen = "whole blood", verified = TRUE),
    liver = list(analyte = "midazolam", units = "nmol", specimen = "tissue", verified = TRUE),
    central = list(analyte = "midazolam", units = "nmol", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "midazolam", units = "nmol", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "midazolam", units = "nmol", specimen = "tissue", verified = TRUE),
    gut_1ohm = list(analyte = "1-OH-midazolam", units = "nmol", specimen = "tissue", verified = TRUE),
    portal_1ohm = list(analyte = "1-OH-midazolam", units = "nmol", specimen = "whole blood", verified = TRUE),
    liver_1ohm = list(analyte = "1-OH-midazolam", units = "nmol", specimen = "tissue", verified = TRUE),
    central_1ohm = list(analyte = "1-OH-midazolam", units = "nmol", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed at baseline. WT carries the two estimated power-function",
        "covariate relationships of Brussee 2018 Table II, both centred on the",
        "cohort median of 16 kg used as the model's reference individual:",
        "CL_H,int = 527 * (WT/16)^0.472, CL_G,int = 5.08 * (WT/16)^0.807,",
        "CL_H,int,M = 235 * (WT/16)^0.651 and Q_cp1 = 14.9 * (WT/16)^0.92.",
        "WT also scales the four fixed distribution volumes linearly from the",
        "76 kg adult reference of Frechen 2013 (exponents k3 and k7 fixed to",
        "1). Finally WT enters the body-surface-area formulas (supplemental",
        "eqs. S3 / S4), which in turn drive liver volume, cardiac output and",
        "intestinal surface area. Note that 16 kg is the covariate-centring",
        "reference of Table II and is NOT the cohort median body weight, which",
        "Table SI reports as 27.4 kg."
      ),
      source_name = "WT"
    ),
    AGE = list(
      description = "Age at the time of dosing",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Age was tested as a covariate on the estimated clearances and was NOT",
        "retained (Brussee 2018 Results). It is nonetheless load-bearing as a",
        "system descriptor: age drives small-intestine volume (Table I,",
        "V_in = 0.0467*AGE + 0.0901 L), plasma albumin (eq. 3), cardiac output",
        "(eq. 6), hematocrit (Table I age bands) and, together with sex,",
        "height (supplemental eqs. S1 / S2) and hence body surface area.",
        "Must be supplied in years; the height polynomials are only valid over",
        "the studied 1-18 year range."
      ),
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "male (SEXF = 0)",
      notes = paste(
        "Sex was tested as a covariate on the estimated clearances and was NOT",
        "retained (Brussee 2018 Results). It enters the model only through the",
        "system physiology: the hepatic fraction of cardiac output (Table I,",
        "0.28 in girls versus 0.255 in boys), the 12-18 year hematocrit band",
        "(0.41 female versus 0.43 male) and the sex-specific height",
        "polynomials (supplemental eqs. S1 male / S2 female). The source data",
        "are reported as a male/female split (Table SI: 148 boys / 116 girls),",
        "so SEXF = 1 - SEXM relative to that reporting."
      ),
      source_name = "SEX"
    )
  )

  population <- list(
    species = "human (children)",
    n_subjects = 264L,
    n_studies = 1L,
    age_range = "1-18 years",
    age_median = "7 years",
    weight_range = "9.1-137.6 kg",
    weight_median = "27.4 kg",
    sex_female_pct = 43.9,
    disease_state = paste(
      "Generally healthy children (American Society of Anesthesiologists",
      "physical status class I or II) undergoing elective surgery, given",
      "midazolam pre-operatively as an oral suspension."
    ),
    dose_range = "3-15 mg midazolam oral suspension, single pre-operative dose (median 10 mg)",
    regions = "United States (Children's Hospital of Philadelphia, PA)",
    notes = paste(
      "865 plasma samples from 264 analysed patients (Table SI). Of 266",
      "enrolled, two 14-year-olds in the sparsely sampled group were excluded",
      "because their recorded body weight was < 12 kg. 31 patients were",
      "densely sampled (median 10 samples each, range 8-11, around 0.25, 0.5,",
      "1, 1.5, 2, 3, 4, 6, 8, 10 and 22 h post-dose) and 233 were sparsely",
      "sampled (median 2 samples each, range 1-3, mostly within the first 4 h).",
      "Both midazolam and 1-OH-midazolam were assayed. Measurements below the",
      "limit of quantification were discarded per the M6 method (4 midazolam,",
      "5 1-OH-midazolam). Measured plasma concentrations were converted to",
      "blood concentrations with the blood:plasma ratio (eq. 1) before fitting,",
      "so the model's predicted concentrations are whole-blood concentrations.",
      "The source PK data are those of Barrett 2013 (reference 17 of the",
      "paper)."
    )
  )

  ini({
    # ---- Absorption -------------------------------------------------------
    # ka could not be estimated from these data and was fixed; Fa was assumed
    # to be 1 for midazolam.
    lka <- fixed(log(4.16))
    label("Midazolam absorption rate constant (1/h)") # Brussee 2018 Table I (K_a = 4.16 1/h) and Model Development ('could not be estimated, and was therefore fixed at 4.16 h-1')
    fa <- fixed(1)
    label("Fraction of the oral dose absorbed from the gut lumen (unitless)") # Brussee 2018 Table I (F_a = 1, ref 29) and eq. 7

    # ---- Midazolam whole-organ intrinsic clearances -----------------------
    # Reference individual is 16 kg (Table II covariate centring), NOT the
    # 27.4 kg cohort median.
    lcl_int_h <- log(527.0)
    label("Midazolam whole-organ intrinsic hepatic clearance at WT = 16 kg (L/h)") # Brussee 2018 Table II (CL_H,int,16kg = 527.0 L/h, RSE 7%)
    e_wt_cl_int_h <- 0.472
    label("Body-weight power exponent on midazolam intrinsic hepatic clearance (unitless)") # Brussee 2018 Table II (k1 = 0.472, RSE 16%)
    lcl_int_g <- log(5.08)
    label("Midazolam whole-organ intrinsic gut wall clearance at WT = 16 kg (L/h)") # Brussee 2018 Table II (CL_G,int,16kg = 5.08 L/h, RSE 10%)
    e_wt_cl_int_g <- 0.807
    label("Body-weight power exponent on midazolam intrinsic gut wall clearance (unitless)") # Brussee 2018 Table II (k2 = 0.807, RSE 10%)

    # ---- Midazolam distribution -------------------------------------------
    # Volumes could not be estimated (oral dosing only) and were fixed from
    # the 76 kg healthy-adult analysis of Frechen 2013, scaled linearly.
    lvc <- fixed(log(20.4))
    label("Midazolam central volume of distribution at WT = 76 kg (L, blood basis)") # Brussee 2018 Table II (V_c,76kg = 20.4 L fix; from Frechen 2013)
    lvp <- fixed(log(55.2))
    label("Midazolam first peripheral volume of distribution at WT = 76 kg (L, blood basis)") # Brussee 2018 Table II (V_p1,76kg = 55.2 L fix; from Frechen 2013)
    lvp2 <- fixed(log(79.1))
    label("Midazolam second peripheral volume of distribution at WT = 76 kg (L, blood basis)") # Brussee 2018 Table II (V_p2,76kg = 79.1 L fix; from Frechen 2013)
    e_wt_vc <- fixed(1)
    label("Body-weight power exponent on all midazolam distribution volumes (unitless)") # Brussee 2018 Table II (k3 = 1 fix)

    lq <- log(14.9)
    label("Midazolam inter-compartmental clearance to peripheral1 at WT = 16 kg (L/h)") # Brussee 2018 Table II (Q_cp1 = 14.9 L/h, RSE 19%)
    e_wt_q <- 0.92
    label("Body-weight power exponent on inter-compartmental clearance to peripheral1 (unitless)") # Brussee 2018 Table II (k4 = 0.92, RSE 21%)
    lq2 <- log(7.5)
    label("Midazolam inter-compartmental clearance to peripheral2 (L/h)") # Brussee 2018 Table II (Q_cp2 = 7.5 L/h, RSE 10%; no covariate identified)

    # ---- 1-OH-midazolam ---------------------------------------------------
    f_m <- fixed(1)
    label("Molar fraction of midazolam metabolised to 1-OH-midazolam (unitless)") # Brussee 2018 Table II (f_M = 1 fix) and Structural Model ('assumed 100%')
    lcl_int_h_1ohm <- log(235.0)
    label("1-OH-midazolam whole-organ intrinsic hepatic clearance at WT = 16 kg (L/h)") # Brussee 2018 Table II (CL_H,int,M,16kg = 235.0 L/h, RSE 6%)
    e_wt_cl_int_h_1ohm <- 0.651
    label("Body-weight power exponent on 1-OH-midazolam intrinsic hepatic clearance (unitless)") # Brussee 2018 Table II (k5 = 0.651, RSE 9%)

    # The 1-OH-midazolam gut wall intrinsic clearance could not be estimated
    # independently ('due to model instability') and was estimated as a
    # multiple of the midazolam gut wall intrinsic clearance of the SAME
    # individual, so it inherits that parameter's IIV and its WT exponent.
    # Named on the pattern of the shipped `ratio_glubile_pgp` ("Ratio of
    # Tel-GLU to telmisartan biliary secretion clearance") in
    # Tsuchitani_2024_telmisartan_pbpk.R. Deliberately NOT `lclrat_1ohm`: that
    # canonical is the ratio of one elimination route to the parent's
    # UNCHANGED-elimination clearance for parallel losses from a single
    # compartment, which is a different quantity from this ratio of the
    # metabolite's own organ intrinsic clearance to the parent's.
    ratio_cl_int_g_1ohm <- 18.4
    label("Ratio of 1-OH-midazolam to midazolam whole-organ intrinsic gut wall clearance (unitless)") # Brussee 2018 Table II (k6 = 18.4, RSE 12%; CL_G,int,M,i = k6 * CL_G,int,i)

    lvc_1ohm <- fixed(log(65.7))
    label("1-OH-midazolam central volume of distribution at WT = 76 kg (L, blood basis)") # Brussee 2018 Table II (V_M,76kg = 65.7 L fix; from Frechen 2013)
    e_wt_vc_1ohm <- fixed(1)
    label("Body-weight power exponent on the 1-OH-midazolam distribution volume (unitless)") # Brussee 2018 Table II (k7 = 1 fix)

    # ---- Inter-individual variability -------------------------------------
    # Table II reports these as VARIANCES (table footnote: 'Inter-individual
    # and residual variability values are shown as variance estimates') on the
    # log scale, per eq. 16 CL_int,i = theta_TV * exp(eta_i).
    etalcl_int_h ~ 0.25 # Brussee 2018 Table II, row 'omega^2 CL_H,int' = 0.25 (RSE 13%, shrinkage 25%)
    etalcl_int_g ~ 1.20 # Brussee 2018 Table II, row 'omega^2 CL_G,int' = 1.20 (RSE 13%, shrinkage 13%)
    etalq ~ 1.05 # Brussee 2018 Table II, row 'omega^2 Q_cp1' = 1.05 (RSE 35%, shrinkage 42%)
    etalq2 ~ 1.06 # Brussee 2018 Table II, row 'omega^2 Q_cp2' = 1.06 (RSE 31%, shrinkage 46%)
    etalcl_int_h_1ohm ~ 0.13 # Brussee 2018 Table II, row 'omega^2 CL_H,int,M' = 0.13 (RSE 18%, shrinkage 31%)

    # ---- Residual unexplained variability ---------------------------------
    # Combined proportional + additive on the linear scale (eq. 17,
    # Y = Cpred*(1+eps1) + eps2). Table II values are variances, so the SDs
    # used here are their square roots. Additive terms are in nmol/L.
    propSd <- sqrt(0.166)
    label("Midazolam proportional residual error (fraction)") # Brussee 2018 Table II (proportional error variance 0.166, RSE 8%; SD = sqrt(0.166) = 0.4074)
    addSd <- fixed(sqrt(0.001))
    label("Midazolam additive residual error (nmol/L)") # Brussee 2018 Table II (additive error variance 0.001 fix; SD = sqrt(0.001) = 0.03162)
    propSd_1ohm <- sqrt(0.292)
    label("1-OH-midazolam proportional residual error (fraction)") # Brussee 2018 Table II (proportional error variance 0.292, RSE 11%; SD = sqrt(0.292) = 0.5404)
    addSd_1ohm <- sqrt(0.528)
    label("1-OH-midazolam additive residual error (nmol/L)") # Brussee 2018 Table II (additive error variance 0.528, RSE 10%; SD = sqrt(0.528) = 0.7266)
  })
  model({
    # =====================================================================
    # SYSTEM PHYSIOLOGY
    # Every quantity in this block is a fixed population value taken from
    # Brussee 2018 Table I or the supplemental material; "For all
    # physiological parameters in Table I, population values were used
    # without interindividual variability or uncertainty."
    # =====================================================================

    # ---- Height (cm) from age, by sex ------------------------------------
    # Simcyp polynomials reproduced in the supplemental material. Valid over
    # the studied 1-18 year range only.
    ht_male <- 1.76179e-5 * AGE^7 - 1.19874e-3 * AGE^6 + 0.0323848 * AGE^5 -
      0.444112 * AGE^4 + 3.2946 * AGE^3 - 13.2191 * AGE^2 +
      33.75 * AGE + 52.62152 # Brussee 2018 supplemental eq. S1
    ht_female <- -1.51027e-6 * AGE^8 + 1.21261e-4 * AGE^7 - 0.0040023 * AGE^6 +
      0.070179 * AGE^5 - 0.708233 * AGE^4 + 4.1872 * AGE^3 -
      14.3393 * AGE^2 + 33.84778 * AGE + 51.535477 # Brussee 2018 supplemental eq. S2
    ht <- SEXF * ht_female + (1 - SEXF) * ht_male

    # ---- Body surface area (m^2), weight-banded ---------------------------
    bsa_lo <- 0.007184 * ht^0.725 * WT^0.425 # Brussee 2018 supplemental eq. S3 (WT < 15 kg)
    bsa_hi <- 0.024265 * ht^0.3964 * WT^0.537 # Brussee 2018 supplemental eq. S4 (WT >= 15 kg)
    bsa <- (WT < 15) * bsa_lo + (WT >= 15) * bsa_hi

    # ---- Tissue volumes (L) ----------------------------------------------
    v_liv <- 0.722 * bsa^1.176 # Brussee 2018 eq. 4 / Table I (liver volume from BSA)
    v_pv <- 0.0052 # Brussee 2018 Table I (portal vein volume fixed at the adult 5.2 mL)
    v_gut <- 0.0467 * AGE + 0.0901 # Brussee 2018 eq. 5 / Table I (small-intestine volume from age)

    # ---- Organ blood flows (L/h) -----------------------------------------
    co <- bsa * (110 + 184.974 * (exp(-0.0378 * AGE) - exp(-0.24477 * AGE))) # Brussee 2018 eq. 6 / Table I (cardiac output)
    q_h <- (0.28 * SEXF + 0.255 * (1 - SEXF)) * co # Brussee 2018 Table I (Q_h = 0.28*CO female, 0.255*CO male)
    q_pv <- 0.75 * q_h # Brussee 2018 Table I (portal vein = 75% of hepatic blood flow)
    q_ha <- 0.25 * q_h # Brussee 2018 Table I (hepatic artery = 25% of hepatic blood flow)
    q_in <- 0.40 * q_h # Brussee 2018 Table I (small intestine)
    q_muc <- 0.80 * q_in # Brussee 2018 Table I (mucosa)
    q_villi <- 0.60 * q_muc # Brussee 2018 Table I (microvilli / villous blood flow)

    # ---- Plasma protein binding ------------------------------------------
    # Pediatric albumin from age (eq. 3), then the McNamara-Alcorn scaling of
    # the adult unbound fraction (eq. 2). See the vignette's Errata: Figure
    # S2A of the supplement plots slightly higher unbound fractions than these
    # printed equations reproduce; the printed equations are used here.
    alb_ped <- 1.1287 * log(AGE) + 33.746 # Brussee 2018 eq. 3 (pediatric plasma albumin, g/L; Johnson 2006)
    alb_adult <- 37.7 # Brussee 2018 Table I ([P]_adult = 37.7 g/L)
    fu_adult <- 0.0303 # Brussee 2018 Table I (midazolam adult fraction unbound in plasma)
    fu_adult_1ohm <- 0.106 # Brussee 2018 Table I (1-OH-midazolam adult fraction unbound in plasma)
    fu_p <- 1 / (1 + (1 - fu_adult) / fu_adult * alb_ped / alb_adult) # Brussee 2018 eq. 2 (midazolam fraction unbound in plasma)
    fu_p_1ohm <- 1 / (1 + (1 - fu_adult_1ohm) / fu_adult_1ohm * alb_ped / alb_adult) # Brussee 2018 eq. 2 (1-OH-midazolam fraction unbound in plasma)

    # ---- Hematocrit and blood:plasma ratio --------------------------------
    # Table I reports hematocrit in five age / sex bands. Band edges are taken
    # at the lower bound of each printed band; the printed 7-12 y and 12-18 y
    # bands overlap at exactly 12 years and age 12 is assigned to the older
    # (sex-split) band.
    hct <- 0.36 * (AGE < 3) +
      0.37 * (AGE >= 3) * (AGE < 7) +
      0.40 * (AGE >= 7) * (AGE < 12) +
      (AGE >= 12) * (0.41 * SEXF + 0.43 * (1 - SEXF)) # Brussee 2018 Table I (Hem bands 1-2y 0.36, 3-6y 0.37, 7-12y 0.40, 12-18y 0.41 female / 0.43 male)
    bp <- 1 + hct * (fu_p - 1) # Brussee 2018 eq. 1 / Table I (B:P = 1 + Hem*(fu*Kp - 1) with Kp = 1)
    bp_1ohm <- 1 + hct * (fu_p_1ohm - 1) # Brussee 2018 Table I (1-OH-midazolam B:P with Kp = 1)

    # Fraction unbound in BLOOD is what the well-stirred liver uses. Eq. 15
    # of the paper is algebraically identical to Q_h*E_H*(B:P) with
    # fu_b = fu_p / (B:P), which confirms this conversion.
    fu_b <- fu_p / bp
    fu_b_1ohm <- fu_p_1ohm / bp_1ohm

    # Fraction unbound in the gut wall, assumed 1 for both analytes.
    fu_g <- 1 # Brussee 2018 Table I (F_u,G = 1) and eq. 9 text

    # ---- Qgut hybrid flow -------------------------------------------------
    # Intestinal radius and length from BSA (eqs. 13-14, in metres), surface
    # area from eq. 12 (in m^2 despite the "dm2" column header of Table I;
    # the 0.66 m^2 adult cut-off confirms the scale).
    r_int <- 0.5 * (0.016 * bsa + 0.0159) # Brussee 2018 eq. 13 (intestinal radius, m)
    h_int <- 2.56 * bsa + 2.95 # Brussee 2018 eq. 14 (intestinal length, m)
    a_int_raw <- 2 * 3.141592653589793 * r_int * (r_int + h_int) # Brussee 2018 eq. 12 (intestinal surface area, m^2)
    a_int <- 0.66 + (a_int_raw - 0.66) * (a_int_raw < 0.66) # Brussee 2018 Model Development ('cut-off at a maximum value of the adult value of 0.66 m2')

    # P_eff,man = 4.4e-4 cm/s = 1.584 cm/h = 0.1584 dm/h; A in dm^2 is
    # 100 * A[m^2], so CL_perm[L/h] = 0.1584 * 100 * A[m^2] = 15.84 * A[m^2].
    cl_perm <- 15.84 * a_int # Brussee 2018 eq. 11 / Table I (CL_perm = P_eff,man * A, P_eff,man = 4.4e-4 cm/s)
    q_gut <- q_villi * cl_perm / (q_villi + cl_perm) # Brussee 2018 eq. 10 (Qgut model; Yang 2007)

    # =====================================================================
    # INDIVIDUAL PARAMETERS
    # =====================================================================
    ka <- exp(lka)
    cl_int_h <- exp(lcl_int_h + etalcl_int_h) * (WT / 16)^e_wt_cl_int_h
    cl_int_g <- exp(lcl_int_g + etalcl_int_g) * (WT / 16)^e_wt_cl_int_g
    cl_int_h_1ohm <- exp(lcl_int_h_1ohm + etalcl_int_h_1ohm) * (WT / 16)^e_wt_cl_int_h_1ohm
    # CL_G,int,M,i = k6 * CL_G,int,i -- the individual (post-eta) parent gut
    # wall intrinsic clearance, so the metabolite inherits its IIV.
    cl_int_g_1ohm <- ratio_cl_int_g_1ohm * cl_int_g

    q_cp1 <- exp(lq + etalq) * (WT / 16)^e_wt_q
    q_cp2 <- exp(lq2 + etalq2)

    vc <- exp(lvc) * (WT / 76)^e_wt_vc
    vp1 <- exp(lvp) * (WT / 76)^e_wt_vc
    vp2 <- exp(lvp2) * (WT / 76)^e_wt_vc
    vc_1ohm <- exp(lvc_1ohm) * (WT / 76)^e_wt_vc_1ohm

    # =====================================================================
    # DERIVED FIRST-PASS QUANTITIES
    # Not used by the ODEs (the compartments below generate the same
    # extraction dynamically) but reported so the vignette can reproduce
    # Figure 4 and eq. 15 directly.
    # =====================================================================
    eg <- fu_g * cl_int_g / (q_gut + fu_g * cl_int_g) # Brussee 2018 eq. 9 (gut wall extraction ratio, Qgut model)
    eh <- fu_b * cl_int_h / (q_h + fu_b * cl_int_h) # Brussee 2018 eq. 8 (hepatic extraction ratio, well-stirred model)
    fg <- 1 - eg # Brussee 2018 eq. 7 text (gut wall bioavailability)
    fh <- 1 - eh # Brussee 2018 eq. 7 text (hepatic bioavailability)
    ftotal <- fa * fg * fh # Brussee 2018 eq. 7 (total oral bioavailability)
    cl_plasma <- q_h * cl_int_h * fu_p / (q_h + fu_p * cl_int_h / bp) # Brussee 2018 eq. 15 (total plasma clearance)

    # =====================================================================
    # ODE SYSTEM (Brussee 2018 Fig. 1)
    # All states hold nmol; all concentrations are BLOOD concentrations in
    # nmol/L, which is the scale the model was fitted on (measured plasma
    # concentrations were converted to blood via eq. 1 before fitting).
    #
    # Blood-flow balance: the portal vein receives q_pv from the systemic
    # circulation and delivers q_pv to the liver; the liver additionally
    # receives q_ha from the systemic circulation and returns
    # q_h = q_pv + q_ha to it, so flow closes at the central compartment.
    #
    # Drug leaves the midazolam gut wall into the portal vein at the Qgut
    # hybrid flow q_gut (permeability-limited, eq. 10), whereas
    # 1-OH-midazolam leaves at the villous blood flow q_villi. That asymmetry
    # is what Fig. 1 prints: E_G has q_gut in its denominator while E_G,M has
    # q_vi. It is mechanistically sensible -- the metabolite is formed inside
    # the enterocyte and does not have to cross the apical membrane -- and the
    # paper reports no permeability term for the metabolite.
    # =====================================================================
    c_gut <- gut / v_gut
    c_pv <- portal / v_pv
    c_liv <- liver / v_liv
    c_b <- central / vc
    c_p1 <- peripheral1 / vp1
    c_p2 <- peripheral2 / vp2
    c_gut_m <- gut_1ohm / v_gut
    c_pv_m <- portal_1ohm / v_pv
    c_liv_m <- liver_1ohm / v_liv
    c_b_m <- central_1ohm / vc_1ohm

    # ---- Midazolam --------------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(gut) <- fa * ka * depot - q_gut * c_gut - fu_g * cl_int_g * c_gut
    d/dt(portal) <- q_gut * c_gut + q_pv * c_b - q_pv * c_pv
    d/dt(liver) <- q_pv * c_pv + q_ha * c_b - q_h * c_liv - fu_b * cl_int_h * c_liv
    d/dt(central) <- q_h * c_liv - q_pv * c_b - q_ha * c_b -
      q_cp1 * (c_b - c_p1) - q_cp2 * (c_b - c_p2)
    d/dt(peripheral1) <- q_cp1 * (c_b - c_p1)
    d/dt(peripheral2) <- q_cp2 * (c_b - c_p2)

    # ---- 1-OH-midazolam ---------------------------------------------------
    d/dt(gut_1ohm) <- f_m * fu_g * cl_int_g * c_gut -
      q_villi * c_gut_m - fu_g * cl_int_g_1ohm * c_gut_m
    d/dt(portal_1ohm) <- q_villi * c_gut_m + q_pv * c_b_m - q_pv * c_pv_m
    d/dt(liver_1ohm) <- q_pv * c_pv_m + q_ha * c_b_m +
      f_m * fu_b * cl_int_h * c_liv -
      q_h * c_liv_m - fu_b_1ohm * cl_int_h_1ohm * c_liv_m
    d/dt(central_1ohm) <- q_h * c_liv_m - q_pv * c_b_m - q_ha * c_b_m

    # ---- Observations -----------------------------------------------------
    # Whole-blood concentrations, matching the scale the residual error was
    # estimated on. Divide by bp / bp_1ohm to recover plasma concentrations.
    Cc <- c_b
    Cc_1ohm <- c_b_m

    Cc ~ add(addSd) + prop(propSd)
    Cc_1ohm ~ add(addSd_1ohm) + prop(propSd_1ohm)
  })
}
