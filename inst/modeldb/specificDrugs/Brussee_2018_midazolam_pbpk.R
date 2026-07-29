Brussee_2018_midazolam_pbpk <- function() {
  description <- paste(
    "PBPK (semi-physiological; well-stirred liver + Qgut gut wall) population PK",
    "model for midazolam and its primary metabolite 1-OH-midazolam in 37 preterm",
    "neonates (gestational age 26-34 weeks, body weight 0.770-2.030 kg at the",
    "time of dosing). Distinguishes first-pass CYP3A-mediated metabolism in the",
    "gut wall (Qgut model) and liver (well-stirred model) from systemic hepatic",
    "elimination of the metabolite. Tissue volumes (V_h, V_pv, V_gw) and hepatic",
    "blood flow Q_h are allometrically scaled from a term-neonate reference",
    "(Bjorkman 2005) by body weight with fixed exponents (1 for volumes, 0.75",
    "for flow); intestinal length scales as 2.736 * WT[g]^0.512 cm (Struijs",
    "2009) so the Qgut hybrid flow varies with body size. Supports oral",
    "administration (depot, full first-pass through gut wall and liver) and IV",
    "(dose directly to central; no first-pass)."
  )
  reference <- paste(
    "Brussee JM, Yu H, Krekels EHJ, de Roos B, Brill MJE, van den Anker JN,",
    "Rostami-Hodjegan A, de Wildt SN, Knibbe CAJ (2018). First-Pass CYP3A-",
    "Mediated Metabolism of Midazolam in the Gut Wall and Liver in Preterm",
    "Neonates. CPT Pharmacometrics Syst Pharmacol 7(6):374-383.",
    "doi:10.1002/psp4.12295.",
    sep = " "
  )
  vignette <- "Brussee_2018_midazolam_pbpk"
  units    <- list(time = "hour", dosing = "microgram", concentration = "microgram/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight at the time of dose administration",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-occasion body weight (the paper records body weight at the day",
        "of dosing, range 0.770-2.030 kg). WT scales the term-neonate reference",
        "physiology (Bjorkman 2005 Table 1, reproduced in Brussee 2018 Table 1):",
        "liver volume V_h = 0.120 * (WT/3.55) L, portal vein V_pv = 0.778 * V_h,",
        "gut wall V_gw = 0.050 * (WT/3.55) L, and hepatic blood flow",
        "Q_h = 13.2 * (WT/3.55)^0.75 L/h with reference WT = 3.55 kg. WT also",
        "enters the empirical intestinal-length formula h_intest =",
        "2.736 * WT[g]^0.512 cm (Struijs 2009) which feeds the Qgut hybrid",
        "flow used in the gut-wall extraction-ratio formula. WT must be",
        "supplied in kg in the data column; the intestinal-length formula",
        "internally converts to grams via WT * 1000."
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human (preterm neonates)",
    n_subjects     = 37L,
    n_studies      = 2L,
    age_range      = "Postnatal age 3-11 days at the time of dosing",
    weight_range   = "Body weight 0.770-2.030 kg at the time of dosing",
    weight_median  = "1.1 kg (typical preterm neonate)",
    ga_range       = "Gestational age at birth 26-34 weeks",
    pna_range      = "Postnatal age 3-11 days",
    birth_weight_range = "Birth weight 0.745-2.135 g",
    sex_female_pct = NA_real_,
    disease_state  = paste(
      "Preterm neonates admitted to the Sophia Children's Hospital neonatal",
      "intensive care unit (Rotterdam, The Netherlands). Drawn from the IV",
      "(de Wildt 2001 Clin Pharmacol Ther 70:525-531) and oral (de Wildt 2002",
      "Br J Clin Pharmacol 53:390-392) midazolam pharmacokinetic studies in",
      "the same NICU cohort."
    ),
    dose_range     = paste(
      "0.1 mg/kg midazolam orally via nasogastric tube (n = 13) or as a",
      "30-minute intravenous infusion (n = 25). Subjects still meeting the",
      "inclusion criteria after >= 72 h received a second 0.1 mg/kg dose via",
      "the other route in a crossover design."
    ),
    regions        = "The Netherlands",
    notes          = paste(
      "Pooled IV + oral midazolam concentration data from de Wildt 2001 and",
      "de Wildt 2002. Plasma midazolam and 1-OH-midazolam were sampled at",
      "0.5, 1, 2, 4, 6, 12, and 24 hours post-dose; observations below LLOQ",
      "(< 1% for midazolam, < 2% for 1-OH-midazolam of 329 and 153 total",
      "samples) were discarded. Brussee 2018 tested birth weight, gestational",
      "age, sex, mechanical respiratory support, body weight per occasion,",
      "postmenstrual age, and postnatal age as covariates on the intrinsic",
      "clearances and volumes; none reached statistical significance and the",
      "final model carries WT only through the physiological scaling of tissue",
      "volumes, flows, and intestinal length."
    )
  )

  ini({
    # Structural parameters from Brussee 2018 Table 2 final estimates.
    # Clearances in L/h, volumes in L (blood basis -- volumes were estimated
    # on blood concentrations per Methods, "the drug and metabolite molar
    # concentrations in blood were calculated to be able to be used in the
    # well-stirred model"). Intrinsic clearances are on unbound-blood basis.

    # Absorption rate constant. Fixed at 10 1/h: "could not be estimated and
    # was assumed to be 10 h^-1, which entails an expected maximum
    # concentration around 12.5 minutes". A sensitivity analysis confirmed
    # values 4.16-25 1/h gave the same intrinsic-clearance estimates.
    lka            <- fixed(log(10));   label("Midazolam absorption rate constant (1/h)")                                   # Brussee 2018 Methods (Eq. 4 paragraph; Table 1 K_a = 10 1/h)

    # Midazolam intrinsic hepatic clearance.
    lcl_int_h      <- log(6.7);         label("Midazolam intrinsic hepatic clearance (L/h, unbound-blood basis)")           # Brussee 2018 Table 2 (CL_H,int = 6.7 L/h, RSE 10%)

    # Midazolam intrinsic gut wall clearance. Table 2 reports 19.6 mL/h
    # (RSE 178%); 19.6 mL/h = 0.0196 L/h.
    lcl_int_g      <- log(0.0196);      label("Midazolam intrinsic gut wall clearance (L/h, unbound basis)")                # Brussee 2018 Table 2 (CL_G,int = 19.6 mL/h = 0.0196 L/h, RSE 178%)

    # Midazolam apparent volume of distribution (empirical central in blood).
    lvc            <- log(3.0);         label("Midazolam apparent volume of distribution (L, blood basis)")                 # Brussee 2018 Table 2 (V = 3.0 L, RSE 11%)

    # 1-OH-midazolam intrinsic hepatic clearance.
    lcl_int_h_1ohm <- log(8.9);         label("1-OH-midazolam intrinsic hepatic clearance (L/h, unbound-blood basis)")      # Brussee 2018 Table 2 (CL_H,int,M = 8.9 L/h, RSE 22%)

    # 1-OH-midazolam apparent volume of distribution.
    lvc_1ohm       <- log(2.7);         label("1-OH-midazolam apparent volume of distribution (L, blood basis)")            # Brussee 2018 Table 2 (V_M = 2.7 L, RSE 43%)

    # IIV. Brussee 2018 Table 2 reports omega^2 variances directly on the
    # log-eta scale (paper Eq. 9 lognormal individual-parameter model).
    # No IIV was reported on CL_G,int (Table 2 shows "-" for that row).
    etalcl_int_h       ~ 0.887                                                                                              # Brussee 2018 Table 2 (omega^2_CL_H,int = 0.887, RSE 26%)
    etalvc             ~ 0.603                                                                                              # Brussee 2018 Table 2 (omega^2_V = 0.603, RSE 26%)
    etalcl_int_h_1ohm  ~ 0.832                                                                                              # Brussee 2018 Table 2 (omega^2_CL_H,int,M = 0.832, RSE 42%)
    etalvc_1ohm        ~ 0.887                                                                                              # Brussee 2018 Table 2 (omega^2_V_M = 0.887, RSE 48%)

    # Residual error -- combined additive + proportional on the linear scale
    # (Brussee 2018 Eq. 10: Y = Cpred*(1+E1) + E2 with E1, E2 ~ N(0, sigma^2)).
    # propSd and addSd are the standard deviations (SD = sqrt(reported
    # variance)). The additive components are FIX in NONMEM at sigma^2 =
    # 0.0001 ("The additive errors were fixed to very small numbers");
    # encoded with fixed() to preserve provenance.
    propSd      <- sqrt(0.201);    label("Midazolam proportional residual error (SD, fraction)")             # Brussee 2018 Table 2 (sigma^2_prop = 0.201, RSE 26%; SD = sqrt(0.201) = 0.4483)
    addSd       <- fixed(sqrt(0.0001));    label("Midazolam additive residual error (SD, microgram/L)")     # Brussee 2018 Table 2 (sigma^2_add = 0.0001 FIX; SD = 0.01)
    propSd_1ohm <- sqrt(0.164);    label("1-OH-midazolam proportional residual error (SD, fraction)")        # Brussee 2018 Table 2 (sigma^2_prop = 0.164, RSE 91%; SD = sqrt(0.164) = 0.4050)
    addSd_1ohm  <- fixed(sqrt(0.0001));    label("1-OH-midazolam additive residual error (SD, microgram/L)")# Brussee 2018 Table 2 (sigma^2_add = 0.0001 FIX; SD = 0.01)
  })
  model({
    # ===== Reference physiology (term-neonate baseline; Bjorkman 2005) =====
    # Reference body weight. Bjorkman 2005 reported volumes / flows for a
    # 3.55 kg term neonate; Brussee 2018 Table 1 carries that reference
    # forward and scales by WT / 3.55.
    wt_ref <- 3.55

    # Allometric tissue volumes (exponent 1 -> linear in WT).
    v_h  <- 0.120 * (WT / wt_ref)       # Liver volume (L); Brussee 2018 Table 1 (V_h,3.55kg = 0.120)
    v_pv <- 0.778 * v_h                 # Portal vein volume (L); Brussee 2018 Table 1 (V_pv = 0.778 * V_h, Aguirre-Reyes 2015)
    v_gw <- 0.050 * (WT / wt_ref)       # Gut wall volume (L); Brussee 2018 Table 1 (V_gw,3.55kg = 0.050)

    # Hepatic blood flow (allometric exponent 0.75).
    q_h    <- 13.2 * (WT / wt_ref)^0.75 # Hepatic blood flow (L/h); Brussee 2018 Table 1 (Q_h,3.55kg = 13.2)

    # Other splanchnic blood flows proportional to Q_h (Brussee 2018 Table 1).
    q_pv    <- 0.75 * q_h               # Portal vein flow to liver (L/h)
    q_ha    <- 0.25 * q_h               # Hepatic artery flow (L/h)
    q_in    <- 0.40 * q_h               # Small intestine flow (L/h)
    q_muc   <- 0.80 * q_in              # Mucosal flow (L/h)
    q_villi <- 0.60 * q_muc             # Villous blood flow (L/h; used in Qgut)

    # Intestinal surface area for the Qgut model. h_intest = 2.736 *
    # WT[g]^0.512 cm (Struijs 2009 ref 27) and r_intest = 1 cm (Ives 2016
    # ref 26). A = 2 * pi * r * h, output in cm^2.
    h_intest <- 2.736 * (WT * 1000)^0.512
    a_intest <- 2 * 3.141592653589793 * 1.0 * h_intest

    # Permeability term in the Qgut model. Effective intestinal permeability
    # for midazolam P_eff = 4.4e-4 cm/s = 1.584 cm/h (Yang 2007 ref 25).
    # CL_perm = P_eff * A (cm/h * cm^2 = cm^3/h) converted to L/h via /1000.
    p_eff   <- 1.584
    cl_perm <- p_eff * a_intest / 1000

    # Qgut hybrid flow (Brussee 2018 Eq. 3; Yang 2007).
    q_gut <- (q_villi * cl_perm) / (q_villi + cl_perm)

    # ===== Drug physicochemical constants (Brussee 2018 Table 1) =====
    # Midazolam fraction unbound in blood (Table 1 reports f_u,B = 0.04094
    # derived via the McNamara-Alcorn 2002 formula from f_u,adult = 0.0303,
    # [P]_pediatric = 27.1 g/L, [P]_adult = 37.0 g/L).
    fu_b      <- 0.04094

    # Midazolam fraction unbound in gut (assumed 1; Yang 2007 ref 25).
    fu_g      <- 1.0

    # Midazolam blood:plasma ratio (Maharaj 2013 ref 30 formula at
    # hematocrit 0.45 and K_p = 1; Table 1).
    bp        <- 0.568

    # 1-OH-midazolam fraction unbound in blood (McNamara-Alcorn from
    # f_u,M,adult = 0.106, Mandema 1992 ref 29).
    fu_b_1ohm <- 0.1394

    # 1-OH-midazolam blood:plasma ratio (Maharaj 2013).
    bp_1ohm   <- 0.613

    # Fraction of midazolam metabolised to 1-OH-midazolam (Brussee 2018
    # Table 1; assumed 1).
    f_m       <- 1.0

    # ===== Individual PK parameters =====
    ka              <- exp(lka)
    cl_int_h        <- exp(lcl_int_h + etalcl_int_h)
    cl_int_g        <- exp(lcl_int_g)
    vc              <- exp(lvc + etalvc)
    cl_int_h_1ohm   <- exp(lcl_int_h_1ohm + etalcl_int_h_1ohm)
    vc_1ohm         <- exp(lvc_1ohm + etalvc_1ohm)

    # ===== Extraction ratios =====
    # Gut wall extraction (Brussee 2018 Eq. 2; Yang 2007 Qgut model).
    eg <- (fu_g * cl_int_g) / (fu_g * cl_int_g + q_gut)

    # Hepatic extraction of parent (Brussee 2018 Eq. 1; well-stirred model).
    eh <- (fu_b * cl_int_h) / (fu_b * cl_int_h + q_h)

    # Hepatic extraction of 1-OH-midazolam (Brussee 2018; same well-stirred
    # formulation applied to the metabolite).
    eh_1ohm <- (fu_b_1ohm * cl_int_h_1ohm) /
               (fu_b_1ohm * cl_int_h_1ohm + q_h)

    # Bioavailabilities (Brussee 2018 Eq. 5 with F_a = 1, Gorski 1994 ref 10).
    fg <- 1 - eg
    fh <- 1 - eh

    # Hepatic blood clearances via well-stirred (Brussee 2018; derived from
    # Q_h and the extraction ratio). Multiplied by central blood concentration
    # below to give mass-per-time elimination rate.
    cl_h_blood      <- q_h * eh
    cl_h_blood_1ohm <- q_h * eh_1ohm

    # ===== ODE system =====
    # All compartments hold mass in dosing units (micrograms). Concentrations
    # in central / central_1ohm are computed as mass / volume = ug/L, which
    # equals ng/mL. vc and vc_1ohm are in blood-concentration terms; observed
    # plasma concentrations are obtained by dividing by the blood:plasma
    # ratio (since B:P = C_blood / C_plasma, so C_plasma = C_blood / B:P).
    #
    # Mass-balance summary for the lumped semi-PBPK representation:
    #
    #   Oral dose enters depot; first-order absorption ka * depot leaves
    #   the lumen. Gut wall and liver are treated as steady-state extractors
    #   (their physiological volumes are sub-millilitre; QSS is reached
    #   essentially instantaneously relative to the central-compartment
    #   timescale of hours). The fraction (1 - eg)(1 - eh) = fg * fh of the
    #   absorbed flux reaches the central blood pool; the complementary
    #   1 - fg * fh is converted to 1-OH-midazolam during the gut-wall +
    #   hepatic first pass. IV doses (event-table cmt = central) bypass
    #   the first pass.
    #
    #   Systemic recirculation of parent through the liver is handled by the
    #   well-stirred plasma-clearance term CL_h_blood * central / vc, every
    #   unit of which is converted to 1-OH-midazolam (f_m = 1).
    #
    #   The metabolite enters its own apparent-volume central pool and is
    #   cleared by hepatic CL_int,H,M via the same well-stirred formulation.

    d/dt(depot)        <- -ka * depot

    d/dt(central)      <- ka * depot * fg * fh -
                          cl_h_blood * central / vc

    d/dt(central_1ohm) <- ka * depot * (1 - fg * fh) * f_m +
                          cl_h_blood * central / vc * f_m -
                          cl_h_blood_1ohm * central_1ohm / vc_1ohm

    # ===== Observations =====
    # Convert blood -> plasma via the blood:plasma ratio.
    Cc      <- central      / vc      / bp
    Cc_1ohm <- central_1ohm / vc_1ohm / bp_1ohm

    Cc      ~ add(addSd)      + prop(propSd)
    Cc_1ohm ~ add(addSd_1ohm) + prop(propSd_1ohm)
  })
}
