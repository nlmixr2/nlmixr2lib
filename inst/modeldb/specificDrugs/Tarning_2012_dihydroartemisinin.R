# Population pharmacokinetic model for oral dihydroartemisinin in
# pregnant and non-pregnant women with uncomplicated Plasmodium
# falciparum or mixed falciparum/vivax malaria on the Thai-Myanmar
# border (Tarning 2012, Antimicrob Agents Chemother 56(4):1997-2007;
# doi:10.1128/AAC.05756-11).

Tarning_2012_dihydroartemisinin <- function() {
  description <- paste(
    "One-compartment population PK model for oral dihydroartemisinin",
    "(parent drug, dosed as a fixed-dose tablet co-formulated with",
    "piperaquine) in 24 pregnant (second / third trimester) and 24",
    "matched non-pregnant women with uncomplicated malaria on the",
    "Thai-Myanmar border (Tarning 2012 AAC). Transit-compartment",
    "absorption with 7 fixed transit compartments (ktr = (n+1)/MTT",
    "with n=7); drug-transit rate is set equal to the absorption rate",
    "from the last transit to central (single estimated ktr).",
    "Allometric scaling of CL/F (exponent 3/4) and V/F (exponent 1) on",
    "body weight centered at the cohort median 48.5 kg. F fixed at 1",
    "with log-normal IIV (CV 30.3%); proportional pregnancy effect on F",
    "(-37.5%) and linear effect of log10 admission parasitaemia on F",
    "(+27.8% per log10 unit centered at 3.98). IIV on V/F (12.8% CV);",
    "between-occasion variability (BOV across 3 dose occasions) on MTT",
    "(50.9% CV) multiplexed by the OCC indicator. Additive residual on",
    "natural-log concentrations (sigma = 0.580), encoded as proportional",
    "residual on the linear-concentration scale per Kloprogge 2018",
    "lumefantrine precedent. Companion file Tarning_2012_piperaquine.R",
    "models the co-administered piperaquine arm.",
    sep = " "
  )
  reference <- paste(
    "Tarning J, Rijken MJ, McGready R, Phyo AP, Hanpithakpong W, Day NPJ,",
    "White NJ, Nosten F, Lindegardh N (2012). Population pharmacokinetics",
    "of dihydroartemisinin and piperaquine in pregnant and nonpregnant",
    "women with uncomplicated malaria.",
    "Antimicrobial Agents and Chemotherapy 56(4):1997-2007.",
    "doi:10.1128/AAC.05756-11.",
    sep = " "
  )
  vignette <- "Tarning_2012_dihydroartemisinin_piperaquine"
  units    <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject at admission. Tarning 2012 Table 1",
        "reports body-weight medians of 51 kg (range 36-58) in pregnant",
        "women and 48 kg (range 37-78) in non-pregnant women. Allometric",
        "scaling of CL/F (exponent 3/4) and V/F (exponent 1) centered at",
        "the model-building reference 48.5 kg (Table 4 footnote: 'Computed",
        "population mean values from NONMEM are calculated for a typical",
        "patient with a body weight of 48.5 kg and an initial logarithmic",
        "parasitaemia of 3.98'). Body-weight allometric scaling was",
        "retained for dihydroartemisinin (delta-OFV = -9.08; Results",
        "'Pharmacokinetics of dihydroartemisinin') but was NOT retained",
        "for the co-administered piperaquine.",
        sep = " "
      ),
      source_name        = "WT"
    ),
    PREG = list(
      description        = "Pregnancy status indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "1 = pregnant (second or third trimester), 0 = non-pregnant.",
        "Time-fixed per subject. Tarning 2012 enrolled 24 pregnant women",
        "(estimated gestational age 13.1-33.4 weeks, median 25.3) and 24",
        "matched non-pregnant women on the Thai-Myanmar border (Table 1).",
        "Pregnancy enters the dihydroartemisinin model as a proportional",
        "decrease in F: F = F_typ * (1 + e_preg_f * PREG) with",
        "e_preg_f = -0.375, i.e. F is 37.5% lower in pregnant than in",
        "non-pregnant women (Results 'Pharmacokinetics of",
        "dihydroartemisinin'). Reference category 0 = non-pregnant.",
        sep = " "
      ),
      source_name        = "PREG"
    ),
    PARA = list(
      description        = "Plasmodium falciparum parasitaemia at admission (asexual parasites/uL)",
      units              = "parasites/uL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Admission-only (time-fixed) parasitaemia. Tarning 2012 Table 1",
        "reports medians of 6,780 (range 96-177,000) parasites/uL in",
        "pregnant women and 9,670 (range 32-136,000) in non-pregnant",
        "women. Table 4 footnote centres the typical patient at",
        "'logarithmic parasitaemia 3.98', which matches log10 of the",
        "pooled-cohort median (log10(6,780) = 3.83 in pregnant women and",
        "log10(9,670) = 3.99 in non-pregnant women bracket 3.98); natural",
        "log would require parasitaemia ~ exp(3.98) ~= 54 parasites/uL,",
        "which is far below the cohort range, so the centering value is",
        "interpreted as log10. Applied as a linear effect on relative",
        "bioavailability F with the log10 transform performed inside",
        "model():",
        "F_para = 1 + e_para_f * (log10(max(PARA, 1)) - 3.98)",
        "with e_para_f = +0.278 per log10 unit (Results: '27.8% linear",
        "increase per unit logarithmic parasitemia'). The max(PARA, 1)",
        "gate keeps F_para finite for PARA <= 0; for downstream",
        "simulations where the absence of parasites should give no",
        "covariate effect, set PARA = 10^3.98 (~ 9,550) instead of 0.",
        "Same log10-inside-model() idiom as the Kloprogge 2014 / 2018",
        "Mahidol-Oxford malaria popPK family.",
        sep = " "
      ),
      source_name        = "PARA"
    ),
    OCC = list(
      description        = paste(
        "Integer-valued dose-occasion indicator for between-occasion",
        "variability on mean transit time (MTT). Values 1, 2, 3 identify",
        "the three consecutive daily doses (0, 24, and 48 h",
        "post-treatment-start)."
      ),
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Decomposed inside model() into binary indicators oc1, oc2, oc3",
        "that multiplex the BOV eta on MTT. Tarning 2012 Results",
        "'Pharmacokinetics of dihydroartemisinin': 'Allowing",
        "interindividual variability in the relative bioavailability",
        "(delta-OFV = -51.0) and intraindividual variability in mean",
        "transit time between dose occasions (delta-OFV = -344)",
        "significantly improved model fit.' Only MTT carries BOV in the",
        "dihydroartemisinin model; F carries IIV instead. Downstream",
        "users supplying multi-dose data should set OCC = 1 for samples",
        "between dose 1 and dose 2, OCC = 2 between dose 2 and dose 3,",
        "and OCC = 3 after dose 3.",
        sep = " "
      ),
      source_name        = "OCC"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 48L,
    n_pregnant      = 24L,
    n_studies       = 1L,
    age_range       = "18-45 years (Table 1; median 25 pregnant, 27.5 non-pregnant)",
    weight_range    = "36-78 kg (Table 1; median 51 pregnant, 48 non-pregnant; allometric reference 48.5 kg)",
    sex_female_pct  = 100,
    disease_state   = paste(
      "Uncomplicated Plasmodium falciparum malaria (or mixed P. falciparum +",
      "P. vivax) in pregnant women in their second or third trimester (GA",
      "13.1-33.4 weeks, median 25.3) and age- and parasitaemia-matched",
      "non-pregnant women. Admission parasitaemia 6,780 (96-177,000)",
      "parasites/uL in pregnant women and 9,670 (32-136,000) in",
      "non-pregnant women (Table 1).",
      sep = " "
    ),
    dose_range      = paste(
      "Holleypharm fixed-dose dihydroartemisinin-piperaquine tablet (40",
      "mg dihydroartemisinin + 320 mg piperaquine tetraphosphate",
      "equivalent to 184 mg base per tablet); once daily by mouth for 3",
      "days at 0, 24, and 48 h. Per-cohort total dihydroartemisinin dose",
      "6.67 (6.12-7.32) mg/kg in pregnant and 6.28 (6.00-7.32) mg/kg in",
      "non-pregnant women (Table 1).",
      sep = " "
    ),
    regions         = "Thai-Myanmar border (Shoklo Malaria Research Unit, Wang Pha Clinic)",
    notes           = paste(
      "Demographics from Tarning 2012 Table 1. 480 dihydroartemisinin",
      "plasma samples per cohort (frequent venipuncture in the first 72",
      "h; Methods 'Drug regimen and blood sampling'). The structural",
      "disposition model was a one-compartment fit (delta-OFV vs",
      "two-compartment = -0.096 with first-order absorption and not",
      "significantly improved with the final transit-absorption + IIV",
      "structure; Results 'Pharmacokinetics of dihydroartemisinin').",
      "Data below the limit of quantification (2.0 ng/mL) were coded as",
      "missing.",
      sep = " "
    )
  )

  ini({
    # Structural population mean parameters come from Tarning 2012
    # Table 4 'Population estimate' column for dihydroartemisinin. The
    # paper reports back-transformed (linear-scale) values; log() is
    # applied here for the nlmixr2 internal scale. Population estimates
    # are centered at a 48.5-kg patient with admission log10(PARA) =
    # 3.98 (Table 4 footnote). Disposition is a 1-compartment model
    # parameterised as CL/F and V/F. The absorption arm is a fixed
    # 7-transit-compartment chain (Methods + Fig 1B) with the
    # drug-transit rate set identical to the absorption rate from the
    # last transit to central (single estimated ktr).
    lcl  <- log(78.0)  ; label("Apparent oral clearance CL/F (L/h, typical 48.5-kg non-pregnant non-malaria-anchor)") # Tarning 2012 Table 4: CL/F = 78.0  (RSE 7.00%; 95% CI 67.4-88.6)
    lvc  <- log(129)   ; label("Apparent central volume of distribution V/F (L, typical 48.5-kg)")                    # Tarning 2012 Table 4: V/F  = 129   (RSE 7.30%; 95% CI 111-149)
    lmtt <- log(0.982) ; label("Mean absorption transit time MTT (h)")                                                # Tarning 2012 Table 4: MTT  = 0.982 (RSE 5.30%; 95% CI 0.890-1.09)

    # Relative bioavailability anchored at 1 (structural, fixed by the
    # source paper; Table 4 reports F as '100 (fixed)'). Log-normal IIV
    # (etalfdepot, CV 30.3%) plus the proportional pregnancy multiplier
    # and the linear log10-parasitaemia multiplier below capture all
    # F-related variability.
    lfdepot <- fixed(log(1)) ; label("Relative bioavailability F (unitless, fixed)")           # Tarning 2012 Table 4: F = 1 (fixed)

    # Allometric exponents. Fixed at the canonical Mahidol-Oxford
    # malaria-popPK values used by the source paper: 3/4 on clearance
    # (CL/F) and 1 on volume (V/F). Reported in Methods 'Pharmacokinetic
    # analysis': 'Body weight was also evaluated as a simultaneous
    # incorporation of an allometric function on all clearance and
    # volume parameters, where clearance values scale to mass to a
    # power of 0.75 and where volume values scale to mass to the first
    # power [e.g., individual clearance value = typical clearance value
    # x (individual body weight / median body weight in the
    # population)^0.75].' Retained in the final dihydroartemisinin
    # covariate model (delta-OFV = -9.08).
    allo_cl <- fixed(3/4) ; label("Allometric exponent on CL/F (unitless, fixed)")             # Tarning 2012 Methods + Results 'Pharmacokinetics of dihydroartemisinin'
    allo_vc <- fixed(1)   ; label("Allometric exponent on V/F (unitless, fixed)")              # Tarning 2012 Methods + Results 'Pharmacokinetics of dihydroartemisinin'

    # Pregnancy proportional effect on F. Tarning 2012 Results
    # 'Pharmacokinetics of dihydroartemisinin': 'a 37.5% lower relative
    # bioavailability in pregnant women than non-pregnant women.'
    # Encoded as multiplicative (1 + e * PREG) with e = -0.375 so the
    # pregnant cohort (PREG = 1) gets F * 0.625.
    e_preg_f <- -0.375 ; label("Pregnancy effect on F: F_pregnant / F_nonpregnant - 1 = -0.375") # Tarning 2012 Table 4: Pregnancy on F = -37.5% (RSE 17.2%; 95% CI 24.0-49.0)

    # log10 admission parasitaemia linear effect on F. Tarning 2012
    # Results 'Pharmacokinetics of dihydroartemisinin': 'an increase in
    # relative bioavailability with increasing initial parasitemia
    # (27.8% linear increase per unit logarithmic parasitemia)'.
    # Centered at the typical patient's log10(PARA) = 3.98 (Table 4
    # footnote). The covariate name e_para_f is used because the active
    # covariate column is the raw PARA count with the log10 transform
    # applied inside model() (same log10-inside-model() idiom as the
    # Kloprogge 2014 / 2018 Mahidol-Oxford precedents).
    e_para_f <- 0.278 ; label("Linear effect of log10 admission parasitaemia on F (per log10 parasites/uL, centered at 3.98)") # Tarning 2012 Table 4: Parasitemia on F = +27.8% (RSE 12.7%; 95% CI 21.3-37.3)

    # IIV. Tarning 2012 Table 4 reports %CV in the 'IIV/BOV' column
    # (footnote: 'Coefficient of variation (%CV) for interindividual
    # variability and between-occasion variability are calculated as
    # SQRT of [exp(estimated variance) - 1]'). The internal log-normal
    # variances are recovered by omega^2 = log((CV/100)^2 + 1):
    #   V CV 12.8% -> log(0.128^2 + 1) = 0.01625
    #   F CV 30.3% -> log(0.303^2 + 1) = 0.08786
    # Only V/F and F carry IIV in the source model; CL/F has a dash in
    # the IIV column of Table 4.
    etalvc     ~ 0.01625  # Tarning 2012 Table 4: IIV on V/F = 12.8% CV (RSE 27.7%; 95% CI 8.70-16.0)
    etalfdepot ~ 0.08786  # Tarning 2012 Table 4: IIV on F   = 30.3% CV (RSE 24.8%; 95% CI 20.9-35.7)

    # Between-occasion variability (BOV) on MTT across 3 dose
    # occasions, multiplexed by the OCC indicator. NONMEM $OMEGA
    # BLOCK(1) SAME pattern: occasion-1 carries the estimated variance
    # and occasions 2/3 fix the same value.
    #   MTT CV 50.9% -> log(0.509^2 + 1) = 0.23069
    etaiov_mtt_1 ~ 0.23069        # Tarning 2012 Table 4: BOV on MTT = 50.9% CV (RSE 15.7%; 95% CI 41.8-59.6)
    etaiov_mtt_2 ~ fix(0.23069)   # NONMEM $OMEGA BLOCK(1) SAME (occasion 2 shares occasion-1 variance)
    etaiov_mtt_3 ~ fix(0.23069)   # NONMEM $OMEGA BLOCK(1) SAME (occasion 3 shares occasion-1 variance)

    # Residual error. Tarning 2012 Methods 'Pharmacokinetic analysis':
    # 'The residual random variability was assumed to be additive since
    # data were transformed into their natural logarithms (i.e.,
    # essentially equivalent to an exponential error model on a linear
    # scale).' Table 4 reports the additive log-scale SD as sigma =
    # 0.580; encoded here as propSd following the Kloprogge 2018
    # lumefantrine precedent (additive-on-log ~= proportional on the
    # linear-concentration scale; exact would be ~lnorm(expSd) with the
    # same magnitude). Sigma is large (0.580 implies geometric SD
    # factor exp(0.580) ~= 1.79); see vignette Assumptions for the
    # log-normal vs proportional approximation note.
    propSd <- 0.580 ; label("Proportional residual SD on linear concentration scale (= SD on log scale)") # Tarning 2012 Table 4: sigma = 0.580 (RSE 5.11%; 95% CI 0.523-0.634)
  })

  model({
    # 1. Decompose OCC into per-occasion binary indicators and build
    #    the multiplexed BOV eta for MTT.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)

    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2 + oc3 * etaiov_mtt_3

    # 2. Individual structural parameters. Dihydroartemisinin final
    #    covariate model retained allometric body-weight scaling on
    #    CL/F (exponent 3/4) and V/F (exponent 1), plus a proportional
    #    pregnancy effect and a linear log10-parasitaemia effect on F.
    #    Allometric centering at the cohort reference 48.5 kg.
    cl  <- exp(lcl)              * (WT / 48.5)^allo_cl
    vc  <- exp(lvc + etalvc)     * (WT / 48.5)^allo_vc
    mtt <- exp(lmtt + iov_mtt)

    # Transit absorption rate constant. Fig 1B: ktr = (n+1)/MTT with
    # n = 7 transit compartments fixed.
    ktr <- (7 + 1) / mtt

    # 3. One-compartment disposition micro-constant.
    kel <- cl / vc

    # 4. ODE system: transit-compartment absorption (depot through 7
    #    transit compartments) feeds the central compartment;
    #    elimination first-order from central.
    d/dt(depot)    <- -ktr * depot
    d/dt(transit1) <-  ktr * depot    - ktr * transit1
    d/dt(transit2) <-  ktr * transit1 - ktr * transit2
    d/dt(transit3) <-  ktr * transit2 - ktr * transit3
    d/dt(transit4) <-  ktr * transit3 - ktr * transit4
    d/dt(transit5) <-  ktr * transit4 - ktr * transit5
    d/dt(transit6) <-  ktr * transit5 - ktr * transit6
    d/dt(transit7) <-  ktr * transit6 - ktr * transit7
    d/dt(central)  <-  ktr * transit7 - kel * central

    # 5. Relative bioavailability F applied to the depot compartment.
    #    Three multiplicative factors:
    #      a. Log-normal IIV: exp(lfdepot + etalfdepot), centered at F
    #         = 1.
    #      b. Proportional pregnancy effect: (1 + e_preg_f * PREG).
    #      c. log10-parasitaemia linear effect: 1 + e_para_f *
    #         (log10(max(PARA, 1)) - 3.98).
    #    The max(PARA, 1) gate keeps the log10 finite for PARA <= 0.
    fpara <- 1 + e_para_f * (log10(max(PARA, 1)) - 3.98)
    fpreg <- 1 + e_preg_f * PREG
    f(depot) <- exp(lfdepot + etalfdepot) * fpreg * fpara

    # 6. Dihydroartemisinin plasma concentration. Dose units are mg,
    #    Vc is L, so central / vc gives mg/L = ug/mL. Multiply by 1000
    #    to obtain ng/mL to match the units declared in metadata and
    #    the LC-MS/MS assay range (LLOQ 2.0 ng/mL).
    Cc <- central / vc * 1000

    Cc ~ prop(propSd)
  })
}
