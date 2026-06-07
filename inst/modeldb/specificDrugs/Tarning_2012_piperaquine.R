# Population pharmacokinetic model for oral piperaquine in pregnant and
# nonpregnant women with uncomplicated Plasmodium falciparum or mixed
# falciparum/vivax malaria on the Thai-Myanmar border
# (Tarning 2012, Antimicrob Agents Chemother 56(4):1997-2007;
# doi:10.1128/AAC.05756-11).

Tarning_2012_piperaquine <- function() {
  description <- paste(
    "Three-compartment population PK model for oral piperaquine in 24",
    "pregnant (second / third trimester) and 24 matched non-pregnant",
    "women with uncomplicated malaria treated with the fixed-dose oral",
    "dihydroartemisinin-piperaquine combination once daily for 3 days",
    "(Tarning 2012 AAC). Transit-compartment absorption with 5 fixed",
    "transit compartments (ktr = (n+1)/MTT with n=5); the drug-transit",
    "rate is set equal to the absorption rate from the last transit to",
    "central (single estimated ktr). F fixed at 1; CL/F and F carry",
    "proportional pregnancy effects (+45.0% on CL/F and +46.8% on F).",
    "IIV on CL/F (21.5% CV) and Vc/F (39.5% CV); between-occasion",
    "variability (BOV across 3 dose occasions) on MTT (45.8% CV) and F",
    "(56.3% CV) multiplexed by the OCC indicator. Additive residual on",
    "natural-log concentrations (sigma = 0.285), encoded as proportional",
    "residual on the linear-concentration scale per Kloprogge 2018",
    "lumefantrine precedent. Companion file Tarning_2012_dihydroartemisinin.R",
    "models the co-administered dihydroartemisinin arm.",
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
        "Pregnancy enters the piperaquine model as proportional effects",
        "on CL/F (+45.0%) and F (+46.8%): CL = CL_typ * (1 + e_preg_cl *",
        "PREG); F = F_typ * (1 + e_preg_f * PREG). Reference category 0",
        "= non-pregnant.",
        sep = " "
      ),
      source_name        = "PREG"
    ),
    OCC = list(
      description        = paste(
        "Integer-valued dose-occasion indicator for between-occasion",
        "variability on mean transit time (MTT) and relative",
        "bioavailability (F). Values 1, 2, 3 identify the three",
        "consecutive daily doses (0, 24, and 48 h post-treatment-start)."
      ),
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Decomposed inside model() into binary indicators oc1, oc2, oc3",
        "that multiplex the BOV etas on MTT and on F (lfdepot).",
        "Tarning 2012 Results 'Pharmacokinetics of piperaquine':",
        "'Allowing intraindividual variability between dose occasions in",
        "mean transit time (delta-OFV = -344) and in relative",
        "bioavailability (delta-OFV = -165) significantly improved the",
        "model fit.' The text additionally notes a time-dependent",
        "increase in F across the three occasions (101%, 130%, and 170%",
        "at doses 1, 2, and 3); this pattern is interpreted as an",
        "emergent feature of the BOV draws in the model-building cohort",
        "rather than a structural per-occasion fixed effect, because",
        "Table 2 lists F as a single fixed value of 1 with BOV 56.3% --",
        "see vignette Errata. Downstream users supplying multi-dose",
        "data should set OCC = 1 for samples between dose 1 and dose 2,",
        "OCC = 2 between dose 2 and dose 3, and OCC = 3 after dose 3.",
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
    weight_range    = "36-78 kg (Table 1; median 51 pregnant, 48 non-pregnant)",
    sex_female_pct  = 100,
    disease_state   = paste(
      "Uncomplicated Plasmodium falciparum malaria (or mixed P. falciparum +",
      "P. vivax) in pregnant women in their second or third trimester (GA",
      "13.1-33.4 weeks, median 25.3) and age- and parasitaemia-matched",
      "non-pregnant women. Admission parasitaemia 6,780 (96-177,000)",
      "parasites/uL in pregnant women and 9,670 (32-136,000) in non-pregnant",
      "women (Table 1).",
      sep = " "
    ),
    dose_range      = paste(
      "Holleypharm fixed-dose dihydroartemisinin-piperaquine tablet (40 mg",
      "dihydroartemisinin + 320 mg piperaquine tetraphosphate equivalent",
      "to 184 mg base per tablet); once daily by mouth for 3 days at 0,",
      "24, and 48 h. Per-cohort total dose of piperaquine base 30.8",
      "(28.3-33.8) mg/kg in pregnant and 29.0 (27.7-33.8) mg/kg in",
      "non-pregnant women (Table 1).",
      sep = " "
    ),
    regions         = "Thai-Myanmar border (Shoklo Malaria Research Unit, Wang Pha Clinic)",
    notes           = paste(
      "Demographics from Tarning 2012 Table 1. 749 piperaquine samples in",
      "the pregnant cohort and 740 in the non-pregnant cohort were",
      "modelled (frequent venipuncture in the first 72 h and sparser",
      "sampling on days 5, 7, 14, 21, 28, 35, 42, 49, 56, 63, 77, and 84;",
      "Methods 'Drug regimen and blood sampling'). Population estimates",
      "in Table 2 are reported with no covariate centering (piperaquine",
      "did not retain an allometric body-weight effect in the final",
      "covariate model; delta-OFV = +4.62 vs base without allometric and",
      "an R-matrix nonpositive-semidefinite warning -- Results",
      "'Pharmacokinetics of piperaquine').",
      sep = " "
    )
  )

  ini({
    # Structural population mean parameters come from Tarning 2012
    # Table 2 'Population estimate' column for piperaquine. The paper
    # reports back-transformed (linear-scale) values; log() is applied
    # here for the nlmixr2 internal scale. Disposition is a 3-compartment
    # model parameterised as CL/F, Vc/F, Q1/F, Vp1/F, Q2/F, Vp2/F. The
    # absorption arm is a fixed 5-transit-compartment chain (Methods +
    # Fig 1A) with the drug-transit rate set identical to the absorption
    # rate from the last transit to central (single estimated ktr).
    lcl <- log(60.2)   ; label("Apparent oral clearance CL/F (L/h, non-pregnant)")            # Tarning 2012 Table 2: CL/F = 60.2  (RSE 10.5%; 95% CI 49.6-74.2)
    lvc <- log(3070)   ; label("Apparent central volume of distribution Vc/F (L)")             # Tarning 2012 Table 2: Vc/F = 3,070 (RSE 12.5%; 95% CI 2,400-3,930)
    lq  <- log(427)    ; label("Apparent inter-compartmental clearance Q1/F (L/h)")            # Tarning 2012 Table 2: Q1/F = 427   (RSE 14.9%; 95% CI 324-575)
    lvp <- log(4440)   ; label("Apparent peripheral-1 volume of distribution Vp1/F (L)")       # Tarning 2012 Table 2: Vp1/F = 4,440 (RSE 17.8%; 95% CI 3,210-6,300)
    lq2 <- log(160)    ; label("Apparent inter-compartmental clearance Q2/F (L/h)")            # Tarning 2012 Table 2: Q2/F = 160   (RSE 10.7%; 95% CI 131-196)
    lvp2 <- log(31400) ; label("Apparent peripheral-2 volume of distribution Vp2/F (L)")       # Tarning 2012 Table 2: Vp2/F = 31,400 (RSE 9.55%; 95% CI 26,500-38,300)
    lmtt <- log(2.04)  ; label("Mean absorption transit time MTT (h)")                         # Tarning 2012 Table 2: MTT = 2.04  (RSE 5.47%; 95% CI 1.85-2.25)

    # Relative bioavailability anchored at 1 (structural, fixed by the
    # source paper; Table 2 reports F as '100 (fixed)'). All structural
    # variability in F is captured via the per-occasion BOV etas
    # (etaiov_fdepot_1..3, CV 56.3% per occasion) plus the proportional
    # pregnancy multiplier below.
    lfdepot <- fixed(log(1)) ; label("Relative bioavailability F (unitless, fixed)")           # Tarning 2012 Table 2: F = 1 (fixed)

    # Pregnancy proportional effects. Tarning 2012 Results
    # 'Pharmacokinetics of piperaquine': '45.0% and 46.8% increase in
    # piperaquine elimination clearance and relative bioavailability,
    # respectively, in pregnant patients compared to non-pregnant
    # patients.' Encoded as multiplicative (1 + e * PREG).
    e_preg_cl <- 0.450 ; label("Pregnancy effect on CL/F: CL_pregnant / CL_nonpregnant - 1 = +0.450") # Tarning 2012 Table 2: Pregnancy on CL = 45.0% (RSE 17.8%; 95% CI 33.8-56.1)
    e_preg_f  <- 0.468 ; label("Pregnancy effect on F: F_pregnant / F_nonpregnant - 1 = +0.468")     # Tarning 2012 Table 2: Pregnancy on F  = 46.8% (RSE 35.6%; 95% CI 18.2-86.0)

    # IIV. Tarning 2012 Table 2 reports %CV in the 'IIV/BOV' column
    # (footnote: 'Coefficient of variation (%CV) for interindividual
    # variability and between-occasion variability are calculated as
    # SQRT of [exp(estimated variance) - 1]'). The internal log-normal
    # variances are recovered by omega^2 = log((CV/100)^2 + 1):
    #   CL CV 21.5% -> log(0.215^2 + 1) = 0.04518
    #   Vc CV 39.5% -> log(0.395^2 + 1) = 0.14496
    # Only CL/F and Vc/F carry IIV in the source model; the other
    # structural parameters (Q1, Vp1, Q2, Vp2) have dashes in the IIV
    # columns of Table 2.
    etalcl ~ 0.04518  # Tarning 2012 Table 2: IIV on CL = 21.5% CV (RSE 27.0%; 95% CI 14.5-26.2)
    etalvc ~ 0.14496  # Tarning 2012 Table 2: IIV on Vc = 39.5% CV (RSE 35.7%; 95% CI 21.7-50.7)

    # Between-occasion variability (BOV) on MTT and F across 3 dose
    # occasions, multiplexed by the OCC indicator. NONMEM $OMEGA
    # BLOCK(1) SAME pattern (Aregbe 2012 alvespimycin, Bender 2009
    # pregabalin precedents): occasion-1 carries the estimated variance
    # and occasions 2/3 fix the same value.
    #   MTT CV 45.8% -> log(0.458^2 + 1) = 0.19050
    #   F   CV 56.3% -> log(0.563^2 + 1) = 0.27536
    etaiov_mtt_1 ~ 0.19050           # Tarning 2012 Table 2: BOV on MTT = 45.8% CV (RSE 24.5%; 95% CI 35.4-56.1)
    etaiov_mtt_2 ~ fix(0.19050)      # NONMEM $OMEGA BLOCK(1) SAME (occasion 2 shares occasion-1 variance)
    etaiov_mtt_3 ~ fix(0.19050)      # NONMEM $OMEGA BLOCK(1) SAME (occasion 3 shares occasion-1 variance)

    etaiov_fdepot_1 ~ 0.27536        # Tarning 2012 Table 2: BOV on F = 56.3% CV (RSE 25.6%; 95% CI 41.8-72.2)
    etaiov_fdepot_2 ~ fix(0.27536)   # NONMEM $OMEGA BLOCK(1) SAME (occasion 2 shares occasion-1 variance)
    etaiov_fdepot_3 ~ fix(0.27536)   # NONMEM $OMEGA BLOCK(1) SAME (occasion 3 shares occasion-1 variance)

    # Residual error. Tarning 2012 Methods 'Pharmacokinetic analysis':
    # 'The residual random variability was assumed to be additive since
    # data were transformed into their natural logarithms (i.e.,
    # essentially equivalent to an exponential error model on a linear
    # scale).' Table 2 reports the additive log-scale SD as sigma =
    # 0.285; encoded here as propSd following the Kloprogge 2018
    # lumefantrine precedent (additive-on-log ~= proportional on the
    # linear-concentration scale for small sigma; exact would be
    # ~lnorm(expSd) with the same magnitude).
    propSd <- 0.285 ; label("Proportional residual SD on linear concentration scale (= SD on log scale)") # Tarning 2012 Table 2: sigma = 0.285 (RSE 4.97%; 95% CI 0.255-0.314)
  })

  model({
    # 1. Decompose OCC into per-occasion binary indicators and build
    #    the multiplexed BOV etas for MTT and F (lfdepot).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)

    iov_mtt    <- oc1 * etaiov_mtt_1    + oc2 * etaiov_mtt_2    + oc3 * etaiov_mtt_3
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 + oc3 * etaiov_fdepot_3

    # 2. Individual structural parameters. Piperaquine final covariate
    #    model retained pregnancy as a proportional effect on CL/F and
    #    on F (Results 'Pharmacokinetics of piperaquine'); body weight
    #    was tested as an allometric scaler but did not improve fit and
    #    was NOT retained, so structural parameters are not scaled by
    #    WT.
    cl  <- exp(lcl + etalcl) * (1 + e_preg_cl * PREG)
    vc  <- exp(lvc + etalvc)
    q   <- exp(lq)
    vp  <- exp(lvp)
    q2  <- exp(lq2)
    vp2 <- exp(lvp2)
    mtt <- exp(lmtt + iov_mtt)

    # Transit absorption rate constant. Fig 1A: ktr = (n+1)/MTT with
    # n = 5 transit compartments fixed.
    ktr <- (5 + 1) / mtt

    # 3. Three-compartment disposition micro-constants.
    kel <- cl  / vc
    k12 <- q   / vc
    k21 <- q   / vp
    k13 <- q2  / vc
    k31 <- q2  / vp2

    # 4. ODE system: transit-compartment absorption (depot through 5
    #    transit compartments) feeds the central compartment; central
    #    exchanges with two peripheral compartments and eliminates
    #    first-order from central.
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(transit3)    <-  ktr * transit2 - ktr * transit3
    d/dt(transit4)    <-  ktr * transit3 - ktr * transit4
    d/dt(transit5)    <-  ktr * transit4 - ktr * transit5
    d/dt(central)     <-  ktr * transit5 - kel * central -
                          k12 * central + k21 * peripheral1 -
                          k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central  - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central  - k31 * peripheral2

    # 5. Relative bioavailability F applied to the depot compartment.
    #    Two multiplicative factors:
    #      a. Per-occasion BOV: exp(lfdepot + iov_fdepot), with lfdepot
    #         fixed at log(1) so the baseline F is 1.
    #      b. Proportional pregnancy effect: (1 + e_preg_f * PREG).
    f(depot) <- exp(lfdepot + iov_fdepot) * (1 + e_preg_f * PREG)

    # 6. Piperaquine plasma concentration. Dose units are mg, Vc is L,
    #    so central / vc gives mg/L = ug/mL. Multiply by 1000 to obtain
    #    ng/mL to match the units declared in metadata and the LC-MS/MS
    #    assay range (LLOQ 1.50 ng/mL).
    Cc <- central / vc * 1000

    Cc ~ prop(propSd)
  })
}
