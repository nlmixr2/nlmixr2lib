# Population pharmacokinetic model for oral lumefantrine in pregnant and
# non-pregnant women with uncomplicated Plasmodium falciparum malaria in
# Uganda after the standard fixed-dose artemether-lumefantrine regimen
# (Kloprogge 2013, CPT: PSP 2:e83; doi:10.1038/psp.2013.59).

Kloprogge_2013_lumefantrine <- function() {
  description <- paste(
    "Population PK model for oral lumefantrine in pregnant and non-pregnant",
    "women with uncomplicated Plasmodium falciparum malaria in Uganda after",
    "the standard fixed-dose oral artemether-lumefantrine treatment",
    "(Kloprogge 2013). Flexible five-compartment transit absorption chain",
    "into a two-compartment disposition model with relative bioavailability",
    "F1 fixed at 1, log-normal IIV on CL / Vp / MTT / F, and covariate",
    "effects of pregnancy on intercompartmental clearance (-36.5%,",
    "categorical) and body temperature on mean absorption transit time",
    "(+16.5% per degC over 36.0-39.8 degC, linear-deviation centered at",
    "the cohort median 36.9 degC).",
    sep = " "
  )
  reference <- paste(
    "Kloprogge F, Piola P, Dhorda M, Muwanga S, Turyakira E, Apinan S,",
    "Lindegardh N, Nosten F, Day NPJ, White NJ, Guerin PJ, Tarning J (2013).",
    "Population Pharmacokinetics of Lumefantrine in Pregnant and Nonpregnant",
    "Women With Uncomplicated Plasmodium falciparum Malaria in Uganda.",
    "CPT: Pharmacometrics & Systems Pharmacology 2:e83.",
    "doi:10.1038/psp.2013.59.",
    sep = " "
  )
  vignette <- "Kloprogge_2013_lumefantrine"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    PREG = list(
      description        = "Pregnancy status indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "1 = pregnant (second or third trimester), 0 = non-pregnant.",
        "Time-fixed per subject. Kloprogge 2013 enrolled 116 pregnant women",
        "(26 with dense venous sampling, 90 with sparse capillary sampling)",
        "and 17 non-pregnant post-partum controls with dense venous",
        "sampling. The paper applies pregnancy as a categorical relative",
        "difference on intercompartmental clearance Q/F:",
        "Q_indiv = TVQ * (1 + e_preg_q * PREG) with e_preg_q = -0.365,",
        "so pregnant women have ~36.5% lower Q/F than non-pregnant women",
        "(reference category 0 = non-pregnant). Sign convention matches",
        "Table 2 'Pregnancy on Q'.",
        sep = " "
      ),
      source_name        = "PREG"
    ),
    BODYTEMP = list(
      description        = "Body temperature at admission",
      units              = "degC",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Admission body temperature, time-fixed per subject. Population",
        "median 36.9 degC (Table 1, all-cohort column); range 36.0-39.8",
        "degC. Kloprogge 2013 applies a linear-deviation effect on the",
        "mean absorption transit time:",
        "MTT_indiv = TVMTT * (1 + e_bodytemp_mtt * (BODYTEMP - 36.9))",
        "with e_bodytemp_mtt = +0.165 per degC. Per the paper's",
        "Discussion, higher fever is plausibly associated with reduced",
        "gut motility and prolonged lumefantrine absorption; the effect",
        "is significant within the observed 36.0-39.8 degC range and",
        "should not be extrapolated outside that window.",
        sep = " "
      ),
      source_name        = "TEMP"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 132L,
    n_studies       = 1L,
    n_pregnant      = 115L,
    n_nonpregnant   = 17L,
    age_range       = "15.0-38.0 years (all-cohort range, Table 1)",
    age_median      = "21.0 years (all-cohort, Table 1)",
    weight_range    = "40.0-83.0 kg (all-cohort range, Table 1)",
    weight_median   = "56.0 kg (all-cohort, Table 1)",
    sex_female_pct  = 100,
    disease_state   = paste(
      "Uncomplicated Plasmodium falciparum malaria; second- or",
      "third-trimester pregnancy compared with non-pregnant post-partum",
      "controls (matched on history of fever, axillary temperature,",
      "smoking, and parasitemia stratum)."
    ),
    dose_range      = paste(
      "Coartem (Novartis): 20 mg artemether + 120 mg lumefantrine per",
      "tablet; four tablets per dose (480 mg lumefantrine) given twice",
      "daily for 3 days (oral; dose times 0, 8, 24, 36, 48, 60 hours)",
      "with 200 mL milk tea to optimise oral bioavailability."
    ),
    regions         = "Uganda (Mbarara National Referral Hospital antenatal clinic)",
    trial_registration = "ClinicalTrials.gov NCT00495508",
    notes           = paste(
      "Demographics from Kloprogge 2013 Table 1. Pregnant cohort spans",
      "venous (n = 26) and capillary (n = 89; one subject excluded for",
      "an unexplained baseline lumefantrine of 7,717 ng/mL) sampling",
      "arms; the model fits both matrices simultaneously via a",
      "matrix-conversion factor (capillary = 0.881 * venous on the",
      "population mean, 95% CI 0.733-1.05, not significant on backward",
      "elimination -- not encoded in this model file, see Errata).",
      "Gestational age range 13.0-39.0 weeks in the pregnant arms;",
      "admission body temperature range 36.0-39.8 degC across the",
      "all-cohort window."
    )
  )

  ini({
    # Structural population mean parameters come from Table 2 'Fixed
    # effect' column of Kloprogge 2013. The paper reports clearance and
    # volume estimates on the linear scale (CL/F = 5.09 L/h, etc.);
    # log() is applied here for the nlmixr2 internal scale.
    lcl  <- log(5.09)  ; label("Apparent elimination clearance CL/F (L/h)")              # Kloprogge 2013 Table 2: CL/F = 5.09 (RSE 7.90%; 95% CI 4.35-5.87)
    lvc  <- log(123)   ; label("Apparent central volume of distribution Vc/F (L)")       # Kloprogge 2013 Table 2: Vc/F = 123 (RSE 8.40%; 95% CI 104-145)
    lq   <- log(1.68)  ; label("Apparent intercompartmental clearance Q/F (L/h)")        # Kloprogge 2013 Table 2: Q/F  = 1.68 (RSE 10.2%; 95% CI 1.35-2.00)
    lvp  <- log(110)   ; label("Apparent peripheral volume of distribution Vp/F (L)")    # Kloprogge 2013 Table 2: Vp/F = 110 (RSE 9.07%; 95% CI 91.7-131)
    lmtt <- log(4.09)  ; label("Mean absorption transit time MTT (h)")                   # Kloprogge 2013 Table 2: MTT  = 4.09 (RSE 5.22%; 95% CI 3.70-4.55)

    # Relative bioavailability is anchored at 1 (structural, fixed by the
    # source paper, with all variability around F captured by IIV /
    # IOV). The paper also fixes the number of transit compartments at
    # 5; that count is a structural feature encoded in model() below,
    # not an ini() parameter.
    lfdepot <- fixed(log(1)) ; label("Relative bioavailability F1 (unitless)")           # Kloprogge 2013 Table 2: F = 'One fixed' (no estimation, no RSE)

    # Covariate effects. Both are linear-deviation forms; see notes on
    # the covariateData entries above and the Discussion / Methods of
    # Kloprogge 2013 (the paper writes the categorical pregnancy effect
    # as a 'relative difference (%) between groups' and the continuous
    # body-temperature effect as the linear form of Equation 1 centered
    # on the population median).
    e_preg_q       <- -0.365 ; label("Pregnancy effect on Q/F: Q_pregnant / Q_nonpregnant - 1 = -0.365") # Kloprogge 2013 Table 2: 'Pregnancy on Q' = -0.365 (RSE 14.3%; 95% CI -0.455 to -0.259)
    e_bodytemp_mtt <-  0.165 ; label("Body-temperature effect on MTT: per-degC fractional increase (1/degC), centered at 36.9 degC") # Kloprogge 2013 Table 2: 'Temperature on MTT' = 0.165 (RSE 44.1%; 95% CI 0.0328-0.329)

    # IIV. Kloprogge 2013 Table 2 reports CV% (footnote: 100 *
    # sqrt(exp(omega^2) - 1)). The internal variances are recovered by
    # omega^2 = log((CV/100)^2 + 1):
    #   CL  CV 17.5% -> log(0.175^2 + 1) = 0.030164
    #   Vp  CV 21.5% -> log(0.215^2 + 1) = 0.045183
    #   MTT CV 36.6% -> log(0.366^2 + 1) = 0.125731
    #   F   CV 44.8% -> log(0.448^2 + 1) = 0.182910
    etalcl     ~ 0.030164  # Kloprogge 2013 Table 2: BSV on CL/F  = 17.5% CV (RSE 20.6%)
    etalvp     ~ 0.045183  # Kloprogge 2013 Table 2: BSV on Vp/F  = 21.5% CV (RSE 65.2%)
    etalmtt    ~ 0.125731  # Kloprogge 2013 Table 2: BSV on MTT   = 36.6% CV (RSE 32.0%)
    etalfdepot ~ 0.182910  # Kloprogge 2013 Table 2: BSV on F     = 44.8% CV (RSE 31.2%)

    # Residual error. Kloprogge 2013 modelled the natural logarithm of
    # the lumefantrine plasma concentration with an additive error on
    # log scale, which is equivalent to proportional error in nlmixr2's
    # linear-concentration space (see references/parameter-names.md
    # 'Residual error'). Table 2 reports the venous variance on the
    # log scale as sigma_venous = 0.0595; the corresponding SD is
    # sqrt(0.0595) = 0.244, encoded here as propSd. The smaller
    # capillary residual variance sigma_capillary = 0.0207 is not
    # reproduced in this single-output model file (see Errata).
    propSd <- sqrt(0.0595) ; label("Proportional residual SD for venous lumefantrine plasma concentration (SD on log scale)") # Kloprogge 2013 Table 2: sigma_venous = 0.0595 (RSE 15.3%; 95% CI 0.0382-0.0914)
  })

  model({
    # Individual PK parameters. Kloprogge 2013 did not retain body
    # weight as a covariate in the final model after backward
    # elimination, so no allometric scaling is applied (the structural
    # parameters are reported as apparent values for the all-cohort
    # population mean, not normalised to a reference weight). Q/F
    # carries the pregnancy effect; MTT carries the body-temperature
    # effect; CL, Vc, Vp are not modified by retained covariates.
    cl  <- exp(lcl  + etalcl)
    vc  <- exp(lvc)
    q   <- exp(lq)              * (1 + e_preg_q * PREG)
    vp  <- exp(lvp  + etalvp)
    mtt <- exp(lmtt + etalmtt)  * (1 + e_bodytemp_mtt * (BODYTEMP - 36.9))

    # Transit-absorption chain rate constant. NN = 5 transit
    # compartments fixed by the paper; with the absorption chain
    # depot -> transit1 -> transit2 -> transit3 -> transit4 ->
    # transit5 -> central (NN + 1 = 6 transitions at rate ktr), the
    # mean transit time is MTT = (NN + 1) / ktr, so ktr = 6 / MTT.
    # This matches the Savic & Karlsson (2007) convention used by
    # Birgersson 2019 artesunate (NN = 3, ktr = 4 / MTT) in the same
    # lab's NONMEM idiom.
    ktr <- 6 / mtt

    # Two-compartment disposition micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system: 5-compartment transit absorption chain feeding into a
    # two-compartment disposition model. The same ktr propagates the
    # dose through the entire depot + transit chain.
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(transit3)    <-  ktr * transit2 - ktr * transit3
    d/dt(transit4)    <-  ktr * transit3 - ktr * transit4
    d/dt(transit5)    <-  ktr * transit4 - ktr * transit5
    d/dt(central)     <-  ktr * transit5 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                    k12 * central - k21 * peripheral1

    # Relative bioavailability applied to the depot (dosing) compartment.
    f(depot) <- exp(lfdepot + etalfdepot)

    # Lumefantrine plasma concentration in ug/mL. Dose units mg, Vc
    # units L -> central / vc has units mg/L = ug/mL. To compare with
    # the paper's day-7 thresholds reported in ng/mL (175 and 280
    # ng/mL), multiply Cc by 1000 in the vignette.
    Cc <- central / vc

    # Proportional residual error on the linear-concentration scale
    # (NONMEM additive-on-log-scale maps to proportional in nlmixr2;
    # see ini() comment on propSd).
    Cc ~ prop(propSd)
  })
}
