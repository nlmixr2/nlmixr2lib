# Population pharmacokinetic model of oral piperaquine in pregnant and
# non-pregnant Sudanese women with uncomplicated Plasmodium falciparum
# malaria after the standard fixed-dose dihydroartemisinin-piperaquine
# regimen (Hoglund 2012, Malaria Journal 11:398;
# doi:10.1186/1475-2875-11-398).

Hoglund_2012_piperaquine <- function() {
  description <- paste(
    "Population PK model for oral piperaquine in pregnant and non-pregnant",
    "Sudanese women with uncomplicated Plasmodium falciparum malaria",
    "(Hoglund 2012). Three-transit-compartment absorption (ka = ktr) into a",
    "three-compartment disposition model. Body weight is the only retained",
    "covariate, applied as an allometric function on all clearances (fixed",
    "exponent 0.75) and volumes (fixed exponent 1.0). Relative bioavailability",
    "F is fixed at 1. The final model retains BSV on CL and F, treats MTT",
    "between-occasion variability as forward-simulation IIV, and uses an",
    "additive residual on the log-transformed observation (proportional in",
    "linear concentration space).",
    sep = " "
  )
  reference <- paste(
    "Hoglund RM, Adam I, Hanpithakpong W, Ashton M, Lindegardh N, Day NPJ,",
    "White NJ, Nosten F, Tarning J (2012).",
    "A population pharmacokinetic model of piperaquine in pregnant and",
    "non-pregnant women with uncomplicated Plasmodium falciparum malaria",
    "in Sudan. Malaria Journal 11:398.",
    "doi:10.1186/1475-2875-11-398.",
    sep = " "
  )
  vignette <- "Hoglund_2012_piperaquine"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight at admission",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling with fixed exponents 0.75 on apparent clearances",
        "(CL/F, Q1/F, Q2/F) and 1.0 on apparent volumes (Vc/F, Vp1/F, Vp2/F).",
        "Hoglund 2012 Methods (page 4): 'Body weight was tried in the model as",
        "a simultaneous incorporation of an allometric function on all",
        "clearance (power of 0.75) and volume parameters (power of 1),",
        "considering the strong biological prior of this covariate",
        "relationship'. The paper does not state the allometric reference",
        "weight explicitly; the cohort medians were 53 kg (non-pregnant) and",
        "59 kg (pregnant) per Table 1, so the pooled-cohort approximate",
        "median 56 kg is used here. Forward-simulating a 53 kg non-pregnant",
        "patient under this reference reproduces the published Table 3",
        "AUC0-90 of ~38000 ng*h/mL to within 3%. Time-fixed at admission",
        "in the source data set.",
        sep = " "
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 24L,
    n_studies       = 1L,
    n_pregnant      = 12L,
    n_nonpregnant   = 12L,
    age_range       = "16.0-43.0 years (all-cohort range, Table 1)",
    age_median      = "21.0 years non-pregnant; 26.0 years pregnant (Table 1)",
    weight_range    = "44.0-81.0 kg (all-cohort range, Table 1)",
    weight_median   = "53.0 kg non-pregnant; 59.0 kg pregnant (Table 1)",
    height_range    = "150-174 cm (all-cohort, Table 1)",
    sex_female_pct  = 100,
    ega_range       = "15.3-40.1 weeks among the pregnant cohort (Table 1)",
    disease_state   = paste(
      "Uncomplicated Plasmodium falciparum malaria; symptomatic women,",
      "12 in their second or third trimester of pregnancy and 12",
      "non-pregnant controls."
    ),
    dose_range      = paste(
      "Dihydroartemisinin-piperaquine fixed-dose combination (Duo Cotecxin,",
      "40 mg dihydroartemisinin + 320 mg piperaquine tetra-phosphate per",
      "tablet) once daily for three days, directly observed, taken with water",
      "under fasting conditions. Tablet count titrated to a daily dose of",
      "20 mg piperaquine tetra-phosphate/kg, equivalent to ~10.5 mg",
      "piperaquine base/kg (Table 1)."
    ),
    regions         = "Sudan (New Halfa Teaching Hospital, New Halfa)",
    notes           = paste(
      "Demographics from Hoglund 2012 Table 1. Sample size n=24 (12+12)",
      "post-randomisation; 14 non-pregnant women were initially recruited",
      "but 2 withdrew consent. PK analysis used 564 post-dose plasma",
      "samples (7 BLQ samples omitted). Sampling pre-dose and at 1.5, 4,",
      "8, 24, 25.5, 28, 32, 48, 49, 50, 52, 56, 60, 72 h after the first",
      "dose and on days 5, 7, 14, 21, 28, 35, 42, 49, 56, 63, 90."
    )
  )

  ini({
    # Structural parameters from Hoglund 2012 Table 2 ("Population estimates"
    # column). Apparent values (relative to F) reported on the linear scale;
    # log() applied here for the nlmixr2 internal log scale. Typical values
    # correspond to the allometric reference weight 56 kg (the pooled-cohort
    # approximate median; the paper does not state the reference, see WT
    # covariate notes).
    lcl  <- log(44.6)
    label("Apparent piperaquine elimination clearance CL/F at WT = 56 kg (L/h)")
    # Hoglund 2012 Table 2: CL/F = 44.6 L/h (RSE 9.90%; 95% CI 37.3-53.8)

    lvc  <- log(1820)
    label("Apparent central volume of distribution Vc/F at WT = 56 kg (L)")
    # Hoglund 2012 Table 2: Vc/F = 1820 L (RSE 11.5%; 95% CI 1450-2240)

    lq   <- log(47.7)
    label("Apparent inter-compartmental clearance to first peripheral compartment Q1/F at WT = 56 kg (L/h)")
    # Hoglund 2012 Table 2: Q1/F = 47.7 L/h (RSE 19.0%; 95% CI 32.4-69.2)

    lvp  <- log(15900)
    label("Apparent first peripheral volume of distribution Vp1/F at WT = 56 kg (L)")
    # Hoglund 2012 Table 2: Vp1/F = 15900 L (RSE 12.3%; 95% CI 12600-20400)

    lq2  <- log(352)
    label("Apparent inter-compartmental clearance to second peripheral compartment Q2/F at WT = 56 kg (L/h)")
    # Hoglund 2012 Table 2: Q2/F = 352 L/h (RSE 11.1%; 95% CI 283-431)

    lvp2 <- log(7520)
    label("Apparent second peripheral volume of distribution Vp2/F at WT = 56 kg (L)")
    # Hoglund 2012 Table 2: Vp2/F = 7520 L (RSE 17.1%; 95% CI 5520-10500)

    lmtt <- log(1.70)
    label("Mean transit time of the 3-compartment transit-absorption chain MTT (h)")
    # Hoglund 2012 Table 2: MTT = 1.70 h (RSE 8.05%; 95% CI 1.45-2.00)

    # Relative bioavailability fixed at 1 (Hoglund 2012 Methods page 4:
    # "Bioavailability was added to the model and the population value was
    # fixed to 100%"). All F variability is captured by IIV around F (BSV
    # and BOV in the source paper; see etalfdepot below).
    lfdepot <- fixed(log(1))
    label("Relative bioavailability F (unitless, fixed at 1)")
    # Hoglund 2012 Table 2: F (%) = 100 fix (no estimation, no RSE)

    # Allometric exponents fixed by the source paper (strong biological prior;
    # exponents are written into the model rather than estimated).
    e_wt_cl <- fixed(0.75)
    label("Allometric WT exponent on all apparent clearance parameters (CL/F, Q1/F, Q2/F)")
    # Hoglund 2012 Methods page 4: clearance power = 0.75 (fixed)

    e_wt_vc <- fixed(1.00)
    label("Allometric WT exponent on all apparent volume parameters (Vc/F, Vp1/F, Vp2/F)")
    # Hoglund 2012 Methods page 4: volume power = 1 (fixed)

    # IIV. Coefficient of variation reported in Hoglund 2012 Table 2 (footnote:
    # "%CV ... calculated as the [exp(estimated variance)-1]^(1/2)", i.e.
    # CV = sqrt(exp(omega^2) - 1), so the internal log-scale variance is
    # recovered as omega^2 = log(CV^2 + 1)).
    #
    #   CL  BSV 22.5%  -> omega^2 = log(0.225^2 + 1) = 0.049389
    #   MTT BOV 60.7%  -> omega^2 = log(0.607^2 + 1) = 0.313762  (BOV, treated as forward-simulation IIV)
    #   F   BSV 34.7%  -> omega^2 = log(0.347^2 + 1) = 0.113707  (BSV only; BOV 64.8% omitted, see vignette Errata)
    etalcl     ~ 0.049389
    # Hoglund 2012 Table 2: BSV on CL/F = 22.5% CV (RSE 31.5%; 95% CI 14.6-29.0)

    etalmtt    ~ 0.313762
    # Hoglund 2012 Table 2: BOV on MTT = 60.7% CV (RSE 22.8%; 95% CI 44.5-76.7).
    # The published variability is between-occasion (per-dose); nlmixr2lib does
    # not natively model BOV, so the value is encoded here as a forward-
    # simulation IIV term on MTT (same approximation used in
    # Kloprogge_2013_lumefantrine.R). See vignette Errata.

    etalfdepot ~ 0.113707
    # Hoglund 2012 Table 2: BSV on F = 34.7% CV (RSE 59.2%; 95% CI 9.52-54.9).
    # The published F additionally carries a BOV of 64.8% CV (omega^2 =
    # log(0.648^2 + 1) = 0.350506), not encoded here. See vignette Errata.

    # Residual error. Hoglund 2012 Methods page 4: "An additive residual error
    # model was assumed since data were transformed into their natural
    # logarithms (i.e. essentially equivalent to an exponential error model
    # on an arithmetic scale)." This NONMEM additive-on-log-scale residual
    # maps to nlmixr2 proportional residual in linear concentration space
    # (see references/parameter-names.md Residual error). Table 2 reports the
    # log-scale variance RUV = 0.0973; the corresponding SD is sqrt(0.0973) =
    # 0.3119, which equals the proportional CV to first order.
    propSd <- sqrt(0.0973)
    label("Proportional residual SD for piperaquine plasma concentration (SD on log scale)")
    # Hoglund 2012 Table 2: RUV = 0.0973 (variance, RSE 5.90%; 95% CI 0.0753-0.120)
  })

  model({
    # Individual PK parameters. Allometric scaling on CL with exponent 0.75
    # and on V with exponent 1, centred on WT = 56 kg (the pooled-cohort
    # approximate median; the paper does not state the reference, see WT
    # covariate notes). The Q1 and Q2 parameters take the clearance exponent
    # (0.75); the peripheral volumes take the volume exponent (1.0).
    cl  <- exp(lcl  + etalcl) * (WT / 56)^e_wt_cl
    vc  <- exp(lvc)           * (WT / 56)^e_wt_vc
    q   <- exp(lq)            * (WT / 56)^e_wt_cl
    vp  <- exp(lvp)           * (WT / 56)^e_wt_vc
    q2  <- exp(lq2)           * (WT / 56)^e_wt_cl
    vp2 <- exp(lvp2)          * (WT / 56)^e_wt_vc

    # Mean transit time and chain rate constant. The Hoglund 2012 final
    # absorption model uses 3 transit compartments with ka = ktr (Hoglund
    # 2012 Methods page 4 and Results page 4: "the absorption rate constant
    # (ka) and the transit-compartment rate constant (ktr) were set equal").
    # With NN = 3 transit compartments and the absorption chain
    # depot -> transit1 -> transit2 -> transit3 -> central (NN + 1 = 4
    # transitions at rate ktr), the mean transit time is MTT = (NN + 1)/ktr,
    # so ktr = 4 / MTT (Savic & Karlsson 2007 convention; same idiom as the
    # Tarning group's Birgersson_2019_artesunate and Kloprogge_2013_lumefantrine
    # implementations).
    mtt <- exp(lmtt + etalmtt)
    ktr <- 4 / mtt

    # Three-compartment disposition micro-constants.
    kel <- cl  / vc
    k12 <- q   / vc
    k21 <- q   / vp
    k13 <- q2  / vc
    k31 <- q2  / vp2

    # ODE system: 3-compartment transit absorption chain feeding into a
    # three-compartment disposition model. The same ktr propagates the dose
    # through the entire depot + transit chain (Figure 1 of Hoglund 2012).
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(transit3)    <-  ktr * transit2 - ktr * transit3
    d/dt(central)     <-  ktr * transit3 - kel * central -
                          k12 * central + k21 * peripheral1 -
                          k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central  - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central  - k31 * peripheral2

    # Relative bioavailability applied to the depot (dosing) compartment.
    # Anchor is fixed at 1; variability is captured by etalfdepot (BSV on F).
    f(depot) <- exp(lfdepot + etalfdepot)

    # Piperaquine plasma concentration. Dose units mg (piperaquine base),
    # Vc units L -> central / vc has units mg/L. Multiplying by 1000 converts
    # to ng/mL for direct comparison against Hoglund 2012 Table 3 (Cmax,
    # day 7, day 28 concentrations all in ng/mL).
    Cc <- 1000 * central / vc

    # Proportional residual error on the linear-concentration scale (NONMEM
    # additive-on-log-scale -> nlmixr2 proportional; see ini() comment on
    # propSd).
    Cc ~ prop(propSd)
  })
}
