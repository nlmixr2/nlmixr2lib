# Population PK model for oral piperaquine in African pregnant women with
# uncomplicated Plasmodium falciparum malaria (PREGACT sub-study,
# NCT00852423; Ding 2024, CPT Pharmacometrics Syst Pharmacol
# 13:1893-1903; doi:10.1002/psp4.13211).

Ding_2024_piperaquine <- function() {
  description <- paste(
    "Population PK model for oral piperaquine in African women in the",
    "second and third trimester of pregnancy with uncomplicated",
    "Plasmodium falciparum malaria (Ding 2024, PREGACT phase 3 sub-study,",
    "n = 755). Two-transit-compartment absorption with kA = kTR feeding a",
    "three-compartment disposition model. Allometric body-weight scaling",
    "on all apparent clearances (fixed exponent 0.75) and apparent volumes",
    "(fixed exponent 1.0) at a reference weight of 70 kg. Relative",
    "bioavailability is anchored at 1 with inter-individual variability,",
    "decreases by 11.9% per log10 unit of baseline parasitaemia above the",
    "cohort median, and increases by a literature-fixed 23.7% with each",
    "successive dose occasion. Neither gestational age nor trimester was a",
    "significant covariate. Predictions are plasma piperaquine base",
    "concentrations in ng/mL.",
    sep = " "
  )
  reference <- paste(
    "Ding J, Hoglund RM, Tagbor H, Tinto H, Valea I, Mwapasa V,",
    "Kalilani-Phiri L, Van Geertruyden JP, Nambozi M, Mulenga M,",
    "Hachizovu S, Ravinetto R, D'Alessandro U, Tarning J (2024).",
    "Population pharmacokinetics of amodiaquine and piperaquine in African",
    "pregnant women with uncomplicated Plasmodium falciparum infections.",
    "CPT: Pharmacometrics & Systems Pharmacology 13(11):1893-1903.",
    "doi:10.1002/psp4.13211.",
    sep = " "
  )
  vignette <- "Ding_2024_antimalarials_pregnancy"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Ding 2024 Methods ("Blood samples",
  # "Concentration quantification": venous plasma piperaquine by LC/MS-MS)
  # and Table 3.
  compartmentData <- list(
    depot       = list(analyte = "piperaquine", units = "mg", specimen = "administration site", verified = TRUE),
    transit1    = list(analyte = "piperaquine", units = "mg", specimen = "administration site", verified = TRUE),
    transit2    = list(analyte = "piperaquine", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "piperaquine", units = "mg", specimen = "plasma",              verified = TRUE),
    peripheral1 = list(analyte = "piperaquine", units = "mg", specimen = "plasma",              verified = TRUE),
    peripheral2 = list(analyte = "piperaquine", units = "mg", specimen = "plasma",              verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at baseline. Ding 2024 Methods ('Covariates model',",
        "Equations 2 and 3): 'body weight was added on all clearance and",
        "volume parameters using an allometric function with a fixed",
        "exponent of 0.75 and 1, respectively'. Equations 2 and 3 write the",
        "normalising constant generically as BWmedian, but the Table 3",
        "footnote states that 'Population estimates are given for a",
        "typical pregnant women weighted 70 kg', so 70 kg is the reference",
        "encoded here. That reading is confirmed by the paper's own",
        "secondary parameters: at the cohort median weight of 54 kg",
        "(Table 1) a 70 kg reference reproduces the published AUC",
        "(33.1 vs 31.6 h*ug/mL) whereas a 54 kg reference gives 27.2, and",
        "by the Discussion, which compares the 69.9 L/h estimate directly",
        "against a literature CL 'for a typical 70-kg adult patient'. The",
        "in-silico simulations (Results, 'Simulations') likewise use 70 kg",
        "hypothetical patients.",
        sep = " "
      ),
      source_name        = "BW"
    ),
    PARA = list(
      description        = "Baseline (enrolment) asexual Plasmodium falciparum parasitaemia",
      units              = "parasites/uL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Admission-only / time-fixed. Enters relative oral bioavailability",
        "with the log10 transform applied inside model(), in the linear",
        "deviation form given in the Table 3 footnote: 'Baseline parasite",
        "counts (log scale) were implemented on relative oral",
        "bioavailability [1 + theta x (log(parasitemia) - 2.83)]'. The",
        "logarithm is base 10 and the centring constant is the cohort",
        "median: log10(680) = 2.83 for the median of 680 parasites/uL",
        "(Table 1), and the Results simulations use 'a baseline",
        "parasitemia of 676 parasites/uL', which is 10^2.83 = 676.1 to",
        "four significant figures. The same max(PARA, 1) gating used by",
        "the sibling malaria models Kloprogge_2014_quinine.R,",
        "Kloprogge_2018_lumefantrine.R and",
        "Tarning_2012_dihydroartemisinin.R is applied so that values below",
        "1 parasite/uL (the cohort minimum in Table 1 is 0) collapse to",
        "log10 = 0 rather than diverging.",
        sep = " "
      ),
      source_name        = "PARA"
    ),
    OCC = list(
      description        = "Dose-occasion counter (1 = first daily dose, 2 = second, 3 = third)",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Ding 2024 did not estimate a dose-occasion effect because 'most",
        "PK samples were collected after the last dose' (Discussion);",
        "instead Methods ('Population PK analysis') states that 'for",
        "piperaquine, the effect of dose occasion on the relative",
        "bioavailability was fixed to 23.7% increase with each dose, as",
        "reported in a large individual-level data population",
        "pharmacokinetic meta-analysis' -- Hoglund 2017 (PLoS Med",
        "14:e1002212), reference 28 of this paper and the source of the",
        "sibling model Hoglund_2017_piperaquine.R. The Table 3 footnote",
        "gives the form: 'Dose occasion was implemented on relative oral",
        "bioavailability [1 + theta x (OCC - 1)]', so OCC = 1 -> F = 1.000,",
        "OCC = 2 -> F = 1.237, OCC = 3 -> F = 1.474. OCC is carried on each",
        "dose event row; the value is used only by f(depot) at dose",
        "administration.",
        sep = " "
      ),
      source_name        = "OCC"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 755L,
    n_studies       = 1L,
    n_observations  = "976 piperaquine plasma concentrations included (982 measured, 6 (0.6%) below the LLOQ omitted; Abstract, Results)",
    age_range       = "20 (15-43) years, median (min-max) (Table 1)",
    weight_range    = "54 (35-115) kg, median (min-max) (Table 1)",
    height_range    = "155 (138-178) cm, median (min-max) (Table 1)",
    ega_range       = "24 (16-36) weeks gestational age, median (min-max) (Table 1)",
    trimester       = "519/755 (69.1%) in the second trimester, remainder in the third (Table 1)",
    sex_female_pct  = 100,
    pregnant_pct    = 100,
    disease_state   = paste(
      "Acute uncomplicated Plasmodium falciparum mono-infection in pregnant",
      "women in the second or third trimester. Median parasitaemia at",
      "enrolment 680 (5-355,400) parasites/uL; median gametocytaemia 0",
      "(0-1200) parasites/uL (Table 1)."
    ),
    dose_range      = paste(
      "Dihydroartemisinin-piperaquine (Sigma-Tau) 3 tablets once daily for",
      "3 consecutive days under direct observation. One tablet contains",
      "40 mg dihydroartemisinin and 320 mg piperaquine tetraphosphate =",
      "171 mg piperaquine base, so the daily dose is 513 mg piperaquine",
      "base (Methods, 'Drug regimen'). Median daily dose 17.8 (8.3-27.4)",
      "mg/kg of piperaquine tetraphosphate (Table 1)."
    ),
    regions         = "Burkina Faso (2 sites), Ghana (3 sites), Malawi (1 site), Zambia (1 site)",
    notes           = paste(
      "PREGACT trial (NCT00852423), a non-inferiority, multi-centre,",
      "randomised, open-label phase 3 trial of 4 artemisinin-based",
      "combination therapies conducted June 2010 to August 2013. Of 763",
      "women randomised to dihydroartemisinin-piperaquine, 8 (1.0%) who",
      "vomited after dosing were excluded from the PK analysis, leaving",
      "755. A venous sample was taken from every woman on day 7 with",
      "additional samples at other clinical visits when possible, so the",
      "data are sparse and essentially uninformative in the absorption",
      "phase. Gestational age and trimester were both non-significant on",
      "every PK parameter (Results), as was age."
    )
  )

  ini({
    # ---- Absorption -------------------------------------------------
    # Ding 2024 Methods ('Population PK analysis'): "Due to uninformative
    # sample collection in absorption phase, absorption relevant
    # parameters such as the absorption rate (ka), mean transit time
    # (MTT), and number of transit compartments were fixed to literature
    # values" -- for piperaquine the literature source is Hoglund 2017
    # (PLoS Med 14:e1002212), reference 28 of the paper, which reports the
    # identical MTT of 2.11 h with 2 transit compartments (see the sibling
    # model Hoglund_2017_piperaquine.R).
    lmtt <- fixed(log(2.11))
    label("Mean transit time of the 2-transit-compartment absorption chain, literature value (MTT, h)")
    # Ding 2024 Table 3: MTT = 2.11 h fixed; number of transit
    # compartments = 2 fixed. No separate ka is reported for piperaquine
    # (contrast Table 2 for amodiaquine, which reports both), so kA = kTR
    # as in Hoglund 2017: the chain depot -> transit1 -> transit2 ->
    # central has (NN + 1) = 3 equal-rate transitions and
    # ktr = 3/MTT (Savic & Karlsson 2007 convention).

    # ---- Disposition (three compartments) ----------------------------
    # Ding 2024 Table 3, "NONMEM Population estimates" column. Values are
    # apparent (relative to F = 1) and reported on the linear scale for a
    # typical pregnant woman weighing 70 kg; log() is applied here for the
    # nlmixr2 internal log scale.
    lcl <- log(69.9)
    label("Apparent piperaquine elimination clearance CL/F at WT = 70 kg (L/h)")
    # Ding 2024 Table 3: CL/F = 69.9 L/h (%RSE 5.1; bootstrap median 69.3,
    # 95% CI 62.5-76.1)

    lvc <- log(4240)
    label("Apparent piperaquine central volume of distribution Vc/F at WT = 70 kg (L)")
    # Ding 2024 Table 3: Vc/F = 4240 L (%RSE 36.0; bootstrap median 4010,
    # 95% CI 1600-7600)

    lq <- log(265)
    label("Apparent inter-compartmental clearance to the first peripheral compartment Q1/F at WT = 70 kg (L/h)")
    # Ding 2024 Table 3: Q1/F = 265 L/h (%RSE 39.9; bootstrap median 243,
    # 95% CI 63-472)

    lvp <- log(3880)
    label("Apparent first peripheral volume of distribution Vp1/F at WT = 70 kg (L)")
    # Ding 2024 Table 3: Vp1/F = 3880 L (%RSE 28.7; bootstrap median 3860,
    # 95% CI 1240-5730)

    lq2 <- log(103)
    label("Apparent inter-compartmental clearance to the second peripheral compartment Q2/F at WT = 70 kg (L/h)")
    # Ding 2024 Table 3: Q2/F = 103 L/h (%RSE 12.8; bootstrap median 100,
    # 95% CI 76-127)

    lvp2 <- log(22900)
    label("Apparent second peripheral volume of distribution Vp2/F at WT = 70 kg (L)")
    # Ding 2024 Table 3: Vp2/F = 22,900 L (%RSE 7.9; bootstrap median
    # 22,800, 95% CI 19,900-26,700)

    # ---- Relative bioavailability ------------------------------------
    # Ding 2024 Methods ('Population PK analysis'): "Relative
    # bioavailability (F) was fixed to unity (100%) for the population to
    # allow the estimation of IIV of the absorption."
    lfdepot <- fixed(log(1))
    label("Relative oral bioavailability of piperaquine at OCC = 1 and median baseline parasitaemia (unitless)")
    # Ding 2024 Table 3: F = 100% fixed

    # ---- Allometric exponents ----------------------------------------
    # Ding 2024 Methods ('Covariates model', Equations 2 and 3): "body
    # weight was added on all clearance and volume parameters using an
    # allometric function with a fixed exponent of 0.75 and 1,
    # respectively". Retained even though it slightly worsened the fit
    # (Results: delta-OFV = 41.4) because pregnant and non-pregnant women
    # are expected to differ systematically in weight.
    e_wt_cl <- fixed(0.75)
    label("Allometric WT exponent on all apparent clearance parameters (CL/F, Q1/F, Q2/F)")
    # Ding 2024 Methods, Equation 2: clearance power = 0.75

    e_wt_vc <- fixed(1.00)
    label("Allometric WT exponent on all apparent volume parameters (Vc/F, Vp1/F, Vp2/F)")
    # Ding 2024 Methods, Equation 3: volume power = 1

    # ---- Covariate effects on relative bioavailability ----------------
    e_para_f <- -0.119
    label("Linear increment in relative oral bioavailability per log10 unit of baseline parasitaemia above 2.83 (per log10 parasites/uL)")
    # Ding 2024 Table 3: Baseline parasites count on F = -11.9% (%RSE
    # 19.3; bootstrap median -12.1, 95% CI -16.2 to -7.3). Table 3
    # footnote gives the form: "Baseline parasite counts (log scale) were
    # implemented on relative oral bioavailability
    # [1 + theta x (log(parasitemia) - 2.83)]". The Results text rounds
    # this to "11.6% decrease in F per one log unit"; the final-estimates
    # table value of -11.9% is used here.

    e_doseocc_f <- fixed(0.237)
    label("Increment in relative oral bioavailability per additional dose occasion, literature value (fraction)")
    # Ding 2024 Table 3: Dose occasion on F = 23.7% fixed, taken from
    # Hoglund 2017 (reference 28) per Methods ('Population PK analysis').

    # ---- Inter-individual variability --------------------------------
    # Ding 2024 Table 3 footnote: "Coefficients of variation for
    # inter-individual variability (IIV) were calculated as
    # 100 x (e^variance - 1)^(1/2)", so the internal log-scale variance is
    # recovered as omega^2 = log(CV^2 + 1).
    #
    #   F     IIV 33.8% -> omega^2 = log(0.338^2 + 1) = 0.1081761
    #   Vc    IIV 112%  -> omega^2 = log(1.120^2 + 1) = 0.8128839
    #
    # Table 3 reports no IIV on CL/F, Q1/F, Vp1/F, Q2/F or Vp2/F, so no
    # eta slots are created for those parameters.
    etalfdepot ~ 0.1081761
    # Ding 2024 Table 3: IIV on F = 33.8% CV (%RSE 9.2; bootstrap median
    # 33.6, 95% CI 27.2-40.2; eta shrinkage 41.0%)

    etalvc ~ 0.8128839
    # Ding 2024 Table 3: IIV on Vc/F = 112% CV (%RSE 29.3; bootstrap
    # median 125, 95% CI 37-308; eta shrinkage 73.1%)

    # ---- Residual unexplained variability -----------------------------
    # Ding 2024 Methods ('Population PK analysis'): "Residual unexplained
    # variability was modeled with an additive error on the
    # log-transformed concentrations, which is essentially equivalent to
    # an exponential residual error on an arithmetic scale." That NONMEM
    # additive-on-log-scale residual maps to a proportional residual in
    # linear concentration space. The Table 2 footnote states "RUV is the
    # residual error variance", so the tabulated number is a variance and
    # the SD is its square root -- the same convention as the sibling
    # Hoglund_2017_piperaquine.R.
    propSd <- sqrt(0.222)
    label("Proportional residual SD for piperaquine plasma concentration (SD on log scale)")
    # Ding 2024 Table 3: RUV = 0.222 (variance; %RSE 9.5; bootstrap median
    # 0.215, 95% CI 0.177-0.256; epsilon shrinkage 17.7%)
  })

  model({
    # Mean transit time and chain rate constant. With NN = 2 transit
    # compartments and kA = kTR, the absorption chain
    # depot -> transit1 -> transit2 -> central has (NN + 1) = 3 equal
    # transitions; ktr = (NN + 1)/MTT = 3/MTT (Savic & Karlsson 2007
    # convention, matching the source of the fixed MTT,
    # Hoglund_2017_piperaquine.R).
    mtt <- exp(lmtt)
    ktr <- 3 / mtt

    # Individual PK parameters. Allometric weight scaling on all apparent
    # clearances (exponent 0.75) and apparent volumes (exponent 1)
    # centred on the 70 kg reference of the Table 3 footnote. IIV is
    # carried only on F and Vc/F (Ding 2024 Table 3).
    cl  <- exp(lcl)             * (WT / 70)^e_wt_cl
    vc  <- exp(lvc  + etalvc)   * (WT / 70)^e_wt_vc
    q   <- exp(lq)              * (WT / 70)^e_wt_cl
    vp  <- exp(lvp)             * (WT / 70)^e_wt_vc
    q2  <- exp(lq2)             * (WT / 70)^e_wt_cl
    vp2 <- exp(lvp2)            * (WT / 70)^e_wt_vc

    # Three-compartment disposition micro-constants (1/h).
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # ODE system (Ding 2024 Figure S2): 2-transit-compartment absorption
    # feeding a 3-compartment disposition model. Compartment amounts are
    # in mg of piperaquine base and volumes are in L, so amount/volume is
    # mg/L and is scaled to ng/mL below.
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(central)     <-  ktr * transit2 - kel * central -
                          k12 * central + k21 * peripheral1 -
                          k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central  - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central  - k31 * peripheral2

    # Relative oral bioavailability. The population anchor
    # lfdepot = log(1) is fixed. Two multiplicative covariate terms apply,
    # both in the linear deviation form of the Ding 2024 Table 3 footnote:
    # baseline parasitaemia centred at log10 = 2.83 (the cohort median of
    # 680 parasites/uL) and the literature-fixed dose-occasion increment.
    # max(PARA, 1) gates parasitaemia values below 1 parasite/uL (the
    # cohort minimum in Table 1 is 0) to log10 = 0, the same convention
    # used by the sibling malaria models.
    f_para   <- 1 + e_para_f    * (log10(max(PARA, 1)) - 2.83)
    f_occ    <- 1 + e_doseocc_f * (OCC - 1)
    f(depot) <- f_para * f_occ * exp(lfdepot + etalfdepot)

    # Plasma piperaquine concentration in ng/mL: amount (mg) / volume (L)
    # is mg/L, multiplied by 1000 to give ng/mL, the units used
    # throughout Ding 2024 Table 3 (day-7 concentration) and Figures 2
    # and 4.
    Cc <- 1000 * central / vc

    # Proportional residual error on the linear-concentration scale, the
    # linear-space equivalent of the paper's additive-on-log-scale error.
    Cc ~ prop(propSd)
  })
}
