# Joint parent-metabolite population PK model for oral amodiaquine and its
# active metabolite desethylamodiaquine in African pregnant women with
# uncomplicated Plasmodium falciparum malaria (PREGACT sub-study,
# NCT00852423; Ding 2024, CPT Pharmacometrics Syst Pharmacol
# 13:1893-1903; doi:10.1002/psp4.13211).

Ding_2024_amodiaquine <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for oral amodiaquine and",
    "its active CYP2C8-derived metabolite desethylamodiaquine in African",
    "women in the second and third trimester of pregnancy with",
    "uncomplicated Plasmodium falciparum malaria (Ding 2024, PREGACT",
    "phase 3 sub-study, n = 771). Two-transit-compartment absorption",
    "feeding a one-compartment amodiaquine disposition model, with",
    "complete bioconversion (molar-corrected) to a two-compartment",
    "desethylamodiaquine disposition model. Allometric body-weight scaling",
    "on all apparent clearances (fixed exponent 0.75) and apparent volumes",
    "(fixed exponent 1.0) at a reference weight of 70 kg. Relative",
    "bioavailability is anchored at 1 with inter-individual variability and",
    "increases linearly by 1.28% per week of gestational age, referenced at",
    "24 weeks. Predictions are plasma amodiaquine and desethylamodiaquine",
    "base concentrations in ng/mL.",
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
  # "Concentration quantification": venous plasma amodiaquine and
  # desethylamodiaquine by LC/MS-MS) and Table 2.
  compartmentData <- list(
    depot            = list(analyte = "amodiaquine",         units = "mg", specimen = "administration site", verified = TRUE),
    transit1         = list(analyte = "amodiaquine",         units = "mg", specimen = "administration site", verified = TRUE),
    transit2         = list(analyte = "amodiaquine",         units = "mg", specimen = "administration site", verified = TRUE),
    central          = list(analyte = "amodiaquine",         units = "mg", specimen = "plasma",              verified = TRUE),
    central_deaq     = list(analyte = "desethylamodiaquine", units = "mg", specimen = "plasma",              verified = TRUE),
    peripheral1_deaq = list(analyte = "desethylamodiaquine", units = "mg", specimen = "plasma",              verified = TRUE)
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
        "normalising constant generically as BWmedian, but the Table 2",
        "footnote states that 'Population estimates are given for a",
        "typical pregnant women weighting 70 kg', so 70 kg is the reference",
        "encoded here. That reading is confirmed by the paper's own",
        "secondary parameters: at the cohort median weight of 55 kg",
        "(Table 1) a 70 kg reference reproduces the published AUC_AQ",
        "(286 vs 283 h*ng/mL) and the published desethylamodiaquine",
        "terminal half-life (14.2 vs 14.1 days), whereas a 55 kg reference",
        "misses both. The in-silico simulations (Results, 'Simulations')",
        "likewise use 70 kg hypothetical patients.",
        sep = " "
      ),
      source_name        = "BW"
    ),
    EGA = list(
      description        = "Maternal estimated gestational age at dosing",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Estimated by fundal height (Ding 2024 Methods, 'Covariates",
        "model'). Enters relative oral bioavailability as the linear",
        "deviation form given in the Table 2 footnote:",
        "F = 1 + theta * (GA - 24), with theta = 0.0128 per week",
        "(Table 2 'Gestational age on F AQ (%) = 1.28'). The reference of",
        "24 weeks is the cohort median gestational age (Table 1) and is",
        "named explicitly in Results: 'Considering pregnant women with",
        "24 weeks of GA as a reference population'. Cohort range 13-36",
        "weeks; the linear form is only supported inside that range and is",
        "NOT extrapolable to the EGA = 0 non-pregnant anchor described in",
        "the covariate register (it would give F = 0.69 there, which the",
        "study did not observe because no non-pregnant women were",
        "enrolled).",
        sep = " "
      ),
      source_name        = "GA"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 771L,
    n_studies       = 1L,
    n_observations  = "1071 amodiaquine and 1087 desethylamodiaquine plasma concentrations (Abstract, Results)",
    age_range       = "22 (15-43) years, median (min-max) (Table 1)",
    weight_range    = "55 (37-104) kg, median (min-max) (Table 1)",
    height_range    = "158 (132-179) cm, median (min-max) (Table 1)",
    ega_range       = "24 (13-36) weeks gestational age, median (min-max) (Table 1)",
    trimester       = "581/771 (75.7%) in the second trimester, remainder in the third (Table 1)",
    sex_female_pct  = 100,
    pregnant_pct    = 100,
    disease_state   = paste(
      "Acute uncomplicated Plasmodium falciparum mono-infection in pregnant",
      "women in the second or third trimester. Median parasitaemia at",
      "enrolment 560 (0-82,292) parasites/uL; median gametocytaemia 0",
      "(0-253) parasites/uL (Table 1)."
    ),
    dose_range      = paste(
      "Artesunate-amodiaquine (Sanofi-Aventis) 2 tablets once daily for 3",
      "consecutive days under direct observation. One tablet contains",
      "100 mg artesunate and 352.64 mg amodiaquine salt = 270 mg",
      "amodiaquine base, so the daily dose is 540 mg amodiaquine base",
      "(Methods, 'Drug regimen'). Median daily dose 12.8 (6.8-19.1) mg/kg",
      "of amodiaquine salt (Table 1)."
    ),
    regions         = "Burkina Faso (2 sites), Ghana (3 sites), Malawi (1 site), Zambia (1 site)",
    notes           = paste(
      "PREGACT trial (NCT00852423), a non-inferiority, multi-centre,",
      "randomised, open-label phase 3 trial of 4 artemisinin-based",
      "combination therapies conducted June 2010 to August 2013. Of 784",
      "women randomised to artesunate-amodiaquine, 13 (1.7%) who vomited",
      "after dosing were excluded from the PK analysis, leaving 771. A",
      "venous sample was taken from every woman on day 7 with additional",
      "samples at other clinical visits when possible, so the data are",
      "sparse and essentially uninformative in the absorption phase. 848",
      "of 1071 (79.2%) amodiaquine and 29 of 1087 (2.7%)",
      "desethylamodiaquine measurements were below the LLOQ and were",
      "omitted (Beal M1), which a categorical VPC supported (Results)."
    )
  )

  ini({
    # ---- Absorption -------------------------------------------------
    # Ding 2024 Methods ('Population PK analysis'): "Due to uninformative
    # sample collection in absorption phase, absorption relevant parameters
    # such as the absorption rate (ka), mean transit time (MTT), and number
    # of transit compartments were fixed to literature values" -- for
    # amodiaquine the literature source is Tarning 2012 (AAC
    # 56:5764-5773), reference 9 of the paper.
    lka <- fixed(log(0.589))
    label("Absorption rate constant from the last transit compartment into amodiaquine central, literature value (1/h)")
    # Ding 2024 Table 2: Ka = 0.589 1/h fixed

    lmtt <- fixed(log(0.236))
    label("Mean transit time of the 2-transit-compartment absorption chain, literature value (h)")
    # Ding 2024 Table 2: MTT = 0.236 h fixed; number of transit
    # compartments = 2 fixed. Because MTT (0.236 h) is SHORTER than the
    # mean time of the single ka step (1/0.589 = 1.70 h), MTT provably
    # excludes the ka step: a mean transit time cannot be smaller than one
    # of the mean times it sums. MTT therefore spans only the two transit
    # compartments, giving ktr = 2/MTT = 8.475 1/h, and ka is the final
    # transfer into central. Mean absorption time = MTT + 1/ka = 1.93 h.
    # This is why the chain rate here is 2/MTT while the sibling
    # Hoglund_2017_piperaquine.R and Ali_2018_amodiaquine.R use 3/MTT:
    # in those models kA = kTR, so all three transfers share one rate and
    # belong to MTT. See the vignette Assumptions and deviations section.

    # ---- Amodiaquine disposition (one compartment) -------------------
    # Ding 2024 Table 2, "NONMEM Estimates" column. Values are apparent
    # (relative to F = 1) and reported on the linear scale for a typical
    # pregnant woman weighing 70 kg; log() is applied here for the
    # nlmixr2 internal log scale.
    lcl <- log(6780)
    label("Apparent amodiaquine elimination clearance CL/F at WT = 70 kg (L/h)")
    # Ding 2024 Table 2: CL/F AQ = 6780 L/h (%RSE 4.9; bootstrap median
    # 6770, 95% CI 6130-7440)

    lvc <- log(272000)
    label("Apparent amodiaquine central volume of distribution Vc/F at WT = 70 kg (L)")
    # Ding 2024 Table 2: Vc/F AQ = 272,000 L (%RSE 8.9; bootstrap median
    # 270,000, 95% CI 228,000-322,000)

    # ---- Desethylamodiaquine disposition (two compartments) ----------
    lcl_deaq <- log(38.3)
    label("Apparent desethylamodiaquine elimination clearance CL/F at WT = 70 kg (L/h)")
    # Ding 2024 Table 2: CL/F DEAQ = 38.3 L/h (%RSE 9.6; bootstrap median
    # 38.0, 95% CI 30.4-43.6)

    lvc_deaq <- log(861)
    label("Apparent desethylamodiaquine central volume of distribution Vc/F at WT = 70 kg (L)")
    # Ding 2024 Table 2: Vc/F DEAQ = 861 L (%RSE 15.7; bootstrap median
    # 897, 95% CI 656-1240)

    lq_deaq <- log(81.7)
    label("Apparent desethylamodiaquine inter-compartmental clearance Q/F at WT = 70 kg (L/h)")
    # Ding 2024 Table 2: Q/F DEAQ = 81.7 L/h (%RSE 7.0; bootstrap median
    # 82.2, 95% CI 72.1-94.5)

    lvp_deaq <- log(13200)
    label("Apparent desethylamodiaquine peripheral volume of distribution Vp/F at WT = 70 kg (L)")
    # Ding 2024 Table 2: Vp/F DEAQ = 13,200 L (%RSE 13.4; bootstrap median
    # 13,400, 95% CI 11,200-17,100)

    # ---- Relative bioavailability ------------------------------------
    # Ding 2024 Methods ('Population PK analysis'): "Relative
    # bioavailability (F) was fixed to unity (100%) for the population to
    # allow the estimation of IIV of the absorption."
    lfdepot <- fixed(log(1))
    label("Relative oral bioavailability of amodiaquine at EGA = 24 weeks (unitless)")
    # Ding 2024 Table 2: F AQ = 100% fixed

    # ---- Allometric exponents ----------------------------------------
    # Ding 2024 Methods ('Covariates model', Equations 2 and 3): "body
    # weight was added on all clearance and volume parameters using an
    # allometric function with a fixed exponent of 0.75 and 1,
    # respectively". Retained even though it slightly worsened the fit
    # (Results: delta-OFV = 35.2) because pregnant and non-pregnant women
    # are expected to differ systematically in weight.
    e_wt_cl <- fixed(0.75)
    label("Allometric WT exponent on all apparent clearance parameters (CL/F AQ, CL/F DEAQ, Q/F DEAQ)")
    # Ding 2024 Methods, Equation 2: clearance power = 0.75

    e_wt_vc <- fixed(1.00)
    label("Allometric WT exponent on all apparent volume parameters (Vc/F AQ, Vc/F DEAQ, Vp/F DEAQ)")
    # Ding 2024 Methods, Equation 3: volume power = 1

    # ---- Gestational-age effect on relative bioavailability ----------
    e_ega_f <- 0.0128
    label("Linear increment in relative oral bioavailability per week of gestational age above 24 weeks (per week)")
    # Ding 2024 Table 2: Gestational age on F AQ = 1.28% per week
    # (%RSE 25.0; bootstrap median 1.29, 95% CI 0.65-1.92). Table 2
    # footnote gives the form: "Gestational age (GA) was implemented on
    # relative oral bioavailability [1 + (theta x (GA - 24))]". The
    # Abstract rounds this to 1.25%/week; Results and Table 2 both give
    # 1.28%/week and the table value is used here.

    # ---- Inter-individual variability --------------------------------
    # Ding 2024 Table 2 footnote: "Coefficients of variation for
    # inter-individual variability (IIV) were calculated as
    # 100 x (e^variance - 1)^(1/2)", so the internal log-scale variance is
    # recovered as omega^2 = log(CV^2 + 1).
    #
    #   F AQ      IIV 30.6% -> omega^2 = log(0.306^2 + 1) = 0.0895079
    #   CL DEAQ   IIV 19.7% -> omega^2 = log(0.197^2 + 1) = 0.0380749
    #   Vc DEAQ   IIV 205%  -> omega^2 = log(2.050^2 + 1) = 1.6491393
    #
    # Table 2 reports no IIV on CL/F AQ, Vc/F AQ, Q/F DEAQ or Vp/F DEAQ,
    # so no eta slots are created for those parameters.
    etalfdepot ~ 0.0895079
    # Ding 2024 Table 2: IIV on F AQ = 30.6% CV (%RSE 8.4; bootstrap
    # median 30.4, 95% CI 24.9-35.4; eta shrinkage 33.7%)

    etalcl_deaq ~ 0.0380749
    # Ding 2024 Table 2: IIV on CL/F DEAQ = 19.7% CV (%RSE 37.0; bootstrap
    # median 19.2, 95% CI 6.2-36.2; eta shrinkage 70.4%)

    etalvc_deaq ~ 1.6491393
    # Ding 2024 Table 2: IIV on Vc/F DEAQ = 205% CV (%RSE 13.4; bootstrap
    # median 191, 95% CI 103-318; eta shrinkage 74.3%)

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
    propSd <- sqrt(0.267)
    label("Proportional residual SD for amodiaquine plasma concentration (SD on log scale)")
    # Ding 2024 Table 2: RUV AQ = 0.267 (variance; %RSE 13.9; bootstrap
    # median 0.262, 95% CI 0.200-0.338; epsilon shrinkage 4.6%)

    propSd_deaq <- sqrt(0.122)
    label("Proportional residual SD for desethylamodiaquine plasma concentration (SD on log scale)")
    # Ding 2024 Table 2: RUV DEAQ = 0.122 (variance; %RSE 10.7; bootstrap
    # median 0.122, 95% CI 0.096-0.148; epsilon shrinkage 22.4%)
  })

  model({
    # Molecular weights of the free bases (g/mol). Ding 2024 Methods
    # ('Population PK analysis') states the parent-metabolite model
    # assumes "complete bioconversion of amodiaquine into
    # desethylamodiaquine" but does not print the conversion factor. The
    # mass flux leaving amodiaquine central is therefore molar-corrected
    # before it enters desethylamodiaquine central, matching the sibling
    # WWARN amodiaquine model Ali_2018_amodiaquine.R (reference 29 of this
    # paper), which states the molar correction explicitly. The paper's
    # own secondary parameters confirm the correction: the printed
    # AUC_DEAQ / AUC_AQ ratio of 45,200 / 283 = 159.7 back-solves a
    # conversion factor of 0.902 (this ratio cancels the dose and F
    # entirely), against 0.921 for the molar reading and 1.0 for a
    # mass-for-mass reading.
    mwAQ        <- 355.85
    mwDEAQ      <- 327.81
    molarFactor <- mwDEAQ / mwAQ

    # Absorption chain rate constants. NN = 2 transit compartments and
    # MTT covers only those two compartments (see the lmtt annotation in
    # ini()), so ktr = NN/MTT; ka governs the final transfer into central.
    ka  <- exp(lka)
    mtt <- exp(lmtt)
    ktr <- 2 / mtt

    # Individual PK parameters. Allometric weight scaling on all apparent
    # clearances (exponent 0.75) and apparent volumes (exponent 1)
    # centred on the 70 kg reference of the Table 2 footnote. IIV is
    # carried only on F, CL/F DEAQ and Vc/F DEAQ (Ding 2024 Table 2).
    cl <- exp(lcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc) * (WT / 70)^e_wt_vc

    cl_deaq <- exp(lcl_deaq + etalcl_deaq) * (WT / 70)^e_wt_cl
    vc_deaq <- exp(lvc_deaq + etalvc_deaq) * (WT / 70)^e_wt_vc
    q_deaq  <- exp(lq_deaq)                * (WT / 70)^e_wt_cl
    vp_deaq <- exp(lvp_deaq)               * (WT / 70)^e_wt_vc

    # Micro-rate constants (1/h).
    kel_aq   <- cl      / vc
    kel_deaq <- cl_deaq / vc_deaq
    k12_deaq <- q_deaq  / vc_deaq
    k21_deaq <- q_deaq  / vp_deaq

    # ODE system (Ding 2024 Figure S1). Compartment amounts are in mg of
    # analyte base and volumes are in L, so amount/volume is mg/L and is
    # scaled to ng/mL below.
    #
    # Absorption: depot -> transit1 -> transit2 at ktr, then transit2 ->
    # amodiaquine central at ka.
    d/dt(depot)    <- -ktr * depot
    d/dt(transit1) <-  ktr * depot    - ktr * transit1
    d/dt(transit2) <-  ktr * transit1 - ka  * transit2

    # Amodiaquine central, one compartment. The entire mass flux leaving
    # amodiaquine central (kel_aq * central) is routed to
    # desethylamodiaquine central under the complete-bioconversion
    # assumption, after the molar correction.
    d/dt(central) <- ka * transit2 - kel_aq * central

    # Desethylamodiaquine central plus one peripheral compartment.
    d/dt(central_deaq)     <-  molarFactor * kel_aq * central -
                               kel_deaq * central_deaq -
                               k12_deaq * central_deaq + k21_deaq * peripheral1_deaq
    d/dt(peripheral1_deaq) <-  k12_deaq * central_deaq - k21_deaq * peripheral1_deaq

    # Relative oral bioavailability on the amodiaquine dose. The
    # population anchor lfdepot = log(1) is fixed; the gestational-age
    # effect is the linear deviation form of the Ding 2024 Table 2
    # footnote, F = 1 + theta * (GA - 24), referenced at the cohort
    # median of 24 weeks. IIV is carried by etalfdepot.
    f_ega    <- 1 + e_ega_f * (EGA - 24)
    f(depot) <- f_ega * exp(lfdepot + etalfdepot)

    # Plasma concentrations in ng/mL: amount (mg) / volume (L) is mg/L,
    # multiplied by 1000 to give ng/mL, the units used throughout Ding
    # 2024 Table 2 (day-7 concentration) and Figures 1 and 3.
    Cc      <- 1000 * central      / vc
    Cc_deaq <- 1000 * central_deaq / vc_deaq

    # Proportional residual error on the linear-concentration scale, the
    # linear-space equivalent of the paper's additive-on-log-scale error.
    Cc      ~ prop(propSd)
    Cc_deaq ~ prop(propSd_deaq)
  })
}
