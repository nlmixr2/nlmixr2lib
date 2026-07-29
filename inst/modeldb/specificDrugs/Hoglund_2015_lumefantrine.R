# Joint parent-metabolite population PK model of oral lumefantrine and its
# major oxidative metabolite desbutyl-lumefantrine in HIV-infected Ugandan
# adults co-administered antiretroviral regimens of efavirenz, nevirapine,
# or lopinavir/ritonavir (Hoglund 2015, Br J Clin Pharmacol 79:636-649;
# doi:10.1111/bcp.12529).

Hoglund_2015_lumefantrine <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for oral lumefantrine and",
    "its major oxidative metabolite desbutyl-lumefantrine in 89 HIV-infected",
    "Ugandan adults receiving artemether-lumefantrine (Coartem) with or",
    "without concomitant antiretroviral therapy (efavirenz, nevirapine, or",
    "lopinavir/ritonavir) (Hoglund 2015). 1-transit-compartment absorption",
    "with ka = ktr feeds a 2-compartment lumefantrine disposition; complete",
    "in-vivo conversion of lumefantrine to a 1-compartment",
    "desbutyl-lumefantrine disposition with stoichiometric molar conversion.",
    "Relative bioavailability F is anchored at 1 (fixed) with log-normal",
    "IIV (47.4 % CV). Three antiretroviral drug-drug interactions are",
    "encoded as linear-deviation effects on parent clearance and on",
    "bioavailability: efavirenz increases LF CL/F by 72.6 %,",
    "lopinavir/ritonavir decreases LF CL/F by 62.1 % and increases",
    "desbutyl-lumefantrine CL/F by 392 %, nevirapine decreases relative",
    "bioavailability by 24.8 %. IIV is retained on LF CL, the mean transit",
    "time, and the relative bioavailability F. NONMEM additive residual",
    "error on log-transformed concentrations is encoded as a proportional",
    "residual in linear concentration space for both parent and metabolite.",
    sep = " "
  )
  reference <- paste(
    "Hoglund RM, Byakika-Kibwika P, Lamorde M, Merry C, Ashton M,",
    "Hanpithakpong W, Day NPJ, White NJ, Abelo A, Tarning J (2015).",
    "Artemether-lumefantrine co-administration with antiretrovirals:",
    "population pharmacokinetics and dosing implications.",
    "British Journal of Clinical Pharmacology 79(4):636-649.",
    "doi:10.1111/bcp.12529.",
    sep = " "
  )
  vignette <- "Hoglund_2015_artemether_lumefantrine"
  units    <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CONMED_EFV = list(
      description        = "Concomitant efavirenz co-administration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "1 = subject is receiving concomitant efavirenz (600 mg once",
        "daily) based antiretroviral therapy, 0 = no concomitant",
        "efavirenz. Time-fixed per subject in the Hoglund 2015 cohort",
        "(Methods: efavirenz arm was treated for at least 1 month before",
        "the artemether-lumefantrine dosing period). Efavirenz is a",
        "CYP3A4 inducer; the indicator enters as a linear-deviation",
        "multiplier on apparent oral lumefantrine clearance: CL/F = TVCL",
        "* (1 + e_efv_cl * CONMED_EFV) with e_efv_cl = +0.726 (Table 2:",
        "EFZ CL/F = +72.6 % +/-17.2 % RSE)."
      ),
      source_name        = "EFZ"
    ),
    CONMED_NVP = list(
      description        = "Concomitant nevirapine co-administration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "1 = subject is receiving concomitant nevirapine (200 mg twice",
        "daily) based antiretroviral therapy, 0 = no concomitant",
        "nevirapine. Time-fixed per subject. Nevirapine is a CYP3A4",
        "inducer and a putative inducer of intestinal P-glycoprotein;",
        "the indicator enters as a linear-deviation multiplier on",
        "relative bioavailability F: f_typ = TVF1 * (1 + e_nvp_fdepot *",
        "CONMED_NVP) with e_nvp_fdepot = -0.248 (Table 2: NEV F =",
        "-24.8 % +/-38.6 % RSE)."
      ),
      source_name        = "NEV"
    ),
    CONMED_LPV = list(
      description        = "Concomitant lopinavir/ritonavir co-administration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "1 = subject is receiving concomitant lopinavir 400 mg /",
        "ritonavir 100 mg twice daily (Aluvia, Abbott); 0 = no",
        "concomitant LPV/r. Time-fixed per subject (Methods: study 1",
        "patients had been on LPV/r for at least 1 month before the",
        "artemether-lumefantrine dose). Lopinavir/ritonavir is a CYP3A4",
        "inhibitor (the ritonavir booster being the principal",
        "perpetrator); the indicator enters as a linear-deviation",
        "multiplier on apparent parent lumefantrine clearance: CL/F =",
        "TVCL * (1 + e_lpv_cl * CONMED_LPV) with e_lpv_cl = -0.621",
        "(Table 2: LOP CL/F = -62.1 % +/-8.48 % RSE). The same",
        "indicator additionally enters as a linear-deviation multiplier",
        "on apparent desbutyl-lumefantrine clearance: CL/F_desbutlum =",
        "TVCL_desbutlum * (1 + e_lpv_cl_desbutlum * CONMED_LPV) with",
        "e_lpv_cl_desbutlum = +3.92 (Table 2: LOP CL/F = +392 % on",
        "desbutyl-lumefantrine, +/-17.6 % RSE)."
      ),
      source_name        = "LOP"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 89L,
    n_studies       = 2L,
    n_study1        = 31L,
    n_study2        = 58L,
    age_range       = "20-70 years (study 1 median 36.5, range 24-51; study 2 median 36, range 20-70; Table 1)",
    weight_range    = "42-91 kg (study 1 median 64, range 45-86; study 2 median 56, range 42-91; Table 1)",
    bmi_range       = "17.0-36.5 kg/m^2 (study 1 median 23.7, range 17.0-34.0; study 2 median 22.3, range 17.3-36.5; Table 1)",
    sex_female_pct  = 73.0,
    disease_state   = paste(
      "HIV-infected adults without active malaria. Exclusion criteria:",
      "anaemia (Hb < 8 g/dL), pregnancy, abnormal liver / renal",
      "function, CYP- or P-gp-inhibiting/inducing comedications, herbal",
      "medicines, prolonged QT interval, intercurrent illness; study 1",
      "additionally excluded viral loads > 400 counts/mL and malaria",
      "parasitaemia."
    ),
    dose_range      = paste(
      "Coartem (Novartis, batch F0660): 20 mg artemether + 120 mg",
      "lumefantrine per tablet; oral 4-tablet dose (80 mg AM + 480 mg",
      "lumefantrine). Study 1 (parallel): single dose with venous",
      "plasma samples at 1, 2, 4, 6, 8, 12, 24, 48, 72 h post-dose;",
      "two arms (no HIV therapy vs LPV/r at least 1 month before",
      "study). Study 2 (3-period crossover): six-dose regimen (4",
      "tablets BID for 3 days), samples after the last dose at 1, 2,",
      "4, 8, 12, 24, 48, 72, 96, 120 h; period 1 antimalarial alone,",
      "1-week washout, HIV regimen alone for 1 month (NVP 200 mg or",
      "EFV 600 mg with NRTI backbone), then six-dose AM-LF",
      "concomitantly with HIV therapy in period 3."
    ),
    regions         = "Uganda (Mulago National Referral Hospital, Kampala)",
    trial_registration = "ClinicalTrials.gov NCT00619944 (study 1) and NCT00620438 (study 2)",
    notes           = paste(
      "Demographics pooled from Hoglund 2015 Table 1. Total observations:",
      "1365 lumefantrine + 341 desbutyl-lumefantrine samples; 99",
      "lumefantrine (7.3 %) and 37 desbutyl-lumefantrine (11 %) samples",
      "below limits of quantification (25 and 1.0 ng/mL respectively)",
      "were omitted. No body weight, BMI, age, sex, or study covariate",
      "was retained in the final model (Results, Lumefantrine and",
      "desbutyl-lumefantrine pharmacokinetics). Companion artemether",
      "extraction from the same cohort:",
      "modellib('Hoglund_2015_artemether')."
    )
  )

  ini({
    # Structural population mean parameters from Hoglund 2015 Table 2,
    # "Population estimates*" column. Values are apparent (relative to
    # F = 1) and given on the linear scale; log() applied here for the
    # nlmixr2 internal log-scale.
    lcl <- log(4.77)
    label("Apparent lumefantrine elimination clearance CL/F (L/h)")
    # Hoglund 2015 Table 2: CL/F = 4.77 L/h (RSE 5.30 %; 95 % CI 4.31-5.30)

    lvc <- log(68.9)
    label("Apparent lumefantrine central volume of distribution Vc/F (L)")
    # Hoglund 2015 Table 2: Vc/F = 68.9 L (RSE 27.1 %; 95 % CI 47.4-117)

    lq  <- log(2.86)
    label("Apparent lumefantrine inter-compartmental clearance Q/F (L/h)")
    # Hoglund 2015 Table 2: Q/F = 2.86 L/h (RSE 19.4 %; 95 % CI 1.72-3.62)

    lvp <- log(111)
    label("Apparent lumefantrine peripheral volume of distribution Vp/F (L)")
    # Hoglund 2015 Table 2: Vp/F = 111 L (RSE 9.14 %; 95 % CI 93.9-132)

    lmtt <- log(6.27)
    label("Mean transit time of the 1-transit-compartment absorption chain MTT (h)")
    # Hoglund 2015 Table 2: MTT = 6.27 h (RSE 21.2 %; 95 % CI 3.75-8.35).
    # Number of transit compartments = 1 (fixed; Table 2: "Number of
    # trans comp = 1 FIX"). With ka = ktr (Results: "k_a and k_tr
    # assumed to be equal"), the absorption chain depot -> transit1 ->
    # central has NN + 1 = 2 equal-rate transitions; ktr = (NN + 1)/MTT
    # = 2/MTT (Savic & Karlsson 2007 convention; same idiom as the
    # companion Hoglund_2018_mefloquine.R, Hoglund_2012_piperaquine.R,
    # Hoglund_2017_piperaquine.R files).

    lfdepot <- fixed(log(1))
    label("Relative bioavailability F (unitless; fixed at 1)")
    # Hoglund 2015 Table 2: Bioavailability = 1 FIX (Methods: "A fixed
    # bioavailability of 100% for the population with an estimated
    # between subject variability was evaluated").

    # Desbutyl-lumefantrine disposition (1-compartment with full
    # LF-to-desbutyl-lumefantrine conversion). Values from Hoglund 2015
    # Table 2 "Desbutyl-Lumefantrine" rows.
    lcl_desbutlum <- log(489)
    label("Apparent desbutyl-lumefantrine elimination clearance CL/F_desbutlum (L/h)")
    # Hoglund 2015 Table 2: CL/F = 489 L/h (RSE 5.98 %; 95 % CI 435-554)

    lvc_desbutlum <- log(22800)
    label("Apparent desbutyl-lumefantrine central volume Vc/F_desbutlum (L)")
    # Hoglund 2015 Table 2: Vc/F = 22 800 L (RSE 7.93 %; 95 % CI 19 600-26 800)

    # Drug-drug interaction effects on lumefantrine clearance and
    # bioavailability. The paper reports each effect as a percentage
    # change relative to the no-comedication reference. They enter the
    # model on the linear-deviation form CL = TVCL * (1 + e * X) so the
    # numeric coefficient is the fractional change (e.g. +72.6 % -> e =
    # +0.726, -62.1 % -> e = -0.621).
    e_efv_cl <- 0.726
    label("Linear-deviation effect of efavirenz on lumefantrine CL/F (fractional)")
    # Hoglund 2015 Table 2: EFZ CL/F = +72.6 % (RSE 17.2 %; 95 % CI 51.5-100)

    e_lpv_cl <- -0.621
    label("Linear-deviation effect of lopinavir/ritonavir on lumefantrine CL/F (fractional)")
    # Hoglund 2015 Table 2: LOP CL/F = -62.1 % (RSE 8.48 %; 95 % CI -72.1 to -51.8)

    e_nvp_fdepot <- -0.248
    label("Linear-deviation effect of nevirapine on relative bioavailability F (fractional)")
    # Hoglund 2015 Table 2: NEV F = -24.8 % (RSE 38.6 %; 95 % CI -42.4 to -4.66)

    e_lpv_cl_desbutlum <- 3.92
    label("Linear-deviation effect of lopinavir/ritonavir on desbutyl-lumefantrine CL/F (fractional)")
    # Hoglund 2015 Table 2: LOP CL/F = +392 % on desbutyl-lumefantrine (RSE 17.6 %; 95 % CI 239-488)

    # IIV. Hoglund 2015 Table 2 footnote: "Coefficients of variation
    # (%CV) for between-subject variability (BSV) were calculated as
    # 100 * (e^variance - 1)^(1/2)", so the internal log-scale variance
    # is recovered as omega^2 = log((CV/100)^2 + 1).
    #   CL  IIV 14.8 % -> omega^2 = log(0.148^2 + 1) = 0.021668
    #   MTT IIV 31.4 % -> omega^2 = log(0.314^2 + 1) = 0.094033
    #   F   IIV 47.4 % -> omega^2 = log(0.474^2 + 1) = 0.202676
    # Only CL, MTT, and F carry IIV in the final model (Results,
    # Lumefantrine and desbutyl-lumefantrine pharmacokinetics:
    # "Between subject variability was retained in the estimates of
    # elimination clearance, the mean transit time and the relative
    # bioavailability of lumefantrine").
    etalcl     ~ 0.021668
    # Hoglund 2015 Table 2: IIV on CL/F = 14.8 % CV (RSE 50.7 %; 95 % CI 5.68-22.7; shrinkage 40.1 %)

    etalmtt    ~ 0.094033
    # Hoglund 2015 Table 2: IIV on MTT = 31.4 % CV (RSE 106 %; 95 % CI 15.1-94.6; shrinkage 46.9 %)

    etalfdepot ~ 0.202676
    # Hoglund 2015 Table 2: IIV on F = 47.4 % CV (RSE 18.3 %; 95 % CI 38.2-57.4; shrinkage 6.21 %)

    # Residual error. Methods, Pharmacokinetic analysis: "All
    # concentrations were converted into their natural logarithms ...
    # An additive residual error on log-transformed data (essentially
    # equivalent to an exponential error on normal scale data) was
    # used." This NONMEM additive-on-log-scale residual maps to a
    # nlmixr2 proportional residual in linear concentration space (see
    # references/parameter-names.md Residual error and the matching
    # comment in Hoglund_2018_mefloquine.R / Hoglund_2017_piperaquine.R).
    # Table 2 reports the log-scale variance as "RUV"; the SD is
    # sqrt(RUV), which equals the proportional CV to first order.
    #   Lumefantrine          RUV = 0.566 -> SD = sqrt(0.566) = 0.7523
    #   Desbutyl-lumefantrine RUV = 0.465 -> SD = sqrt(0.465) = 0.6819
    propSd <- sqrt(0.566)
    label("Proportional residual SD for lumefantrine plasma concentration (SD on log scale)")
    # Hoglund 2015 Table 2: RUV (LF) = 0.566 (variance on log scale, RSE 7.83 %; 95 % CI 0.479-0.643; epsilon shrinkage 4.05 %)

    propSd_desbutlum <- sqrt(0.465)
    label("Proportional residual SD for desbutyl-lumefantrine plasma concentration (SD on log scale)")
    # Hoglund 2015 Table 2: RUV (desbutyl-lumefantrine) = 0.465 (variance on log scale, RSE 13.3 %; 95 % CI 0.357-0.591; epsilon shrinkage 6.05 %)
  })

  model({
    # Molecular weights (g/mol). Lumefantrine C30H32Cl3NO = 528.94 g/mol;
    # desbutyl-lumefantrine C26H24Cl3NO = 472.83 g/mol. Used at the
    # systemic LF -> desbutyl-lumefantrine formation step to convert
    # mass-rate of lumefantrine elimination into mass-rate of metabolite
    # formation under the source paper's assumption of complete in-vivo
    # 1:1 molar conversion (Methods: "Both artemether and lumefantrine
    # were assumed to be fully converted into their respective
    # metabolite").
    mw_lf        <- 528.94
    mw_desbutlum <- 472.83

    # Individual PK parameters with DDI covariate effects. CL_LF carries
    # both the EFV and LPV/r linear-deviation effects (mutually exclusive
    # in the trial design, so at most one of CONMED_EFV / CONMED_LPV is
    # 1 at any time; the multiplicative form is the canonical encoding
    # for independently estimated categorical effects). CL_desbutlum
    # carries only the LPV/r effect (Hoglund 2015 retained no EFV or NVP
    # effect on desbutyl-lumefantrine clearance). Vc, Q, Vp do not carry
    # IIV in the final model (Table 2 IIV columns "-").
    cl <- exp(lcl + etalcl) *
            (1 + e_efv_cl * CONMED_EFV) *
            (1 + e_lpv_cl * CONMED_LPV)
    vc <- exp(lvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    cl_desbutlum <- exp(lcl_desbutlum) *
                      (1 + e_lpv_cl_desbutlum * CONMED_LPV)
    vc_desbutlum <- exp(lvc_desbutlum)

    # Typical relative bioavailability with the NVP linear-deviation
    # effect; log-normal IIV applied multiplicatively via f(depot)
    # below. Anchor lfdepot = log(1) is fixed (paper convention).
    fdepot_typ <- exp(lfdepot) * (1 + e_nvp_fdepot * CONMED_NVP)

    # Mean transit time and chain rate constant. With NN = 1 transit
    # compartment and ka = ktr, the absorption chain depot -> transit1
    # -> central has NN + 1 = 2 equal-rate transitions;
    # ktr = (NN + 1)/MTT = 2/MTT (Savic & Karlsson 2007).
    mtt <- exp(lmtt + etalmtt)
    ktr <- 2 / mtt

    # Disposition micro-constants.
    kel          <- cl           / vc
    k12          <- q            / vc
    k21          <- q            / vp
    kel_desbutlum <- cl_desbutlum / vc_desbutlum

    # ODE system: 1-transit-compartment absorption feeding a
    # 2-compartment lumefantrine disposition; full conversion of
    # lumefantrine to a 1-compartment desbutyl-lumefantrine
    # disposition with stoichiometric molar conversion.
    d/dt(depot)            <- -ktr * depot
    d/dt(transit1)         <-  ktr * depot    - ktr * transit1
    d/dt(central)          <-  ktr * transit1 - kel * central -
                                k12 * central + k21 * peripheral1
    d/dt(peripheral1)      <-  k12 * central  - k21 * peripheral1
    d/dt(central_desbutlum) <- kel * central  * (mw_desbutlum / mw_lf) -
                                kel_desbutlum * central_desbutlum

    # Relative bioavailability applied to the depot. NVP effect via
    # fdepot_typ; log-normal IIV via etalfdepot.
    f(depot) <- fdepot_typ * exp(etalfdepot)

    # Plasma concentrations in ng/mL. Dose in mg, Vc in L gives
    # central / vc in mg/L = ug/mL; multiplying by 1000 yields ng/mL,
    # the unit Hoglund 2015 reports throughout (e.g. LLOQ 25 ng/mL for
    # lumefantrine and 1.0 ng/mL for desbutyl-lumefantrine).
    Cc            <- 1000 * central            / vc
    Cc_desbutlum  <- 1000 * central_desbutlum  / vc_desbutlum

    # Proportional residual error on the linear-concentration scale
    # (NONMEM additive-on-log-scale -> nlmixr2 proportional; see ini()
    # comment on propSd / propSd_desbutlum).
    Cc           ~ prop(propSd)
    Cc_desbutlum ~ prop(propSd_desbutlum)
  })
}
