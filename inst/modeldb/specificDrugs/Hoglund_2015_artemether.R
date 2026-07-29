# Joint parent-metabolite population PK model of oral artemether and its
# active metabolite dihydroartemisinin (DHA), with enzymatic auto-induction
# of artemether elimination, in HIV-infected Ugandan adults co-administered
# antiretroviral regimens of efavirenz, nevirapine, or lopinavir/ritonavir
# (Hoglund 2015, Br J Clin Pharmacol 79:636-649; doi:10.1111/bcp.12529).

Hoglund_2015_artemether <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for oral artemether and",
    "its active metabolite dihydroartemisinin (DHA) in 89 HIV-infected",
    "Ugandan adults receiving artemether-lumefantrine (Coartem) with or",
    "without concomitant antiretroviral therapy (efavirenz, nevirapine, or",
    "lopinavir/ritonavir) (Hoglund 2015). 3-transit-compartment absorption",
    "with ka = ktr feeds a 1-compartment artemether disposition; complete",
    "in-vivo conversion of artemether to a 1-compartment DHA disposition",
    "with stoichiometric molar conversion. Enzymatic auto-induction of the",
    "first-pass demethylation of artemether is modelled as a Hill function",
    "of time-since-first-dose with fixed maximum maturation (100 % increase",
    "in apparent CL) and fixed maturation half-time (62 h, literature",
    "value) and an estimated Hill coefficient (0.445). Relative",
    "bioavailability F is anchored at 1 (fixed) with log-normal IIV",
    "(58.6 % CV). Three antiretroviral drug-drug interactions are encoded",
    "as linear-deviation effects: lopinavir/ritonavir increases AM CL/F by",
    "32.8 % and DHA CL/F by 143 %; efavirenz and nevirapine decrease",
    "relative bioavailability by 71.5 % and 66.3 % respectively; nevirapine",
    "additionally decreases DHA CL/F by 44.5 %. IIV is retained on AM CL,",
    "DHA CL, the mean transit time, and F. NONMEM additive residual error",
    "on log-transformed concentrations is encoded as a proportional",
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
        "(Methods: HIV regimen administered alone for 1 month before",
        "the second artemether-lumefantrine dosing period). Efavirenz",
        "enters this model as a linear-deviation multiplier on relative",
        "bioavailability F: f_typ = TVF1 * (1 + e_efv_fdepot *",
        "CONMED_EFV) with e_efv_fdepot = -0.715 (Table 3: EFZ F =",
        "-71.5 % +/-5.94 % RSE; mechanism: induction of intestinal",
        "P-glycoprotein and/or intestinal CYP3A4)."
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
        "nevirapine. Time-fixed per subject. Nevirapine enters this",
        "model as TWO independent linear-deviation effects: (i) on",
        "relative bioavailability F: f_typ = TVF1 * (1 + e_nvp_fdepot *",
        "CONMED_NVP) with e_nvp_fdepot = -0.663 (Table 3: NEV F =",
        "-66.3 % +/-7.41 % RSE); and (ii) on apparent dihydroartemisinin",
        "clearance: CL_dha/F = TVCL_dha * (1 + e_nvp_cl_dha *",
        "CONMED_NVP) with e_nvp_cl_dha = -0.445 (Table 3: NEV CL/F on",
        "DHA = -44.5 % +/-14.5 % RSE; the negative direction is",
        "physiologically unexpected since NVP has not been reported to",
        "affect the UGT system that clears DHA, and is flagged as such",
        "in the source paper Discussion)."
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
        "concomitant LPV/r. Time-fixed per subject. Lopinavir/ritonavir",
        "enters this model as TWO independent linear-deviation",
        "multipliers: (i) on apparent artemether clearance: CL/F = TVCL",
        "* (1 + e_lpv_cl * CONMED_LPV) with e_lpv_cl = +0.328 (Table 3:",
        "LOP CL/F = +32.8 % +/-21.5 % RSE); and (ii) on apparent",
        "dihydroartemisinin clearance: CL_dha/F = TVCL_dha * (1 +",
        "e_lpv_cl_dha * CONMED_LPV) with e_lpv_cl_dha = +1.43 (Table 3:",
        "LOP CL/F on DHA = +143 % +/-19.7 % RSE; mechanism: lopinavir",
        "induction of CYP2B6/2C9/2C19 and possibly other unknown",
        "pathways)."
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
      "lumefantrine per tablet; oral 4-tablet dose (80 mg artemether +",
      "480 mg lumefantrine). Study 1 (parallel): single dose; study 2",
      "(crossover): standard 6-dose regimen (4 tablets BID for 3",
      "days). Samples for AM/DHA: 1, 2, 4, 6, 8, 12, 24, 48, 72 h",
      "post-dose (study 1); after the last dose of the 6-dose regimen",
      "at 1, 2, 4, 8, 12, 24, 48, 72, 96, 120 h (study 2)."
    ),
    regions         = "Uganda (Mulago National Referral Hospital, Kampala)",
    trial_registration = "ClinicalTrials.gov NCT00619944 (study 1) and NCT00620438 (study 2)",
    notes           = paste(
      "Demographics pooled from Hoglund 2015 Table 1. Total observations:",
      "934 AM + DHA plasma samples; 277 AM (30 %) and 181 DHA (19 %)",
      "below LLOQ (1.4 ng/mL for both) were omitted (the M3 method was",
      "attempted but produced unrealistic estimates). No body weight,",
      "BMI, age, sex, or study covariate was retained in the final model",
      "(Results, Artemether and dihydroartemisinin pharmacokinetics).",
      "Companion lumefantrine extraction from the same cohort:",
      "modellib('Hoglund_2015_lumefantrine')."
    )
  )

  ini({
    # Structural population mean parameters from Hoglund 2015 Table 3,
    # "Population estimates*" column. Values are apparent (relative to
    # F = 1) and given on the linear scale; log() applied here for the
    # nlmixr2 internal log-scale.
    lcl <- log(317)
    label("Apparent baseline artemether elimination clearance CL/F (L/h)")
    # Hoglund 2015 Table 3: CL/F = 317 L/h (RSE 8.25 %; 95 % CI 270-374).
    # Baseline CL/F before auto-induction maturation; effective CL/F
    # rises with time-since-first-dose according to the Hill maturation
    # function defined by (mat_mmax, mat_mf50, mat_hill) below.

    lvc <- log(1090)
    label("Apparent artemether central volume of distribution Vc/F (L)")
    # Hoglund 2015 Table 3: Vc/F = 1090 L (RSE 8.69 %; 95 % CI 917-1291)

    lmtt <- log(0.970)
    label("Mean transit time of the 3-transit-compartment absorption chain MTT (h)")
    # Hoglund 2015 Table 3: MTT = 0.970 h (RSE 6.58 %; 95 % CI 0.853-1.10).
    # Number of transit compartments = 3 (fixed; Table 3: "Number of
    # trans comp = 3 FIX"). With ka = ktr (Results, Artemether and
    # dihydroartemisinin pharmacokinetics: "k_a and k_tr assumed to be
    # equal"), the absorption chain depot -> transit1 -> transit2 ->
    # transit3 -> central has NN + 1 = 4 equal-rate transitions;
    # ktr = (NN + 1)/MTT = 4/MTT (Savic & Karlsson 2007 convention).

    lfdepot <- fixed(log(1))
    label("Relative bioavailability F (unitless; fixed at 1)")
    # Hoglund 2015 Table 3: Bioavailability = 1 FIX (Methods: "A fixed
    # bioavailability of 100 % for the population with an estimated
    # between subject variability was evaluated").

    # Auto-induction (maturation) parameters. Methods: "Enzymatic
    # auto-induction of artemether metabolism was evaluated ... Auto-
    # induction was modelled using a maturation-model with a maturation
    # half-life fixed to a mean of literature values (i.e. 62 h) and
    # the maximum maturation were assumed to be 100 %." The functional
    # form is the Hill maturation
    #   maturation(t) = t^Hill / (MF50^Hill + t^Hill)
    #   CL_AM(t) = CL_AM_baseline * (1 + MMAX * maturation(t))
    # with t = time since first dose. MF50 is the time at which the
    # maturation factor reaches 50 % (Table 3 caption: "MF50 is the
    # time in which the maturation has reached 50 %"). MMAX is the
    # maximum fractional increment in CL_AM at full induction (Table 3
    # caption: "MMAX is the maximum maturation in the maturation
    # model"). Both MMAX and MF50 are FIXED; only Hill is estimated.
    mat_mmax <- fixed(1)
    label("Maximum auto-induction increment in artemether CL/F (fraction; fixed at 100 %)")
    # Hoglund 2015 Table 3: MMAX = 1 FIX

    mat_mf50 <- fixed(62)
    label("Auto-induction half-time on artemether CL/F: time at which maturation reaches 50 % (h; fixed at literature value)")
    # Hoglund 2015 Table 3: MF50 = 62 FIX (literature value, Methods)

    mat_hill <- 0.445
    label("Hill coefficient governing the slope of the artemether auto-induction maturation function (unitless)")
    # Hoglund 2015 Table 3: Hill = 0.445 (RSE 43.9 % based on 1000
    # separate bootstrap runs with a reduced dataset; the main bootstrap
    # did not provide a CI for this parameter)

    # Dihydroartemisinin disposition (1-compartment with full AM-to-DHA
    # conversion). Values from Hoglund 2015 Table 3 "Dihydroartemisinin"
    # rows.
    lcl_dha <- log(160)
    label("Apparent dihydroartemisinin elimination clearance CL/F_dha (L/h)")
    # Hoglund 2015 Table 3: CL/F = 160 L/h (RSE 4.93 %; 95 % CI 145-174)

    lvc_dha <- log(14.9)
    label("Apparent dihydroartemisinin central volume Vc/F_dha (L)")
    # Hoglund 2015 Table 3: Vc/F = 14.9 L (RSE 39.8 %; 95 % CI 4.22-27.9)

    # Drug-drug interaction effects. Each effect enters on the linear-
    # deviation form parameter = TV * (1 + e * X) with the numeric
    # coefficient equal to the fractional change relative to the
    # no-comedication reference.
    e_lpv_cl <- 0.328
    label("Linear-deviation effect of lopinavir/ritonavir on artemether CL/F (fractional)")
    # Hoglund 2015 Table 3: LOP CL/F = +32.8 % (RSE 21.5 %; 95 % CI 21.0-47.0)

    e_efv_fdepot <- -0.715
    label("Linear-deviation effect of efavirenz on relative bioavailability F (fractional)")
    # Hoglund 2015 Table 3: EFZ F = -71.5 % (RSE 5.94 %; 95 % CI -79.3 to -62.0)

    e_nvp_fdepot <- -0.663
    label("Linear-deviation effect of nevirapine on relative bioavailability F (fractional)")
    # Hoglund 2015 Table 3: NEV F = -66.3 % (RSE 7.41 %; 95 % CI -75.3 to -55.9)

    e_lpv_cl_dha <- 1.43
    label("Linear-deviation effect of lopinavir/ritonavir on dihydroartemisinin CL/F (fractional)")
    # Hoglund 2015 Table 3: LOP CL/F on DHA = +143 % (RSE 19.7 %; 95 % CI 96.2-207)

    e_nvp_cl_dha <- -0.445
    label("Linear-deviation effect of nevirapine on dihydroartemisinin CL/F (fractional)")
    # Hoglund 2015 Table 3: NEV CL/F on DHA = -44.5 % (RSE 14.5 %; 95 % CI -56.6 to -31.2)

    # IIV. Hoglund 2015 Table 3 footnote: "Coefficients of variation
    # (%CV) for between-subject variability (BSV) were calculated as
    # 100 * (e^variance - 1)^(1/2)", so the internal log-scale variance
    # is recovered as omega^2 = log((CV/100)^2 + 1).
    #   CL_AM IIV  9.8 % -> omega^2 = log(0.098^2 + 1) = 0.009558
    #   MTT   IIV 51.6 % -> omega^2 = log(0.516^2 + 1) = 0.236065
    #   F     IIV 58.6 % -> omega^2 = log(0.586^2 + 1) = 0.295201
    #   CL_DHA IIV 53.9 %-> omega^2 = log(0.539^2 + 1) = 0.255046
    # IIV is retained on AM CL, DHA CL, MTT, and F (Results: "Between
    # subject variability was retained in the estimates of elimination
    # clearance of artemether, elimination clearance of dihydroartemisinin,
    # mean transit time and relative bioavailability").
    etalcl     ~ 0.009558
    # Hoglund 2015 Table 3: IIV on AM CL/F = 9.8 % CV (RSE 60.6 %; 95 % CI 2.75-16.1; shrinkage 50.6 %)

    etalmtt    ~ 0.236065
    # Hoglund 2015 Table 3: IIV on MTT = 51.6 % CV (RSE 18.6 %; 95 % CI 41.1-61.6; shrinkage 19.0 %)

    etalfdepot ~ 0.295201
    # Hoglund 2015 Table 3: IIV on F = 58.6 % CV (RSE 27.6 %; 95 % CI 39.9-76.6; shrinkage 10.7 %)

    etalcl_dha ~ 0.255046
    # Hoglund 2015 Table 3: IIV on DHA CL/F = 53.9 % CV (RSE 31.5 %; 95 % CI 33.4-70.4; shrinkage 18.7 %)

    # Residual error. Methods, Pharmacokinetic analysis: "All
    # concentrations were converted into their natural logarithms ...
    # An additive residual error on log-transformed data (essentially
    # equivalent to an exponential error on normal scale data) was
    # used." This NONMEM additive-on-log-scale residual maps to a
    # nlmixr2 proportional residual in linear concentration space (see
    # references/parameter-names.md Residual error and the matching
    # comment in the companion Hoglund_2015_lumefantrine.R). Table 3
    # reports the log-scale variance as "RUV"; the SD is sqrt(RUV),
    # which equals the proportional CV to first order.
    #   Artemether         RUV = 0.724 -> SD = sqrt(0.724) = 0.8509
    #   Dihydroartemisinin RUV = 0.707 -> SD = sqrt(0.707) = 0.8408
    propSd <- sqrt(0.724)
    label("Proportional residual SD for artemether plasma concentration (SD on log scale)")
    # Hoglund 2015 Table 3: RUV (AM) = 0.724 (variance on log scale, RSE 5.60 %; 95 % CI 0.644-0.803; epsilon shrinkage 7.25 %)

    propSd_dha <- sqrt(0.707)
    label("Proportional residual SD for dihydroartemisinin plasma concentration (SD on log scale)")
    # Hoglund 2015 Table 3: RUV (DHA) = 0.707 (variance on log scale, RSE 4.40 %; 95 % CI 0.645-0.764; epsilon shrinkage 8.28 %)
  })

  model({
    # Molecular weights (g/mol). Artemether C16H26O5 = 298.37 g/mol;
    # dihydroartemisinin C15H24O5 = 284.35 g/mol. Used at the systemic
    # AM -> DHA formation step to convert mass-rate of artemether
    # elimination into mass-rate of DHA formation under the source
    # paper's assumption of complete in-vivo 1:1 molar conversion
    # (Methods: "Both artemether and lumefantrine were assumed to be
    # fully converted into their respective metabolite").
    mw_am  <- 298.37
    mw_dha <- 284.35

    # Enzymatic auto-induction maturation factor on AM CL/F. The Hill
    # function increases from 0 at t = 0 (no induction) to 1 as
    # t -> infinity (full induction); MMAX scales the asymptotic
    # increment in apparent CL/F. With MMAX = 1, AM CL/F doubles at
    # full induction (Methods: "the maximum maturation were assumed to
    # be 100 %"). t is the model time since the first dose.
    maturation_am <- t^mat_hill / (mat_mf50^mat_hill + t^mat_hill)
    ind_factor    <- 1 + mat_mmax * maturation_am

    # Individual PK parameters with auto-induction and DDI covariate
    # effects. CL_AM carries the LPV/r linear-deviation effect; CL_DHA
    # carries the LPV/r and NVP linear-deviation effects. Vc_AM, Vc_DHA
    # do not carry IIV in the final model (Table 3 IIV columns "-").
    cl     <- exp(lcl     + etalcl)     * ind_factor *
                (1 + e_lpv_cl     * CONMED_LPV)
    vc     <- exp(lvc)
    cl_dha <- exp(lcl_dha + etalcl_dha) *
                (1 + e_lpv_cl_dha * CONMED_LPV) *
                (1 + e_nvp_cl_dha * CONMED_NVP)
    vc_dha <- exp(lvc_dha)

    # Typical relative bioavailability with EFV and NVP linear-deviation
    # effects; log-normal IIV applied multiplicatively via f(depot)
    # below. Anchor lfdepot = log(1) is fixed (paper convention).
    fdepot_typ <- exp(lfdepot) *
                    (1 + e_efv_fdepot * CONMED_EFV) *
                    (1 + e_nvp_fdepot * CONMED_NVP)

    # Mean transit time and chain rate constant. With NN = 3 transit
    # compartments and ka = ktr, the absorption chain depot -> transit1
    # -> transit2 -> transit3 -> central has NN + 1 = 4 equal-rate
    # transitions; ktr = (NN + 1)/MTT = 4/MTT (Savic & Karlsson 2007).
    mtt <- exp(lmtt + etalmtt)
    ktr <- 4 / mtt

    # Disposition micro-constants.
    kel     <- cl     / vc
    kel_dha <- cl_dha / vc_dha

    # ODE system: 3-transit-compartment absorption feeding a
    # 1-compartment artemether disposition; full conversion of
    # artemether to a 1-compartment dihydroartemisinin disposition
    # with stoichiometric molar conversion.
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(transit3)    <-  ktr * transit2 - ktr * transit3
    d/dt(central)     <-  ktr * transit3 - kel * central
    d/dt(central_dha) <-  kel * central  * (mw_dha / mw_am) -
                          kel_dha * central_dha

    # Relative bioavailability applied to the depot. EFV / NVP effects
    # via fdepot_typ; log-normal IIV via etalfdepot.
    f(depot) <- fdepot_typ * exp(etalfdepot)

    # Plasma concentrations in ng/mL. Dose in mg, Vc in L gives
    # central / vc in mg/L = ug/mL; multiplying by 1000 yields ng/mL,
    # the unit Hoglund 2015 reports throughout (e.g. LLOQ 1.4 ng/mL
    # for AM and DHA).
    Cc     <- 1000 * central     / vc
    Cc_dha <- 1000 * central_dha / vc_dha

    # Proportional residual error on the linear-concentration scale
    # (NONMEM additive-on-log-scale -> nlmixr2 proportional; see ini()
    # comment on propSd / propSd_dha).
    Cc     ~ prop(propSd)
    Cc_dha ~ prop(propSd_dha)
  })
}
