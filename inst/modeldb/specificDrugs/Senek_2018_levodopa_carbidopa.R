Senek_2018_levodopa_carbidopa <- function() {
  description <- "Joint population PK model for levodopa and carbidopa after a single dose of dispersible levodopa/carbidopa (4:1) microtablets in healthy adults and in Parkinson's disease patients with motor fluctuations: levodopa two-compartment and carbidopa one-compartment disposition, each analyte absorbed through two parallel transit-compartment chains (fast and slow) that reproduce double-peak profiles, an estimated pre-dose concentration in patients decaying at the terminal slope, allometric body-weight scaling, a hyperbolic carbidopa-dose inhibition of levodopa CL/F, Hoehn and Yahr stage on both CL/F, age on carbidopa CL/F and levodopa relative bioavailability, and sex on the carbidopa slow-chain mean transit time (Senek 2018)"
  reference <- paste(
    "Senek M, Nyholm D, Nielsen EI.",
    "Population pharmacokinetics of levodopa/carbidopa microtablets in healthy",
    "subjects and Parkinson's disease patients.",
    "Eur J Clin Pharmacol. 2018;74(10):1299-1307.",
    "doi:10.1007/s00228-018-2497-2"
  )
  vignette <- "Senek_2018_levodopa_carbidopa"
  units <- list(time = "min", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling of every disposition parameter with reference weight 70 kg and exponents fixed to 0.75 (CL/F, Q/F) and 1 (VC/F, VP/F) (Methods, 'Covariate model', first display equation).",
      source_name = "WT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Centred linear effect at the pooled-cohort mean of 48.7 years (Table 1 footnote) on carbidopa CL/F and on levodopa relative bioavailability: P = TVP * (1 + coeff * (AGE - 48.7)) (Methods, 'Covariate model', second display equation).",
      source_name = "Age"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "Centred at the pooled-cohort female fraction 0.389 (Table 1 footnote 'sex, 0.389'; 14 of 36 subjects female), so neither sex is the reference",
      notes = "Centred linear effect on the carbidopa slow-chain mean transit time: MTT2_CD = TVMTT2 * (1 + 0.178 * (SEXF - 0.389)). The 0.389 centring equals the female fraction, which pins the coding as 1 = female; the paper's 'MTT2,CD ... 19% (22 min) longer for women' reproduces under this coding.",
      source_name = "Sex"
    ),
    DIS_PD_HY = list(
      description = "Hoehn and Yahr stage of Parkinson's disease, extended to 0 for healthy subjects",
      units = "(ordinal stage 0-5; 0 = healthy)",
      type = "categorical",
      reference_category = "n/a -- centred at the pooled-cohort mean of 1.53 (Table 1 footnote 'HY, 1.53')",
      notes = "Centred linear effect on levodopa CL/F (coefficient -0.0616) and carbidopa CL/F (coefficient +0.0141): CL/F = TVCL * (1 + coeff * (DIS_PD_HY - 1.53)). Healthy subjects carry 0 ('Hoehn and Yahr score (HY, disease stage, set to 0 for healthy subjects)', Methods 'Covariate model'). Study 2 patients were HY 2-5 (mean 3.2).",
      source_name = "HY"
    ),
    DOSE_CARBIDOPA_MG = list(
      description = "Total carbidopa dose administered over the study period",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters levodopa CL/F hyperbolically: CL/F_LD = TVCL / (1 + DOSE_CARBIDOPA_MG / INTER), INTER = 86.7 mg (Methods third display equation; Table 2 'INTERCL/F-CDAMT'). The TOTAL carbidopa dose over the study period, not the per-administration dose: 25 mg in Study 1, 27.5-102.5 mg in Study 2, and 75 mg (= 18.75 + 5 x 11.25) for the six-dose Study 3 regimen (Results, 'the total CD dose administered during the study period was used to describe the CD/LD interaction'). Time-fixed per subject.",
      source_name = "CDAMT"
    ),
    DIS_HEALTHY = list(
      description = "Healthy-subject indicator, 1 = healthy volunteer, 0 = Parkinson's disease patient",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Parkinson's disease patient)",
      notes = "Selects the population-specific residual-error magnitudes (Table 2 'Healthy subjects' / 'Patients' rows, estimated separately per analyte) and switches off the pre-dose concentration term, which the paper introduces for the residual levodopa and carbidopa that Study 2 patients carried over despite overnight washout (Methods, 'Base model'). Study 1 subjects were drug-naive healthy volunteers.",
      source_name = "(study association: Study 1 healthy vs Study 2 patients)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "levodopa", units = "mg", specimen = "administration site", verified = TRUE),
    transit_fast1 = list(analyte = "levodopa", units = "mg", specimen = "administration site", verified = TRUE),
    transit_fast2 = list(analyte = "levodopa", units = "mg", specimen = "administration site", verified = TRUE),
    transit_fast3 = list(analyte = "levodopa", units = "mg", specimen = "administration site", verified = TRUE),
    transit_fast4 = list(analyte = "levodopa", units = "mg", specimen = "administration site", verified = TRUE),
    transit_fast5 = list(analyte = "levodopa", units = "mg", specimen = "administration site", verified = TRUE),
    depot2 = list(analyte = "levodopa", units = "mg", specimen = "administration site", verified = TRUE),
    transit_slow1 = list(analyte = "levodopa", units = "mg", specimen = "administration site", verified = TRUE),
    transit_slow2 = list(analyte = "levodopa", units = "mg", specimen = "administration site", verified = TRUE),
    transit_slow3 = list(analyte = "levodopa", units = "mg", specimen = "administration site", verified = TRUE),
    transit_slow4 = list(analyte = "levodopa", units = "mg", specimen = "administration site", verified = TRUE),
    transit_slow5 = list(analyte = "levodopa", units = "mg", specimen = "administration site", verified = TRUE),
    transit_slow6 = list(analyte = "levodopa", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "levodopa", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "levodopa", units = "mg", specimen = "plasma", verified = TRUE),
    depot_carbidopa = list(analyte = "carbidopa", units = "mg", specimen = "administration site", verified = TRUE),
    transit_fast1_carbidopa = list(
      analyte = "carbidopa",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit_fast2_carbidopa = list(
      analyte = "carbidopa",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit_fast3_carbidopa = list(
      analyte = "carbidopa",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    depot2_carbidopa = list(analyte = "carbidopa", units = "mg", specimen = "administration site", verified = TRUE),
    transit_slow1_carbidopa = list(
      analyte = "carbidopa",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit_slow2_carbidopa = list(
      analyte = "carbidopa",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit_slow3_carbidopa = list(
      analyte = "carbidopa",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit_slow4_carbidopa = list(
      analyte = "carbidopa",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit_slow5_carbidopa = list(
      analyte = "carbidopa",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit_slow6_carbidopa = list(
      analyte = "carbidopa",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit_slow7_carbidopa = list(
      analyte = "carbidopa",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit_slow8_carbidopa = list(
      analyte = "carbidopa",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit_slow9_carbidopa = list(
      analyte = "carbidopa",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit_slow10_carbidopa = list(
      analyte = "carbidopa",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central_carbidopa = list(analyte = "carbidopa", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 36L,
    n_studies = 2L,
    age_range = "19-82 years",
    age_mean = "48.7 years (healthy 26.0, patients 71.4)",
    weight_range = "55-96 kg",
    weight_mean = "71.7 kg (healthy), 75.4 kg (patients)",
    sex_female_pct = 38.9,
    disease_state = "18 healthy adult volunteers (Study 1) and 18 patients with Parkinson's disease experiencing motor fluctuations (Study 2; Hoehn and Yahr stage 2-5, mean 3.2; 9.6 years since diagnosis, 9.5 years on levodopa)",
    dose_range = "Single dose of levodopa/carbidopa 5/1.25 mg microtablets dispersed in 100 mL water: 100/25 mg in Study 1; 150% of the individual levodopa-equivalent morning dose, 110/27.5 to 410/102.5 mg, after overnight washout in Study 2",
    regions = "Sweden",
    notes = "Pooled from two single-centre, open-label, single-dose microtablet studies (Nyholm 2012 Clin Neuropharmacol, Senek 2017 Eur J Clin Pharmacol). Demographics from Table 1: Study 1 9 M / 9 F, age 26.0 +/- 6.2 [19-46] years, WT 71.7 +/- 11.3 [59-95] kg; Study 2 13 M / 5 F, age 71.4 +/- 6.3 [61-82] years, WT 75.4 +/- 11.0 [55-96] kg. One of 19 Study 2 patients was excluded for an extreme triple-peak profile. Study 3 (10 healthy volunteers, 4 M / 6 F, six doses 75/18.75 mg then 45/11.25 mg every 2.4 h; Nyholm 2013 Acta Neurol Scand) was used only for external evaluation of the levodopa model (Fig. 2). Levodopa and carbidopa below the limit of detection (5.6%) were handled with the M6 method."
  )

  ini({
    # Levodopa disposition (Table 2, 'Levodopa' column = NONMEM final estimates).
    # Typical CL/F is for 0 mg carbidopa and the mean Hoehn and Yahr stage 1.53 (Table 2 footnote dagger).
    lcl <- log(1.31); label("Levodopa CL/F at 70 kg, 0 mg carbidopa and Hoehn and Yahr 1.53 (L/min)") # Table 2 'CL/F (L/min)*(WT/70)^0.75' levodopa 1.31 (RSE 10.7%)
    lvc <- log(45.4); label("Levodopa central volume VC/F at 70 kg (L)") # Table 2 'VC/F (L)*(WT/70)' levodopa 45.4
    lq <- log(0.667); label("Levodopa inter-compartmental clearance Q/F at 70 kg (L/min)") # Table 2 'Q/F (L/min)*(WT/70)^0.75' 0.667
    lvp <- log(44.9); label("Levodopa peripheral volume VP/F at 70 kg (L)") # Table 2 'VP/F (L)*(WT/70)' 44.9

    # Levodopa absorption: two parallel transit chains (5 and 6 transit compartments, Results 'Base model')
    lmtt1 <- log(16.1); label("Levodopa fast-chain mean transit time MTT1 (min)") # Table 2 'MTT1 (min)' levodopa 16.1
    lmtt2 <- log(86.2); label("Levodopa slow-chain mean transit time MTT2 (min)") # Table 2 'MTT2 (min)' levodopa 86.2
    logitffo <- log(0.749 / (1 - 0.749)); label("Logit of the levodopa dose fraction entering the fast chain, fa1 = 0.749 (unitless)") # Table 2 'fa1' levodopa 0.749, estimated on the logit scale (footnote double-dagger)
    lfdepot <- fixed(log(1)); label("Levodopa relative bioavailability Frel (unitless)") # Table 2 'FREL' levodopa '1 FIX'
    lc0 <- log(0.120); label("Levodopa pre-dose concentration in patients (ug/mL)") # Table 2 'Pre-dose concentration (ug/mL)' levodopa 0.120; Results 'estimated to 0.12'

    # Levodopa covariate effects
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on levodopa CL/F and Q/F (unitless)") # Methods 'Covariate model': exponent held at 0.75 for CL/F and Q/F
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on levodopa VC/F and VP/F (unitless)") # Methods 'Covariate model': exponent held at 1 for VC/F and VP/F
    led50 <- log(86.7); label("Total carbidopa dose that halves levodopa CL/F, INTER (mg)") # Table 2 'INTERCL/F-CDAMT (mg)' 86.7; Results 'estimated to 86.7 mg'
    e_dis_pd_hy_cl <- fixed(-0.0616); label("Fractional change in levodopa CL/F per Hoehn and Yahr stage from 1.53 (unitless)") # Table 2 'CL/F-HY' levodopa '-0.0616 FIX' (AALASSO coefficient)
    e_age_fdepot <- fixed(0.00172); label("Fractional change in levodopa Frel per year of age from 48.7 years (1/year)") # Table 2 'FREL-AGE' levodopa '0.00172 FIX' (AALASSO coefficient)

    # Carbidopa disposition (Table 2, 'Carbidopa' column)
    lcl_carbidopa <- log(1.05); label("Carbidopa CL/F at 70 kg, age 48.7 years and Hoehn and Yahr 1.53 (L/min)") # Table 2 'CL/F (L/min)*(WT/70)^0.75' carbidopa 1.05
    lvc_carbidopa <- log(168); label("Carbidopa volume VC/F at 70 kg (L)") # Table 2 'VC/F (L)*(WT/70)' carbidopa 168

    # Carbidopa absorption: two parallel transit chains (3 and 10 transit compartments, Results 'Base model')
    lmtt1_carbidopa <- log(34.6); label("Carbidopa fast-chain mean transit time MTT1 (min)") # Table 2 'MTT1 (min)' carbidopa 34.6
    lmtt2_carbidopa <- log(121); label("Carbidopa slow-chain mean transit time MTT2 at the mean sex fraction 0.389 (min)") # Table 2 'MTT2 (min)' carbidopa 121
    logitffo_carbidopa <- log(0.570 / (1 - 0.570)); label("Logit of the carbidopa dose fraction entering the fast chain, fa1 = 0.570 (unitless)") # Table 2 'fa1' carbidopa 0.570, logit scale (footnote double-dagger)
    lfdepot_carbidopa <- fixed(log(1)); label("Carbidopa relative bioavailability Frel (unitless)") # Table 2 'FREL' carbidopa '1 FIX'
    lc0_carbidopa <- log(0.0256); label("Carbidopa pre-dose concentration in patients (ug/mL)") # Table 2 'Pre-dose concentration (ug/mL)' carbidopa 0.0256

    # Carbidopa covariate effects
    e_wt_cl_carbidopa <- fixed(0.75); label("Allometric exponent of body weight on carbidopa CL/F (unitless)") # Methods 'Covariate model': exponent held at 0.75 for CL/F
    e_wt_vc_carbidopa <- fixed(1); label("Allometric exponent of body weight on carbidopa VC/F (unitless)") # Methods 'Covariate model': exponent held at 1 for VC/F
    e_age_cl_carbidopa <- fixed(-0.0119); label("Fractional change in carbidopa CL/F per year of age from 48.7 years (1/year)") # Table 2 'CL/F-AGE' carbidopa '-0.0119 FIX' (AALASSO coefficient)
    e_dis_pd_hy_cl_carbidopa <- fixed(0.0141); label("Fractional change in carbidopa CL/F per Hoehn and Yahr stage from 1.53 (unitless)") # Table 2 'CL/F-HY' carbidopa '0.0141 FIX' (AALASSO coefficient)
    e_sexf_mtt2_carbidopa <- fixed(0.178); label("Fractional change in carbidopa MTT2 for female sex, centred at 0.389 (unitless)") # Table 2 'MTT2-SEX' carbidopa '0.178 FIX' (AALASSO coefficient)

    # IIV. Table 2 prints IIV as CV% = 100 * sqrt(omega^2): the footnote states the RSEs are
    # on the SD scale (PsN sumo convention), and the unfixed-coefficient column prints the
    # levodopa pre-dose IIV as 1.80 where the main column prints 180. Variance = (CV/100)^2.
    # fa1 IIV is on the logit scale (footnote double-dagger). No covariances reported.
    etalcl ~ 0.0225 # Table 2 'IIV CL/F' levodopa 15.0; 0.150^2
    etalvc ~ 0.156025 # Table 2 'IIV VC/F' levodopa 39.5; 0.395^2
    etalmtt1 ~ 0.119025 # Table 2 'IIV MTT1' levodopa 34.5; 0.345^2
    etalmtt2 ~ 0.023104 # Table 2 'IIV MTT2' levodopa 15.2; 0.152^2
    etalogitffo ~ 0.992016 # Table 2 'IIV fa1' levodopa 99.6 on the logit scale; 0.996^2
    etalfdepot ~ 0.080089 # Table 2 'IIV Frel' levodopa 28.3; 0.283^2
    etalc0 ~ 3.24 # Table 2 'IIV Pre-dose concentration' levodopa 180; 1.80^2
    etalcl_carbidopa ~ 0.043681 # Table 2 'IIV CL/F' carbidopa 20.9; 0.209^2
    etalmtt1_carbidopa ~ 0.101761 # Table 2 'IIV MTT1' carbidopa 31.9; 0.319^2
    etalmtt2_carbidopa ~ 0.073984 # Table 2 'IIV MTT2' carbidopa 27.2; 0.272^2
    etalogitffo_carbidopa ~ 0.659344 # Table 2 'IIV fa1' carbidopa 81.2 on the logit scale; 0.812^2
    etalfdepot_carbidopa ~ 0.231361 # Table 2 'IIV Frel' carbidopa 48.1; 0.481^2
    etalc0_carbidopa ~ 0.802816 # Table 2 'IIV Pre-dose concentration' carbidopa 89.6; 0.896^2

    # Residual error: combined proportional + additive, estimated separately per analyte and per population
    propSd_healthy <- 0.150; label("Levodopa proportional residual SD, healthy subjects (fraction)") # Table 2 'Proportional error (%)' healthy levodopa 15.0
    propSd_patient <- 0.0702; label("Levodopa proportional residual SD, patients (fraction)") # Table 2 'Proportional error (%)' patients levodopa 7.02
    addSd_healthy <- 0.00406; label("Levodopa additive residual SD, healthy subjects (ug/mL)") # Table 2 'Additive error (ug/mL)' healthy levodopa 0.00406
    addSd_patient <- fixed(0); label("Levodopa additive residual SD, patients (ug/mL)") # Table 2 'Additive error' patients levodopa '-'; Results: estimated to a low value and set to 0
    propSd_carbidopa_healthy <- 0.0770; label("Carbidopa proportional residual SD, healthy subjects (fraction)") # Table 2 'Proportional error (%)' healthy carbidopa 7.70
    propSd_carbidopa_patient <- 0.0386; label("Carbidopa proportional residual SD, patients (fraction)") # Table 2 'Proportional error (%)' patients carbidopa 3.86
    addSd_carbidopa_healthy <- 0.00415; label("Carbidopa additive residual SD, healthy subjects (ug/mL)") # Table 2 'Additive error (ug/mL)' healthy carbidopa 0.00415
    addSd_carbidopa_patient <- 0.00684; label("Carbidopa additive residual SD, patients (ug/mL)") # Table 2 'Additive error (ug/mL)' patients carbidopa 0.00684
  })

  model({
    # Covariate centring constants are the pooled Study 1 + 2 means (Table 1 footnote).
    # Linear covariate form: P = TVP * (1 + coeff * (COV - mean)) (Methods, second display equation).
    wtn <- WT / 70

    # Levodopa individual parameters. Carbidopa dose inhibits levodopa CL/F as
    # CL = TVCL / (1 + CDAMT / INTER) (Methods, third display equation), written here
    # with INTER as the half-maximal dose ed50.
    ed50 <- exp(led50)
    cl <- exp(lcl + etalcl) * wtn^e_wt_cl *
      (1 + e_dis_pd_hy_cl * (DIS_PD_HY - 1.53)) *
      ed50 / (ed50 + DOSE_CARBIDOPA_MG)
    vc <- exp(lvc + etalvc) * wtn^e_wt_vc
    q <- exp(lq) * wtn^e_wt_cl
    vp <- exp(lvp) * wtn^e_wt_vc
    mtt1 <- exp(lmtt1 + etalmtt1)
    mtt2 <- exp(lmtt2 + etalmtt2)
    logitffo_ind <- logitffo + etalogitffo
    ffo <- expit(logitffo_ind)
    fdepot <- exp(lfdepot + etalfdepot) * (1 + e_age_fdepot * (AGE - 48.7))
    c0 <- exp(lc0 + etalc0)

    # Carbidopa individual parameters
    cl_carbidopa <- exp(lcl_carbidopa + etalcl_carbidopa) * wtn^e_wt_cl_carbidopa *
      (1 + e_age_cl_carbidopa * (AGE - 48.7)) *
      (1 + e_dis_pd_hy_cl_carbidopa * (DIS_PD_HY - 1.53))
    vc_carbidopa <- exp(lvc_carbidopa) * wtn^e_wt_vc_carbidopa
    mtt1_carbidopa <- exp(lmtt1_carbidopa + etalmtt1_carbidopa)
    mtt2_carbidopa <- exp(lmtt2_carbidopa + etalmtt2_carbidopa) *
      (1 + e_sexf_mtt2_carbidopa * (SEXF - 0.389))
    logitffo_carbidopa_ind <- logitffo_carbidopa + etalogitffo_carbidopa
    ffo_carbidopa <- expit(logitffo_carbidopa_ind)
    fdepot_carbidopa <- exp(lfdepot_carbidopa + etalfdepot_carbidopa)
    c0_carbidopa <- exp(lc0_carbidopa + etalc0_carbidopa)

    # Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    kel_carbidopa <- cl_carbidopa / vc_carbidopa

    # Transit rate constants, Savic 2007 (the paper's reference 11): a chain of n
    # transit compartments between the dosing compartment and central has n + 1
    # first-order transfers, so ktr = (n + 1) / MTT.
    ktr1 <- (5 + 1) / mtt1
    ktr2 <- (6 + 1) / mtt2
    ktr1_carbidopa <- (3 + 1) / mtt1_carbidopa
    ktr2_carbidopa <- (10 + 1) / mtt2_carbidopa

    # Levodopa terminal slope lambda2 of the two-compartment system
    lambda2 <- 0.5 * ((kel + k12 + k21) - sqrt((kel + k12 + k21)^2 - 4 * kel * k21))

    # Levodopa ODEs. The dose is split fa1 / (1 - fa1) between depot (fast chain)
    # and depot2 (slow chain) by f(); dose the same amount into both.
    d/dt(depot) <- -ktr1 * depot
    d/dt(transit_fast1) <- ktr1 * (depot - transit_fast1)
    d/dt(transit_fast2) <- ktr1 * (transit_fast1 - transit_fast2)
    d/dt(transit_fast3) <- ktr1 * (transit_fast2 - transit_fast3)
    d/dt(transit_fast4) <- ktr1 * (transit_fast3 - transit_fast4)
    d/dt(transit_fast5) <- ktr1 * (transit_fast4 - transit_fast5)
    d/dt(depot2) <- -ktr2 * depot2
    d/dt(transit_slow1) <- ktr2 * (depot2 - transit_slow1)
    d/dt(transit_slow2) <- ktr2 * (transit_slow1 - transit_slow2)
    d/dt(transit_slow3) <- ktr2 * (transit_slow2 - transit_slow3)
    d/dt(transit_slow4) <- ktr2 * (transit_slow3 - transit_slow4)
    d/dt(transit_slow5) <- ktr2 * (transit_slow4 - transit_slow5)
    d/dt(transit_slow6) <- ktr2 * (transit_slow5 - transit_slow6)
    d/dt(central) <- ktr1 * transit_fast5 + ktr2 * transit_slow6 -
      kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Carbidopa ODEs, same dosing convention with depot_carbidopa / depot2_carbidopa
    d/dt(depot_carbidopa) <- -ktr1_carbidopa * depot_carbidopa
    d/dt(transit_fast1_carbidopa) <- ktr1_carbidopa * (depot_carbidopa - transit_fast1_carbidopa)
    d/dt(transit_fast2_carbidopa) <- ktr1_carbidopa * (transit_fast1_carbidopa - transit_fast2_carbidopa)
    d/dt(transit_fast3_carbidopa) <- ktr1_carbidopa * (transit_fast2_carbidopa - transit_fast3_carbidopa)
    d/dt(depot2_carbidopa) <- -ktr2_carbidopa * depot2_carbidopa
    d/dt(transit_slow1_carbidopa) <- ktr2_carbidopa * (depot2_carbidopa - transit_slow1_carbidopa)
    d/dt(transit_slow2_carbidopa) <- ktr2_carbidopa * (transit_slow1_carbidopa - transit_slow2_carbidopa)
    d/dt(transit_slow3_carbidopa) <- ktr2_carbidopa * (transit_slow2_carbidopa - transit_slow3_carbidopa)
    d/dt(transit_slow4_carbidopa) <- ktr2_carbidopa * (transit_slow3_carbidopa - transit_slow4_carbidopa)
    d/dt(transit_slow5_carbidopa) <- ktr2_carbidopa * (transit_slow4_carbidopa - transit_slow5_carbidopa)
    d/dt(transit_slow6_carbidopa) <- ktr2_carbidopa * (transit_slow5_carbidopa - transit_slow6_carbidopa)
    d/dt(transit_slow7_carbidopa) <- ktr2_carbidopa * (transit_slow6_carbidopa - transit_slow7_carbidopa)
    d/dt(transit_slow8_carbidopa) <- ktr2_carbidopa * (transit_slow7_carbidopa - transit_slow8_carbidopa)
    d/dt(transit_slow9_carbidopa) <- ktr2_carbidopa * (transit_slow8_carbidopa - transit_slow9_carbidopa)
    d/dt(transit_slow10_carbidopa) <- ktr2_carbidopa * (transit_slow9_carbidopa - transit_slow10_carbidopa)
    d/dt(central_carbidopa) <- ktr1_carbidopa * transit_fast3_carbidopa +
      ktr2_carbidopa * transit_slow10_carbidopa - kel_carbidopa * central_carbidopa

    f(depot) <- fdepot * ffo
    f(depot2) <- fdepot * (1 - ffo)
    f(depot_carbidopa) <- fdepot_carbidopa * ffo_carbidopa
    f(depot2_carbidopa) <- fdepot_carbidopa * (1 - ffo_carbidopa)

    # Pre-dose concentration (patients only) as residual drug from prior dosing, decaying
    # at the individual terminal slope (Methods 'Base model'). Placing the residual
    # amount on the slow eigenvector of the two-compartment system makes it decay
    # mono-exponentially at lambda2; for one-compartment carbidopa lambda = kel.
    central(0) <- (1 - DIS_HEALTHY) * c0 * vc
    peripheral1(0) <- (1 - DIS_HEALTHY) * c0 * vc * k12 / (k21 - lambda2)
    central_carbidopa(0) <- (1 - DIS_HEALTHY) * c0_carbidopa * vc_carbidopa

    # Observations: amounts in mg over volumes in L give mg/L = ug/mL
    Cc <- central / vc
    Cc_carbidopa <- central_carbidopa / vc_carbidopa

    # Population-specific combined error, Var = (prop * C)^2 + add^2
    propSd_ld <- propSd_healthy * DIS_HEALTHY + propSd_patient * (1 - DIS_HEALTHY)
    addSd_ld <- addSd_healthy * DIS_HEALTHY + addSd_patient * (1 - DIS_HEALTHY)
    sdCc <- sqrt((propSd_ld * Cc)^2 + addSd_ld^2)
    propSd_cd <- propSd_carbidopa_healthy * DIS_HEALTHY + propSd_carbidopa_patient * (1 - DIS_HEALTHY)
    addSd_cd <- addSd_carbidopa_healthy * DIS_HEALTHY + addSd_carbidopa_patient * (1 - DIS_HEALTHY)
    sdCc_carbidopa <- sqrt((propSd_cd * Cc_carbidopa)^2 + addSd_cd^2)

    Cc ~ add(sdCc)
    Cc_carbidopa ~ add(sdCc_carbidopa)
  })
}
