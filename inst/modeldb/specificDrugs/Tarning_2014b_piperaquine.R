# Population pharmacokinetic-pharmacodynamic model of oral piperaquine in
# patients with Plasmodium vivax malaria on the Thailand-Myanmar border
# (Tarning 2014, CPT Pharmacometrics Syst Pharmacol 3:e132;
# doi:10.1038/psp.2014.29).
#
# NOTE: this is NOT the same paper as Tarning_2014_piperaquine.R, which is
# Tarning 2014 Antimicrob Agents Chemother 58(4):2052-2058
# (doi:10.1128/AAC.02318-13), a food-effect study in P. falciparum malaria.

Tarning_2014b_piperaquine <- function() {
  description <- paste(
    "Population PK-PD model for oral piperaquine given as",
    "dihydroartemisinin-piperaquine in 116 patients (aged 3-57 years) with",
    "Plasmodium vivax malaria in Thailand (Tarning 2014 CPT Pharmacometrics",
    "Syst Pharmacol). Three transit-absorption compartments with the",
    "absorption rate constant equal to the transit rate constant",
    "(ktr = ka = 4 / MTT) feeding a three-compartment disposition model;",
    "fixed allometric scaling of all clearances (power 0.75) and volumes",
    "(power 1) on body weight centred at 50 kg; relative bioavailability",
    "fixed to 1 with IIV. Venous (Cc) and capillary (Ccap = 1.41 x Cc)",
    "plasma were fitted simultaneously with separate log-scale residual",
    "errors. Coupled interval-censored time-to-event model for recurrent",
    "P. vivax malaria: a constant baseline hazard (8.94 relapses/year)",
    "multiplied by a relapse surge function with peaks at 3, 6 and 9 weeks",
    "after the first dose (amplitude +123%) and by a sigmoid inhibitory",
    "effect of venous piperaquine (PC50 = 6.92 ng/mL, shape 9.28, maximal",
    "inhibition 1), exposed as the survivor function `sur`. Time 0 must be",
    "the first dose, because the relapse surges are timed from it.",
    sep = " "
  )
  reference <- paste(
    "Tarning J, Thana P, Phyo AP, Lwin KM, Hanpithakpong W, Ashley EA,",
    "Day NPJ, Nosten F, White NJ (2014). Population Pharmacokinetics and",
    "Antimalarial Pharmacodynamics of Piperaquine in Patients With",
    "Plasmodium vivax Malaria in Thailand. CPT: Pharmacometrics & Systems",
    "Pharmacology 3(8):e132. doi:10.1038/psp.2014.29.",
    sep = " "
  )
  vignette <- "Tarning_2014b_piperaquine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "piperaquine", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "piperaquine", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "piperaquine", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "piperaquine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "piperaquine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "piperaquine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "piperaquine", units = "mg", specimen = "plasma", verified = TRUE),
    cumhaz = list(
      analyte = "Cumulative hazard of recurrent P. vivax malaria",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight at admission",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed per subject at admission. Applied as a fixed allometric",
        "function on every clearance (exponent 0.75) and every volume",
        "(exponent 1) parameter, centred at 50 kg. The Supplementary",
        "Methods state that body weight was 'centered on the median weight",
        "of the population' (Table 1 PK cohort median 50 kg, range 9-67),",
        "and the supplementary NONMEM control streams (s03, s04) code",
        "TVCL = THETA(1) * (WT/50.0)**0.75 and TVV2 = THETA(2) * (WT/50.0).",
        "Body weight was the only covariate retained (Results, dOFV = 31.0).",
        sep = " "
      ),
      source_name = "WT"
    )
  )

  # Screened by Tarning 2014 and not retained. The Supplementary Methods list
  # the admission covariates tested by stepwise addition (p < 0.05) and
  # elimination (p < 0.01) on all population parameters: age, gender,
  # height, body weight, body temperature, days of fever before admission,
  # hematocrit, initial parasitaemia and gametocytaemia. Body temperature on
  # Vp2/F and age on Q1/F were selected forward but removed backward
  # (Results), and sex was selected forward on the baseline hazard but
  # removed backward (Results, 'Pharmacodynamics of piperaquine'). Days of
  # fever and gametocytaemia have no canonical column and are recorded in
  # population$notes instead.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at admission",
      units = "years",
      type = "continuous",
      notes = "Selected on Q1/F in the forward step (dOFV = 5.46) but not retained in the backward deletion. Table 1 PK cohort median 23.5 years (range 3-57)."
    ),
    BODYTEMP = list(
      description = "Body temperature at admission",
      units = "degC",
      type = "continuous",
      notes = "Selected on Vp2/F in the forward step (dOFV = 11.5) but not retained in the backward deletion. Table 1 PK cohort median 37 degC (range 36-40.4)."
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Selected on the baseline relapse hazard in the forward step (P < 0.05) but removed in the backward elimination, giving a covariate-free pharmacodynamic model. Table 1 PK cohort 90 male / 26 female."
    ),
    HT = list(
      description = "Height at admission",
      units = "cm",
      type = "continuous",
      notes = "Screened; not retained. Table 1 PK cohort median 160 cm (range 88-176)."
    ),
    HCT = list(
      description = "Hematocrit at admission",
      units = "%",
      type = "continuous",
      notes = "Screened; not retained. Table 1 PK cohort median 39% (range 22-50)."
    ),
    PARA = list(
      description = "Initial P. vivax parasitaemia",
      units = "parasites/uL",
      type = "continuous",
      notes = "Screened; not retained. Table 1 PK cohort median 7,185 /uL (range 80-58,714)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 116L,
    n_subjects_pkpd = 62L,
    n_subjects_efficacy = 241L,
    n_studies = 1L,
    age_range = "3-57 years; median 23.5 years (Table 1, PK cohort)",
    weight_range = "9-67 kg; median 50 kg (Table 1, PK cohort)",
    sex_female_pct = 22.4,
    disease_state = "Uncomplicated Plasmodium vivax mono-infection (Karen and Burmese patients)",
    dose_range = paste(
      "Oral dihydroartemisinin-piperaquine (Duo-cotecxin; 40 mg",
      "dihydroartemisinin + 320 mg piperaquine phosphate per tablet),",
      "2.2 mg/kg/day dihydroartemisinin and 17.8 mg/kg/day piperaquine",
      "phosphate once daily for 3 days with milk. Table 1 reports a total",
      "piperaquine BASE dose of 30.1 (27.2-36.1) mg/kg in the PK cohort.",
      sep = " "
    ),
    regions = "Thailand-Myanmar border (Shoklo Malaria Research Unit clinics, Mae Sot)",
    notes = paste(
      "PK sub-study nested in a 500-patient efficacy trial of DHA-PQ vs",
      "chloroquine; 250 patients in the DHA-PQ arm, 241 after exclusions.",
      "791 piperaquine plasma concentrations (435 venous, 356 capillary)",
      "from 116 patients: 62 sampled randomly over 69 days plus 54 who gave",
      "a single sample at recurrent malaria. Below-LLOQ data (< 0.5%) were",
      "omitted. The PK-PD time-to-event model was fitted to the 62 patients",
      "with both PK and PD data and evaluated against all 241. Covariates",
      "screened also included days of fever before admission (median 2,",
      "range 1-8) and gametocytaemia (65 yes / 51 no), neither retained.",
      "Eta shrinkage was up to 43.5% and epsilon shrinkage 15.0%.",
      sep = " "
    )
  )

  ini({
    # ---- Absorption -------------------------------------------------
    # 'A transit-compartment model with a fixed number of transit
    # compartments (n = 3, transit rate constant was set to be identical to
    # absorption rate constant)' (Results); Figure 1 caption:
    # KTR = (n + 1) / MTT. Control stream s03: KTR = 4/MTT.
    lmtt <- log(1.93); label("Mean absorption transit time MTT (h)") # Table 2: MTT = 1.93 h (RSE 9.16%)

    # ---- Disposition (typical 50-kg patient, venous matrix) ---------
    lcl <- log(54.1); label("Apparent elimination clearance CL/F (L/h)") # Table 2: CL/F = 54.1 L/h (RSE 5.56%)
    lvc <- log(2770); label("Apparent central volume of distribution Vc/F (L)") # Table 2: VC/F = 2,770 L (RSE 10.4%)
    lq <- log(663); label("Apparent intercompartmental clearance to peripheral1, Q1/F (L/h)") # Table 2: Q1/F = 663 L/h (RSE 19.7%)
    lvp <- log(10600); label("Apparent first peripheral volume of distribution Vp1/F (L)") # Table 2: VP1/F = 10,600 L (RSE 10.7%)
    lq2 <- log(72.7); label("Apparent intercompartmental clearance to peripheral2, Q2/F (L/h)") # Table 2: Q2/F = 72.7 L/h (RSE 12.1%)
    lvp2 <- log(25900); label("Apparent second peripheral volume of distribution Vp2/F (L)") # Table 2: VP2/F = 25,900 L (RSE 12.8%)

    # Relative bioavailability fixed to 1 for the population, carrying IIV
    # (Results: 'A relative bioavailability parameter (i.e., fixed to 100%
    # for the population) was implemented to allow quantification of the
    # interindividual variability in the absorption').
    lfdepot <- fixed(log(1)); label("Relative oral bioavailability F (unitless)") # Table 2: F = 1 fixed

    # ---- Allometric scaling -----------------------------------------
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL/F, Q1/F and Q2/F (unitless)") # Supplementary Methods and control stream s03: exponent 0.75 on clearances, fixed
    e_wt_vc <- fixed(1); label("Allometric exponent on Vc/F, Vp1/F and Vp2/F (unitless)") # Supplementary Methods and control stream s03: exponent 1.0 on volumes, fixed

    # ---- Capillary / venous matrix ----------------------------------
    # Control stream s03: CPC = CPV * CV, i.e. capillary = venous x TF, the
    # direction of the library's `cfcap` (Ccap = Cc * cfcap). Results: '41%
    # higher capillary concentrations compared with venous concentrations'.
    cfcap <- 1.41; label("Capillary:venous transformation factor for piperaquine plasma concentrations (unitless)") # Table 2: TF = 1.41 (RSE 2.11%)

    # ---- Pharmacodynamics: time to recurrent P. vivax malaria -------
    # Constant baseline hazard reported per year, converted to the model's
    # hourly time base: 8.94 / (365.25 * 24) = 1.0199e-3 /h (the control
    # stream s04 initial estimate is 0.001 /h).
    lhaz_base <- log(8.94 / (365.25 * 24)); label("Baseline hazard of recurrent P. vivax malaria (1/h)") # Table 2: theta_BHZ = 8.94 relapses/year (RSE 22.1%)

    # Multiple surge function (Supplementary Methods Eqs. 4-5; control
    # stream s04, SU):
    #   Relapse(t) = prod_k [1 + SA / (((t - k * PT)^2 / SW^2)^SHP + 1)],
    #   k = 1, 2, 3, PT = 504 h (3 weeks) fixed.
    lsurge_amp <- log(1.23); label("Peak fractional amplitude of each relapse surge in the hazard (unitless)") # Table 2: theta_SA = 123% (RSE 23.7%)
    lsurge_sw <- log(49.5); label("Half-width of each relapse surge at half its peak amplitude (h)") # Table 2: theta_SW = 49.5 h (RSE 54.3%)
    surge_shp <- fixed(4); label("Steepness (shape) factor of the relapse surge function (unitless)") # Table 2: theta_SHP = 4 fixed
    surge_period <- fixed(504); label("Interval between relapse surges, first surge at one interval after the first dose (h)") # Supplementary Methods: PT 'fixed peak time (i.e. 3 week intervals)'; control stream s04 peaks at 504, 1008, 1512 h

    # Inhibitory sigmoid Emax effect of the predicted VENOUS concentration
    # on the hazard, maximal inhibition 1 (Supplementary Methods Eq. 6;
    # control stream s04: 1 - (CP+DE)^SHP / ((CP+DE)^SHP + IC50^SHP)).
    lic50 <- log(6.92); label("Venous piperaquine concentration giving a 50% reduction of the relapse hazard, PC50 (ng/mL)") # Table 2: PC50 = 6.92 ng/mL (RSE 15.5%)
    lhill <- log(9.28); label("Shape factor of the inhibitory Emax function on the relapse hazard (unitless)") # Table 2: Gamma = 9.28 (RSE 25.6%)

    # ---- Inter-individual variability -------------------------------
    # Table 2 reports %CV = 100 * sqrt(exp(omega^2) - 1); the control
    # streams carry the variances, which reproduce every %CV exactly
    # (e.g. 100 * sqrt(exp(1.07) - 1) = 138%). IIV on CL/F, Vc/F, Q2/F and
    # TF was removed (Results, dOFV = 2.40).
    etalq ~ 1.07 # Table 2: IIV Q1/F 138% CV (RSE 32.2%); control stream s03 OMEGA(3) = 1.07
    etalvp ~ 0.125 # Table 2: IIV VP1/F 36.5% CV (RSE 40.2%); control stream s03 OMEGA(4) = 0.125
    etalvp2 ~ 0.186 # Table 2: IIV VP2/F 45.2% CV (RSE 71.0%); control stream s03 OMEGA(6) = 0.186
    etalmtt ~ 0.26 # Table 2: IIV MTT 54.5% CV (RSE 38.4%); control stream s03 OMEGA(7) = 0.26
    etalfdepot ~ 0.196 # Table 2: IIV F 46.5% CV (RSE 18.9%); control stream s03 OMEGA(8) = 0.196

    # ---- Residual error ---------------------------------------------
    # Concentrations were natural-log transformed and fitted with an
    # additive error on the log scale (control stream s03:
    # Y = LOG(pred) + W * EPS(1), SIGMA 1 FIX), so W is the log-scale SD.
    # Table 2: sigma_1 = capillary, sigma_2 = venous.
    expSd <- 0.321; label("Log-scale residual SD for venous plasma piperaquine (unitless)") # Table 2: sigma_2 = 0.321 (RSE 5.79%)
    expSd_Ccap <- 0.270; label("Log-scale residual SD for capillary plasma piperaquine (unitless)") # Table 2: sigma_1 = 0.270 (RSE 5.59%)

    # Placeholder residual on the survivor-function output; the source
    # likelihood is interval-censored time-to-event, not a continuous
    # endpoint. Same convention as Chotsiri_2019_piperaquine.R.
    addSd <- 0.001; label("Placeholder additive residual error on the relapse-free survival output `sur` (unitless); not from the source")
  })

  model({
    # ---- Individual parameters (control stream s03 $PK) -------------
    mtt <- exp(lmtt + etalmtt)
    cl <- exp(lcl) * (WT / 50)^e_wt_cl
    vc <- exp(lvc) * (WT / 50)^e_wt_vc
    q <- exp(lq + etalq) * (WT / 50)^e_wt_cl
    vp <- exp(lvp + etalvp) * (WT / 50)^e_wt_vc
    q2 <- exp(lq2) * (WT / 50)^e_wt_cl
    vp2 <- exp(lvp2 + etalvp2) * (WT / 50)^e_wt_vc

    # Three transit compartments plus the dosing compartment give four
    # first-order transfers to the central compartment: KTR = 4 / MTT.
    ktr <- 4 / mtt
    ka <- ktr

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # ---- ODE system -------------------------------------------------
    # Dose in mg of piperaquine BASE into `depot`.
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(central) <- ka * transit3 - kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    f(depot) <- exp(lfdepot + etalfdepot)

    # Venous concentration in ng/mL (mg / L * 1000; control stream
    # S2 = V2/1000) and capillary concentration CPC = CPV * TF.
    Cc <- 1000 * central / vc
    Ccap <- Cc * cfcap

    # ---- Time to recurrent P. vivax malaria (control stream s04) ----
    # Surges are timed from t = 0, which the dataset must place at the
    # first dose. The 1e-6 ng/mL offset is the source's guard (DE).
    surge_amp <- exp(lsurge_amp)
    surge_sw <- exp(lsurge_sw)
    surge <- (1 + surge_amp / (((t - surge_period)^2 / surge_sw^2)^surge_shp + 1)) *
      (1 + surge_amp / (((t - 2 * surge_period)^2 / surge_sw^2)^surge_shp + 1)) *
      (1 + surge_amp / (((t - 3 * surge_period)^2 / surge_sw^2)^surge_shp + 1))
    ic50 <- exp(lic50)
    hill <- exp(lhill)
    pqeff <- 1 - (Cc + 1e-6)^hill / ((Cc + 1e-6)^hill + ic50^hill)
    hazard <- exp(lhaz_base) * surge * pqeff

    d/dt(cumhaz) <- hazard
    cumhaz(0) <- 0
    sur <- exp(-cumhaz)

    Cc ~ lnorm(expSd)
    Ccap ~ lnorm(expSd_Ccap)
    sur ~ add(addSd)
  })
}
