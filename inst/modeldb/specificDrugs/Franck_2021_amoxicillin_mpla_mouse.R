Franck_2021_amoxicillin_mpla_mouse <- function() {
  description <- paste(
    "Preclinical (mouse). Sequential PK/PD disease-treatment and survival",
    "(time-to-event) model for oral amoxicillin with or without intraperitoneal",
    "monophosphoryl lipid A (MPLA) in mice with Streptococcus pneumoniae",
    "serotype 1 pneumonia (Franck 2021). The PK layer is the two-compartment",
    "submodel with the MPLA x dose clearance interaction, fixed at its typical",
    "values. Separate first-order effect compartments link serum amoxicillin to",
    "lung and spleen. Lung bacteria grow with a delayed onset",
    "(kg * (1 - exp(-klag * t))), are removed by treatment-unrelated killing",
    "and natural death (kkill, multiplied by 1.40 under MPLA) and by a steep",
    "sigmoidal amoxicillin Emax kill (Hill 20). Bacteria reach the spleen by a",
    "Savic gamma-kernel transit (n = 23, MTT = 40.8 h) scaled by the current",
    "lung burden, and are killed there by a power amoxicillin effect",
    "(kAMX * Ce^5.06) and a first-order MPLA kill. Survival over 14 days",
    "follows a surge-function hazard reduced exponentially by the model-predicted",
    "time of serum amoxicillin above the MIC (h) and by MPLA coadministration.",
    "Model time zero is the time of infection; treatment is given at 12 h.",
    sep = " "
  )
  reference <- paste(
    "Franck S, Michelet R, Casilag F, Sirard JC, Wicha SG, Kloft C.",
    "A Model-Based Pharmacokinetic/Pharmacodynamic Analysis of the Combination",
    "of Amoxicillin and Monophosphoryl Lipid A Against S. pneumoniae in Mice.",
    "Pharmaceutics. 2021;13(4):469. doi:10.3390/pharmaceutics13040469.",
    "PK/PD parameters from Table 1, survival parameters from Table 2, equations",
    "from Supplementary Section S2 (Eqs. 1-7) and Figure 1.",
    sep = " "
  )
  vignette <- "Franck_2021_amoxicillin_mpla_pneumonia"
  units <- list(time = "h", dosing = "ug", concentration = "ug/mL")

  covariateData <- list(
    DOSE_AMOXICILLIN_UG = list(
      description = "Administered single oral amoxicillin dose in ug (absolute amount per mouse)",
      units = "ug",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters only through the MPLA pharmacokinetic interaction on clearance",
        "(Supplementary Eq. 1), so it has no effect when CONMED_MPLA = 0. The",
        "paper's arithmetic reproduces exactly with the mg/kg dose multiplied",
        "by a 25 g body weight (1.2 mg/kg = 30 ug).",
        sep = " "
      ),
      source_name = "DOSE"
    ),
    CONMED_MPLA = list(
      description = "Monophosphoryl lipid A on board (1 = at or after the 2.0 mg/kg IP MPLA dose, 0 = before it or no MPLA)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Time-varying: set 0 before the MPLA dose (given with the amoxicillin",
        "dose 12 h after infection) and 1 from that time on for MPLA-treated",
        "mice. Applying it from the time of infection instead gives a 2.28",
        "rather than the published 1.71 log10 CFU/lung reduction at 36 h, so the",
        "published model switches the effect on at treatment. Acts on the lung",
        "kkill, the spleen MPLA kill, the amoxicillin clearance and the survival",
        "hazard.",
        sep = " "
      ),
      source_name = "MPLA"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "amoxicillin", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "amoxicillin", units = "ug", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "amoxicillin", units = "ug", specimen = "serum", verified = TRUE),
    effect_lung = list(analyte = "amoxicillin", units = "ug/mL", specimen = "not applicable", verified = TRUE),
    effect_spleen = list(analyte = "amoxicillin", units = "ug/mL", specimen = "not applicable", verified = TRUE),
    bacteria_lung = list(analyte = "Streptococcus pneumoniae", units = "CFU", specimen = "tissue", verified = TRUE),
    bacteria_spleen = list(analyte = "Streptococcus pneumoniae", units = "CFU", specimen = "tissue", verified = TRUE),
    t_above_mic = list(analyte = "amoxicillin", units = "h", specimen = "not applicable", verified = TRUE),
    cumhaz = list(analyte = "death", units = NA_character_, specimen = "not applicable", verified = TRUE)
  )

  # Per-paper states and outputs: two latent organ effect compartments, the
  # organ-specific bacterial burdens and their log10 observations, and the
  # time-above-MIC accumulator (precedent: Lallemand_2023_benzylpenicillin_horse.R,
  # Assmus_2025_benznidazole_mouse.R).
  paper_specific_compartments <- c(
    "effect_lung",
    "effect_spleen",
    "bacteria_lung",
    "bacteria_spleen",
    "log_cfu_lung",
    "log_cfu_spleen",
    "t_above_mic"
  )

  population <- list(
    species = "mouse (RjOrl:Swiss / CD-1 and Balb/cJRj, female, S. pneumoniae serotype 1 pneumonia model)",
    n_subjects = 936,
    n_studies = 3,
    age_range = "6-8 weeks",
    weight_range = "~25 g",
    sex_female_pct = 100,
    disease_state = paste(
      "Pneumonia after intranasal infection with 1-4 x 10^6 CFU Streptococcus",
      "pneumoniae serotype 1 (clinical isolate E1586, amoxicillin MIC 0.016 mg/L);",
      "treatment given 12 h after infection.",
      sep = " "
    ),
    dose_range = paste(
      "Amoxicillin single oral gavage 0.4 or 14 mg/kg (PK study) and 0.2, 0.4 or",
      "1.2 mg/kg (PD and survival studies), with or without monophosphoryl",
      "lipid A 2.0 mg/kg intraperitoneally; untreated and MPLA-alone groups.",
      sep = " "
    ),
    regions = "France (Institut Pasteur de Lille)",
    notes = paste(
      "Supplementary Section S1: 106 RjOrl:Swiss mice in the PK study (serum",
      "amoxicillin), 634 RjOrl:Swiss and Balb/cJRj mice in the PD study",
      "(lung and spleen CFU at -12 to 36 h relative to treatment, one organ",
      "harvest per mouse) and 196 mice in the survival study (monitored every",
      "24 h for 14 days after infection). Pooled studies; mouse type was",
      "tested and not retained. No IIV could be estimated for the PD or",
      "survival data (one observation per mouse).",
      sep = " "
    )
  )

  ini({
    # --- Pharmacokinetics: fixed at the PK submodel typical values (Table 1,
    # footnote *: 'Fixed parameter estimates of developed PK submodel (Table S1)').
    lka <- fixed(log(5.04))
    label("First-order absorption rate constant ka (1/h)") # Table 1: ka = 5.04 1/h *
    ltlag <- fixed(log(0.125))
    label("Absorption lag time tlag (h)") # Table 1: tlag = 0.125 h *
    lvc <- fixed(log(15.4))
    label("Apparent central volume of distribution Vc/F (mL)") # Table 1: Vc/F = 15.4 mL *
    lvp <- fixed(log(50.7))
    label("Apparent peripheral volume of distribution Vp/F (mL)") # Table 1: Vp/F = 50.7 mL *
    lq <- fixed(log(71.9))
    label("Apparent intercompartmental clearance Q/F (mL/h)") # Table 1: Q/F = 71.9 mL/h *
    lcl <- fixed(log(124))
    label("Apparent amoxicillin clearance without MPLA, CL/F (mL/h)") # Table 1: CL_AMX/F = 124 mL/h *
    lfdepot <- fixed(log(1))
    label("Oral bioavailability F (fraction)") # Table 1 abbreviations: 'F: Bioavailability of AMX fixed to 1'
    e_dose_mpla_cl <- fixed(-0.145)
    label("Additive change in CL/F per ug amoxicillin dose when MPLA is coadministered (mL/h/ug)") # Table 1: FC_AMX+MPLA = -0.145 mL/h/ug *; Supplementary Eq. 1

    # --- Effect compartments (Supplementary Eq. 4)
    lke0_lung <- log(0.125)
    label("Serum-to-lung effect-compartment rate constant ke0,lung (1/h)") # Table 1: ke0,lung = 0.125 1/h (RSE 19.7%)
    lke0_spleen <- log(0.0435)
    label("Serum-to-spleen effect-compartment rate constant ke0,spleen (1/h)") # Table 1: ke0,spleen = 0.0435 1/h (RSE 17.7%)

    # --- Bacterial disease submodel (Supplementary Eqs. 2-3)
    bl_log_cfu_lung <- fixed(6.12)
    label("Initial lung bacterial burden at infection (log10 CFU/lung)") # Table 1: N_bacteria,t=0 = 6.12 log10(CFU/lung) ** (fixed; 'at -12 h', i.e. at infection)
    lkg <- log(0.477)
    label("First-order bacterial growth rate constant in lung kg (1/h)") # Table 1: kg = 0.477 1/h (RSE 7.00%)
    lklag <- log(0.0595)
    label("First-order rate constant for the delayed onset of lung growth klag (1/h)") # Table 1: klag = 0.0595 1/h (RSE 46.2%)
    lkkill_lung <- log(0.274)
    label("Treatment-unrelated killing and natural death rate constant in lung kkill,lung (1/h)") # Table 1: kkill,lung = 0.274 1/h (RSE 20.3%)
    lntr <- log(23.0)
    label("Number of lung-to-spleen transit compartments n (unitless)") # Table 1: n = 23.0 (RSE 13.1%)
    lmtt <- log(40.8)
    label("Mean lung-to-spleen transit time MTT (h)") # Table 1: MTT = 40.8 h (RSE 4.30%); Supplementary Section S3 text says 42.0 h

    # --- Disease and treatment submodel
    e_conmed_mpla_kkill_lung <- 1.40
    label("Multiplicative factor on kkill,lung under MPLA (unitless ratio)") # Table 1: MPLA_lung = 1.40 (RSE 6.10%); 'fractional change of kkill,lung in presence of MPLA', Supplementary Section S2 'proportionality factor'
    lemax <- log(0.255)
    label("Maximum amoxicillin kill rate in lung Emax (1/h)") # Table 1: Emax = 0.255 1/h (RSE 6.00%)
    lec50 <- log(0.00109)
    label("Lung effect-compartment amoxicillin concentration for half-maximal kill EC50 (ug/mL)") # Table 1: EC50 = 0.00109 ug/mL (RSE 29.4%; bootstrap CI 0.000134-0.00146); Results text '0.0109' is a typo
    hill_lung <- fixed(20)
    label("Hill factor of the lung amoxicillin Emax kill (unitless)") # Table 1: H_lung = 20 ** (fixed after log-likelihood profiling, 95% CI 1.96-105)
    lkmpla_spleen <- log(3.71)
    label("First-order MPLA kill rate constant in spleen kMPLA,spleen (1/h)") # Table 1: kMPLA,spleen = 3.71 1/h (RSE 27.5%)
    lkamx_spleen <- log(10^13.7)
    label("Amoxicillin power-model kill coefficient in spleen kAMX (1/h per (ug/mL)^Hspleen)") # Table 1: kAMX = 13.7 [log10(h-1)], i.e. log10(kAMX) = 13.7
    hill_spleen <- 5.06
    label("Exponent of the spleen amoxicillin power kill model H_spleen (unitless)") # Table 1: H_spleen = 5.06 (RSE 23.9%)

    # --- Survival (time-to-event) model, time since infection in h (Table 2; Supplementary Eqs. 6-7)
    lsa_haz <- log(0.0404)
    label("Surge amplitude of the baseline hazard SA (1/h)") # Table 2: SA = 0.0404 1/h (RSE 20.9%)
    lsw_haz <- log(35.7)
    label("Surge width at half-maximum intensity SW (h)") # Table 2: SW = 35.7 h (RSE 17.3%)
    lgam_haz <- log(2.24)
    label("Shape parameter of the surge peak gamma (unitless)") # Table 2: gamma = 2.24 (RSE 26.0%)
    lpt_haz <- log(89.2)
    label("Peak time of the baseline hazard PT after infection (h)") # Table 2: PT = 89.2 h (RSE 4.80%)
    e_tmic_haz <- -0.926
    label("Log hazard ratio per hour of serum amoxicillin above the MIC (1/h)") # Table 2: beta T>MIC = -0.926 (RSE 16.7%)
    e_conmed_mpla_haz <- -1.32
    label("Log hazard ratio for MPLA coadministration (unitless)") # Table 2: beta MPLA_TTE = -1.32 (RSE 16.6%)
    mic <- fixed(0.016)
    label("Amoxicillin MIC of S. pneumoniae serotype 1 isolate E1586 (ug/mL)") # Section 2.1 / Supplementary Section S1: MIC_AMX = 0.016 mg/L

    # --- Residual error (additive on the log10 scale; Table 1 footnote ***, SD scale)
    addSd_log_cfu_lung <- 1.12
    label("Additive residual error on log10 lung CFU (log10 CFU/lung)") # Table 1: RUV lung = 1.12 (RSE 3.90%)
    addSd_log_cfu_spleen <- 1.81
    label("Additive residual error on log10 spleen CFU (log10 CFU/spleen)") # Table 1: RUV spleen = 1.81 (RSE 4.40%)
    addSd_sur <- fixed(0.001)
    label("Placeholder additive residual error on the survival-probability output (unitless); not from the source") # not in the source -- placeholder so the forward-simulation output has an error model
  })
  model({
    # 1. Pharmacokinetics (typical values; no IIV in the PK/PD model).
    ka <- exp(lka)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q <- exp(lq)
    cl <- exp(lcl) + e_dose_mpla_cl * DOSE_AMOXICILLIN_UG * CONMED_MPLA
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    f(depot) <- exp(lfdepot)
    alag(depot) <- exp(ltlag)

    Cc <- central / vc

    # 2. Effect compartments in lung and spleen (Supplementary Eq. 4).
    ke0_lung <- exp(lke0_lung)
    ke0_spleen <- exp(lke0_spleen)
    d/dt(effect_lung) <- ke0_lung * (Cc - effect_lung)
    d/dt(effect_spleen) <- ke0_spleen * (Cc - effect_spleen)

    # 3. Lung bacteria (Supplementary Eq. 2 plus the treatment terms of
    #    Figure 1). `t` is the time since infection.
    kg <- exp(lkg)
    klag <- exp(lklag)
    kkill_lung <- exp(lkkill_lung) * e_conmed_mpla_kkill_lung^CONMED_MPLA
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    kill_amx_lung <- emax * effect_lung^hill_lung / (ec50^hill_lung + effect_lung^hill_lung)
    d/dt(bacteria_lung) <- kg * (1 - exp(-klag * t)) * bacteria_lung -
      kkill_lung * bacteria_lung - kill_amx_lung * bacteria_lung
    bacteria_lung(0) <- 10^bl_log_cfu_lung

    # 4. Lung-to-spleen transit: Savic analytical transit input (Supplementary
    #    Eq. 3, ktr = (n + 1) / MTT, n! by the Stirling approximation) with the
    #    current lung burden in place of a dose. The transit does not deplete
    #    the lung, and the spleen has no outflow other than drug killing.
    ntr <- exp(lntr)
    mtt <- exp(lmtt)
    ktr <- (ntr + 1) / mtt
    lnfac_ntr <- log(sqrt(2 * pi)) + (ntr + 0.5) * log(ntr) - ntr
    transit_in_spleen <- bacteria_lung * ktr * (ktr * t)^ntr * exp(-ktr * t - lnfac_ntr)

    kmpla_spleen <- exp(lkmpla_spleen)
    kamx_spleen <- exp(lkamx_spleen)
    kill_amx_spleen <- kamx_spleen * effect_spleen^hill_spleen
    d/dt(bacteria_spleen) <- transit_in_spleen - kill_amx_spleen * bacteria_spleen -
      kmpla_spleen * CONMED_MPLA * bacteria_spleen

    # 5. Survival. Time of serum amoxicillin above the MIC (h), then a surge
    #    baseline hazard (Supplementary Eq. 6) whose covariate effects are
    #    time-constant in the source (Eq. 7), so S(t) = exp(-H0(t) * HR).
    d/dt(t_above_mic) <- (Cc >= mic)
    sa_haz <- exp(lsa_haz)
    sw_haz <- exp(lsw_haz)
    gam_haz <- exp(lgam_haz)
    pt_haz <- exp(lpt_haz)
    hazard0 <- sa_haz / (((t - pt_haz)^2 / sw_haz^2)^gam_haz + 1)
    d/dt(cumhaz) <- hazard0
    hr_haz <- exp(e_tmic_haz * t_above_mic + e_conmed_mpla_haz * CONMED_MPLA)
    hazard <- hazard0 * hr_haz
    sur <- exp(-cumhaz * hr_haz)

    # 6. Observations (log10 CFU per organ; additive error on the log10 scale).
    log_cfu_lung <- log10(bacteria_lung)
    log_cfu_spleen <- log10(bacteria_spleen)
    log_cfu_lung ~ add(addSd_log_cfu_lung)
    log_cfu_spleen ~ add(addSd_log_cfu_spleen)
    sur ~ add(addSd_sur)
  })
}
