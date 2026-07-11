Hien_2017_cipargamin <- function() {
  description <- paste(
    "Population PK / PD model of the spiroindolone antimalarial cipargamin",
    "(KAE609) in Vietnamese adults with acute uncomplicated Plasmodium",
    "falciparum malaria (Hien 2017). PK is a flexible transit-absorption",
    "chain (NN = 3 transit compartments fixed) into a one-compartment",
    "disposition model with fixed allometric body-weight scaling on CL/F and",
    "V/F. Bioavailability is anchored at F = 100% with between-subject",
    "variability. PD is a two-population parasite-clearance model: the",
    "asexual parasite pool at enrolment is split by a sensitive fraction",
    "Fsen and its complement (refractory subpopulation); sensitive parasites",
    "grow at a fixed Kgrow (= ln(10)/48 = 0.0479 /h corresponding to 10-fold",
    "multiplication per 48-h intraerythrocytic cycle) and are killed via an",
    "Emax / EC50 saturable-effect term driven by cipargamin plasma",
    "concentration; refractory parasites become active (join the sensitive",
    "pool) at first-order rate Kact and are then subject to the same drug",
    "kill. The typical Emax is dose-dependent per Hien 2017 Table 3",
    "equation 7: Emax_i = TVEmax * (dose/10)^COVdose_Emax. Structural PK",
    "estimates are the Table 3 population values referenced to a 59 kg",
    "typical patient receiving 10 mg cipargamin. Reproduces the individually",
    "predicted parasite clearance profiles in Hien 2017 Fig 3.",
    sep = " "
  )
  reference <- paste(
    "Hien TT, White NJ, Thuy-Nhien NT, Hoa NT, Thuan PD, Tarning J, Nosten F,",
    "Magnusson B, Jain JP, Hamed K.",
    "Estimation of the In Vivo MIC of Cipargamin in Uncomplicated Plasmodium",
    "falciparum Malaria.",
    "Antimicrob Agents Chemother. 2017;61(2):e01940-16.",
    "doi:10.1128/AAC.01940-16.",
    sep = " "
  )
  vignette <- "Hien_2017_cipargamin"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  paper_specific_compartments <- c(
    "parasite_sensitive",
    "parasite_refractory"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight at enrolment",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at enrolment. Fixed allometric scaling on CL/F (exponent",
        "0.75) and V/F (exponent 1.0) referenced at the population typical",
        "weight of 59 kg (Table 3 footnote a: 'Population mean parameter",
        "estimates were calculated for a typical patient with a body weight",
        "of 59 kg receiving a study drug dose of 10 mg')."
      ),
      source_name        = "WT"
    ),
    PARA = list(
      description        = "Baseline asexual Plasmodium falciparum parasitaemia at enrolment",
      units              = "parasites/uL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject; measured by microscopy and validated qPCR",
        "of the P. falciparum 18S rRNA gene (limit of detection 22",
        "parasites/mL = 0.022 parasites/uL). Enrolment inclusion required",
        "asexual parasite count 5,000-50,000 / uL (Methods 'Patients').",
        "Used inside model() to seed the initial-condition split into",
        "sensitive and refractory sub-populations:",
        "parasite_sensitive(0) = fsen * PARA and",
        "parasite_refractory(0) = (1 - fsen) * PARA. Both ODE states carry",
        "units of parasites/uL, matching the source paper's Fig 3",
        "individually predicted clearance-profile scale."
      ),
      source_name        = "PARA"
    ),
    DOSE_CIPARGAMIN_MG = list(
      description        = "Administered cipargamin dose",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-subject cipargamin dose (single oral administration on day 1;",
        "values 10, 15, 20, 21 (one patient administered 21 mg in error;",
        "analysed in the 20-mg cohort per Table 1 footnote a), or 30 mg).",
        "Enters the PD sub-model through the dose-dependent typical Emax",
        "equation Emax_i = TVEmax * (DOSE_CIPARGAMIN_MG / 10)^e_dose_emax",
        "(Hien 2017 Table 3 equation 7). Constant within a subject in the",
        "single-dose study design; must still be carried as a data column",
        "because the PD term evaluates at simulation time."
      ),
      source_name        = "DOSE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 25L,
    n_studies      = 1L,
    age_range      = "20-52 years",
    age_median     = "31 years (median across all cohorts; Table 1)",
    weight_range   = "42-78 kg",
    weight_median  = "59 kg (typical value used for parameter estimates; Table 3 footnote a)",
    sex_female_pct = 0,
    race_ethnicity = c(Asian = 100),
    disease_state  = paste(
      "Acute uncomplicated Plasmodium falciparum monoinfection confirmed",
      "by microscopy (asexual parasite count 5,000-50,000 / uL); axillary",
      "temperature >= 37.5 degC or oral / tympanic / rectal temperature",
      ">= 38 degC at screening or within the previous 24 h (Methods",
      "'Patients')."
    ),
    dose_range     = paste(
      "Single oral doses of cipargamin (formerly KAE609): 30 mg (n = 6),",
      "20 mg (n = 5; includes one patient in the 30-mg cohort who received",
      "21 mg in error per Table 1 footnote a), 15 mg (n = 7), or 10 mg",
      "(n = 7). Capsules supplied at doses of 1, 10, and 25 mg; combined",
      "under direct observation on day 1."
    ),
    regions        = "Vietnam (single-centre, Hospital for Tropical Diseases, Ho Chi Minh City)",
    trial_registration = "ClinicalTrials.gov NCT01836458",
    notes          = paste(
      "Adaptive single-dose de-escalation design (Methods 'Study design').",
      "Rescue treatment (dihydroartemisinin-piperaquine + primaquine) was",
      "given at day 42 or earlier upon early treatment failure or rising",
      "parasitaemia. PK-PD analysis pooled uncensored microscopy and qPCR",
      "parasite densities; DNA density estimates were censored when",
      "gametocyte-stage Pfs25 mRNA > 10% of the total measurement to",
      "reduce confounding by gametocytemia (Methods 'Pharmacodynamic",
      "modeling'). Supplement Text S1 (0.6 MB) was not on disk at",
      "extraction time; canonical allometric exponents 0.75 (CL) and 1.0",
      "(V) applied per Methods 'Pharmacokinetic modeling' ('Body weight",
      "was incorporated as a fixed allometric function on the clearance",
      "and volume parameters'). See vignette Assumptions and deviations."
    )
  )

  ini({
    # =========================================================================
    # PK: transit-absorption + one-compartment disposition (Hien 2017 Table 3)
    # -------------------------------------------------------------------------
    # All values are population typical values referenced to a 59 kg patient
    # (Table 3 footnote a). IIV %CV values are converted to log-normal
    # variance via omega^2 = log(1 + CV^2) per the paper's own reporting
    # formula (Table 3 footnote a): "The coefficient of variation (CV) for
    # interindividual variability was calculated as 100 * sqrt(exp(estimate)
    # - 1)."
    # =========================================================================
    lcl     <- log(1.72);   label("Apparent cipargamin clearance CL/F for a 59 kg patient (L/h)")             # Hien 2017 Table 3 row 'CL/F (liters/h) = 1.72 (95% CI 1.55-1.94; %RSE 5.69)'
    lvc     <- log(40.6);   label("Apparent cipargamin central volume V/F for a 59 kg patient (L)")           # Hien 2017 Table 3 row 'V/F (liters) = 40.6 (95% CI 37.3-44.7; %RSE 4.71)'
    lmtt    <- log(0.867);  label("Mean transit time MTT through the NN = 3 transit compartments (h)")        # Hien 2017 Table 3 row 'MTT (h) = 0.867 (95% CI 0.682-1.11; %RSE 12.6)'
    lka     <- log(1.65);   label("First-order absorption rate constant Ka from the last transit to central (1/h)")   # Hien 2017 Table 3 row 'Ka (1/h) = 1.65 (95% CI 1.04-2.81; %RSE 25.2)'
    lfdepot <- fixed(log(1));  label("Reference relative oral bioavailability F (fraction; fixed at 100%)")   # Hien 2017 Table 3 row 'F (%) = 100 fix' (paper fixes typical F to 100% and estimates only its IIV)
    nn_fix  <- fixed(3);    label("Number of Savic-style transit compartments (integer, unitless)")           # Hien 2017 Table 3 row 'No. trans comp = 3 fix' (paper fixed the number of transit compartments at 3 during estimation)

    # Fixed allometric exponents (Methods 'Pharmacokinetic modeling'
    # paragraph: 'Body weight was incorporated as a fixed allometric function
    # on the clearance and volume parameters.') Uses the canonical Anderson-
    # Holford 0.75 exponent on CL/F and 1.0 on V/F; the paper's supplement
    # (Text S1) with the explicit allometric equations was not on disk at
    # extraction time. See vignette Assumptions and deviations.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on cipargamin CL/F with body weight (unitless; fixed)") # Methods 'Pharmacokinetic modeling' + canonical Anderson-Holford 0.75
    e_wt_vc <- fixed(1.0);  label("Allometric exponent on cipargamin V/F with body weight (unitless; fixed)")  # Methods 'Pharmacokinetic modeling' + canonical Anderson-Holford 1.0

    # =========================================================================
    # PD: two-population parasite clearance model (Hien 2017 Table 3)
    # -------------------------------------------------------------------------
    # Kgrow is fixed to ln(10) / 48 h per Table 3 row 'K_grow (1/h) = 0.0479
    # fix' (parasite multiplication rate fixed to 10-fold multiplication per
    # 48-h intraerythrocytic cycle, references 14 and 15 in the paper).
    # Emax is dose-dependent through equation 7: Emax_i = TVEmax *
    # (dose/10)^COVdose_Emax with dose in mg and the reference dose 10 mg
    # (the smallest cohort). e_dose_emax is estimated on the linear scale
    # per Table 3.
    # =========================================================================
    lkgrow      <- fixed(log(0.0479));  label("Fixed asexual-parasite growth rate constant Kgrow (1/h); = ln(10)/48 for 10-fold multiplication per 48-h cycle") # Hien 2017 Table 3 row 'K_grow (1/h) = 0.0479 fix' (paper: 'K_grow was fixed to 10 per life cycle (i.e., 48 h)')
    lemax       <- log(0.564);          label("Typical maximum parasite-kill rate at DOSE = 10 mg, TVEmax (1/h)")                                                # Hien 2017 Table 3 row 'E_max (1/h) = 0.564 (95% CI 0.383-0.710; %RSE 12.9)'
    lec50       <- log(0.354);          label("Cipargamin concentration giving half-maximal parasite-kill rate EC50 (ng/mL)")                                    # Hien 2017 Table 3 row 'EC50 (ng/ml) = 0.354 (95% CI 0.222-0.466; %RSE 16.7)'
    lfsen       <- log(0.991);          label("Population fraction of asexual parasites that are fully drug-sensitive at enrolment Fsen (fraction)")             # Hien 2017 Table 3 row 'F_sen (%) = 99.1 (95% CI 98.8-99.4; %RSE 0.132)'
    lkact       <- log(0.0987);         label("First-order activation rate of refractory to sensitive parasites Kact (1/h)")                                     # Hien 2017 Table 3 row 'K_act (1/h) = 0.0987 (95% CI 0.0592-0.123; %RSE 14.7)'
    e_dose_emax <- 0.0463;              label("Power exponent of administered dose on typical Emax; Emax_i = TVEmax * (DOSE_CIPARGAMIN_MG/10)^e_dose_emax (unitless)") # Hien 2017 Table 3 row 'COVdose_E_max = 0.0463 (95% CI 0.0318-0.0604; %RSE 14.9)' with reference dose 10 mg per equation 7 in the supplemental text

    # =========================================================================
    # Between-subject variability (BSV / IIV). All values converted from the
    # Table 3 %CV values via omega^2 = log(1 + CV^2). Fsen omega uses the
    # log-normal convention per the paper's Table 3 footnote a; the eta on
    # lfsen can occasionally shift an individual fsen slightly above 1 (with
    # the corresponding refractory fraction slightly below zero), and the
    # vignette's Assumptions and deviations documents this limitation of
    # the log-normal-on-fraction encoding.
    # =========================================================================
    etalcl     ~ 0.03362  # Hien 2017 Table 3 %CV for IIV on CL/F = 18.5%; omega^2 = log(1 + 0.185^2) = 0.03362
    etalmtt    ~ 0.36011  # Hien 2017 Table 3 %CV for IIV on MTT = 65.2%; omega^2 = log(1 + 0.652^2) = 0.36011
    etalka     ~ 1.41146  # Hien 2017 Table 3 %CV for IIV on Ka = 176%; omega^2 = log(1 + 1.76^2) = 1.41146
    etalfdepot ~ 0.06617  # Hien 2017 Table 3 %CV for IIV on F (bioavailability) = 26.2%; omega^2 = log(1 + 0.262^2) = 0.06617

    etalemax   ~ 0.32418  # Hien 2017 Table 3 %CV for IIV on E_max = 62.2%; omega^2 = log(1 + 0.622^2) = 0.32418
    etalfsen   ~ 0.51240  # Hien 2017 Table 3 %CV for IIV on F_sen = 81.8%; omega^2 = log(1 + 0.818^2) = 0.51240
    etalkact   ~ 0.15672  # Hien 2017 Table 3 %CV for IIV on K_act = 41.5%; omega^2 = log(1 + 0.415^2) = 0.15672

    # =========================================================================
    # Residual error. Table 3 reports the PK residual as sigma = 17.5% CV,
    # implemented here as a proportional residual on the linear plasma
    # concentration scale. The PD residual is sigma = 109% CV, implemented as
    # a proportional residual on the linear total-parasitaemia (parasites/uL)
    # scale; this is a large residual variance consistent with a wide
    # dynamic range of parasite densities across observations. Neither the
    # main text nor Table 3 explicitly documents whether the NONMEM ERROR
    # block used log-additive residuals; the sigma-as-%CV reporting is
    # taken to mean linear-space proportional per the standard nlmixr2
    # interpretation. See vignette Assumptions and deviations.
    # =========================================================================
    propSd                  <- 0.175; label("Proportional residual error on plasma cipargamin concentration (fraction)")         # Hien 2017 Table 3 PK model row 'sigma (%CV) = 17.5 (95% CI 14.3-21.1; %RSE 10.3)'
    propSd_parasitemia_total <- 1.09; label("Proportional residual error on total parasitaemia (fraction; parasitemia_total)")   # Hien 2017 Table 3 PK-PD model row 'sigma (%CV) = 109 (95% CI 94.8-178; %RSE 17.0)'
  })

  model({
    # =========================================================================
    # 1. Individual PK parameters. Allometric weight scaling applies to CL/F
    #    and V/F only (Methods 'Pharmacokinetic modeling'); MTT, Ka, and F
    #    are body-weight-independent.
    # =========================================================================
    cl     <- exp(lcl + etalcl)         * (WT / 59)^e_wt_cl
    vc     <- exp(lvc)                  * (WT / 59)^e_wt_vc
    mtt    <- exp(lmtt + etalmtt)
    ka     <- exp(lka + etalka)
    fdepot <- exp(lfdepot + etalfdepot)

    # Savic transit-absorption-chain rate. NN = 3 transit compartments fixed
    # per Table 3; MTT is the mean transit time through those NN transits
    # (following the convention used in Bienczak_2016_nevirapine.R and
    # Ali_2018_amodiaquine.R). The final transit-to-central step is at
    # rate ka rather than ktr.
    ktr <- nn_fix / mtt

    # Elimination rate constant.
    kel <- cl / vc

    # =========================================================================
    # 2. Individual PD parameters. Emax is dose-dependent per equation 7;
    #    fsen may occasionally exceed 1 for extreme etas (documented as a
    #    limitation of the log-normal-on-fraction encoding).
    # =========================================================================
    emax_typical <- exp(lemax + etalemax)
    emax         <- emax_typical * (DOSE_CIPARGAMIN_MG / 10)^e_dose_emax
    ec50         <- exp(lec50)
    fsen         <- exp(lfsen + etalfsen)
    kact         <- exp(lkact + etalkact)
    kgrow        <- exp(lkgrow)

    # =========================================================================
    # 3. ODE system.
    #    3a. PK: depot -> transit1 -> transit2 -> transit3 -> central
    #        (three ktr transitions from depot through transit2, then ka
    #        into central; convention matches Bienczak_2016_nevirapine.R).
    #    3b. PD: two parasite subpopulations with drug-kill on the sensitive
    #        pool and activation of the refractory pool into the sensitive
    #        pool.
    # =========================================================================

    # Plasma cipargamin concentration (ng/mL): amt (mg) / vc (L) = mg/L =
    # 1000 * ng/mL, i.e., divide by vc and multiply by 1e6 / 1e3 = 1000.
    Cc <- 1000 * central / vc

    # Drug-driven kill rate (h^-1); Emax / EC50 saturable function.
    kkill <- emax * Cc / (ec50 + Cc)

    d/dt(depot)               <- -ktr * depot
    d/dt(transit1)            <-  ktr * depot    - ktr * transit1
    d/dt(transit2)            <-  ktr * transit1 - ktr * transit2
    d/dt(transit3)            <-  ktr * transit2 - ka  * transit3
    d/dt(central)             <-  ka  * transit3 - kel * central
    d/dt(parasite_sensitive)  <-  kgrow * parasite_sensitive -
                                  kkill * parasite_sensitive +
                                  kact  * parasite_refractory
    d/dt(parasite_refractory) <- -kact  * parasite_refractory

    # Initial conditions for the parasite pools: total baseline parasitaemia
    # PARA (parasites/uL) is measured at enrolment and split into the
    # sensitive fraction (fsen * PARA) and its refractory complement
    # ((1 - fsen) * PARA). Both ODE states carry units of parasites/uL.
    parasite_sensitive(0)  <- fsen       * PARA
    parasite_refractory(0) <- (1 - fsen) * PARA

    # =========================================================================
    # 4. Bioavailability applied to the depot (dosing) compartment. Reference
    #    F is anchored at 1 (Table 3 row 'F (%) = 100 fix'); IIV enters via
    #    etalfdepot.
    # =========================================================================
    f(depot) <- fdepot

    # =========================================================================
    # 5. Observation and residual error.
    #    Total observable parasitaemia is the sum of the sensitive and
    #    refractory subpopulations (Fig 2 legend: 'the observed total
    #    parasitemia is the sum of sensitive and refractory parasites').
    # =========================================================================
    parasitemia_total <- parasite_sensitive + parasite_refractory

    Cc                ~ prop(propSd)
    parasitemia_total ~ prop(propSd_parasitemia_total)
  })
}
