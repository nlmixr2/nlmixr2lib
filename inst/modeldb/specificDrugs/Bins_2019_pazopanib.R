Bins_2019_pazopanib <- function() {
  description <- paste(
    "Two-compartment population PK model for oral pazopanib in adult cancer patients",
    "(Bins 2019; the Yu 2017 base model extended with a CYP3A4*22 effect on clearance).",
    "Each dose is split over two parallel first-order absorption routes: a fast route with",
    "no lag time and a slow route that starts after a lag time. Relative bioavailability",
    "falls with the dose in an Emax manner (reference 200 mg) and decays exponentially",
    "over the treatment course to an asymptote. Clearance is 36% lower in CYP3A4*22",
    "carriers and is scaled by the relative amount of a metabolic enzyme pool induced",
    "by co-administered ifosfamide (Kerbusch 2000 autoinduction model, fixed). IIV on",
    "the fast absorption rate, clearance, peripheral volume and relative bioavailability;",
    "inter-occasion variability on relative bioavailability; combined additive and",
    "proportional residual error."
  )
  reference <- paste(
    "Bins S, Huitema ADR, Laven P, el Bouazzaoui S, Yu H, van Erp N, van Herpen C,",
    "Hamberg P, Gelderblom H, Steeghs N, Sleijfer S, van Schaik RHN, Mathijssen RHJ,",
    "Koolen SLW. Impact of CYP3A4*22 on Pazopanib Pharmacokinetics in Cancer Patients.",
    "Clin Pharmacokinet. 2019;58(5):651-658. doi:10.1007/s40262-018-0719-5.",
    "Base model: Yu H, van Erp N, Bins S, et al. Development of a Pharmacokinetic Model",
    "to Describe the Complex Pharmacokinetics of Pazopanib in Cancer Patients.",
    "Clin Pharmacokinet. 2017;56(3):293-303. doi:10.1007/s40262-016-0443-y.",
    "Ifosfamide sub-model: Kerbusch T, et al. Br J Clin Pharmacol. 2000;49(6):555-561.",
    "doi:10.1046/j.1365-2125.2000.00217.x."
  )
  vignette <- "Bins_2019_pazopanib"
  dosing <- c("depot", "depot2", "central_ifosfamide")
  units <- list(
    time = "h",
    dosing = "mg (pazopanib, depot and depot2); umol (ifosfamide, central_ifosfamide)",
    concentration = "mg/L (pazopanib Cc); umol/L (ifosfamide Cc_ifosfamide)"
  )

  compartmentData <- list(
    depot = list(analyte = "pazopanib", units = "mg", specimen = "administration site", verified = TRUE),
    depot2 = list(analyte = "pazopanib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "pazopanib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "pazopanib", units = "mg", specimen = "tissue", verified = TRUE),
    central_ifosfamide = list(analyte = "ifosfamide", units = "umol", specimen = "plasma", verified = TRUE),
    enz_pool = list(
      analyte = "ifosfamide-inducible metabolic enzyme (relative amount)",
      units = "(fraction of baseline)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    SNP_CYP3A4_RS35599367 = list(
      description = "CYP3A4*22 (rs35599367, 15389C>T) carrier indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = CYP3A4 *1/*1 wild type",
      notes = "1 = carrier (all 11 carriers in Bins 2019 were heterozygous *1/*22; Table 1). Control stream sets FL1 = 1 when CYP3 is 1 or 2, so a homozygote takes the same effect.",
      source_name = "CYP3"
    ),
    DOSE_PAZOPANIB_MG = list(
      description = "Pazopanib amount of the current administration (mg)",
      units = "mg",
      type = "continuous",
      reference_category = "200 mg (relative bioavailability = 1)",
      notes = "Drives the dose-dependent relative bioavailability 1 - Emax * (DOSE - 200) / (ED50 + DOSE - 200). Per administration: 400 for 400 mg twice daily. Doses in the source data ranged 200-1000 mg once daily (Yu 2017 Table 1 also lists 1200 mg).",
      source_name = "DOSE"
    ),
    T_FIRSTDOSE = list(
      description = "Time since the first pazopanib dose of the treatment course (h)",
      units = "h",
      type = "continuous",
      reference_category = "0 at the first dose",
      notes = "The control stream uses the NONMEM TIME data item in TVF1 = 1 - DCRP + DCRP * EXP(-LAMBDA * TIME); the decay rate is printed per day and divided by 24 in $PK, so TIME is in hours since the start of treatment. When a simulation starts at the first dose, set T_FIRSTDOSE equal to the model time.",
      source_name = "TIME"
    ),
    OCC = list(
      description = "Pharmacokinetic sampling occasion (1, 2 or 3)",
      units = "(count)",
      type = "categorical",
      reference_category = "n/a",
      notes = "Control stream PKDAY: IOV1 = ETA(9), ETA(10), ETA(11) for PKDAY 1, 2, 3. Up to three sampling occasions per patient (Yu 2017 Table 1). Records with OCC outside 1-3 receive no IOV.",
      source_name = "PKDAY"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 94L,
    n_studies = 3L,
    n_observations = 761L,
    age_range = "median 57 years (IQR 48-64)",
    sex_female_pct = 29,
    race_ethnicity = c(Caucasian = 97, Asian = 2, Black = 1),
    disease_state = "Adults with advanced solid tumours (tumour types not recorded) treated with pazopanib alone or combined with docetaxel or ifosfamide.",
    dose_range = "Pazopanib 200-1000 mg once daily as starting dose (200 mg n = 18, 400 mg n = 38, 600 mg n = 3, 800 mg n = 20, 1000 mg n = 15).",
    regions = "Netherlands",
    genotypes = "CYP3A4*22: 83 *1/*1, 11 *1/*22, 0 *22/*22 (minor allele frequency 6%).",
    notes = "Bins 2019 Table 2: 13 patients from a pazopanib TDM study (NTR3293), 40 from a phase I pazopanib + docetaxel study (NTR2309) and 41 from a phase I pazopanib + ifosfamide study (NTR2063); 94 of the 96 patients of the Yu 2017 model were genotyped."
  )

  ini({
    # Structural parameters: Bins 2019 Table 4 'Final model estimates'.
    lka_fast <- log(0.39); label("Fast absorption rate constant KaF (1/h)") # Table 4 'Fast absorption rate constant' 0.39 h-1 (RSE 23%); control stream THETA(3) KA1
    lka_slow <- log(0.12); label("Slow absorption rate constant KaS (1/h)") # Table 4 'Slow absorption rate constant' 0.12 h-1 (RSE 28%); control stream THETA(4) KA2
    logitfrel <- -0.575364; label("Logit of the fraction of the bioavailable dose absorbed by the fast route (unitless)") # Table 4 'Fraction of fast absorption' 36% (RSE 36%); logit(0.36) = log(0.36 / 0.64); control stream TVFR = THETA(10) is the odds, F1 = TVFR / (1 + TVFR) * TVF
    ltlag <- log(0.97); label("Lag time of the slow absorption route (h)") # Table 4 'Lag time between fast and slow absorption' 0.97 h (RSE 6%); control stream THETA(5) ALAG3 on the slow-route dose compartment
    lcl <- log(0.28); label("Clearance in CYP3A4 wild type, before ifosfamide induction (L/h)") # Table 4 'Clearance' 0.28 L/h (RSE 23%); control stream THETA(1)
    lvc <- log(2.4); label("Central volume of distribution (L)") # Table 4 'Volume of distribution of central compartment' 2.4 L (RSE 34%); control stream THETA(2)
    lq <- log(0.97); label("Inter-compartmental clearance (L/h)") # Table 4 'Inter-compartment clearance' 0.97 L/h (RSE 29%); control stream THETA(6)
    lvp <- log(24.3); label("Peripheral volume of distribution (L)") # Table 4 'Volume of distribution of peripheral compartment' 24.3 L (RSE 28%); control stream THETA(7)

    # CYP3A4*22: Ppop = theta_x * theta_y^pg (Bins 2019 Methods equation).
    # Table 4 prints the carrier clearance (0.18 L/h, RSE 18%) rather than the
    # multiplier THETA(13); theta_y = 0.18 / 0.28.
    e_snp_cyp3a4_rs35599367_cl <- 0.642857; label("Multiplicative effect of CYP3A4*22 carriage on clearance (unitless)") # Table 4 'Clearance CYP3A4*22' 0.18 L/h / 'Clearance' 0.28 L/h = 0.643 (36% lower; abstract '35% lower clearance'); control stream THETA(13) CYP3 as CL * THETA(13)**FL1

    # Relative bioavailability (Yu 2017 Eq. 3-5; control stream TVB, TVF1).
    lfdepot <- fixed(log(1)); label("Relative bioavailability at the 200 mg reference dose at the first dose (unitless)") # Yu 2017 Table 2 'rF200' 1 FIX; control stream TVB = 1 at DOSE = 200
    lemax <- fixed(log(1)); label("Maximum fractional decrease of relative bioavailability with dose, Emax (unitless)") # Yu 2017 Table 2 'Emax' 1 FIX; control stream THETA(9) VMAXDOSE (1 FIX)
    led50 <- log(466); label("Dose increment above 200 mg giving half the 200 mg relative bioavailability, ED50 (mg)") # Table 4 'Dose level with half of bioavailability at dose 200 mg' 466 mg (RSE 23%); control stream THETA(8) KMDOSE
    lfdepot_exp_inf <- log(0.49); label("Long-term asymptote of the time-dependent relative-bioavailability factor, 1 - Ad (unitless)") # Table 4 'Magnitude of decrease in relative bioavailability in time' 51% (RSE 27%) -> 1 - 0.51 = 0.49; control stream TVF1 = 1 - DCRP + DCRP * EXP(-LAMBDA * TIME), DCRP = THETA(12)
    lfdepot_exp_kdes <- log(0.00666667); label("First-order decay rate of the time-dependent relative-bioavailability factor (1/h)") # Table 4 'First-order decay constant' 0.16 day-1 (RSE 44%) -> 0.16 / 24 = 0.00666667 h-1; control stream LAMBDA = THETA(11) / 24

    # Ifosfamide sub-model (Kerbusch 2000 Table 1), fixed: Yu 2017 Section
    # 2.2.1 'a metabolic enzyme turn-over model that has been developed to
    # describe ifosfamide auto-induction [9] was used and fixed'. The source
    # data carried each patient's ifosfamide EBEs (ICLI, IKENZ, IIC50, IVI);
    # here the Kerbusch population values and IIV stand in for them.
    lcl_ifosfamide <- fixed(log(2.94)); label("Initial (pre-induction) ifosfamide clearance (L/h)") # Kerbusch 2000 Table 1 CL 2.94 L/h; control stream ICLI
    lvc_ifosfamide <- fixed(log(43.5)); label("Ifosfamide volume of distribution (L)") # Kerbusch 2000 Table 1 V 43.5 L; control stream IVI
    lkdeg <- fixed(log(0.0546)); label("Enzyme degradation rate constant K_enz,out (1/h)") # Kerbusch 2000 Table 1 K_enz,out 0.0546 h-1; control stream IKENZ
    lic50 <- fixed(log(30.7)); label("Ifosfamide concentration giving 50% inhibition of enzyme degradation (umol/L)") # Kerbusch 2000 Table 1 IC50 30.7 uM; control stream IIC50

    # Between-subject variability: Bins 2019 Table 4 final-model CV%.
    # omega^2 = log(1 + CV^2) (see the vignette for the scale check against the
    # published simulated trough IQRs). $OMEGA 2, 4, 5 and 8 (V2, KA2, ALAG3,
    # TVFR) are 0 FIX in the control stream and carry no eta here.
    etalka_fast ~ 1.103610 # Table 4 BSV 'Fast absorption rate constant' 142% (RSE 20%); log(1 + 1.42^2)
    etalcl ~ 0.070414 # Table 4 BSV 'Clearance' 27% (RSE 28%); log(1 + 0.27^2)
    etalvp ~ 0.664844 # Table 4 BSV 'Volume of distribution of peripheral compartment' 98% (RSE 18%); log(1 + 0.98^2)
    etalfdepot ~ 0.134652 # Table 4 BSV 'Relative bioavailability' 38% (RSE 15%); log(1 + 0.38^2)

    # Within-subject (inter-occasion) variability on relative bioavailability,
    # one eta per sampling occasion (control stream $OMEGA BLOCK(1) ... SAME x3).
    etaiov_fdepot_1 ~ 0.450679 # Table 4 WSV 'Relative bioavailability' 76% (RSE 22%); log(1 + 0.76^2)
    etaiov_fdepot_2 ~ fixed(0.450679) # same variance as occasion 1 ($OMEGA BLOCK(1) SAME)
    etaiov_fdepot_3 ~ fixed(0.450679) # same variance as occasion 1 ($OMEGA BLOCK(1) SAME)

    # Ifosfamide IIV (Kerbusch 2000 Table 1), fixed; omega^2 = log(1 + CV^2).
    etalcl_ifosfamide ~ fixed(0.058328) # Kerbusch 2000 Table 1 IIV CL 24.5%; log(1 + 0.245^2)
    etalvc_ifosfamide ~ fixed(0.053299) # Kerbusch 2000 Table 1 IIV V 23.4%; log(1 + 0.234^2)
    etalkdeg ~ fixed(0.050242) # Kerbusch 2000 Table 1 IIV K_enz,out 22.7%; log(1 + 0.227^2)

    # Residual error: Y = IPRED * (1 + ERR(1)) + ERR(2).
    propSd <- 0.08; label("Proportional residual error (fraction)") # Table 4 'Proportional residual error' 8.0% (RSE 11%)
    addSd <- 2.91; label("Additive residual error (mg/L)") # Table 4 'Additive residual error' 2.91 mg/L (RSE 27%)
  })

  model({
    # Inter-occasion variability on relative bioavailability (control stream
    # IOV1 by PKDAY).
    occ1 <- (OCC == 1)
    occ2 <- (OCC == 2)
    occ3 <- (OCC == 3)
    iov_fdepot <- occ1 * etaiov_fdepot_1 + occ2 * etaiov_fdepot_2 + occ3 * etaiov_fdepot_3

    # Pazopanib disposition parameters.
    ka_fast <- exp(lka_fast + etalka_fast)
    ka_slow <- exp(lka_slow)
    tlag <- exp(ltlag)
    cl <- exp(lcl + etalcl) * e_snp_cyp3a4_rs35599367_cl^SNP_CYP3A4_RS35599367
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp + etalvp)

    # Relative bioavailability (control stream TVB, TVF1, TVF; Yu 2017
    # Eq. 3-5). Dose part: 1 at 200 mg, Emax decrease above it.
    emax <- exp(lemax)
    ed50 <- exp(led50)
    dose_excess <- DOSE_PAZOPANIB_MG - 200
    fdepot_dose <- exp(lfdepot) * (1 - emax * dose_excess / (ed50 + dose_excess))
    # Time part: decays from 1 at the first dose to fdepot_exp_inf.
    fdepot_exp_inf <- exp(lfdepot_exp_inf)
    fdepot_exp_kdes <- exp(lfdepot_exp_kdes)
    fdepot_time <- fdepot_exp_inf + (1 - fdepot_exp_inf) * exp(-fdepot_exp_kdes * T_FIRSTDOSE)
    fdepot <- fdepot_dose * fdepot_time * exp(etalfdepot + iov_fdepot)
    # Split over the fast and slow routes (control stream F1, F3).
    frel <- expit(logitfrel)

    # Ifosfamide (Kerbusch 2000) and the shared enzyme pool.
    cl_ifosfamide <- exp(lcl_ifosfamide + etalcl_ifosfamide)
    vc_ifosfamide <- exp(lvc_ifosfamide + etalvc_ifosfamide)
    kdeg <- exp(lkdeg + etalkdeg)
    ic50 <- exp(lic50)
    Cc_ifosfamide <- central_ifosfamide / vc_ifosfamide

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka_fast * depot
    d/dt(depot2) <- -ka_slow * depot2
    d/dt(central) <- ka_fast * depot + ka_slow * depot2 - kel * enz_pool * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(central_ifosfamide) <- -cl_ifosfamide / vc_ifosfamide * enz_pool * central_ifosfamide
    d/dt(enz_pool) <- kdeg - kdeg * enz_pool * (1 - Cc_ifosfamide / (Cc_ifosfamide + ic50))
    enz_pool(0) <- 1

    # Every pazopanib dose is given to both depot and depot2 (control stream
    # compartments 1 and 3).
    f(depot) <- frel * fdepot
    f(depot2) <- (1 - frel) * fdepot
    alag(depot2) <- tlag

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
