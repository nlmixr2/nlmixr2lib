Abuhelwa_2022_ketamine_hr <- function() {
  description <- "Population PK/PD model for the acute effect of racemic ketamine on heart rate (HR) in adults with treatment-refractory depression (Abuhelwa 2022). HR is an immediate (no-delay) Emax function of plasma ketamine rising from a baseline of 72.7 beats/min towards a maximum attainable HR fixed to 220 beats/min, with EC50 7580 ng/mL; between-subject variability is on the baseline and on EC50 (Table 3 and Results; the supplementary control stream instead places the second eta on the maximum attainable HR); norketamine has no effect. Ketamine and norketamine PK (two-compartment each, first-order s.c./i.m. absorption with a common ka and 64.4 percent bioavailability, all ketamine clearance forming norketamine, allometric scaling to 70 kg) are the final PK model of Abuhelwa_2022_ketamine, held fixed."
  reference <- "Abuhelwa AY, Somogyi AA, Loo CK, Glue P, Barratt DT, Foster DJR. Population Pharmacokinetics and Pharmacodynamics of the Therapeutic and Adverse Effects of Ketamine in Patients With Treatment-Refractory Depression. Clin Pharmacol Ther. 2022;112(3):720-729. doi:10.1002/cpt.2640."
  vignette <- "Abuhelwa_2022_ketamine"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL", hr = "beats/min")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling referenced to 70 kg on every clearance (exponent 0.75: CL, Q, CL/fm, Q/fm) and every volume (exponent 1: V1, V2, V1/fm, V2/fm) of ketamine and norketamine (Methods 'General modeling approach'; supplementary NONMEM $PK, WTS = 70).",
      source_name = "WEIGHT"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "ketamine", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ketamine", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ketamine", units = "ug", specimen = "plasma", verified = TRUE),
    central_nk = list(
      analyte = "norketamine",
      units = "ug (ketamine mass equivalents)",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1_nk = list(
      analyte = "norketamine",
      units = "ug (ketamine mass equivalents)",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 21L,
    n_studies = 2L,
    age_range = "29-69 years",
    age_median = "50 years (mean 49.5)",
    weight_range = "not reported",
    sex_female_pct = 33.3,
    disease_state = "Treatment-refractory major depression (MADRS >= 20 before each treatment)",
    dose_range = "Racemic ketamine 0.1-0.5 mg/kg single doses in a within-subject dose-titration design, by i.v., s.c. or i.m. injection; doses at least 1 week apart",
    regions = "Australia (University of New South Wales, Sydney)",
    genotype = "CYP2B6 *6/*6 : *1/*6 or *5 : *1/*1 = 2 : 8 : 11",
    n_observations = "278 ketamine and 278 norketamine plasma concentrations (Table S1); 14 percent of concentrations were below the 0.25 ng/mL LLOQ and were excluded (M1 method)",
    notes = "Demographics from Abuhelwa 2022 Table 1 (7 female, 14 male). Data are from 61 ketamine treatment sessions of the Loo 2016 multiple-route dose-titration crossover study, including the Lai 2014 i.v. pilot. Body weight was used as a covariate but its distribution is not reported. Heart rate was measured before each treatment and at 5, 10, 30, 60 and 240 min after it (322 readings in 19 subjects, Table S1)."
  )

  ini({
    # PK block: the final PK model of Table 2, held FIXED while the PD
    # parameters were estimated (sequential PPP&D approach; every PK $THETA
    # and $OMEGA carries FIX in the supplementary PD control streams).
    # Ketamine (Table 2; supplementary NONMEM $THETA). Typical values for 70 kg.
    lka <- fixed(log(6.53)); label("First-order absorption rate constant, common to s.c. and i.m. (1/h)") # Table 2 'KA ketamine' = 6.53 h-1 (RSE 19.1%); one THETA for both routes
    lfdepot <- fixed(log(0.644)); label("Absolute bioavailability of s.c. and i.m. ketamine (fraction)") # Table 2 'F extravascular' = 0.644 (RSE 11.5%); F2 = F3 in the supplementary code
    lcl <- fixed(log(69.6)); label("Ketamine clearance at 70 kg (L/h)") # Table 2 'CL ketamine' = 69.6 L/h (RSE 10.5%)
    lvc <- fixed(log(79.3)); label("Ketamine central volume of distribution at 70 kg (L)") # Table 2 'V1 ketamine' = 79.3 L (RSE 18.7%)
    lq <- fixed(log(121)); label("Ketamine intercompartmental clearance at 70 kg (L/h)") # Table 2 'Q ketamine' = 121 L/h (RSE 14.3%)
    lvp <- fixed(log(87.4)); label("Ketamine peripheral volume of distribution at 70 kg (L)") # Table 2 'V2 ketamine' = 87.4 L (RSE 14.8%)

    # Norketamine apparent parameters (divided by fm, which is assumed to be 1;
    # Methods 'Pharmacokinetic modeling').
    lcl_nk <- fixed(log(30.3)); label("Norketamine apparent clearance CL/fm at 70 kg (L/h)") # Table 2 'CL norketamine/fm' = 30.3 L/h (RSE 31.2%)
    lvc_nk <- fixed(log(73.6)); label("Norketamine apparent central volume V1/fm at 70 kg (L)") # Table 2 'V1 norketamine/fm' = 73.6 L (RSE 12.4%)
    lq_nk <- fixed(log(54.2)); label("Norketamine apparent intercompartmental clearance Q/fm at 70 kg (L/h)") # Table 2 'Q norketamine/fm' = 54.2 L/h (RSE 11.3%)
    lvp_nk <- fixed(log(339)); label("Norketamine apparent peripheral volume V2/fm at 70 kg (L)") # Table 2 'V2 norketamine/fm' = 339 L (RSE 34.2%)

    # Allometric exponents, fixed (Methods 'General modeling approach').
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent of WT on ketamine CL and Q (unitless)") # Methods: 'clearances scaled to the exponent of 0.75'
    e_wt_vc_vp <- fixed(1); label("Allometric exponent of WT on ketamine V1 and V2 (unitless)") # Methods: 'volume to 1'
    e_wt_cl_nk <- fixed(0.75); label("Allometric exponent of WT on norketamine CL/fm (unitless)") # Methods: 'clearances scaled to the exponent of 0.75'
    e_wt_q_nk <- fixed(0.75); label("Allometric exponent of WT on norketamine Q/fm (unitless)") # Methods: 'clearances scaled to the exponent of 0.75'
    e_wt_vc_nk <- fixed(1); label("Allometric exponent of WT on norketamine V1/fm (unitless)") # Methods: 'volume to 1'
    e_wt_vp_nk <- fixed(1); label("Allometric exponent of WT on norketamine V2/fm (unitless)") # Methods: 'volume to 1'

    # IIV: Table 2 BSV %CV = sqrt(omega2) x 100, confirmed by the supplementary
    # $OMEGA values (e.g. sqrt(0.0341) = 0.185).
    etalcl ~ fixed(0.0341) # Table 2 BSV 18.5%CV; supplementary $OMEGA PPV_CLK = 0.0341
    etalvc ~ fixed(0.219) # Table 2 BSV 46.8%CV; supplementary $OMEGA PPV_V1K = 0.219
    etalvp ~ fixed(0.182) # Table 2 BSV 42.7%CV; supplementary $OMEGA PPV_V2K = 0.182
    etalcl_nk ~ fixed(0.424) # Table 2 BSV 65.1%CV; supplementary $OMEGA PPV_CLNK = 0.424
    etalvc_nk ~ fixed(0.217) # Table 2 BSV 46.6%CV; supplementary $OMEGA PPV_V1NK = 0.217

    # Residual error: Y = IPRED + sqrt((IPRED * prop)^2 + add^2) * EPS, SIGMA 1 FIX.
    propSd <- fixed(0.294); label("Proportional residual error, ketamine (fraction)") # Table 2 'RUVPROP-ketamine' = 0.294 (RSE 9.6%)
    addSd <- fixed(0.00001); label("Additive residual error, ketamine (ng/mL)") # Table 2 'RUVADD-ketamine' = 0.00001 Fix
    propSd_nk <- fixed(0.257); label("Proportional residual error, norketamine (fraction)") # Table 2 'RUVPROP-norketamine' = 0.257 (RSE 18.6%)
    addSd_nk <- fixed(3.88); label("Additive residual error, norketamine (ng/mL)") # Table 2 'RUVADD-norketamine' = 3.88 ng/mL (RSE 14.7%)

    # HR PD model (Table 3 'HR' column; supplementary 'Heart rate PKPD model'
    # control stream). EMAX is the maximum attainable heart rate (a LEVEL):
    # EFF = EBASE + (EMAX - EBASE) * C / (EC50 + C).
    lrbase_hr <- log(72.7); label("Baseline heart rate (beats/min)") # Table 3 HR 'BASE' = 72.7 beats/min (bootstrap 67-78)
    lec50_hr <- log(7580); label("Plasma ketamine concentration at half-maximal HR increase (ng/mL)") # Table 3 HR 'EC50' = 7580 ng/mL (bootstrap 2198-8990)
    lrmax_hr <- fixed(log(220)); label("Maximum attainable heart rate (beats/min)") # Table 3 HR 'Emax' = 220 Fix; Results: estimates above 220 bpm were implausible
    lhill_hr <- fixed(log(1)); label("Hill coefficient of the HR concentration-effect curve (unitless)") # Table 3 HR 'HILL' = 1 Fix

    # Table 3 PPV %CV = sqrt(omega2) x 100 (the convention the PK omegas confirm).
    # The second HR random effect is on EC50, as Table 3 ('PPV_EC50' row) and
    # the Results ('Both BP and HR PD models supported BSV on the baseline and
    # EC50 parameters') state. The supplementary HR control stream puts it on
    # EMAX instead (TVEMAX = THETA(19) * EXP(ETA(7))), but with omega2 = 2.43
    # that placement gives a 95th-percentile HR of about 155 beats/min at the
    # end of a 0.5 mg/kg 40-min infusion, against about 120 beats/min in the
    # paper's own simulated 90% band (Figures S6 and S7); EC50 placement gives
    # about 122.
    etalrbase_hr ~ 0.023104 # Table 3 HR 'PPV_BASE' = 15.2%CV; omega2 = 0.152^2
    etalec50_hr ~ 2.4336 # Table 3 HR 'PPV_EC50' = 156.0%CV; omega2 = 1.56^2

    addSd_hr <- 9.15; label("Additive residual error, HR (beats/min)") # Table 3 HR 'RUV_ADDPD' = 9.15 beats/min (bootstrap 6.25-12.7)
    propSd_hr <- fixed(0.0001); label("Proportional residual error, HR (fraction)") # Table 3 HR 'RUV_PROPPD' = 0.0001 Fix
  })

  model({
    # Individual parameters with allometric scaling to 70 kg (supplementary $PK).
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q <- exp(lq) * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp

    cl_nk <- exp(lcl_nk + etalcl_nk) * (WT / 70)^e_wt_cl_nk
    vc_nk <- exp(lvc_nk + etalvc_nk) * (WT / 70)^e_wt_vc_nk
    q_nk <- exp(lq_nk) * (WT / 70)^e_wt_q_nk
    vp_nk <- exp(lvp_nk) * (WT / 70)^e_wt_vp_nk

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    kel_nk <- cl_nk / vc_nk
    k12_nk <- q_nk / vc_nk
    k21_nk <- q_nk / vp_nk

    # i.v. doses go to central; s.c. and i.m. doses go to depot (the NONMEM
    # code keeps separate s.c. and i.m. depots but gives them the same KA and
    # F, so one depot is equivalent). All ketamine clearance forms
    # norketamine in ketamine mass units (fm = 1, no molecular-weight
    # conversion), exactly as in the supplementary $DES.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(central_nk) <- kel * central - kel_nk * central_nk - k12_nk * central_nk + k21_nk * peripheral1_nk
    d/dt(peripheral1_nk) <- k12_nk * central_nk - k21_nk * peripheral1_nk

    f(depot) <- exp(lfdepot)

    # ug / L = ng/mL
    Cc <- central / vc
    Cc_nk <- central_nk / vc_nk

    Cc ~ add(addSd) + prop(propSd)
    Cc_nk ~ add(addSd_nk) + prop(propSd_nk)

    # Immediate Emax towards a maximum attainable HR (supplementary $ERROR:
    # EFF = EBASE + (EMAX - EBASE) * C / (EC50 + C)).
    rbase_hr <- exp(lrbase_hr + etalrbase_hr)
    rmax_hr <- exp(lrmax_hr)
    ec50_hr <- exp(lec50_hr + etalec50_hr)
    hill_hr <- exp(lhill_hr)
    cdrv_hr <- max(Cc, 0)
    hr <- rbase_hr + (rmax_hr - rbase_hr) * cdrv_hr^hill_hr / (ec50_hr^hill_hr + cdrv_hr^hill_hr)

    hr ~ add(addSd_hr) + prop(propSd_hr)
  })
}
