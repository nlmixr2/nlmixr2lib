Tadayasu_2014_linagliptin <- function() {
  description <- "Two-compartment target-mediated drug disposition population PK model for linagliptin with quasi-equilibrium concentration-dependent binding to DPP-4 in both the central and peripheral compartments, coupled with an occupancy-based DPP-4-inhibition pharmacodynamic model (DPP-4 inhibition = Emax * Cbound/BMAX in the central compartment), in Japanese patients with type 2 diabetes mellitus (Tadayasu 2014 Table 3)."
  reference <- paste(
    "Tadayasu Y, Sarashina A, Tsuda Y, Tatami S, Friedrich C, Retlich S,",
    "Staab A, Takano M. Population pharmacokinetic/pharmacodynamic analysis",
    "of the DPP-4 inhibitor linagliptin in Japanese patients with type 2",
    "diabetes mellitus. J Pharm Pharm Sci. 2013;16(5):708-721.",
    "doi:10.18433/j3s304."
  )
  vignette <- "Tadayasu_2014_linagliptin"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL (linagliptin total plasma concentration; converted from the internal nmol/L scale via MW 472.54 g/mol); % of baseline (DPP-4 activity)")

  covariateData <- list()

  population <- list(
    species         = "human",
    n_subjects      = 36L,
    n_subjects_pk   = 36L,
    n_subjects_pd   = 36L,
    n_studies       = 1L,
    studies         = "Phase II trial of linagliptin in Japanese patients with T2DM (Horie 2011 [ref 14]); 28-day randomised double-blind placebo-controlled multiple-dose study with 0.5, 2.5, or 10 mg QD arms. The 0.5 mg arm failed the quasi-equilibrium assumption and was excluded from the final model.",
    age_range       = "40-69 years (range across both modelled dose groups; Tadayasu 2014 Methods)",
    age_median      = "60.2 years in the 2.5 mg arm; 59.1 years in the 10 mg arm (Tadayasu 2014 Methods)",
    bmi_range       = "18.4-34.4 kg/m^2 (range across both modelled dose groups)",
    bmi_median      = "26.0 kg/m^2 in the 2.5 mg arm; 23.8 kg/m^2 in the 10 mg arm",
    sex_female_pct  = 25,
    race_ethnicity  = c(Asian = 100),
    disease_state   = "Type 2 diabetes mellitus (T2DM); baseline HbA1c 7.1 +/- 0.5 % in the 2.5 mg arm and 7.2 +/- 0.9 % in the 10 mg arm; baseline fasting plasma glucose 154.7 +/- 25.1 mg/dL and 158.4 +/- 28.6 mg/dL respectively. Patients with hepatic, renal, neurological, cardiovascular, gastrointestinal, metabolic, hormonal disorders, hyperlipidemia, or hypertension were excluded.",
    dose_range      = "2.5 mg or 10 mg PO once daily for 28 days (final model excludes the 0.5 mg arm; Table 2)",
    regions         = "Japan (5 clinical sites)",
    co_medication   = "Individual antidiabetic treatment discontinued 14 days prior to first drug administration; no concomitant antidiabetic drugs during the trial.",
    n_observations  = "1018 linagliptin plasma concentrations + 1058 DPP-4 inhibition observations across 36 patients in the 2.5 mg and 10 mg arms used for analysis (Tadayasu 2014 Table 2).",
    notes           = "Demographics summarised in Tadayasu 2014 Methods (Patient population); modelled population includes only the 2.5 mg and 10 mg arms. The model structure (TMDD-QE with binding in both compartments) was carried over from the prior non-Japanese analysis of Friedrich 2011 [ref 18 in Tadayasu 2014]; all parameter values were re-estimated in the Japanese cohort."
  )

  ini({
    # ---- Structural PK parameters (Tadayasu 2014 Table 3) ----
    # F1 = 1 is set as the structural reference (Table 3 row F1, footnote a / paper Methods);
    # the IIV on F1 (etalfdepot) captures relative bioavailability dispersion around 1.
    lfdepot <- fixed(log(1));     label("Reference relative bioavailability (F1 = 1; fixed structural anchor)")        # Table 3 row F1
    lka     <- log(1.63);         label("First-order absorption rate constant Ka (1/h)")                                # Table 3 row KA
    lcl     <- log(121);          label("Apparent clearance of free linagliptin CL/F1 (L/h)")                            # Table 3 row CL/F1
    lvc     <- log(633);          label("Apparent central volume of distribution V2/F1 (L)")                             # Table 3 row V2/F1
    lq      <- log(73.0);         label("Apparent inter-compartmental clearance Q3/F1 (L/h)")                            # Table 3 row Q3/F1
    lvp     <- log(683);          label("Apparent peripheral volume of distribution V3/F1 (L)")                          # Table 3 row V3/F1
    lbmax   <- log(6.07);         label("Concentration of DPP-4 binding sites in the central compartment BMAX (nmol/L)") # Table 3 row BMAX
    lkd     <- log(0.108);        label("Linagliptin-DPP-4 dissociation constant KD (nmol/L)")                           # Table 3 row KD
    lamax2  <- log(534);          label("Apparent amount of DPP-4 binding sites in the peripheral compartment AMAX2/F1 (nmol)") # Table 3 row AMAX2/F1
    emax    <- 0.925;             label("Maximum fractional DPP-4 inhibition (Emax/100)")                                # Table 3 row EMAX (92.5 %)

    # ---- IIV (Tadayasu 2014 Table 3; CV % converted to log-normal variance via omega^2 = log(1 + CV^2)) ----
    # Final model uses diagonal OMEGA (the BMAX-CL block correlation R = 0.837 mentioned in Results was tested
    # but not implemented because it overestimated the variability).
    etalfdepot ~ 0.19727  # omega^2 = log(1 + 0.467^2) = 0.19727 (CV 46.7 %)  Table 3 row IIV in F1
    etalka     ~ 0.43288  # omega^2 = log(1 + 0.736^2) = 0.43288 (CV 73.6 %)  Table 3 row IIV in KA
    etalcl     ~ 0.38754  # omega^2 = log(1 + 0.688^2) = 0.38754 (CV 68.8 %)  Table 3 row IIV in CL
    etalbmax   ~ 0.01996  # omega^2 = log(1 + 0.142^2) = 0.01996 (CV 14.2 %)  Table 3 row IIV in BMAX

    # ---- Residual error (Tadayasu 2014 Table 3) ----
    # PK: paper coded as additive on log-transformed concentrations, which the authors equate
    # to a proportional error with CV 27.0 % on the untransformed scale (Methods, Pharmacokinetic model).
    propSd <- 0.270;        label("Proportional residual error for linagliptin plasma concentration (fraction)")  # Table 3 row 'Proportional residual variability PK'

    # PD: paper specifies the custom residual model Y = Yhat + (100 - Yhat) * epsilon with epsilon ~ N(0, sigma^2)
    # on DPP-4 inhibition (%). This is algebraically identical to a proportional residual on DPP-4 activity
    # (= 100 - inhibition), and is encoded here as prop() on the Dpp4Act output (% of baseline activity).
    propSd_Dpp4Act <- 0.201; label("Proportional residual error for DPP-4 activity (= 100 - inhibition); fraction") # Table 3 row 'Proportional residual variability PD'
  })

  model({
    # ---- Constants ----
    # Linagliptin (BI 1356) molecular weight: small molecule, xanthine-based DPP-4 inhibitor,
    # CAS 668270-12-0, C25H28N8O2 -> 472.54 g/mol. Same value used in Retlich_2015_linagliptin.
    mw_linag    <- 472.54
    nmol_per_mg <- 1e6 / mw_linag

    # ---- Individual PK parameters ----
    f1     <- exp(lfdepot + etalfdepot)
    ka     <- exp(lka     + etalka)
    cl     <- exp(lcl     + etalcl)
    vc     <- exp(lvc)
    q      <- exp(lq)
    vp     <- exp(lvp)
    bmax_c <- exp(lbmax   + etalbmax)
    kd     <- exp(lkd)
    amax2  <- exp(lamax2)

    # Apparent peripheral binding-site concentration (BMAX2 = AMAX2 / V3, nmol/L).
    bmax_p <- amax2 / vp

    # ---- Quasi-equilibrium free-vs-total drug algebra (Gibiansky 2008) ----
    # For each compartment with binding-site concentration B and dissociation constant KD:
    #   C_total = C_free + B * C_free / (KD + C_free)
    #   complex^2 - (C_total + B + KD) * complex + B * C_total = 0
    #   complex = 0.5 * ((C_total + B + KD) - sqrt((C_total + B + KD)^2 - 4 * B * C_total))
    #   C_free  = C_total - complex
    # Same algebra is used in the upstream Retlich_2015_linagliptin model.
    ctot_c    <- central / vc
    discr_c   <- (ctot_c + bmax_c + kd)^2 - 4 * bmax_c * ctot_c
    complex_c <- 0.5 * ((ctot_c + bmax_c + kd) - sqrt(discr_c))
    cfree_c   <- ctot_c - complex_c

    ctot_p    <- peripheral1 / vp
    discr_p   <- (ctot_p + bmax_p + kd)^2 - 4 * bmax_p * ctot_p
    complex_p <- 0.5 * ((ctot_p + bmax_p + kd) - sqrt(discr_p))
    cfree_p   <- ctot_p - complex_p

    # ---- ODE system (Tadayasu 2014 Figure 1) ----
    # Linear (non-specific) clearance and inter-compartmental flux act on the FREE drug
    # concentration in each compartment; only unbound linagliptin crosses membranes.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - cl * cfree_c - q * cfree_c + q * cfree_p
    d/dt(peripheral1) <-               q * cfree_c - q * cfree_p

    # ---- Bioavailability (mg dose to nmol amount, multiplied by relative F) ----
    f(depot) <- f1 * nmol_per_mg

    # ---- Observation 1: linagliptin total plasma concentration ----
    # Internal kinetics use nmol/L (matches Tadayasu 2014 Table 3 units of BMAX, KD, AMAX2,
    # and reported concentrations). The observed Cc is converted to ng/mL so that the
    # user-facing concentration is dimensionally consistent with the mg dose unit.
    # Conversion: c_nmolL * MW_linag (g/mol) / 1000 = ng/mL
    #   1 nmol/L * MW (g/mol) = MW * 1e-9 g/L = MW ng/L = MW / 1000 ng/mL.
    Cc <- ctot_c * mw_linag / 1000
    Cc ~ prop(propSd)

    # ---- Observation 2: DPP-4 activity and inhibition (Tadayasu 2014 Eq. for occupancy model) ----
    # Occupancy = fraction of central-compartment DPP-4 molecules bound by linagliptin
    #           = complex_c / BMAX  (equivalently cfree_c / (KD + cfree_c)).
    # DPP-4 inhibition (%) = 100 * Emax * occupancy.
    # The paper's residual model Y = Yhat + (100 - Yhat) * eps on inhibition is encoded as
    # a proportional residual on DPP-4 activity (= 100 - inhibition), which is algebraically
    # equivalent for a symmetric Gaussian residual.
    occupancy <- complex_c / bmax_c
    Dpp4Inh   <- 100 * emax * occupancy
    Dpp4Act   <- 100 - Dpp4Inh
    Dpp4Act ~ prop(propSd_Dpp4Act)
  })
}
