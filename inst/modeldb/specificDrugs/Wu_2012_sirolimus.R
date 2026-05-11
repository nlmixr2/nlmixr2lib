Wu_2012_sirolimus <- function() {
  description <- "Two-compartment population PK model for oral sirolimus with saturable Michaelis-Menten absorption in patients with advanced cancer (Wu 2012). Hematocrit power covariate on apparent oral clearance."
  reference <- "Wu K, Cohen EEW, House LK, Ramirez J, Zhang W, Ratain MJ, Bies RR. Nonlinear Population Pharmacokinetics of Sirolimus in Patients With Advanced Cancer. CPT Pharmacometrics Syst Pharmacol. 2012;1(11):e17. doi:10.1038/psp.2012.18"
  vignette <- "Wu_2012_sirolimus"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    HCT = list(
      description        = "Hematocrit (packed red-blood-cell volume fraction)",
      units              = "%",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline hematocrit; range 11.8-44.7% (median 35.1%) per Wu 2012 Table 2 footnote to Table 1. Negative power-form effect on apparent oral clearance: CL/F = th_CL * (35.1 / HCT)^th_HCT; higher hematocrit -> lower apparent CL, attributed to RBC sequestration of sirolimus (>95% of whole-blood drug resides in red cells).",
      source_name        = "HCT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 76L,
    n_observations = 563L,
    n_studies      = 4L,
    age_range      = "22-83 years",
    age_median     = "mean 57.7 years",
    weight_range   = "32.8-154.6 kg",
    weight_median  = "mean 79.76 kg",
    sex_female_pct = 48.7,
    race_ethnicity = "Not reported (single-centre cohort at The University of Chicago).",
    disease_state  = "Adults with advanced solid tumors enrolled in phase I trials of oral sirolimus.",
    dose_range     = "1-60 mg/week oral sirolimus, plus one trial of 4 mg once daily. Liquid formulation in trials 1, 3, and 4; tablet in trial 2. Drug formulation was not a significant covariate on absorption parameters.",
    regions        = "USA (single centre).",
    notes          = "Wu 2012 Table 2 baseline demographics: hematocrit 35.0% (11.8-44.7), hemoglobin 11.83 g/dL (7.5-15.7), creatinine 0.85 mg/dL (0.3-1.4), GFR 88.69 mL/min/BSA (42->120), albumin 3.80 g/dL (2.5-4.6). Hematocrit was the only significant covariate retained in backward elimination; in forward addition gender shrank V1/F to 0.75x the male value but did not survive backward elimination, and body weight modestly increased V2/F (exponent 0.676) without statistically significant OFV change. 16 of 808 raw samples (2.0%) were BLQ; M3 method was tested but not retained."
  )

  ini({
    # Structural parameters - Wu 2012 Table 1 final NONMEM estimates.
    # All clearances and volumes are apparent (X/F): bioavailability is folded
    # in because no IV data were collected. Full dose enters the depot
    # (intestinal-lumen compartment A1) and absorption proceeds via the
    # Michaelis-Menten term; F = 1 in the structural code is the standard
    # apparent-parameter convention, NOT a literature value.
    lcl <- log(12.9);  label("Apparent oral clearance CL1/F at reference hematocrit 35.1% (L/h)")  # Table 1 NONMEM theta1 = 12.9 L/h (16.3% RSE)
    lvc <- log(53.4);  label("Apparent central volume of distribution V1/F (L)")                    # Table 1 NONMEM V1/F = 53.4 L (38.0% RSE)
    lq  <- log(29.0);  label("Apparent intercompartmental clearance CL2/F (L/h)")                   # Table 1 NONMEM CL2/F = 29.0 L/h (8.17% RSE)
    lvp <- log(611);   label("Apparent peripheral volume of distribution V2/F (L)")                 # Table 1 NONMEM V2/F = 611 L (11.3% RSE)

    # Michaelis-Menten saturable absorption from intestinal lumen to central
    # compartment. The paper prints Vm units as "ug/l.h" in Table 1 but the
    # model equation d(A1)/dt = -Vm*A1/(Km+A1) requires mass/time (A1 and Km
    # are amounts in mg per the same Table 1 caption and equation text).
    # Reading Vm as 4.56 mg/h gives sensible absorption kinetics across the
    # 1-60 mg dose range (saturation onset near Km = 13.8 mg).
    lvmax <- log(4.56); label("Maximum oral absorption rate Vm (mg/h)")  # Table 1 NONMEM Vm = 4.56 (37.7% RSE); units inferred from equation mass balance (paper unit print is inconsistent)
    lkm   <- log(13.8); label("Amount in intestinal lumen at 50% of Vm (mg)")  # Table 1 NONMEM Km = 13.8 mg (50.3% RSE)

    # Covariate effect: CL1/F = th1 * (35.1 / HCT)^th2. Reference hematocrit is
    # the population median 35.1% (Table 1 footnote). Positive exponent means
    # CL DECREASES with INCREASING hematocrit (Discussion: HCT 35.1 -> 44.7
    # drops CL from 12.9 to 12.4 L/h, ~4% change).
    e_hct_cl <- 0.14; label("Power exponent on (35.1/HCT) for CL1/F (unitless)")  # Table 1 NONMEM theta2 = 0.14 (55.4% RSE)

    # IIV - Wu 2012 Methods "Eq. 1: P_ij = PTV_j * exp(eta_ij)" (log-normal).
    # Table 1 reports IIV as a percentage; converted via omega^2 = log(1 + CV^2).
    # IIV on CL1/F and V1/F print identically at 52.4%; the table reports
    # diagonal entries only with no off-diagonal/correlation, so the etas are
    # encoded as independent (no $OMEGA BLOCK reported).
    etalcl ~ 0.2425   # 52.4% CV; omega^2 = log(1 + 0.524^2) = 0.2425
    etalvc ~ 0.2425   # 52.4% CV; omega^2 = log(1 + 0.524^2) = 0.2425
    etalq  ~ 0.4035   # 70.5% CV; omega^2 = log(1 + 0.705^2) = 0.4035
    etalvp ~ 0.0366   # 19.3% CV; omega^2 = log(1 + 0.193^2) = 0.0366

    # Residual error - Wu 2012 Methods: Cobs = Cpred*(1 + eps1) + eps2; Table 1
    # reports proportional SD = 2.17% and additive SD = 0.5 ng/mL. Proportional
    # is poorly identified (97.0% RSE) but is retained per the paper's final
    # model.
    propSd <- 0.0217; label("Proportional residual error (fraction)")        # Table 1 NONMEM proportional = 2.17% (97.0% RSE)
    addSd  <- 0.5;    label("Additive residual error (ng/mL)")                # Table 1 NONMEM additive = 0.5 ng/mL (35.5% RSE)
  })

  model({
    # Individual PK parameters (apparent values; F folded into CL/V).
    # Hematocrit centered on the population median 35.1%.
    cl   <- exp(lcl + etalcl) * (35.1 / HCT)^e_hct_cl
    vc   <- exp(lvc + etalvc)
    q    <- exp(lq + etalq)
    vp   <- exp(lvp + etalvp)
    vmax <- exp(lvmax)
    km   <- exp(lkm)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Wu 2012 ODE system (Results, page 2). depot = A1 (intestinal lumen,
    # amount in mg), central = A2 (apparent central, mg), peripheral1 = A3
    # (apparent peripheral, mg). Initial dose enters depot.
    d/dt(depot)       <- -vmax * depot / (km + depot)
    d/dt(central)     <-  vmax * depot / (km + depot) - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                                k12 * central - k21 * peripheral1

    # Concentration: A2 (mg) / V1 (L) = mg/L; * 1000 -> ng/mL (paper unit
    # used in Figure 4 and all reported observations).
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
