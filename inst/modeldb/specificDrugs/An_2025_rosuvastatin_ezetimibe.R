An_2025_rosuvastatin_ezetimibe <- function() {
  description <- "Joint population PK/PD model of co-administered rosuvastatin and ezetimibe with enterohepatic recirculation, fitted to a two-part open-label multiple-dose crossover drug-interaction study in 50 healthy Korean male volunteers (An 2025 Table 2). Rosuvastatin has first-order absorption and two-compartment disposition. Total ezetimibe (unchanged ezetimibe plus its phenolic glucuronide, which together are the measured analyte) has a four-compartment structure: a gastrointestinal / absorption compartment, central and peripheral compartments, and a gallbladder reservoir. Drug moves from the ezetimibe central compartment into the gallbladder continuously at kbm, and is released back into the GASTROINTESTINAL compartment (not the central compartment) during three 0.75 h post-prandial windows, from which it is re-absorbed at ka - this GI-linked return is the structure the authors found superior on goodness-of-fit. The two PK models drive one shared LDL-cholesterol indirect-response compartment through independent, multiplicative (Bliss-independent) inhibition of LDL-C production, with Imax fixed at 1 and the Hill coefficient fixed at 1; there is no PK or PD interaction term. The meal gate is anchored to TIME AFTER DOSE, matching the protocol's standard meals 4, 10 and 24 h after an administration, so a single-dose event table reproduces the paper's Equations 3-6 indicator exactly. No covariate was significant on any PK or PD parameter; the screened covariates are recorded in covariatesDataExcluded. Concentrations are in ng/mL and LDL-C in mg/dL, so model() scales amount/volume by 1000 to convert mg/L to ng/mL."
  reference <- paste(
    "An H, Shin D.",
    "Population pharmacokinetics and pharmacodynamics with enterohepatic",
    "recirculation of co-medication of rosuvastatin and ezetimibe.",
    "Drug Des Devel Ther. 2025;19:4775-4787.",
    "doi:10.2147/DDDT.S522863.",
    sep = " "
  )
  vignette <- "An_2025_rosuvastatin_ezetimibe"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # An 2025, "Covariate Analysis and Model Evaluation" and both Results
  # sections: age, weight, serum creatinine, albumin, ALP, ALT, AST and GGT
  # were screened by univariate testing followed by stepwise forward addition
  # and backward elimination, on both the PK and the PD parameters. None was
  # retained -- "No covariate effect was significant for the estimated PK
  # parameters, possibly because the volunteers were all young, healthy males"
  # and "Due to the homogeneity of the study population, no significant
  # covariates were identified". The final model therefore has no covariates.
  covariateData <- list()

  covariatesDataExcluded <- list(
    AGE   = list(description = "Age at screening", units = "years", type = "continuous",
                 notes = "Screened (centered on the mean) in An 2025 but not retained in the final PK or PD model."),
    WT    = list(description = "Body weight at screening", units = "kg", type = "continuous",
                 notes = "Screened (centered on the mean) in An 2025 but not retained; enrolment required weight within +/-20% of ideal body weight, so the range was narrow."),
    CREAT = list(description = "Serum creatinine at screening", units = "mg/dL", type = "continuous",
                 notes = "Screened in An 2025 but not retained; subjects with Cockcroft-Gault creatinine clearance below 80 mL/min were excluded."),
    ALB   = list(description = "Serum albumin at screening", units = "g/dL", type = "continuous",
                 notes = "Screened in An 2025 but not retained. An 2025 Table 1 reports albumin in g/dL (US convention), not the canonical g/L."),
    ALP   = list(description = "Serum alkaline phosphatase at screening", units = "U/L", type = "continuous",
                 notes = "Screened in An 2025 but not retained."),
    ALT   = list(description = "Serum alanine aminotransferase at screening", units = "U/L", type = "continuous",
                 notes = "Screened in An 2025 but not retained."),
    AST   = list(description = "Serum aspartate aminotransferase at screening", units = "U/L", type = "continuous",
                 notes = "Screened in An 2025 but not retained."),
    GGT   = list(description = "Serum gamma-glutamyltransferase at screening", units = "U/L", type = "continuous",
                 notes = "Screened in An 2025 but not retained.")
  )

  compartmentData <- list(
    depot_rosuvastatin       = list(analyte = "rosuvastatin", units = "mg", specimen = "administration site", verified = TRUE),
    central_rosuvastatin     = list(analyte = "rosuvastatin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_rosuvastatin = list(analyte = "rosuvastatin", units = "mg", specimen = "plasma", verified = TRUE),
    depot_ezetimibe          = list(analyte = "total ezetimibe (ezetimibe + ezetimibe phenolic glucuronide)", units = "mg", specimen = "administration site", verified = TRUE),
    central_ezetimibe        = list(analyte = "total ezetimibe (ezetimibe + ezetimibe phenolic glucuronide)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_ezetimibe    = list(analyte = "total ezetimibe (ezetimibe + ezetimibe phenolic glucuronide)", units = "mg", specimen = "plasma", verified = TRUE),
    gallbladder_ezetimibe    = list(analyte = "total ezetimibe (ezetimibe + ezetimibe phenolic glucuronide)", units = "mg", specimen = "bile", verified = TRUE),
    ldl                      = list(analyte = "low-density lipoprotein cholesterol", units = "mg/dL", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species       = "human",
    n_subjects    = 50L,
    n_studies     = 1L,
    age_range     = "19-45 years by protocol; median 24 (range 19-33) in Part A and 24 (19-37) in Part B (An 2025 Table 1).",
    weight_range  = "Median 69.1 kg (55.0-88.0) in Part A and 68.8 kg (44.8-84.0) in Part B (An 2025 Table 1).",
    sex_female_pct = 0,
    race_ethnicity = "Korean.",
    disease_state = "Healthy male volunteers with no clinically significant medical history; creatinine clearance at or above 80 mL/min by Cockcroft-Gault.",
    dose_range    = "Rosuvastatin 20 mg and/or ezetimibe 10 mg once daily for 7 days, as monotherapy or co-therapy.",
    regions       = "Republic of Korea (Gachon University Gil Medical Center).",
    notes         = "An 2025 Methods and Table 1. Two-part, open-label, multiple-dose, two-treatment, two-period, two-sequence crossover drug-interaction study (ClinicalTrials.gov NCT02289430, IRB GCIRB2014-324). 56 subjects enrolled (28 per part); 50 (25 per part) contributed to the population analysis, giving 25 rosuvastatin monotherapy, 25 ezetimibe monotherapy and 50 co-therapy concentration-time profiles. Part A compared rosuvastatin 20 mg with rosuvastatin 20 mg + ezetimibe 10 mg; Part B compared ezetimibe 10 mg with the same combination; 14-day washout between periods. Steady-state PK sampling on day 7 at 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 8, 12, 24, 48 and 72 h post-dose. Standard meals were given 4, 10 and 24 h after the final dose, and these define the gallbladder-emptying windows. LDL-C was measured pre-dose on day 1 (baseline) and 24 h after the final dose. Concentrations below the LLOQ (8.38% of rosuvastatin and 1.25% of ezetimibe samples, all after 48 h) were treated as missing. Estimation used the SAEM algorithm in Monolix 2024R1."
  )

  ini({
    # =====================================================================
    # Rosuvastatin: first-order absorption, two-compartment disposition.
    # All values are APPARENT (oral data only; the paper writes Cl_R, Vc_R
    # etc. without an explicit /F but the study has no intravenous arm).
    # An 2025 Table 2, "PK Model for Rosuvastatin" block.
    # =====================================================================
    lka_rosuvastatin <- log(0.21);   label("Rosuvastatin first-order absorption rate constant ka,R (1/h)")            # An 2025 Table 2: ka,R = 0.21 1/h, RSE 9.3%
    lcl_rosuvastatin <- log(92.27);  label("Apparent rosuvastatin clearance Cl,R/F (L/h)")                            # An 2025 Table 2: Cl,R = 92.27 L/h, RSE 8.5%
    lvc_rosuvastatin <- log(222.23); label("Apparent rosuvastatin central volume Vc,R/F (L)")                         # An 2025 Table 2: Vc,R = 222.23 L, RSE 13.0%
    lq_rosuvastatin  <- log(24.16);  label("Apparent rosuvastatin inter-compartmental clearance Q,R/F (L/h)")         # An 2025 Table 2: Q,R = 24.16 L/h, RSE 12.9%
    lvp_rosuvastatin <- log(650.71); label("Apparent rosuvastatin peripheral volume Vp,R/F (L)")                      # An 2025 Table 2: Vp,R = 650.71 L, RSE 24.1%

    # =====================================================================
    # Total ezetimibe: four-compartment model with enterohepatic
    # recirculation (An 2025 Equations 3-6 and Table 2, "PK Model for
    # Ezetimibe" block). "Total ezetimibe" is the measured analyte: the sum
    # of unchanged ezetimibe and ezetimibe phenolic glucuronide, the latter
    # accounting for 80-90% of the total and having comparable potency.
    # =====================================================================
    lka_ezetimibe <- log(0.64);   label("Total ezetimibe first-order absorption rate constant ka,E (1/h)")            # An 2025 Table 2: ka,E = 0.64 1/h, RSE 6.5%
    lcl_ezetimibe <- log(20.04);  label("Apparent total ezetimibe clearance Cl,E/F (L/h)")                            # An 2025 Table 2: Cl,E = 20.04 L/h, RSE 6.7%
    lvc_ezetimibe <- log(31.98);  label("Apparent total ezetimibe central volume Vc,E/F (L)")                         # An 2025 Table 2: Vc,E = 31.98 L, RSE 9.8%
    lq_ezetimibe  <- log(44.53);  label("Apparent total ezetimibe inter-compartmental clearance Q,E/F (L/h)")         # An 2025 Table 2: Q,E = 44.53 L/h, RSE 8.9%
    lvp_ezetimibe <- log(363.06); label("Apparent total ezetimibe peripheral volume Vp,E/F (L)")                      # An 2025 Table 2: Vp,E = 363.06 L, RSE 9.5%

    # Enterohepatic recirculation. kbm is the paper's kb,E (central ->
    # gallbladder, continuous and independent of food); kehc is the paper's
    # ke,E (gallbladder -> gastrointestinal compartment, gated by the meal
    # switch GBE). ke,E was FIXED as the reciprocal of the 0.75 h bile-release
    # duration rather than estimated -- Table 2 prints "Fix" in its RSE cell.
    lkbm_ezetimibe  <- log(0.013);      label("Total ezetimibe central-to-gallbladder biliary transfer rate constant kb,E (1/h)")  # An 2025 Table 2: kb,E = 0.013 1/h, RSE 5.3%
    lkehc_ezetimibe <- fixed(log(1.33)); label("Total ezetimibe gallbladder release rate constant ke,E (1/h)")                     # An 2025 Table 2: ke,E = 1.33 1/h, "Fix"; Results text "fixed as the inverse of duration, 1.33 h-1 = (0.75 h)-1"

    # Meal gate. An 2025 defines GBE = 1 for 0.75 h after each of three meals
    # given 4, 10 and 24 h after an administration, and 0 otherwise. These are
    # protocol constants, not estimated quantities, so they are fixed. Unlike
    # the Keunecke 2020 regorafenib model, whose tmeal values are clock times
    # since midnight, these are TIMES AFTER DOSE -- model() reads them against
    # tad(), so no midnight-origin event table is required.
    tmeal1 <- fixed(4);    label("Time after dose of the first meal triggering gallbladder emptying (h)")   # An 2025 Methods: "A standard meal was provided at 4, 10, and 24 h after the final dose"
    tmeal2 <- fixed(10);   label("Time after dose of the second meal triggering gallbladder emptying (h)")  # An 2025 Methods: meals at 4, 10 and 24 h after the final dose
    tmeal3 <- fixed(24);   label("Time after dose of the third meal triggering gallbladder emptying (h)")   # An 2025 Methods: meals at 4, 10 and 24 h after the final dose
    dge    <- fixed(0.75); label("Duration of each gallbladder-emptying window (h)")                        # An 2025 Results: "the duration of bile release in each EHC cycle was set to 0.75 h"

    # =====================================================================
    # LDL-C indirect response (An 2025 Equation 1 and Table 2, "PD Model for
    # LDL" block). Both drugs inhibit LDL-C production multiplicatively and
    # independently; Imax is fixed at 1 ("The model assumes Imax = 1 (full
    # inhibitory effect), which was fixed rather than estimated") and the Hill
    # coefficient at 1 ("Hill coefficient ... is simply fixed at 1").
    # kout is NOT tabulated: the paper reports the baseline and kin, and kout
    # follows from the drug-free steady-state identity baseline = kin / kout,
    # which model() applies per individual.
    # =====================================================================
    lbase <- log(92.2); label("Typical baseline (drug-free) LDL-C (mg/dL)")                                  # An 2025 Table 2: Baseline LDL = 92.2 mg/dL, RSE 3.9%
    lkin  <- log(1.9);  label("Zero-order LDL-C production rate kin (mg/dL/h)")                              # An 2025 Table 2: kin = 1.9 mg/dL*h, RSE 22.2%
    lic50_rosuvastatin <- log(4.6);  label("Rosuvastatin concentration giving 50% inhibition of LDL-C production IC50,R (ng/mL)")     # An 2025 Table 2: IC50,R = 4.6 ng/mL, RSE 7.8%
    lic50_ezetimibe    <- log(36.9); label("Total ezetimibe concentration giving 50% inhibition of LDL-C production IC50,E (ng/mL)")  # An 2025 Table 2: IC50,E = 36.9 ng/mL, RSE 11.8%

    # =====================================================================
    # Inter-individual variability. Methods: "All parameters were assumed to
    # follow a log-normal distribution, and an exponential error model was
    # used for interindividual variability", i.e. P_i = P * exp(eta_i).
    # An 2025 Table 2 tabulates the random effects under a column headed
    # "CV%", and the table's own abbreviation footnote defines "CV,
    # coefficient of variation", so each entry is read as a coefficient of
    # variation and converted to the log-scale variance ini() expects via
    # omega^2 = log(1 + CV^2). The alternative reading (the printed number is
    # 100 * omega, i.e. the log-scale SD) changes every value by under 5%
    # except omega Vp,R, where 108.25 gives 0.776 as a CV against 1.172 as an
    # SD; see the vignette Errata.
    #
    # Two rows of Table 2 are mislabelled in the source: the fifth
    # rosuvastatin row and the fifth ezetimibe row both print "omega Vc"
    # a second time. They are the peripheral-volume random effects -- the
    # central-volume rows appear immediately above them (49.39 and 43.62), the
    # parameter blocks list Vc then Q then Vp in that order, and the ezetimibe
    # correlation block names Vp,E as a correlated random effect, which
    # requires an omega Vp,E to exist.
    # =====================================================================

    # An 2025 Table 2: omega ka,R = 29.21 CV% (RSE 22.3%) and omega Cl,R =
    # 43.52 CV% (RSE 14.4%), with correlation ka,R - Cl,R = -0.87 (RSE 13.9%).
    # The off-diagonal is corr * sd1 * sd2 on the log scale.
    etalka_rosuvastatin + etalcl_rosuvastatin ~ c(
      0.0818771,
      -0.1036777, 0.1734482
    )
    # An 2025 Table 2: omega Vc,R = 49.39 CV% (RSE 17.2%); log(1 + 0.4939^2)
    etalvc_rosuvastatin ~ 0.2182815
    # An 2025 Table 2: omega Q,R = 31.95 CV% (RSE 34.3%); log(1 + 0.3195^2)
    etalq_rosuvastatin ~ 0.0971995
    # An 2025 Table 2: omega Vp,R = 108.25 CV% (RSE 21.2%); log(1 + 1.0825^2)
    etalvp_rosuvastatin ~ 0.7755592

    # An 2025 Table 2: omega Cl,E = 33.87 CV% (RSE 14.6%), omega Q,E = 41.8
    # CV% (RSE 17.1%), omega Vp,E = 46.13 CV% (RSE 16.1%), with correlations
    # Cl,E - Q,E = 0.74 (RSE 15.7%), Cl,E - Vp,E = 0.58 (RSE 25.5%) and
    # Q,E - Vp,E = 0.74 (RSE 17.2%). The resulting 3x3 correlation matrix has
    # determinant 0.204 and is positive definite.
    etalcl_ezetimibe + etalq_ezetimibe + etalvp_ezetimibe ~ c(
      0.1086012,
      0.0978603, 0.1610332,
      0.0839547, 0.1304335, 0.1929298
    )
    # An 2025 Table 2: omega ka,E = 29.27 CV% (RSE 17.5%); log(1 + 0.2927^2)
    etalka_ezetimibe ~ 0.0822003
    # An 2025 Table 2: omega Vc,E = 43.62 CV% (RSE 18.4%); log(1 + 0.4362^2)
    etalvc_ezetimibe ~ 0.1741805
    # An 2025 Table 2: omega kb,E = 11.93 CV% (RSE 34.0%); log(1 + 0.1193^2)
    etalkbm_ezetimibe ~ 0.0141322

    # An 2025 Table 2: omega baseline LDL = 26.8 CV% (RSE 10.8%); log(1 + 0.268^2)
    etalbase ~ 0.0693619
    # An 2025 Table 2: omega kin = 45.2 CV% (RSE 48.8%); log(1 + 0.452^2)
    etalkin ~ 0.1859018
    # An 2025 Table 2: omega IC50,R = 38.4 CV% (RSE 22.6%); log(1 + 0.384^2)
    etalic50_rosuvastatin ~ 0.1375473
    # An 2025 Table 2: omega IC50,E = 29.6 CV% (RSE 59.4%); log(1 + 0.296^2)
    etalic50_ezetimibe ~ 0.0839881

    # =====================================================================
    # Residual error. An 2025 Equation 2 is the Monolix "combined1" form
    # C_ij = f_ij + (a + b * f_ij) * e_ij, so a is an additive SD in the
    # units of the observation and b is a proportional SD (fraction).
    # =====================================================================
    addSd_rosuvastatin  <- 0.17; label("Rosuvastatin additive residual error a,R (ng/mL)")            # An 2025 Table 2: a,R = 0.17, RSE 22.3%
    propSd_rosuvastatin <- 0.20; label("Rosuvastatin proportional residual error b,R (fraction)")     # An 2025 Table 2: b,R = 0.20, RSE 3.7%
    # Ezetimibe carries no proportional term: "The residual variability was
    # described using an additive error model with zero proportional error
    # (b = 0) in Equation 2", and Table 2 lists a,E alone.
    addSd_ezetimibe     <- 0.33; label("Total ezetimibe additive residual error a,E (ng/mL)")         # An 2025 Table 2: a,E = 0.33, RSE 2.7%
    addSd_ldl              <- 0.73; label("LDL-C additive residual error a,LDL (mg/dL)")                 # An 2025 Table 2: a,LDL = 0.73, RSE 0.5%
    propSd_ldl             <- 0.12; label("LDL-C proportional residual error b,LDL (fraction)")          # An 2025 Table 2: b,LDL = 0.12, RSE 0.03%
  })

  model({
    # 1. Individual parameters. Rosuvastatin (An 2025 Table 2).
    ka_rosuvastatin <- exp(lka_rosuvastatin + etalka_rosuvastatin)
    cl_rosuvastatin <- exp(lcl_rosuvastatin + etalcl_rosuvastatin)
    vc_rosuvastatin <- exp(lvc_rosuvastatin + etalvc_rosuvastatin)
    q_rosuvastatin  <- exp(lq_rosuvastatin + etalq_rosuvastatin)
    vp_rosuvastatin <- exp(lvp_rosuvastatin + etalvp_rosuvastatin)

    # 2. Individual parameters. Total ezetimibe (An 2025 Table 2).
    ka_ezetimibe   <- exp(lka_ezetimibe + etalka_ezetimibe)
    cl_ezetimibe   <- exp(lcl_ezetimibe + etalcl_ezetimibe)
    vc_ezetimibe   <- exp(lvc_ezetimibe + etalvc_ezetimibe)
    q_ezetimibe    <- exp(lq_ezetimibe + etalq_ezetimibe)
    vp_ezetimibe   <- exp(lvp_ezetimibe + etalvp_ezetimibe)
    kbm_ezetimibe  <- exp(lkbm_ezetimibe + etalkbm_ezetimibe)
    kehc_ezetimibe <- exp(lkehc_ezetimibe)

    # 3. Individual PD parameters. kout is derived, not estimated: An 2025
    #    reports the baseline and kin, and the drug-free steady state of
    #    Equation 1 (all C = 0, d(LDL-C)/dt = 0) gives baseline = kin / kout.
    #    Deriving it per individual keeps each subject's simulated profile
    #    starting at that subject's own baseline.
    base <- exp(lbase + etalbase)
    kin  <- exp(lkin + etalkin)
    kout <- kin / base
    ic50_rosuvastatin <- exp(lic50_rosuvastatin + etalic50_rosuvastatin)
    ic50_ezetimibe    <- exp(lic50_ezetimibe + etalic50_ezetimibe)

    # 4. Gallbladder-emptying switch, the GBE indicator of An 2025 Equations 3
    #    and 6: 1 while bile is being released, 0 otherwise. The paper places
    #    the three release periods 4, 10 and 24 h AFTER an administration, so
    #    the gate is read against time after dose rather than absolute time.
    #    With once-daily dosing tad() never reaches tmeal3 within an interval,
    #    so the 4 h and 10 h meals recur after every dose and the 24 h window
    #    opens only once the final dose is no longer followed by another -
    #    exactly the protocol. For a single-dose event table tad() equals t
    #    and the gate reduces to the paper's literal {4, 10, 24} h indicator.
    tafterdose <- tad()
    gbe <- (tafterdose >= tmeal1) * (tafterdose < tmeal1 + dge) +
           (tafterdose >= tmeal2) * (tafterdose < tmeal2 + dge) +
           (tafterdose >= tmeal3) * (tafterdose < tmeal3 + dge)

    # 5. Rosuvastatin PK: first-order absorption into a two-compartment
    #    disposition (An 2025 Figure 2, left half).
    d/dt(depot_rosuvastatin)       <- -ka_rosuvastatin * depot_rosuvastatin
    d/dt(central_rosuvastatin)     <-  ka_rosuvastatin * depot_rosuvastatin -
                                        (cl_rosuvastatin / vc_rosuvastatin) * central_rosuvastatin -
                                        (q_rosuvastatin / vc_rosuvastatin) * central_rosuvastatin +
                                        (q_rosuvastatin / vp_rosuvastatin) * peripheral1_rosuvastatin
    d/dt(peripheral1_rosuvastatin) <-  (q_rosuvastatin / vc_rosuvastatin) * central_rosuvastatin -
                                        (q_rosuvastatin / vp_rosuvastatin) * peripheral1_rosuvastatin

    # 6. Total ezetimibe PK, An 2025 Equations 3-6 verbatim. depot_ezetimibe
    #    is the paper's gastrointestinal (GI) tract compartment A_GI,E: it
    #    receives the oral dose, empties into central at ka,E, and is
    #    RE-FILLED by the gallbladder during each emptying window. The
    #    paper's k23,E / k32,E are written here as the equivalent
    #    Q,E / Vc,E and Q,E / Vp,E, which is the parameterisation Table 2
    #    reports.
    d/dt(depot_ezetimibe)       <- -ka_ezetimibe * depot_ezetimibe +
                                     gbe * kehc_ezetimibe * gallbladder_ezetimibe
    d/dt(central_ezetimibe)     <-  ka_ezetimibe * depot_ezetimibe -
                                     (cl_ezetimibe / vc_ezetimibe) * central_ezetimibe -
                                     (q_ezetimibe / vc_ezetimibe) * central_ezetimibe -
                                     kbm_ezetimibe * central_ezetimibe +
                                     (q_ezetimibe / vp_ezetimibe) * peripheral1_ezetimibe
    d/dt(peripheral1_ezetimibe) <-  (q_ezetimibe / vc_ezetimibe) * central_ezetimibe -
                                     (q_ezetimibe / vp_ezetimibe) * peripheral1_ezetimibe
    d/dt(gallbladder_ezetimibe) <-  kbm_ezetimibe * central_ezetimibe -
                                     gbe * kehc_ezetimibe * gallbladder_ezetimibe

    # 7. Plasma concentrations. Amounts are in mg and volumes in L, so the
    #    ratio is mg/L; the factor 1000 converts to the ng/mL in which An 2025
    #    reports concentrations and both IC50 values.
    Cc_rosuvastatin <- 1000 * central_rosuvastatin / vc_rosuvastatin
    Cc_ezetimibe    <- 1000 * central_ezetimibe / vc_ezetimibe

    # 8. LDL-C indirect response, An 2025 Equation 1. The two inhibition
    #    terms multiply, which is the Bliss-independence form the Discussion
    #    names; Imax and the Hill coefficient are both fixed at 1, so each
    #    term is the plain 1 - C/(IC50 + C). The state starts at the
    #    individual drug-free baseline.
    d/dt(ldl) <- kin *
                   (1 - Cc_rosuvastatin / (ic50_rosuvastatin + Cc_rosuvastatin)) *
                   (1 - Cc_ezetimibe / (ic50_ezetimibe + Cc_ezetimibe)) -
                 kout * ldl
    ldl(0) <- base

    # 9. Observations, An 2025 Equation 2.
    Cc_rosuvastatin ~ add(addSd_rosuvastatin) + prop(propSd_rosuvastatin)
    Cc_ezetimibe    ~ add(addSd_ezetimibe)
    ldl             ~ add(addSd_ldl) + prop(propSd_ldl)
  })
}
