"Ait-Oudhia_2016_sunitinib" <- function() {
  description <- paste0(
    "Joint population PK/PD model for sunitinib and its equipotent active ",
    "metabolite SU12662 in adults with advanced hepatocellular carcinoma (HCC) ",
    "receiving 37.5 mg sunitinib PO QD. Parent drug and metabolite each ",
    "follow a 2-compartment oral PK structure with first-order absorption; ",
    "each oral sunitinib dose deposits Dose into the parent depot and ",
    "fM * Dose (fM = 0.21 fixed, Houk 2009) into the SU12662 depot. The ",
    "active free (unbound) drug concentration ACub = (1 - fb_D) * Cc + ",
    "(1 - fb_M) * Cc_su12662 (fb_D = 0.9, fb_M = 0.95 fixed, free fractions ",
    "0.1 and 0.05) inhibits the zero-order production rate of plasma sVEGFR2, ",
    "captured with an indirect-response model dsVEGFR2/dt = kin / (1 + alpha ",
    "* INH) - kout * sVEGFR2 with INH = ACub / (kd + ACub) (kd = 4 ug/L ",
    "fixed, Mendel 2003) and kin = R0 * kout. Tumor volume follows a first-",
    "order growth dTG/dt = kg * (1 - H(t)) * TG with kg derived from baseline ",
    "tumor volume by kg = ln(2) / (114 * TG0^0.14) (Taouli 2005) and H(t) = ",
    "Imax * dsVEGFR2 / (dsVEGFR2 + dIC50) with Imax = 1 fixed and dsVEGFR2 = ",
    "R0 - sVEGFR2(t). The paper reports a significant covariate effect of ",
    "the DCE-MRI volume-transfer constant Ktrans on dIC50 (power coefficient ",
    "2.12) but the cohort-median Ktrans required to centre that effect is ",
    "not reported in the paper or supplements on disk; the effect is omitted ",
    "from model() and documented in the vignette. A Cox-style time-to-tumor ",
    "progression hazard h(t) = b0 * exp(b1 * dAUC24h_sVEGFR2) is described ",
    "in the paper but evaluated post-simulation in the vignette, not encoded ",
    "as an ODE."
  )
  reference <- paste(
    "Ait-Oudhia S, Mager DE, Pokuri V, Tomaszewski G, Groman A, Zagst P,",
    "Fetterly G, Iyer R.",
    "Bridging Sunitinib Exposure to Time-to-Tumor Progression in",
    "Hepatocellular Carcinoma Patients With Mathematical Modeling of an",
    "Angiogenic Biomarker.",
    "CPT Pharmacometrics Syst Pharmacol. 2016;5(6):297-304.",
    "doi:10.1002/psp4.12084.",
    sep = " "
  )
  vignette <- "Ait-Oudhia_2016_sunitinib"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    TUM_VOL = list(
      description        = "Baseline HCC tumor volume measured by DCE-MRI at the start of sunitinib treatment, used as the per-subject initial condition of the tumor compartment and to derive the growth rate constant kg = ln(2) / (114 * TG0^0.14) per Taouli 2005.",
      units              = "mm^3",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-subject time-fixed. The paper reports baseline tumor volumes for the n = 8 DCE-MRI subset; the validation vignette uses a representative typical baseline that reproduces Figure 2d. Use the subject's pre-treatment tumor volume on data ingestion; the column is consumed both as the d/dt(tumor) initial condition and inside the kg formula so a finite, positive value is required.",
      source_name        = "TG0"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 16L,
    n_studies      = 1L,
    age_range      = "adults with confirmed advanced HCC; demographics reported in Supplementary Table S1 (not on disk)",
    weight_range   = "not extractable from main text",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = "Advanced hepatocellular carcinoma (HCC) with 1-4 lesions, ECOG performance status 0/1/2, life expectancy >= 12 weeks. All patients received transarterial chemoembolisation (TACE) with doxorubicin 30 mg on Day 8 of Cycle 1.",
    dose_range     = "Sunitinib 37.5 mg PO QD, Days 1-7 + Days 15-35 of each 6-week cycle (4 weeks on / 2 weeks off, with a TACE-related break between Days 8 and 14). Repeat cycles until disease progression or unacceptable toxicity.",
    regions        = "Single-centre phase II pilot study at Roswell Park Cancer Institute (Buffalo, NY, USA).",
    notes          = "Single-arm open-label phase II pilot (n = 16 total; n = 8 had repeated DCE-MRI). PK / sVEGFR2 sampling at 24 h post-dose on Days 8, 10, and 35. Tumor volume measured by DCE-MRI on Days 0, 8, 10, and 35. Median observed TTP 7 months (8 months per Discussion). PK model parameters ka_D, Q_D, V2_D, ka_M, Q_M, V2_M and their IIVs (where reported as fixed) were inherited from the Houk 2009 sunitinib popPK meta-analysis (Ref 30 of the paper) under a MAP-Bayesian framework, and kout was inherited from the Lindauer 2010 sVEGFR2 popPD model (Ref 33). Detailed baseline demographics (age, weight, sex, race) are in Supplementary Table S1, which is not on disk; populate when the supplement becomes available."
  )

  ini({
    # ----------------------------------------------------------------------
    # Sunitinib (parent) PK -- Table 1.
    # The model is 2-compartment with first-order oral absorption. All
    # values in apparent units (V/F, CL/F). MAP-Bayesian re-fit of the
    # Houk 2009 priors against the HCC trough-concentration dataset.
    # ka_D, Q_D, V2_D are FIXED from Houk 2009 (Table 1 footnote "a").
    # ----------------------------------------------------------------------
    lvc <- log(1777);              label("Sunitinib apparent central volume V1_D/Fcentral (L)") # Table 1
    lcl <- log(30.3);               label("Sunitinib apparent clearance CL_D/Fcentral (L/h)")    # Table 1
    lka <- fixed(log(0.195));       label("Sunitinib first-order absorption rate ka_D (1/h, fixed Houk 2009)")    # Table 1 footnote a
    lq  <- fixed(log(6.37));        label("Sunitinib apparent inter-compartmental clearance Q_D/Fperipheral (L/h, fixed Houk 2009)")  # Table 1 footnote a
    lvp <- fixed(log(588));         label("Sunitinib apparent peripheral volume V2_D/Fperipheral (L, fixed Houk 2009)")                # Table 1 footnote a

    # ----------------------------------------------------------------------
    # SU12662 (metabolite) PK -- Table 1.
    # The metabolite is treated as receiving fM * Dose into its own depot
    # at each oral sunitinib dose (fM = 0.21 fixed). ka_M, Q_M, V2_M are
    # FIXED from Houk 2009.
    # ----------------------------------------------------------------------
    lvc_su12662 <- log(1840);       label("SU12662 apparent central volume V1_M/Fcentral (L)")   # Table 1
    lcl_su12662 <- log(19.72);      label("SU12662 apparent clearance CL_M/Fcentral (L/h)")      # Table 1
    lka_su12662 <- fixed(log(0.487)); label("SU12662 first-order absorption rate ka_M (1/h, fixed Houk 2009)")  # Table 1 footnote a
    lq_su12662  <- fixed(log(27.7));  label("SU12662 apparent inter-compartmental clearance Q_M/Fperipheral (L/h, fixed Houk 2009)")  # Table 1 footnote a
    lvp_su12662 <- fixed(log(345));   label("SU12662 apparent peripheral volume V2_M/Fperipheral (L, fixed Houk 2009)")                # Table 1 footnote a

    # ----------------------------------------------------------------------
    # sVEGFR2 indirect-response biomarker -- Table 2.
    # kout fixed from Lindauer 2010 (Table 2 footnote "a"); kin derived as
    # R0 * kout so the baseline reproduces R0 at steady state.
    # alpha is the dimensionless intrinsic activity scaling the saturable
    # inhibition term INH = ACub / (kd + ACub). The Table 2 row reports
    # alpha with units (ug/L)^-1, which is inconsistent with the equation
    # printed in the paper (Eq 7: kin / (1 + alpha * INH) requires alpha
    # dimensionless because INH is dimensionless). The value 0.77 is used
    # as-is per the paper text; the unit annotation is treated as a Table 2
    # typo (see vignette Assumptions and deviations).
    # ----------------------------------------------------------------------
    lr0_svegfr2    <- log(18.3);        label("sVEGFR2 baseline plasma concentration R0 (ug/L)") # Table 2
    lalpha_svegfr2 <- log(0.77);        label("sVEGFR2 intrinsic activity alpha (dimensionless)") # Table 2
    lkout_svegfr2  <- fixed(log(0.175)); label("sVEGFR2 first-order elimination rate kout (1/day, fixed Lindauer 2010)") # Table 2 footnote a

    # ----------------------------------------------------------------------
    # Tumor growth inhibition -- Table 2.
    # dIC50 is the dsVEGFR2 (= R0 - sVEGFR2) concentration producing 50% of
    # Imax (Imax = 1 fixed). The paper reports a significant Ktrans^beta
    # covariate effect on dIC50 (beta = 2.12, Table 2) but the cohort-
    # median Ktrans needed to centre that power covariate is not reported
    # in the paper or supplements on disk; the effect is omitted here and
    # documented in the vignette Assumptions and deviations section.
    # ----------------------------------------------------------------------
    ldic50 <- log(1.83); label("Tumor growth inhibition dIC50 (ug/L)") # Table 2

    # ----------------------------------------------------------------------
    # Inter-individual variability -- Table 1 and Table 2.
    # IIVs marked with footnote "a" were fixed from the corresponding
    # priors (Houk 2009 for PK, Lindauer 2010 for sVEGFR2 kout); the
    # remaining IIVs were estimated by the paper's MAP-Bayesian re-fit.
    # Variances are computed from the reported %CV via
    # omega^2 = log(CV^2 + 1).
    # ----------------------------------------------------------------------
    # Sunitinib PK IIVs -- all fixed from Houk 2009 (Table 1 footnote a)
    etalvc ~ fixed(0.1821); label("IIV variance on log V1_D/F (44.7% CV, fixed Houk 2009)")     # Table 1 footnote a
    etalcl ~ fixed(0.1342); label("IIV variance on log CL_D/F (37.9% CV, fixed Houk 2009)")     # Table 1 footnote a
    etalka ~ fixed(0.5066); label("IIV variance on log ka_D (81.2% CV, fixed Houk 2009)")       # Table 1 footnote a

    # SU12662 PK IIVs -- V1_M and CL_M BSVs fixed; ka_M BSV estimated
    # (Table 1 reports the "a" footnote on V1_M and CL_M BSV cells but
    # not on the ka_M BSV cell). Q_M and V2_M carry no IIV in Table 1.
    etalvc_su12662 ~ fixed(0.3534); label("IIV variance on log V1_M/F (65.1% CV, fixed Houk 2009)")  # Table 1 footnote a
    etalcl_su12662 ~ fixed(0.2410); label("IIV variance on log CL_M/F (52.2% CV, fixed Houk 2009)")  # Table 1 footnote a
    etalka_su12662 ~ 0.5843;        label("IIV variance on log ka_M (89.1% CV; estimated)")          # Table 1

    # sVEGFR2 IIVs -- kout BSV is reported with RSE on the BSV column
    # only (Table 2 footnote "a" is on the kout point estimate, not the
    # BSV). The BSV 67% CV is treated as estimated.
    etalr0_svegfr2    ~ 0.1118;  label("IIV variance on log sVEGFR2 baseline R0 (34.4% CV)") # Table 2
    etalalpha_svegfr2 ~ 0.04315; label("IIV variance on log sVEGFR2 alpha (21% CV)")          # Table 2
    etalkout_svegfr2  ~ 0.3708;  label("IIV variance on log sVEGFR2 kout (67% CV)")           # Table 2

    # Tumor growth IIV
    etaldic50 ~ 0.1218; label("IIV variance on log dIC50 (36% CV)") # Table 2

    # ----------------------------------------------------------------------
    # Residual error -- Table 1 (e_sunitinib, e_metabolite) and Table 2
    # (e_vegfR2 proportional; e_TG additive). The paper uses proportional
    # residuals for sunitinib, SU12662, and sVEGFR2 (Eq 13) and an
    # additive residual for tumor growth (Eq 14).
    # ----------------------------------------------------------------------
    propSd            <- 0.31;  label("Proportional residual error on sunitinib (fraction)") # Table 1
    propSd_su12662 <- 0.16;  label("Proportional residual error on SU12662 (fraction)")   # Table 1
    propSd_svegfr2    <- 0.26;  label("Proportional residual error on sVEGFR2 (fraction)")   # Table 2
    addSd_tumor       <- 15.5;  label("Additive residual error on tumor volume (mm^3)")       # Table 2
  })

  model({
    # ------------------------------------------------------------------
    # Fixed structural constants from the paper (not estimated).
    # ------------------------------------------------------------------
    fM   <- 0.21   # Fraction of sunitinib dose metabolised to SU12662 (Methods + Table 1 footnote)
    fb_D <- 0.9    # Sunitinib plasma protein bound fraction (Methods, Refs 15-17)
    fb_M <- 0.95   # SU12662 plasma protein bound fraction (Methods, Refs 15-17)
    kd   <- 4      # sVEGFR2 inhibition Kd in ug/L (= 4 ng/mL; Methods, Ref 31 Mendel 2003)
    Imax <- 1      # Maximum drug inhibitory effect on tumor growth (Methods + Table 2 footnote)
    tvdt_a <- 114  # Coefficient of the Taouli 2005 TVDT relationship: TVDT = 114 * TG0^0.14 (h)
    tvdt_b <- 0.14 # Exponent of the Taouli 2005 TVDT relationship

    # ------------------------------------------------------------------
    # Individual parameters -- typical * exp(eta).
    # ------------------------------------------------------------------
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl)
    ka <- exp(lka + etalka)
    q  <- exp(lq)
    vp <- exp(lvp)

    vc_su12662 <- exp(lvc_su12662 + etalvc_su12662)
    cl_su12662 <- exp(lcl_su12662 + etalcl_su12662)
    ka_su12662 <- exp(lka_su12662 + etalka_su12662)
    q_su12662  <- exp(lq_su12662)
    vp_su12662 <- exp(lvp_su12662)

    r0_svegfr2    <- exp(lr0_svegfr2    + etalr0_svegfr2)
    alpha_svegfr2 <- exp(lalpha_svegfr2 + etalalpha_svegfr2)
    kout_svegfr2  <- exp(lkout_svegfr2  + etalkout_svegfr2)  # 1/day in the paper
    kin_svegfr2   <- r0_svegfr2 * kout_svegfr2

    dic50 <- exp(ldic50 + etaldic50)

    # ------------------------------------------------------------------
    # PK micro-constants (all rate constants 1/h to match the units list).
    # The sVEGFR2 kout is reported in 1/day; convert to 1/h to share a
    # common time scale across the ODE system. The Taouli 2005 TVDT
    # constants are in hours.
    # ------------------------------------------------------------------
    kel        <- cl / vc
    k12        <- q  / vc
    k21        <- q  / vp
    kel_su12662 <- cl_su12662 / vc_su12662
    k12_su12662 <- q_su12662  / vc_su12662
    k21_su12662 <- q_su12662  / vp_su12662

    kout_svegfr2_h <- kout_svegfr2 / 24  # 1/day -> 1/h
    kin_svegfr2_h  <- kin_svegfr2  / 24  # ug/L/day -> ug/L/h, matches kout_svegfr2_h

    # ------------------------------------------------------------------
    # Plasma concentrations and active free unbound concentration.
    # AMT is in mg, apparent V in L, so `central/vc` is in mg/L. The
    # paper reports concentrations in ug/L throughout (R0 = 18.3 ug/L,
    # kd = 4 ug/L, dIC50 = 1.83 ug/L); scale Cc by 1000 (mg/L -> ug/L)
    # so the ACub equation and the sVEGFR2 indirect-response constants
    # share a single unit basis with the paper.
    # ------------------------------------------------------------------
    Cc          <- 1000 * central          / vc          # sunitinib plasma (ug/L)
    Cc_su12662  <- 1000 * central_su12662  / vc_su12662  # SU12662 plasma (ug/L)
    ACub        <- (1 - fb_D) * Cc + (1 - fb_M) * Cc_su12662   # active unbound concentration (ug/L)
    INH         <- ACub / (kd + ACub)

    # ------------------------------------------------------------------
    # Tumor growth rate from per-subject baseline tumor volume (Taouli
    # 2005). kg is computed at every integration step using TUM_VOL as
    # the fixed per-subject baseline; this matches the paper's
    # parameterisation kg = ln(2) / (114 * TG0^0.14) where TG0 is the
    # observed baseline tumor volume in mm^3.
    # ------------------------------------------------------------------
    tg0_safe <- TUM_VOL
    kg <- log(2) / (tvdt_a * tg0_safe^tvdt_b)

    # ------------------------------------------------------------------
    # Inhibition function H(t) (Methods + Eq 8).
    # ------------------------------------------------------------------
    dsVEGFR2 <- r0_svegfr2 - svegfr2
    Hfun     <- Imax * dsVEGFR2 / (dsVEGFR2 + dic50)

    # ------------------------------------------------------------------
    # ODEs.
    # The metabolite depot receives fM * Dose at every sunitinib dose
    # event via the dose-multiplier f(depot_su12662) <- fM applied to a
    # paired dose record on depot_su12662 in the event table.
    # ------------------------------------------------------------------
    d/dt(depot)                <- -ka * depot
    d/dt(central)              <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1)          <-  k12 * central - k21 * peripheral1

    d/dt(depot_su12662)        <- -ka_su12662 * depot_su12662
    d/dt(central_su12662)      <-  ka_su12662 * depot_su12662 - kel_su12662 * central_su12662 - k12_su12662 * central_su12662 + k21_su12662 * peripheral1_su12662
    d/dt(peripheral1_su12662)  <-  k12_su12662 * central_su12662 - k21_su12662 * peripheral1_su12662

    d/dt(svegfr2)              <-  kin_svegfr2_h / (1 + alpha_svegfr2 * INH) - kout_svegfr2_h * svegfr2
    d/dt(tumor)                <-  kg * (1 - Hfun) * tumor

    # Bioavailability of the metabolite depot is the metabolised
    # fraction fM. Sunitinib depot inherits the unmodified dose.
    f(depot_su12662) <- fM

    # ------------------------------------------------------------------
    # Initial conditions.
    # ------------------------------------------------------------------
    svegfr2(0) <- r0_svegfr2
    tumor(0)   <- TUM_VOL

    # ------------------------------------------------------------------
    # Observations and residual error (Eq 13 + Eq 14).
    # ------------------------------------------------------------------
    Cc         ~ prop(propSd)
    Cc_su12662 ~ prop(propSd_su12662)
    svegfr2    ~ prop(propSd_svegfr2)
    tumor      ~ add(addSd_tumor)
  })
}
