Bouhaddou_2020_iadademstat_mouse <- function() {
  description <- "QSP. Preclinical (mouse, athymic nude, NCI-H510A small-cell lung cancer xenograft). Semimechanistic in vivo PK/PD model of the covalent LSD1 (KDM1A) inhibitor iadademstat (ORY-1001) from Bouhaddou 2020. An oral two-compartment PK model with dose-dependent (saturable) bioavailability F = D / (D50 + D) drives the unbound plasma concentration (fu = 0.75), assumed equal to the free intracellular tumour concentration. That concentration irreversibly inactivates LSD1 (kinact / Ki), the drug-bound LSD1 is degraded by Michaelis-Menten kinetics, percent target engagement (TE) suppresses GRP mRNA through a sigmoidal inhibition with a self-regulating production rate, and the GRP level switches tumour cells between a proliferating and a quiescent population (reversible cytostasis). Every PD parameter is carried over unchanged from the in vitro NCI-H510A model (modellib('Bouhaddou_2020_iadademstat_invitro')); only the intrinsic growth constant kP was re-estimated on drug-free xenograft growth. The output tumor_vol (P + Q, mm^3) is the Figure 5 endpoint. Deterministic: the parameters were fitted by least squares (MATLAB lsqnonlin / fmincon multistart), so the model carries no between-animal variability and no residual error."
  reference <- paste(
    "Bouhaddou M, Yu LJ, Lunardi S, Stamatelos SK, Mack F, Gallo JM,",
    "Birtwistle MR, Walz AC. (2020). Predicting In Vivo Efficacy from In Vitro",
    "Data: Quantitative Systems Pharmacology Modeling for an Epigenetic Modifier",
    "Drug in Cancer. Clin Transl Sci 13(2):419-429. doi:10.1111/cts.12727.",
    "Parameter values from Supplementary Table S2 and the deposited MATLAB code",
    "(Supplementary zip, vivo_model/rateconstants_vivo_BEST.txt, RunModelVivo.m)."
  )
  vignette <- "Bouhaddou_2020_iadademstat"

  # The four PD states are paper-mechanistic with no canonical analogue in
  # inst/references/compartment-names.md: drug-bound (inactivated) LSD1, the
  # GRP mRNA biomarker, and the proliferating / quiescent tumour-cell pools of
  # the epigenetic switch (named as in Mazzocco_2015_temozolomide.R).
  paper_specific_compartments <- c("lsd1b", "grp", "prolif", "quiesc")

  dosing <- "depot"

  units <- list(
    time = "h",
    dosing = "ng",
    concentration = "ng/mL for Cc; nM for the unbound tumour concentration roc and the LSD1 states; mm^3 for tumor_vol"
  )

  # Issue #482: what each ODE state holds, in what units, in what matrix.
  # verified = TRUE -- checked against Bouhaddou 2020 Eqs. 1-20, Table S1 and
  # the deposited RunModelVivo.m.
  compartmentData <- list(
    depot = list(
      analyte = "iadademstat",
      units = "ng",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "iadademstat",
      units = "ng",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "iadademstat",
      units = "ng",
      specimen = "tissue",
      verified = TRUE
    ),
    lsd1b = list(
      analyte = "LSD1 (KDM1A) covalently bound to iadademstat",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    grp = list(
      analyte = "gastrin-releasing peptide (GRP) mRNA, relative to the untreated level",
      units = "unitless (1 = untreated)",
      specimen = "tumor",
      verified = TRUE
    ),
    prolif = list(
      analyte = "proliferating NCI-H510A tumour cells",
      units = "mm^3",
      specimen = "not applicable",
      verified = TRUE
    ),
    quiesc = list(
      analyte = "quiescent NCI-H510A tumour cells",
      units = "mm^3",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  # No subject-level covariates: a single typical mouse. The dose is entered
  # as an absolute amount (ng) because the saturable bioavailability is a
  # function of that amount; the authors converted ug/kg with a fixed
  # 0.025 kg mouse weight (RunDosingRegVivo.m: dose = dose_ugkg * 1E3 * 0.025).
  covariateData <- list()

  population <- list(
    species = "mouse (athymic nude, female, NCI-H510A xenograft; PK in tumour-free mice)",
    n_subjects = "PK: 3 mice per time point at each of 3 single oral doses (20, 40 and 10000 ug/kg; 42 plasma samples); efficacy: 13-15 xenograft-bearing mice per arm in 5 arms",
    n_studies = 2L,
    age_range = "7-8 weeks at subcutaneous implantation of 5 million NCI-H510A cells",
    weight_range = "0.025 kg assumed by the authors for the ug/kg to ng dose conversion",
    sex_female_pct = 100,
    disease_state = "Subcutaneous NCI-H510A small-cell lung cancer (SCLC) xenograft",
    dose_range = "Oral iadademstat 20 or 40 ug/kg once daily 5-on/2-off, 200 ug/kg weekly, or 400 ug/kg on days 7, 16 and 23 (Figure 5); single oral 20, 40 and 10000 ug/kg for PK (Figure 4)",
    regions = "Roche Innovation Center New York (USA) and Oryzon Genomics (Spain)",
    notes = "The PD parameters come from in vitro NCI-H510A experiments (target engagement, GRP mRNA, cell growth and viability); only the drug-free growth constant kP was fitted to xenograft data. Fraction unbound 0.75 was measured in mouse plasma."
  )

  ini({
    # ==================================================================
    # Oral two-compartment PK -- Bouhaddou 2020 Eqs. 1-5, values from
    # Table S2 ('Value in vivo') and rateconstants_vivo_BEST.txt.
    # ==================================================================
    lka <- log(0.88205)
    label("First-order absorption rate constant ka (1/h)") # Table S2 ka = 0.88205
    lcl <- log(117.66)
    label("Clearance CL (mL/h)") # Table S2 'kCL' = 117.66; used as CL = CLml / V in RunModelVivo.m
    lvc <- log(577.49)
    label("Central volume of distribution V (mL)") # Table S2 V = 577.49 mL
    lk12 <- log(0.47004)
    label("Central-to-peripheral rate constant k12 (1/h)") # Table S2 k12 = 0.47004
    lk21 <- log(0.25469)
    label("Peripheral-to-central rate constant k21 (1/h)") # Table S2 k21 = 0.25469
    ld50 <- log(1123.9)
    label("Dose giving half-maximal bioavailability D50 (ng)") # Table S2 D50 = 1123.9 (a.u.; ng per RunModelVivo.m)
    fu <- fixed(0.75)
    label("Fraction unbound in mouse plasma fu (unitless), measured") # Table S2 fu = 0.75 (Experiment); Results p. 424
    mw <- fixed(230.4)
    label("Molecular weight of iadademstat (g/mol), ng/mL to nM conversion") # rateconstants_vivo_BEST.txt RO_molweight = 230.4

    # ==================================================================
    # Target engagement -- Eqs. 6-9, Table S2 (in vitro values, 'Same' in
    # vivo).
    # ==================================================================
    ki <- fixed(19)
    label("Inhibitory constant Ki (nM), measured") # Table S2 Ki = 19 (Experiment)
    lkinact <- log(0.87052)
    label("Maximal inactivation rate constant kinact (1/h)") # Table S2 kinact = 0.87052
    lvmax_lsd1b <- log(0.03515)
    label("Maximal degradation rate of drug-bound LSD1 Vmax (nM/h)") # Table S2 vmax = 0.03515
    lkm_lsd1b <- log(0.9028)
    label("Michaelis-Menten constant of bound-LSD1 degradation Km (nM)") # Table S2 km = 0.9028
    llsd1tot <- log(3.0833)
    label("Total LSD1 in the cell LSD1TOTAL (nM)") # Table S2 LSD1TOTAL = 3.0833

    # ==================================================================
    # GRP mRNA biomarker -- Eqs. 10-12, Table S2.
    # ==================================================================
    lte50 <- log(26.513)
    label("Target engagement giving half-maximal GRP suppression K50 (%)") # Table S2 K50 = 26.513
    hill_te <- 2.0819
    label("Hill coefficient of TE on GRP production n (unitless)") # Table S2 n = 2.0819
    lkdeg_grp <- log(0.025567)
    label("First-order GRP mRNA degradation rate constant kdeg (1/h)") # Table S2 kdeg = 0.025567
    b_grp <- 0.083048
    label("Intercept b of the GRP-dependent production rate kmax = -m*GRP + b (1/h)") # Table S2 b = 0.083048

    # ==================================================================
    # Tumour-cell epigenetic switch -- Eqs. 13-20, Table S2. kP is the only
    # parameter re-estimated in vivo.
    # ==================================================================
    lkp <- log(0.0037)
    label("Intrinsic growth rate constant of proliferating cells kP (1/h)") # Table S2 kP in vivo = 0.0037
    lk50p <- log(100000)
    label("Proliferating-cell level giving half-maximal growth k50P (mm^3)") # Table S2 k50P = 100000 ('Same'); rateconstants_vivo_BEST.txt
    lkmaxpq <- log(0.5932)
    label("Maximal proliferating-to-quiescent switching rate kmaxPQ (1/h)") # Table S2 kmaxPQ = 0.5932
    lkmaxqp <- log(4.4908)
    label("Maximal quiescent-to-proliferating switching rate kmaxQP (1/h)") # Table S2 kmaxQP = 4.4908
    k50pq <- 0.8836
    label("Biomarker level giving half-maximal P-to-Q switching k50PQ (unitless)") # Table S2 k50PQ = 0.8836
    k50qp <- 0.99897
    label("Biomarker level giving half-maximal Q-to-P switching k50QP (unitless)") # Table S2 k50QP = 0.99897
    hill_pq <- 3.8313
    label("Hill coefficient of the P-to-Q switch nPQ (unitless)") # Table S2 nPQ = 3.8313
    hill_qp <- 37.455
    label("Hill coefficient of the Q-to-P switch nQP (unitless)") # Table S2 nQP = 37.455
    bs <- 0.033328
    label("Basal bias term of Q-to-P switching bs (unitless)") # Table S2 bs = 0.033328

    prolif0 <- fixed(70)
    label("Initial proliferating tumour volume P(0) (mm^3)") # Table S1 P in vivo = 70 mm^3
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    kel <- cl / vc
    d50 <- exp(ld50)

    # Eq. 1 prints dX/dt = -ka * D/(D50 + D) * D, which does not balance
    # against Eq. 2 (ka * X). The deposited RunModelVivo.m resolves it: the
    # depot is first-order (dX = -ka*X) and the administered amount D is
    # scaled once, at dosing, by the fraction D/(D50 + D) -- i.e. a
    # dose-dependent bioavailability. podo(depot) is that dose amount (ng).
    fdepot <- podo(depot) / (d50 + podo(depot))

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central +
      k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    f(depot) <- fdepot

    # Eq. 5: plasma concentration (ng/mL).
    Cc <- central / vc

    # Eq. 4: pharmacologically active (unbound, intracellular) concentration,
    # converted to nM with the molecular weight exactly as RunModelVivo.m
    # does (fac_ng2nM = (1000 / MW) / V * fu). The code integrates this as an
    # ODE equal to fac_ng2nM * dQc/dt from zero, which is identical to the
    # algebraic form here.
    roc <- Cc * fu * 1000 / mw

    # Target engagement, Eqs. 7-9.
    kinact <- exp(lkinact)
    vmax_lsd1b <- exp(lvmax_lsd1b)
    km_lsd1b <- exp(lkm_lsd1b)
    lsd1tot <- exp(llsd1tot)
    lsd1u <- lsd1tot - lsd1b
    d/dt(lsd1b) <- kinact * roc / (ki + roc) * lsd1u -
      vmax_lsd1b / (km_lsd1b + lsd1b) * lsd1b
    TE <- lsd1b / lsd1tot * 100

    # GRP mRNA biomarker, Eqs. 10-12. With TE = 0 the steady state is
    # GRP = 1 (the untreated level, Table S1).
    te50 <- exp(lte50)
    kdeg_grp <- exp(lkdeg_grp)
    m_grp <- b_grp - kdeg_grp
    kmax_grp <- -m_grp * grp + b_grp
    d/dt(grp) <- kmax_grp * te50^hill_te / (te50^hill_te + TE^hill_te) -
      kdeg_grp * grp
    grp(0) <- 1

    # Epigenetic switch, Eqs. 13-20. BM = 1 - GRP (Eq. 18). GRP cannot rise
    # above 1 in this model, so BM >= 0; abs() (as in the authors' in vitro
    # RunModelVitro.m) only guards against a round-off BM of -1e-17 being
    # raised to the non-integer power nPQ.
    kp <- exp(lkp)
    k50p <- exp(lk50p)
    kmaxpq <- exp(lkmaxpq)
    kmaxqp <- exp(lkmaxqp)
    bm_grp <- abs(1 - grp)
    bmpq <- bm_grp
    bmqp <- 1 - bm_grp
    vp <- kp * (1 - prolif / (k50p + prolif)) * prolif
    vpq <- kmaxpq * bmpq^hill_pq / (k50pq^hill_pq + bmpq^hill_pq) * prolif
    vqp <- kmaxqp * (bmqp^hill_qp / (k50qp^hill_qp + bmqp^hill_qp) + bs) * quiesc
    d/dt(prolif) <- vp - vpq + vqp
    d/dt(quiesc) <- vpq - vqp
    prolif(0) <- prolif0

    # Tumour volume (mm^3), the Figure 5b endpoint (RunDosingRegVivo.m:
    # sum(yout(:,4:5), 2)).
    tumor_vol <- prolif + quiesc
  })
}
