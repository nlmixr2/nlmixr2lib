Gibbs_2017_evolocumab <- function() {
  description <- "Target-mediated drug disposition (quasi-steady-state approximation) PK/PD model for evolocumab, an anti-PCSK9 monoclonal antibody, in healthy subjects and statin-treated patients with hypercholesterolemia (Gibbs 2017): one-compartment evolocumab disposition with first-order SC absorption, linear clearance plus PCSK9-complex-mediated (target) elimination, coupled to a PCSK9 turnover pool and to a type-4 indirect-response model in which unbound PCSK9 inhibits LDL-cholesterol elimination."
  reference <- "Gibbs JP, Doshi S, Kuchimanchi M, Grover A, Emery MG, Dodds MG, Gibbs MA, Somaratne R, Wasserman SM, Blom D. Impact of target-mediated elimination on the dose and regimen of evolocumab, a human monoclonal antibody against proprotein convertase subtilisin/kexin type 9 (PCSK9). J Clin Pharmacol. 2017;57(5):616-626. doi:10.1002/jcph.840"
  vignette <- "Gibbs_2017_evolocumab"

  # Time in days; doses entered in mg and converted to nmol on the dose record
  # (see f(depot) in model()); all drug and target concentrations are molar (nM),
  # matching the paper's equations 1-7 and Tables 2-3. LDL-C is in mg/dL.
  units <- list(time = "day", dosing = "mg", concentration = "nM")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE: each entry was checked against Gibbs
  # 2017 equations 1-4 and 7 and the Figure 1 schematic.
  #
  # NOTE the mixed amount/concentration bases, which are the paper's own:
  # equations 1-2 carry evolocumab as an AMOUNT (Adepot, TDA) while equations 3
  # and 7 carry PCSK9 and LDL-C as CONCENTRATIONS (TLC, LDL). TDA is converted
  # to the concentration TDC by TDA/V (equation 4).
  compartmentData <- list(
    depot        = list(analyte = "evolocumab", units = "nmol", specimen = "administration site", verified = TRUE),
    central      = list(analyte = "total (unbound + PCSK9-bound) evolocumab", units = "nmol", specimen = "serum", verified = TRUE),
    total_target = list(analyte = "total (unbound + evolocumab-bound) PCSK9", units = "nM", specimen = "serum", verified = TRUE),
    ldl          = list(analyte = "low-density lipoprotein cholesterol", units = "mg/dL", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    DIS_HEALTHY = list(
      description        = "Healthy-participant cohort indicator, 1 = healthy subject (phase 1a study 20080397), 0 = hypercholesterolemic patient on stable statin therapy (phase 1b study 20080398)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (statin-treated hypercholesterolemic patient; the population in which BASE_PCSK9 = 5.27 nM was estimated)",
      notes              = "The only covariate retained in the Gibbs 2017 final model. Multiplicative effect on the baseline PCSK9 turnover set point: BASE_PCSK9 * theta1^DIS_HEALTHY with theta1 = 0.637 (Table 2). The paper reports the two resulting typical values explicitly in the Results, which over-determines the coding and fixes the direction: 5.27 nM (379 ng/mL) in statin-treated patients and 5.27 * 0.637 = 3.36 nM (242 ng/mL) in healthy subjects. Because BASE_PCSK9 also sets ksyn (= kdeg * BASE_PCSK9) and the LDL-C baseline balance (kin), this single covariate propagates into both the target-mediated clearance of evolocumab and the achievable LDL-C reduction. Time-fixed per subject.",
      source_name        = "theta 1 (Gibbs 2017 Table 2 footnote: 'fold change in baseline PCSK9 for healthy subjects vs statin-treated patients')"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 101L,
    n_observations   = 4910L,
    n_studies        = 2L,
    n_healthy        = 44L,
    n_patients       = 57L,
    age_mean         = "45.5 years",
    weight_mean      = "81.4 kg",
    sex_female_pct   = 26.7,
    race_ethnicity   = "83.2% white (pooled across both studies; the paper reports only the pooled white percentage)",
    disease_state    = "Pooled: healthy subjects (phase 1a study 20080397, single ascending SC dose) and patients with hypercholesterolemia on stable low- to moderate-intensity statin therapy (phase 1b study 20080398, multiple ascending SC dose). Two additional phase 1b cohorts received high-intensity statin therapy or carried a heterozygous familial hypercholesterolemia diagnosis and were dosed 140 mg Q2W x 3.",
    dose_range       = "Phase 1a: single SC 7, 21, 70, 210, or 420 mg. Phase 1b: SC 14 mg QW x 6, 35 mg QW x 6, 140 mg Q2W x 3, 280 mg Q2W x 3, or 420 mg QM x 2.",
    pcsk9_baseline   = "5.27 nM (379 ng/mL) typical in statin-treated patients; 3.36 nM (242 ng/mL) typical in healthy subjects. Baseline PCSK9 was 1.58-fold higher in the phase 1b study than the phase 1a study.",
    ldlc_baseline    = "116 mg/dL typical (BASE_LDL-C, Table 3)",
    notes            = "A total of 73 participants received evolocumab and 28 received placebo across the two studies; the model dataset comprised the 101 individuals with evaluable PK/PD data. Baseline demographics are given in the paper's Table S1 (Supporting Information), which is not in the on-disk source set; the demographic values recorded here are the pooled summaries stated in the Results narrative. Estimation used NONMEM 7.2 with SAEM followed by importance sampling, and the M3 method for below-limit-of-quantification unbound evolocumab and unbound PCSK9 observations."
  )

  ini({
    # ================= Evolocumab PK (Gibbs 2017 Table 2) =================
    # One-compartment disposition, first-order SC absorption. Reference
    # subject for BASE_PCSK9 is the statin-treated patient (DIS_HEALTHY = 0).
    lka     <- log(0.245);        label("First-order SC absorption rate constant ka (1/day)")                  # Table 2: ka = 0.245 1/day (RSE 12%)
    lcl     <- log(0.256);        label("Linear (non-target-mediated) clearance CL (L/day)")                    # Table 2: CL = 0.256 L/day (RSE 17%)
    lvc     <- log(2.66);         label("Central volume of distribution V (L)")                                 # Table 2: V = 2.66 L (RSE 5.9%)
    lfdepot <- fixed(log(0.72));  label("Subcutaneous bioavailability F (fraction)")                            # Methods: "Bioavailability (F) was fixed to 72% (Amgen internal data)"

    # ================= PCSK9 target turnover and binding (Table 2) =================
    lkdeg         <- log(2.12);   label("First-order PCSK9 degradation rate constant kdeg (1/day)")             # Table 2: kdeg = 2.12 1/day (RSE 7.3%); Discussion: "elimination half-life of approximately 8 hours"
    lrbase_target <- log(5.27);   label("Baseline total PCSK9 concentration BASE_PCSK9 in statin-treated patients (nM)") # Table 2: BASE_PCSK9 = 5.27 nM (RSE 2.7%); Results: "5.27 nM (or 379 ng/mL) in hypercholesterolemic patients stably treated with statins"
    lkss          <- log(0.253);  label("Quasi-steady-state constant kss = (kint + koff)/kon (nM)")             # Table 2: kss = 0.253 nM (RSE 23%)
    lkint         <- log(0.0529); label("Evolocumab-PCSK9 complex elimination (internalisation) rate constant kint (1/day)") # Table 2: kint = 0.0529 1/day (RSE 3.5%)

    # ---- Covariate effect (Table 2, theta 1) ----
    e_dis_healthy_rbase_target <- 0.637; label("Fold change in baseline PCSK9 for healthy subjects vs statin-treated patients (unitless)") # Table 2: theta1 = 0.637 (RSE 4.1%); Results: 5.27 * 0.637 = 3.36 nM (242 ng/mL) in healthy subjects

    # ================= LDL-C indirect response (Gibbs 2017 Table 3) =================
    lkout      <- log(0.305);      label("Maximal first-order LDL-C elimination rate constant kout, i.e. in the absence of PCSK9 (1/day)") # Table 3: kout = 0.305 1/day (RSE 4.0%); Discussion: "half-life of approximately 2.3 days ... in the absence of PCSK9"
    lrbase_ldl <- log(116);        label("Baseline LDL-C concentration BASE_LDL-C (mg/dL)")                      # Table 3: BASE_LDL-C = 116 mg/dL (RSE 2.9%)
    limax      <- fixed(log(1.00)); label("Maximal proportional inhibition of kout by PCSK9, Imax (unitless)")   # Table 3: Imax = 1.00 (fix); Results: "Initial attempts to estimate Imax suggested that it was near the boundary of 1.00, so its value was fixed to 1.00"
    lic50      <- log(1.46);       label("Unbound PCSK9 concentration giving half-maximal inhibition of kout, IC50 (nM)") # Table 3: IC50 = 1.46 nM (RSE 10%); Abstract: "IC50 ... was 1.46 nM"

    # ================= Between-subject variability =================
    # Gibbs 2017 Tables 2-3 report the random-effects column as VARIANCES: the
    # table footnote defines "CV%, coefficient of variation calculated as
    # sqrt(var) * 100", and every printed CV% reproduces as sqrt(estimate):
    #   sqrt(0.807)  = 0.898 -> 90%  (ka)
    #   sqrt(0.158)  = 0.397 -> 40%  (V)
    #   sqrt(0.924)  = 0.961 -> 96%  (CL)
    #   sqrt(0.0348) = 0.187 -> 19%  (BASE_PCSK9)
    #   sqrt(1.17)   = 1.082 -> 108% (kss)
    #   sqrt(0.0465) = 0.216 -> 22%  (BASE_LDL-C)
    #   sqrt(0.481)  = 0.694 -> 69%  (IC50)
    # so the estimates are carried through unchanged (no CV% -> omega^2
    # back-transform is needed or appropriate here).
    #
    # Full 3x3 block on ka, V, CL (Table 2 "Covariance" rows). The parenthesised
    # percentages on those rows are correlations, and they reproduce from the
    # covariances, confirming the block is a covariance (not correlation) matrix:
    #   0.253 / sqrt(0.807 * 0.158) =  0.709 ->  70%
    #  -0.671 / sqrt(0.807 * 0.924) = -0.777 -> -78%
    #  -0.183 / sqrt(0.158 * 0.924) = -0.479 -> -48%
    # Lower-triangle order below is var(ka); cov(ka,V), var(V); cov(ka,CL), cov(V,CL), var(CL).
    etalka + etalvc + etalcl ~ c(0.807,
                                 0.253, 0.158,
                                -0.671, -0.183, 0.924)   # Table 2: variances 0.807 / 0.158 / 0.924 and covariances 0.253 / -0.671 / -0.183
    etalrbase_target ~ 0.0348   # Table 2: BASE_PCSK9 intersubject variance 0.0348 (~19% CV)
    etalkss          ~ 1.17     # Table 2: kss intersubject variance 1.17 (~108% CV)
    etalrbase_ldl    ~ 0.0465   # Table 3: BASE_LDL-C intersubject variance 0.0465 (~22% CV)
    etalic50         ~ 0.481    # Table 3: IC50 intersubject variance 0.481 (~69% CV)
    # kdeg, kint, theta1, kout and Imax carry no between-subject variability:
    # Tables 2-3 report 0.00 with footnote a, "Intersubject random variance was
    # fixed at 0", and Methods states "BSV for kdeg and kint were fixed at 0%".

    # ================= Residual unexplained variability =================
    # Tables 2-3 report sigma^2 (variances); nlmixr2 expects standard deviations,
    # so each value below is sqrt(sigma^2), and each reproduces the paper's own
    # parenthesised summary.
    propSd            <- 0.240;  label("Proportional residual error on unbound evolocumab (fraction)")  # Table 2: sigma^2 (evolocumab) = 0.0576 -> sqrt = 0.240; Results: "residual variabilities for unbound evolocumab ... were 24%"
    propSd_freePcsk9  <- 0.3069; label("Proportional residual error on unbound PCSK9 (fraction)")        # Table 2: sigma^2 (PCSK9) = 0.0942 -> sqrt = 0.30692; Results: "... and unbound PCSK9 were 24% and 31%"
    propSd_ldl        <- 0.1140; label("Proportional residual error on LDL-C (fraction)")                # Table 3: sigma^2 (proportional error) = 0.0130 -> sqrt = 0.11402; Results: "proportional and additive residual variability of 11% and 7.8 mg/dL"
    addSd_ldl         <- 7.7524; label("Additive residual error on LDL-C (mg/dL)")                       # Table 3: sigma^2 (additive error) = 60.1 -> sqrt = 7.7524, and the table's own parenthesised 7.8 is flagged by footnote b as "the standard deviation ... calculated as sqrt(var)"
  })

  model({
    # ---- Physical constants ----------------------------------------------
    # Evolocumab molecular weight. NOT reported in Gibbs 2017: the paper works
    # entirely in molar units and never prints a mg <-> nmol bridge, but doses
    # are administered in mg, so one is required. 141,800 g/mol is the value in
    # the FDA-approved Repatha (evolocumab) prescribing information, and is the
    # same constant already used by the sibling registry model
    # inst/modeldb/specificDrugs/Kuchimanchi_2018_evolocumab.R (MW_EVO).
    # Non-paper provenance; see the vignette's Assumptions and deviations
    # section, which also shows that the paper's reported simulation endpoints
    # move by < 1 percentage point when this constant is varied to 150,000.
    MW_EVO      <- 141800        # g/mol
    nmol_per_mg <- 1e6 / MW_EVO  # 7.052 nmol per mg

    # ---- Individual PK parameters (Table 2) -------------------------------
    ka     <- exp(lka + etalka)
    cl     <- exp(lcl + etalcl)
    vc     <- exp(lvc + etalvc)
    fdepot <- exp(lfdepot)
    kdeg   <- exp(lkdeg)
    kint   <- exp(lkint)
    kss    <- exp(lkss + etalkss)

    # Baseline PCSK9 set point. theta1 multiplies BASE_PCSK9 for healthy
    # subjects; the statin-treated patient (DIS_HEALTHY = 0) is the reference.
    rbase_target <- exp(lrbase_target + etalrbase_target) *
                    e_dis_healthy_rbase_target^DIS_HEALTHY

    # Zero-order PCSK9 production. Methods: ksyn and kdeg are "related to the
    # baseline PCSK9 (BASE_PCSK9) level as ksyn/kdeg", i.e. ksyn = kdeg * BASE.
    ksyn <- kdeg * rbase_target

    # ---- Individual PD parameters (Table 3) -------------------------------
    kout      <- exp(lkout)
    imax      <- exp(limax)
    ic50      <- exp(lic50 + etalic50)
    rbase_ldl <- exp(lrbase_ldl + etalrbase_ldl)

    # Zero-order LDL-C production, set so that the LDL-C state holds at
    # BASE_LDL-C in the drug-free state, where unbound PCSK9 = BASE_PCSK9:
    #   0 = kin - kout * (1 - Imax * BASE_PCSK9/(IC50 + BASE_PCSK9)) * BASE_LDL-C
    kin <- kout * (1 - imax * rbase_target / (ic50 + rbase_target)) * rbase_ldl

    # ---- Quasi-steady-state TMDD algebra (Gibbs 2017 equations 4-6) -------
    # tdc: total drug concentration = TDA / V                     (equation 4)
    # fdc: free (unbound) drug concentration                      (equation 5)
    # flc: free (unbound) ligand [PCSK9] concentration            (equation 6)
    # Equation 5 is the positive root of the one-to-one QSS binding quadratic
    #   FDC^2 + FDC * (TLC - TDC + kss) - kss * TDC = 0,
    # so FDC = 0.5 * [(TDC - TLC - kss) + sqrt((TDC - TLC - kss)^2 + 4*kss*TDC)].
    # The drug-target complex concentration is TDC - FDC, hence equation 6.
    tdc <- central / vc
    tlc <- total_target
    qss <- tdc - tlc - kss
    fdc <- 0.5 * (qss + sqrt(qss * qss + 4 * kss * tdc))
    flc <- tlc - (tdc - fdc)

    # ---- ODE system (Gibbs 2017 equations 1-3 and 7, Figure 1) ------------
    # Equation 2 is written by the authors as "- k * FDC * V - kint * TLC * FDC * V/(kss + FDC)";
    # k is the first-order elimination rate constant kel = CL/V, so k * FDC * V = CL * FDC.
    d/dt(depot)        <- -ka * depot
    d/dt(central)      <-  ka * depot -
                           cl * fdc -
                           kint * tlc * fdc * vc / (kss + fdc)
    d/dt(total_target) <-  ksyn - kdeg * tlc -
                           (kint - kdeg) * fdc * tlc / (kss + fdc)
    d/dt(ldl)          <-  kin - kout * (1 - imax * flc / (ic50 + flc)) * ldl

    # Initial conditions: both turnover pools start at their drug-free baselines.
    total_target(0) <- rbase_target
    ldl(0)          <- rbase_ldl

    # SC bioavailability, combined with the mg -> nmol conversion so that dose
    # records are entered in mg while the compartments are molar.
    f(depot) <- fdepot * nmol_per_mg

    # ---- Observations and residual error ---------------------------------
    # The three measured analytes are unbound evolocumab, unbound PCSK9 (both
    # nM; the assays deliberately do not detect the drug-target complex, see
    # Methods "PK/PD Sampling and Assays"), and LDL-C (mg/dL).
    Cc         <- fdc
    freePcsk9  <- flc

    Cc        ~ prop(propSd)
    freePcsk9 ~ prop(propSd_freePcsk9)
    ldl       ~ add(addSd_ldl) + prop(propSd_ldl)
  })
}
