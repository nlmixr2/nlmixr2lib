Valenzuela_2011_oxaliplatin <- function() {
  description <- "Population PK/PD model for hyperthermic intraperitoneal oxaliplatin (HIO) and induced neutropenia in 30 adults with peritoneal carcinomatosis after cytoreductive surgery (Valenzuela 2011). PK: peritoneum-as-depot first-order absorption (parameterized in the paper as peritoneum-to-plasma clearance Qa and peritoneum volume Va = vd, with ka = Qa/Va as a secondary parameter) feeding an open two-compartment plasma disposition; bioavailability F was fixed to 1 so Cl/F, Vc/F, Q2/F, Vp/F are apparent. PD: Friberg semi-mechanistic myelosuppression chain (one proliferating compartment plus three transit compartments feeding circulating ANC) with a linear drug effect Edrug = alpha * Cc on the proliferation rate and a (Circ0/Circ)^gamma feedback amplification; MTT was fixed at 118 h and the circulating-cell elimination rate constant kCirc was fixed at 0.07 per h (both from Friberg 2002). No subject covariates were retained in the final model; ten demographic and biochemistry covariates were screened graphically and showed no correlation with PK/PD parameters."
  reference <- paste(
    "Valenzuela B, Nalda-Molina R, Bretcha-Boix P, Escudero-Ortiz V,",
    "Duart MJ, Carbonell V, Sureda M, Rebollo JP, Farre J, Brugarolas A,",
    "Perez-Ruixo JJ. Pharmacokinetic and Pharmacodynamic Analysis of",
    "Hyperthermic Intraperitoneal Oxaliplatin-Induced Neutropenia in",
    "Subjects with Peritoneal Carcinomatosis. AAPS J. 2011 Mar;13(1):72-82.",
    "doi:10.1208/s12248-010-9249-2 (PMID 21210260).",
    "Erratum: AAPS J. 2011 Jun;13(2):318 (referenced in the PubMed",
    "'Erratum in' field of PMID 21210260; the erratum text was not",
    "available on disk at extraction time, and a public-web search",
    "for its full content was inconclusive -- operator confirmation",
    "is recommended; see vignette Errata).",
    "Friberg semi-mechanistic myelosuppression backbone:",
    "Friberg LE, Henningsson A, Maas H, Nguyen L, Karlsson MO.",
    "Model of chemotherapy-induced myelosuppression with parameter",
    "consistency across drugs. J Clin Oncol 2002;20(24):4713-4721.",
    "doi:10.1200/JCO.2002.02.140 (PMID 12488418; reference 30 in",
    "Valenzuela 2011, source of fixed MTT = 118 h and fixed kCirc = 0.07/h).",
    sep = " "
  )
  vignette <- "Valenzuela_2011_oxaliplatin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L",
                anc = "10^9 cells/L")

  covariateData <- list()

  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Age",
      units              = "year",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened graphically in Valenzuela 2011 (Results, Pharmacodynamics paragraph: 'The exploratory graphical analysis of the effect of age, sex, body weight, serum creatinine, albumin, serum ALT, serum AST, total bilirubin, hemoglobin, and hematocrit did not suggest any correlation between these covariates and PK/PD parameters'). No formal covariate analysis was attempted because of the small sample size (N=30). Baseline distribution (Table I): mean 57.9, SD 10.5, range 32.0-75.0 year.",
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Sex (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Screened graphically; no correlation. Baseline distribution (Table I): 60% female (18/30).",
      source_name        = "SEX"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened graphically; no correlation. Baseline distribution (Table I): mean 69.3, SD 12.1, range 42.0-90.0 kg.",
      source_name        = "WT"
    ),
    CREAT = list(
      description        = "Serum creatinine (reported indirectly via Cockcroft-Gault CRCL)",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "The paper reports Cockcroft-Gault creatinine clearance (Table I: mean 85.3, SD 32.9, range 23.2-150.0 mL/min, truncated at 150) rather than the raw serum creatinine. Screened against PK/PD parameters; no correlation.",
      source_name        = "SCR"
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened graphically; no correlation. Baseline distribution (Table I): mean 33.3, SD 10.6, range 14.9-50.5 g/L.",
      source_name        = "ALB"
    ),
    ALT = list(
      description        = "Serum alanine aminotransferase",
      units              = "IU/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened graphically; no correlation. Baseline distribution (Table I): mean 45.1, SD 20.1, range 19.0-100 IU/L.",
      source_name        = "ALT"
    ),
    AST = list(
      description        = "Serum aspartate aminotransferase",
      units              = "IU/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened graphically; no correlation. Baseline distribution (Table I): mean 38.9, SD 19.4, range 10.0-83.0 IU/L.",
      source_name        = "AST"
    ),
    TBILI = list(
      description        = "Total serum bilirubin",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened graphically; no correlation. Baseline distribution (Table I): mean 0.6, SD 0.3, range 0.2-1.6 umol/L (as reported; the units column header in Table I reads 'umol/L' but the magnitudes are more consistent with mg/dL -- confirm the column header against the published PDF before relying on the absolute scale).",
      source_name        = "TBIL"
    ),
    HGB = list(
      description        = "Hemoglobin",
      units              = "g/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened graphically; no correlation. Baseline distribution (Table I): mean 11.3, SD 1.5, range 6.4-13.0 g/dL.",
      source_name        = "HGB"
    ),
    HCT = list(
      description        = "Hematocrit",
      units              = "(fraction)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened against PK/PD parameters in Valenzuela 2011 (Results, Pharmacodynamics paragraph) but no baseline distribution is tabulated in the paper. No correlation detected.",
      source_name        = "HCT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 30L,
    n_studies      = 1L,
    age_range      = "32-75 years",
    age_mean       = "57.9 (SD 10.5) years",
    weight_range   = "42-90 kg",
    weight_mean    = "69.3 (SD 12.1) kg",
    bsa_mean       = "1.7 (SD 0.2) m^2 (range 1.4-2.0)",
    sex_female_pct = 60,
    race_ethnicity = NA_character_,
    disease_state  = "Adults with peritoneal carcinomatosis (primary tumour: ovarian n=10, colorectal n=9, appendiceal n=5, gastric n=3, endometrial n=2, primary papillary n=1) undergoing cytoreductive surgery followed by hyperthermic intraperitoneal oxaliplatin (HIO). World Health Organization performance status 0-2; life expectancy >=3 months. Eligibility required normal hepatic and renal function (bilirubin <=1.5x ULN, AST/ALT <=2.5x ULN, serum creatinine <=1.5x ULN) and acceptable bone marrow function (WBC >3.5e9/L, neutrophils >1.5e9/L, platelets >100e9/L).",
    dose_range     = "Single hyperthermic intraperitoneal oxaliplatin dose of 360 mg/m^2 administered in 4% icodextrin perfusate (2.5-6.0 L volume) at perfusate temperature 42-43 degC; HIO duration 30-60 min (mean 40 min). All patients additionally received intraperitoneal 5-FU 15 mg/kg over a 1-h infusion on each of the 5 postoperative days; 5-FU exposure was assumed negligible for the neutropenia model (low intrinsic neutropenic effect; 81.5% of 5-FU plasma concentrations below 0.04 mg/L LOQ).",
    regions        = "USP Hospital San Jaime (Torrevieja, Spain); enrollment 2006-2009",
    n_obs_peritoneum = 140L,
    n_obs_plasma     = 338L,
    n_obs_anc        = 678L,
    notes          = "Single-arm Phase 1-2 safety / tolerability / PK / PD study; baseline characteristics in Table I. Peritoneal cancer index (Sugarbaker): mean 8.6, range 0.0-39.0. Complete cytoreduction (CC0) achieved in 66.7%. Cockcroft-Gault creatinine clearance: mean 85.3 mL/min (SD 32.9, range 23.2-150.0; values >150 truncated to 150). Liver metastases present in 16.7%."
  )

  ini({
    # ------------------------------------------------------------------
    # Pharmacokinetic parameters - Valenzuela 2011 Table II (page 76),
    # 'Original dataset' column. Bioavailability F was fixed to 1, so
    # the reported Cl/F, Vc/F, Q2/F, Vp/F are apparent values.
    #
    # The paper parameterises peritoneal absorption as peritoneum-to-
    # plasma clearance Qa (L/h) and peritoneum volume Va (L), with the
    # first-order absorption rate constant ka = Qa/Va computed as a
    # secondary parameter (Methods, "Pharmacokinetic and Pharmacodynamic
    # Model" paragraph; Eq 1: dA/dt = -Qa/Va * A). nlmixr2lib has no
    # canonical (lqa, lva) pair, so the equivalent (lka, lvd)
    # canonical-name parameterisation is used here. The peritoneum
    # volume Va maps onto `lvd` (canonical "paper-mechanistic volume
    # term where the standard vc shape does not apply") and the
    # peritoneum-to-plasma clearance Qa is derived inside model() as
    # `qa <- ka * vd` for downstream use. The IIVs are encoded as a
    # correlated (etalka, etalvd) block (see the IIV block below) so
    # that the marginal distributions of the derived (Qa, Va) reproduce
    # the paper's independent omega_Qa and omega_Va exactly.
    # ------------------------------------------------------------------
    lka <- log(0.3242);  label("First-order absorption rate ka = Qa/Va from peritoneum to plasma (1/h)")  # Table II Qa=2.70 / Va=8.33; secondary param per Methods
    lvd <- log(8.33);    label("Peritoneum volume of distribution Va (L)")                                # Table II Va
    lcl <- log(1.61);    label("Apparent systemic clearance Cl/F from plasma (L/h)")                      # Table II Cl/F
    lq  <- log(77.0);    label("Apparent intercompartmental clearance Q2/F (L/h)")                        # Table II Q2/F
    lvc <- log(19.2);    label("Apparent central volume of distribution Vc/F (L)")                        # Table II Vc/F
    lvp <- log(72.8);    label("Apparent peripheral volume of distribution Vp/F (L)")                     # Table II Vp/F

    # ------------------------------------------------------------------
    # Pharmacodynamic parameters - Valenzuela 2011 Table II, system-
    # and drug-related blocks.
    #
    # MTT is fixed to 118 h from Friberg 2002 (the paper's reference 30)
    # because the ANC data lacked information to estimate it (Methods,
    # "Pharmacokinetic and Pharmacodynamic Model" paragraph).
    #
    # kCirc is fixed to 0.07 per h from Friberg 2002 ("the population
    # mean half life of neutrophils previously determined, 0.07 h^-1
    # (30)") and is kept as a separate parameter from ktr (Eq 8:
    # dCirc/dt = ktr * Transit3 - kCirc * Circ). Figure 1 schematic
    # labels k_circ as "(= k_tr)" -- that annotation is the standard
    # Friberg cartoon and contradicts the text's explicit "0.07 h^-1"
    # statement; we follow the text + Eq 8, with steady-state initial
    # conditions set in model() to maintain a self-consistent baseline.
    # See vignette Errata for the documented figure-vs-text discrepancy.
    # ------------------------------------------------------------------
    lcirc0 <- log(7.05);        label("Baseline absolute neutrophil count Circ0 (10^9 cells/L)")          # Table II Circ0
    lmtt   <- fixed(log(118));  label("Mean transit time MTT through the proliferation -> circulation chain (h)")  # Table II 'MTT 118 Fixed'; ref 30 (Friberg 2002)
    gamma  <- 0.135;            label("Feedback exponent gamma on (Circ0/Circ) (unitless)")              # Table II gamma
    lalpha <- log(0.182);       label("Linear drug-effect slope alpha on the proliferation rate (L/mg); Edrug = alpha * Cc")  # Table II alpha; Eq 9
    kcirc  <- fixed(0.07);      label("Circulating-cell elimination rate constant kCirc (1/h); fixed from Friberg 2002")  # text + ref 30 (Friberg 2002); see comments above re: Fig 1 vs text

    # ------------------------------------------------------------------
    # Inter-individual variability - Valenzuela 2011 Table II 'Inter-
    # individual variability (CV %)' rows. Each CV% is converted to a
    # log-normal variance via omega^2 = log(CV^2 + 1).
    #
    # The (etalka, etalvd) block is a reparameterisation of the paper's
    # independent etas on (Qa, Va):
    #   etalqa = etalka + etalvd  and  etalva = etalvd
    # If etalqa and etalva are independent in the paper, then the
    # equivalent (etalka, etalvd) joint distribution is:
    #   var(etalka)         = var(etalqa) + var(etalva)
    #   cov(etalka, etalvd) = -var(etalva)
    #   var(etalvd)         =  var(etalva)
    # This block reproduces the paper's marginal omega_Qa = 34.1 % and
    # omega_Va = 17.7 % and zero (Qa, Va) covariance after the
    # parameterisation change. omega_Q2 is reported as '-' in Table II (no
    # IIV estimated on Q2).
    # ------------------------------------------------------------------
    etalka + etalvd ~ c(0.14091, -0.030857, 0.030857)
                                                                                              # omega_Qa=34.1% -> var(lqa)=log(1+0.341^2)=0.11005; omega_Va=17.7% -> var(lva)=log(1+0.177^2)=0.030857
    etalcl    ~ 0.55015                                                                       # omega_Cl/F = 85.6 % -> log(1 + 0.856^2) = 0.55015
    etalvc    ~ 0.28903                                                                       # omega_Vc/F = 57.9 % -> log(1 + 0.579^2) = 0.28903
    etalvp    ~ 0.053744                                                                      # omega_Vp/F = 23.5 % -> log(1 + 0.235^2) = 0.053744

    etalcirc0 ~ 0.16467                                                                       # omega_Circ0 = 42.3 % -> log(1 + 0.423^2) = 0.16467
    etalmtt   ~ 0.10215                                                                       # omega_MTT   = 32.8 % -> log(1 + 0.328^2) = 0.10215 (typical fixed; IIV estimated)
    etalalpha ~ 1.09504                                                                       # omega_alpha = 141  % -> log(1 + 1.41^2)  = 1.09504

    # ------------------------------------------------------------------
    # Residual error - Valenzuela 2011 Table II 'Residual variability
    # (CV %)' rows. The paper Methods state that residual variability
    # was modelled as an additive error on the natural-log-transformed
    # observations and predictions; the resulting magnitudes are
    # reported as approximate coefficients of variation (Statistical
    # Model paragraph). This LTBS additive-on-log encoding is
    # numerically equivalent to a proportional residual on the linear
    # scale for small sigma, and is encoded here with `prop()` to match
    # the Friberg_2002_paclitaxel.R precedent.
    # ------------------------------------------------------------------
    propSd       <- 0.147; label("Proportional residual error on plasma oxaliplatin Cc (fraction)")             # Table II sigma_1 (plasma)
    propSd_Cperi <- 0.165; label("Proportional residual error on peritoneal oxaliplatin Cperi (fraction)")      # Table II sigma_2 (peritoneum)
    propSd_ANC   <- 0.497; label("Proportional residual error on absolute neutrophil count ANC (fraction)")     # Table II sigma   (ANC)
  })

  model({
    # 1. Individual PK parameters
    ka <- exp(lka + etalka)
    vd <- exp(lvd + etalvd)
    cl <- exp(lcl + etalcl)
    q  <- exp(lq)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)

    # 2. Individual PD parameters
    circ0 <- exp(lcirc0 + etalcirc0)
    mtt   <- exp(lmtt   + etalmtt)
    alpha <- exp(lalpha + etalalpha)

    # 3. Derived rate constants and baseline steady-state amounts
    #    ktr  = (n_transit + 1) / MTT for a chain of 3 transit
    #           compartments (n = 3), so ktr = 4 / MTT (Eq 10).
    #    kprol = ktr at steady state (Methods text: "dProl/dt = 0
    #           and therefore kProl = ktr").
    #    Steady-state precursor amounts: with kCirc != ktr (paper's
    #    fixed kCirc = 0.07 vs ktr = 4/118 = 0.0339), the steady-state
    #    requirement dCirc/dt = 0 implies Transit3 = (kCirc/ktr) *
    #    Circ0, and the upstream-equality d/dt(Transit_i) = 0 chains
    #    that ratio back through Transit2, Transit1, and Prol.
    ktr     <- 4 / mtt
    prol_ss <- circ0 * kcirc / ktr

    # 4. Observed concentrations
    Cperi <- depot   / vd
    Cc    <- central / vc

    # 5. Drug-effect on the proliferation rate (Eq 9) and feedback
    #    amplification term (Eq 4). Edrug = alpha * Cc; values > 1
    #    indicate the proliferation rate is being driven negative
    #    (net killing of precursor cells) -- the paper reports 50th /
    #    95th peak-Edrug percentiles of 0.44 and 5.64.
    edrug <- alpha * Cc
    feed  <- (circ0 / circ)^gamma

    # 6. ODE system (Eq 1-3 PK, Eq 4-8 PD)
    d/dt(depot)       <- -ka * depot                                                          # Eq 1: dA/dt = -ka * A
    d/dt(central)     <-  ka * depot   - (cl + q) / vc * central + q / vp * peripheral1       # Eq 2
    d/dt(peripheral1) <-  q  / vc * central - q / vp * peripheral1                            # Eq 3
    d/dt(precursor1)  <-  ktr * precursor1 * feed * (1 - edrug) - ktr * precursor1            # Eq 4: Prol with self-replication, feedback, drug effect (kProl = ktr)
    d/dt(precursor2)  <-  ktr * precursor1 - ktr * precursor2                                 # Eq 5: Transit1
    d/dt(precursor3)  <-  ktr * precursor2 - ktr * precursor3                                 # Eq 6: Transit2
    d/dt(precursor4)  <-  ktr * precursor3 - ktr * precursor4                                 # Eq 7: Transit3
    d/dt(circ)        <-  ktr * precursor4 - kcirc * circ                                     # Eq 8: Circ with kCirc separate from ktr

    # 7. Initial conditions: circulating cells at observed baseline;
    #    precursor pools at the (kCirc / ktr) * Circ0 steady-state
    #    amount so the chain is fully balanced at t = 0.
    precursor1(0) <- prol_ss
    precursor2(0) <- prol_ss
    precursor3(0) <- prol_ss
    precursor4(0) <- prol_ss
    circ(0)       <- circ0

    # 8. Observation model. Outputs:
    #      Cc    = plasma oxaliplatin concentration (mg/L)
    #      Cperi = peritoneal oxaliplatin concentration (mg/L)
    #      ANC   = absolute neutrophil count (10^9 cells/L)
    ANC <- circ
    Cc    ~ prop(propSd)
    Cperi ~ prop(propSd_Cperi)
    ANC   ~ prop(propSd_ANC)
  })
}
