Troy_2020_arylsulfataseA <- function() {
  description <- "Four-compartment population PK model for recombinant human arylsulfatase A (TAK-611, formerly SHP611) following intrathecal administration in children with metachromatic leukodystrophy (Troy 2020): a two-compartment CNS subsystem (intrathecal dose enters csf, exchanging with a putative brain-tissue compartment cns_tissue via Q_CSF) draining through a hypothetical transit compartment into a one-compartment serum disposition. CSF and CNS distribution volumes are proportional to age-predicted physiologic volumes (Matsuzawa 2001 regressions); serum CL and V_central are allometrically scaled on body weight."
  reference <- "Troy S, Wasilewski M, Beusmans J, Godfrey CJ. Pharmacokinetic Modeling of Intrathecally Administered Recombinant Human Arylsulfatase A (TAK-611) in Children With Metachromatic Leukodystrophy. Clin Pharmacol Ther. 2020;107(6):1394-1404. doi:10.1002/cpt.1752"
  vignette <- "Troy_2020_arylsulfataseA"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL"
  )

  # cns_tissue is the compartment peripheral to csf that Troy 2020 introduced
  # to fit the CSF trough data: "The inclusion of a compartment peripheral to
  # CSF, potentially representing the central nervous system (CNS) brain
  # tissue, greatly improved the prediction of CSF data" (Results). Its volume
  # is scaled to the age-predicted white-matter + gray-matter volume, so it is
  # a lumped brain-parenchyma state rather than any single canonical
  # brain_<region>, and it is peripheral to csf rather than to central (so
  # peripheral1 would misstate the topology). Same role and same name as the
  # lumped CNS-tissue compartment in Luu_2017_nusinersen.R.
  paper_specific_compartments <- c("cns_tissue")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Amounts are in mg throughout (the intrathecal dose unit);
  # concentrations are formed in model() with an explicit mg/L -> ng/mL factor
  # of 1000, mirroring the source control stream's S1 = VCSF/1000 and
  # S2 = VCENT/1000 scaling (Troy 2020 Data S1 $PK).
  compartmentData <- list(
    csf        = list(analyte = "arylsulfatase A (TAK-611)", units = "mg", specimen = "CSF", verified = TRUE),
    central    = list(analyte = "arylsulfatase A (TAK-611)", units = "mg", specimen = "serum", verified = TRUE),
    transit1   = list(analyte = "arylsulfatase A (TAK-611)", units = "mg", specimen = "administration site", verified = TRUE),
    cns_tissue = list(analyte = "arylsulfatase A (TAK-611)", units = "mg", specimen = "CSF", verified = TRUE)
  )

  covariateData <- list(
    AGE = list(
      description        = "Chronological age at the first intrathecal dose",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Troy 2020 Data S1 ($PK) carries age in MONTHS and advances it over the study:",
        "AGEC = AGE + TAFD/24/30. This model file keeps the register-canonical AGE in YEARS",
        "and converts inside model() (agec <- 12 * AGE + t / (24 * 30)), so the event table",
        "must supply AGE in years and must place the first dose at t = 0 for agec to track",
        "time after first dose. agec drives the Matsuzawa 2001 age regressions for physiologic",
        "CSF, white-matter and gray-matter volumes, which in turn set V_CSF and V_CNS.",
        "Cohort baseline age: median 36.5 months (range 19.0-107), i.e. 3.04 years (1.58-8.92)."
      ),
      source_name        = "AGE"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline body weight, allometrically scaling the systemic components only:",
        "TV_CL = theta_CL * (WT/15)^0.75 and TV_Vcentral = theta_Vcentral * (WT/15)",
        "(Troy 2020 Methods display equations; Data S1 $PK CLscale / Vscale).",
        "Reference weight 15 kg. Time-varying weight was tested and REJECTED",
        "(OFV increased by 3.649; Troy 2020 Results, Covariate analysis), so this column is",
        "baseline-only. Cohort median baseline weight 14.1 kg (range 10.5-24.8)."
      ),
      source_name        = "WT"
    )
  )

  # Covariates Troy 2020 screened but did NOT retain in the final model. Kept
  # here for provenance; none is referenced in model().
  covariatesDataExcluded <- list(
    FORM_arylsulfataseA_processB = list(
      description = "Indicator for TAK-611 manufactured by the revised process B (raised mannose-6-phosphate and sialic acid) versus the original process A",
      units       = "(binary)",
      type        = "binary",
      reference_category = "0 (process A)",
      notes       = paste(
        "Troy 2020 Results, Covariate analysis: relative bioavailability of process-B material",
        "was estimated at 0.948 (95% CI 0.546, 1.350), 'not significantly different from that of",
        "TAK-611 produced using process A (reference value of 1)'. Not carried in the final model",
        "(absent from Table 1 and from the Data S1 control stream), so no bioavailability term is",
        "encoded. Cohort 4 (n = 6) received process B; cohorts 1-3 (n = 18) received process A."
      )
    ),
    WT_TIMEVARYING = list(
      description = "Time-varying body weight, linearly interpolated between quarterly measurements",
      units       = "kg",
      type        = "continuous",
      reference_category = NULL,
      notes       = paste(
        "Troy 2020 Results, Covariate analysis: replacing baseline weight with continuously",
        "interpolated weight 'did not improve the fit of the data (objective function value",
        "increased by 3.649)'. The final model uses baseline WT only."
      )
    ),
    ADA_TITER = list(
      description = "Anti-drug antibody titer in serum and in CSF",
      units       = "titer",
      type        = "continuous",
      reference_category = NULL,
      notes       = paste(
        "Troy 2020 Results, Covariate analysis and Figure S1: ADA titer and neutralizing activity",
        "were assessed GRAPHICALLY ONLY and were never entered as formal covariates ('Although ADA",
        "titers were not included formally as covariates in the model...', Discussion). Five subjects",
        "had week-38 serum titers > 1,000 and four of those showed reduced serum TAK-611; no",
        "neutralizing ADAs were detected in CSF. No coefficient is published, so nothing is encoded.",
        "The Data S1 control stream reads SERADAI / CSFADAI but never references them in $PK."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 24L,
    n_studies      = 1L,
    studies        = "NCT01510028, phase I/II multicenter open-label dose-escalation",
    age_range      = "19.0-107 months (1.58-8.92 years) at baseline",
    age_median     = "36.5 months (3.04 years)",
    weight_range   = "10.5-24.8 kg",
    weight_median  = "14.1 kg",
    disease_state  = "Children with metachromatic leukodystrophy (late-infantile onset; first symptoms at or before 30 months of age, ambulatory at screening)",
    dose_range     = "10, 30 or 100 mg intrathecal every other week for 38 weeks (up to 20 doses); cohort 4 received 100 mg of process-B material",
    administration_routes = "Intrathecal, via a surgically implanted intrathecal drug delivery device (or lumbar puncture when the device was unusable)",
    regions        = "Multicenter (North America and Europe)",
    notes          = paste(
      "Four dose cohorts of n = 6 (Troy 2020 Results). Model built on 321 CSF samples",
      "(median 11 per subject, range 5-19) and 387 serum samples (median 11.5, range 1-17);",
      "60 CSF (18.7%) and 117 serum (30.2%) samples were below the limit of quantification and",
      "were EXCLUDED rather than censored. CSF samples are predose troughs; serum was sampled",
      "over 48 h at week 0 and week 38. NONMEM 7.3 ADVAN13, FOCE with interaction, log-transformed",
      "concentrations. Sex and race distributions are not reported in the paper."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Structural fixed effects: Troy 2020 Table 1 ("Fixed effects").
    # The control stream (Data S1 $PK) mu-references theta3-theta6 on the log
    # scale (MU_1..MU_4 with EXP(MU_n + ETA)), and Table 1 reports the
    # back-transformed values; those back-transformed values are what is
    # wrapped in log() here.
    #
    # theta1 and theta2 are DIMENSIONLESS PROPORTIONALITY FACTORS, not
    # volumes. Data S1 $PK: TVVCSF = THETA(V_CSF) * CSFVOL and
    # TVVCNS = THETA(V_CNS) * CNSTVOL, where CSFVOL and CNSTVOL are the
    # age-predicted physiologic volumes computed in model(). Reading 0.183
    # and 1.69 as litres would build a completely different model; the
    # resulting volumes at the cohort median age of 36.5 months are
    # V_CSF = 0.0268 L and V_CNS = 1.566 L, matching the individual-level
    # medians of 0.0272 L and 1.58 L in Table 2.
    lvcsf <- log(0.183);   label("Proportionality of V_CSF to the age-predicted physiologic CSF volume (unitless multiplier)")                    # Troy 2020 Table 1 theta1 = 0.183 (RSE 76.1%; 95% CI 0.0889, 1.72)
    lvcns <- log(1.69);    label("Proportionality of V_CNS to the age-predicted white-matter + gray-matter volume (unitless multiplier)")         # Troy 2020 Table 1 theta2 = 1.69 (RSE 121%; 95% CI 0.421, 3.97)
    lqcsf <- log(0.00280); label("Intercompartmental clearance between csf and cns_tissue (L/h)")                                                 # Troy 2020 Table 1 theta3 = 0.00280 L/h (RSE 22.3%; 95% CI 0.000949, 0.0222)
    lktr  <- log(0.581);   label("First-order transit rate constant from csf to serum (1/h)")                                                     # Troy 2020 Table 1 theta4 K_trans = 0.581 /h (RSE 40.2%; 95% CI 0.0917, 0.809)
    lcl   <- log(3.01);    label("Systemic clearance from the central compartment at the reference weight of 15 kg (L/h)")                        # Troy 2020 Table 1 theta5 = 3.01 L/h (RSE 19.3%; 95% CI 2.10, 3.40)
    lvc   <- log(69.1);    label("Central (serum) volume of distribution at the reference weight of 15 kg (L)")                                   # Troy 2020 Table 1 theta6 V_central = 69.1 L (RSE 5.55%; 95% CI 0.963, 81.0)

    # ---------------------------------------------------------------------
    # Allometric exponents on the SYSTEMIC components, reference weight 15 kg.
    # Troy 2020 Methods prints both display equations:
    #   TV_Vcentral = theta_Vcentral * (WT / 15 kg)
    #   TV_CL       = theta_CL       * (WT / 15 kg)^0.75
    # The exponents are structural constants of the allometric model (Data S1
    # $PK: CLscale = 0.75 * LOG(WTNORM), Vscale = 1 * LOG(WTNORM)); neither
    # appears in Table 1 and neither carries an RSE or CI, so both are fixed.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on (WT/15) for systemic clearance (unitless)")                                             # Troy 2020 Methods display equation 2; Data S1 $PK CLscale
    e_wt_vc <- fixed(1);    label("Allometric exponent on (WT/15) for the central volume (unitless)")                                             # Troy 2020 Methods display equation 1; Data S1 $PK Vscale

    # ---------------------------------------------------------------------
    # Interindividual variability: Troy 2020 Table 1 ("Interindividual
    # variability random effects"), diagonal Omega (Data S1 $PROB: "IIV=Diag").
    # The Estimate column holds the VARIANCE: the paper's CV column is
    # sqrt(omega^2) for every one of the six rows (e.g. sqrt(0.409) = 63.9%,
    # sqrt(0.176) = 42.0%, sqrt(1.05e-08) = 0.0102%), so no CV -> variance
    # back-transformation is applied.
    etalqcsf ~ 0.0607    # Troy 2020 Table 1 omega(1,1) = 0.0607 (CV 24.6%, shrinkage 78.4%)
    etalktr  ~ 0.176     # Troy 2020 Table 1 omega(2,2) = 0.176  (CV 42.0%, shrinkage 18.7%)
    etalcl   ~ 1.05e-08  # Troy 2020 Table 1 omega(3,3) = 1.05e-08 (CV 0.0102%, shrinkage 100%) -- estimated to the boundary; retained as published
    etalvc   ~ 0.162     # Troy 2020 Table 1 omega(4,4) = 0.162  (CV 40.3%, shrinkage 22.3%)
    etalvcsf ~ 0.409     # Troy 2020 Table 1 omega(5,5) = 0.409  (CV 63.9%, shrinkage 11.6%)
    etalvcns ~ 0.425     # Troy 2020 Table 1 omega(6,6) = 0.425  (CV 65.2%, shrinkage 56.7%)

    # ---------------------------------------------------------------------
    # Residual error. Troy 2020 Methods: "Concentrations were log-transformed
    # before analysis and an additive residual error model on the log-scale was
    # used", implemented in Data S1 $ERROR as
    #   IPRED = LOG(F + DEL); Y = IPRED + INDCSF*EPS(1) + INDCENT*EPS(2)
    # which is exactly nlmixr2's lnorm() error, hence expSd rather than propSd
    # (the proportional approximation to a log-additive residual is not usable
    # at these magnitudes -- the CSF SD is ~98%).
    # Table 1 reports sigma as the VARIANCE in the Estimate column and the SD in
    # the CV column (footnote b, "Reported as SD"): sqrt(0.963) = 0.981 and
    # sqrt(0.456) = 0.675, matching the printed SDs exactly.
    expSd      <- sqrt(0.456); label("Serum residual SD on the log scale (log ng/mL)")                                                            # Troy 2020 Table 1 sigma(2,2) = 0.456 variance; SD 0.675 (RSE 6.09%, shrinkage 5.01%)
    expSd_Ccsf <- sqrt(0.963); label("CSF residual SD on the log scale (log ng/mL)")                                                              # Troy 2020 Table 1 sigma(1,1) = 0.963 variance; SD 0.981 (RSE 11.5%, shrinkage 5.14%)
  })

  model({
    # -------------------------------------------------------------------
    # 1. Derived covariate terms
    # -------------------------------------------------------------------
    # Age in MONTHS, advancing over the study. Data S1 $PK:
    #   AGEC = AGE + TAFD/24/30
    # with AGE in months and TAFD in hours. AGE is supplied here in the
    # register-canonical unit of years, so it is multiplied by 12; t is time
    # after the first dose in hours, so t/(24*30) converts to months.
    agec <- 12 * AGE + t / (24 * 30)

    # Age-predicted physiologic brain and CSF volumes (L). These are the
    # Matsuzawa 2001 regressions for healthy infants and children aged
    # 1-120 months, transcribed verbatim from Troy 2020 Data S1 $PK, which
    # cites them as "Cereb. Cortex (2001) 11 (4): 335-342,
    # doi: 10.1093/cercor/11.4.335". Troy 2020 Results credits the same
    # source (reference 20) for "data on volumetric changes in CSF and white
    # and gray matter in children aged 1-120 months". The published
    # regressions give volumes in mL, hence the /1000 to litres.
    zgm    <- (agec / 100)^(-0.5) - 1.462
    gmvol  <- (653 - 33.54 * zgm) / 1000
    z1     <- log(agec / 100) + 0.7596
    wmvol  <- (298.28 + 72.74 * z1) / 1000
    csfvol <- (149.88 + 13.90 * z1) / 1000
    cnsvol <- wmvol + gmvol

    # -------------------------------------------------------------------
    # 2. Individual parameters
    # -------------------------------------------------------------------
    # CNS volumes are the estimated proportionality factor times the
    # age-predicted physiologic volume, so vcsf and vcns are in litres even
    # though lvcsf and lvcns are dimensionless (Data S1 $PK: TVVCSF, TVVCNS).
    vcsf <- exp(lvcsf + etalvcsf) * csfvol
    vcns <- exp(lvcns + etalvcns) * cnsvol
    qcsf <- exp(lqcsf + etalqcsf)
    ktr  <- exp(lktr  + etalktr)
    # Systemic components carry the weight allometry, reference weight 15 kg.
    cl   <- exp(lcl   + etalcl)   * (WT / 15)^e_wt_cl
    vc   <- exp(lvc   + etalvc)   * (WT / 15)^e_wt_vc

    # -------------------------------------------------------------------
    # 3. Micro-constants (Data S1 $PK reparameterisation to rate constants)
    # -------------------------------------------------------------------
    k14 <- qcsf / vcsf   # csf -> cns_tissue
    k41 <- qcsf / vcns   # cns_tissue -> csf
    kel <- cl / vc       # systemic elimination out of central

    # -------------------------------------------------------------------
    # 4. ODE system (Troy 2020 Figure 3; Data S1 $DES). States are declared in
    #    the order of the source control stream's $MODEL block so the mapping
    #    is 1:1: A(1) csf, A(2) central, A(3) transit, A(4) cns. The
    #    intrathecal dose enters csf ($MODEL: "COMP = (CSF DEFDOSE)"), so
    #    dosing records must carry cmt = "csf".
    #    The transit chain uses one rate constant in both directions
    #    (K13 = K32 = KTRANS), so it is a pure delay and adds no clearance.
    # -------------------------------------------------------------------
    d/dt(csf)        <- -k14 * csf + k41 * cns_tissue - ktr * csf
    d/dt(central)    <-  ktr * transit1 - kel * central
    d/dt(transit1)   <-  ktr * csf - ktr * transit1
    d/dt(cns_tissue) <-  k14 * csf - k41 * cns_tissue

    # -------------------------------------------------------------------
    # 5. Observations. States hold mg and volumes are in L, so mg/L is scaled
    #    by 1000 to give ng/mL -- the source control stream folds the same
    #    factor into its scaling terms (Data S1 $PK: S1 = VCSF/1000,
    #    S2 = VCENT/1000; $ERROR: CCSF = A(1)/S1, CSERUM = A(2)/S2).
    #    Ccns is derived for simulation only; TAK-611 was never measured in
    #    brain tissue.
    # -------------------------------------------------------------------
    Cc   <- 1000 * central    / vc
    Ccsf <- 1000 * csf        / vcsf
    Ccns <- 1000 * cns_tissue / vcns

    Cc   ~ lnorm(expSd)
    Ccsf ~ lnorm(expSd_Ccsf)
  })
}
