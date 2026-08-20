Dias_2024_quetiapine_rat <- function() {
  description <- paste0(
    "Preclinical (rat, Wistar). Semimechanistic pharmacodynamic (PD-only) ",
    "precursor-pool model for medial prefrontal cortex (mPFC) extracellular ",
    "dopamine after a single 5 mg/kg intravenous dose of quetiapine given ",
    "either as a solution (FQ) or as lipid core nanocapsules (QLNC), in ",
    "naive rats and in poly(i:c) schizophrenia phenotyped rats (SPR). ",
    "Dopamine precursors are synthesised into a pool at a zero-order rate ",
    "and released into the extracellular dopamine compartment by a ",
    "first-order rate constant that quetiapine stimulates linearly; ",
    "dopamine is removed by a reuptake process that a negative-feedback ",
    "modulator exacerbates, giving the fast return to baseline. Quetiapine ",
    "acts through a Sheiner effect compartment driven by the unbound brain ",
    "concentration, and co-administered nanoparticles competitively blunt ",
    "the drug effect while decaying first-order. Disease severity enters as ",
    "the continuous prepulse-inhibition score on both the dopamine baseline ",
    "and the drug-effect slope. This model has NO PK layer: the unbound ",
    "quetiapine brain concentration must be supplied by the user as the ",
    "time-varying covariate CU_QTP_BRAIN, generated from the companion ",
    "semimechanistic popPK model (Carreno et al. 2020, ",
    "doi:10.1124/jpet.120.000109)."
  )
  reference <- paste(
    "Dias BB, Carreno F, Helfer VE, Olivo LB, Staudt KJ, Paese K,",
    "Barreto F, Meyer FS, Herrmann AP, Guterres SS, Rates SMK,",
    "de Araujo BV, Troconiz IF, Dalla Costa T.",
    "Pharmacokinetic/pharmacodynamic modeling of cortical dopamine",
    "concentrations after quetiapine lipid core nanocapsules administration",
    "to schizophrenia phenotyped rats.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(4):638-648.",
    "doi:10.1002/psp4.13107.",
    "Upstream PK driver (not implemented here):",
    "Carreno F et al. J Pharmacol Exp Ther. 2020;375(1):49-58,",
    "doi:10.1124/jpet.120.000109.",
    sep = " "
  )
  vignette <- "Dias_2024_quetiapine_rat"
  units <- list(
    time = "h",
    dosing = "n/a (PD-only; unbound brain quetiapine supplied externally via CU_QTP_BRAIN)",
    concentration = "ng/mL"
  )

  # Mechanistic states named by the source paper (Figure 2 / NM-TRAN $MODEL)
  # that do not map onto a canonical compartment role.
  paper_specific_compartments <- c(
    "auc_dopamine", "modulator", "nano"
  )

  covariateData <- list(
    CU_QTP_BRAIN = list(
      description = paste0(
        "Unbound quetiapine concentration in medial prefrontal cortex ",
        "interstitial fluid (ng/mL), supplied externally as a time-varying ",
        "covariate. This is the driver of the whole PD system."
      ),
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Time-varying. The paper contains no PK layer: Methods state that ",
        "'typical unbound brain QTP concentration profiles for each ",
        "experimental group were generated from the semimechanistic PK ",
        "model developed previously (Figure S2) and were used as the ",
        "drivers of the DA response'. In the NM-TRAN control stream ",
        "(Supplementary Material S1) the value arrives as the data item ",
        "CBRAIN and is assigned CB = CBRAIN in $PK, then carried into $DES ",
        "as CONC. The interpolation block in $PK (LCP / LTM / SLOP) is ",
        "inert because LCP and LTM are never updated after the first ",
        "record, so CONC collapses to the record's CBRAIN value; rxode2's ",
        "covariate interpolation supplies the same driver. Profiles for all ",
        "four experimental groups are plotted in Figure 4a and Figure S2. ",
        "Source of the profiles: Carreno et al. 2020 JPET, ",
        "doi:10.1124/jpet.120.000109."
      ),
      source_name = "CBRAIN"
    ),
    SCORE_PPI = list(
      description = paste0(
        "Prepulse inhibition of the acoustic startle response (%), the ",
        "rodent sensorimotor-gating measure used as the continuous ",
        "disease-severity covariate for the schizophrenia phenotype."
      ),
      units = "%",
      type = "continuous",
      reference_category = "33.4 (PPImd, the median PPI across the whole animal study population)",
      notes = paste0(
        "Centred at PPImd = 33.4% and entered linearly on both the dopamine ",
        "baseline (Eq 7) and the quetiapine effect slope (Eq 6). Group ",
        "medians reported in Supplementary Material S1: 49.7% for naive ",
        "animals and 23.6% for SPR. Higher PPI (less disease) gives a ",
        "higher dopamine baseline and a larger drug effect. Sex and body ",
        "weight were screened and not retained (p > 0.05)."
      ),
      source_name = "PPI"
    ),
    FORM_QTP_QLNC = list(
      description = paste0(
        "Quetiapine formulation indicator: 1 = lipid core nanocapsules ",
        "(QLNC), 0 = quetiapine solution (FQ)."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (FQ, quetiapine solution)",
      notes = paste0(
        "Gates the nanoparticle state: the NM-TRAN control stream sets ",
        "A_0(6) = DNANO and DADT(6) = -KELNANO * A(6) only when NANO == 1, ",
        "and uses the blunting term 1/(1 + NP) in the drug effect only when ",
        "NANO == 1. Because NP0 is multiplied by this indicator, nano(0) is ",
        "zero for FQ, the nano ODE is identically zero, and the blunting ",
        "term collapses to 1 -- so the single unified expression used in ",
        "model() reproduces both NM-TRAN branches exactly."
      ),
      source_name = "NANO"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male).",
      units = "(binary)",
      type = "binary",
      notes = paste0(
        "Screened in the stepwise covariate analysis and not retained ",
        "(Results: 'sex and weight did not impact model parameters ",
        "(p > 0.05)'). Both sexes were studied in every experimental group."
      )
    ),
    WT = list(
      description = "Body weight (kg).",
      units = "kg",
      type = "continuous",
      notes = paste0(
        "Screened in the stepwise covariate analysis and not retained ",
        "(Results: 'sex and weight did not impact model parameters ",
        "(p > 0.05)')."
      )
    )
  )

  compartmentData <- list(
    precursor1 = list(
      analyte = "dopamine precursor pool",
      units = "ng/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    dopamine = list(
      analyte = "dopamine",
      units = "ng/mL",
      specimen = "brain ISF",
      verified = TRUE
    ),
    auc_dopamine = list(
      analyte = "dopamine",
      units = "ng/mL*h",
      specimen = "not applicable",
      verified = TRUE
    ),
    modulator = list(
      analyte = "dopamine reuptake negative-feedback modulator",
      units = "unitless (fraction of baseline; 1 at steady state)",
      specimen = "not applicable",
      verified = TRUE
    ),
    effect = list(
      analyte = "quetiapine",
      units = "ng/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    nano = list(
      analyte = "quetiapine lipid core nanocapsules",
      units = "ng/mL",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "rat (Wistar; naive and prenatal poly(i:c) schizophrenia phenotyped)",
    n_subjects = 49L,
    n_studies = 1L,
    n_observations = 924L,
    sex_female_pct = 49.0,
    disease_state = paste0(
      "Naive rats and schizophrenia phenotyped rats (SPR). SPR offspring ",
      "were generated by a single 4 mg/kg i.v. poly(i:c) bolus to pregnant ",
      "dams at gestational day 15 (naive dams received 2 mL/kg saline); the ",
      "phenotype was confirmed in adult offspring (PND75) by the prepulse ",
      "inhibition of the acoustic startle response test."
    ),
    dose_range = paste0(
      "Single 5 mg/kg intravenous bolus of quetiapine via the lateral caudal ",
      "vein, given either as a solution (FQ, 5 mg/mL) or as quetiapine lipid ",
      "core nanocapsules (QLNC, 1 mg/mL)."
    ),
    regions = "Brazil (Federal University of Rio Grande do Sul, Porto Alegre)",
    notes = paste0(
      "Eight groups by disease status, sex and treatment (Methods, ",
      "'Animals'): FQ naive male n=7, FQ naive female n=6, FQ SPR male n=5, ",
      "FQ SPR female n=6, QLNC naive male n=7, QLNC naive female n=6, ",
      "QLNC SPR male n=6, QLNC SPR female n=6 (50 animals dosed). Pooled ",
      "into four analysis groups: FQ-naive n=13, FQ-SPR n=11, QLNC-naive ",
      "n=13, QLNC-SPR n=12 (Figure 1 caption), i.e. 49 animals contributing ",
      "924 dopamine microdialysate observations; one animal was excluded ",
      "for a PPI value discrepant from its group (75.38%). Dopamine in mPFC ",
      "extracellular fluid was sampled by intracerebral microdialysis ",
      "(CMA 12 probe, 3 mm PAES membrane, 20 kDa cutoff, artificial CSF at ",
      "1 uL/min): four 20-min baseline dialysates, then 15 further 20-min ",
      "dialysates to 280 min after dosing. No correction for in vitro probe ",
      "recovery was applied, so dopamine is reported as measured in the ",
      "dialysate. The sex_female_pct of 49.0 is computed from the group ",
      "sizes above (24 of 49 analysed animals are female)."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # All population values are FINAL estimates from Supplementary Material
    # S1, Table S1 ("Population parameters estimated for the
    # semi-mechanistic PK/PD model and bootstrap analysis"), cross-checked
    # against the NM-TRAN $PK theta assignments in the same supplement.
    # NM-TRAN theta mapping:
    #   THETA(1)=Pool0  THETA(2)=Kin   THETA(3)=DA0   THETA(4)=EQTP
    #   THETA(5)=residual SD          THETA(6)=Kmod  THETA(7)=Keo
    #   THETA(8)=KNP    THETA(9)=NP0  THETA(10)=theta_DA0,SPR
    #   THETA(11)=theta_EQTP,SPR
    # ---------------------------------------------------------------------

    # Structural PD parameters. All are estimated on the natural scale in
    # NONMEM with exponential IIV (THETA * EXP(ETA)), which is encoded here
    # as a log-scale typical value with additive eta.

    lrbase_precursor1 <- log(1.16)
    label("Baseline dopamine precursor pool concentration Pool0 (ng/mL)")
    # Table S1: Pool0 1.16 ng/mL (29% RSE); bootstrap median 1.16 (0.82-1.47); no IIV

    lkpin <- log(0.313)
    label("Zero-order dopamine precursor synthesis rate Kin (ng/mL/h)")
    # Table S1: Kin 0.313 ng/mL*h (22% RSE); bootstrap median 0.311 (0.225-0.428)

    lrbase_dopamine <- log(0.288)
    label("Baseline mPFC extracellular dopamine concentration DA0 (ng/mL)")
    # Table S1: DA0 0.288 ng/mL (3% RSE); bootstrap median 0.288 (0.274-0.304)

    leqtp <- log(37.7)
    label("Linear slope of the quetiapine effect on dopamine release EQTP (mL/ng)")
    # Table S1: EQTP 37.7 (66% RSE); bootstrap median 39.3 (33.4-59.4)

    lke0 <- log(0.418)
    label("Effect-compartment rate constant Keo (1/h)")
    # Table S1: Keo 0.418 /h (50% RSE); bootstrap median 0.436 (0.371-0.661)

    lkmod <- log(0.564)
    label("Negative-feedback modulator turnover rate constant Kmod (1/h)")
    # Table S1: Kmod 0.564 /h (44% RSE); bootstrap median 0.566 (0.352-0.967)

    lknp <- log(3.34)
    label("First-order nanoparticle degradation rate constant KNP (1/h)")
    # Table S1: KNP 3.34 /h (17% RSE); bootstrap median 3.22 (2.54-4.32); no IIV

    lnp0 <- log(4.7e4)
    label("Initial nanoparticle effect amount NP0 (ng/mL), QLNC groups only")
    # Table S1: NP0 4.7e4 ng/mL (144% RSE); bootstrap median 5.0e4 (6.4e3-8.2e5)

    # Covariate effects. Both are linear on the natural scale about the
    # population median PPI (PPImd = 33.4%), per Results Eqs (6) and (7)
    # and NM-TRAN $PK lines TVBASE / TVSLPE.

    e_score_ppi_rbase_dopamine <- 0.0095
    label("PPI effect on the dopamine baseline DA0 (per % PPI, linear about PPImd)")
    # Table S1: theta_DA0,SPR 0.0095 (20% RSE); bootstrap median 0.0093 (0.0056-0.0123)
    # Discussion restates the whole expression: DA0 = 0.288 * [1 + 0.0095 * (PPI - PPImd)]

    e_score_ppi_eqtp <- 0.0243
    label("PPI effect on the quetiapine effect slope EQTP (per % PPI, linear about PPImd)")
    # Table S1: theta_EQTP,SPR 0.0243 (47% RSE); bootstrap median 0.0265 (0.0077-0.0407)

    # IIV. Results: "IIV was found to be significant on the following model
    # parameters: Kin, EQTP, Ke0, Kmod, DA0, and NP0". Pool0 and KNP carry no
    # IIV (Table S1 shows "-"), matching $OMEGA fixing those etas to zero.
    #
    # Table S1 footnote defines the reported IIV column as
    #   "inter-individual variability expressed as CV (%) calculated as
    #    omega2 x 100, where omega2 is the standard deviation of the variance
    #    of the random effect"
    # i.e. the tabulated percentage is 100 * omega (the SD on the log scale),
    # so the nlmixr2 variance is (CV% / 100)^2. Shrinkage from the same table
    # is noted per line.

    etalkpin ~ 0.36^2            # Table S1: Kin IIV 36% (27% RSE), shrinkage 29%
    etalrbase_dopamine ~ 0.16^2  # Table S1: DA0 IIV 16% (16% RSE), shrinkage 8%
    etaleqtp ~ 0.51^2            # Table S1: EQTP IIV 51% (45% RSE), shrinkage 43%
    etalke0 ~ 0.54^2             # Table S1: Keo IIV 54% (49% RSE), shrinkage 38%
    etalkmod ~ 0.62^2            # Table S1: Kmod IIV 62% (25% RSE), shrinkage 58%
    etalnp0 ~ 1.51^2             # Table S1: NP0 IIV 151% (29% RSE), shrinkage 38%

    # Residual error. NM-TRAN $ERROR is Y = IPRED + W * EPS(1) with
    # W = THETA(5), i.e. a pure additive error whose SD is estimated as a
    # theta; Methods confirm "residual variability was described with an
    # additive error model".
    addSd <- 0.0759
    label("Additive residual error on dopamine (ng/mL)")
    # Table S1: residual error 0.0759 ng/mL (5% RSE), shrinkage 7%;
    # bootstrap median 0.0759 (0.0690-0.0821)
  })

  model({
    # =====================================================================
    # Translation of the NM-TRAN control stream in Supplementary Material
    # S1. NONMEM A(i) compartments map to named states as follows:
    #   A(1) = POOL       -> precursor1    (dopamine precursor pool)
    #   A(2) = DOPAMINE   -> dopamine      (mPFC extracellular dopamine)
    #   A(3) = OUTPUT_DA  -> auc_dopamine  (running integral of dopamine)
    #   A(4) = MODULATOR  -> modulator     (negative-feedback modulator)
    #   A(5) = EFFECT     -> effect        (Sheiner effect compartment)
    #   A(6) = NANO       -> nano          (nanoparticle effect at target)
    # =====================================================================

    # ---- Constants ------------------------------------------------------
    # PPImd, the median prepulse inhibition across the whole animal study
    # population, used to centre both covariate effects.
    ppiMedian <- 33.4  # %; NM-TRAN $PK "PPIMEAN=33.4"; Results text after Eq (7)

    # ---- Individual parameters -------------------------------------------
    rbase_precursor1 <- exp(lrbase_precursor1)          # BASEPOOL (no IIV)
    kpin <- exp(lkpin + etalkpin)                       # RFORM = Kin
    knp <- exp(lknp)                                    # KELNANO (no IIV)

    # Covariate effects are applied to the typical value exactly as in
    # $PK: TVBASE = THETA(3) * (1 + THETA(10) * (PPI - 33.4)) then
    # BASE = TVBASE * EXP(ETA(3)); likewise TVSLPE / SLPE.
    rbase_dopamine <- exp(lrbase_dopamine + etalrbase_dopamine) *
      (1 + e_score_ppi_rbase_dopamine * (SCORE_PPI - ppiMedian))   # BASE = DA0
    eqtp <- exp(leqtp + etaleqtp) *
      (1 + e_score_ppi_eqtp * (SCORE_PPI - ppiMedian))             # SLPE = EQTP

    ke0 <- exp(lke0 + etalke0)                          # KEO
    kmod <- exp(lkmod + etalkmod)                       # KMOD

    # NP0 is gated by the formulation indicator: $PK sets DNANO = 0 unless
    # NANO == 1, so the solution (FQ) groups start with an empty nano state.
    np0 <- exp(lnp0 + etalnp0) * FORM_QTP_QLNC          # DNANO

    # ---- Derived rate constants ------------------------------------------
    # Results: "the parameters Kin and Kout are derived from the expressions
    # Krel x Pool0, and Krel x Pool0 / DA0", i.e. Krel = Kin / Pool0 and
    # Kout = Kin / DA0. NM-TRAN: KBASE = RFORM / BASEPOOL, KEL = RFORM / BASE.
    kpout <- kpin / rbase_precursor1                    # KBASE = Krel (1/h)
    kout <- kpin / rbase_dopamine                       # KEL = Kout (1/h)

    # ---- Initial conditions ----------------------------------------------
    # NM-TRAN: A_0(1) = BASEPOOL, A_0(2) = BASE, A_0(4) = 1, A_0(6) = DNANO.
    # A(3) and A(5) start at zero (NONMEM default). The system is at exact
    # steady state in the absence of drug: with deff = 1 and modulator = 1,
    # d/dt(precursor1) = Kin - (Kin/Pool0)*Pool0 = 0 and
    # d/dt(dopamine) = Kin - (Kin/DA0)*DA0 = 0.
    precursor1(0) <- rbase_precursor1
    dopamine(0) <- rbase_dopamine
    modulator(0) <- 1
    nano(0) <- np0

    # ---- Drug effect ------------------------------------------------------
    # Results: f(CQTP) = EQTP * Cu,brain,e and g(NP) = 1 / (1 + NP), so the
    # bracketed release multiplier of Eqs (1) and (2) is
    #   1 + f(CQTP) * g(NP) = 1 + EQTP * Ce / (1 + NP).
    # NM-TRAN writes this as two branches -- DEFF = 1 + SLPE * A(5) when
    # NANO == 0 and DEFF = 1 + SLPE * A(5) / (1 + ANTAG) when NANO == 1 --
    # which the single expression below reproduces exactly, because nano is
    # identically zero for the FQ groups.
    deff <- 1 + eqtp * effect / (1 + nano)

    # ---- ODE system --------------------------------------------------------
    # Eq (1) / DADT(1)
    d/dt(precursor1) <- kpin - kpout * deff * precursor1
    # Eq (2) / DADT(2)
    d/dt(dopamine) <- kpout * deff * precursor1 - kout * dopamine * modulator
    # DADT(3): running integral of dopamine. NM-TRAN forms the predicted
    # microdialysate concentration as IPRED = A(3) / TIN, the MEAN dopamine
    # concentration over each 20-min collection interval (Methods: "DA
    # microdialysate concentrations were described by the integral over each
    # collection interval; therefore, no assumptions regarding collection
    # times were made"), resetting A(3) each interval. Carrying the running
    # integral as a state lets a user recover exactly that interval mean
    # without reset events, by differencing across collection boundaries:
    #   mean DA over (t1, t2] = (auc_dopamine(t2) - auc_dopamine(t1)) /
    #                           (t2 - t1)
    # The validation vignette does precisely this. The residual error below
    # is applied to the instantaneous dopamine concentration.
    d/dt(auc_dopamine) <- dopamine
    # Eq (4) / DADT(4)
    d/dt(modulator) <- kmod * (dopamine / rbase_dopamine - modulator)
    # DADT(5) = KEO * CONC - A(5): the effect compartment equilibrates
    # towards ke0 * Cu,brain with an implicit first-order rate constant of
    # 1 /h (the bare "- A(5)" term). This is the form the authors ran and it
    # is NOT the printed Eq (3), which reads
    #   dCu,brain,e/dt = Keo * (Cu,brain - Cu,brain,e).
    # The control stream is authoritative here: Figure 4a plots Cu,brain and
    # Cu,brain,e together for all four groups and the plateau ratio
    # Ce/Cu,brain equals Keo in every panel (0.450 -> ~0.40, 0.500 -> ~0.39,
    # 0.568 -> ~0.56, 1.01 -> ~1.02, with QLNC-SPR crossing above Cu,brain),
    # which only this form reproduces; under printed Eq (3) the ratio would
    # approach 1.0 in every panel. Resolved with the operator 2026-08-05.
    # See the vignette Errata for the full argument.
    d/dt(effect) <- ke0 * CU_QTP_BRAIN - effect
    # Eq (5) / DADT(6)
    d/dt(nano) <- -knp * nano

    # ---- Observation and error ---------------------------------------------
    dopamine ~ add(addSd)
  })
}
