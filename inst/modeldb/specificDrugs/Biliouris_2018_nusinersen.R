Biliouris_2018_nusinersen <- function() {
  description <- paste(
    "Semimechanistic nine-compartment population PK model of nusinersen",
    "(antisense oligonucleotide, Biogen ISIS 396443) intended for extrapolation",
    "of nonhuman-primate (cynomolgus monkey) fits to paediatric patients with",
    "spinal muscular atrophy, following intrathecal lumbar-puncture bolus",
    "administration. Structure (Biliouris 2018): CSF is the sampled dosing site,",
    "coupled bidirectionally to three anatomically distinct spinal-cord segments",
    "(cervical, thoracic, lumbar), to the brain (with a deeper brain-tissue",
    "redistribution compartment), and to the pons; CSF also drains one-way into",
    "the plasma / central compartment, which has its own peripheral distribution",
    "and first-order elimination. Volumes scale linearly with body weight",
    "relative to a 2.8 kg reference; rate constants scale as (WT/2.8)^(-0.08);",
    "the CSF physiological volume V_CSF is overridden by a stepwise",
    "paediatric-age function (120, 130, 135, 140 mL for age classes <0.25,",
    "0.25-0.5, 0.5-1, 1-2 years) rather than by WT. Residual error is",
    "proportional per observed matrix; twelve exponential-IIV terms are carried",
    "on the elimination and transfer rate constants and on the CSF / plasma",
    "volumes. All THETA / OMEGA / SIGMA values are held FIXED from the paper's",
    "final NONMEM simulation control stream on file."
  )
  reference <- paste(
    "Biliouris K, Gaitonde P, Yin W, Norris DA, Wang Y, Henry S, Fey R,",
    "Nestorov I, Schmidt S, Rogge M, Lesko LJ, Trame MN.",
    "A Semi-Mechanistic Population Pharmacokinetic Model of Nusinersen: An",
    "Antisense Oligonucleotide for the Treatment of Spinal Muscular Atrophy.",
    "CPT Pharmacometrics Syst Pharmacol. 2018;7(9):581-592.",
    "doi:10.1002/psp4.12323",
    sep = " "
  )
  vignette <- "Biliouris_2018_nusinersen"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  paper_specific_compartments <- c(
    "spinal_cord_cervical",
    "spinal_cord_lumbar",
    "spinal_cord_thoracic",
    "pons"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    csf                  = list(analyte = "nusinersen", units = "mg", specimen = "CSF", verified = FALSE),
    central              = list(analyte = "nusinersen", units = "mg", specimen = "plasma", verified = FALSE),
    spinal_cord_cervical = list(analyte = "nusinersen", units = "mg", specimen = "administration site", verified = FALSE),
    brain                = list(analyte = "nusinersen", units = "mg", specimen = "tissue", verified = FALSE),
    peripheral1          = list(analyte = "nusinersen", units = "mg", specimen = "plasma", verified = FALSE),
    spinal_cord_lumbar   = list(analyte = "nusinersen", units = "mg", specimen = "administration site", verified = FALSE),
    brain_deep           = list(analyte = "nusinersen", units = "mg", specimen = "tissue", verified = FALSE),
    spinal_cord_thoracic = list(analyte = "nusinersen", units = "mg", specimen = "administration site", verified = FALSE),
    pons                 = list(analyte = "nusinersen", units = "mg", specimen = "administration site", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight (used with LOCF for time-varying body weight in the paediatric extrapolation)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reference weight 2.8 kg per Biliouris 2018 simulation control stream",
        "($PK: `(WT/2.8)`). Volumes are scaled linearly with (WT/2.8) EXCEPT",
        "V_CSF, whose WT-scaling term is commented out in the control stream",
        "(`TVV1 = THETA(1) ; * (WT/2.8)`) so V_CSF is set purely by the",
        "age-based physiological table (see AGE covariate). Rate constants are",
        "scaled as (WT/2.8)^(-0.08). LOCF is applied to WT in the source",
        "dataset (`Description: All St3 IDs with LOCF WT`)."
      ),
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Age in years, used to select the physiological CSF volume V_CSF for paediatric extrapolation",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "V_CSF is set by a stepwise physiological table (Biliouris 2018",
        "simulation control stream $PK): AGE < 0.25 y -> 120 mL; 0.25 <= AGE",
        "< 0.5 -> 130 mL; 0.5 <= AGE < 1 -> 135 mL; 1 <= AGE < 2 -> 140 mL;",
        "otherwise the THETA1 default 150 mL (adult reference). AGE is used",
        "only for this V_CSF selection; it does NOT enter any rate constant."
      ),
      source_name        = "AGE"
    )
  )

  population <- list(
    species        = "cynomolgus monkey (juvenile) with paediatric-human physiological scaling",
    n_subjects     = NA,
    n_studies      = NA,
    age_range      = "monkey fit + paediatric-human simulation window of 0-2 years for the V_CSF table",
    weight_range   = "reference body weight 2.8 kg (monkey / newborn-scale)",
    disease_state  = "Spinal muscular atrophy (SMA) target indication; monkey $PROBLEM tag 'SMA MONKEY'",
    dose_range     = "Intrathecal lumbar-puncture bolus (simulated at 1-12 mg per paediatric-clinic protocol)",
    administration_routes = "Intrathecal bolus into CSF",
    regions        = "Preclinical / translational (cynomolgus-monkey fit extrapolated to paediatric SMA patients)",
    notes          = paste(
      "The on-disk source is the paper's simulation NONMEM control stream",
      "(`$PROBLEM SMA MONKEY`, `$SIMULATION (12345678) ONLYSIM",
      "SUBPROBLEM=1000`) with all $THETA / $OMEGA / $SIGMA values held",
      "FIXED. Compartment layout (paper: `$MODEL COMP=(CSF) COMP=(PLASMA)",
      "COMP=(CSPCORD) COMP=(BRAIN) COMP=(PERIPH1) COMP=(LSPCORD)",
      "COMP=(DTBRAIN) COMP=(TSPCORD) COMP=(PONS)`). Full baseline",
      "demographics for the fitted monkey cohort are not on disk with the",
      "supplementary control stream; refer to Biliouris 2018 main text",
      "Methods and Table 1 for n_subjects and age ranges."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters (Biliouris 2018 simulation NMcode SMA
    # `$THETA` block, all FIX). Compartment indices in the paper are:
    #   1 CSF, 2 PLASMA, 3 CSPCORD (cervical SC), 4 BRAIN,
    #   5 PERIPH1, 6 LSPCORD (lumbar SC), 7 DTBRAIN (deep brain),
    #   8 TSPCORD (thoracic SC), 9 PONS.
    # ------------------------------------------------------------------

    # ---- Physiological volumes (mL) ----
    lv_csf    <- fixed(log(150));   label("CSF adult reference volume V1 (mL; overridden by AGE for < 2 y)")  # $THETA(1) 150 FIX ; (1. V1) mL
    lv_plasma <- fixed(log(937));   label("Plasma volume V2 (mL) at reference weight 2.8 kg")                # $THETA(2) 937 FIX ; (2. V2) mL
    lv_cerv   <- fixed(log(1.91));  label("Cervical spinal-cord volume V3 (mL) at reference weight 2.8 kg")  # $THETA(3) 1.91 FIX ; (3. V3) mL
    lv_brain  <- fixed(log(53.8));  label("Brain volume V4 (mL) at reference weight 2.8 kg")                 # $THETA(4) 53.8 FIX ; (4. V4) mL
    lv_lumb   <- fixed(log(1.08));  label("Lumbar spinal-cord volume V6 (mL) at reference weight 2.8 kg")    # $THETA(20) 1.08 FIX ; (20. V6) mL
    lv_thor   <- fixed(log(1.52));  label("Thoracic spinal-cord volume V8 (mL) at reference weight 2.8 kg")  # $THETA(25) 1.52 FIX ; (25. V8) mL
    lv_pons   <- fixed(log(2.11));  label("Pons volume V9 (mL) at reference weight 2.8 kg")                  # $THETA(28) 2.11 FIX ; (28. V9) mL

    # ---- Rate constants (1/h) ----
    lk_csf_cerv     <- fixed(log(0.00171));  label("CSF -> cervical spinal-cord rate K13 (1/h)")             # $THETA(5) 0.00171 FIX ; (5. K13) 1/h; comment: Rate const. CSF->SPCORD
    lk_cerv_csf     <- fixed(log(0.0001));   label("Cervical spinal-cord -> CSF rate K31 (1/h)")             # $THETA(6) 0.0001 FIX ; (6. K31) 1/h; comment: Rate const. SPCORD->CSF
    lk_csf_brain    <- fixed(log(0.006));    label("CSF -> brain rate K14 (1/h)")                            # $THETA(7) 0.006 FIX ; (7. K14) 1/h; comment: Rate const. CSF->BRAIN
    lk_brain_csf    <- fixed(log(0.0004));   label("Brain -> CSF rate K41 (1/h)")                            # $THETA(8) 0.0004 FIX ; (8. K41) 1/h; comment: Rate const. BRAIN->CSF
    lk_csf_plasma   <- fixed(log(0.0891));   label("CSF -> plasma rate K12 (1/h)")                           # $THETA(9) 0.0891 FIX ; (9. K12) 1/h; comment: Rate const. CSF->PLASMA
    lkel            <- fixed(log(0.206));    label("Plasma elimination rate K20 (1/h)")                      # $THETA(10) 0.206 FIX ; (10. K20) 1/h; comment: Elim rate const PLASMA->out
    lk_plasma_periph <- fixed(log(0.00818)); label("Plasma -> peripheral rate K25 (1/h)")                    # $THETA(11) 0.00818 FIX ; (11. K25) 1/h
    lk_periph_plasma <- fixed(log(0.0001));  label("Peripheral -> plasma rate K52 (1/h)")                    # $THETA(12) 0.0001 FIX ; (12. K52) 1/h
    lk_csf_lumb     <- fixed(log(0.00286));  label("CSF -> lumbar spinal-cord rate K16 (1/h)")               # $THETA(21) 0.00286 FIX ; (21. K16) 1/h
    lk_lumb_csf     <- fixed(log(0.0003));   label("Lumbar spinal-cord -> CSF rate K61 (1/h)")               # $THETA(22) 0.0003 FIX ; (22. K61) 1/h
    lk_brain_deep   <- fixed(log(0.00257));  label("Brain -> deep-brain-tissue rate K47 (1/h)")              # $THETA(23) 0.00257 FIX ; (23. K47) 1/h
    lk_deep_brain   <- fixed(log(0.0001));   label("Deep-brain-tissue -> brain rate K74 (1/h)")              # $THETA(24) 0.0001 FIX ; (24. K74) 1/h
    lk_csf_thor     <- fixed(log(0.0021));   label("CSF -> thoracic spinal-cord rate K18 (1/h)")             # $THETA(26) 0.0021 FIX ; (26.K18) 1/h
    lk_thor_csf     <- fixed(log(0.00045));  label("Thoracic spinal-cord -> CSF rate K81 (1/h)")             # $THETA(27) 0.00045 FIX ; (27.K81) 1/h
    lk_csf_pons     <- fixed(log(0.00157));  label("CSF -> pons rate K19 (1/h)")                             # $THETA(29) 0.00157 FIX ; (29.K19) 1/h
    lk_pons_csf     <- fixed(log(0.0002));   label("Pons -> CSF rate K91 (1/h)")                             # $THETA(30) 0.0002 FIX ; (30.K91) 1/h

    # ---- Body-weight scaling exponent on all rate constants ----
    # From simulation-file description ("scale exp -0.08"); rate constants
    # are multiplied by (WT/2.8)^e_wt_k, e_wt_k = -0.08.
    e_wt_k <- fixed(-0.08); label("Body-weight scaling exponent applied to all rate constants (unitless)")   # Description header: "scale exp -0.08"; $PK: `* 1/((WT/2.8)**0.08)`

    # ------------------------------------------------------------------
    # Inter-individual variability -- $OMEGA FIX block (variances on log
    # scale, exponential IIV `P = TVP * EXP(ETA)`). Mapping from the
    # paper's ETA(i) index to model parameter (from `$PK` block):
    #   ETA(1)  -> K20     ETA(2)  -> K12    ETA(3)  -> V1
    #   ETA(4)  -> V2      ETA(5)  -> K13    ETA(6)  -> K14
    #   ETA(7)  -> K41     ETA(8)  -> K16    ETA(9)  -> K18
    #   ETA(10) -> K81     ETA(11) -> K19    ETA(12) -> K91
    # ETA on K31, K25, K52, V3, V4, V6, K61, K47, K74, V8, V9 is
    # commented out in $PK and therefore not carried in the model.
    # ------------------------------------------------------------------
    etalkel            ~ fixed(0.421)   # $OMEGA 1. -> K20  (ETA(1))
    etalk_csf_plasma   ~ fixed(1.3)     # $OMEGA 2. -> K12  (ETA(2))
    etalv_csf          ~ fixed(0.854)   # $OMEGA 3. -> V1   (ETA(3))
    etalv_plasma       ~ fixed(0.832)   # $OMEGA 4. -> V2   (ETA(4))
    etalk_csf_cerv     ~ fixed(0.453)   # $OMEGA 5. -> K13  (ETA(5))
    etalk_csf_brain    ~ fixed(13.6)    # $OMEGA 6. -> K14  (ETA(6)) -- verbatim from paper; extreme magnitude noted in vignette Errata
    etalk_brain_csf    ~ fixed(0.153)   # $OMEGA 7. -> K41  (ETA(7))
    etalk_csf_lumb     ~ fixed(0.124)   # $OMEGA 8. -> K16  (ETA(8))
    etalk_csf_thor     ~ fixed(0.28)    # $OMEGA 9. -> K18  (ETA(9))
    etalk_thor_csf     ~ fixed(0.129)   # $OMEGA 10. -> K81 (ETA(10))
    etalk_csf_pons     ~ fixed(0.654)   # $OMEGA 11. -> K19 (ETA(11))
    etalk_pons_csf     ~ fixed(0.289)   # $OMEGA 12. -> K91 (ETA(12))

    # ------------------------------------------------------------------
    # Residual error -- proportional per observed matrix.
    # NONMEM formulation: `W = SQRT(THETA(k)^2 * IPRED^2)`, `Y = IPRED +
    # W * ERR(1)` with `$SIGMA 1 FIX`. So `W = |THETA(k)| * |IPRED|` and
    # THETA(k) is directly the proportional-SD fraction on the linear
    # scale. THETA(13..19) map to CSF, plasma, cervical SC, brain,
    # lumbar SC, thoracic SC, pons.
    # ------------------------------------------------------------------
    propSd_Ccsf   <- fixed(0.714);  label("CSF proportional residual SD (fraction)")            # $THETA(13) 0.714 FIX ; (13. ERR1) -- CSF (CMT 1)
    propSd        <- fixed(0.412);  label("Plasma proportional residual SD (fraction)")         # $THETA(14) 0.412 FIX ; (14. ERR2) -- PLASMA (CMT 2)
    propSd_Ccerv  <- fixed(0.101);  label("Cervical spinal-cord proportional residual SD (fraction)")  # $THETA(15) 0.101 FIX ; (15. ERR3) -- CSPCORD (CMT 3)
    propSd_Cbrain <- fixed(1.77);   label("Brain proportional residual SD (fraction)")          # $THETA(16) 1.77 FIX ; (16. ERR4) -- BRAIN (CMT 4)
    propSd_Clumb  <- fixed(0.312);  label("Lumbar spinal-cord proportional residual SD (fraction)")    # $THETA(17) 0.312 FIX ; (17. ERR6) -- LSPCORD (CMT 6)
    propSd_Cthor  <- fixed(0.0239); label("Thoracic spinal-cord proportional residual SD (fraction)")  # $THETA(18) 0.0239 FIX ; (18. ERR8) -- TSPCORD (CMT 8)
    propSd_Cpons  <- fixed(0.0103); label("Pons proportional residual SD (fraction)")           # $THETA(19) 0.0103 FIX ; (19. ERR9) -- PONS (CMT 9)
  })

  model({
    # ---- Weight-scaling factors ----
    # Volumes:  V = TVV * (WT / 2.8)
    # Rate consts: K = TVK * (WT / 2.8)^e_wt_k, with e_wt_k = -0.08.
    wt_ratio <- WT / 2.8
    wt_kfac  <- wt_ratio^e_wt_k

    # ---- V_CSF: age-based physiological override (mL) ----
    # AGE < 0.25         -> 120
    # 0.25 <= AGE < 0.5  -> 130
    # 0.5  <= AGE < 1    -> 135
    # 1    <= AGE < 2    -> 140
    # otherwise           -> exp(lv_csf) = 150 (adult reference).
    # Not WT-scaled: paper's `TVV1 = THETA(1)   ; * (WT/2.8)` -- the
    # (WT/2.8) factor is commented out. Individual variability is applied
    # multiplicatively via etalv_csf.
    tv_v_csf <- ifelse(AGE < 0.25, 120,
                ifelse(AGE < 0.5,  130,
                ifelse(AGE < 1,    135,
                ifelse(AGE < 2,    140, exp(lv_csf)))))
    v_csf <- tv_v_csf * exp(etalv_csf)

    # ---- Other volumes: linearly WT-scaled from 2.8 kg reference ----
    v_plasma <- exp(lv_plasma + etalv_plasma) * wt_ratio
    v_cerv   <- exp(lv_cerv)                  * wt_ratio
    v_brain  <- exp(lv_brain)                 * wt_ratio
    v_lumb   <- exp(lv_lumb)                  * wt_ratio
    v_thor   <- exp(lv_thor)                  * wt_ratio
    v_pons   <- exp(lv_pons)                  * wt_ratio

    # ---- Rate constants: scaled by (WT/2.8)^(-0.08) ----
    k_csf_cerv      <- exp(lk_csf_cerv     + etalk_csf_cerv)     * wt_kfac
    k_cerv_csf      <- exp(lk_cerv_csf)                          * wt_kfac
    k_csf_brain     <- exp(lk_csf_brain    + etalk_csf_brain)    * wt_kfac
    k_brain_csf     <- exp(lk_brain_csf    + etalk_brain_csf)    * wt_kfac
    k_csf_plasma    <- exp(lk_csf_plasma   + etalk_csf_plasma)   * wt_kfac
    kel             <- exp(lkel            + etalkel)            * wt_kfac
    k_plasma_periph <- exp(lk_plasma_periph)                     * wt_kfac
    k_periph_plasma <- exp(lk_periph_plasma)                     * wt_kfac
    k_csf_lumb      <- exp(lk_csf_lumb     + etalk_csf_lumb)     * wt_kfac
    k_lumb_csf      <- exp(lk_lumb_csf)                          * wt_kfac
    k_brain_deep    <- exp(lk_brain_deep)                        * wt_kfac
    k_deep_brain    <- exp(lk_deep_brain)                        * wt_kfac
    k_csf_thor      <- exp(lk_csf_thor     + etalk_csf_thor)     * wt_kfac
    k_thor_csf      <- exp(lk_thor_csf     + etalk_thor_csf)     * wt_kfac
    k_csf_pons      <- exp(lk_csf_pons     + etalk_csf_pons)     * wt_kfac
    k_pons_csf      <- exp(lk_pons_csf     + etalk_pons_csf)     * wt_kfac

    # ---- ODE system (paper $DES). ODE order below matches the paper's
    # compartment numbering 1..9 so the paper's `A(i)` maps positionally
    # to the same slot in rxode2:
    #   1 csf, 2 central (plasma), 3 spinal_cord_cervical, 4 brain,
    #   5 peripheral1, 6 spinal_cord_lumbar, 7 brain_deep,
    #   8 spinal_cord_thoracic, 9 pons.
    d/dt(csf) <- -k_csf_cerv * csf - k_csf_brain * csf - k_csf_plasma * csf -
                  k_csf_lumb * csf - k_csf_thor * csf - k_csf_pons  * csf +
                  k_cerv_csf * spinal_cord_cervical +
                  k_brain_csf * brain +
                  k_lumb_csf * spinal_cord_lumbar +
                  k_thor_csf * spinal_cord_thoracic +
                  k_pons_csf * pons
    d/dt(central) <- k_csf_plasma * csf - kel * central -
                     k_plasma_periph * central + k_periph_plasma * peripheral1
    d/dt(spinal_cord_cervical) <-  k_csf_cerv  * csf - k_cerv_csf * spinal_cord_cervical
    d/dt(brain) <- k_csf_brain * csf - k_brain_csf * brain -
                   k_brain_deep * brain + k_deep_brain * brain_deep
    d/dt(peripheral1) <-  k_plasma_periph * central - k_periph_plasma * peripheral1
    d/dt(spinal_cord_lumbar) <-  k_csf_lumb  * csf - k_lumb_csf * spinal_cord_lumbar
    d/dt(brain_deep) <-  k_brain_deep * brain - k_deep_brain * brain_deep
    d/dt(spinal_cord_thoracic) <-  k_csf_thor  * csf - k_thor_csf * spinal_cord_thoracic
    d/dt(pons) <-  k_csf_pons  * csf - k_pons_csf * pons

    # ---- Bioavailability: doses input in mg become ng-in-state so V (mL)
    # yields concentrations directly in ng/mL. Intrathecal bolus dose
    # into CSF only.
    f(csf) <- 1e6

    # ---- Observations: amount / volume for each sampled tissue ----
    Ccsf   <- csf                    / v_csf
    Cc     <- central                / v_plasma
    Ccerv  <- spinal_cord_cervical   / v_cerv
    Cbrain <- brain                  / v_brain
    Clumb  <- spinal_cord_lumbar     / v_lumb
    Cthor  <- spinal_cord_thoracic   / v_thor
    Cpons  <- pons                   / v_pons

    # ---- Residual error (proportional per output) ----
    Ccsf   ~ prop(propSd_Ccsf)
    Cc     ~ prop(propSd)
    Ccerv  ~ prop(propSd_Ccerv)
    Cbrain ~ prop(propSd_Cbrain)
    Clumb  ~ prop(propSd_Clumb)
    Cthor  ~ prop(propSd_Cthor)
    Cpons  ~ prop(propSd_Cpons)
  })
}
