Henin_2009_capecitabine <- function() {
  description <- "Longitudinal Markov-proportional-odds model for hand-and-foot syndrome (HFS) toxicity grades 0-2 in cancer patients receiving capecitabine. Capecitabine exposure is described with a kinetic-pharmacodynamic (K-PD) one-compartment delay; the per-week effective drug rate drives a sigmoid Emax that shifts the cumulative log-odds for the next HFS grade conditional on the previous grade. Baseline-Cockcroft-Gault creatinine clearance is the only structural covariate."
  reference <- paste(
    "Henin E, You B, VanCutsem E, Hoff PM, Cassidy J, Twelves C,",
    "Zuideveld KP, Sirzen F, Dartois C, Freyer G, Tod M, Girard P (2009).",
    "A dynamic model of hand-and-foot syndrome in patients receiving capecitabine.",
    "Clin Pharmacol Ther 85(4):418-425.",
    "doi:10.1038/clpt.2008.220.",
    "DDMORE Foundation Model Repository: DDMODEL00000214.",
    sep = " "
  )
  vignette <- "Henin_2009_capecitabine"
  ddmore_id <- "DDMODEL00000214"
  replicate_of <- NULL
  units <- list(
    time     = "week",
    dosing   = "mg",
    exposure = "K-PD effective rate K * central (mg/week)",
    outcome  = "HFS grade (0, 1, 2)"
    # units$concentration intentionally omitted: this is a categorical-likelihood
    # K-PD model with no observed plasma concentration. checkModelConventions()
    # reports an info-level note about the missing key; see the vignette's
    # Assumptions and deviations section.
  )

  covariateData <- list(
    CRCL = list(
      description        = "Baseline creatinine clearance (Cockcroft-Gault), not BSA-normalized",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column 'CLCR' (Cockcroft-Gault) maps to the canonical CRCL.",
        "Henin 2009 reports raw Cockcroft-Gault in mL/min (NOT 1.73 m^2-normalized);",
        "value used as-is. Time-fixed at baseline (BCLCR = CLCR at TIME=0",
        "in the source NONMEM .mod, with CLCR constant per subject in the",
        "shipped simulated dataset). Reference value 75.5 mL/min (median",
        "Cockcroft-Gault CrCl in the model-development cohort)."
      ),
      source_name        = "CLCR"
    )
  )

  population <- list(
    n_subjects     = 595,
    n_studies      = 2,
    age_range      = "Adult cancer patients",
    weight_range   = "Adult cancer patients",
    disease_state  = "Patients with metastatic colorectal cancer or advanced/metastatic breast cancer receiving capecitabine in two phase III trials",
    dose_range     = "Capecitabine 1250 mg/m^2 BID for 14 days followed by a 7-day rest (3-week cycles); typical doses 4000-4300 mg/day for adults of ~1.6-1.7 m^2 BSA",
    regions        = "International phase III trials",
    notes          = paste(
      "Demographics summarized from the DDMORE Foundation Model Repository",
      "bundle (DDMODEL00000214) and the linked publication abstract; the",
      "publication PDF is not on disk for this extraction. NOBS = 18,445",
      "HFS-grade observations across 595 subjects (Output_real_HFSmodel.lst).",
      "Per-subject demographics in the shipped Simulated_GHFS_HFSmodel.csv",
      "show baseline CrCl roughly 60-150 mL/min and adult body weights",
      "(~50-100 kg)."
    )
  )

  ini({
    # K-PD compartment (capecitabine effect-rate kinetics)
    lk     <- log(0.102)  ; label("KPD elimination rate constant K (1/week)")                                                          # Output_real_HFSmodel.lst FINAL THETA TH 1 (TVK)

    # Markov-state-conditional cumulative log-odds intercepts for P(grade <= 0)
    b00    <- 4.14        ; label("Cumulative log-odds intercept B00 for grade <= 0 given previous grade 0 (logit units)")             # Output_real_HFSmodel.lst FINAL THETA TH 2 (B00)
    b10    <- 0.855       ; label("Cumulative log-odds intercept B10 for grade <= 0 given previous grade 1 (logit units)")             # Output_real_HFSmodel.lst FINAL THETA TH 3 (B10)
    b20    <- 1.47        ; label("Cumulative log-odds intercept B20 for grade <= 0 given previous grade 2 (logit units)")             # Output_real_HFSmodel.lst FINAL THETA TH 4 (B20)

    # Markov-state-conditional Emax of capecitabine effect on cumulative log-odds
    lemax0 <- log(3.17)   ; label("Maximum log-odds reduction EMAX0 by drug effect, previous grade 0 (logit units)")                   # Output_real_HFSmodel.lst FINAL THETA TH 5 (EMAX0)
    lemax1 <- log(6.65)   ; label("Maximum log-odds reduction EMAX1 by drug effect, previous grade 1 (logit units)")                   # Output_real_HFSmodel.lst FINAL THETA TH 10 (EMAX1)
    lemax2 <- log(8.92)   ; label("Maximum log-odds reduction EMAX2 by drug effect, previous grade 2 (logit units)")                   # Output_real_HFSmodel.lst FINAL THETA TH 11 (EMAX2)

    # Concentration-of-effect (in K*central units, i.e., effective mg/week) at half-Emax
    led50  <- log(12900)  ; label("Effective rate ED50 at half-Emax (mg/week, in K*central units)")                                    # Output_real_HFSmodel.lst FINAL THETA TH 6 (ED50)

    # Markov-state-conditional log-cumulative-logit increments for P(grade <= 1) - P(grade <= 0)
    lb01   <- log(0.626)  ; label("Log of cumulative-logit increment B01 between grade <= 0 and grade <= 1, previous grade 0 (logit units, log scale)")  # Output_real_HFSmodel.lst FINAL THETA TH 7 (B01)
    lb11   <- log(7.24)   ; label("Log of cumulative-logit increment B11 between grade <= 0 and grade <= 1, previous grade 1 (logit units, log scale)")  # Output_real_HFSmodel.lst FINAL THETA TH 8 (B11)
    lb21   <- log(0.330)  ; label("Log of cumulative-logit increment B21 between grade <= 0 and grade <= 1, previous grade 2 (logit units, log scale)")  # Output_real_HFSmodel.lst FINAL THETA TH 9 (B21)

    # Covariate effect: baseline Cockcroft-Gault CrCl on the cumulative-logit intercept (additive)
    e_crcl_b <- 0.00650   ; label("Additive effect of (CRCL - 75.5) on cumulative-logit intercept (logit units per mL/min)")           # Output_real_HFSmodel.lst FINAL THETA TH 12 (TCLCR)

    # Block IIV (NONMEM $OMEGA BLOCK(2) FIXED): etalk on K, etab00 shared additive shift on the cumulative-logit
    # intercept (added to b00, b10, b20 alike). The "etab00" name is associated with b00 to satisfy the eta+typical
    # naming convention; mechanistically the same draw shifts every Markov-state intercept (the NONMEM source
    # statement "IIV=ETA(2)+(BCLCR-75.5)*TCLCR" is added to A0 and A1 regardless of SWM1).
    etalk + etab00 ~ fixed(c(0.802, 0.735, 1.50))                                                                                       # Output_real_HFSmodel.lst FINAL OMEGA BLOCK(2) FIXED: var(ETA1)=0.802, cov=0.735, var(ETA2)=1.50
  })

  model({
    # 1. K-PD compartment: dose into central, first-order elimination at rate k.
    #    central holds the "delayed exposure" amount (mg); central * k is the
    #    instantaneous K-PD effective rate (mg/week) that drives the Emax.
    k <- exp(lk + etalk)
    d/dt(central) <- -k * central

    # 2. Per-Markov-state Emax and cumulative-logit increments (back-transformed from log)
    emax0 <- exp(lemax0)
    emax1 <- exp(lemax1)
    emax2 <- exp(lemax2)
    ed50  <- exp(led50)
    bb01  <- exp(lb01)
    bb11  <- exp(lb11)
    bb21  <- exp(lb21)

    # 3. K-PD effective rate (mg/week)
    effrate <- k * central

    # 4. Per-Markov-state Emax effect (logit units): shrinks the cumulative-logit intercept.
    eff_s0 <- emax0 * effrate / (effrate + ed50)
    eff_s1 <- emax1 * effrate / (effrate + ed50)
    eff_s2 <- emax2 * effrate / (effrate + ed50)

    # 5. Random + covariate intercept shift, added to all three state-conditional log-odds.
    iiv <- etab00 + e_crcl_b * (CRCL - 75.5)

    # 6. Cumulative log-odds A0 (P(grade <= 0)) and A1 (P(grade <= 1)) per Markov state.
    a0_s0 <- b00 - eff_s0 + iiv
    a1_s0 <- a0_s0 + bb01
    a0_s1 <- b10 - eff_s1 + iiv
    a1_s1 <- a0_s1 + bb11
    a0_s2 <- b20 - eff_s2 + iiv
    a1_s2 <- a0_s2 + bb21

    # 7. Cumulative probabilities and per-grade probabilities, exposed per Markov state
    #    (consumer simulates the categorical draw and feeds the realized grade as the
    #    next Markov state). Naming: P<grade>_s<previous>.
    pc0_s0 <- exp(a0_s0) / (1 + exp(a0_s0))
    pc1_s0 <- exp(a1_s0) / (1 + exp(a1_s0))
    P0_s0 <- pc0_s0
    P1_s0 <- pc1_s0 - pc0_s0
    P2_s0 <- 1 - pc1_s0

    pc0_s1 <- exp(a0_s1) / (1 + exp(a0_s1))
    pc1_s1 <- exp(a1_s1) / (1 + exp(a1_s1))
    P0_s1 <- pc0_s1
    P1_s1 <- pc1_s1 - pc0_s1
    P2_s1 <- 1 - pc1_s1

    pc0_s2 <- exp(a0_s2) / (1 + exp(a0_s2))
    pc1_s2 <- exp(a1_s2) / (1 + exp(a1_s2))
    P0_s2 <- pc0_s2
    P1_s2 <- pc1_s2 - pc0_s2
    P2_s2 <- 1 - pc1_s2
  })
}
