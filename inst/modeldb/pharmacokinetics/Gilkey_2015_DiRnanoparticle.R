Gilkey_2015_DiRnanoparticle <- function() {
  description <- "Preclinical (mouse, BALB/c, 4-6 weeks). PBPK model for fluorescently labeled (DiR) block-copolymer nanoparticles in mice, developed as a surrogate model for dexamethasone-encapsulated nanoparticles in pediatric acute lymphoblastic leukemia therapy. Five compartments: plasma, liver, spleen, kidneys, and a virtual 'other' compartment introduced to close the mass balance for ~50% of injected dose that experimental imaging could not account for in the four sampled organs. Single 100 uL IV bolus of 5 ug/mL DiR-NPs; the model treats plasma initial concentration as 5 ug/mL per paper convention (rather than 0.5 ug dose distributed into 1.7 mL plasma volume) -- see vignette Assumptions section."
  reference <- "Gilkey MJ, Krishnan V, Scheetz L, Jia X, Rajasekaran AK, Dhurjati PS. Physiologically based pharmacokinetic modeling of fluorescently labeled block copolymer nanoparticles for controlled drug delivery in leukemia therapy. CPT Pharmacometrics Syst Pharmacol. 2015;4(3):e13. doi:10.1002/psp4.13"
  vignette <- "Gilkey_2015_DiRnanoparticle"
  units <- list(
    time = "min",
    dosing = "ug",
    concentration = "ug/mL"
  )

  covariateData <- list()

  population <- list(
    species        = "mouse (BALB/c, female, 4-6 weeks old)",
    n_subjects     = 3,
    n_studies      = 1,
    age_range      = "4-6 weeks",
    weight_range   = "not reported",
    sex_female_pct = 100,
    disease_state  = "Healthy adult female BALB/c mice (no induced pathology). The fluorescently labeled nanoparticles are a surrogate for dexamethasone-encapsulated block-copolymer nanoparticles being developed for pediatric acute lymphoblastic leukemia (ALL) therapy.",
    dose_range     = "100 uL IV bolus of 5 ug/mL DiR-encapsulated nanoparticles (0.5 ug total mass); single dose, submandibular blood sampling.",
    regions        = "Single-center preclinical (University of Delaware).",
    notes          = "Sample size n = 3 mice per time point; destructively sampled organs (liver, spleen, kidneys; heart, lungs, intestine, gonads, bladder, brain were also harvested but showed no fluorescence). Source data are reused from Krishnan et al. 2013 Mol Pharm 10:2199-2210 (reference 17 of the source paper).",
    scope_note     = "Mechanistic PBPK simulator: no IIV, no residual error -- intended for typical-value simulation. Equations 6-10 of the source paper; parameter values from Table 1. See vignette for steady-state, mass-balance, and paper Figure 1-3 replication checks."
  )

  ini({
    # Fitted parameters (paper Table 1; flagged with footnote 'a' = determined to fit the data).
    # Fixed-physiological parameters (V_P, V_L, V_S, V_K, Q_S, Q_K, Q_O) are declared directly
    # in model() and traced to refs 6, 7 of the source paper (Birnbaum 1994 ILSI; Davies & Morris 1993).

    # Volume of the virtual 'other' compartment (mL). Fitted to close mass balance.
    lvoth <- log(1.01);  label("Volume of 'other' compartment (mL)")                  # Table 1, V_O

    # Liver plasma flow (mL/min). Re-estimated because few compartments were retained
    # and the literature value gave overshoots in initial liver concentration
    # (paper Discussion, "To create better agreement ... Q_L was estimated").
    lql   <- log(0.75);  label("Plasma flow through liver Q_L (mL/min)")              # Table 1, Q_L

    # Kidney clearance (mL/min). Sensitive parameter per paper sensitivity analysis.
    lkk   <- log(2.74);  label("Renal elimination clearance K_K (mL/min)")            # Table 1, K_K

    # Tissue distribution ratios (dimensionless). All fitted.
    lrl   <- log(7.87);  label("Distribution ratio in liver R_L (unitless)")          # Table 1, R_L
    lrs   <- log(1.17);  label("Distribution ratio in spleen R_S (unitless)")         # Table 1, R_S
    lrk   <- log(1.85);  label("Distribution ratio in kidneys R_K (unitless)")        # Table 1, R_K
    lro   <- log(10.13); label("Distribution ratio in 'other' R_O (unitless)")        # Table 1, R_O

    # 'Other'-to-tissue distribution ratios (dimensionless). All fitted.
    lrlo  <- log(14.90); label("Distribution ratio 'other'-to-liver R_LO (unitless)") # Table 1, R_LO
    lrso  <- log(8.43);  label("Distribution ratio 'other'-to-spleen R_SO (unitless)")# Table 1, R_SO
    lrko  <- log(2.50);  label("Distribution ratio 'other'-to-kidneys R_KO (unitless)") # Table 1, R_KO
  })

  model({
    # === Fitted parameters: exp() to linear scale ===
    voth <- exp(lvoth)
    ql   <- exp(lql)
    kk   <- exp(lkk)
    rl   <- exp(lrl)
    rs   <- exp(lrs)
    rk   <- exp(lrk)
    ro   <- exp(lro)
    rlo  <- exp(lrlo)
    rso  <- exp(lrso)
    rko  <- exp(lrko)

    # === Fixed-physiological parameters (paper Table 1, no 'a' footnote) ===
    # Volumes (mL) from refs 6 and 7 (Birnbaum 1994 ILSI; Davies & Morris 1993).
    vp   <- 1.70   # Plasma volume
    vliv <- 1.30   # Liver volume
    vspl <- 0.10   # Spleen volume
    vkid <- 0.34   # Kidneys volume

    # Plasma flow rates (mL/min) through spleen, kidneys, and 'other' compartment.
    qs   <- 0.09   # Spleen plasma flow Q_S
    qk   <- 1.30   # Kidneys plasma flow Q_K
    qo   <- 0.78   # 'Other' plasma flow Q_O

    # === Concentration aliases (used in the ODE right-hand sides) ===
    # States are mass amounts (ug); concentrations (ug/mL) are amount / volume.
    cp   <- plasma / vp
    cliv <- liv    / vliv
    cspl <- spl    / vspl
    ckid <- kid    / vkid
    coth <- oth    / voth

    # === Mass-balance ODEs (paper equations 6 through 10) ===
    # Equations are written in the paper as time derivatives of CONCENTRATION
    # multiplied by 1/V_i. Multiplying both sides by V_i gives the
    # time derivative of AMOUNT, which is what rxode2 integrates by default.

    # Plasma (paper Eq 6)
    d/dt(plasma) <- ql * (cliv / rl) +
                    qs * (cspl / rs) +
                    qk * (ckid / rk) +
                    qo * (coth / ro) -
                    cp * (ql + qs + qk + qo)

    # Liver (paper Eq 7). Net flow from plasma is (Q_L - Q_S) because
    # spleen returns to liver in this venous-mixing arrangement
    # (mesenteric / splanchnic style, although the paper does not explicitly
    # call it portal). Spleen and 'other' deliver to liver via their
    # respective distribution ratios.
    d/dt(liv)    <- cp * (ql - qs) +
                    qs * (cspl / rs) +
                    qo * (coth / rlo) -
                    (cliv / rl) * (ql + qo)

    # Spleen (paper Eq 8)
    d/dt(spl)    <- cp * qs +
                    qo * (coth / rso) -
                    (cspl / rs) * (qs + qo)

    # Kidneys (paper Eq 9). Includes the renal clearance K_K in the outflow.
    d/dt(kid)    <- cp * qk +
                    qo * (coth / rko) -
                    (ckid / rk) * (qk + qo + kk)

    # 'Other' / virtual compartment (paper Eq 10)
    d/dt(oth)    <- qo * (cp - coth / ro)

    # === Initial conditions ===
    # Paper treats t = 0 as a step input with plasma concentration = injection
    # concentration = 5 ug/mL (Results, "at t = 0 all of the injected
    # nanoparticles are contained in the plasma (5 ug/mL) and undergo an
    # exponential decay as expected"). Initial amount in plasma is therefore
    # CP(0) * V_P = 5 * 1.70 = 8.5 ug, not the actually injected mass
    # (100 uL * 5 ug/mL = 0.5 ug). See vignette Assumptions for the
    # mass-balance implication.
    plasma(0) <- 5.0 * 1.70

    # === Observation variables ===
    # Canonical central observation
    Cc       <- cp

    # Per-tissue concentrations as additional outputs for figure replication
    Cliver   <- cliv
    Cspleen  <- cspl
    Ckidneys <- ckid
    Cother   <- coth
  })
}
