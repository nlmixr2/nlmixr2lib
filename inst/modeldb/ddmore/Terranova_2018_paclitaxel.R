Terranova_2018_paclitaxel <- function() {
  description <- "Preclinical xenograft-mouse Dynamic Energy Budget tumor-growth-inhibition (DEB-TGI) PK/PD model with paclitaxel-induced tumor kill and tumor-driven plus drug-driven host cachexia (Terranova 2018; DDMODEL00000274 paclitaxel scenario). Two-compartment paclitaxel PK (rate constants K10/K12/K21 and central volume V1 fixed from upstream popPK) drives a Simeoni-style tumor inhibition arm (proliferating tumor VU1 plus three damaged-cell transit compartments VU2/VU3/VU4) coupled to a host energy budget (structural body component Z, enzyme density EN). Body weight is W = density_V * (1 + xi * EN) * Z; tumor weight is Wu = density_Vu * (VU1 + VU2 + VU3 + VU4). The host-tumor coupling makes the structural-body dynamics piecewise: three switch branches (SWITCH1 / SWITCH2 thresholds with a delta_Vmax cap) determine whether the tumor draws from host enzymes preferentially, from the structural body component, or hits the catabolic body-loss cap."
  reference <- paste(
    "Terranova N., Tosca E. M., Borella E., Pesenti E., Rocchetti M., Magni P. (2018).",
    "Modeling tumor growth inhibition and toxicity outcome after administration of",
    "anticancer agents in xenograft mice: A Dynamic Energy Budget (DEB) approach.",
    "J Theor Biol 450:1-14.",
    "doi:10.1016/j.jtbi.2018.04.012.",
    "DDMORE Foundation Model Repository: DDMODEL00000274 (paclitaxel scenario).",
    sep = " "
  )
  vignette <- "Terranova_2018_paclitaxel"
  paper_specific_compartments <- c("bodyZ", "bodyEn", "tumor1", "tumor2", "tumor3", "tumor4")

  units <- list(
    time          = "day",
    dosing        = "amu",
    concentration = "amu/avu",
    bodyWeight    = "g",
    tumorWeight   = "g",
    notes         = "Paclitaxel dosing and concentration units are not declared by the bundle; amu = arbitrary mass units, avu = arbitrary volume units. Vc, K10, IC50, k2, and the AMT column in Simulated_DEB_TGI_data.csv (3e+07 per dose) all share an implicit unit set carried verbatim from DDMODEL00000274. Body weight and tumor weight are in grams as confirmed by the bundle's Output_simulated_*.pdf axis labels (mouse body weight ~22 g; tumor weight ~0.2-0.5 g)."
  )

  ddmore_id    <- "DDMODEL00000274"
  replicate_of <- NULL

  covariateData <- list()

  population <- list(
    n_subjects     = 2L,
    n_studies      = 1L,
    species        = "mouse (xenograft tumor model)",
    age_range      = NA_character_,
    weight_range   = "approximately 21 g at study start (W_initial = 21.2 g per the .ctl)",
    sex_female_pct = NA_real_,
    disease_state  = "Preclinical xenograft-mouse oncology efficacy / cachexia model. Terranova 2018 develops the Dynamic Energy Budget (DEB) framework as a generalisation of the Simeoni 2004 TGI model that simultaneously predicts tumor growth and host body-weight loss (cachexia). The framework was fit across multiple anticancer agents in xenograft mice; the DDMODEL00000274 bundle implements the paclitaxel scenario only ('Among the drugs considered in the paper, only PACLITAXEL has been used' -- Model_Accomodations.txt).",
    dose_range     = "Paclitaxel intravenous bolus AMT = 3e+07 (source units; consistent with K10/K12/K21 in 1/day and V1 = 813.1) at days 8, 12, 16, 20, 24 in the bundle's Simulated_DEB_TGI_data.csv (treated subject ID = 1; control subject ID = 2 receives no doses).",
    regions        = NA_character_,
    notes          = "DDMORE bundle 274 ships only a single Simulated_DEB_TGI_data.csv with two virtual subjects (one treated, one control) and a Simulated_*.pdf rendering of the typical-value trajectory; no Output_real_*.lst (real-data fit listing) is shipped, and the .ctl is configured for $SIM ONLYSIM with $OMEGA 0 FIX and $SIGMA 1 FIX for both DVID arms. The .ctl $THETA values are therefore the simulation-truth (publication-derived) point estimates rather than re-fitted final estimates; see the Errata section of the validation vignette for the full bundle-versus-publication caveat list. The Terranova 2018 publication (J Theor Biol 450:1-14) was paywalled and not available on disk at extraction time, so the simulation-truth values were not cross-checked against published tables."
  )

  ini({
    # ----------------------------------------------------------------------
    # Parameter values come from the .ctl $THETA block in
    # Executable_Terranova_2017_oncology_TGI_HM.ctl. The bundle ships no
    # Output_real_*.lst (real-data fit) and no Output_simulated_*.lst
    # (NONMEM listing) -- only an Output_simulated_*.pdf rendering of the
    # typical trajectory and a Simulated_DEB_TGI_data.csv event table.
    # The .ctl is configured for $SIM (12345) (54321) ONLYSIM with
    # $OMEGA 0 FIX and $SIGMA 1.0 FIX, so the $THETA values are the
    # publication-derived simulation truth, not refitted final estimates.
    # See ddmore-source.md and the validation vignette Errata for the
    # full caveat list. The Terranova 2018 publication is paywalled and
    # was not available on disk for cross-checking the values against
    # published tables at extraction time.
    # ----------------------------------------------------------------------

    # ----- Estimated DEB-TGI parameters (.ctl $THETA non-FIX entries) -----
    lmu          <- log(0.0223)  ; label("Body-weight reduction rate constant from tumor (mu, 1/day)")                            # .ctl $THETA(1)  = 0.0223  ; mu_POP
    lmu_u        <- log(13.3)    ; label("Cachexia coupling parameter coupling tumor mass to host energy budget (mu_u, unitless)") # .ctl $THETA(2)  = 13.3    ; mu_u_POP
    lgu          <- log(11.7)    ; label("Tumor energy-budget threshold for cachexia onset (gu, unitless)")                       # .ctl $THETA(3)  = 11.7    ; gu_POP
    ldelta_vmax  <- log(0.185)   ; label("Maximum body-weight loss rate cap (delta_Vmax, g/day)")                                 # .ctl $THETA(4)  = 0.185   ; delta_Vmax_POP
    lw_initial   <- log(21.2)    ; label("Initial mouse body weight (g)")                                                         # .ctl $THETA(5)  = 21.2    ; W_initial_POP
    lvu1_initial <- log(0.0023)  ; label("Initial proliferating tumor mass (g)")                                                  # .ctl $THETA(6)  = 0.0023  ; Vu1_initial_POP
    lic50        <- log(0.461)   ; label("Paclitaxel concentration giving 50% inhibition of tumor growth (IC50; same units as central / V1)")  # .ctl $THETA(7)  = 0.461  ; IC50_POP
    lk1          <- log(0.462)   ; label("Damaged-tumor-cell transit rate constant (k1, 1/day)")                                  # .ctl $THETA(8)  = 0.462   ; k1_POP
    lk2          <- log(6.53e-4) ; label("Linear paclitaxel cell-kill coefficient (k2, 1/((conc-unit)*day))")                    # .ctl $THETA(9)  = 6.53e-4 ; k2_POP

    # Residual-error scale coefficients (b_W, b_Wu in the .ctl $ERROR block).
    # In the source the residual structure is
    #     Y = IPRED + b * sqrt(IPRED) * EPS,  EPS ~ N(0, sigma=1)
    # i.e. additive on the linear scale with SD = b * sqrt(IPRED) (Poisson-like).
    # nlmixr2 / rxode2 do not have a direct sqrt-style residual term, so we
    # encode the b coefficients here with plain `add()` residual error and
    # call out the simplification in the vignette Errata. Forward-simulation
    # validation in the vignette uses `rxode2::zeroRe()` and compares typical
    # trajectories against the bundle, so the residual form does not affect
    # the F.2 self-consistency check.
    addSd_bodyWeight  <- 0.101 ; label("Additive residual-error coefficient on body weight (b_W; original form: SD = b_W * sqrt(W))")    # .ctl $THETA(10) = 0.101  ; b_W
    addSd_tumorWeight <- 0.134 ; label("Additive residual-error coefficient on tumor weight (b_Wu; original form: SD = b_Wu * sqrt(Wu))") # .ctl $THETA(11) = 0.134  ; b_Wu

    # ----- Parameters fixed from upstream literature (.ctl $THETA FIX) -----
    lk10        <- fixed(log(20.832))  ; label("Paclitaxel central-elimination rate constant (k10, 1/day) -- FIXED from upstream popPK")               # .ctl $THETA(12) = 20.832  FIX ; K10_POP
    lk12        <- fixed(log(0.144))   ; label("Paclitaxel central-to-peripheral rate constant (k12, 1/day) -- FIXED from upstream popPK")              # .ctl $THETA(13) = 0.144   FIX ; K12_POP
    lk21        <- fixed(log(2.011))   ; label("Paclitaxel peripheral-to-central rate constant (k21, 1/day) -- FIXED from upstream popPK")              # .ctl $THETA(14) = 2.011   FIX ; K21_POP
    lvc         <- fixed(log(813.1))   ; label("Paclitaxel central-compartment volume Vc (source-unit volume; consistent with k10/k12/k21 and AMT) -- FIXED from upstream popPK")  # .ctl $THETA(15) = 813.1   FIX ; V1_POP
    len_initial <- fixed(log(1.3))     ; label("Initial enzyme density EN (DEB host-energy state; unitless) -- FIXED")                                   # .ctl $THETA(16) = 1.3     FIX ; En_initial_POP
    lxi         <- fixed(log(0.184))   ; label("DEB host-body composition coupling xi between EN and structural body component Z (unitless) -- FIXED")  # .ctl $THETA(17) = 0.184   FIX ; xi_POP
    lni         <- fixed(log(1.2242))  ; label("Tumor proliferation rate scaling ni (DEB; unitless) -- FIXED")                                           # .ctl $THETA(18) = 1.2242  FIX ; ni_POP
    lgr         <- fixed(log(12.2))    ; label("Tumor energy-budget rate constant gr (DEB; unitless) -- FIXED")                                          # .ctl $THETA(19) = 12.2    FIX ; gr_POP
    lv1inf      <- fixed(log(22.6))    ; label("Asymptotic tumor-volume scale V1inf (g; DEB) -- FIXED")                                                  # .ctl $THETA(20) = 22.6    FIX ; V1inf_POP
    lrho_b      <- fixed(log(1.0))     ; label("Basal proliferation factor rho_b (set to 1 = unrestricted untreated growth) -- FIXED")                   # .ctl $THETA(21) = 1.0     FIX ; rho_b_POP
  })

  model({
    # ------------------------------------------------------------------
    # 1. Source-fixed mechanistic constants (.ctl $PK lines 91-95).
    #    DENSITY_V, DENSITY_VU set tissue-to-volume conversion (both 1
    #    in the .ctl). OMEG = 0.75 is the DEB allometric exponent.
    # ------------------------------------------------------------------
    density_V  <- 1
    density_Vu <- 1
    omeg       <- 0.75

    # ------------------------------------------------------------------
    # 2. Typical-value parameters (no IIV; .ctl has $OMEGA 0 FIX).
    # ------------------------------------------------------------------
    mu          <- exp(lmu)
    mu_u        <- exp(lmu_u)
    gu          <- exp(lgu)
    delta_vmax  <- exp(ldelta_vmax)
    w_initial   <- exp(lw_initial)
    vu1_initial <- exp(lvu1_initial)
    ic50        <- exp(lic50)
    k1          <- exp(lk1)
    k2          <- exp(lk2)

    k10         <- exp(lk10)
    k12         <- exp(lk12)
    k21         <- exp(lk21)
    vc          <- exp(lvc)
    en_initial  <- exp(len_initial)
    xi          <- exp(lxi)
    ni          <- exp(lni)
    gr          <- exp(lgr)
    v1inf       <- exp(lv1inf)
    rho_b       <- exp(lrho_b)

    # ------------------------------------------------------------------
    # 3. Derived constants used inside the DEB-TGI dynamics.
    #    M is the tumor-volume-normalised proliferation rate (.ctl $PK
    #    L97). z_initial closes the body-weight identity W = (1 + xi*EN)*Z
    #    at t=0 so the user-visible bodyWeight starts at w_initial (.ctl
    #    $PK L99).
    # ------------------------------------------------------------------
    M         <- ni / (v1inf^(1/3) * gr)
    z_initial <- w_initial / (1 + en_initial * xi)

    # ------------------------------------------------------------------
    # 4. Drug-dependent quantities used inside the DES branches.
    #    Cc_pacl = paclitaxel central-compartment concentration (.ctl
    #    $DES L123). rho is the inhibition-modulated proliferation
    #    factor (Hill / Emax form, .ctl L124). ku is the cachexia
    #    coupling fraction (.ctl L125).
    # ------------------------------------------------------------------
    Cc_pacl <- central / vc
    rho     <- rho_b * (1 - Cc_pacl / (ic50 + Cc_pacl))
    ku      <- (mu_u * tumor1) / (bodyZ + mu_u * tumor1)

    # ------------------------------------------------------------------
    # 5. SWITCH thresholds selecting between three structural-body
    #    dynamics regimes (.ctl $DES L126-127).
    # ------------------------------------------------------------------
    switch1 <- ((1 - ku) * ni * bodyEn * bodyZ^(2/3) - gr * M * bodyZ) /
               (gr + (1 - ku) * bodyEn)
    switch2 <- ((1 - ku) * ni * bodyEn * bodyZ^(2/3) - gr * M * bodyZ) /
               ((1 - ku) * (bodyEn + omeg * gr))

    # ------------------------------------------------------------------
    # 6. Piecewise dev_VU1 (proliferating-tumor flux) and dev_Z
    #    (structural-body flux). Three branches, mutually exclusive in
    #    the .ctl IF / IF / IF block (.ctl $DES L129-146):
    #      branch A  (switch1 >= 0)
    #      branch B  (switch1 < 0  AND  switch2 >= -delta_vmax)
    #      branch C  (switch1 < 0  AND  switch2 <= -delta_vmax)
    #    The 1.0e-5 floor in the denominators is from the source as a
    #    numerical guard against divide-by-zero (preserved verbatim).
    #    Math note: it can be shown that switch1 < 0  =>  switch2 < 0
    #    (denom1 > denom2 > 0 with shared numerator), so the three
    #    cases above are exhaustive.
    # ------------------------------------------------------------------
    if (switch1 >= 0) {
      dev_VU1 <- ((ni * bodyZ^(2/3) + M * bodyZ) * gr * ku * bodyEn) /
                 ((gr * gu) + (1 - ku) * gu * bodyEn + 1.0e-5) -
                 mu * tumor1 -
                 k2 * tumor1 * Cc_pacl
      dev_Z   <- ((1 - ku) * ni * bodyEn * bodyZ^(2/3) - gr * M * bodyZ) /
                 (gr + (1 - ku) * bodyEn + 1.0e-5)
    } else if (switch2 >= -delta_vmax) {
      dev_VU1 <- (gr * M * ku * bodyZ) / (gu * (1 - ku) + 1.0e-5) -
                 mu * tumor1 -
                 k2 * Cc_pacl * tumor1
      dev_Z   <- ((1 - ku) * ni * bodyEn * bodyZ^(2/3) - gr * M * bodyZ) /
                 ((1 - ku) * (bodyEn + omeg * gr) + 1.0e-5)
    } else {
      dev_VU1 <- (ku / (gu + 1.0e-5)) *
                 (bodyEn * ni * bodyZ^(2/3) +
                  delta_vmax * bodyEn +
                  delta_vmax * omeg * gr) -
                 mu * tumor1 -
                 k2 * Cc_pacl * tumor1
      dev_Z   <- -delta_vmax
    }

    # ------------------------------------------------------------------
    # 7. ODE system -- direct port of the .ctl $DES DADT lines.
    #    Compartment ordering matches the source so the bundle's
    #    Simulated_DEB_TGI_data.csv (CMT = 1 doses) maps onto
    #    rxode2's positional compartment numbering: central = 1,
    #    peripheral1 = 2, bodyZ = 3, bodyEn = 4, tumor1..tumor4 = 5..8.
    #    Drug PK is the standard 2-compartment first-order-elimination
    #    form (DADT(1) and DADT(3) in the .ctl). The cachexia / TGI
    #    states use the dev_Z / dev_VU1 piecewise terms above.
    # ------------------------------------------------------------------
    d/dt(central)     <- k21 * peripheral1 - (k10 + k12) * central                       # .ctl DADT(1)
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1                               # .ctl DADT(3)
    d/dt(bodyZ)       <- dev_Z                                                           # .ctl DADT(2)
    d/dt(bodyEn)      <- (ni / bodyZ^(1/3)) *                                            # .ctl DADT(4)
                          (rho * (v1inf / (tumor1 + bodyZ))^(2/3) - bodyEn)
    d/dt(tumor1)      <- dev_VU1                                                         # .ctl DADT(5)
    d/dt(tumor2)      <- k2 * Cc_pacl * tumor1 - k1 * tumor2                             # .ctl DADT(6)
    d/dt(tumor3)      <- k1 * tumor2 - k1 * tumor3                                       # .ctl DADT(7)
    d/dt(tumor4)      <- k1 * tumor3 - k1 * tumor4                                       # .ctl DADT(8)

    # ------------------------------------------------------------------
    # 8. Initial conditions (.ctl $PK A_0 lines 103-110 after the
    #    "cut and replace" swap that puts Q1 in slot 1 and Z in slot 2).
    # ------------------------------------------------------------------
    central(0)     <- 0
    peripheral1(0) <- 0
    bodyZ(0)       <- z_initial
    bodyEn(0)      <- en_initial
    tumor1(0)      <- vu1_initial
    tumor2(0)      <- 0
    tumor3(0)      <- 0
    tumor4(0)      <- 0

    # ------------------------------------------------------------------
    # 9. Observations (.ctl $ERROR L182-185):
    #    bodyWeight  = density_V  * (1 + xi * EN) * Z              -> DVID = 1
    #    tumorWeight = density_Vu * (VU1 + VU2 + VU3 + VU4)        -> DVID = 2
    # ------------------------------------------------------------------
    bodyWeight  <- density_V  * (1 + xi * bodyEn) * bodyZ
    tumorWeight <- density_Vu * (tumor1 + tumor2 + tumor3 + tumor4)

    bodyWeight  ~ add(addSd_bodyWeight)
    tumorWeight ~ add(addSd_tumorWeight)
  })
}
