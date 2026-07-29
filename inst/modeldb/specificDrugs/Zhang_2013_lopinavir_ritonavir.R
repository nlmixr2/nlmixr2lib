Zhang_2013_lopinavir_ritonavir <- function() {
  description <- paste(
    "Integrated population PK model for lopinavir (1-compartment, first-order",
    "absorption) and ritonavir (2-compartment, transit-chain absorption with",
    "N=2 transit compartments) co-administered to HIV-infected adults (n=21)",
    "and children (n=74; 35 of whom received rifampicin-based antitubercular",
    "treatment). Ritonavir plasma concentration inhibits lopinavir apparent",
    "clearance via a sigmoidal Emax DDI (Emax=0.82, EC50=0.098 mg/L, Hill=2.8).",
    "Rifampicin coadministration increases apparent clearance and reduces",
    "relative bioavailability of both drugs, with separate magnitudes for",
    "adults vs children. Ritonavir dose (mg/kg) drives a linear increase in",
    "relative bioavailability of both drugs; lopinavir-on-adults is the only",
    "dose-effect cell not supported by the data. Diurnal variation is encoded",
    "as a step function with overnight reduction in apparent clearance",
    "(adults 51%, children 27%) and increased bioavailability at the evening",
    "lopinavir dose for adults (+19%). Allometric scaling on apparent CL/Q",
    "(exponent 0.75) and apparent V/Vp (exponent 1) with reference body",
    "weight 65 kg (Zhang 2013)."
  )
  reference <- paste(
    "Zhang C, Denti P, Decloedt EH, Ren Y, Karlsson MO, McIlleron H.",
    "Model-based evaluation of the pharmacokinetic differences between",
    "adults and children for lopinavir and ritonavir in combination with",
    "rifampicin. Br J Clin Pharmacol. 2013;76(5):741-751.",
    "doi:10.1111/bcp.12101."
  )
  vignette <- "Zhang_2013_lopinavir_ritonavir"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling with reference 65 kg (Zhang 2013 Methods: 65 kg is the median body weight in the adult cohort). Exponents fixed at 0.75 for apparent CL and Q, 1 for apparent V and Vp.",
      source_name        = "WT"
    ),
    CHILD = list(
      description        = "Pediatric vs adult indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (adult)",
      notes              = "Zhang 2013 cohort: adults are 21 HIV-infected adults aged 26-58 years (median 36); children are 74 HIV-infected children aged 6 months to 4.5 years (median 21 months). The paper stratifies typical-value CL (lopinavir and ritonavir), ka (lopinavir), MTT (ritonavir), rifampicin effects on CL and F, and the linear ritonavir-dose effect on F by adult vs child without fitting an age-continuous maturation term (the Anderson-Holford maturation function was tested and not supported by the data).",
      source_name        = "CHILD"
    ),
    CONMED_RIF = list(
      description        = "Concomitant rifampicin-based antitubercular treatment",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no rifampicin)",
      notes              = "Zhang 2013 cohort received rifampicin 600 mg daily (adults) or 10 mg/kg daily (children) as part of standard first-line rifampicin-isoniazid-based antituberculosis therapy. Pharmacokinetic sampling occurred after at least 2 weeks of rifampicin co-administration, so the effect is at the chronic post-induction equilibrium. Encoded as a binary on/off flag because the source paper fits a single rifampicin-period multiplicative effect on CL and a paper-cell-specific bioavailability anchor; no time-decaying induction trajectory is modeled (Zhang 2013 Methods).",
      source_name        = "RIF"
    ),
    DOSE_RTV_MGKG = list(
      description        = "Per-administration ritonavir dose per kg body weight",
      units              = "mg/kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Zhang 2013 Methods (paper formula F = BIO * (1 + SLP * (DoseRTV - DoseRTV_STD))): the per-administration ritonavir dose in mg/kg is supplied as a covariate column on every record. The per-population reference dose DoseRTV_STD = 1.5 mg/kg (adult median ritonavir dose without rifampicin co-administration) for adults and 2.9 mg/kg (children's median ritonavir dose without rifampicin co-administration) for children. Cohort-observed values: 1.5 mg/kg (adult standard dose), 2.9 mg/kg (child standard dose), 14 mg/kg (super-boosted children), 6 mg/kg (doubled-dose children).",
      source_name        = "DoseRTV"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 95L,
    n_studies      = 3L,
    n_observations = 1226L,
    age_range      = "6 months-4.5 years (children); 26-58 years (adults)",
    age_median     = "21 months (children); 36 years (adults)",
    weight_range   = "5.0-17.0 kg (children); 43.0-110.0 kg (adults)",
    weight_median  = "10.2 kg (children); 64.5 kg (adults)",
    sex_female_pct = 60,
    race_ethnicity = "Predominantly South African (Western Cape paediatric and adult HIV cohorts); no further race / ethnicity stratification reported in Zhang 2013 Table 1.",
    disease_state  = "HIV infection. 35 of 74 children additionally had active tuberculosis on rifampicin-based antitubercular treatment ('super-boosted' LPV/r with LPV:RTV 1:1 ratio in 15, doubled doses of LPV/r in 20). 21 adults were HIV-infected volunteers without tuberculosis, evaluated at four sequential treatment occasions across LPV/r 400/100 mg without rifampicin, 400/100 mg with rifampicin, 600/150 mg with rifampicin, and 800/200 mg with rifampicin.",
    dose_range     = "Children: LPV/r oral solution 230/57.5 mg/m^2 12-hourly (standard), super-boosted (LPV:RTV 1:1) or doubled 12-hourly doses when on antituberculosis treatment. Median LPV dose 11.6 mg/kg (no TB) or up to ~24 mg/kg (TB doubled). Adults: LPV/r tablet 400/100, 600/150, or 800/200 mg 12-hourly.",
    regions        = "South Africa (Western Cape; Cape Town adult cohort plus paediatric cohorts including Cape Town and Eastern Cape sites).",
    notes          = paste(
      "Integrated joint substrate-and-perpetrator popPK model:",
      "lopinavir is the primary substrate (1-cmt + first-order absorption);",
      "ritonavir is both the booster (sigmoidal Emax inhibition of lopinavir",
      "apparent CL) and a second sampled substrate (2-cmt + N=2 transit-chain",
      "absorption). LLOQ: 0.05 mg/L (LPV), 0.025 mg/L (RTV). Number of",
      "transit compartments NN=2 carried over from the Zhang 2012 adults-only",
      "precursor (doi:10.1111/j.1365-2125.2011.04154.x Table 2 NN=2.03);",
      "Zhang 2013 itself does not list NN. Diurnal variation parameters",
      "(evening effect on CL and F) are encoded in model() via rxode2's",
      "model-time variable t modulo 24 (EVENING_PERIOD = (t %% 24) >= 12),",
      "assuming the simulation event table is aligned to t = 0 = morning dose;",
      "see vignette Assumptions and deviations for the alignment requirement.",
      "Inter-occasion variability (IOV) on bioavailability, CL, ka, and MTT",
      "reported in Zhang 2013 Table 2 is NOT carried in this ini() block",
      "(IIV is retained). Vignette documents this simplification. Bootstrap",
      "n=50 samples; small sample size noted in the paper as a precision",
      "caveat."
    )
  )

  ini({
    # =========================================================================
    # Structural pharmacokinetic parameters -- Zhang 2013 Table 2.
    #
    # Lopinavir is the primary substrate; the canonical names lcl / lvc / lka
    # apply to it. Ritonavir is the second sampled drug AND the perpetrator of
    # the sigmoidal-Emax DDI on lopinavir CL; its parameters carry the _rtv
    # metabolite-suffix.
    #
    # Apparent CL and Q are allometrically scaled with body weight using a
    # fixed exponent of 0.75 (canonical for hepatic and renal clearance);
    # apparent V and Vp use a fixed exponent of 1 (canonical for distribution
    # volume); reference body weight is 65 kg (Zhang 2013 Methods: median
    # adult body weight in the cohort).
    # =========================================================================

    # --- Lopinavir (substrate) ------------------------------------------------
    lcl    <- log(23)
    label("Lopinavir apparent oral CL/F adult typical value (L/h, at 65 kg)")  # Table 2 col 1 LPV adults CL = 23 (20, 26) L/h
    lvc    <- log(57)
    label("Lopinavir apparent V/F (L, at 65 kg; shared adult/child)")          # Table 2 col 2 LPV V = 57 (53, 61) L (single shared value)
    lka    <- log(1.1)
    label("Lopinavir first-order absorption rate ka adult typical value (1/h)") # Table 2 col 1 LPV adults ka = 1.1 (0.9, 1.3) /h

    # Pediatric vs adult shifts on lopinavir CL and ka (log-additive on the
    # log-transformed parameter). Encoded as log(child / adult) so adding to
    # lcl / lka gives the child typical value exactly.
    e_child_cl <- log(15 / 23)
    label("Log shift on lopinavir CL for CHILD = 1 (= log(15/23))")            # Table 2 col 2 LPV children CL = 15 (13, 17) L/h
    e_child_ka <- log(0.38 / 1.1)
    label("Log shift on lopinavir ka for CHILD = 1 (= log(0.38/1.1))")         # Table 2 col 2 LPV children ka = 0.38 (0.30, 0.46) /h

    # --- Ritonavir (perpetrator + second sampled substrate) ------------------
    lcl_rtv    <- log(22)
    label("Ritonavir apparent oral CL/F adult typical value (L/h, at 65 kg)") # Table 2 col 4 RTV adults CL = 22 (19, 25) L/h
    lvc_rtv    <- log(39)
    label("Ritonavir apparent central V/F (L, at 65 kg; shared adult/child)") # Table 2 col 4/5 RTV V = 39 (32, 46) L (single shared value)
    lq_rtv     <- log(30)
    label("Ritonavir apparent inter-compartmental Q/F (L/h, at 65 kg)")       # Table 2 col 4 RTV Q = 30 (26, 34) L/h
    lvp_rtv    <- log(53)
    label("Ritonavir apparent peripheral Vp/F (L, at 65 kg)")                 # Table 2 col 4 RTV Vp = 53 (47, 59) L
    lka_rtv    <- log(2.3)
    label("Ritonavir first-order rate from last transit to central, ka (1/h)") # Table 2 col 4 RTV ka = 2.3 (1.9, 2.7) /h (post-transit rate in Savic chain)
    lmtt_rtv   <- log(1.1)
    label("Ritonavir mean transit time adult typical value (h)")              # Table 2 col 4 RTV MTT = 1.1 (0.94, 1.3) h

    e_child_cl_rtv <- log(13 / 22)
    label("Log shift on ritonavir CL for CHILD = 1 (= log(13/22))")            # Table 2 col 5 RTV children CL = 13 (10, 16) L/h
    e_child_mtt_rtv <- log(2.2 / 1.1)
    label("Log shift on ritonavir MTT for CHILD = 1 (= log(2.2/1.1))")         # Table 2 col 5 RTV children MTT = 2.2 (2.0, 2.4) h

    # Allometric exponents fixed at canonical values per the paper's Methods
    # ("Similar formulas were used for inter-compartmental clearance (Q) and
    # peripheral volume (Vp)" -- allometric scaling on CL/Q with exponent
    # 0.75 and on V/Vp with exponent 1).
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent on apparent CL and Q (unitless, fixed)")        # Zhang 2013 Methods: allometric scaling
    e_wt_vc <- fixed(1.00)
    label("Allometric exponent on apparent V and Vp (unitless, fixed)")        # Zhang 2013 Methods: allometric scaling

    # =========================================================================
    # Rifampicin effect on apparent CL of lopinavir and ritonavir.
    # Table 2 'RIF on CL (+)' row: adult and child shifts differ; encoded as
    # per-population multiplicative factors.
    # =========================================================================
    e_rif_cl_adult     <- 0.58
    label("RIF multiplicative shift on lopinavir CL, adults (unitless)")       # Table 2 col 1 LPV adults RIF on CL = +58% (55%, 61%)
    e_rif_cl_child     <- 0.48
    label("RIF multiplicative shift on lopinavir CL, children (unitless)")     # Table 2 col 2 LPV children RIF on CL = +48% (46%, 50%)
    e_rif_cl_rtv_adult <- 0.34
    label("RIF multiplicative shift on ritonavir CL, adults (unitless)")       # Table 2 col 4 RTV adults RIF on CL = +34% (32%, 36%)
    e_rif_cl_rtv_child <- 0.22
    label("RIF multiplicative shift on ritonavir CL, children (unitless)")     # Table 2 col 5 RTV children RIF on CL = +22% (20%, 24%)

    # =========================================================================
    # Relative bioavailability anchors at the per-population median ritonavir
    # dose (DOSE_RTV_STD), parameterised as log-additive shifts from the
    # adult-no-RIF reference (which is fixed at 1 per the paper's footnote
    # 'Bioavailability has been fixed to 1 for adults without rifampicin co-
    # treatment'). The four cells (adult/child x RIF/no-RIF) are reproduced
    # exactly via the indicator-driven log-additive decomposition. Reported
    # F values come from Zhang 2013 Table 2 'Bioavailability no RIF' /
    # 'Bioavailability with RIF' rows.
    # =========================================================================

    # Lopinavir bioavailability anchors (adult-no-RIF = 1 FIX = reference).
    lfdepot                <- fixed(log(1))
    label("LPV F adult no-RIF reference, FIX = 1 (unitless)")                  # Table 2 col 1 'Bioavailability no RIF' = 1 FIX
    e_child_fdepot         <- log(0.79 / 1)
    label("Log shift on LPV F for CHILD = 1 at no-RIF (= log(0.79))")          # Table 2 col 2 LPV children no-RIF F = 0.79 (0.75, 0.83) (at median RTV dose 2.9 mg/kg)
    e_rif_fdepot           <- log(0.75 / 1)
    label("Log shift on LPV F for CONMED_RIF = 1 in adults (= log(0.75))")     # Table 2 col 1 LPV adults with-RIF F = 0.75 (0.74, 0.76)
    e_child_rif_fdepot     <- log(0.33 / (0.79 * 0.75))
    label("Log interaction shift on LPV F for CHILD * CONMED_RIF (=log(0.33/(0.79*0.75)))") # Table 2 col 2 LPV children with-RIF F = 0.33 (0.31, 0.35); interaction so that exp(0+e_child+e_rif+e_child_rif)=0.33

    # Ritonavir bioavailability anchors (adult-no-RIF = 1 FIX = reference).
    lfdepot_rtv            <- fixed(log(1))
    label("RTV F adult no-RIF reference, FIX = 1 (unitless)")                  # Table 2 col 4 'Bioavailability no RIF' = 1 FIX
    e_child_fdepot_rtv     <- log(0.25 / 1)
    label("Log shift on RTV F for CHILD = 1 at no-RIF (= log(0.25))")          # Table 2 col 5 RTV children no-RIF F = 0.25 (0.24, 0.26) (at median RTV dose 2.9 mg/kg)
    e_rif_fdepot_rtv       <- log(0.48 / 1)
    label("Log shift on RTV F for CONMED_RIF = 1 in adults (= log(0.48))")     # Table 2 col 4 RTV adults with-RIF F = 0.48 (0.47, 0.49)
    e_child_rif_fdepot_rtv <- log(0.021 / (0.25 * 0.48))
    label("Log interaction shift on RTV F for CHILD * CONMED_RIF (=log(0.021/(0.25*0.48)))") # Table 2 col 5 RTV children with-RIF F = 0.021 (0.019, 0.023)

    # =========================================================================
    # Linear ritonavir-dose effect on relative bioavailability of both drugs.
    # Zhang 2013 Methods: F = BIO * (1 + SLP * (DoseRTV - DoseRTV_STD)) with
    # DoseRTV_STD = 1.5 mg/kg (adults) or 2.9 mg/kg (children). The SLP values
    # come from Zhang 2013 Table 2 'Slope of RTV dose effect on bioavailability'
    # row. For lopinavir adults the dose-on-F effect was not supported by the
    # data and is therefore zero in this row of the table; encoded as a
    # CHILD-gated child-only effect on LPV F.
    # =========================================================================
    e_dose_rtv_fdepot_child     <- 0.019
    label("Per-mg/kg RTV dose slope on LPV F, children only (1/(mg/kg))")      # Table 2 col 2 LPV children RTV-dose slope = 0.019 (0.017, 0.021)
    e_dose_rtv_fdepot_rtv_adult <- 0.46
    label("Per-mg/kg RTV dose slope on RTV F, adults (1/(mg/kg))")             # Table 2 col 4 RTV adults RTV-dose slope = 0.46 (0.44, 0.48)
    e_dose_rtv_fdepot_rtv_child <- 0.026
    label("Per-mg/kg RTV dose slope on RTV F, children (1/(mg/kg))")           # Table 2 col 5 RTV children RTV-dose slope = 0.026 (0.024, 0.028)

    # =========================================================================
    # Diurnal (overnight) effect on apparent CL of both drugs.
    # Zhang 2013 Table 2 'Evening effect on CL (-)' row: CL is reduced
    # overnight (post-evening-dose interval) by 51% in adults and 27% in
    # children, for BOTH lopinavir and ritonavir (the two drugs share the
    # diurnal magnitude per the paper's column layout, indicated by the
    # repeated value across the LPV and RTV columns).
    # =========================================================================
    e_evening_cl_adult <- -0.51
    label("Overnight (PM) multiplicative shift on CL, adults (unitless)")      # Table 2 col 1/4 'Evening effect on CL (-)' = 51% (48%, 54%) (both LPV and RTV)
    e_evening_cl_child <- -0.27
    label("Overnight (PM) multiplicative shift on CL, children (unitless)")    # Table 2 col 2/5 'Evening effect on CL (-)' = 27% (25%, 29%) (both LPV and RTV)

    # =========================================================================
    # Diurnal effect on relative bioavailability at the evening lopinavir dose
    # (adults only; not supported for children or for ritonavir in this paper).
    # =========================================================================
    e_evening_fdepot_adult <- 0.19
    label("Evening-dose multiplicative shift on LPV F, adults only (unitless)") # Table 2 col 1 'Evening effect on F (+)' = 19% (18%, 21%); blank in LPV-child / RTV-adult / RTV-child columns

    # =========================================================================
    # Lopinavir-ritonavir sigmoidal-Emax DDI on lopinavir apparent CL.
    # Zhang 2013 Methods equation (paper Methods 'using a sigmoid relationship
    # as reported before [10, 11]'):
    #
    #   CL_LPV(t) = CL0_LPV * (1 - Emax * C_RTV^Hill / (EC50^Hill + C_RTV^Hill))
    #
    # =========================================================================
    emax  <- 0.82
    label("Maximum fractional inhibition of LPV CL by ritonavir (unitless)")   # Table 2 col 3 'Emax' = 0.82 (0.81, 0.83)
    ec50  <- 0.098
    label("RTV concentration producing 50% of Emax on LPV CL (mg/L)")          # Table 2 col 3 'EC50' = 0.098 (0.093, 0.10) mg/L
    hill  <- 2.8
    label("Hill exponent of the LPV CL-vs-RTV inhibition (unitless)")          # Table 2 col 3 'Hill' = 2.8 (2.7, 2.9)

    # =========================================================================
    # Inter-individual variability (IIV). Zhang 2013 Table 2 reports IIV on
    # apparent CL (LPV 22%, RTV 26%), apparent V (LPV 22%; RTV blank, i.e. no
    # IIV term reported), and F (LPV 26%, RTV 65%) with a correlation between
    # LPV and RTV IIV F of 82%. Within-drug correlations are not reported
    # (treated as zero). Inter-occasion variability terms reported in the same
    # table (IOV CL / ka / F / MTT) are NOT carried in this ini() block; the
    # vignette documents this simplification.
    #
    # CV% -> omega^2 conversion: omega^2 = log(1 + CV^2).
    # =========================================================================
    etalcl     ~ log(1 + 0.22^2)
    # Table 2 col 2 'IIV CL (%CV)' LPV = 22% (20%, 24%)
    etalvc     ~ log(1 + 0.22^2)
    # Table 2 col 2 'IIV V (%CV)' LPV = 22% (19%, 25%)
    etalcl_rtv ~ log(1 + 0.26^2)
    # Table 2 col 4 'IIV CL (%CV)' RTV = 26% (24%, 28%)

    # Correlated IIV block on LPV F and RTV F. Covariance is
    # rho * sqrt(var_F_lpv * var_F_rtv).
    etalfdepot + etalfdepot_rtv ~ c(
      log(1 + 0.26^2),
      0.82 * sqrt(log(1 + 0.26^2) * log(1 + 0.65^2)),
      log(1 + 0.65^2)
    )
    # Table 2 col 2 'IIV F (%CV)' LPV = 26% (24%, 28%); col 4 'IIV F (%CV)' RTV = 65% (62%, 68%);
    # 'Correlation between lopinavir and ritonavir IIV F' = 82% (76%, 88%)

    # =========================================================================
    # Residual unexplained variability. Zhang 2013 Table 2 reports a combined
    # additive + proportional structure:
    #   - Additive error LPV = 0.054 mg/L (CI 0.051-0.057). Ritonavir is not
    #     reported with a separate additive component in the table; per the
    #     joint-LPV/RTV typesetting of the additive row (values cluster in
    #     the lopinavir columns only), the additive is implemented on the
    #     lopinavir output only and the ritonavir output uses pure
    #     proportional error (see vignette Assumptions and deviations).
    #   - Proportional error LPV = 13% (12%, 14%); RTV = 21% (19%, 23%).
    # =========================================================================
    addSd      <- 0.054
    label("LPV additive residual error (mg/L)")                                # Table 2 col 1 'Additive error (mg/L)' = 0.054 (0.051, 0.057)
    propSd     <- 0.13
    label("LPV proportional residual error (fraction)")                        # Table 2 col 2 'Proportional error (%)' LPV = 13% (12%, 14%)
    propSd_rtv <- 0.21
    label("RTV proportional residual error (fraction)")                        # Table 2 col 4/5 'Proportional error (%)' RTV = 21% (19%, 23%)
  })

  model({
    # =========================================================================
    # 1. Derived covariate terms.
    #
    # Diurnal (overnight) period indicator. rxode2's `t` is the integration
    # time. Modulo 24 h gives the wall-clock hour with the simulation-time
    # convention that t = 0 corresponds to the morning dose. The overnight
    # period (post-evening-dose interval) is t mod 24 in [12, 24). See the
    # Bienczak 2016 nevirapine model (inst/modeldb/specificDrugs/
    # Bienczak_2016_nevirapine.R, lines 225-228) for the same pattern in a
    # continuous-cosine form. Step-function form here matches Zhang 2013's
    # parameter table (separate per-dose / per-interval shifts).
    # =========================================================================
    evening_period <- (t - 24 * floor(t / 24)) >= 12

    # Per-population standard ritonavir dose for the dose-on-F deviation. The
    # Zhang 2013 Methods identifies 1.5 mg/kg as the adult cohort's median
    # ritonavir dose without rifampicin co-administration and 2.9 mg/kg as
    # the children's median.
    dose_rtv_std <- 1.5 * (1 - CHILD) + 2.9 * CHILD

    # =========================================================================
    # 2. Individual pharmacokinetic parameters.
    #
    # Apparent CL of both drugs picks up the population stratum, allometric
    # body-weight scaling, the per-population rifampicin multiplicative shift,
    # the per-population overnight (diurnal) multiplicative shift, and the
    # individual log-additive IIV. The ritonavir-on-lopinavir sigmoidal Emax
    # inhibition is applied directly to lopinavir CL below.
    # =========================================================================

    # --- Ritonavir (compute first so its central concentration is available
    # for the lopinavir CL inhibition term below) -----------------------------
    cl_rtv_rif_factor <- 1 + e_rif_cl_rtv_adult * (1 - CHILD) * CONMED_RIF +
                              e_rif_cl_rtv_child * CHILD * CONMED_RIF
    cl_diurnal        <- 1 + evening_period *
                              (e_evening_cl_adult * (1 - CHILD) +
                               e_evening_cl_child * CHILD)

    cl_rtv <- exp(lcl_rtv + etalcl_rtv + e_child_cl_rtv * CHILD) *
              (WT / 65)^e_wt_cl *
              cl_rtv_rif_factor *
              cl_diurnal

    vc_rtv <- exp(lvc_rtv) * (WT / 65)^e_wt_vc
    vp_rtv <- exp(lvp_rtv) * (WT / 65)^e_wt_vc
    q_rtv  <- exp(lq_rtv)  * (WT / 65)^e_wt_cl
    ka_rtv <- exp(lka_rtv)

    # Number of transit compartments fixed at NN = 2 (carried over from the
    # adult-only Zhang 2012 precursor doi:10.1111/j.1365-2125.2011.04154.x
    # Table 2 NN = 2.03). The transit-chain rate ktr_rtv = NN / MTT per the
    # Savic 2007 parameterisation.
    mtt_rtv <- exp(lmtt_rtv + e_child_mtt_rtv * CHILD)
    ktr_rtv <- 2 / mtt_rtv

    # Ritonavir plasma concentration used both as the model observation and as
    # the inhibition driver on lopinavir CL.
    crtv <- central_rtv / vc_rtv

    # --- Lopinavir (the substrate; CL is reduced by ritonavir Cp) ------------
    inhib            <- emax * crtv^hill / (ec50^hill + crtv^hill)
    cl_rif_factor    <- 1 + e_rif_cl_adult * (1 - CHILD) * CONMED_RIF +
                             e_rif_cl_child * CHILD * CONMED_RIF

    cl <- exp(lcl + etalcl + e_child_cl * CHILD) *
          (WT / 65)^e_wt_cl *
          cl_rif_factor *
          cl_diurnal *
          (1 - inhib)
    vc <- exp(lvc + etalvc) * (WT / 65)^e_wt_vc
    ka <- exp(lka + e_child_ka * CHILD)

    # Micro-constants
    kel <- cl  / vc
    k12 <- q_rtv / vc_rtv
    k21 <- q_rtv / vp_rtv

    # =========================================================================
    # 3. ODE system.
    #
    # Lopinavir: depot -> central -> elimination (1-compartment).
    # Ritonavir: depot -> transit1_rtv -> transit2_rtv -> central_rtv ->
    #            peripheral1_rtv recycle + elimination (2-compartment with
    #            Savic 2007 N=2 transit-chain absorption).
    # =========================================================================
    d/dt(depot)        <- -ka * depot
    d/dt(central)      <-  ka * depot - kel * central

    d/dt(depot_rtv)        <- -ktr_rtv * depot_rtv
    d/dt(transit1_rtv)     <-  ktr_rtv * depot_rtv      - ktr_rtv * transit1_rtv
    d/dt(transit2_rtv)     <-  ktr_rtv * transit1_rtv   - ka_rtv  * transit2_rtv
    d/dt(central_rtv)      <-  ka_rtv * transit2_rtv    -
                                (cl_rtv / vc_rtv) * central_rtv  -
                                k12 * central_rtv + k21 * peripheral1_rtv
    d/dt(peripheral1_rtv)  <-  k12 * central_rtv - k21 * peripheral1_rtv

    # =========================================================================
    # 4. Relative bioavailability anchors per (CHILD, CONMED_RIF, EVENING)
    # cell, plus the linear ritonavir-dose deviation term.
    #
    # Zhang 2013 Methods (paper formula): F = BIO + SLP * (DoseRTV - DoseRTV_STD)
    # is an ADDITIVE (not multiplicative) deviation on the bioavailability
    # anchor. Verified against the paper's Table 3 (computed for children on
    # the three dose strategies with the median doses actually received):
    #
    #   Super-boosted children (RIF, 14 mg/kg RTV):
    #     LPV F = 0.33 + 0.019 * (14 - 2.9) = 0.541  (Table 3 = 0.53)
    #     RTV F = 0.021 + 0.026 * (14 - 2.9) = 0.310 (Table 3 = 0.31)
    #   Doubled-dose children (RIF, 6 mg/kg RTV):
    #     LPV F = 0.33 + 0.019 * (6 - 2.9) = 0.389   (Table 3 = 0.38)
    #     RTV F = 0.021 + 0.026 * (6 - 2.9) = 0.102  (Table 3 = 0.10)
    #
    # The multiplicative reading would give 0.40 / 0.027 / 0.349 / 0.023, all
    # far off from Table 3. The additive reading matches.
    #
    # Evening-dose effect on F (adult LPV only) is applied as a multiplicative
    # scalar on the composite (base + dose_delta) per the paper's '+19%'
    # description.
    #
    # Log-normal IIV (etalfdepot / etalfdepot_rtv) acts multiplicatively on the
    # composite typical-value F.
    # =========================================================================
    f_lpv_base       <- exp(lfdepot + e_child_fdepot * CHILD +
                            e_rif_fdepot * CONMED_RIF +
                            e_child_rif_fdepot * CHILD * CONMED_RIF)
    f_lpv_dose_delta <- e_dose_rtv_fdepot_child * CHILD *
                              (DOSE_RTV_MGKG - dose_rtv_std)
    f_lpv_evening    <- 1 + e_evening_fdepot_adult * (1 - CHILD) * evening_period
    f_lpv            <- (f_lpv_base + f_lpv_dose_delta) *
                              f_lpv_evening *
                              exp(etalfdepot)

    f_rtv_base       <- exp(lfdepot_rtv + e_child_fdepot_rtv * CHILD +
                            e_rif_fdepot_rtv * CONMED_RIF +
                            e_child_rif_fdepot_rtv * CHILD * CONMED_RIF)
    f_rtv_dose_delta <- (e_dose_rtv_fdepot_rtv_adult * (1 - CHILD) +
                         e_dose_rtv_fdepot_rtv_child * CHILD) *
                              (DOSE_RTV_MGKG - dose_rtv_std)
    f_rtv            <- (f_rtv_base + f_rtv_dose_delta) * exp(etalfdepot_rtv)

    f(depot)     <- f_lpv
    f(depot_rtv) <- f_rtv

    # =========================================================================
    # 5. Observation variables and residual error.
    #
    # Lopinavir plasma concentration Cc uses the combined additive +
    # proportional residual structure; ritonavir plasma concentration Cc_rtv
    # uses pure proportional residual error per Table 2 (the additive row in
    # Table 2 reports a single value clustered with the lopinavir columns).
    # =========================================================================
    Cc      <- central     / vc
    Cc_rtv  <- central_rtv / vc_rtv

    Cc      ~ add(addSd) + prop(propSd)
    Cc_rtv  ~ prop(propSd_rtv)
  })
}
