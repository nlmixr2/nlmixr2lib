Snelder_2020_ASP8232 <- function() {
  description <- "Integrated multiple-target-mediated drug disposition (TMDD) population PK-PD model for the small-molecule vascular adhesion protein-1 (VAP-1 / AOC3) inhibitor ASP8232 in adults (Snelder 2020). First-order oral absorption with lag time, three-compartment disposition (central + two peripheral compartments; V3 fixed equal to V2). ASP8232 binds under a quasi-steady-state assumption to soluble VAP-1 (sVAP-1) in the central compartment and to membrane-bound VAP-1 (mVAP-1) in the central and both peripheral compartments; sVAP-1 turnover and drug-target complex elimination are assumed negligible so all VAP-1 pool sizes are held constant per subject. A single dissociation constant KD applies to every binding site. Three observed analytes: (1) total ASP8232 in the central compartment (free + drug-sVAP-1 complex, matching the LC-MS assay), (2) total sVAP-1 in the central compartment (matching the ELISA assay), and (3) VAP-1 plasma activity (a power function of free sVAP-1). Pooled fit across four studies -- two phase 1 (first-in-human 8232-CL-0001 and renal-impairment / T2DM-CKD 8232-CL-0002) and two phase 2 (VIDI in diabetic macular oedema, ALBUM in diabetic kidney disease). Between-population PK differences are explained by baseline eGFR (CKD-EPI) effects on clearance (sigmoid Emax, Hill coefficient fixed at 10) and relative bioavailability (power). Female subjects have 12.5% higher VAP-1 concentrations. Full omega block on CL, sVAP-1c, and the activity slope SL. Log-additive residual errors with a multiplicative factor of 1.88 on the phase-2 (VIDI + ALBUM) residuals for ASP8232 PK and VAP-1 activity."
  reference <- "Snelder N, Hoefman S, Garcia-Hernandez A, Onkels H, Larsson TE, Bergmann KR. Population pharmacokinetics and pharmacodynamics of a novel vascular adhesion protein-1 inhibitor using a multiple-target mediated drug disposition model. J Pharmacokinet Pharmacodyn. 2021;48(1):39-53. doi:10.1007/s10928-020-09717-w. PMID:32930923."
  vignette <- "Snelder_2020_ASP8232"
  units <- list(
    time          = "h",
    dosing        = "nmol",
    concentration = "nmol/L",
    dosing_notes  = "Amounts are carried internally in nmol so C = A / V yields nmol/L (equivalent to nM), directly comparable to KD, sVAP-1c, and mVAP-1 (all in nM per the paper). Convert an mg dose to nmol by multiplying by 1000/444, i.e. dividing by the ASP8232 free-base molecular weight of 444 g/mol (Snelder 2020 Main modeling assumption 6). Example: a 40 mg oral dose corresponds to 40 * 1000 / 444 = 90.09 nmol."
  )

  covariateData <- list(
    CRCL = list(
      description        = "Baseline estimated glomerular filtration rate (CKD-EPI equation), BSA-normalised.",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject (baseline value). CKD-EPI equation per Snelder 2020 Methods. Two structural covariate effects: (a) sigmoid Emax on CL with EC50 = 77 mL/min/1.73 m^2 and Hill exponent fixed at 10 (nearly a step function around EC50); (b) centred power on relative bioavailability F1 with reference 90 mL/min/1.73 m^2 (see Assumptions and deviations in the vignette; the paper reports the exponent -0.257 without stating the centring value, so a rounded normal-range value of 90 is used).",
      source_name        = "eGFR (CKD-EPI, mL/min/1.73 m^2)"
    ),
    SEXF = list(
      description        = "Female sex indicator (1 = female, 0 = male).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = male",
      notes              = "Snelder 2020 Table 4 reports the covariate as a fractional female-vs-male increase of 0.125 on all VAP-1 concentrations (both sVAP-1 and every mVAP-1 pool). The paper's source column is 'Sex' with 68.9% male in the pooled cohort (Table 3); this file uses the canonical SEXF (female = 1) with the sign flipped as needed.",
      source_name        = "Sex (Male / Female)"
    ),
    STUDY_ASP8232_PHASE2 = list(
      description        = "Phase-2 study cohort indicator (1 = VIDI or ALBUM; 0 = 8232-CL-0001 or 8232-CL-0002).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = phase 1 study",
      notes              = "Selects the log-additive residual-error magnitude on ASP8232 total plasma concentration and on VAP-1 plasma activity between the phase 1 studies (reference) and the phase 2 studies. Snelder 2020 Table 4 estimates a single multiplicative 'Factor res error phase 2 studies' = 1.88 that scales BOTH the ASP8232-PK residual SD and the VAP-1-activity residual SD in the phase 2 subjects. The VAP-1-concentration residual is NOT affected because that assay was run only in the phase 2 studies (VIDI + ALBUM). Set once per subject from the trial identifier; time-fixed.",
      source_name        = "Study (Phase 1 vs Phase 2)"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 363L,
    n_studies        = 4L,
    study_names      = c(
      "8232-CL-0001 (Phase 1 first-in-human single- and multiple-ascending-oral-dose in healthy volunteers; unpublished data)",
      "8232-CL-0002 / NCT02218099 (Phase 1 renal-impairment PK/PD in healthy volunteers, mild/moderate/severe renal impairment, and DKD)",
      "VIDI / 8232-CL-3001 / NCT02302079 (Phase 2 diabetic macular oedema; 40 mg PO QD, +/- 0.3 mg intravitreal ranibizumab)",
      "ALBUM / 8232-CL-0004 / NCT02358096 (Phase 2 diabetic kidney disease; 40 mg PO QD)"
    ),
    age_range        = "adults (individual age ranges not tabulated in the pooled paper).",
    weight_range     = "49.1 to 158 kg (pooled Table 2, N = 363; medians per study 73.6 - 91.7 kg).",
    weight_median    = "83.6 kg (pooled Table 2 median).",
    sex_female_pct   = 31.1,
    sex_notes        = "68.9% male in the pooled cohort (Table 1 and Table 3). Per-study male fraction: 8232-CL-0001 80.4%, 8232-CL-0002 63.6%, VIDI 50.0%, ALBUM 77.5%.",
    egfr_range       = "14.2 to 136 mL/min/1.73 m^2 (pooled Table 2, N = 363; healthy volunteers median 104, ALBUM median 44).",
    egfr_reference   = "90 mL/min/1.73 m^2 (nominal 'normal renal function' assumed for the F1 power-model centering; the paper does not report the centering value explicitly).",
    disease_state    = "Healthy volunteers plus adult patients with renal impairment, type 2 diabetes mellitus with CKD, diabetic kidney disease, and diabetic macular oedema.",
    dose_range       = "0.1-100 mg single oral (8232-CL-0001); 0.2-200 mg QD PO with loading doses of 1-600 mg (8232-CL-0001 multiple-dose); 200 mg single oral or 150 mg QD PO with 250 mg loading (8232-CL-0002); 40 mg PO QD (VIDI, ALBUM). Higher single-dose cohorts of 300-1000 mg (24 subjects) and multiple-dose 800 mg (12 subjects) from 8232-CL-0001 were excluded during model development as outside the considered clinically relevant exposure range.",
    regions          = "Multinational (VIDI + ALBUM were multi-centre; phase 1 study locations not specified in the pooled paper).",
    data_records     = "3498 ASP8232 plasma concentration records; 5893 VAP-1 plasma activity records; 1714 VAP-1 plasma concentration records (Table 1). 1.8%, 2.3%, and 1.4% respectively identified as outliers (|CWRES| > 3) and excluded during model development.",
    notes            = "See Snelder 2020 Tables 1-3 for the study-by-study demographic and baseline distributions. Study 8232-CL-0001 is unpublished (the paper reports demographics only). Molar-unit conversions in the model use the paper's assumed molecular weights of 444 g/mol for ASP8232 (free base) and 84,622 g/mol for the VAP-1 monomer (UniProt Q16853) with mVAP-1 vs sVAP-1 difference assumed negligible."
  )

  ini({
    # ------------------------------------------------------------------
    # ASP8232 disposition -- Snelder 2020 Table 4. All L and L/h are
    # apparent (i.e. per unit F1) because F1 is fixed at 1.
    # ------------------------------------------------------------------
    lka   <- log(3.12);   label("First-order absorption rate ka (1/h)")                        # Snelder 2020 Table 4 (k_a = 3.12 1/h, RSE 12.3%)
    ltlag <- log(0.311);  label("Absorption lag time (h)")                                     # Snelder 2020 Table 4 (LAG = 0.311 h, RSE 8.39%)
    lcl   <- log(17.6);   label("Apparent clearance CL/F1 (L/h) at eGFR << EC50")              # Snelder 2020 Table 4 (CL/F1 = 17.6 L/h, RSE 4.52%); intercept of the sigmoid-Emax model on CL
    lvc   <- log(210);    label("Apparent central volume V2/F1 (L)")                           # Snelder 2020 Table 4 (V2/F1 = 210 L, RSE 7.76%)
    lq    <- log(37.6);   label("Apparent inter-compartmental clearance to peripheral1 Q/F1 (L/h)")  # Snelder 2020 Table 4 (Q/F1 = 37.6 L/h, RSE 12.1%)
    lq2   <- log(80.5);   label("Apparent inter-compartmental clearance to peripheral2 Q2/F1 (L/h)") # Snelder 2020 Table 4 (Q2/F1 = 80.5 L/h, RSE 15.5%)
    lvp2  <- log(26.7);   label("Apparent peripheral2 volume V4/F1 (L)")                       # Snelder 2020 Table 4 (V4/F1 = 26.7 L, RSE 17.4%)
    # V3/V2 = 1 fixed (Table 4 footnote a,b); vp is set = vc inside model()
    # F1 = 1 fixed (apparent-parameter convention); covariate on F1 rescales it
    lfdepot <- fixed(log(1)); label("Reference relative bioavailability F1")     # Snelder 2020: only apparent parameters (CL/F1, V2/F1, ...) are estimated; F1 anchored at 1

    # ------------------------------------------------------------------
    # TMDD binding -- concentrations in nmol/L, KD in nmol/L. Every
    # binding site (sVAP-1c and mVAP-1 in central and both peripheral
    # compartments) shares KD per assumption 4 of the paper's "Main
    # modeling assumptions".
    # ------------------------------------------------------------------
    lkd     <- log(0.929); label("Dissociation constant KD (nmol/L)")                          # Snelder 2020 Table 4 (K_D = 0.929 nM, RSE 8.18%)
    lsvap1c <- log(5.52);  label("Central-compartment soluble VAP-1 total concentration sVAP-1c (nmol/L)")  # Snelder 2020 Table 4 (sVAP-1c = 5.52 nM, RSE 1.97%); NONMEM Bmax in supplement per paper
    f_mvap1c  <- 2.13;               label("Factor for mVAP-1 in central compartment (mVAP-1c = f_mvap1c * sVAP-1c)")           # Snelder 2020 Table 4 (Factor for mVAP-1c = 2.13, RSE 15.9%); paper's Bmax2 / Bmax
    f_mvap1p1 <- 52;                 label("Factor for mVAP-1 in peripheral 1 compartment (mVAP-1p1 = f_mvap1p1 * sVAP-1c)")    # Snelder 2020 Table 4 (Factor for mVAP-1p1 = 52, RSE 11%); paper's Bmax3 / Bmax
    f_mvap1p2 <- fixed(1);           label("Factor for mVAP-1 in peripheral 2 compartment (mVAP-1p2 = f_mvap1p2 * sVAP-1c) []")  # Snelder 2020 Table 4 footnote a,d (Factor for mVAP-1p2 fixed at 1)

    # ------------------------------------------------------------------
    # VAP-1 activity model -- power function on free (unbound) sVAP-1c.
    # ------------------------------------------------------------------
    lsl     <- log(851);   label("Slope for VAP-1 plasma activity power model SL (1/h)")       # Snelder 2020 Table 4 (SL = 851 1/h, RSE 2.99%)
    pow_act <- 0.851;                label("Power exponent for VAP-1 plasma activity power model (unitless)")   # Snelder 2020 Table 4 (POW = 0.851, RSE 1.76%)

    # ------------------------------------------------------------------
    # Covariate effects. All source-traced to Table 4 of Snelder 2020.
    # ------------------------------------------------------------------
    emax_egfr_cl  <- 1.3;                     label("Emax of the sigmoid-Emax effect of baseline eGFR on CL (fractional)")     # Snelder 2020 Table 4 (Emax eGFR on CL/F1 = 1.3, RSE 15.5%)
    lec50_egfr_cl <- log(77);                 label("EC50 of the sigmoid-Emax effect of baseline eGFR on CL (mL/min/1.73 m^2)")  # Snelder 2020 Table 4 (EC50 eGFR on CL/F1 = 77 mL/min/1.73 m^2, RSE 3.95%)
    hill_egfr_cl  <- fixed(10);               label("Hill exponent of the sigmoid-Emax effect of baseline eGFR on CL []")    # Snelder 2020 Table 4 (POW eGFR on CL/F1 = 10 fixed; nearly-step-function around EC50)
    e_egfr_fdepot <- -0.257;                  label("Power exponent of baseline eGFR on F1 (unitless; centred at 90 mL/min/1.73 m^2)")  # Snelder 2020 Table 4 (eGFR on F1 = -0.257, RSE 32.8%)
    e_sexf_vap1   <- 0.125;                   label("Fractional increase in all VAP-1 concentrations for females (unitless)")   # Snelder 2020 Table 4 (SEX on VAP-1 concentrations = 0.125, RSE 28.4%; direction confirmed by paper text 'VAP-1 concentrations were found to be 12.5% higher for females')

    # ------------------------------------------------------------------
    # Study-level residual-error factor for phase-2 studies (VIDI +
    # ALBUM). Multiplies the log-additive SD on ASP8232 PK and on
    # VAP-1 activity. VAP-1 concentration residual is unaffected
    # because sVAP-1 was measured only in the two phase 2 studies.
    # ------------------------------------------------------------------
    factor_res_phase2 <- 1.88;                label("Multiplicative factor on log-additive residual SD for phase 2 studies (unitless)")  # Snelder 2020 Table 4 (Factor res error phase 2 studies = 1.88, RSE 10.4%)

    # ------------------------------------------------------------------
    # Inter-individual variability -- full 3x3 omega block on CL,
    # sVAP-1c, and SL. Variances and covariances read directly from
    # Snelder 2020 Table 4. Order (nlmixr2 lower-triangular c()):
    #   var(CL/F1)       = 0.128    (RSE 15.4%)
    #   cov(CL/F1,sVAP1c)= 0.0213   (RSE 42.4%)
    #   var(sVAP1c)      = 0.0735   (RSE 8.59%)
    #   cov(CL/F1,SL)    = -0.0301  (RSE 29.3%)
    #   cov(sVAP1c,SL)   = -0.0222  (RSE 23.7%)
    #   var(SL)          = 0.0574   (RSE 13.7%)
    # Encoded on the log parameters as etalcl, etalsvap1c, etalsl per
    # the log-normal transform convention.
    # ------------------------------------------------------------------
    etalcl + etalsvap1c + etalsl ~ c(0.128,
                                     0.0213, 0.0735,
                                     -0.0301, -0.0222, 0.0574)

    # ------------------------------------------------------------------
    # Residual error -- log-additive on all three analytes. The paper
    # reports the variances sigma^2 log(PK), sigma^2 log(VAP-1
    # concentration), sigma^2 log(VAP-1 activity); SDs are the square
    # root of these. Encoded via `~ lnorm(sd)` in model() so
    # log(Y_obs) = log(Y_pred) + eps, eps ~ N(0, sd), matching the
    # paper's "log-transformed data were modeled with additive residual
    # error" convention.
    # ------------------------------------------------------------------
    expSd                <- sqrt(0.115);      label("Log-additive residual SD on ASP8232 total plasma concentration Cc, phase 1 studies (unitless)")  # Snelder 2020 Table 4 (sigma^2 log(PK) = 0.115, RSE 8.33%; SD = sqrt of the reported variance)
    expSd_svap1          <- sqrt(0.0351);     label("Log-additive residual SD on total sVAP-1 plasma concentration (unitless; phase 2 studies only)")  # Snelder 2020 Table 4 (sigma^2 log(VAP-1 concentration) = 0.0351, RSE 5.92%; SD = sqrt of the reported variance)
    expSd_vap1activity   <- sqrt(0.0696);     label("Log-additive residual SD on VAP-1 plasma activity, phase 1 studies (unitless)")                   # Snelder 2020 Table 4 (sigma^2 log(VAP-1 activity) = 0.0696, RSE 8.67%; SD = sqrt of the reported variance)
  })

  model({
    # ------------------------------------------------------------------
    # Molecular weights (Snelder 2020 Main modeling assumption 6):
    #   ASP8232 free base : 444 g/mol
    #   VAP-1 monomer     : 84,622 g/mol  (UniProt Q16853; mVAP-1 vs sVAP-1
    #                                       difference assumed negligible)
    # Amounts are carried internally in nmol so C = A / V yields nmol/L
    # (= nM), directly comparable to KD, sVAP-1c, and mVAP-1 in nM.
    # Dose amounts supplied to rxode2::et() are expected in nmol; the
    # vignette converts a mg dose via `nmol = mg * 1000 / 444`.
    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    # 1. Individual PK parameters (log-normal random effects; sigmoid-
    #    Emax and power covariate models on CL and F1).
    # ------------------------------------------------------------------
    ka     <- exp(lka)
    tlag   <- exp(ltlag)
    egfr_ratio_cl <- (CRCL^hill_egfr_cl) / (exp(lec50_egfr_cl)^hill_egfr_cl + CRCL^hill_egfr_cl)
    cl     <- exp(lcl + etalcl) * (1 + emax_egfr_cl * egfr_ratio_cl)
    vc     <- exp(lvc)
    vp     <- vc                                # V3 = V2 (Table 4 footnote a,b)
    q      <- exp(lq)
    q2     <- exp(lq2)
    vp2    <- exp(lvp2)
    fdepot <- exp(lfdepot) * (CRCL / 90)^e_egfr_fdepot

    # ------------------------------------------------------------------
    # 2. Individual VAP-1 pool concentrations (nmol/L). sVAP-1c
    #    carries the only VAP-1 IIV; the mVAP-1 pools are proportional
    #    factors on sVAP-1c and inherit the 12.5% female shift.
    # ------------------------------------------------------------------
    sex_vap1_mult <- 1 + e_sexf_vap1 * SEXF
    svap1c        <- exp(lsvap1c + etalsvap1c) * sex_vap1_mult
    mvap1c        <- svap1c * f_mvap1c
    mvap1p1       <- svap1c * f_mvap1p1
    mvap1p2       <- svap1c * f_mvap1p2
    tvap1c        <- svap1c + mvap1c            # total VAP-1 in central (Snelder 2020 Eq 5 tVAP-1c)
    kd            <- exp(lkd)

    # ------------------------------------------------------------------
    # 3. Micro-constants derived from apparent parameters. All rate
    #    constants act on the free-fraction concentration in the
    #    originating compartment (Snelder 2020 Eqs 6-9).
    # ------------------------------------------------------------------
    kel <- cl / vc
    k23 <- q  / vc
    k32 <- q  / vp
    k24 <- q2 / vc
    k42 <- q2 / vp2

    # ------------------------------------------------------------------
    # 4. Quasi-steady-state free-fractions (Snelder 2020 Eqs 2-5).
    #    General closed form for the QSS binding equilibrium
    #    KD = [D_free] * [T_free] / [complex]:
    #      phi = ([D] - [T] - KD + sqrt(([D] - [T] - KD)^2 + 4*KD*[D])) / (2*[D])
    #    Applied per compartment with:
    #      central       : [D] = asp_c    , [T] = tVAP-1c     -> phi_drug_c
    #      peripheral1   : [D] = asp_p1   , [T] = mVAP-1p1    -> phi_drug_p1
    #      peripheral2   : [D] = asp_p2   , [T] = mVAP-1p2    -> phi_drug_p2
    #    plus the free-fraction of tVAP-1 in central (Eq 5):
    #      [D] = tVAP-1c , [T] = asp_c    -> phi_tvap1_c
    #    A `+ 1e-12` guard on the denominators prevents divide-by-zero
    #    at t = 0 (asp = 0 in the peripheral compartments) without
    #    affecting the mid-simulation values (concentrations are many
    #    orders of magnitude larger during dosing).
    # ------------------------------------------------------------------
    asp_c  <- central     / vc
    asp_p1 <- peripheral1 / vp
    asp_p2 <- peripheral2 / vp2

    disc_dc   <- asp_c  - tvap1c   - kd
    phi_drug_c  <- (disc_dc  + sqrt(disc_dc  * disc_dc  + 4 * kd * asp_c )) / (2 * (asp_c  + 1e-12))

    disc_dp1  <- asp_p1 - mvap1p1  - kd
    phi_drug_p1 <- (disc_dp1 + sqrt(disc_dp1 * disc_dp1 + 4 * kd * asp_p1)) / (2 * (asp_p1 + 1e-12))

    disc_dp2  <- asp_p2 - mvap1p2  - kd
    phi_drug_p2 <- (disc_dp2 + sqrt(disc_dp2 * disc_dp2 + 4 * kd * asp_p2)) / (2 * (asp_p2 + 1e-12))

    disc_tc   <- tvap1c - asp_c    - kd
    phi_tvap1_c <- (disc_tc  + sqrt(disc_tc  * disc_tc  + 4 * kd * tvap1c)) / (2 * (tvap1c + 1e-12))

    # ------------------------------------------------------------------
    # 5. ODE system in nmol. All disposition rates act on the free
    #    drug in the originating compartment.
    # ------------------------------------------------------------------
    rate_elim   <- phi_drug_c  * kel * central
    rate_c_to_p1 <- phi_drug_c  * k23 * central
    rate_p1_to_c <- phi_drug_p1 * k32 * peripheral1
    rate_c_to_p2 <- phi_drug_c  * k24 * central
    rate_p2_to_c <- phi_drug_p2 * k42 * peripheral2

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <- ka * depot - rate_elim - rate_c_to_p1 + rate_p1_to_c - rate_c_to_p2 + rate_p2_to_c
    d/dt(peripheral1) <- rate_c_to_p1 - rate_p1_to_c
    d/dt(peripheral2) <- rate_c_to_p2 - rate_p2_to_c

    # Absorption lag time and relative bioavailability. Dose in nmol.
    alag(depot) <- tlag
    f(depot)    <- fdepot

    # ------------------------------------------------------------------
    # 6. Observation predictions (Snelder 2020 Eqs 10-12; linear scale;
    #    log-normal residual error applied via `~ lnorm(sd)`).
    #
    # ASP8232_total_c (nmol/L) : Eq 10 -- LC-MS measures free + drug-
    #                                     sVAP1c complex only (not
    #                                     drug-mVAP-1c complex, which
    #                                     is membrane-bound and outside
    #                                     the plasma matrix).
    # sVAP1_total_c   (nmol/L) : Eq 11 -- ELISA measures total sVAP-1
    #                                     (free + drug-bound), which
    #                                     equals sVAP-1c since VAP-1
    #                                     turnover is neglected.
    # VAP1_activity   (arb.  ) : Eq 12 -- power function of free
    #                                     sVAP-1 concentration.
    # ------------------------------------------------------------------
    sl_i         <- exp(lsl + etalsl)
    Cc           <- phi_drug_c * asp_c + (1 - phi_drug_c) * asp_c * svap1c / tvap1c
    svap1        <- svap1c
    vap1activity <- sl_i * (phi_tvap1_c * svap1c)^pow_act

    # Study-dependent residual SD on ASP8232 PK and VAP-1 activity
    # (multiplicative factor for phase 2). VAP-1 concentration is
    # measured only in phase 2 studies; its SD is not scaled.
    sd_cc           <- expSd              * (1 + STUDY_ASP8232_PHASE2 * (factor_res_phase2 - 1))
    sd_vap1activity <- expSd_vap1activity * (1 + STUDY_ASP8232_PHASE2 * (factor_res_phase2 - 1))

    Cc           ~ lnorm(sd_cc)
    svap1        ~ lnorm(expSd_svap1)
    vap1activity ~ lnorm(sd_vap1activity)
  })
}
