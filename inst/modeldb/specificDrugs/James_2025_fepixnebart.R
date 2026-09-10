James_2025_fepixnebart <- function() {
  description <- "Simultaneous population PK/PD model for fepixnebart (LY3016859, a humanized IgG4 monoclonal antibody against epiregulin and TGF-alpha) and its soluble target epiregulin in adults with chronic pain (James 2025): two-compartment IV disposition with parallel linear and Michaelis-Menten elimination, an indirect-response epiregulin turnover compartment in which fepixnebart inhibits epiregulin degradation through a sigmoid Emax (Emax fixed to 1) relationship, estimated allometric weight exponents on CL, Q, Vc and Vp, sex and glomerular filtration rate on CL, and sex and pain indication on Vc. The drug-effect fraction doubles as the predicted soluble target engagement."
  reference <- "James DE, Bailey J, van der Walt J-S, Winkler J, Schoemaker R. Population pharmacokinetics and pharmacodynamics of fepixnebart (LY3016859) and epiregulin in patients with chronic pain. Clin Pharmacokinet. 2025;64(5):757-766. doi:10.1007/s40262-025-01506-3"
  vignette <- "James_2025_fepixnebart"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # `epiregulin` is the soluble-target turnover state of the indirect-response
  # arm (James 2025 Eq. 3). It is declared paper-specific rather than promoted
  # to a canonical compartment: a single paper is not enough to found a new
  # canonical compartment name, and a second soluble-EGFR-ligand turnover model
  # is the trigger to revisit. The existing canonicals do not fit -- `target`
  # is the FREE target of a TMDD/receptor-binding model and `total_target` is
  # specifically the quasi-steady-state TMDD total-target state, whereas this
  # state is an indirect-response biomarker pool carried in concentration units.
  paper_specific_compartments <- c("epiregulin")

  compartmentData <- list(
    central     = list(analyte = "fepixnebart", units = "mg",    specimen = "serum",          verified = TRUE),
    peripheral1 = list(analyte = "fepixnebart", units = "mg",    specimen = "not applicable", verified = TRUE),
    epiregulin  = list(analyte = "epiregulin",  units = "pg/mL", specimen = "serum",          verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on CL, Q, Vc and Vp with a reference individual of 70 kg (James 2025 Table 3 footnote b and Sect. 2.4). All four exponents were ESTIMATED rather than fixed at their theoretical 0.75 / 1 values: estimating the Vc and Vp exponents dropped the objective function by 31.18 points (Sect. 3.3). Cohort mean 90.4 kg, range 47-148 kg (Table 1).",
      source_name        = "Body weight"
    ),
    SEXF = list(
      description        = "Sex: 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Female sex raises both CL and Vc. James 2025 Table 3 footnotes c and d give the log-scale increments directly: CLSEX = 0.153 and VCSEX = 0.162 for female patients, 0 for male patients; the Table 3 body reports the corresponding fold changes exp(0.153) = 1.17 on CL and exp(0.162) = 1.18 on Vc. Cohort 53.1% female (Table 1).",
      source_name        = "SEX"
    ),
    CRCL = list(
      description        = "Body-surface-area-adjusted glomerular filtration rate estimated with the six-variable Modification of Diet in Renal Disease (MDRD-6) study equation",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "MDRD-6 estimated GFR, BSA-normalized (James 2025 Sect. 2.5). Enters the log-scale CL as 0.21 * log(CRCL / 88) / (log(119) - log(65)), i.e. a log-linear term rescaled so the coefficient spans the 95% percentile range of the observed GFR: 88 is the cohort median and 65-119 mL/min/1.73 m^2 the 95% range (Table 3 footnote c). The Table 3 body reports the resulting fold change across that range, 1.23. Cohort mean 88.6, range 52.0-138 mL/min/1.73 m^2 (Table 1).",
      source_name        = "GFR"
    ),
    DIS_DPN = list(
      description        = "Pain indication: 1 = painful diabetic peripheral neuropathic pain (DPNP), 0 = other indication",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (chronic low back pain, the reference indication)",
      notes              = "One of the two non-reference levels of the three-level pain indication in James 2025 (CLBP reference, DPNP, OA). Table 3 footnote d: VCSTUDY = -0.102 for DPNP, additive on the log scale, giving the tabulated fold change exp(-0.102) = 0.903 on Vc. n = 124 of 386 (Table 1). Paired with DIS_OA; both are 0 for the CLBP reference group.",
      source_name        = "Disease (study)"
    ),
    DIS_OA = list(
      description        = "Pain indication: 1 = osteoarthritis knee pain, 0 = other indication",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (chronic low back pain, the reference indication)",
      notes              = "The second non-reference level of the three-level pain indication. Table 3 footnote d: VCSTUDY = 0.00735 for OA, giving the tabulated fold change exp(0.00735) = 1.01 on Vc -- an essentially null effect whose 408.2% RSE the paper reports without comment; it is retained because the paper retained the whole three-level term. n = 113 of 386 (Table 1). Paired with DIS_DPN.",
      source_name        = "Disease (study)"
    )
  )

  # Covariates the paper SCREENED but did not retain in the final model. They
  # are recorded here for provenance only and are deliberately absent from
  # model(): James 2025 Table 2 reports the single-covariate objective-function
  # drops, none of which met the p < 0.001 backward-deletion criterion.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline",
      units       = "years",
      type        = "continuous",
      notes       = "Screened on CL (dOFV -3.33, p = 0.0680), EC50 (dOFV -0.48, p = 0.4888) and Vc (dOFV -2.54, p = 0.1108); not retained (James 2025 Table 2, runs 301 / 307 / 312). Cohort mean 59.4 years, range 20-84."
    ),
    ADA_POS = list(
      description = "Time-varying treatment-emergent anti-drug antibody status: 1 = ADA-positive at the current time, 0 = ADA-negative",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a TIME-VARYING covariate on CL (dOFV -2.73, p = 0.0983) and Vc (dOFV -0.02, p = 0.9004); not retained (James 2025 Table 2, runs 304 / 315). 70 of 255 fepixnebart-treated participants developed treatment-emergent ADAs, 69 of them neutralizing, highest titer 1:2560 (Sect. 3.2) -- so the negative finding is not for want of ADA-positive subjects."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 386L,
    n_studies      = 3L,
    age_range      = "20-84 years",
    age_mean       = "59.4 years",
    weight_range   = "47-148 kg",
    weight_mean    = "90.4 kg",
    sex_female_pct = 53.1,
    disease_state  = "Chronic pain: chronic low back pain (n = 149), painful diabetic peripheral neuropathy (n = 124), and osteoarthritis knee pain (n = 113).",
    renal_function = "MDRD-6 estimated GFR mean 88.6 mL/min/1.73 m^2, range 52.0-138 (Table 1).",
    dose_range     = "750 mg intravenous loading dose followed by three 500 mg intravenous doses every 2 weeks (4 infusions in total), each infused over 1 h.",
    n_observations_drug       = 2444L,
    n_observations_epiregulin = 2436L,
    notes          = "Pooled from three 26-week phase 2 proof-of-concept, randomized, double-blind, placebo-controlled studies (NCT04456686 osteoarthritis, NCT04476108 DPNP, NCT04529096 CLBP), each with an 8-week double-blind treatment period and an 18-week follow-up. Participants were randomized 2:1 to fepixnebart or placebo, so the 386-subject dataset includes placebo recipients who contribute epiregulin baseline data; 255 received fepixnebart. Demographics from James 2025 Table 1 (mean (range))."
  )

  ini({
    # ---- Structural PK (James 2025 Table 3, "Estimate" column) ------------
    # The paper tabulates the clearances in mL/h; they are recorded here in
    # L/h so that they are dimensionally consistent with the volumes (L),
    # the mg dosing units and the mg/L concentration units.
    lcl   <- log(6.72 / 1000);  label("Linear clearance (L/h) for a 70 kg male at the median GFR")   # James 2025 Table 3: CL = 6.72 mL/h
    lvc   <- log(2.42);         label("Central volume of distribution (L) for a 70 kg male with CLBP") # James 2025 Table 3: Vc = 2.42 L
    lq    <- log(12.7 / 1000);  label("Intercompartmental clearance (L/h) for a 70 kg individual")   # James 2025 Table 3: Q  = 12.7 mL/h
    lvp   <- log(2.06);         label("Peripheral volume of distribution (L) for a 70 kg individual") # James 2025 Table 3: Vp = 2.06 L

    # Michaelis-Menten (target-mediated) elimination from the central
    # compartment. Vmax is tabulated in ug/h and recorded here in mg/h.
    lvmax <- log(41.5 / 1000);  label("Maximum Michaelis-Menten elimination rate (mg/h)")            # James 2025 Table 3: Vmax = 41.5 ug/h
    lkm   <- log(0.966);        label("Fepixnebart concentration at half-maximal Vmax (mg/L)")       # James 2025 Table 3: Km = 0.966 mg/L

    # ---- Epiregulin indirect-response model (James 2025 Eq. 3, Table 3) ---
    lec50  <- log(3.42);        label("Fepixnebart concentration at 50% effect on epiregulin (mg/L)") # James 2025 Table 3: EC50 = 3.42 mg/L
    lrbase <- log(278);         label("Baseline soluble epiregulin concentration (pg/mL)")            # James 2025 Table 3: Baseline epiregulin = 278 pg/mL
    lkdeg  <- log(0.0234);      label("Epiregulin degradation rate in the absence of fepixnebart (1/h)") # James 2025 Table 3: Kdeg = 0.0234 /h
    lhill  <- log(0.723);       label("Hill coefficient of the fepixnebart effect on epiregulin (unitless)") # James 2025 Table 3: Hill factor = 0.723
    # "Preliminary assessment suggested that suppression was very pronounced,
    # and the maximal effect of fepixnebart on epiregulin (Emax) could
    # therefore be fixed to one" (James 2025 Sect. 2.4).
    lemax  <- fixed(log(1));    label("Maximal fractional inhibition of epiregulin degradation (unitless)") # James 2025 Sect. 2.4: Emax fixed to 1

    # ---- Allometric exponents on body weight (all ESTIMATED) --------------
    # Table 3 footnote b; Sect. 2.4 describes the a-priori allometric form
    # PAR_i = theta1 * (WT_i / 70)^theta2 * exp(eta_i), and Sect. 3.3 records
    # that the exponents were estimated rather than fixed at 0.75 / 1.
    e_wt_cl <- 1.06;            label("Allometric exponent of body weight on CL (unitless)")   # James 2025 Table 3: exponent of WT on CL = 1.06
    e_wt_q  <- 2.81;            label("Allometric exponent of body weight on Q (unitless)")    # James 2025 Table 3: exponent of WT on Q  = 2.81
    e_wt_vc <- 0.623;           label("Allometric exponent of body weight on Vc (unitless)")   # James 2025 Table 3: exponent of WT on Vc = 0.623
    e_wt_vp <- 1.06;            label("Allometric exponent of body weight on Vp (unitless)")   # James 2025 Table 3: exponent of WT on Vp = 1.06

    # ---- Covariate effects, all additive on the log scale -----------------
    # Table 3 footnote c writes CL out in full; the 0.21 coefficient is
    # divided by (log(119) - log(65)) inside model() so that it expresses the
    # effect across the 95% percentile range of GFR rather than per log unit.
    e_crcl_cl <- 0.21;          label("Effect of GFR on log CL across the 95% GFR percentile range (unitless)") # James 2025 Table 3 footnote c
    e_sexf_cl <- 0.153;         label("Effect of female sex on log CL (CLSEX; unitless)")      # James 2025 Table 3 footnote c: CLSEX = 0.153 for female
    e_sexf_vc <- 0.162;         label("Effect of female sex on log Vc (VCSEX; unitless)")      # James 2025 Table 3 footnote d: VCSEX = 0.162 for female
    e_dpn_vc  <- -0.102;        label("Effect of DPNP indication on log Vc (VCSTUDY; unitless)") # James 2025 Table 3 footnote d: VCSTUDY = -0.102 for DPNP
    e_oa_vc   <- 0.00735;       label("Effect of osteoarthritis indication on log Vc (VCSTUDY; unitless)") # James 2025 Table 3 footnote d: VCSTUDY = 0.00735 for OA

    # ---- IIV -------------------------------------------------------------
    # James 2025 Table 3 reports interindividual variability as a percentage
    # in the "IIV" column, read here as the coefficient of variation of a
    # log-normal random effect, so omega^2 = log(1 + CV^2). At these
    # magnitudes the alternative reading (the percentage IS the log-scale SD)
    # changes each omega by at most 0.02 -- see the vignette Errata. The
    # paper reports no OMEGA off-diagonals, so the block is diagonal, and
    # "IIV was estimated on all parameters except for Q" (Sect. 3.3).
    etalcl    ~ log(1 + 0.212^2); label("IIV on linear clearance")                      # James 2025 Table 3: IIV on CL = 21.2%
    etalvc    ~ log(1 + 0.178^2); label("IIV on central volume")                        # James 2025 Table 3: IIV on Vc = 17.8%
    etalvp    ~ log(1 + 0.164^2); label("IIV on peripheral volume")                     # James 2025 Table 3: IIV on Vp = 16.4%
    etalvmax  ~ log(1 + 0.270^2); label("IIV on maximum Michaelis-Menten elimination rate") # James 2025 Table 3: IIV on Vmax = 27.0%
    etalkm    ~ log(1 + 0.266^2); label("IIV on Michaelis-Menten constant")             # James 2025 Table 3: IIV on Km = 26.6%
    etalec50  ~ log(1 + 0.195^2); label("IIV on EC50")                                  # James 2025 Table 3: IIV on EC50 = 19.5%
    etalrbase ~ log(1 + 0.249^2); label("IIV on baseline epiregulin")                   # James 2025 Table 3: IIV on baseline epiregulin = 24.9%
    etalkdeg  ~ log(1 + 0.262^2); label("IIV on epiregulin degradation rate")           # James 2025 Table 3: IIV on Kdeg = 26.2%
    etalhill  ~ log(1 + 0.143^2); label("IIV on the Hill coefficient")                  # James 2025 Table 3: IIV on Hill factor = 14.3%

    # ---- Residual unexplained variability ---------------------------------
    # Table 3 legend: "RUV proportional residual unexplained variability
    # (residual error)"; both endpoints are proportional-only.
    propSd     <- 0.138; label("Proportional residual error on fepixnebart concentration (fraction)") # James 2025 Table 3: RUV LY3016859 = 13.8%
    propSd_Epi <- 0.181; label("Proportional residual error on epiregulin concentration (fraction)")  # James 2025 Table 3: RUV Epiregulin = 18.1%
  })

  model({
    # ---- 1. Covariate terms, all additive on the log scale ----------------
    # James 2025 Table 3 footnote c:
    #   CL = exp(log(6.72) + 0.21 * log(GFR/88)/(log(119) - log(65))
    #            + CLSEX + 1.06 * log(WT/70))
    # 88 is the median GFR and 65-119 mL/min/1.73 m^2 the 95% percentile
    # range of the GFR data, so the 0.21 coefficient is the log fold change
    # across that range rather than per log unit of GFR.
    cov_cl <- e_crcl_cl * log(CRCL / 88) / (log(119) - log(65)) +
              e_sexf_cl * SEXF
    # Table 3 footnote d:
    #   Vc = exp(log(2.42) + VCSEX + VCSTUDY + 0.623 * log(WT/70))
    # with VCSTUDY = 0 / -0.102 / 0.00735 for CLBP / DPNP / OA, so chronic
    # low back pain is the reference indication and contributes nothing.
    cov_vc <- e_sexf_vc * SEXF + e_dpn_vc * DIS_DPN + e_oa_vc * DIS_OA

    # ---- 2. Individual parameters ----------------------------------------
    cl    <- exp(lcl + etalcl + cov_cl) * (WT / 70)^e_wt_cl   # L/h
    vc    <- exp(lvc + etalvc + cov_vc) * (WT / 70)^e_wt_vc   # L
    q     <- exp(lq)                    * (WT / 70)^e_wt_q    # L/h; no IIV (Sect. 3.3)
    vp    <- exp(lvp + etalvp)          * (WT / 70)^e_wt_vp   # L
    vmax  <- exp(lvmax + etalvmax)                            # mg/h
    km    <- exp(lkm   + etalkm)                              # mg/L
    ec50  <- exp(lec50 + etalec50)                            # mg/L
    hill  <- exp(lhill + etalhill)                            # unitless
    emax  <- exp(lemax)                                       # unitless, fixed at 1
    kdeg  <- exp(lkdeg + etalkdeg)                            # 1/h
    rbase <- exp(lrbase + etalrbase)                          # pg/mL
    # Epiregulin turnover starts at steady state, so Ksyn = baseline * Kdeg.
    ksyn  <- rbase * kdeg                                     # pg/mL/h

    # ---- 3. Micro-constants ----------------------------------------------
    k12 <- q / vc
    k21 <- q / vp

    # ---- 4. ODE system (James 2025 Eqs. 1-2) ------------------------------
    # dLY_C/dt = K21*LY_P - K12*LY_C - Vmax*(LY_C/Vc)/(Km + LY_C/Vc)
    #            - CL*(LY_C/Vc)
    # dLY_P/dt = -K21*LY_P + K12*LY_C
    Cc <- central / vc                                        # mg/L
    d/dt(central)     <- k21 * peripheral1 - k12 * central -
                         vmax * Cc / (km + Cc) - cl * Cc
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ---- 5. Epiregulin indirect response (James 2025 Eqs. 3-4) ------------
    # dEpiregulin/dt = Ksyn - Epiregulin * Kdeg
    #                  * (1 - Emax*(LY_C/Vc)^Hill / (EC50^Hill + (LY_C/Vc)^Hill))
    # The Hill term is simultaneously the paper's predicted soluble target
    # engagement (Eq. 4): "If the fepixnebart-induced elimination of
    # epiregulin is assumed to be the same as the target engagement process,
    # then the Hill equation describes the fepixnebart target engagement."
    te <- emax * Cc^hill / (ec50^hill + Cc^hill)              # fraction engaged
    d/dt(epiregulin) <- ksyn - epiregulin * kdeg * (1 - te)
    epiregulin(0)    <- rbase

    # ---- 6. Observations --------------------------------------------------
    Epi <- epiregulin                                         # pg/mL

    Cc  ~ prop(propSd)
    Epi ~ prop(propSd_Epi)
  })
}
