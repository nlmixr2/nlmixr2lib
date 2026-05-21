Svensson_2017_bedaquiline <- function() {
  description <- "Pharmacodynamic exposure-response model for the mycobacterial load (MBL, n bacteria per sample inoculum) in adult patients with drug-resistant pulmonary tuberculosis treated with bedaquiline plus an optimized background regimen. The latent MBL state declines mono-exponentially with a half-life HL that is prolonged by 28.1% in patients with pre-XDR or XDR tuberculosis and shortened by individual bedaquiline weekly-average plasma concentration CAV via an Emax model with the maximum fractional effect on HL fixed at -100% (EC50 1.42 mg/L). The per-subject starting MBL_0 is informed by the baseline mean Time-to-Positivity in MGIT liquid culture (TTP_MGIT_BASE) via a power-form covariate (exponent -3.69 around the cohort median 6.8 days). Inter-individual variability on log HL uses a Box-Cox-transformed eta distribution (Petersson 2009 form, shape 0.66, variance 0.33); inter-occasion variability in sputum sampling on log MBL (variance 3.71) is folded into the residual log-scale error. The Svensson 2017 source's full 3-component model (longitudinal MBL plus per-sample probability of bacterial presence plus MGIT-tube logistic-growth-driven time-to-event for observed TTP) is reduced here to the MBL component, with the latent MBL state treated directly as the observable; the probability-of-presence and tube-growth-driven TTP machinery are measurement-model artifacts of how MBL was inferred from TTP data and are dropped (see vignette Assumptions and deviations). Bedaquiline CAV is supplied as a time-varying covariate column from any popPK source; the upstream popPK paper (Svensson 2016 CPT PSP, reference 21) is shipped in nlmixr2lib as modellib('Svensson_2016_bedaquiline')."
  reference <- paste(
    "Svensson E. M., Karlsson M. O. (2017).",
    "Modelling of mycobacterial load reveals bedaquiline's exposure-response",
    "relationship in patients with drug-resistant TB.",
    "Journal of Antimicrobial Chemotherapy 72(12):3398-3405.",
    "doi:10.1093/jac/dkx317.",
    "Bedaquiline weekly-average plasma concentration CAV is computed from",
    "the upstream population PK model of Svensson E. M., Dosne A.-G.,",
    "Karlsson M. O. (2016).",
    "Population pharmacokinetics of bedaquiline and metabolite M2 in",
    "drug-resistant tuberculosis patients - the effect of time-varying",
    "weight and albumin. CPT Pharmacometrics Syst Pharmacol 5(12):682-691.",
    "doi:10.1002/psp4.12147;",
    "see modellib('Svensson_2016_bedaquiline').",
    "Box-Cox-transformed IIV distribution follows Petersson K. J. F.,",
    "Hanze E., Savic R. M. et al. (2009).",
    "Semiparametric distributions with estimated shape parameters.",
    "Pharm Res 26(9):2174-2185.",
    "doi:10.1007/s11095-009-9931-1.",
    sep = " "
  )
  vignette <- "Svensson_2017_bedaquiline"
  units <- list(time = "week", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CAV = list(
      description        = "Bedaquiline weekly-average plasma concentration (Cav,W in the source). Time-varying covariate updated weekly to reflect the loading-vs-maintenance dose schedule and the slow tissue accumulation of bedaquiline over the treatment period. Set to 0 for placebo periods.",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used in the Emax effect on MBL half-life: `bdq_factor = 1 - emax_bdq * CAV / (ec50_bdq + CAV)`. Set to 0 for placebo subjects; the Emax term then evaluates to 1 and HL is unchanged from the typical placebo half-life of 0.81 weeks. Svensson 2017 obtains individual CAV values as empirical-Bayes secondary metrics from the upstream Svensson 2016 (CPT PSP) population PK model for bedaquiline; users can supply CAV from any popPK source (e.g., by post-processing a `modellib('Svensson_2016_bedaquiline')` simulation into a 7-day rolling mean of central-compartment concentration in mg/L) or as a constant scalar representative of a target exposure scenario. The bedaquiline maximum effect on HL is fixed at -100% (emax_bdq = 1 FIX in `ini()`), so as CAV grows without bound, bdq_factor approaches 0 and HL approaches 0 -- biologically the model extrapolates poorly above the observed exposure range and Svensson 2017 cautions against simulation at markedly higher bedaquiline doses (Discussion paragraph 7).",
      source_name        = "Cav,W"
    ),
    DIS_TB_XDR = list(
      description        = "Pre-XDR or XDR tuberculosis drug-resistance indicator. 1 = subject's Mycobacterium tuberculosis isolate is classified as pre-extensively-drug-resistant (pre-XDR; resistant to isoniazid + rifampicin plus a second-line fluoroquinolone OR an injectable, but not both) or extensively-drug-resistant (XDR; resistant to isoniazid + rifampicin plus a second-line fluoroquinolone AND an injectable); 0 = multidrug-resistant (MDR; resistant to isoniazid + rifampicin only), drug-susceptible, or unclassified TB. Time-fixed per subject.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (MDR / susceptible / missing-treated-as-MDR; the typical-value reference). Per Svensson 2017 Results paragraph 4, the 19% of subjects with missing TB-type information were assigned to the reference (MDR) group after testing that the missing group did not differ significantly from MDR.",
      notes              = "Multiplicative effect on the MBL half-life: `xdr_factor = 1 + e_xdr_hl * DIS_TB_XDR` with `e_xdr_hl = 0.281` (Svensson 2017 Table 2 (pre-)XDR effect on half-life MBL = 28.1% (95% CI 9.1-51.5%)). Pre-XDR/XDR patients have a 28.1% longer MBL half-life than MDR/susceptible patients; per Svensson 2017 Table 3 this translates into a 2-4 weeks longer median time-to-sputum-culture-conversion and a notably lower SCC rate at week 20 under matched bedaquiline exposure.",
      source_name        = "(pre-)XDR"
    ),
    TTP_MGIT_BASE = list(
      description        = "Baseline (pre-treatment) mean time-to-positivity in the mycobacterial growth indicator tube (MGIT) liquid culture system, in days. Computed in the source as the mean of three replicate spot-sputum-sample TTP values collected the day before the start of treatment. Time-fixed per subject.",
      units              = "days",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form covariate on the starting mycobacterial load: `mbl0_i = exp(lmbl0) * (TTP_MGIT_BASE / 6.8)^e_ttp_mbl0` with reference 6.8 days (Svensson 2017 Table 1 cohort median, n = 191) and `e_ttp_mbl0 = -3.69` (Svensson 2017 Table 2 baseline TTP effect on MBL_0). The exponent is negative because longer baseline TTP corresponds to fewer viable bacteria in the inoculum and a lower starting MBL. Per Svensson 2017 Discussion paragraph 3, the estimated effect predicts a four-times longer median time-to-sputum-culture-conversion in patients with the lowest baseline TTP versus the highest.",
      source_name        = "mTTP0"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 189L,
    n_studies      = 1L,
    age_range      = "18-63 years (Table 1; trial enrolment criterion 18-65 years)",
    age_median     = "33 years",
    weight_range   = "35-83 kg (Table 1)",
    weight_median  = "54 kg",
    sex_female_pct = 34.5,
    race_ethnicity = c(
      Caucasian = 10.2,
      Black     = 39.8,
      Hispanic  = 13.6,
      Asian     = 7.3,
      Other     = 28.6,
      Missing   = 0.5
    ),
    disease_state  = "Adult patients with newly diagnosed multidrug-resistant pulmonary tuberculosis (MDR-TB), with pre-XDR (25.1%) and XDR (5.31%) strata pooled into the DIS_TB_XDR = 1 group and drug-susceptible (3.86%) + MDR (46.4%) + missing (18.9%, assigned to MDR) pooled into the DIS_TB_XDR = 0 reference group. Lung cavitation prevalence 91.7%. HIV co-infection allowed if CD4+ T-cell count >= 300 cells/mm3 and patient not on antiretroviral therapy; HIV-positive rate 14.6%.",
    dose_range     = "Open-label oral bedaquiline 400 mg once daily for the first 2 weeks (loading), then 200 mg three times weekly for 22 weeks in the active arm; placebo equivalent in the placebo arm; both arms received an optimized background regimen of five second-line anti-TB drugs (kanamycin, ofloxacin, ethionamide, pyrazinamide, terizidone, with a few predefined substitutions allowed) for 18-24 months under directly observed therapy. The PD model is fitted to TTP observations through week 8 (stage 1) or week 24 (stage 2) and does not contain a dosing compartment itself; bedaquiline exposure enters via the time-varying CAV covariate.",
    regions        = "Multicenter international (TMC207-C208 phase IIb registration trial; specific regional breakdown not reported in the main text).",
    notes          = "Baseline demographics from Svensson 2017 Table 1 (full enrolled cohort n = 206 with TTP data; n = 189 included in the final model fit after exclusions described in the Supplementary data). Trial: TMC207-C208 (ClinicalTrials.gov NCT00449644), Janssen Pharmaceuticals, shared with the authors through the PreDiCT-TB consortium (http://www.predict-tb.eu). The cohort was originally enrolled in two stages: stage 1 (47 subjects; 23 on bedaquiline) received the randomized intervention for 8 weeks, stage 2 (159 subjects; 79 on bedaquiline) received the randomized intervention for 24 weeks; both stages contributed to the model dataset. Baseline TTP values from before initiation of treatment were used as a covariate (TTP_MGIT_BASE) and not as observations. The dataset used for model building included 5833 TTP observations (56.6% positive) from 189 individuals (98 in the placebo arm and 91 in the bedaquiline arm)."
  )

  ini({
    # ========================================================================
    # MBL in patients (Svensson 2017 Table 2 'MBL in patients' subblock).
    # ------------------------------------------------------------------------
    # Time unit: weeks (matches the source's TAST and HL units).
    # ========================================================================

    # Typical starting mycobacterial load per sample inoculum at the start
    # of treatment.
    lmbl0 <- log(2.14e3)        ; label("Typical starting MBL_0 (n bacteria/inoculum, log scale)") # Svensson 2017 Table 2 MBL_0 = 2.14e3 (95% CI 1.39e3-3.46e3)

    # Baseline-TTP power effect on starting MBL_0.  Power form per Svensson
    # 2017 Equation (1): `MBL_0_i = MBL_0 * (mTTP_0_i / mTTP_0_p)^COVTTP`
    # with the reference baseline TTP `mTTP_0_p` = 6.8 days (Table 1
    # cohort median, n = 191).  Negative exponent because longer baseline
    # TTP -> fewer baseline bacteria in the inoculum.
    e_ttp_mbl0 <- -3.69         ; label("Power exponent of baseline TTP on starting MBL_0 (unitless)") # Svensson 2017 Table 2 baseline TTP effect on MBL_0 = -3.69 (95% CI -4.15 to -3.30)

    # Typical half-life of MBL decline in the placebo arm (no bedaquiline,
    # MDR-TB).  Estimated point estimate from the Box-Cox-IIV model fit.
    lhl <- log(0.81)            ; label("Typical half-life of MBL decline (weeks, log scale)") # Svensson 2017 Table 2 half-life MBL = 0.81 weeks (95% CI 0.71-0.93)

    # Pre-XDR / XDR effect on MBL half-life (additive fractional).  Pooled
    # pre-XDR + XDR stratum has a 28.1% longer half-life than the MDR /
    # susceptible reference.
    e_xdr_hl <- 0.281           ; label("Multiplicative pre-XDR/XDR effect on MBL half-life (unitless fraction)") # Svensson 2017 Table 2 (pre-)XDR effect on half-life MBL = 28.1% (95% CI 9.1-51.5)

    # Bedaquiline EC50 for the Emax effect on half-life.  In mg/L total
    # plasma concentration (no albumin adjustment - the paper tested
    # adjustment for albumin-driven free fraction and found it did not
    # improve the fit; Discussion paragraph 4).
    lec50_bdq <- log(1.42)      ; label("Bedaquiline EC50 for Emax effect on MBL half-life (mg/L, log scale)") # Svensson 2017 Table 2 EC50 bedaquiline = 1.42 mg/L (95% CI 1.00-2.05)

    # Bedaquiline maximum fractional effect on HL is fixed at -100% by the
    # source (because the limited range of observed exposures did not
    # support estimating Emax < 1).  Wrapped in fixed() so the constraint
    # is load-bearing provenance.
    emax_bdq <- fixed(1)        ; label("Bedaquiline maximum fractional effect on HL (fixed at 1 = -100%)") # Svensson 2017 Table 2 bedaquiline maximal effect on half-life MBL = -1 FIX

    # Box-Cox shape parameter for the IIV transform on log HL (Petersson
    # 2009 form, `eta_BC = (exp(etalhl)^lambda_bc - 1) / lambda_bc`; see
    # model() comment for the per-eta application).  At lambda_bc -> 0
    # this collapses to a plain log-normal IIV; at lambda_bc = 0.66 it
    # produces a moderately right-skewed half-life distribution.
    lambda_bc <- 0.66           ; label("Box-Cox shape parameter for IIV transform on log HL (Petersson 2009)") # Svensson 2017 Table 2 Box-Cox transformation IIV half-life MBL = 0.66 (95% CI 0.34-1.05)

    # IIV on (Box-Cox-transformed) log HL.  Variance, on the original
    # (pre-transform) eta scale.
    etalhl ~ 0.33               # Svensson 2017 Table 2 IIV half-life MBL (variance) = 0.33 (95% CI 0.25-0.45)

    # Residual log-scale error on the observed MBL.  Folds the source's
    # IOV_sputum (between-sample lognormal variability in the sputum
    # sampling procedure; variance 3.71 in the source) into a single
    # log-normal residual on the latent MBL state, because the extraction
    # treats MBL as directly observable (see description and vignette
    # Assumptions and deviations for the deviation from the source's full
    # 3-component measurement model).
    expSd_MBL <- sqrt(3.71)     ; label("Log-scale residual SD on observed MBL (folds IOV_sputum from source)") # Svensson 2017 Table 2 IOV sputum sampling MBL (variance) = 3.71 (95% CI 3.29-4.38); SD = sqrt(3.71) = 1.926
  })

  model({
    # ========================================================================
    # 1. Box-Cox transform on the IIV eta for log HL.  Source form
    #    (Petersson 2009, reference 29 of Svensson 2017; see also the
    #    NM-TRAN translation used by Schoemaker_2018_levetiracetam and the
    #    .ctl line `TETA1 = (EXP(ETA(1))**SHP1 - 1) / SHP1`).
    # ========================================================================
    eta_bc <- (exp(etalhl)^lambda_bc - 1) / lambda_bc

    # ========================================================================
    # 2. Bedaquiline Emax effect on MBL half-life.  Source Emax model
    #    (Svensson 2017 Results 'Covariate and PK effects'): the maximum
    #    fractional shortening of HL is fixed to -100% and the EC50 is
    #    1.42 mg/L of bedaquiline total plasma concentration.  CAV is the
    #    weekly-average bedaquiline concentration supplied as a
    #    time-varying covariate column; for placebo subjects CAV = 0 and
    #    bdq_factor = 1.
    # ========================================================================
    ec50_bdq   <- exp(lec50_bdq)
    bdq_factor <- 1 - emax_bdq * CAV / (ec50_bdq + CAV)

    # ========================================================================
    # 3. Pre-XDR / XDR multiplicative effect on MBL half-life.
    # ========================================================================
    xdr_factor <- 1 + e_xdr_hl * DIS_TB_XDR

    # ========================================================================
    # 4. Individual MBL half-life and first-order decay rate.  The Box-Cox
    #    IIV multiplies the typical HL via exp(eta_bc); the XDR and
    #    bedaquiline factors then scale HL additionally.
    # ========================================================================
    hl_i    <- exp(lhl + eta_bc) * xdr_factor * bdq_factor
    kel_mbl <- log(2) / hl_i

    # ========================================================================
    # 5. Individual starting MBL_0 (n bacteria/inoculum) from baseline TTP.
    # ========================================================================
    mbl0_i <- exp(lmbl0) * (TTP_MGIT_BASE / 6.8)^e_ttp_mbl0

    # ========================================================================
    # 6. ODE for the latent mycobacterial load state.  Mono-exponential
    #    decline per Svensson 2017 Equation (1).  Initial condition at
    #    time 0 (start of treatment, TAST = 0) is the per-subject
    #    baseline-TTP-informed MBL_0.
    # ========================================================================
    d/dt(mbl) <- -kel_mbl * mbl
    mbl(0)    <- mbl0_i

    # ========================================================================
    # 7. Observation.  The Svensson 2017 source does not observe MBL
    #    directly - it is the latent driver of TTP / probability-of-presence
    #    observations.  This extraction treats MBL as the observable, with
    #    the source's per-sample IOV_sputum (lognormal between-sample
    #    variability on MBL with variance 3.71) folded in as the residual
    #    log-scale error.  See the vignette's Assumptions and deviations
    #    section for the deviation from the source's full measurement model.
    # ========================================================================
    MBL <- mbl
    MBL ~ lnorm(expSd_MBL)
  })
}
