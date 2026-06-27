Netterberg_2018_breast_cancer_FN_biomarkers <- function() {
  description <- "Joint turnover model for serum interleukin-6 (IL-6) and C-reactive protein (CRP) in early breast cancer patients (n=49) receiving adjuvant FEC and docetaxel chemotherapy. Each biomarker follows a first-order indirect-response form (d/dt(BioM) = Rin_BioM * (1 + g(t)) - kout_BioM * BioM) with the Netterberg 2018 Eq. 1 Lorentzian surge function g(t) = SA / ((t - PT)/SW)^4 + 1 driving the post-dose elevation. CRP production is additionally stimulated by the relative change from baseline of IL-6 (RCFB_IL6) through a linear Slope term. Surge activation in each cycle is gated by the binary covariates MIX_ELEV_IL6 / MIX_ELEV_CRP, whose probabilities Pelevation,IL-6 (63.4%) and Pelevation,CRP (44.3%) are the NONMEM mixture-model frequencies the paper estimated."
  reference <- paste(
    "Netterberg I, Karlsson MO, Nielsen EI, Quartino AL, Lindman H, Friberg LE.",
    "The risk of febrile neutropenia in breast cancer patients following adjuvant",
    "chemotherapy is predicted by the time course of interleukin-6 and C-reactive",
    "protein by modelling.",
    "Br J Clin Pharmacol. 2018;84(3):490-500.",
    "doi:10.1111/bcp.13477.",
    sep = " "
  )
  vignette <- "Netterberg_2018_breast_cancer_FN"
  units <- list(
    time = "hour",
    dosing = "n/a (no drug doses; the biomarker surge is gated by binary cycle-level covariates MIX_ELEV_IL6 / MIX_ELEV_CRP)",
    concentration = "IL-6 in pg/mL; CRP in mg/L"
  )

  covariateData <- list(
    MIX_ELEV_IL6 = list(
      description = "Binary cycle-level indicator: 1 = this cycle has an elevated IL-6 production surge (the IL-6 surge function g_IL6(t) is active); 0 = no elevated IL-6 production in this cycle (the surge function is zeroed out). Time-fixed within a single chemotherapy cycle; resampled between cycles 1 and 4.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Implements the NONMEM mixture-model gate (Netterberg 2018 Methods 'Characterization of the IL-6 and CRP time courses'). For simulation, draw as a Bernoulli random variable with probability Pelevation_IL6 = 0.634 (Table 2) independently in each cycle. The paper reports 16 mixture subpopulations across the 2 cycles x 2 biomarkers indicator combinations; the marginal IL-6-elevation probability is 0.634 and is identical for cycles 1 and 4 (per Methods).",
      source_name = "mixture component (NONMEM MIXTURE/PMIX)"
    ),
    MIX_ELEV_CRP = list(
      description = "Binary cycle-level indicator: 1 = this cycle has an elevated CRP production surge (the CRP surge function g_CRP(t) is active); 0 = no elevated CRP production (the surge function is zeroed out, but the CRP production may still be driven by the IL-6 RCFB regulation through Slope). Time-fixed within a single chemotherapy cycle; resampled between cycles 1 and 4.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Implements the NONMEM mixture-model gate (Netterberg 2018 Methods). For simulation, draw as a Bernoulli random variable with probability Pelevation_CRP = 0.443 (Table 2) independently in each cycle. Note that an MIX_ELEV_CRP = 0 patient may still show CRP elevation through the Slope * RCFB_IL6 coupling when MIX_ELEV_IL6 = 1; the paper states the 'actual number of CRP elevations were higher than for IL-6' for this reason.",
      source_name = "mixture component (NONMEM MIXTURE/PMIX)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age at enrolment (used in the prior-to-chemotherapy TTE submodel; not used in the biomarker model).",
      units = "years",
      type = "continuous",
      notes = "Screened in the cycle-level cohort summary (Netterberg 2018 Table 1, median 54 years, range 31-73) but not retained as a covariate on biomarker dynamics. Use the Netterberg_2018_breast_cancer_FN_tte_prechemo model for the age effect on the FN hazard."
    ),
    WT = list(
      description = "Body weight at enrolment (cycle-level cohort summary; not a covariate on biomarker dynamics in the final model).",
      units = "kg",
      type = "continuous",
      notes = "Screened in the cycle-level cohort summary (Netterberg 2018 Table 1, median 70 kg, range 54-111) but not retained as a covariate on biomarker dynamics."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 49L,
    n_studies = 1L,
    age_range = "31-73 years; median 54 years (Netterberg 2018 Table 1)",
    weight_range = "54-111 kg; median 70 kg (Netterberg 2018 Table 1)",
    sex_female_pct = 100,
    disease_state = "Early breast cancer; receiving adjuvant chemotherapy.",
    dose_range = "Most (n=39) received three cycles of FEC (epirubicin 75 mg/m^2 + 5-FU 600 mg/m^2 + cyclophosphamide 600 mg/m^2) followed by three cycles of docetaxel 80 mg/m^2 Q3W. Six patients received the reverse order; small numbers received variant regimens; trastuzumab was added per local routine care when applicable.",
    regions = "Sweden (Uppsala University Hospital, ethics approval Dnr 2006/353).",
    notes = "Biomarker dataset: 445 IL-6 and 482 CRP measurements collected primarily in cycles 1 and 4 (five sampling occasions per cycle, with the first 10 enrolled patients sampled on seven occasions in cycle 1). No drug concentrations were measured (no PK in the analysis). 11 patients developed FN; 12 FN episodes overall (6 in cycle 1, 6 in cycle 4). One patient with an anomalously high baseline IL-6 in cycle 1 (~28x the cohort median) was excluded for both biomarkers in cycle 1. Source: Netterberg 2018 Methods 'Patients, treatment and data'."
  )

  ini({
    # ----- Baseline biomarker concentrations (Table 2) -----
    # Log-transformed so the random effects are log-normal; back-transformed values match Table 2.
    lbl_il6 <- log(2.50);  label("Typical baseline IL-6 concentration IL-60 (pg/mL)")   # Netterberg 2018 Table 2 IL-60 = 2.50 (RSE 9.2%)
    lbl_crp <- log(1.88);  label("Typical baseline CRP concentration CRP0 (mg/L)")      # Netterberg 2018 Table 2 CRP0 = 1.88 (RSE 12%)

    # ----- First-order turnover rate constants (Table 2) -----
    lkout_il6 <- log(0.0141); label("First-order IL-6 turnover rate constant kout_IL-6 (1/h); half-life ~49 h")  # Netterberg 2018 Table 2 kout,IL-6 = 0.0141 (RSE 25%)
    lkout_crp <- log(0.0224); label("First-order CRP turnover rate constant kout_CRP (1/h); half-life ~31 h")    # Netterberg 2018 Table 2 kout,CRP = 0.0224 (RSE 13%)

    # ----- Surge function parameters (Eq. 1: g(t) = SA / (((t - PT)/SW)^4 + 1)) -----
    # SA is dimensionless (relative increase in Rin); PT is peak time (h); SW is width (h).
    lsa_il6 <- log(7.99);  label("IL-6 surge amplitude SA_IL-6 (unitless relative increase)")  # Netterberg 2018 Table 2 SA_IL-6 = 7.99 (RSE 16%)
    lsa_crp <- log(4.40);  label("CRP surge amplitude SA_CRP (unitless relative increase)")    # Netterberg 2018 Table 2 SA_CRP = 4.40 (RSE 21%)
    lsw_il6 <- log(32.4);  label("IL-6 surge width SW_IL-6 (h)")                                  # Netterberg 2018 Table 2 SW_IL-6 = 32.4 (RSE 11%)
    lsw_crp <- log(53.8);  label("CRP surge width SW_CRP (h)")                                    # Netterberg 2018 Table 2 SW_CRP = 53.8 (RSE 17%)
    lpt_il6 <- log(137.0); label("IL-6 surge peak time PT_IL-6 (h after cycle start)")            # Netterberg 2018 Table 2 PT_IL-6 = 137 (RSE 9.7%)
    # CRP peak time is constrained to be the IL-6 peak time + an estimated additional time (PT_CRP+ = 50.3 h);
    # therefore PT_CRP = 137 + 50.3 = 187.3 h (paper text rounds to 187).
    lpt_crp_plus <- log(50.3); label("Time added to PT_IL-6 to get the CRP peak time, PT_CRP+ (h)")  # Netterberg 2018 Table 2 PT_CRP+ = 50.3 (RSE 32%)

    # ----- IL-6 -> CRP regulation slope -----
    # The CRP production is stimulated by the relative change from baseline of IL-6
    # (RCFB_IL6(t)). Per Netterberg 2018 Results: "The model improved when the CRP
    # production was stimulated by a change in IL-6 [RCFBIL6(t)] using a linear
    # function (OFV dropped 61 units)." Coupled-equation form: see model() below.
    slope_il6_crp <- 1.05; label("Linear slope relating RCFB_IL-6(t) to CRP production (unitless)")  # Netterberg 2018 Table 2 Slope = 1.05 (RSE 18%); units RCFBIL-6(t)^(-1) per Table 2 header

    # ----- Mixture-model probabilities (Table 2; Table S1 for the 16 subpopulations) -----
    # The paper reports Pelevation,IL-6 = 63.4% (RSE 10%) and Pelevation,CRP = 44.3%
    # (RSE 20%); these are not estimated parameters in this rxode2/nlmixr2 translation
    # but the marginal probabilities the user feeds into the simulation as Bernoulli
    # draws for MIX_ELEV_IL6 / MIX_ELEV_CRP per cycle. They are recorded in `covariateData`
    # notes for the gating covariates (see top of file).

    # ----- Inter-individual variability (IIV) -- log-normal; Table 2 reports %CV -----
    # %CV -> log-scale variance: var = log(CV^2 + 1).
    # IL-60: 68.0% CV -> var = log(0.680^2 + 1) = 0.3858
    # CRP0:  80.5% CV -> var = log(0.805^2 + 1) = 0.4862
    # kout,IL-6: 130%  CV -> var = log(1.30^2 + 1) = 0.9883
    etalbl_il6  ~ 0.3858  # Netterberg 2018 Table 2 IIV(IL-60) = 68.0% CV
    etalbl_crp  ~ 0.4862  # Netterberg 2018 Table 2 IIV(CRP0)  = 80.5% CV
    etalkout_il6 ~ 0.9883 # Netterberg 2018 Table 2 IIV(kout,IL-6) = 130%  CV

    # ----- Interoccasion variability (IOV; encoded as additional log-normal noise on the surge params) -----
    # Table 2 IOV(SA_CRP) = 61.4% CV  -> var = log(0.614^2 + 1) = 0.3198
    # Table 2 IOV(SW_CRP) = 83.8% CV  -> var = log(0.838^2 + 1) = 0.5219
    # Table 2 IOV(PT_IL-6) = 59.7% CV -> var = log(0.597^2 + 1) = 0.3036
    # Table 2 IOV(PT_CRP+) = 81.3% CV -> var = log(0.813^2 + 1) = 0.4920
    # rxode2/nlmixr2 does not have a first-class IOV construct; in this single-occasion
    # simulation model the IOV is encoded as additional eta terms applied per simulation
    # run. Users simulating multiple cycles per subject should resample these etas at
    # cycle change.
    etalsa_crp     ~ 0.3198  # Netterberg 2018 Table 2 IOV(SA_CRP)  = 61.4% CV
    etalsw_crp     ~ 0.5219  # Netterberg 2018 Table 2 IOV(SW_CRP)  = 83.8% CV
    etalpt_il6     ~ 0.3036  # Netterberg 2018 Table 2 IOV(PT_IL-6) = 59.7% CV
    etalpt_crp_plus ~ 0.4920 # Netterberg 2018 Table 2 IOV(PT_CRP+) = 81.3% CV

    # ----- Residual error -----
    # IL-6: proportional 54.7% (Table 2 'Proportional error (IL-6) = 54.7%'). No log transform.
    # CRP: "For CRP log-transformed concentrations were used in the analysis. Proportional
    # (or additive on the log-scale) and combined ... residual error models were investigated."
    # The reported 53.0% proportional error on CRP is "additive on the log-scale"
    # ~ proportional on the linear scale (standard NONMEM convention).
    propSd_Cc_il6 <- 0.547; label("Proportional residual error on IL-6 (fraction)")  # Netterberg 2018 Table 2 Proportional error (IL-6) = 54.7% (RSE 4.7%)
    propSd_Cc_crp <- 0.530; label("Proportional residual error on CRP (fraction; additive on log-scale ~ proportional on linear scale)")  # Netterberg 2018 Table 2 Proportional error (CRP) = 53.0% (RSE 4.1%)
  })

  model({
    # ----- Individual parameters (log-normal IIV) -----
    bl_il6     <- exp(lbl_il6   + etalbl_il6)
    bl_crp     <- exp(lbl_crp   + etalbl_crp)
    kout_il6_i <- exp(lkout_il6 + etalkout_il6)
    kout_crp_i <- exp(lkout_crp)  # no IIV per Table 2

    # ----- Surge parameters (with IOV applied as additional log-normal noise) -----
    sa_il6_i      <- exp(lsa_il6)                                # SA_IL-6 has no IOV per Table 2
    sa_crp_i      <- exp(lsa_crp      + etalsa_crp)              # IOV(SA_CRP) = 61.4% CV
    sw_il6_i      <- exp(lsw_il6)                                # SW_IL-6 has no IOV per Table 2
    sw_crp_i      <- exp(lsw_crp      + etalsw_crp)              # IOV(SW_CRP) = 83.8% CV
    pt_il6_i      <- exp(lpt_il6      + etalpt_il6)              # IOV(PT_IL-6) = 59.7% CV
    pt_crp_plus_i <- exp(lpt_crp_plus + etalpt_crp_plus)         # IOV(PT_CRP+) = 81.3% CV
    pt_crp_i      <- pt_il6_i + pt_crp_plus_i                    # paper constraint: PT_CRP = PT_IL-6 + PT_CRP+

    # ----- Surge functions (Netterberg 2018 Eq. 1, exponent fixed to 4) -----
    # Gated by the cycle-level MIX_ELEV_IL6 / MIX_ELEV_CRP binary covariates: when the
    # covariate is 0 the surge function is zeroed (matching the paper's
    # mixture-component "Equation 1 set to 0 for the subpopulations that had
    # no elevated biomarker production").
    g_il6 <- MIX_ELEV_IL6 * sa_il6_i / (((t - pt_il6_i) / sw_il6_i)^4 + 1)
    g_crp <- MIX_ELEV_CRP * sa_crp_i / (((t - pt_crp_i) / sw_crp_i)^4 + 1)

    # ----- Steady-state Rin from BioM0 = Rin / kout (Netterberg 2018 Eq. 2b) -----
    rin_il6 <- kout_il6_i * bl_il6
    rin_crp <- kout_crp_i * bl_crp

    # ----- Turnover ODEs (Netterberg 2018 Eq. 2a) -----
    # IL-6: standard surge.
    d/dt(il6) <- rin_il6 * (1 + g_il6) - kout_il6_i * il6
    il6(0)    <- bl_il6

    # CRP: surge plus linear regulation by RCFB_IL-6(t) = (il6 - bl_il6) / bl_il6 with
    # the Slope coefficient (Netterberg 2018 Results "The model improved when the CRP
    # production was stimulated by a change in IL-6 [RCFBIL6(t)] using a linear
    # function (OFV dropped 61 units)"). The exact differential-equation listing is
    # in Supplementary Material 1 (not included with this on-disk extraction);
    # this is the standard linear-regulation form consistent with the Slope units
    # `RCFBIL-6(t)^(-1)` reported in Table 2.
    rcfb_il6 <- (il6 - bl_il6) / bl_il6
    d/dt(crp) <- rin_crp * (1 + g_crp + slope_il6_crp * rcfb_il6) - kout_crp_i * crp
    crp(0)    <- bl_crp

    # ----- Observations -----
    # IL-6 in pg/mL (linear; proportional residual)
    # CRP in mg/L  (proportional residual; paper analysed log-transformed CRP with
    # "additive on log-scale" error, which is proportional on the linear scale)
    Cc_il6 <- il6
    Cc_crp <- crp
    Cc_il6 ~ prop(propSd_Cc_il6)
    Cc_crp ~ prop(propSd_Cc_crp)
  })
}
