Barnett_2018_coproporphyrin_I <- function() {
  description <- "Semi-mechanistic turnover model for the endogenous OATP1B-substrate biomarker coproporphyrin I (CPI) in healthy adult males (Barnett 2018), with simultaneous plasma and urine outputs and competitive rifampicin inhibition of biliary CPI clearance. CPI is produced at a zero-order synthesis rate ksyn, distributed in volume Vcpi, and eliminated via biliary clearance CLb,CPI (the dominant route, ~88% of total CL under baseline conditions) and renal clearance CLr,CPI. Rifampicin inhibits CLb,CPI competitively through KiCPI driven by the instantaneous plasma rifampicin concentration; CLr,CPI and ksyn are unaffected. A binary RIF-coadministration covariate additionally captures a paper-reported ~50% reduction in Vcpi during the rifampicin phase (Barnett 2018 Table 1)."
  reference <- paste(
    "Barnett S, Ogungbenro K, Menochet K, Shen H, Lai Y, Humphreys WG, Galetin A.",
    "Gaining Mechanistic Insight Into Coproporphyrin I as Endogenous Biomarker",
    "for OATP1B-Mediated Drug-Drug Interactions Using Population Pharmacokinetic",
    "Modeling and Simulation.",
    "Clin Pharmacol Ther. 2018;104(3):564-574.",
    "doi:10.1002/cpt.983.",
    "The rifampicin perpetrator PK is parameterised in",
    "modellib('Barnett_2018_rifampicin'); supply its central-compartment",
    "concentration as CP_RIF_UM (umol/L) to drive the competitive OATP1B",
    "inhibition term in this model.",
    sep = " "
  )
  vignette <- "Barnett_2018_coproporphyrin_I"
  units <- list(time = "hour", dosing = "none", concentration = "nmol/L")

  covariateData <- list(
    OCC = list(
      description        = "Integer-valued occasion / period indicator (Barnett 2018 study design: OCC1 = rifampicin-only period, OCC2 = rosuvastatin-only period, OCC3 = combined rifampicin + rosuvastatin period).",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Time-varying within subject; constant within an occasion. CPI plasma + urine samples were collected on all three occasions in the source clinical study (Lai et al. 2016 n=12 healthy-male SLCO1B1-wildtype cohort that supplied the Barnett 2018 fit). Decomposed inside model() into binary indicators oc1, oc2, oc3 that multiplex the per-occasion IOV etas on log-ksyn and log-CLb,CPI; Barnett 2018 Table 1 reports a single shared IOV variance per parameter across occasions, mirroring the Wilkins_2008_rifampicin OMEGA BLOCK(1) SAME idiom.",
      source_name        = "OCC"
    ),
    CONMED_RIF = list(
      description        = "Concomitant single-dose rifampicin co-administration indicator (1 = within a 600 mg rifampicin co-administration period in the Barnett 2018 study design; 0 = baseline or pre-RIF period).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no rifampicin co-administration; baseline phase)",
      notes              = "Time-varying within subject. Used here as the binary period-level covariate that captures Barnett 2018 Table 1's paper-reported reduction of Vcpi from 6.59 L (baseline) to 3.4 L (RIF phase). Distinct semantics from the multi-day CYP3A4-induction use of CONMED_RIF in Svensson_2014_bedaquiline (where the indicator switches on at day 3 of daily dosing): here rifampicin acts acutely as a competitive OATP1B inhibitor (single 600 mg oral dose), so CONMED_RIF = 1 from the time of the rifampicin dose through the end of the RIF-phase plasma / urine sampling window of that occasion, with no induction lag.",
      source_name        = "CONMED_RIF"
    ),
    CP_RIF_UM = list(
      description        = "Instantaneous rifampicin plasma concentration as a time-varying perpetrator covariate driving competitive OATP1B inhibition of biliary CPI clearance (Barnett 2018 Eq. 4).",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Set to 0 outside the rifampicin co-administration window so the OATP1B inhibition term in the d/dt(Ccpi) equation collapses to the baseline form (Eq. 3). The PK trajectory of rifampicin is parameterised in modellib('Barnett_2018_rifampicin'); users typically simulate the rifampicin model first and then feed its central-compartment concentration (after MW conversion to umol/L; rifampicin MW 822.94 g/mol) as the CP_RIF_UM column on the CPI event table. Reference peak: a single 600 mg oral rifampicin dose in the Barnett 2018 cohort produces a typical Cmax of ~29 umol/L after MTT-delayed absorption (computed from the modellib('Barnett_2018_rifampicin') typical-value parameters; see the validation vignette for the worked simulation).",
      source_name        = "CRIF"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 12L,
    n_studies        = 1L,
    n_observations   = 522L,
    age_range        = "Healthy adult males; demographic distribution not tabulated by Barnett 2018 (source dataset Lai et al. 2016, 12 healthy male subjects, SLCO1B1 c.521 T>C wildtype only; no OATP1B1*5 / *15 carriers).",
    weight_range     = "(not extracted; Barnett 2018 Methods do not tabulate per-subject weights for the n=12 cohort.)",
    sex_female_pct   = 0,
    disease_state    = "Healthy adult male volunteers in a three-occasion (7-day washout between OCC1-2 and OCC2-3) drug-drug-interaction crossover study; CPI was monitored as a candidate endogenous biomarker of OATP1B-mediated DDIs.",
    dose_range       = "Endogenous biomarker (no exogenous dose); rifampicin was co-administered as a 600 mg oral dose on OCC1 and OCC3 to perturb the CPI biliary clearance via OATP1B inhibition.",
    regions          = "(not extracted; the underlying clinical study Lai et al. 2016 region was not explicitly stated in Barnett 2018 Methods.)",
    notes            = "Demographics inferred from Barnett 2018 Methods (Clinical data section) and the cited source clinical study Lai Y et al., Pharmacol Res Perspect 2016;4(3):e00207. The 420 CPI plasma samples (OCC1 = 144, OCC2 = 144, OCC3 = 132) plus 102 CPI urine samples (pretreatment = 34, posttreatment = 68) were fit simultaneously in NONMEM using FOCE. Identifiability analysis was performed in DAISY (Bellu 2007); the structural CPI model was shown to be globally identifiable given the available plasma + urine data."
  )

  ini({
    # Structural parameters -- Barnett 2018 Table 1, CPI rows
    # (column 'Estimate (SE%) Population'). SE% are given in parentheses
    # next to each value below. ksyn is the zero-order synthesis rate
    # of CPI; Vcpi is the apparent volume of distribution; CLb,CPI and
    # CLr,CPI are the biliary and renal clearance components.
    lksyn   <- log(12.7);  label("CPI synthesis rate ksyn (nmol/h)")                       # Table 1, CPI row 'k syn' = 12.7 (SE 6%). The unit is reported in the paper as 'nM/h' but the steady-state baseline = ksyn / (CLb + CLr) only resolves dimensionally with ksyn in nmol/h (amount/time); see vignette Errata for the unit discussion.
    lclb    <- log(12.3);  label("CPI biliary clearance CLb,CPI (L/h, baseline)")          # Table 1, CPI row 'CL b,CPI' = 12.3 (SE 10%)
    lclr    <- log(1.64);  label("CPI renal clearance CLr,CPI (L/h)")                      # Table 1, CPI row 'CL R,CPI' = 1.64 (SE 6%)
    lvc     <- log(6.59);  label("CPI volume of distribution Vcpi (L, baseline phase)")    # Table 1, CPI row 'V CPI (L)' = 6.59 (SE 12%)
    lki     <- log(1.15);  label("Total rifampicin OATP1B inhibition constant KiCPI (umol/L)")  # Table 1, CPI row 'Ki CPI (uM) c' = 1.15 (SE 9%); footnote c notes the unbound Ki = 0.13 umol/L after correction by RIF fu = 0.11.

    # Binary RIF-coadministration covariate effect on Vcpi: the paper
    # reports Vcpi reducing from 6.59 L (baseline) to 3.4 L during the
    # rifampicin phase, which encodes as the multiplicative factor
    # (1 + e_rif_vc * CONMED_RIF) on Vcpi with e_rif_vc = (3.4 / 6.59) - 1.
    e_rif_vc <- -0.4841;   label("Multiplicative reduction in Vcpi during the rifampicin phase (unitless)")  # Table 1, CPI row 'V CPI (L) (RIF)' = 3.4 (SE 13%); derived as 3.4/6.59 - 1 = -0.4841.

    # IIV -- log-normal variances computed from Table 1 IIV CV% column
    # via the standard omega^2 = log(1 + CV^2) conversion.
    etalksyn ~ 0.00084   # Table 1, CPI IIV ksyn   = 2.9%  CV (RSE 318.7%, very poorly estimated); log(1 + 0.029^2) = 0.000841
    etalclb  ~ 0.02024   # Table 1, CPI IIV CLb    = 14.3% CV (RSE 63.8%); log(1 + 0.143^2) = 0.02024
    etalclr  ~ 0.00898   # Table 1, CPI IIV CLr    =  9.5% CV (RSE 31.9%); log(1 + 0.095^2) = 0.00898
    etalvc   ~ 0.11679   # Table 1, CPI IIV Vcpi   = 35.2% CV (RSE 24.7%); log(1 + 0.352^2) = 0.11679
    etalki   ~ 0.03473   # Table 1, CPI IIV KiCPI  = 18.8% CV (RSE 40.6%); log(1 + 0.188^2) = 0.03473

    # IOV -- Barnett 2018 reports a single shared IOV variance on ksyn
    # and on CLb,CPI across the three study occasions. First occasion's
    # eta estimated; remainder fixed to the same value (NONMEM
    # OMEGA BLOCK(1) SAME idiom).
    etaiov_ksyn_1 ~ 0.00302         # Table 1, CPI IOV ksyn = 5.5%  CV (RSE 45.6%); log(1 + 0.055^2) = 0.00302
    etaiov_ksyn_2 ~ fix(0.00302)
    etaiov_ksyn_3 ~ fix(0.00302)
    etaiov_clb_1  ~ 0.03294         # Table 1, CPI IOV CLb  = 18.3% CV (RSE 38%);   log(1 + 0.183^2) = 0.03294
    etaiov_clb_2  ~ fix(0.03294)
    etaiov_clb_3  ~ fix(0.03294)

    # Residual error -- combined additive + proportional for both
    # plasma and urine, with the plasma additive component fixed and
    # the urine additive component estimated.
    propSd        <- 0.139;        label("Proportional residual error, CPI plasma (fraction)")  # Table 1, CPI row 'r prop (%) - plasma' = 13.9 (SE 10%) -> fraction 0.139
    addSd         <- fixed(0.001); label("Additive residual error, CPI plasma (nmol/L, fixed)")  # Table 1, CPI row 'r add (nM) - plasma' = 0.001 FIXED
    propSd_Ucpi   <- 0.342;        label("Proportional residual error, CPI urine (fraction)")    # Table 1, CPI row 'r prop (%) - urine'  = 34.2 (SE 6%) -> fraction 0.342
    addSd_Ucpi    <- 2.69;         label("Additive residual error, CPI urine (nmol/L)")          # Table 1, CPI row 'r add (nM) - urine'   = 2.69 (SE 28%)
  })

  model({
    # Decompose the integer-valued OCC column into binary indicators
    # for IOV multiplexing on log-ksyn and log-CLb,CPI.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)

    iov_ksyn <- oc1 * etaiov_ksyn_1 + oc2 * etaiov_ksyn_2 + oc3 * etaiov_ksyn_3
    iov_clb  <- oc1 * etaiov_clb_1  + oc2 * etaiov_clb_2  + oc3 * etaiov_clb_3

    # Individual typical-value parameters with log-normal IIV / IOV.
    ksyn <- exp(lksyn + etalksyn + iov_ksyn)
    cl_b <- exp(lclb  + etalclb  + iov_clb)
    cl_r <- exp(lclr  + etalclr)
    vc   <- exp(lvc   + etalvc) * (1 + e_rif_vc * CONMED_RIF)
    ki   <- exp(lki   + etalki)

    # Rifampicin competitive inhibition of biliary CPI clearance
    # (Barnett 2018 Eq. 4). CP_RIF_UM is the instantaneous rifampicin
    # plasma concentration in umol/L; setting it to 0 outside the RIF
    # window collapses cl_b_eff to cl_b (baseline Eq. 3).
    cl_b_eff <- cl_b / (1 + CP_RIF_UM / ki)

    # ODE system. central is the CPI amount in plasma (nmol);
    # urine is the cumulative CPI amount eliminated via the renal
    # route (nmol). Ccpi is the plasma concentration in nmol/L.
    Ccpi <- central / vc

    d/dt(central) <- ksyn - (cl_b_eff + cl_r) * Ccpi
    d/dt(urine)   <- cl_r * Ccpi

    # Steady-state baseline initial conditions from Eq. 3 with the
    # baseline (no-RIF) clearance form: Css = ksyn / (CLb + CLr).
    # central(0) uses the CURRENT vc (including any CONMED_RIF
    # multiplicative effect) so that Ccpi(0) = central(0) / vc =
    # Css_baseline regardless of CONMED_RIF at t = 0. This anchors
    # the simulation start at the pre-perturbation typical baseline
    # of ~0.91 nmol/L; the subsequent dynamics evolve as the
    # competitive OATP1B inhibition (driven by CP_RIF_UM) and the
    # multiplicative Vcpi shift (driven by CONMED_RIF) take effect.
    Css_baseline <- ksyn / (cl_b + cl_r)
    central(0)   <- Css_baseline * vc
    urine(0)     <- 0

    # Outputs and residual error.
    # Cc = Ccpi (paper notation Ccpi); the canonical observation name
    # is reused so PKNCA-style validation works without renaming.
    Cc   <- Ccpi
    Ucpi <- urine

    Cc   ~ add(addSd)      + prop(propSd)
    Ucpi ~ add(addSd_Ucpi) + prop(propSd_Ucpi)
  })
}
