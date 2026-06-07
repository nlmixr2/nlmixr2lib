Keizer_2011_E7820_human <- function() {
  description <- paste(
    "Population PK/PD model for the alpha2-integrin inhibitor E7820 in patients with",
    "advanced solid tumors or lymphoma (Keizer 2011 clinical column). One-compartment oral PK",
    "with first-order absorption (PK structure and parameter values inherited from an earlier",
    "phase I popPK analysis of the same study and reproduced in Keizer 2011 Table II; the",
    "absorption model was simplified from the original turnover-absorption form to a first-order",
    "form to ease multi-dose simulations). PD is an indirect-response (turnover) model for",
    "alpha2-integrin expression on platelets, with an Emax inhibition function (Emax fixed at 1,",
    "Hill exponent gamma fixed at 1) acting on the input rate kin. BSV is reported on baseline",
    "integrin expression and on drug sensitivity (IC50). No tumor-growth submodel is included in",
    "the clinical analysis (Keizer 2011 Figure 3 caption: 'The clinical model had the same",
    "structure, but did not incorporate a sub-model for tumor size'). Parameter values from",
    "Keizer 2011 Tables II (clinical PK) and III (clinical integrin PD)."
  )
  reference <- paste(
    "Keizer RJ, Funahashi Y, Semba T, Wanders J, Beijnen JH, Schellens JHM, Huitema ADR.",
    "Evaluation of alpha2-integrin expression as a biomarker for tumor growth inhibition for",
    "the investigational integrin inhibitor E7820 in preclinical and clinical studies.",
    "AAPS J. 2011;13(2):230-239.",
    "doi:10.1208/s12248-011-9260-2.",
    "PK structure (clinical column of Table II) is annotated 'Estimates(5)' in the source,",
    "pointing to the earlier phase I popPK report on the same study; values are reproduced inline",
    "in Table II and used here as fixed inputs to the PD analysis (paper Methods: 'The empirical",
    "Bayesian estimates ... and typical parameter values for ka obtained from this PK model and",
    "the observed data were used to drive the PD model for alpha2-integrin expression').",
    sep = " "
  )
  vignette <- "Keizer_2011_E7820"
  units <- list(
    time          = "day",
    dosing        = "mg",
    concentration = "ng/mL (E7820 plasma concentration); integrin expression in MESF (molecules of equivalent soluble fluorochrome)"
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 36L,
    n_studies      = 1L,
    age_range      = "40.0-82.0 years",
    age_median     = "64.8 years (mean)",
    weight_range   = "43.2-113.6 kg",
    weight_median  = "70.9 kg (mean)",
    height_range   = "150.5-182.9 cm",
    bsa_range      = "1.375-2.402 m^2 (mean 1.805)",
    sex_female_pct = round(16 / 36 * 100, 1),
    race_ethnicity = c(Caucasian = 86.1, Black = 2.8, Hispanic = 11.1),
    disease_state  = "advanced solid tumors or lymphoma (phase I oncology dose-escalation)",
    dose_range     = "10, 20, 40, 70, 100, or 200 mg E7820 once daily for 28 days, followed by 7-day washout, repeated up to 9 cycles",
    regions        = "(not reported in the modelling paper)",
    notes          = paste(
      "Table I of Keizer 2011 lists 36 patients enrolled across the six dose levels (3, 4, 3, 3,",
      "17, 6 at 10, 20, 40, 70, 100, 200 mg respectively). The PD analysis used 462 alpha2-integrin",
      "measurements at 209 unique timepoints from 29 of the 36 patients. Blood samples for",
      "integrin expression were collected predose and 6 h post-dose on cycle 1 day 1, predose and",
      "24 h post-dose on cycle 1 day 28, and predose on day 28 of subsequent cycles.",
      "Integrin expression on platelets was measured by flow cytometry with FITC-conjugated",
      "anti-integrin alpha2 antibody and reported in MESF (molecules of equivalent soluble",
      "fluorochrome) units."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Stage 1: E7820 PK (clinical column, Keizer 2011 Table II).
    # Values from an earlier phase I popPK report ('Estimates(5)'),
    # reproduced inline in Table II and used as fixed inputs to the PD
    # analysis in this paper (Methods: EBEs and typical ka were carried
    # over and the absorption model simplified to first-order). All
    # encoded as fixed() because they were not re-estimated jointly with
    # the PD model in this paper.
    # Reported in h^-1 units; converted to day^-1 for a single time axis.
    # ------------------------------------------------------------------
    lka  <- fixed(log(21.336))   ; label("Absorption rate ka (1/day; = 0.889 1/h)")                       # Keizer 2011 Table II clinical ka = 0.889 h^-1 (RSE 7%; carried from ref (5))
    lcl  <- fixed(log(149.76))   ; label("Apparent clearance CL/F (L/day; = 6.24 L/h)")                   # Keizer 2011 Table II clinical CL = 6.24 L/h (RSE 5%; carried from ref (5))
    lvc  <- fixed(log(66.0))     ; label("Apparent central volume V/F (L)")                               # Keizer 2011 Table II clinical V = 66.0 L (RSE 8%; carried from ref (5))

    # ------------------------------------------------------------------
    # Stage 2: alpha2-integrin indirect-response model
    # (Keizer 2011 Table III clinical column).
    # Baseline I_base = 8350 MESF. kin reported as 825.6 MESF/day;
    # kout reported as 0.099 /day. Cross-check at steady state:
    # kin / I_base = 825.6 / 8350 = 0.0988 ~= 0.099 -- consistent.
    # The model below makes kin a derived quantity from rbase and kout.
    # ------------------------------------------------------------------
    lrbase   <- log(8350)    ; label("Baseline alpha2-integrin expression I_base (MESF)")               # Keizer 2011 Table III clinical I_base = 8350 MESF (RSE 1%)
    lkout    <- log(0.099)   ; label("Turnover rate kout for integrin pool (1/day)")                    # Keizer 2011 Table III clinical kout = 0.099 day^-1
    emax_int <- fixed(1)     ; label("Maximal Emax of E7820 on integrin input rate (fixed at 1)")       # Keizer 2011 Table III clinical Emax,C = 1 (fixed; same as preclinical model)
    lic50    <- log(2840)    ; label("E7820 plasma conc at 50% maximal integrin inhibition IC50 (ng/mL)") # Keizer 2011 Table III clinical IC50 = 2840 ng/mL (RSE 32%)
    lhill_int <- fixed(log(1)) ; label("Hill exponent on E7820->integrin (unitless; fixed at 1)")         # Keizer 2011 Table III clinical gamma = 1 (fixed; same as preclinical model)

    # ------------------------------------------------------------------
    # Inter-individual variability.
    #   PK BSV from Keizer 2011 Table II (CV%, with CL~V correlation):
    #     CL  CV 51%  -> var = log(1 + 0.51^2) = 0.23146
    #     V   CV 48%  -> var = log(1 + 0.48^2) = 0.20733
    #     rho CL~V = 62%  -> cov = 0.62 * sqrt(var_CL * var_V) = 0.13580
    #   PD BSV from Keizer 2011 Table III (CV%):
    #     I_base CV 40%   -> var = log(1 + 0.40^2) = 0.14842
    #     IC50   CV 110%  -> var = log(1 + 1.10^2) = 0.79299
    # PK BSV values are carried from the earlier popPK report (ref (5))
    # reproduced in Table II; encoded here as a fixed block-diagonal
    # variance via the standard ~ c(var, cov, var) syntax. They are
    # *not* re-estimated in this paper; if a downstream re-fit is done,
    # treat these as starting values for the block.
    # ------------------------------------------------------------------
    etalcl + etalvc ~ c(0.23146, 0.13580, 0.20733)  # Keizer 2011 Table II clinical omega_CL = 51%, omega_V = 48%, rho CL~V = 62%
    etalrbase ~ 0.14842  # Keizer 2011 Table III clinical omega_I_base = 40% (RSE 19%)
    etalic50  ~ 0.79299  # Keizer 2011 Table III clinical omega_eff = 110% (RSE 32%)

    # ------------------------------------------------------------------
    # Residual error. Paper reports exponential residual error on integrin
    # (Keizer 2011 Table III clinical column: sigma_exp = 11.5%). PK in the
    # PD analysis was deterministic (EBEs used; PK residual error
    # disregarded), so no Cc residual error is encoded.
    # ------------------------------------------------------------------
    expSd_integrin <- 0.115  ; label("Exponential residual SD on integrin expression (fraction)")       # Keizer 2011 Table III clinical sigma_exp = 11.5% (RSE 23%)
  })

  model({
    # ----- Individual structural parameters ------------------------------
    ka <- exp(lka)  # no BSV reported on ka (Table II clinical column reports a point estimate only)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)

    rbase    <- exp(lrbase + etalrbase)
    kout     <- exp(lkout)
    ic50     <- exp(lic50 + etalic50)
    hill_int <- exp(lhill_int)

    # Steady-state synthesis rate: at baseline d(integrin)/dt = 0 implies
    # kin = kout * rbase = 0.099 * 8350 = 826.7 MESF/day (matches Table III).
    kin <- kout * rbase

    # ----- PK ODEs --------------------------------------------------------
    # Dose in mg, V in L: central state in mg, central/vc in mg/L = ug/mL.
    # Convert to ng/mL for comparison with IC50 (in ng/mL).
    kel <- cl / vc
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    Cc <- central / vc * 1000  # ng/mL

    # ----- Integrin turnover ODE -----------------------------------------
    # Keizer 2011 (Integrin Expression in Mice section, also used here):
    #   dI/dt = kin * (1 - Emax * Cp^gamma / (IC50^gamma + Cp^gamma)) - kout * I
    # With Emax = 1 and gamma = 1 (both fixed), the inhibition collapses to
    #   inh = Cp / (IC50 + Cp).
    inh_int <- emax_int * Cc^hill_int / (ic50^hill_int + Cc^hill_int)
    d/dt(integrin) <- kin * (1 - inh_int) - kout * integrin
    integrin(0)    <- rbase

    # ----- Observations --------------------------------------------------
    integrin ~ lnorm(expSd_integrin)
  })
}
