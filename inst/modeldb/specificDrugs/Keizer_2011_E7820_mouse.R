Keizer_2011_E7820_mouse <- function() {
  description <- paste(
    "Preclinical (mouse, female nude with subcutaneous KP-1 pancreatic-carcinoma xenograft).",
    "Sequential PK/PD/tumor-growth model for the alpha2-integrin inhibitor E7820 (Keizer 2011).",
    "Stage 1: one-compartment oral PK with first-order absorption, per-kg parameterisation.",
    "Stage 2: indirect-response (turnover) model for alpha2-integrin expression on platelets,",
    "with an Emax inhibition function (Emax fixed at 1, Hill exponent gamma fixed at 1) acting",
    "on the input rate kin. Stage 3: exponential tumor growth on diameter with an initial-slow-growth",
    "term (1 - exp(-beta*t)) gating the growth rate, and a sigmoidal Emax inhibition driven by relative",
    "alpha2-integrin inhibition ((I_base - integrin)/I_base) with Hill coefficient fixed at 5.",
    "Parameters from Keizer 2011 Tables II (preclinical PK), III (preclinical integrin PD), and IV",
    "(tumor growth)."
  )
  reference <- paste(
    "Keizer RJ, Funahashi Y, Semba T, Wanders J, Beijnen JH, Schellens JHM, Huitema ADR.",
    "Evaluation of alpha2-integrin expression as a biomarker for tumor growth inhibition for",
    "the investigational integrin inhibitor E7820 in preclinical and clinical studies.",
    "AAPS J. 2011;13(2):230-239.",
    "doi:10.1208/s12248-011-9260-2.",
    sep = " "
  )
  vignette <- "Keizer_2011_E7820"
  units <- list(
    time          = "day",
    dosing        = "mg/kg",
    concentration = "ng/mL (E7820 plasma concentration); integrin expression in arbitrary % units (platelet flow-cytometry signal); tumor size in mm (longest caliper diameter)"
  )

  covariateData <- list()

  population <- list(
    species        = "mouse (female nude / female KSN Slc; PK was elucidated in KSN Slc mice, tumor-growth experiments used 7-week-old female nude mice transplanted subcutaneously with KP-1 human pancreatic-carcinoma cells)",
    n_subjects     = 42L,
    n_studies      = 2L,
    age_range      = "6-8 weeks at start of PK studies; 7 weeks at tumor implantation",
    weight_range   = "(not reported in the modelling paper; per-kg parameterisation throughout)",
    sex_female_pct = 100,
    race_ethnicity = NA,
    disease_state  = "subcutaneous KP-1 human pancreatic-carcinoma xenograft (5e6 cells/head); tumor-growth experiments only",
    dose_range     = "PK: single IV 25 mg/kg, single oral 25-100 mg/kg, repeated oral 50 mg/kg (12 h interval). Tumor: vehicle, 12.5, 25, 50, 100, or 200 mg/kg oral gavage twice daily for 21 days starting 1 week after implantation.",
    regions        = "preclinical (in-vivo xenograft)",
    notes          = paste(
      "PK cohort: roughly 3 mice IV 25 mg/kg + 4 mice each at oral 25/50/100 mg/kg + 4 mice repeated 50 mg/kg.",
      "Tumor cohort: 7-week old female nude mice transplanted subcutaneously with KP-1 cells; treatment with",
      "E7820 or vehicle started 1 week after transplantation. 119 alpha2-integrin measurements and 210 tumor",
      "size measurements were available from the tumor experiment. Tumor diameter measured twice weekly by",
      "caliper; platelet alpha2-integrin expression measured weekly by flow cytometry with FITC-conjugated",
      "anti-integrin alpha2 antibody. Body weight was not used as a covariate -- all PK parameters are",
      "reported on a per-kg basis."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Stage 1: E7820 PK (preclinical column, Keizer 2011 Table II).
    # Reported in h^-1 units; converted to day^-1 here so the PK and PD
    # are on a single time axis (the PD turnover model is in days).
    # ------------------------------------------------------------------
    # CL, V, and F are reported in Keizer 2011 Methods as noncompartmental
    # mean values from a separate PK experiment, then carried as fixed
    # inputs into the PD fit. Encoded as fixed() to preserve provenance.
    lka     <- fixed(log(72))      ; label("Absorption rate ka (1/day; = 3 1/h)")                       # Keizer 2011 Table II preclinical ka = 3 h^-1 (fixed; "an absorption rate (ka) of 3 h^-1 was used throughout")
    lcl     <- fixed(log(11.304))  ; label("Apparent clearance CL/F (L/day/kg; = 0.471 L/h/kg)")        # Keizer 2011 Table II preclinical CL = 0.471 L/h/kg (NCA, fixed in PD fit)
    lvc     <- fixed(log(0.684))   ; label("Apparent central volume V/F (L/kg)")                        # Keizer 2011 Table II preclinical V = 0.684 L/kg (NCA, fixed in PD fit)
    lfdepot <- fixed(log(0.693))   ; label("Oral bioavailability F (unitless)")                         # Keizer 2011 Table II preclinical F = 69.3% (NCA, fixed in PD fit)

    # ------------------------------------------------------------------
    # Stage 2: alpha2-integrin indirect-response model
    # (Keizer 2011 Table III preclinical column).
    # Baseline I_base = 22.5 (% units on the platelet flow-cytometry scale).
    # kin (input rate) reported as 3.21 %/day; kout (output rate) reported
    # as 0.143 /day. Cross-check at steady state: kin / I_base = 3.21/22.5
    # = 0.1427 ~= 0.143 -- consistent. The model below makes kin a derived
    # quantity from rbase and kout (parameter-names.md canonical baseline
    # form: rbase + kout drive the turnover; kin solved from steady state).
    # ------------------------------------------------------------------
    lrbase   <- log(22.5)    ; label("Baseline alpha2-integrin expression I_base (% units)")            # Keizer 2011 Table III preclinical I_base = 22.5% (RSE 1%)
    lkout    <- log(0.143)   ; label("Turnover rate kout for integrin pool (1/day)")                    # Keizer 2011 Table III preclinical kout = 0.143 day^-1
    emax_int <- fixed(1)     ; label("Maximal Emax of E7820 on integrin input rate (fixed at 1)")       # Keizer 2011 Table III preclinical Emax,C = 1 (fixed; "the shape parameter gamma in the Hill equation was fixed to 1")
    lic50    <- log(656)     ; label("E7820 plasma conc at 50% maximal integrin inhibition IC50 (ng/mL)") # Keizer 2011 Table III preclinical IC50 = 656 ng/mL (RSE 18%)
    lhill_int <- fixed(log(1)) ; label("Hill exponent on E7820->integrin (unitless; fixed at 1)")         # Keizer 2011 Table III preclinical gamma = 1 (fixed; sigmoidal Emax collapsed to non-sigmoidal Emax)

    # ------------------------------------------------------------------
    # Stage 3: Tumor growth model (Keizer 2011 Table IV; preclinical only).
    # Exponential growth with an initial-slow-growth gate (1 - exp(-beta*t))
    # and a sigmoidal Emax kill term driven by relative integrin inhibition.
    # ------------------------------------------------------------------
    lrbase_tumor <- log(5.14)   ; label("Baseline tumor diameter T_base (mm)")                            # Keizer 2011 Table IV T_base = 5.14 mm (RSE 1%)
    lalpha       <- log(0.0903) ; label("Maximal exponential tumor growth rate alpha (1/day)")            # Keizer 2011 Table IV alpha = 0.0903 day^-1 (RSE 103%)
    lbeta        <- log(0.0391) ; label("Rate of decline of initial resistance to growth beta (1/day)")   # Keizer 2011 Table IV beta = 0.0391 day^-1 (RSE 130%)
    lemax_tumor  <- log(0.0472) ; label("Max effect of integrin inhibition on tumor growth Emax,I (1/day)") # Keizer 2011 Table IV Emax,I = 0.0472 day^-1 (RSE 5%)
    lii50        <- log(0.114)  ; label("Relative integrin inhibition at 50% max effect II50 (fraction)") # Keizer 2011 Table IV II50 = 11.4% (RSE 6%)
    lhill_tumor  <- fixed(log(5)) ; label("Hill exponent integrin-inhibition -> tumor growth (unitless; fixed at 5)") # Keizer 2011 Table IV gamma = 5 (fixed)

    # ------------------------------------------------------------------
    # Inter-individual variability. Source reports BSV as CV%;
    # variance = log(CV^2 + 1).
    #   Integrin model: BSV on I_base only (CV 31% -> var 0.09177).
    #   Tumor model: BSV on T_base (CV 6% -> var 0.003594) and II50 (CV 20% -> var 0.03922).
    # ------------------------------------------------------------------
    etalrbase       ~ 0.09177  # Keizer 2011 Table III preclinical omega_I_base = 31% (RSE 40%)
    etalrbase_tumor ~ 0.003594 # Keizer 2011 Table IV omega_T_base = 6% (RSE 13%)
    etalii50        ~ 0.03922  # Keizer 2011 Table IV omega_II50 = 20% (RSE 50%)

    # ------------------------------------------------------------------
    # Residual error. Paper reports exponential residual error model for
    # both the integrin and tumor outputs (sigma_exp in CV-percent form).
    # For small CV the exponential / log-normal form is numerically close
    # to a proportional model; encoded here with lnorm() residual to match
    # the source's exponential specification.
    # ------------------------------------------------------------------
    expSd_integrin   <- 0.067  ; label("Exponential residual SD on integrin expression (fraction)")     # Keizer 2011 Table III preclinical sigma_exp = 6.7% (RSE 8%)
    expSd_tumor_size <- 0.06   ; label("Exponential residual SD on tumor diameter (fraction)")          # Keizer 2011 Table IV sigma_exp = 6% (RSE 6%)
  })

  model({
    # ----- Individual structural parameters -------------------------------
    # ka, hill exponents, and Emax,C are fixed (no IIV).
    ka       <- exp(lka)
    cl       <- exp(lcl)
    vc       <- exp(lvc)
    fdepot   <- exp(lfdepot)

    rbase       <- exp(lrbase + etalrbase)
    kout        <- exp(lkout)
    ic50        <- exp(lic50)
    hill_int    <- exp(lhill_int)

    rbase_tumor <- exp(lrbase_tumor + etalrbase_tumor)
    alpha_g     <- exp(lalpha)
    beta_g      <- exp(lbeta)
    emax_tumor  <- exp(lemax_tumor)
    ii50        <- exp(lii50 + etalii50)
    hill_tumor  <- exp(lhill_tumor)

    # Steady-state synthesis rate of the integrin turnover model:
    # at baseline d(integrin)/dt = 0 implies kin = kout * rbase.
    kin <- kout * rbase

    # ----- PK ODEs (per-kg amounts; central in mg/kg, vc in L/kg) ---------
    kel <- cl / vc
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central
    f(depot)      <- fdepot

    # Plasma concentration in ng/mL: state in mg/kg, vc in L/kg gives
    # central/vc in mg/L = ug/mL; convert to ng/mL for the IC50 scale.
    Cc <- central / vc * 1000

    # ----- Integrin turnover ODE -----------------------------------------
    # Keizer 2011 Equation (Integrin Expression in Mice section):
    #   dI/dt = kin * (1 - Emax * Cp^gamma / (IC50^gamma + Cp^gamma)) - kout * I
    # With Emax = 1 and gamma = 1 (both fixed), the inhibition collapses to
    #   inh = Cp / (IC50 + Cp).
    inh_int <- emax_int * Cc^hill_int / (ic50^hill_int + Cc^hill_int)
    d/dt(integrin) <- kin * (1 - inh_int) - kout * integrin
    integrin(0)    <- rbase

    # ----- Tumor growth ODE ----------------------------------------------
    # Keizer 2011 (Tumor Growth in Mice section):
    #   dT/dt = alpha * (1 - exp(-beta*t)) * T - effect * T
    #   effect = Emax,T * I_rel^gamma / (II50^gamma + I_rel^gamma)
    #   I_rel = (I_base - I) / I_base  (relative integrin inhibition, fraction)
    # The (1 - exp(-beta*t)) factor gates growth to zero at study start
    # (capturing the initial slow tumor growth on vehicle) and approaches
    # full growth rate alpha as t increases.
    inh_rel <- (rbase - integrin) / rbase
    tge     <- emax_tumor * inh_rel^hill_tumor /
               (ii50^hill_tumor + inh_rel^hill_tumor)
    d/dt(tumor_size) <- alpha_g * (1 - exp(-beta_g * t)) * tumor_size -
                       tge * tumor_size
    tumor_size(0)    <- rbase_tumor

    # ----- Observations --------------------------------------------------
    # Two PD outputs (integrin expression and tumor diameter) with
    # exponential / log-normal residual error per the paper's
    # "exponential residual error model" specification. The plasma
    # concentration Cc is available as a third output for diagnostic /
    # vignette use but the paper reports no residual-error model for Cc
    # in the preclinical sequential PK fit (PK was fit by NCA / pooled,
    # not population, in mice).
    integrin   ~ lnorm(expSd_integrin)
    tumor_size ~ lnorm(expSd_tumor_size)
  })
}
