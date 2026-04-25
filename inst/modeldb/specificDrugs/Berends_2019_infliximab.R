Berends_2019_infliximab <- function() {
  description <- "Two-compartment TMDD-QSS population PK/target-dynamics model of infliximab and free TNF in adults with moderate-to-severe ulcerative colitis (Berends 2019)"
  reference <- "Berends SE, Strik AS, Van Selm S, Lowenberg M, Ponsioen CY, D'Haens GR, Mathot RAA. Tumor necrosis factor-mediated disposition of infliximab in ulcerative colitis patients. J Pharmacokinet Pharmacodyn. 2019;46(6):543-551. doi:10.1007/s10928-019-09652-5"
  vignette <- "Berends_2019_infliximab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL using median-normalization (ALB/38)^e_alb_cl per Berends 2019 covariate equation P = P_TV * (COV/COV_median)^theta. Reference is the cohort-median albumin of 38 g/L (Table 1).",
      source_name        = "ALB"
    ),
    ADA_POS = list(
      description        = "Anti-drug antibody (antibodies-to-infliximab) positivity",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ADA-negative)",
      notes              = "Multiplicative effect on CL using P = P_TV * theta^ADA per Berends 2019 categorical covariate equation; CL is multiplied by 2.15 for ADA-positive subjects. ADA-status was detected in 7/20 patients during follow-up. Source paper text uses 'antibodies-to-infliximab'; renamed to canonical ADA_POS per inst/references/covariate-columns.md.",
      source_name        = "ADA"
    )
  )

  population <- list(
    n_subjects     = 20L,
    n_studies      = 1L,
    age_range      = "19-69 years (median 36)",
    weight_range   = "47-90 kg (median 70)",
    sex_female_pct = 35,
    race_ethnicity = "Not reported (single-center Netherlands cohort)",
    disease_state  = "Adults with moderate-to-severe ulcerative colitis (95% endoscopic Mayo score 3, 95% corticosteroid-refractory, 65% pancolitis, 35% hospitalized).",
    dose_range     = "5 mg/kg IV infliximab induction at weeks 0, 2, and 6 (one subject also received an additional dose at day 5).",
    regions        = "Single-center prospective cohort, Amsterdam, Netherlands.",
    albumin_median = "38 g/L (range 23-45)",
    crp_median     = "25.3 mg/L (range 0.6-196.2)",
    concomitant_thiopurines_pct = 55,
    sccai_median   = "10 (range 1-15)",
    notes          = "Anti-TNF naive patients with active disease at baseline. The dataset comprises 214 IFX serum and 214 TNF serum measurements; antibodies-to-infliximab were detected in 7/20 patients during follow-up. Quantified IFX with a homogenous mobility shift assay (LLOQ 0.6 ug/mL, CV 12%); free TNF with an ultrasensitive Singulex immunoassay (LLOQ 10 fg/mL, CV 15%)."
  )

  ini({
    # ---- Structural PK (Table 2 final-model estimates) -----------------------
    lcl <- log(0.404); label("Population clearance CL (L/day)")                     # Berends 2019 Table 2
    lvc <- log(3.18);  label("Population central volume of distribution Vc (L)")    # Berends 2019 Table 2
    lvp <- log(1.64);  label("Population peripheral volume of distribution Vp (L)") # Berends 2019 Table 2
    lq  <- log(0.344); label("Population intercompartmental clearance Q (L/day)")   # Berends 2019 Table 2

    # ---- Covariate effects (Table 2 final-model estimates) -------------------
    # ADA on CL: categorical, P = P_TV * theta^ADA (Methods "Categorical
    # covariates"); CL is multiplied by e_ada_cl when ADA = 1.
    # Albumin on CL: continuous, P = P_TV * (ALB/38)^e_alb_cl (Methods
    # "Continuous covariates"); 38 g/L is the cohort median (Table 1).
    e_ada_cl <-  2.15;  label("Multiplicative factor on CL for ADA-positive subjects (unitless)") # Berends 2019 Table 2 (ADA-CL)
    e_alb_cl <- -1.13;  label("Power exponent on (ALB/38) for CL (unitless)")                     # Berends 2019 Table 2 (Alb-CL)

    # ---- TMDD-QSS target parameters (Table 2 / Eq. 8-13) ---------------------
    # Paper reports Kss in nM and Bmax in pM. Converted here to IFX-mass-
    # equivalent ug/mL using IFX MW = 149 kDa so that all drug-target QSS
    # arithmetic runs on a single ug/mL scale consistent with mg dosing and
    # L volumes (Cc = central/vc is natively mg/L = ug/mL):
    #   Kss [ug/mL] = Kss [nM] * MW_IFX_kDa / 1000
    #               = 13.6 * 149/1000 = 2.0264 ug/mL
    #   Bmax [ug/mL] = Bmax [nM] * MW_IFX_kDa / 1000
    #                = 3.80e-4 * 149/1000 = 5.662e-5 ug/mL
    # Internalization rate is paper's ke(P) (= kint in QSS notation). kdeg
    # fixed per sensitivity analysis (Suppl. Table 1; Results "Simulations").
    lkss  <- log(2.0264);   label("Steady-state dissociation constant Kss in IFX-equivalent ug/mL (paper: 13.6 nM)")           # Berends 2019 Results "Final model" + Abstract (text 13.6 nM); Table 2 bootstrap 13.7 nM; Table 2 "Estimate" displays integer-rounded 14
    lBmax <- log(5.662e-5); label("Baseline TNF concentration Bmax in IFX-equivalent ug/mL (paper: 0.38 pM = 19.8 pg/mL)")     # Berends 2019 Table 2 (Bmax 0.38 pM = 19.8 pg/mL)
    lkint <- log(0.984);    label("Complex internalization rate ke(P) = kint (1/day)")                                          # Berends 2019 Table 2 (ke(P))
    lkdeg <- fixed(log(5.12)); label("TNF degradation rate kdeg (1/day; FIXED per sensitivity analysis)")                       # Berends 2019 Table 2 (kdeg fixed = 5.12; Results "Simulations")

    # ---- Inter-individual variability ---------------------------------------
    # CV%-style values from Table 2 converted to log-normal variance via
    # omega^2 = log(CV^2 + 1):
    #   29.2 % -> 0.08185 (CL)
    #   22.7 % -> 0.05022 (Vc)
    #   74.2 % -> 0.43847 (Vp)
    #   39.2 % -> 0.14302 (Bmax)
    # CL-Vc are reported as correlated; the off-diagonal in Table 2
    # ("Cov. CL-Vc (%)" = 12.3, 95% CI 0-21.7) is interpreted as the
    # correlation coefficient (rho = 0.123) since reading 12.3 directly as a
    # covariance gives a correlation > 1 (mathematically impossible) and
    # reading it as omega-scale covariance similarly violates Cauchy-Schwarz.
    # With rho = 0.123, the off-diagonal covariance is
    # 0.123 * sqrt(0.08185 * 0.05022) = 0.007884.
    etalcl + etalvc ~ c(0.08185,
                        0.007884, 0.05022)  # Berends 2019 Table 2: IIV CL 29.2%, IIV Vc 22.7%, Cov 12.3% (interpreted as correlation 0.123)
    etalvp   ~ 0.43847    # Berends 2019 Table 2 IIV Vp 74.2 % CV
    etalBmax ~ 0.14302    # Berends 2019 Table 2 IIV Bmax 39.2 % CV

    # ---- Residual error (Table 2 final-model estimates) ---------------------
    # Both errors reported as proportional. IFX residual applies to Cc (ug/mL);
    # TNF residual applies to total TNF concentration (Rtot, IFX-eq ug/mL).
    # Proportional errors are unit-agnostic so they apply directly across the
    # ug/mL <-> pg/mL <-> nM display conversions.
    CcpropSd   <- 0.210; label("Proportional residual error for IFX serum concentration (fraction)") # Berends 2019 Table 2
    RtotpropSd <- 0.406; label("Proportional residual error for total TNF concentration (fraction)") # Berends 2019 Table 2
  })

  model({
    # Molecular weights used for input/output unit conversions only (not in
    # the QSS arithmetic, which runs entirely in IFX-equivalent ug/mL).
    # Source: paper Methods "Model development": "IFX and TNF concentrations
    # were converted to nanomolar using their molecular weights of 149 kDa
    # (IFX) and 52 kDa (TNF)".
    MW_IFX_kDa <- 149  # infliximab MW, kDa
    MW_TNF_kDa <- 52   # TNF MW, kDa

    # ---- Individual PK parameters ------------------------------------------
    # Covariate forms per Berends 2019 Methods:
    #   ADA on CL: CL = CL_TV * e_ada_cl^ADA (categorical)
    #   ALB on CL: CL = CL_TV * (ALB / 38)^e_alb_cl (median-normalized power)
    cl <- exp(lcl + etalcl) * (e_ada_cl^ADA_POS) * (ALB / 38)^e_alb_cl
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q  <- exp(lq)

    # ---- Individual TMDD-QSS parameters ------------------------------------
    kss  <- exp(lkss)                  # ug/mL, IFX-mass-equivalent
    Bmax <- exp(lBmax + etalBmax)      # ug/mL, IFX-mass-equivalent (baseline TNF)
    kint <- exp(lkint)                 # 1/day (= ke(P) in paper)
    kdeg <- exp(lkdeg)                 # 1/day (FIXED at 5.12)

    # TNF synthesis rate from steady-state constraint on Rtot at baseline:
    # dRtot/dt|t=0 = 0  ->  ksyn = kdeg * Bmax  (ug/mL/day)
    ksyn <- kdeg * Bmax

    # ---- Initial condition for endogenous TNF state ------------------------
    total_target(0) <- Bmax            # ug/mL, IFX-mass-equivalent

    # ---- QSS algebraic relations (Eq. 9-10) --------------------------------
    # Total drug central concentration ctot = Atot/Vc (ug/mL).
    # Free drug Cfree from the QSS quadratic (Eq. 10):
    #   C = 0.5 * [(Ctot - Rtot - Kss) + sqrt((Ctot - Rtot - Kss)^2 + 4*Kss*Ctot)]
    # Complex concentration: P = Rtot * Cfree / (Kss + Cfree).
    ctot    <- central / vc
    disc    <- ctot - total_target - kss
    cfree   <- 0.5 * (disc + sqrt(disc * disc + 4 * kss * ctot))
    complex <- total_target * cfree / (kss + cfree)

    # ---- ODE system (Eq. 11-13) --------------------------------------------
    # Drug central (mg/day): linear elimination on free drug, distribution on
    # free drug, and TMDD internalization on the bound complex.
    #   cl[L/day] * cfree[mg/L] = mg/day; kint[1/day] * vc[L] * complex[mg/L] = mg/day
    # Drug peripheral (mg/day): standard 2-compartment exchange on free drug.
    # Total TNF (ug/mL/day): zero-order synthesis, first-order degradation,
    # and complex-mediated internalization (kint - kdeg) per Eq. 13.
    d/dt(central)      <- -cl * cfree - q * cfree + (q / vp) * peripheral1 -
                          kint * vc * complex
    d/dt(peripheral1)  <-  q * cfree - (q / vp) * peripheral1
    d/dt(total_target) <-  ksyn - kdeg * total_target -
                          (kint - kdeg) * (total_target * cfree / (kss + cfree))

    # ---- Observation variables ---------------------------------------------
    # Cc = total IFX serum concentration in ug/mL (the bioassay measures total
    # drug per Methods "Serum measurements" / Eq. 9 C_tot = Ac/Vc).
    Cc <- ctot

    # Rtot = total TNF serum concentration in IFX-equivalent ug/mL (matches
    # Eq. 13's state). Display conversions:
    #   Rtot [nM]    = Rtot [ug/mL] * 1000 / MW_IFX_kDa
    #   Rtot [pg/mL] = Rtot [nM]    * MW_TNF_kDa * 1000
    #                = Rtot [ug/mL] * MW_TNF_kDa / MW_IFX_kDa * 1e6
    # Spot-check: Bmax = 5.662e-5 * 52/149 * 1e6 = 19.76 pg/mL (paper 19.8).
    Rtot       <- total_target
    Cc_nM      <- ctot         * 1000 / MW_IFX_kDa                          # ug/mL -> nM (drug)
    Rtot_nM    <- total_target * 1000 / MW_IFX_kDa                          # ug/mL -> nM (target on the IFX-equivalent molar scale, identical to the paper's nM TNF since the QSS imposes 1:1 stoichiometry)
    Rtot_pgml  <- total_target * MW_TNF_kDa / MW_IFX_kDa * 1e6              # ug/mL -> pg/mL via TNF MW
    free_TNF_pgml <- (total_target - complex) * MW_TNF_kDa / MW_IFX_kDa * 1e6

    # ---- Residual-error models ---------------------------------------------
    Cc   ~ prop(CcpropSd)
    Rtot ~ prop(RtotpropSd)
  })
}
