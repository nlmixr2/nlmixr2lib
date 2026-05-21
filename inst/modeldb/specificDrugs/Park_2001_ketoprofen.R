Park_2001_ketoprofen <- function() {
  description <- "One-compartment oral PK plus Holford-Sheiner effect-compartment for synovial fluid disposition of ketoprofen in adults with arthritis at steady state on 100 mg oral twice-daily dosing (Park 2001 Tables 2-3, Eq. 1; effect-compartment elimination rate keo = 0.16 1/h, peak synovial:plasma ratio 0.77 with 3.1 h time lag)."
  reference <- "Park JY, Sohn JH, Yoon YR, Shon JH, Cha IJ, Seo SS, Choi JS, Shin JG. Disposition Kinetics of Ketoprofen into Synovial Fluid Following Systemic Administration: Population Pharmacokinetic Analysis. Kor J Clin Pharmacol Ther. 2001;9(1):97-107. doi:10.12793/jkscpt.2001.9.1.97"
  vignette <- "Park_2001_ketoprofen"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 17L,                                # Park 2001 Methods + Table 1
    n_studies      = 1L,                                 # Single-centre study at Inje University Pusan Paik Hospital
    age_range      = "22-63 years",                      # Park 2001 Table 1
    age_median     = "44.2 years (mean, SD 13.3)",       # Park 2001 Table 1
    weight_range   = "46-75 kg",                         # Park 2001 Table 1
    weight_median  = "62.2 kg (mean, SD 8.8)",           # Park 2001 Table 1
    sex_female_pct = 52.94,                              # Park 2001 Table 1: 8 male / 9 female
    race_ethnicity = NULL,                               # Not reported (Korean single-centre study)
    disease_state  = "Arthritis: 7 rheumatoid arthritis, 10 osteoarthritis; all enrolled with normal blood chemistry and urinalysis screens.",
    dose_range     = "100 mg oral ketoprofen twice daily for at least 4 days to reach steady state.",
    regions        = "Korea (single centre, Inje University Pusan Paik Hospital).",
    n_observations = "12 of 17 patients had full plasma sampling at 0 (pre-dose), 0.5, 1, 2, 3, 4, 5, 6, 8, and 12 h post-dose; all 17 patients had 1-4 synovial fluid samples drawn at varied steady-state times (Methods).",
    notes          = "Plasma PK parameters (ka, kel, Vd/F, CL/F) were estimated individually with WinNONLIN from the 12 patients with full plasma sampling and reported as arithmetic mean +/- SD across individuals (Park 2001 Table 2). The synovial-fluid effect-compartment elimination rate constant keo was the only population parameter estimated with NONMEM V level 1.1 using all 17 patients' synovial concentrations (Park 2001 Table 3, Eq. 1-3). The Holford-Sheiner steady-state assumption k1e0 = keo (Methods, Fig. 1) collapses synovial influx and efflux into a single rate constant."
  )

  ini({
    # Structural PK parameters - reference body weight 62.2 kg (Table 1 mean)
    # Reported values in Park 2001 Table 2 are per-kg apparent quantities; absolute
    # values used here are the reported means multiplied by the mean body weight.
    lka  <- log(0.80);  label("Absorption rate constant ka (1/h)")           # Park 2001 Table 2: ka = 0.80 +/- 0.69 1/h (WinNONLIN individual fits, n=12)
    lcl  <- log(12.44); label("Apparent oral clearance CL/F (L/h)")          # Park 2001 Table 2: CL/F = 0.20 +/- 0.13 L/kg/h * 62.2 kg = 12.44 L/h
    lvc  <- log(24.26); label("Apparent oral central volume Vd/F (L)")       # Park 2001 Table 2: Vd/F = 0.39 +/- 0.25 L/kg * 62.2 kg = 24.26 L
    lkeo <- log(0.16);  label("Synovial fluid equilibration rate keo (1/h)") # Park 2001 Table 3: keo = 0.16 1/h (95% CI 0.14-0.19); NONMEM population estimate

    # IIV - log-normal; omega^2 = log(CV^2 + 1).
    # ka, CL/F, Vd/F CVs are derived from the Table 2 SD/mean ratios across the
    # 12 WinNONLIN individual fits; keo CV is the reported NONMEM CV(keo) = 24.5%
    # from Table 3, originally fit with a proportional (1+eta) error model in Park
    # 2001 Eq. 2 and converted here to the equivalent log-normal variance.
    etalka  ~ 0.5562   # CV(ka)  = 0.69/0.80 = 86.3 percent -> omega^2 = log(1 + 0.8625^2) = 0.5562 (Park 2001 Table 2)
    etalcl  ~ 0.3525   # CV(CL)  = 0.13/0.20 = 65.0 percent -> omega^2 = log(1 + 0.65^2)   = 0.3525 (Park 2001 Table 2)
    etalvc  ~ 0.3443   # CV(Vc)  = 0.25/0.39 = 64.1 percent -> omega^2 = log(1 + 0.641^2)  = 0.3443 (Park 2001 Table 2)
    etalkeo ~ 0.0583   # CV(keo) = 24.5 percent             -> omega^2 = log(1 + 0.245^2) = 0.0583 (Park 2001 Table 3; paper used (1+eta) proportional eta model)

    # Residual error
    # Park 2001 NONMEM fit (Eq. 3) reports CV_sigma = 40% on the synovial fluid
    # observations; this maps directly to propSd_Csf. The paper does NOT report a
    # population residual error for plasma -- plasma PK was estimated individually
    # by WinNONLIN with no population residual-error term. propSd is assigned a
    # plausible value (0.15) so a popPK simulation produces reasonable plasma
    # noise; see vignette Assumptions and deviations section.
    propSd     <- 0.15; label("Proportional residual error on plasma Cc (fraction)")          # not reported in Park 2001; assumed 0.15 for plasma popPK simulation
    propSd_Csf <- 0.40; label("Proportional residual error on synovial fluid Csf (fraction)") # Park 2001 Table 3: CV_sigma = 40% (proportional residual-error model, Eq. 3)
  })
  model({
    # Individual PK parameters
    ka  <- exp(lka  + etalka)
    cl  <- exp(lcl  + etalcl)
    vc  <- exp(lvc  + etalvc)
    keo <- exp(lkeo + etalkeo)

    # Micro-rate constant
    kel <- cl / vc

    # One-compartment oral PK
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central
    Cc            <- central / vc

    # Effect compartment representing synovial fluid.
    # Park 2001 Eq. 1: dCsf/dt = k1e * Cp - keo * Csf.
    # Holford-Sheiner steady-state mass-balance assumption k1e = keo
    # (Methods, after Eq. 1) gives the standard effect-compartment form:
    d/dt(effect) <- keo * (Cc - effect)
    Csf          <- effect

    # Observations (Park 2001 Eq. 3 for synovial; plasma error assumed -- see ini())
    Cc  ~ prop(propSd)
    Csf ~ prop(propSd_Csf)
  })
}
