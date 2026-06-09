Fauchet_2013_zidovudine <- function() {
  description <- paste(
    "One-compartment population PK model for oral zidovudine (ZDV) and",
    "its glucuronide metabolite 3'-azido-3'-deoxy-5'-glucuronylthymidine",
    "(G-ZDV) in HIV-1-infected children, infants, and adolescents",
    "(Fauchet 2013, retrospective Paris-area therapeutic-drug-monitoring",
    "cohort, n = 247, age 0.5-18 years). First-order absorption with a",
    "fixed ka = 2.86 1/h (inherited from Panhard 2007) delivers ZDV into",
    "a one-compartment central compartment with apparent total clearance",
    "CL_p/F and apparent volume V/F. The metabolite is described by a",
    "single G-ZDV state driven by a lumped metabolic formation rate",
    "constant CL_m/V_m and a first-order metabolite elimination rate",
    "constant k_el. The metabolite distribution volume V_m is not",
    "identifiable from plasma data alone and is set structurally to 1 L",
    "(same convention used by Lee 2016 for raltegravir glucuronide).",
    "Body weight enters as an estimated power-allometric covariate on",
    "CL_p/F (exponent 0.858) and on V/F (exponent 0.534), centered on",
    "the cohort median 32.2 kg; age, sex, dosage form, and antiretroviral",
    "cotreatments (3TC, ddI, ABC, LPV, RTO, NFV, NVP, EFV) were all",
    "tested and none was retained at p < 0.01.",
    sep = " "
  )
  reference <- paste(
    "Fauchet F, Treluyer JM, Frange P, Urien S, Foissac F, Bouazza N,",
    "Benaboud S, Blanche S, Hirt D. (2013).",
    "Population pharmacokinetics study of recommended zidovudine doses",
    "in HIV-1-infected children.",
    "Antimicrob Agents Chemother 57(10):4801-4808.",
    "doi:10.1128/AAC.00911-13.",
    sep = " "
  )
  vignette <- "Fauchet_2013_zidovudine"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Estimated power-allometric exponents on CL_p/F (e_wt_cl = 0.858, 95% CI 0.67-1.10) and on V/F (e_wt_vc = 0.534, 95% CI 0.22-0.81) centered on the cohort median 32.2 kg (Table 3). The paper tested the canonical fixed allometric exponents (0.75 on CL and 1.0 on V) and reports those did not improve the fit. Source column is BW (kg); rename to canonical WT before passing the dataset to rxSolve.",
      source_name        = "BW"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 247L,
    n_studies      = 1L,
    age_range      = "0.5 to 18 years",
    age_median     = "10.6 years",
    weight_range   = "6.1 to 84 kg",
    weight_median  = "32.2 kg",
    sex_female_pct = 52,
    race_ethnicity = NA_character_,
    disease_state  = "HIV-1 infection on combination antiretroviral therapy. ZDV was administered as oral syrup (33.7%) or tablet (66.3%) twice daily (94.4%) or three times daily (5.6%) for therapy with or without concomitant nucleoside reverse-transcriptase inhibitors (3TC 51%, ddI 10%, ABC 11%), protease inhibitors (LPV 39%, RTO 41%, NFV 13%), and non-nucleoside reverse-transcriptase inhibitors (NVP 18%, EFV 11%).",
    dose_range     = "Oral zidovudine, BID (94.4%) or TID (5.6%). Per-dose range 18-300 mg. Total daily dose median 600 mg (range 36-600).",
    regions        = "Greater Paris, France (retrospective therapeutic-drug-monitoring cohort 1998-2012)",
    notes          = "Demographics from Table 2 and Results 'Demographic data' paragraph. The cohort comprised 119 boys (48%) and 128 girls (52%). 782 plasma ZDV and 554 plasma G-ZDV concentrations were analysed (mean 3.6 samples per child, range 1-19). 46 ZDV (5.9%) and 13 G-ZDV (2.3%) samples were below the limit of quantification (LOQ = 0.05 mg/L for both); BQL values were replaced by LOQ/2 in the final dataset (Methods 'Population pharmacokinetic analysis' paragraph and Results 'Population pharmacokinetics' paragraph 2). The authors used this model to compare WHO and FDA dose recommendations and proposed revised doses for the 20-40 kg weight bands (Tables 1 and 4)."
  )

  ini({
    # Final parameter estimates from Fauchet 2013 Antimicrob Agents Chemother
    # 57:4801-4808, Table 3 (page 4803). NONMEM v6.2 FOCE, simultaneous
    # parent + metabolite fit. The published model parameterises everything
    # on the molar scale (the authors divide each plasma mg/L by the molar
    # mass 267.2 g/mol for ZDV and 443.3 g/mol for G-ZDV, and each mg dose
    # by 267.2 g/mol -- Methods 'Population pharmacokinetic analysis' page
    # 4802). The encoded model below keeps mg dosing and mg/L observations
    # at the boundary and applies the molar mass ratio explicitly inside
    # the metabolite ODE; the structural rate constants, allometric
    # exponents and IIV variances are the published molar-frame values.

    # Structural PK -- ZDV parent (one-compartment with first-order absorption).
    lka <- fixed(log(2.86)); label("First-order absorption rate constant ka of ZDV (1/h, FIXED at the Panhard 2007 estimate)")  # Table 3 footnote b: ka FIXED at 2.86 1/h (Panhard 2007 reference 28)
    lcl <- log(89.7);        label("Apparent ZDV total elimination clearance CL_p/F at WT = 32.2 kg (L/h)")                     # Table 3: CL_p/F = 89.7 L/h, RSE 7.1%, 95% CI 77-102
    lvc <- log(229);         label("Apparent ZDV central volume of distribution V/F at WT = 32.2 kg (L)")                       # Table 3: V/F   = 229 L,   RSE 12.4%, 95% CI 181-291

    # Structural PK -- G-ZDV metabolite (one-compartment downstream of ZDV).
    # The metabolite distribution volume V_m is not identifiable from plasma
    # alone -- the paper only estimates the lumped rate constant CL_m/V_m.
    # Inside model() the metabolite state central_gluc holds the metabolite
    # amount in a structurally-fixed V_m_gluc = 1 L frame (same convention
    # as Lee 2016 raltegravir glucuronide), so the simulated Cc_gluc in mg/L
    # is numerically equal to the central_gluc state.
    lkmet_gluc <- log(12.6); label("G-ZDV metabolic formation rate constant Cl_m/V_m in the molar frame (1/h)")                  # Table 3: Cl_m/V_m = 12.6 1/h, RSE 20.7%, 95% CI 8.2-18
    lkel_gluc  <- log(2.27); label("G-ZDV elimination rate constant k_el (1/h)")                                                  # Table 3: k_el     = 2.27 1/h, RSE 16.2%, 95% CI 1.6-3.1

    # Body-weight allometric exponents (estimated, NOT canonical 0.75 / 1).
    e_wt_cl <- 0.858; label("Body-weight power exponent on CL_p/F centered at WT = 32.2 kg (unitless)")                            # Table 3: beta_BW^(CL/F) = 0.858, RSE 11.1%, 95% CI 0.67-1.1
    e_wt_vc <- 0.534; label("Body-weight power exponent on V/F centered at WT = 32.2 kg (unitless)")                               # Table 3: beta_BW^(V/F)  = 0.534, RSE 24.5%, 95% CI 0.22-0.81

    # IIV. Table 3 footnote a defines the header omega explicitly as
    # "coefficient of variation for between-subject variability", so the
    # tabulated values 0.701, 0.807, 0.352 are CV-fractions, not log-
    # scale variances. nlmixr2 expects the log-scale variance, recovered
    # from CV via the lognormal identity omega^2 = log(1 + CV^2):
    #   var(etalcl)       = log(1 + 0.701^2) = 0.3997 (CV 70.1% -> log var)
    #   var(etalvc)       = log(1 + 0.807^2) = 0.5015 (CV 80.7% -> log var)
    #   var(etalkmet_gluc)= log(1 + 0.352^2) = 0.1168 (CV 35.2% -> log var)
    # The CL_p/F vs V/F correlation 0.733 is FIXED per Table 3 (the Corr
    # row has no RSE column). The derived log-scale covariance is
    #   cov(etalcl, etalvc) = 0.733 * sqrt(0.3997 * 0.5015) = 0.3282.
    # Table 3 omega values (CVs): CL_p/F = 0.701 (RSE 17.9%, 95% CI 0.572-0.820); V/F = 0.807 (RSE 22.9%, 95% CI 0.597-0.963); Corr(CL_p, V) = 0.733 (FIXED).
    # Log-scale: var(etalcl) = log(1 + 0.701^2) = 0.3997; var(etalvc) = log(1 + 0.807^2) = 0.5015; cov(etalcl,etalvc) = 0.733 * sqrt(0.3997*0.5015) = 0.3282.
    etalcl + etalvc ~ c(0.3997,
                        0.3282, 0.5015)
    # Table 3 omega_Cl_m/V_m (CV) = 0.352 (RSE 20.7%, 95% CI 0.213-0.473) -> log var = log(1 + 0.352^2) = 0.1168.
    etalkmet_gluc ~ 0.1168

    # Residual error. ZDV uses a combined error model with the additive
    # component held FIXED at LOQ/2 = 0.025 mg/L (Results 'Population
    # pharmacokinetics' paragraph 3 -- the paper text states "0.025
    # mg/liter, i.e., 0.1 mmol-1" but the molar conversion 0.025 mg/L /
    # 267.2 g/mol = 9.36e-5 mmol/L = 0.094 umol/L indicates the molar unit
    # in Table 3 is umol/L not mmol/L; the published mass value 0.025 mg/L
    # is used here unchanged). G-ZDV uses a proportional-only error model.
    # Paper Table 3 reports the proportional residuals as sigma values
    # (0.56 and 0.69, RSE ~ 8%) interpreted as SDs in the linear-
    # concentration space (see vignette source trace for rationale).
    addSd       <- fixed(0.025); label("ZDV additive residual SD (mg/L, FIXED at LOQ/2)")                                         # Table 3 footnote b: additive component FIXED at LOQ/2 = 0.025 mg/L = 0.1 umol/L molar
    propSd      <- 0.56;         label("ZDV proportional residual SD (fraction)")                                                  # Table 3: sigma_ZDV    = 0.56, RSE 8.5%, 95% CI 0.514-0.614
    propSd_gluc <- 0.69;         label("G-ZDV proportional residual SD (fraction)")                                                # Table 3: sigma_G-ZDV  = 0.69, RSE 7.0%, 95% CI 0.643-0.742
  })

  model({
    # Constants. Molar masses 267.2 and 443.3 g/mol come from Fauchet 2013
    # Methods 'Population pharmacokinetic analysis' paragraph 1. The
    # molar-mass ratio shows up in the metabolite ODE because the paper's
    # CL_m/V_m rate constant is in molar units; converting it to a
    # mass-units formation rate of metabolite (mg G-ZDV per L per hour)
    # multiplies by MW_GZDV / MW_ZDV (= 1.659 mg G-ZDV per mg ZDV
    # consumed, stoichiometry 1:1).
    mw_zdv  <- 267.2
    mw_gzdv <- 443.3
    # Lumped metabolite volume (not identifiable in the published model;
    # set structurally to 1 L so the central_gluc state value in mg
    # equals the metabolite concentration in mg/L numerically).
    vc_gluc <- 1

    # Individual parameters. WT allometric covariate is applied as the
    # estimated power exponents on (WT / 32.2)^e_wt_cl and (WT / 32.2)^e_wt_vc
    # per Methods 'Population pharmacokinetic analysis' final paragraph.
    ka        <- exp(lka)
    cl        <- exp(lcl + etalcl) * (WT / 32.2) ^ e_wt_cl
    vc        <- exp(lvc + etalvc) * (WT / 32.2) ^ e_wt_vc
    kmet_gluc <- exp(lkmet_gluc + etalkmet_gluc)
    kel_gluc  <- exp(lkel_gluc)

    # Parent elimination micro-constant. The paper's CL_p/F is total
    # apparent ZDV elimination clearance (includes both metabolic and
    # renal routes); the metabolite formation rate kmet_gluc is an
    # independently parameterised rate constant (Karlsson-Sheiner-style
    # non-mass-balance parent-metabolite model), not derived from CL_p/F.
    kel <- cl / vc

    # ODE system.
    # Parent ZDV: depot (mg) -> central (mg); kel = CL_p/F / V/F.
    # Metabolite G-ZDV: formed from parent at rate (MW_GZDV/MW_ZDV) * kmet
    # * C_parent (mg per L per hour) into a V_m = 1 L compartment, then
    # eliminated first-order at rate kel_gluc.
    d/dt(depot)        <- -ka * depot
    d/dt(central)      <-  ka * depot - kel * central
    d/dt(central_gluc) <-  (mw_gzdv / mw_zdv) * kmet_gluc * (central / vc) * vc_gluc -
                           kel_gluc * central_gluc

    # Observations (mg/L). Parent: central amount / parent volume.
    # Metabolite: central_gluc state / lumped V_m_gluc = 1 L.
    Cc      <- central / vc
    Cc_gluc <- central_gluc / vc_gluc

    Cc      ~ add(addSd) + prop(propSd)
    Cc_gluc ~ prop(propSd_gluc)
  })
}
