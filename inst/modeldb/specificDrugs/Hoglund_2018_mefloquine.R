# Population pharmacokinetic model of mefloquine in Burmese adults with
# acute uncomplicated Plasmodium falciparum malaria treated with the
# standard 3-day artesunate-mefloquine combination (Hoglund 2018,
# Malaria Journal 17:322; doi:10.1186/s12936-018-2466-3).

Hoglund_2018_mefloquine <- function() {
  description <- paste(
    "Population PK model for oral mefloquine in Burmese adults with",
    "uncomplicated Plasmodium falciparum malaria treated with the",
    "standard 3-day artesunate-mefloquine combination (Hoglund 2018).",
    "One-transit-compartment absorption with ka = ktr feeds a",
    "two-compartment disposition model. No covariates were retained in",
    "the final model: body-weight allometric scaling (fixed exponents",
    "0.75 / 1.0), sex, admission parasitaemia, and validated molecular",
    "markers of mefloquine and artemisinin resistance (pfmdr1, pfcrt,",
    "atp6, pfk13) were tested in a step-wise covariate search but did",
    "not significantly improve the model. Relative bioavailability F is",
    "implicitly 1 (the paper tested adding an estimated F with IIV and",
    "excluded it from the final model). NONMEM additive residual error",
    "on the log-transformed observation is encoded here as a",
    "proportional residual in the linear concentration space.",
    sep = " "
  )
  reference <- paste(
    "Hoglund RM, Ruengweerayut R, Na-Bangchang K (2018).",
    "Population pharmacokinetics of mefloquine given as a 3-day",
    "artesunate-mefloquine in patients with acute uncomplicated",
    "Plasmodium falciparum malaria in a multidrug-resistant area along",
    "the Thai-Myanmar border. Malaria Journal 17:322.",
    "doi:10.1186/s12936-018-2466-3.",
    sep = " "
  )
  vignette <- "Hoglund_2018_mefloquine"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  population <- list(
    species         = "human",
    n_subjects      = 129L,
    n_studies       = 1L,
    n_cured         = 93L,
    n_recrudescent  = 36L,
    age_range       = "16-50 years (Table 1)",
    age_median      = "25 years (Table 1)",
    weight_range    = "39-73.5 kg (Table 1)",
    weight_median   = "52.5 kg (Table 1)",
    sex_female_pct  = 49.6,
    disease_state   = paste(
      "Acute uncomplicated Plasmodium falciparum malaria in adults",
      "(Burmese migrant population, all aged over 15 years) in a",
      "multidrug-resistant transmission area. Pregnant and",
      "breast-feeding women were excluded; severe or complicated",
      "malaria and severe malnutrition were exclusion criteria",
      "(Methods, Patients and treatment)."
    ),
    dose_range      = paste(
      "Standard 3-day artesunate-mefloquine combination. Mefloquine",
      "was administered on days 0 and 1 only: day 0 mefloquine",
      "15 mg/kg (given as a fixed 750 mg dose = 3 tablets of 250 mg",
      "mefloquine; Atlantic Pharmaceutical Company, Thailand) plus",
      "artesunate 4 mg/kg (200 mg); day 1 mefloquine 10 mg/kg (given",
      "as a fixed 500 mg dose = 2 tablets of 250 mg mefloquine) plus",
      "artesunate 4 mg/kg (200 mg); day 2 artesunate 4 mg/kg plus",
      "primaquine 0.6 mg/kg (no mefloquine). All doses observed; if",
      "patients vomited within 30 min, the dose was repeated."
    ),
    regions         = "Thai-Myanmar border (Mae Tao clinic for migrant workers, Tak Province, Thailand; March 2008 to February 2009)",
    parasitaemia_range = "1260-84,000 parasites/uL at admission (median 5320; Table 1)",
    notes           = paste(
      "Demographics from Hoglund 2018 Table 1. Total observations 653",
      "post-dose mefloquine measurements; whole-blood matrix",
      "quantified by HPLC (LOQ 2 ng/mL, LOD 0.5 ng/mL). Baseline",
      "mefloquine concentrations were present in 5 patients (3.88%)",
      "and were ignored. 36 of 129 patients (27.9%) had recrudescent",
      "infection on days 21-40 of the 42-day follow-up; 93 had",
      "successful treatment. Genotyping (pfmdr1, pfcrt, atp6, pfk13",
      "propeller domain codons 440-680) was performed but the",
      "molecular-resistance markers were not retained as PK",
      "covariates (Results, Pharmacokinetic model)."
    )
  )

  ini({
    # Structural parameters from Hoglund 2018 Table 2 ("Parameter value
    # (%RSE)" column). Apparent values relative to F = 1 reported on
    # the linear scale; log() applied here for the nlmixr2 internal
    # log-scale. No covariates were retained in the final model: body
    # weight allometric scaling (fixed exponents 0.75 / 1.0), sex,
    # admission parasitaemia, and validated molecular markers of
    # mefloquine and artemisinin resistance (pfmdr1, pfcrt, atp6,
    # pfk13) were all tested in a step-wise covariate search and
    # excluded (Results, Pharmacokinetic model).
    lcl  <- log(2.77)
    label("Apparent mefloquine elimination clearance CL/F (L/h)")
    # Hoglund 2018 Table 2: CL/F = 2.77 L/h (RSE 4.81%; 95% CI 2.52-3.04)

    lvc  <- log(359)
    label("Apparent central volume of distribution Vc/F (L)")
    # Hoglund 2018 Table 2: Vc/F = 359 L (RSE 3.38%; 95% CI 335-384)

    lq   <- log(11.7)
    label("Apparent inter-compartmental clearance Q/F (L/h)")
    # Hoglund 2018 Table 2: Q/F = 11.7 L/h (RSE 8.37%; 95% CI 10.1-13.9)

    lvp  <- log(474)
    label("Apparent peripheral volume of distribution Vp/F (L)")
    # Hoglund 2018 Table 2: Vp/F = 474 L (RSE 8.02%; 95% CI 406-554)

    lmtt <- log(3.89)
    label("Mean transit time of the 1-transit-compartment absorption chain MTT (h)")
    # Hoglund 2018 Table 2: MTT = 3.89 h (RSE 5.18%; 95% CI 3.54-4.32).
    # Number of transit compartments = 1 (fixed; Results, Pharmacokinetic
    # model: "The absorption phase was described by 1 transit
    # compartments"). The final model holds ka = ktr (Results,
    # Pharmacokinetic model: "Estimating the transit rate constant and
    # the absorption rate constant separately resulted in a significantly
    # improved model. However, it also resulted in a model with a low
    # precision in the estimation of the apparent volume of distribution
    # of the central compartment (a relative standard error of 46.9%)").
    # The absorption chain depot -> transit1 -> central has
    # NN + 1 = 2 equal-rate transitions; ktr = (NN + 1) / MTT = 2 / MTT
    # (Savic & Karlsson 2007 convention; same idiom as the companion
    # Hoglund_2012_piperaquine.R and Hoglund_2017_piperaquine.R files).

    # IIV. Hoglund 2018 Table 2 "IIV %CV (%RSE)" column. The reported
    # %CV is the coefficient of variation of the log-normal
    # individual-parameter distribution; the internal log-scale
    # variance is recovered as omega^2 = log(CV^2 + 1) (Methods, Eq 1:
    # P_i = theta_P * exp(eta_i,P), exponential IIV form).
    #
    #   CL  IIV 38.0%   -> omega^2 = log(0.380^2 + 1) = 0.134880
    #   Vc  IIV "-"     no IIV reported in Table 2; omitted
    #   MTT IIV 43.6%   -> omega^2 = log(0.436^2 + 1) = 0.174034
    #   Q   IIV 55.0%   -> omega^2 = log(0.550^2 + 1) = 0.264285
    #   Vp  IIV 63.1%   -> omega^2 = log(0.631^2 + 1) = 0.335158
    etalcl  ~ 0.134880
    # Hoglund 2018 Table 2: IIV on CL/F = 38.0% CV (RSE 14.8%; 95% CI 27.1-48.0; shrinkage 33.7%)

    etalmtt ~ 0.174034
    # Hoglund 2018 Table 2: IIV on MTT = 43.6% CV (RSE 11.5%; 95% CI 33.7-53.0; shrinkage 18.5%)

    etalq   ~ 0.264285
    # Hoglund 2018 Table 2: IIV on Q/F = 55.0% CV (RSE 27.3%; 95% CI 22.4-80.5; shrinkage 50.8%)

    etalvp  ~ 0.335158
    # Hoglund 2018 Table 2: IIV on Vp/F = 63.1% CV (RSE 15.5%; 95% CI 43.9-81.7; shrinkage 43.2%)

    # Residual error. Methods, Pharmacokinetic analysis: "The natural
    # logarithm of quantified mefloquine concentrations was analysed".
    # The NONMEM additive-on-log-scale residual maps to a nlmixr2
    # proportional residual in linear concentration space (see
    # references/parameter-names.md Residual error and the matching
    # comment in Hoglund_2012_piperaquine.R / Hoglund_2017_piperaquine.R).
    # Table 2 reports the log-scale variance sigma^2 = 0.0902; the SD
    # is sqrt(0.0902) = 0.300333, which equals the proportional CV to
    # first order.
    propSd <- sqrt(0.0902)
    label("Proportional residual SD for mefloquine whole-blood concentration (SD on log scale)")
    # Hoglund 2018 Table 2: sigma^2 = 0.0902 (variance on log scale, RSE 18.7%; 95% CI 0.0591-0.124; shrinkage 17.9%)
  })

  model({
    # Individual PK parameters. Hoglund 2018 final model retained no
    # covariates; population log-parameter + IIV only. Vc carries no
    # IIV (Table 2 IIV column is "-" for Vc/F).
    cl  <- exp(lcl  + etalcl)
    vc  <- exp(lvc)
    q   <- exp(lq   + etalq)
    vp  <- exp(lvp  + etalvp)

    # Mean transit time and chain rate constant. Hoglund 2018 final
    # absorption model uses NN = 1 transit compartment with ka = ktr.
    # The absorption chain depot -> transit1 -> central has
    # NN + 1 = 2 equal-rate transitions; ktr = (NN + 1) / MTT = 2 / MTT
    # (Savic & Karlsson 2007 convention).
    mtt <- exp(lmtt + etalmtt)
    ktr <- 2 / mtt

    # Two-compartment disposition micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system: 1-transit-compartment absorption feeding a
    # two-compartment disposition model. The same ktr propagates the
    # dose through depot and transit1 into central.
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1
    d/dt(central)     <-  ktr * transit1 - kel * central -
                          k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central  - k21 * peripheral1

    # Mefloquine whole-blood concentration. Dose in mg of mefloquine
    # base (tablet count of 250 mg mefloquine), Vc in L -> central / vc
    # has units mg/L. Multiplying by 1000 converts to ng/mL for direct
    # comparison against Hoglund 2018 Table 3 (Cmax in ng/mL) and the
    # paper's reported LOQ of 2 ng/mL. Predictions correspond to whole
    # blood (Methods, Drug quantification).
    Cc <- 1000 * central / vc

    # Proportional residual error on the linear-concentration scale
    # (NONMEM additive-on-log-scale -> nlmixr2 proportional; see
    # ini() comment on propSd).
    Cc ~ prop(propSd)
  })
}
