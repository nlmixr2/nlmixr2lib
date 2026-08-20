Kado_2020_benzathine_benzylpenicillin_g <- function() {
  description <- "One-compartment population PK model for penicillin released from benzathine benzylpenicillin G (Bicillin L-A) with three parallel absorption pathways (slow and fast via a transit compartment, plus an immediate pathway that bypasses the transit compartment) and route-specific structural parameters for intramuscular (IM) and subcutaneous (SC) administration, developed from 311 dried-blood-spot penicillin concentrations in a randomized crossover of 15 healthy adult male volunteers each receiving 1.2 MIU IM and 1.2 MIU SC into the dorsogluteal region (Kado 2020)."
  reference   <- paste(
    "Kado JH, Salman S, Henderson R, Hand R, Wyber R, Page-Sharp M, Batty K,",
    "Carapetis J, Manning L. Subcutaneous administration of benzathine",
    "benzylpenicillin G has favourable pharmacokinetic characteristics for the",
    "prevention of rheumatic heart disease compared with intramuscular",
    "injection: a randomized, crossover, population pharmacokinetic study in",
    "healthy adult volunteers. J Antimicrob Chemother. 2020;75(10):2986-2993.",
    "doi:10.1093/jac/dkaa282",
    sep = " "
  )
  vignette    <- "Kado_2020_benzathine_benzylpenicillin_g"
  units       <- list(time = "day", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    depot1   = list(analyte = "benzathine benzylpenicillin g", units = "mg", specimen = "administration site", verified = FALSE),
    transit1 = list(analyte = "benzathine benzylpenicillin g", units = "mg", specimen = "administration site", verified = FALSE),
    depot2   = list(analyte = "benzathine benzylpenicillin g", units = "mg", specimen = "administration site", verified = FALSE),
    transit2 = list(analyte = "benzathine benzylpenicillin g", units = "mg", specimen = "administration site", verified = FALSE),
    depot3   = list(analyte = "benzathine benzylpenicillin g", units = "mg", specimen = "administration site", verified = FALSE),
    central  = list(analyte = "benzathine benzylpenicillin g", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    ROUTE_SC = list(
      description        = "Route-of-administration indicator: 1 = SC injection, 0 = IM injection.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (IM injection; the reference route for structural parameters).",
      notes              = "Per-dose-record covariate. IM is the reference route (Kado 2020 Table 1 IM structural parameters are the paper's reference); SC structural parameters were estimated as multiplicative factors on the IM values ('Relative structural model parameters for SC administration' rows). ROUTE_SC switches the three absorption half-lives (t1/2,abs-1/-2/-3), the transit half-life (t1/2,tr), the two dose-split ratios (RAT-transit, RAT-slowfast), and applies the SC relative bioavailability F_SC = 0.957. IIVs are also route-specific per Kado 2020 Methods 'PK modelling' section.",
      source_name        = "ROUTE_SC"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 15L,
    n_studies      = 1L,
    age_range      = "18.3-46.6 years (median 24.6)",
    age_median     = "24.6 years",
    weight_range   = "not reported (BMI 19.2-25.8 kg/m^2; entry criterion 18.5-26.0 kg/m^2)",
    weight_median  = "not reported (assumed ~70 kg for reference simulations per Kado 2020 Methods 'Simulations')",
    sex_female_pct = 0,
    race_ethnicity = c(Caucasian = 66.7, Asian = 26.7, African = 6.7),
    disease_state  = "healthy adult male volunteers (nonsmokers, no chronic illness, no allergy to penicillin/cephalosporins, no prescription/OTC/herbal medications for at least 7 days prior)",
    dose_range     = "1.2 MIU (1016.6 mg / 2.3 mL Bicillin L-A) single dose IM followed by 1.2 MIU single dose SC, or vice-versa, with a 10-week washout between doses",
    regions        = "Australia (Perth, Western Australia; Linear Clinical Research volunteer database)",
    notes          = "Randomized crossover pharmacokinetic study; each participant received both routes. 311 valid dried-blood-spot penicillin concentrations were included in the final analysis (5.3% BLQ). 53 paired plasma/DBS samples at 12 h and 14 days confirmed DBS validity (Spearman r = 0.99). Sampling: baseline, 2, 6, 12, 24, 48 h then 3, 5, 7, 14, 21, 28, 42 days post-dose. Baseline demographics: Table (Results, para 1). Final model estimates: Table 1. Re-parameterized route comparison: Table 2. Simulation-derived Cmax/Cmin/T>20/T>10: Table 3."
  )

  ini({
    # ==========================================================================
    # Elimination (common to both routes)
    # ==========================================================================
    # Paper reports kel = 1.32 h^-1 for a 70 kg reference individual, fixed
    # from the previously published Broderick / Kado analysis (Kado 2020
    # Methods 'PK modelling' -- 'Initial modelling using the previously
    # published structure (including fixing the elimination rate)').
    # Converted to 1/day (paper's absorption half-lives are in days).
    lkel <- fixed(log(1.32 * 24));  label("First-order elimination rate constant for penicillin (1/day)")  # Table 1: kel = 1.32 h^-1 fixed = 31.68 1/day
    lvc  <- log(39.6);              label("Apparent central volume of distribution V/F (L, 70 kg reference)")     # Table 1: V/F = 39.6 L

    # ==========================================================================
    # Absorption HALF-LIVES for IM administration (paper's reference route)
    # Converted from t1/2 to rate constants inside model(): ka = log(2) / t1/2.
    # ==========================================================================
    lthalf_abs1_im <- log(10.2);    label("Slow absorption half-life for IM route, t1/2,abs-1,IM (day)")            # Table 1: t1/2,abs-1,IM = 10.2 days
    lthalf_abs2_im <- log(0.97);    label("Fast absorption half-life for IM route, t1/2,abs-2,IM (day)")            # Table 1: t1/2,abs-2,IM = 0.97 days
    lthalf_abs3_im <- log(0.368);   label("Immediate absorption half-life for IM route, t1/2,abs-3,IM (day)")       # Table 1: t1/2,abs-3,IM = 0.368 days
    lthalf_tr_im   <- log(0.978);   label("Transit half-life for IM route, t1/2,tr,IM (day)")                       # Table 1: t1/2,tr,IM = 0.978 days

    # ==========================================================================
    # Dose-split ratios for IM (paper's RAT-transit, RAT-slowfast)
    # RAT-transit  = (dose via transit path) / (dose via immediate path)
    # RAT-slowfast = (dose via slow path)    / (dose via fast path), within transit
    # Fractions derived inside model():
    #   f_immediate = 1 / (RAT_transit + 1)
    #   f_transit   = RAT_transit / (RAT_transit + 1)
    #   f_slow      = f_transit * RAT_slowfast / (RAT_slowfast + 1)
    #   f_fast      = f_transit / (RAT_slowfast + 1)
    # These fractions reproduce Table 2 %dose splits for IM
    # (71.3 / 24.9 / 3.8 %) within rounding.
    # ==========================================================================
    lrat_transit_im  <- log(25.8);  label("Log dose-split ratio RAT-transit for IM route (transit vs immediate paths, unitless)")  # Table 1: RAT-IM = 25.8
    lrat_slowfast_im <- log(2.96);  label("Log dose-split ratio RAT-slowfast for IM route (slow vs fast within transit path, unitless)")  # Table 1: RAT2-IM = 2.96

    # ==========================================================================
    # SC-relative log multiplicative factors on the IM structural parameters
    # (Kado 2020 Table 1 "Relative structural model parameters for SC
    # administration" section).
    # log(SC value) = log(IM value) + <log-relative factor>.
    # Encoded as log parameters so that the SC IIVs (which the paper
    # estimates on the log-scale absolute SC values) attach directly to
    # the eta<X> convention -- i.e. etalthalf_abs1_screl is the IIV on
    # log(SC absolute t1/2,abs-1), and is equivalent to IIV on the SC
    # absolute value because log(IM * relative) = log(IM) + log(relative).
    # ==========================================================================
    lthalf_abs1_screl  <- log(2.12);    label("Log SC-relative multiplier on t1/2,abs-1 (log(SC / IM), unitless)")   # Table 1: t1/2,abs-1,SC relative = 2.12
    lthalf_abs2_screl  <- log(0.997);   label("Log SC-relative multiplier on t1/2,abs-2 (log(SC / IM), unitless)")   # Table 1: t1/2,abs-2,SC relative = 0.997
    lthalf_abs3_screl  <- log(0.719);   label("Log SC-relative multiplier on t1/2,abs-3 (log(SC / IM), unitless)")   # Table 1: t1/2,abs-3,SC relative = 0.719
    lthalf_tr_screl    <- log(0.937);   label("Log SC-relative multiplier on t1/2,tr (log(SC / IM), unitless)")      # Table 1: t1/2,tr,SC relative = 0.937
    lrat_transit_screl <- log(1.965);   label("Log SC-relative multiplier on RAT-transit (log(SC / IM), unitless)")  # Table 1: RAT-SC relative = 1.965
    lrat_slowfast_screl <- log(2.53);   label("Log SC-relative multiplier on RAT-slowfast (log(SC / IM), unitless)") # Table 1: RAT2-SC relative = 2.53
    lfsc               <- log(0.957);   label("Log SC relative bioavailability F_SC (log(SC / IM), unitless)")        # Table 1: F_SC_relative = 0.957

    # ==========================================================================
    # Inter-individual variability
    # ==========================================================================
    # Paper Table 1 IIVs are reported as "100% * sqrt(variability estimate)"
    # -- i.e. the reported percentages are CV%. For log-normal individual
    # parameters the corresponding omega^2 = log(1 + CV^2).
    #
    # Kado 2020 estimated separate IIVs for IM vs SC administration on the
    # absorption parameters (Methods 'PK modelling'), plus a common IIV on
    # V/F. Because each subject received both routes in the crossover, we
    # carry both an IM and an SC eta per subject; the ROUTE_SC indicator
    # switches which eta contributes inside model() (see algebraic form
    # eta_im * (1 - ROUTE_SC) + eta_sc * ROUTE_SC).
    #
    # The paper's Table 1 also reports a "r(t1/2,abs-1, t1/2,abs-1)"
    # correlation of 0.661 for IM and 0.840 for SC. The apparent
    # self-correlation is a printing error (a parameter cannot correlate
    # with itself); the pair the correlation actually spans is not
    # unambiguously stated in the paper text and is not captured here.
    # Simulations that require the between-parameter correlation should
    # override the omega block; see vignette Errata.
    etalvc                  ~ 0.01432   # Table 1: IIV on V/F = 12 CV%; omega^2 = log(1 + 0.12^2) = 0.01432

    etalthalf_abs1_im       ~ 0.12876   # Table 1: IIV on t1/2,abs-1,IM = 37 CV%; omega^2 = log(1 + 0.37^2) = 0.12876
    etalthalf_abs2_im       ~ 0.05188   # Table 1: IIV on t1/2,abs-2,IM = 23 CV%; omega^2 = log(1 + 0.23^2) = 0.05188
    etalrat_transit_im      ~ 0.05610   # Table 1: IIV on RAT-IM      = 24 CV%; omega^2 = log(1 + 0.24^2) = 0.05610
    etalrat_slowfast_im     ~ 0.10917   # Table 1: IIV on RAT2-IM     = 34 CV%; omega^2 = log(1 + 0.34^2) = 0.10917

    etalthalf_abs1_screl    ~ 0.23044   # Table 1: IIV on t1/2,abs-1,SC = 51 CV%; omega^2 = log(1 + 0.51^2) = 0.23044 (attached to the log SC-relative multiplier; log-space variance is identical to the variance on log(SC absolute))
    etalthalf_abs2_screl    ~ 0.23044   # Table 1: IIV on t1/2,abs-2,SC = 51 CV%; omega^2 = log(1 + 0.51^2) = 0.23044
    etalrat_transit_screl   ~ 0.25489   # Table 1: IIV on RAT-SC       = 54 CV%; omega^2 = log(1 + 0.54^2) = 0.25489
    etalrat_slowfast_screl  ~ 0.25489   # Table 1: IIV on RAT2-SC      = 54 CV%; omega^2 = log(1 + 0.54^2) = 0.25489

    # ==========================================================================
    # Residual variability
    # ==========================================================================
    # Table 1 residual variability RV = 22%. Proportional error model on
    # linear concentration scale (Kado 2020 Table 1 caption: "IIV and RV
    # are presented as 100% * sqrt(variability estimate)").
    propSd <- 0.22;                  label("Proportional residual error (fraction)")   # Table 1: RV = 22 CV%
  })

  model({
    # Convert paper's t1/2 half-lives and dose-split ratios into individual
    # subject values, switching between IM and SC parameters via ROUTE_SC.
    #
    # For each parameter the encoding
    #     x_i = exp(l<x>_im
    #             + l<x>_screl        * ROUTE_SC
    #             + eta<l<x>_im>      * (1 - ROUTE_SC)
    #             + eta<l<x>_screl>   * ROUTE_SC)
    # collapses to
    #     ROUTE_SC = 0 -> x_i = exp(lX_im + etalX_im)                     (IM)
    #     ROUTE_SC = 1 -> x_i = exp(lX_im + lX_screl + etalX_screl)       (SC)
    # The SC log-space variance attaches to the SC-relative multiplier
    # (etalX_screl), which is equivalent to attaching it to the SC absolute
    # value because log(IM * relative) = log(IM) + log(relative).
    thalf_abs1 <- exp(lthalf_abs1_im +
                      lthalf_abs1_screl    * ROUTE_SC +
                      etalthalf_abs1_im    * (1 - ROUTE_SC) +
                      etalthalf_abs1_screl * ROUTE_SC)
    thalf_abs2 <- exp(lthalf_abs2_im +
                      lthalf_abs2_screl    * ROUTE_SC +
                      etalthalf_abs2_im    * (1 - ROUTE_SC) +
                      etalthalf_abs2_screl * ROUTE_SC)
    # No IIV on t1/2,abs-3 or t1/2,tr in Table 1 -- immediate absorption
    # and the transit half-life are fixed-effect kinetic components.
    thalf_abs3 <- exp(lthalf_abs3_im + lthalf_abs3_screl * ROUTE_SC)
    thalf_tr   <- exp(lthalf_tr_im   + lthalf_tr_screl   * ROUTE_SC)

    rat_transit  <- exp(lrat_transit_im +
                        lrat_transit_screl    * ROUTE_SC +
                        etalrat_transit_im    * (1 - ROUTE_SC) +
                        etalrat_transit_screl * ROUTE_SC)
    rat_slowfast <- exp(lrat_slowfast_im +
                        lrat_slowfast_screl    * ROUTE_SC +
                        etalrat_slowfast_im    * (1 - ROUTE_SC) +
                        etalrat_slowfast_screl * ROUTE_SC)

    # ---- Absorption + transit rate constants ---------------------------------
    ka_slow <- log(2) / thalf_abs1
    ka_fast <- log(2) / thalf_abs2
    ka_imm  <- log(2) / thalf_abs3
    ktr     <- log(2) / thalf_tr

    # ---- Dose-fraction split per pathway -------------------------------------
    # Kado 2020 Methods: RAT-transit governs the transit-vs-immediate split;
    # RAT-slowfast governs the slow-vs-fast split within the transit path.
    f_immediate <- 1 / (rat_transit + 1)
    f_transit   <- rat_transit / (rat_transit + 1)
    f_slow      <- f_transit * rat_slowfast / (rat_slowfast + 1)
    f_fast      <- f_transit /                (rat_slowfast + 1)

    # ---- SC relative bioavailability -----------------------------------------
    # F_SC applies uniformly to all three absorption pathways; IM has F = 1
    # (fixed reference). Encoded as exp(lfsc * ROUTE_SC) which is 1 when
    # ROUTE_SC = 0 (IM) and F_SC = exp(lfsc) = 0.957 when ROUTE_SC = 1.
    f_route <- exp(lfsc * ROUTE_SC)

    # ---- Disposition ---------------------------------------------------------
    kel <- exp(lkel)
    vc  <- exp(lvc + etalvc)

    # ---- ODE system ----------------------------------------------------------
    # Three parallel absorption pathways feeding one central compartment:
    #   Slow path:      depot1 -> transit1 -> central   (via ktr, ka_slow)
    #   Fast path:      depot2 -> transit2 -> central   (via ktr, ka_fast)
    #   Immediate path: depot3 -> central                (via ka_imm; bypasses transit)
    # The dose is deposited into all three depots simultaneously, with
    # bioavailability fractions f_slow / f_fast / f_immediate (derived
    # from the paper's RAT-transit and RAT-slowfast) partitioning the
    # total dose across pathways.
    d/dt(depot1)   <- -ktr * depot1
    d/dt(transit1) <-  ktr * depot1  - ka_slow * transit1
    d/dt(depot2)   <- -ktr * depot2
    d/dt(transit2) <-  ktr * depot2  - ka_fast * transit2
    d/dt(depot3)   <- -ka_imm * depot3
    d/dt(central)  <- (ka_slow * transit1 +
                       ka_fast * transit2 +
                       ka_imm  * depot3 -
                       kel     * central)

    # ---- Bioavailability per depot -------------------------------------------
    # Users dose all three depots (depot1 = slow, depot2 = fast,
    # depot3 = immediate) with the same total dose amount per injection;
    # the f() values partition the dose so that the total absorbed mass
    # equals dose * f_route (= dose for IM; = 0.957 * dose for SC).
    f(depot1) <- f_slow      * f_route
    f(depot2) <- f_fast      * f_route
    f(depot3) <- f_immediate * f_route

    # ---- Observation ---------------------------------------------------------
    # central in mg, vc in L -> mg/L; multiply by 1000 to report ng/mL to
    # match Kado 2020 Table 3 concentration units and the 20 / 10 ng/mL
    # PK/PD targets discussed throughout the paper.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
