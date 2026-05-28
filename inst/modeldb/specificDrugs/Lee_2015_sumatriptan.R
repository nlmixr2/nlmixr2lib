Lee_2015_sumatriptan <- function() {
  description <- "One-compartment population PK model for oral sumatriptan in healthy Korean male volunteers (Lee 2015): two parallel absorption routes (first-order absorption with lag time, and a transit-compartment chain with the Savic 2007 analytical input form) into a single central compartment with linear elimination. Captures the multiple-peaks absorption phenomenon reported in oral sumatriptan."
  reference <- paste(
    "Lee J, Lim M, Seong SJ, Park S-M, Gwon M-R, Han S, Lee SM, Kim W, Yoon Y-R, Yoo H-D.",
    "Population pharmacokinetic analysis of the multiple peaks phenomenon in sumatriptan.",
    "Transl Clin Pharmacol. 2015;23(2):66-74.",
    "doi:10.12793/tcp.2015.23.2.66.",
    sep = " "
  )
  vignette <- "Lee_2015_sumatriptan"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 26L,
    n_studies      = 1L,
    age_range      = "22-28 years",
    age_median     = "23.9 years (mean)",
    weight_range   = "51-84 kg",
    weight_median  = "66.7 kg (mean)",
    sex_female_pct = 0,
    race_ethnicity = "Korean (single-center Korean cohort)",
    disease_state  = "Healthy adult male volunteers",
    dose_range     = "Single 50 mg oral dose of sumatriptan succinate (50 mg expressed as sumatriptan free base)",
    regions        = "Republic of Korea (Kyungpook National University Hospital Clinical Trial Center, Daegu)",
    notes          = "Retrospective re-analysis of the reference-formulation arm of a single-center, randomized, open-label, two-period, single-dose crossover bioequivalence study (364 plasma concentrations). Demographics per Table 1 of Lee 2015. Only data from the reference formulation entered the popPK analysis."
  )

  # Implementation notes (see vignette 'Assumptions and deviations' for the
  # full justification):
  # * M5 base model: Lee 2015 develops a base structural model (M5, one-cmt
  #   + parallel first-order + Savic transit absorption) and selects M6 (M5
  #   plus a creatinine-clearance covariate on the transit-fraction f) as
  #   the final published model. The supplement's NONMEM code only covers
  #   M5; the M6 covariate functional form and reference CrCL are not
  #   stated anywhere on disk. Per the skill's "Covariate encoding
  #   ambiguous" stop-and-ask trigger and the operator's directive, this
  #   extraction encodes M5 (no covariate) and documents the M6 covariate
  #   in the vignette's Assumptions and deviations / Errata section.
  # * Savic analytical input rate: the NONMEM code uses Stirling's
  #   approximation log(N!) = log(sqrt(2*pi)) + (N+0.5)*log(N) - N to
  #   evaluate the gamma kernel for non-integer N. rxode2's built-in
  #   transit(n, mtt, bio) computes the same kernel using lgamma(n+1)
  #   directly; the numerical difference at N = 11 is ~0.03% in
  #   log(N!) and is well below simulation precision. The packaged
  #   model therefore uses transit() rather than transcribing the
  #   Stirling approximation.
  # * Two depots, two dose records: depot1 receives (1 - fr) * AMT via
  #   standard bolus dosing with first-order absorption (rate ka1) and
  #   lag time alag1; depot2 has its bolus suppressed (f(depot2) = 0)
  #   and is fed by the analytical Savic input rate, then absorbs into
  #   central at rate ka2. This matches the supplement's $MODEL block
  #   (COMP=(DEPOT1), COMP=(DEPOT2), COMP=(CENTRAL, DEFOBS)) and the $PK
  #   assignments F1 = 1-FR and F2 = 0.
  # * Dosing convention required by rxode2's transit(): rxode2's
  #   transit(n, mtt, bio) is compartment-aware: it fires only when the
  #   dose event targets the same compartment whose d/dt() contains the
  #   transit() call (otherwise podo() / tad() return zero and the gamma
  #   kernel is silent). Lee 2015 routes one user-facing dose through
  #   two parallel arms, so any user data set MUST include two dose
  #   events per administration: one to `depot` (which carries the
  #   first-order arm's bolus, scaled by f(depot) = 1 - fr) and one to
  #   `depot2` (which is suppressed via f(depot2) = 0 but whose
  #   presence makes transit() read podo() = AMT for the gamma chain).
  #   Both records carry the same AMT; the model's bioavailability
  #   settings perform the (1 - fr) / fr split internally. See the
  #   validation vignette for the canonical event-table construction.

  ini({
    # Structural parameters - all from Lee 2015 Table 3 ("Final parameter
    # estimates and bootstrap results"). The supplement's NONMEM code
    # ($PK block, page 74) parameterises each structural quantity as
    # THETA(i) * EXP(ETA(i)), so each typical value is log-transformed
    # here (per nlmixr2lib parameter-naming convention).
    lcl    <- log(418);  label("Apparent clearance CL/F (L/h)")                                # Table 3 final estimate CL/F = 418 L/h (RSE 4%)
    lvc    <- log(56.9); label("Apparent central volume of distribution V/F (L)")              # Table 3 final estimate V/F = 56.9 L (RSE 35.4%)
    lka1   <- log(0.62); label("First-order absorption rate constant from depot1 (1/h)")       # Table 3 final estimate ka1 = 0.62 1/h (RSE 9.13%)
    lka2   <- log(0.29); label("Absorption rate constant from final transit cmt -> central (1/h)") # Table 3 final estimate ka2 = 0.29 1/h (RSE 6.89%)
    lmtt   <- log(1.94); label("Mean transit time of the transit-compartment chain (h)")       # Table 3 final estimate MTT = 1.94 h (RSE 9.89%)
    ln     <- log(11);   label("Number of transit compartments before depot2 (continuous)")     # Table 3 final estimate n = 11 (RSE 23.2%); estimated as a continuous quantity via the Stirling/gamma analytical input form
    lalag1 <- log(0.24); label("Lag time for first-order absorption from depot1 (h)")          # Table 3 final estimate ALAG1 = 0.24 h (RSE 1.32%)
    lfr    <- log(0.56); label("Fraction of dose routed through the transit-compartment arm (unitless)") # Table 3 final estimate f = 0.56 (RSE 6.18%); the remaining (1 - fr) goes via depot1 first-order absorption

    # Inter-individual variability. Lee 2015 reports BSV as CV% on a
    # log-normal scale (NONMEM exponential ETA pattern); the internal
    # variance is omega^2 = log(1 + CV^2). No BSV reported on ka1, n, or
    # ALAG1 in Table 3 ("-" in the BSV column), so those etas are omitted.
    # No $OMEGA BLOCK is declared in the supplement's NONMEM code (page
    # 74), so the etas are independent (diagonal omega).
    etalcl ~ 0.0337   # log(1 + 0.185^2) = log(1.0342) -- Table 3 BSV(CL/F) = 18.5% CV (RSE 25.9%, shrinkage 6.41%)
    etalvc ~ 0.4072   # log(1 + 0.709^2) = log(1.5027) -- Table 3 BSV(V/F) = 70.9% CV (RSE 28.9%, shrinkage 2.16%)
    etalka2 ~ 0.0587  # log(1 + 0.246^2) = log(1.0605) -- Table 3 BSV(ka2) = 24.6% CV (RSE 39.3%, shrinkage 14.2%)
    etalmtt ~ 0.1193  # log(1 + 0.356^2) = log(1.1267) -- Table 3 BSV(MTT) = 35.6% CV (RSE 29.6%, shrinkage 7.1%)
    etalfr  ~ 0.0205  # log(1 + 0.144^2) = log(1.0207) -- Table 3 BSV(f)  = 14.4% CV (RSE 40.7%, shrinkage 18.6%)

    # Residual error - combined additive plus proportional in linear
    # concentration space. The supplement's $ERROR block writes
    #   W = SQRT(THETA(9)^2 + THETA(10)^2 * IPRED^2)
    #   Y = IPRED + W * ERR(1)
    # so THETA(9) is the additive SD and THETA(10) is the proportional
    # SD (both on the IPRED scale, which is ng/mL after the NONMEM
    # S3 = V3/1000 rescaling).
    propSd <- 0.21;  label("Proportional residual error SD (fraction)") # Table 3 final estimate proportional error = 0.21 (RSE 7.3%)
    addSd  <- 0.30;  label("Additive residual error SD (ng/mL)")        # Table 3 final estimate additive error = 0.30 ng/mL (RSE 15.2%)
  })

  model({
    # Individual structural parameters (back-transformed to the linear
    # scale, with eta where applicable).
    cl    <- exp(lcl   + etalcl)
    vc    <- exp(lvc   + etalvc)
    ka1   <- exp(lka1)
    ka2   <- exp(lka2  + etalka2)
    mtt   <- exp(lmtt  + etalmtt)
    n     <- exp(ln)
    alag1 <- exp(lalag1)
    fr    <- exp(lfr   + etalfr)

    kel <- cl / vc

    # ODE system. Lee 2015 Figure 2 schematic + supplement $DES block.
    # depot1: first-order absorption arm, dose record entry point. After
    #   the lag time alag1, drug absorbs at rate ka1 into central.
    # depot2: transit-compartment-chain arm. The analytical Savic 2007
    #   input rate transit(n, mtt, fr) carries fr * AMT through an
    #   n-compartment chain with rate ktr = (n+1)/mtt; depot2 then
    #   absorbs into central at rate ka2. f(depot2) = 0 suppresses the
    #   bolus deposition for any data record that targets depot2
    #   directly so that all transit-arm input arrives via the
    #   analytical kernel.
    # central: linear elimination at rate kel = cl/vc.
    d/dt(depot)   <- -ka1 * depot
    d/dt(depot2)  <-  transit(n, mtt, fr) - ka2 * depot2
    d/dt(central) <-  ka1 * depot + ka2 * depot2 - kel * central

    # Dose-routing bookkeeping. The user dose record targets depot; with
    # f(depot) = 1 - fr only the first-order fraction (1 - fr) of AMT
    # enters depot as a bolus. f(depot2) = 0 suppresses any bolus
    # deposition into depot2 so the transit-arm input arrives solely
    # through transit(). alag(depot) applies the lag only to the
    # first-order arm (Lee 2015 supplement: ALAG1 on COMP 1 only).
    f(depot)    <- 1 - fr
    f(depot2)   <- 0
    alag(depot) <- alag1

    # Plasma concentration: central in mg, vc in L gives mg/L = ug/mL;
    # multiply by 1000 to match Lee 2015's bioanalytical units of
    # ng/mL (NONMEM supplement's S3 = V3/1000 rescaling).
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
