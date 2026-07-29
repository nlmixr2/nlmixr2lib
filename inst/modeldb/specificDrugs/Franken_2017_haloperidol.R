Franken_2017_haloperidol <- function() {
  description <- paste(
    "One-compartment population PK model for haloperidol in 28 terminally",
    "ill adult palliative-care patients (Franken 2017). Two parallel",
    "first-order absorption routes (oral and subcutaneous) with",
    "route-specific absorption rate constants fixed from literature",
    "(Ka oral = 0.236 1/h, Ka SC = 20 1/h derived from intramuscular Tmax",
    "= 20 min). Oral bioavailability F = 0.861 is estimated; SC F is",
    "assumed to be 1. IIV is included on F, CL, and Vd; the IIV on F and",
    "CL was 99% correlated and is encoded with correlation fixed to unity",
    "(BLOCK pattern). Residual variability is additive on log-transformed",
    "concentrations (LTBS). Covariate analysis (body weight, age, sex,",
    "primary diagnosis, plasma creatinine, urea, bilirubin, GGT, ALP,",
    "ALT, AST, CRP, albumin, concomitant CYP2D6 / CYP3A inducers and",
    "inhibitors, time-to-death) did not retain any covariate in the",
    "final model."
  )
  reference <- paste(
    "Franken LG, Mathot RAA, Masman AD, Baar FPM, Tibboel D,",
    "van Gelder T, Koch BCP, de Winter BCM.",
    "Population pharmacokinetics of haloperidol in terminally ill",
    "adult patients. Eur J Clin Pharmacol. 2017;73(10):1271-1277.",
    "doi:10.1007/s00228-017-2283-6.",
    sep = " "
  )
  vignette <- "Franken_2017_haloperidol"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ug/L"
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 28,
    n_studies      = 1,
    age_range      = "43-93 years",
    age_median     = "69.5 years",
    weight_range   = "35-108 kg",
    weight_median  = "67 kg",
    sex_female_pct = 46.4,
    race_ethnicity = c(Caucasian = 92.9, AfroCaribbean = 7.1),
    disease_state  = paste(
      "Terminally ill adult palliative-care patients with advanced",
      "malignancy (neoplasm in 100%; epithelial-tissue primary in 89.3%)",
      "receiving haloperidol for the treatment of delirium. Survival",
      "prognosis at enrolment was 2 days to 3 months; patients followed",
      "until time of death."
    ),
    dose_range     = paste(
      "Oral 0.5-2 mg/day (tablets or liquid) or subcutaneous bolus",
      "0.5-5 mg/day, dosed per Dutch national palliative guidelines."
    ),
    regions        = "The Netherlands (Laurens Cadenza palliative care centre, Rotterdam)",
    notes          = paste(
      "Demographics from Franken 2017 Table 1. NONMEM 7.2 + PsN 4.4.8",
      "with FOCE-I; ADVAN5 subroutine used to handle the two parallel",
      "absorption routes. 87 sparse plasma samples (median 3 per",
      "subject, range 1-9) by venous puncture or indwelling catheter.",
      "Assay: LC-MS/MS with LLOQ 0.5 ug/L (range 0.5-125 ug/L). 14.6%",
      "BQL after excluding samples > 200 h after last dose; M1",
      "(discard) handling used in the final model. Body weight was",
      "missing for ~35% of subjects and was imputed at the population",
      "median (67 kg) during covariate testing, but allometric scaling",
      "was not retained in the final model. Final model validated by",
      "500-run bootstrap and NPDE analysis."
    )
  )

  ini({
    # ================================================================
    # Absorption rate constants - both routes fixed from literature.
    # Methods 'Structural model' / Table 2 footnote a: Ka could not be
    # estimated because of limited absorption-phase data; oral Ka =
    # 0.236 1/h from a prior haloperidol PK study, and SC Ka = 20 1/h
    # was derived from the intramuscular Tmax of 20 min (= ln(2)/k
    # back-calculation) used as a reference because no SC literature
    # exists for the intravenous formulation used here.
    # ================================================================
    lka_oral <- fixed(log(0.236))
    label("Ka oral route (1/h, FIXED literature)")              # Table 2, footnote a: Ka_oral = 0.236 1/h fixed
    lka_sc   <- fixed(log(20))
    label("Ka subcutaneous route (1/h, FIXED literature)")      # Table 2, footnote a: Ka_SC = 20 1/h fixed from IM Tmax = 20 min

    # ================================================================
    # Oral bioavailability - estimated. SC F is structurally fixed at
    # 1 ('Population pharmacokinetic method': the bioavailability of
    # subcutaneous haloperidol was assumed to be 100%).
    # ================================================================
    lfdepot <- log(0.861)
    label("Oral bioavailability F (log scale)")                 # Table 2: F = 0.861 (RSE 18%)

    # ================================================================
    # Structural disposition - Franken 2017 Table 2 Final model.
    # ================================================================
    lcl <- log(29.3)
    label("Clearance CL (L/h)")                                 # Table 2: CL = 29.3 L/h (RSE 11%)
    lvc <- log(1260)
    label("Volume of distribution Vd (L)")                      # Table 2: Vd = 1260 L (RSE 19%)

    # ================================================================
    # Inter-individual variability.
    # IIV reported as CV% in Table 2; converted to omega^2 via
    # omega^2 = log(1 + CV^2).
    # 'Structural model': the IIV on CL and F showed a 99% correlation
    # and were fixed to unity with the addition of an extra theta.
    # Encoded here as a 2x2 BLOCK with the off-diagonal equal to
    # sqrt(var_F * var_CL), which forces rho = 1; this is the
    # mathematically equivalent representation of the published
    # shared-eta-with-scaling-theta form.
    # ================================================================
    etalfdepot + etalcl ~ c(
      log(1 + 0.55^2),
      sqrt(log(1 + 0.55^2) * log(1 + 0.43^2)),
      log(1 + 0.43^2)
    )
    # Table 2 IIV: F = 55% CV (RSE 43%, shrinkage 37%), CL = 43% CV
    # (RSE 34%, shrinkage 29%); rho fixed to 1.

    etalvc ~ log(1 + 0.70^2)
    # Table 2 IIV: Vd = 70% CV (RSE 21%, shrinkage 31%).

    # ================================================================
    # Residual variability.
    # 'Structural model': additive residual error on log-transformed
    # concentrations (LTBS). Table 2 reports 0.258 as the NONMEM $SIGMA
    # variance on the log scale; the SD supplied to lnorm() is sqrt.
    # ================================================================
    expSd <- sqrt(0.258)
    label("Residual SD on log scale (LTBS)")                    # Table 2: residual = 0.258 (RSE 22%, shrinkage 19%); variance on log scale -> SD = sqrt(0.258)
  })

  model({
    # ----------------------------------------------------------------
    # Route-specific absorption rate constants (both fixed).
    # ----------------------------------------------------------------
    ka_oral <- exp(lka_oral)
    ka_sc   <- exp(lka_sc)

    # ----------------------------------------------------------------
    # Oral bioavailability with IIV (shared eta with CL; rho = 1 in
    # the OMEGA block above). SC bioavailability is the rxode2 default
    # of 1 because no f(depot) assignment is made for the SC depot.
    # ----------------------------------------------------------------
    f_oral <- exp(lfdepot + etalfdepot)

    # ----------------------------------------------------------------
    # Individual structural parameters.
    # ----------------------------------------------------------------
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    # ----------------------------------------------------------------
    # ODE system. Two parallel depots feed a single central
    # compartment with route-specific first-order absorption.
    # Convention: depot (cmt = 1) is the SC depot, depot2 (cmt = 2)
    # is the oral depot; central (cmt = 3) is the sampling
    # compartment. (Matches the dual-depot convention used in
    # Franken_2015_morphine.R: SC -> depot, oral -> depot2.)
    # ----------------------------------------------------------------
    d/dt(depot)   <- -ka_sc   * depot
    d/dt(depot2)  <- -ka_oral * depot2
    d/dt(central) <-  ka_sc   * depot + ka_oral * depot2 - kel * central

    # Oral bioavailability assignment. SC depot stays at rxode2
    # default F = 1 (per paper assumption).
    f(depot2) <- f_oral

    # ----------------------------------------------------------------
    # Observation. Internal states are mg of haloperidol; plasma
    # concentrations are reported in ug/L (paper LLOQ 0.5 ug/L).
    # 1 mg/L = 1000 ug/L.
    # ----------------------------------------------------------------
    Cc <- central / vc * 1000
    Cc ~ lnorm(expSd)
  })
}
