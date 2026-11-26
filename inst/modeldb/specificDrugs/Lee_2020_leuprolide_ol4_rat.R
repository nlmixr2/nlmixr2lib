Lee_2020_leuprolide_ol4_rat <- function() {
  description <- paste(
    "Preclinical (rat).",
    "Level A IVIVC-PK model for the OL4 (RG502:RG502H 25:75) leuprolide acetate",
    "electrosprayed PLGA microsphere formulation",
    "after a single subcutaneous dose in normal male Wistar rats.",
    "The dose is split at time zero across a non-capsuled section that is",
    "already in the release compartment (depot), a diffusive-release section",
    "emptying into it through a 2-stage transit chain at ktr-d, and",
    "an erosive-release section emptying into it through a 2-stage",
    "transit chain at ktr-e. Drug in the release compartment is absorbed at",
    "the solution-group ka into a one-compartment disposition model with",
    "apparent CL and Vd carried from the authors' earlier solution study.",
    "Release parameters from Table 2 of Lee 2020."
  )
  reference <- paste(
    "Lee DS, Kang DW, Choi GW, Choi HG, Cho HY.",
    "Development of Level A In Vitro-Vivo Correlation for Electrosprayed",
    "Microspheres Containing Leuprolide: Physicochemical, Pharmacokinetic,",
    "and Pharmacodynamic Evaluation.",
    "Pharmaceutics. 2020;12(1):36. doi:10.3390/pharmaceutics12010036"
  )
  vignette <- "Lee_2020_leuprolide"

  units <- list(time = "h", dosing = "ng", concentration = "ng/mL")

  # Equation 1 starts each section with its share of the dose (NS = dose*E1,
  # DS1 = dose*E2, ES1 = dose*E3), so the event table carries three
  # simultaneous dose records of the full amount and the f() fractions below
  # split them.
  dosing <- c("depot", "transit_fast1", "transit_slow1")

  covariateData <- list()

  compartmentData <- list(
    depot = list(
      analyte = "leuprolide",
      units = "ng",
      specimen = "administration site",
      verified = TRUE
    ),
    transit_fast1 = list(
      analyte = "leuprolide",
      units = "ng",
      specimen = "administration site",
      verified = TRUE
    ),
    transit_fast2 = list(
      analyte = "leuprolide",
      units = "ng",
      specimen = "administration site",
      verified = TRUE
    ),
    transit_slow1 = list(
      analyte = "leuprolide",
      units = "ng",
      specimen = "administration site",
      verified = TRUE
    ),
    transit_slow2 = list(
      analyte = "leuprolide",
      units = "ng",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "leuprolide",
      units = "ng",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "rat (Wistar)",
    n_subjects = 5L,
    n_studies = 1L,
    sex_female_pct = 0,
    disease_state = "healthy (normal male rats)",
    dose_range = "single subcutaneous dose equivalent to 0.1 mg/kg leuprolide acetate (25 ug per rat implied by CL/F x AUCinf)",
    regions = "Korea",
    notes = paste(
      "OL4 was never dosed in vivo; its release parameters come only",
      "from the in vitro release test (Figure 4, n = 3 tubes) time-scaled to",
      "in vivo by the Levy-plot correlation. The rat context is that of the",
      "Lucrin depot group whose fit the correlation was built on (n_subjects).",
      "Release parameters estimated by fitting the IVIVC-predicted in vivo",
      "release profile of this formulation (its observed in vitro release",
      "time-scaled through the Levy-plot Hill correlation of Equation 3) with",
      "the absorption step excluded; no plasma data informed them.",
      "Body weight is not reported; Table 3's CL/F x AUCinf products all",
      "equal 25,000 ng, the per-rat dose the vignette uses."
    )
  )

  ini({
    # ---------------- Disposition (carried from Lee 2018) ------------------
    # Methods 2.9.1: 'ka, CL, and Vd were all estimated by one-compartmental
    # analysis from the leuprolide solution-administered group of our previous
    # study [10]', and are applied unchanged to every formulation (Methods
    # 2.9.3), so they are not estimated in this paper. Apparent values: only
    # subcutaneous data informed them.
    lka <- fixed(log(16.67))
    label("Absorption rate constant ka from the release compartment (1/h)")
    # Methods 2.9.1: '16.67 h-1 of ka'

    lcl <- fixed(log(514.46))
    label("Apparent clearance CL/F (mL/h)")
    # Methods 2.9.1: '514.46 mL/h of CL'

    lvc <- fixed(log(487.40))
    label("Apparent volume of distribution Vd/F (mL)")
    # Methods 2.9.1: '487.40 mL of Vd'

    # ---------------- Release compartment ----------------------------------
    # Equation 1 constrains E1 + E2 + E3 = 1, so the three section fractions
    # are encoded by stick-breaking: logitfrel is E1 itself and logitfrel2 is
    # E2 conditional on the dose not being in the non-capsuled section. The
    # back-transform in model() reproduces the printed E1 / E2 / E3 exactly.
    logitfrel <- log(0.210 / (1 - 0.210))
    label("Non-capsuled (initial burst) fraction E1, logit scale (unitless)")
    # Table 2, column 'OL4' row 'E1' = 0.210

    logitfrel2 <- log((0.780 / (1 - 0.210)) / (1 - 0.780 / (1 - 0.210)))
    label("Diffusive-release fraction E2 conditional on not being non-capsuled, logit scale (unitless)")
    # Table 2, column 'OL4' row 'E2' = 0.780; conditional value 0.780 / (1 - 0.210) = 0.98734, so
    # the erosive remainder is (1 - 0.210) * (1 - 0.98734) = 0.010 = row 'E3' = 0.010

    lktr_fast <- log(0.047)
    label("Diffusive-release transit rate constant ktr-d (1/h)")
    # Table 2, column 'OL4' row 'ktr-d (h-1)' = 0.047

    lktr_slow <- log(0.0000002)
    label("Erosive-release transit rate constant ktr-e (1/h)")
    # Table 2, column 'OL4' row 'ktr-e (h-1)' = 0.0000002

    # ---------------- Residual unexplained variability ---------------------
    # The model was fitted to mean profiles in Berkeley Madonna and no residual
    # error or between-animal variability is reported.
    propSd <- fixed(0)
    label("Proportional residual SD on plasma leuprolide (fraction; 0, not reported in the source)")
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    ktr_fast <- exp(lktr_fast)
    ktr_slow <- exp(lktr_slow)
    kel <- cl / vc

    # Section fractions (Equation 1: E1 + E2 + E3 = 1). frel = E1 (NS),
    # frel2 = E2 (DS) and frel3 = E3 (ES).
    frel <- expit(logitfrel)
    frel2 <- (1 - frel) * expit(logitfrel2)
    frel3 <- (1 - frel) * (1 - expit(logitfrel2))

    # Equation 1. NS is the release compartment (depot); DS1..DSn and
    # ES1..ESn are the diffusive and erosive transit chains (transit_fast<i>
    # and transit_slow<i>). Equation 1 as printed writes ktr-d inside the ES
    # chain but ktr-e for its outflow into NS, which would not conserve mass;
    # Methods 2.9.1 and Figure 2 state that ES releases at ktr-e, so the whole
    # ES chain runs at ktr-e here.
    # Table 2, column 'OL4' row 'Number of transit (DS)' = 2.020, implemented as 2 stages
    # (nearest integer); see the rounding note below.
    # Table 2, column 'OL4' row 'Number of transit (ES)' = 2.068, implemented as 2 stages
    # (nearest integer); see the rounding note below.
    # Table 2 prints non-integer transit-stage counts for OL1-OL5 although
    # Equation 1 builds each chain from discrete compartments DS1..DSn and
    # ES1..ESn (Berkeley Madonna arrays). Rounding each printed count to the
    # nearest integer reproduces every curve of Figure 6A (for example OL3,
    # 2.759 -> 3 stages, gives the 0.024 ng/mL trough near 8-10 h and the
    # 0.095 ng/mL plateau near 56 h that the figure shows), whereas truncating
    # it or treating the count as a continuous gamma-distributed delay does
    # not. See the vignette, section 'Transit-stage counts'.
    d/dt(depot) <- ktr_fast * transit_fast2 + ktr_slow * transit_slow2 - ka * depot
    d/dt(transit_fast1) <- -ktr_fast * transit_fast1
    d/dt(transit_fast2) <- ktr_fast * transit_fast1 - ktr_fast * transit_fast2
    d/dt(transit_slow1) <- -ktr_slow * transit_slow1
    d/dt(transit_slow2) <- ktr_slow * transit_slow1 - ktr_slow * transit_slow2
    d/dt(central) <- ka * depot - kel * central

    # Initial values of NS, DS1 and ES1 are dose*E1, dose*E2 and dose*E3
    # (Methods 2.9.1); each receives the full dose amount on its own record.
    f(depot) <- frel
    f(transit_fast1) <- frel2
    f(transit_slow1) <- frel3

    Cc <- central / vc

    Cc ~ prop(propSd)
  })
}
