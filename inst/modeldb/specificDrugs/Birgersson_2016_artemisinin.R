Birgersson_2016_artemisinin <- function() {
  description <- "One-compartment population PK model for oral artemisinin in 15 healthy male Vietnamese volunteers, with a seven-compartment transit-absorption chain (number of transits fixed at 7). The published final model carries inter-occasion variability on apparent clearance and mean transit-time and inter-individual variability on relative bioavailability; for forward simulation in nlmixr2lib the IOV terms are mapped onto etas (etalcl, etalmtt) so a single-occasion simulation reproduces the population variability. The published full covariate analysis found no clinically significant effect (>20 %) of formulation, dose level (160 vs 500 mg), or concomitant piperaquine on the structural PK parameters, so no covariates are carried in the model."
  reference <- paste(
    "Birgersson S, Van Toi P, Truong NT, Dung NT, Ashton M, Hien TT,",
    "Abelo A, Tarning J (2016).",
    "Population pharmacokinetic properties of artemisinin in healthy",
    "male Vietnamese volunteers.",
    "Malaria Journal 15:90.",
    "doi:10.1186/s12936-016-1134-8.",
    sep = " "
  )
  vignette <- "Birgersson_2016_artemisinin"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 15L,
    n_studies      = 1L,
    n_observations = 786L,
    age_range      = "19-41 years",
    age_median     = "23 years",
    weight_range   = "43-80 kg",
    weight_median  = "58 kg",
    height_range   = "158-180 cm",
    height_median  = "167 cm",
    sex_female_pct = 0,
    race_ethnicity = c(Asian = 100),
    disease_state  = "healthy adult male volunteers",
    dose_range     = paste(
      "Four-way single-dose oral crossover with three-week washout:",
      "T1 = 160 mg artemisinin (micronized, 2 x 80 mg capsules);",
      "T2 = 160 mg artemisinin (reference Vietnamese low-dose formulation,",
      "2 x 80 mg capsules);",
      "T3 = 500 mg artemisinin (reference Vietnamese dose-strength formulation,",
      "2 x 250 mg capsules);",
      "T4 = 160 mg artemisinin (micronized) + 720 mg piperaquine phosphate",
      "(2 x [80 mg artemisinin + 360 mg piperaquine phosphate] tablets)."
    ),
    regions        = "Vietnam (Hospital for Tropical Diseases, Ho Chi Minh City)",
    notes          = paste(
      "Demographics from Birgersson 2016 Table 1.",
      "All 15 subjects completed all four treatment periods.",
      "Blood samples (n = 786) were drawn pre-dose and at 0.25, 0.5, 1, 1.5,",
      "2, 2.5, 3, 4, 5, 6, 7, 8, 10, and 12 h after each dose; about 5 %",
      "of samples (all within 30 min of dosing) were below the LC-MS/MS limit",
      "of quantification (1.03 ng/mL) and omitted from the fit."
    )
  )

  ini({
    # Structural population parameters (Birgersson 2016 Table 2, 'Final model'
    # column). Tested allometric body-weight scaling did not improve the fit
    # and was not retained; the typical values below are therefore the
    # whole-population estimates without size normalization.

    lcl     <- log(417)                ; label("Apparent oral clearance CL/F (L/h)")                          # Birgersson 2016 Table 2: CL/F = 417 L/h (RSE 9.32 %, 95 % CI 350-501)
    lvc     <- log(1210)               ; label("Apparent central volume of distribution V/F (L)")            # Birgersson 2016 Table 2: V/F = 1210 L (RSE 9.02 %, 95 % CI 1030-1450)
    lmtt    <- log(0.787)              ; label("Mean transit time of the absorption chain MTT (h)")          # Birgersson 2016 Table 2: MTT = 0.787 h (RSE 5.97 %, 95 % CI 0.702-0.891)
    nn_fix  <- fixed(7)                ; label("Number of transit compartments (integer, unitless)")         # Birgersson 2016 Table 2: 'Nr. trans comp' = 7 fix; carried in model() as a structural constant
    lfdepot <- fixed(log(1))           ; label("Reference relative oral bioavailability F (unitless)")        # Birgersson 2016 Table 2: F (%) = 100 (fixed); structural anchor of unity allowing IIV on F to be identified

    # Inter-individual and inter-occasion variability (Birgersson 2016 Table 2,
    # 'IIV/IOV* CV % (RSE %)' column). The paper retained IOV on CL/F and MTT
    # and IIV on F/F (the * marks IOV in the table). To allow single-occasion
    # forward simulation in nlmixr2lib, IOV is mapped onto inter-individual
    # etas; this is documented as a deviation in the vignette and is the
    # conventional simplification when carrying an IOV-bearing popPK model
    # into a library intended for clinical-trial simulation. Variances on the
    # log scale are computed as omega^2 = log(CV^2 + 1) per the paper's own
    # Table 2 footnote ('CV % = 100 * sqrt(exp(variance) - 1)').
    #   CL/F : omega^2 = log(0.171^2 + 1) = 0.02883
    #   MTT  : omega^2 = log(0.539^2 + 1) = 0.25530
    #   F    : omega^2 = log(0.343^2 + 1) = 0.11122
    etalcl     ~ 0.02883                                                                                      # Birgersson 2016 Table 2: 17.1 % IOV* on CL/F (RSE 34.3 %); encoded as IIV for single-occasion simulation
    etalmtt    ~ 0.25530                                                                                      # Birgersson 2016 Table 2: 53.9 % IOV* on MTT  (RSE 20.3 %); encoded as IIV for single-occasion simulation
    etalfdepot ~ 0.11122                                                                                      # Birgersson 2016 Table 2: 34.3 % IIV on F/F (RSE 52.3 %)

    # Residual variability. The paper modeled the residual as 'an additive
    # error model on the log-transformed drug concentrations, being
    # essentially equivalent to an exponential residual error on an
    # arithmetic scale' and reports sigma = 51.6 % (CV) in Table 2, i.e.,
    # the SD on the log scale is 0.516. Per the verification-checklist
    # convention 'NONMEM additive-on-log-scale == proportional in nlmixr2
    # linear space', propSd carries the SD on the log scale directly.
    propSd <- 0.516 ; label("Proportional residual error SD (on log scale)")                                  # Birgersson 2016 Table 2: sigma = 51.6 % (RSE 5.84 %, 95 % CI 44.9-58.1)
  })

  model({
    # Individual PK parameters. No covariates retained in the final model
    # (Birgersson 2016 Results: 'Bodyweight, as an allometric function, did
    # not improve the model fit and was not retained in the final model.
    # ... No other tested covariates were significant in the stepwise
    # covariate approach.').
    cl  <- exp(lcl  + etalcl)
    vc  <- exp(lvc)
    mtt <- exp(lmtt + etalmtt)

    # Transit-absorption chain rate (Savic 2007 parameterization, identical
    # to the convention used by Birgersson 2019 artesunate). With nn_fix = 7
    # transit compartments between depot and central, the chain has
    # nn_fix + 1 = 8 first-order steps and ktr = (nn_fix + 1) / MTT.
    ktr <- (nn_fix + 1) / mtt

    # Elimination micro-constant from the single (central) disposition
    # compartment.
    kel <- cl / vc

    # ODE system: depot -> transit1 -> ... -> transit7 -> central ->
    # elimination. Same first-order rate ktr for every absorption step.
    d/dt(depot)    <- -ktr * depot
    d/dt(transit1) <-  ktr * depot    - ktr * transit1
    d/dt(transit2) <-  ktr * transit1 - ktr * transit2
    d/dt(transit3) <-  ktr * transit2 - ktr * transit3
    d/dt(transit4) <-  ktr * transit3 - ktr * transit4
    d/dt(transit5) <-  ktr * transit4 - ktr * transit5
    d/dt(transit6) <-  ktr * transit5 - ktr * transit6
    d/dt(transit7) <-  ktr * transit6 - ktr * transit7
    d/dt(central)  <-  ktr * transit7 - kel * central

    # Bioavailability on the depot compartment. lfdepot is fixed at log(1)
    # so exp(lfdepot) = 1 at the population mean; etalfdepot carries the
    # 34.3 % IIV on F/F.
    f(depot) <- exp(lfdepot + etalfdepot)

    # Plasma concentration. Doses are in mg and V/F is in L, so
    # central / vc has units of mg/L = 1000 ng/mL. The multiplicative 1000
    # converts to the reported ng/mL units (Birgersson 2016 Table 3 and the
    # LC-MS/MS assay range of 1.03-762 ng/mL described in the Drug analysis
    # subsection).
    Cc <- central / vc * 1000

    # Residual error. NONMEM 'additive on log scale' maps to proportional
    # error in nlmixr2's linear space (see references/verification-checklist.md
    # and Birgersson 2019 artesunate for the sibling convention).
    Cc ~ prop(propSd)
  })
}
