# Joint parent-metabolite population PK model of oral artemether and its
# active metabolite dihydroartemisinin (DHA) in pregnant women with
# uncomplicated Plasmodium falciparum malaria in Uganda
# (Tarning 2012, Malaria Journal 11:293; doi:10.1186/1475-2875-11-293).

Tarning_2012_artemether <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for oral artemether and its",
    "active metabolite dihydroartemisinin (DHA) in 21 pregnant women (2nd or 3rd",
    "trimester) with uncomplicated Plasmodium falciparum malaria in Uganda",
    "after the standard fixed-dose oral artemether-lumefantrine regimen",
    "(Tarning 2012). Absorption is flexible: zero-order dissolution into a depot",
    "of duration DUR feeds a 6-compartment transit chain at rate ktr; the same",
    "ktr empties transit6 into central (ka set equal to ktr). Disposition is",
    "1-compartment for both artemether and DHA with complete in-vivo conversion",
    "of artemether to DHA. Relative bioavailability F is fixed at 1 with",
    "log-normal IIV; no statistically significant covariates were retained in",
    "the final model (Methods / Results); a single combined additive residual",
    "on log-transformed plasma concentrations is shared by both species.",
    sep = " "
  )
  reference <- paste(
    "Tarning J, Kloprogge F, Piola P, Dhorda M, Muwanga S, Turyakira E,",
    "Nuengchamnong N, Nosten F, Day NPJ, White NJ, Guerin PJ, Lindegardh N",
    "(2012). Population pharmacokinetics of Artemether and dihydroartemisinin",
    "in pregnant women with uncomplicated Plasmodium falciparum malaria in",
    "Uganda. Malaria Journal 11:293. doi:10.1186/1475-2875-11-293.",
    sep = " "
  )
  vignette <- "Tarning_2012_artemether"
  units    <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  population <- list(
    species         = "human",
    n_subjects      = 21L,
    n_studies       = 1L,
    n_observations  = 316L,
    age_range       = "16-35 years (median 21; Table 1)",
    weight_range    = "49-88 kg (median 55; Table 1)",
    sex_female_pct  = 100,
    disease_state   = paste(
      "Uncomplicated Plasmodium falciparum malaria in the second or third",
      "trimester of pregnancy (estimated gestational age 13-36 weeks,",
      "median 27 weeks; Table 1)."
    ),
    dose_range      = paste(
      "Coartem (Novartis): 20 mg artemether + 120 mg lumefantrine per tablet;",
      "four tablets per dose (80 mg artemether) given twice daily for 3 days",
      "(oral; dose times 0, 8, 24, 36, 48, 60 hours) with 200 mL milk tea to",
      "optimise oral bioavailability of lumefantrine. The PK analysis used",
      "dense post-last-dose samples only (Methods)."
    ),
    regions         = "Uganda (Mbarara National Referral Hospital antenatal clinic)",
    trial_registration = "ClinicalTrials.gov NCT00495508",
    gestational_age_range = "13-36 weeks (median 27 weeks; Table 1)",
    notes           = paste(
      "Demographics from Tarning 2012 Table 1. No statistically significant",
      "covariates retained in the final model; a full-covariate model with",
      "estimated gestational age on CL_ARM, V_ARM, CL_DHA, V_DHA, and MTT",
      "found the EGA effect distributed -7.0% to +5.5% per week (95% CI;",
      "Results). Companion lumefantrine model from the same Mbarara cohort:",
      "modellib('Kloprogge_2013_lumefantrine')."
    )
  )

  ini({
    # Structural population mean parameters from Tarning 2012 Table 2
    # "Population estimate" column. The paper modelled the natural logarithm
    # of the molar plasma concentration in NONMEM; CL/F and V/F estimates
    # are reported in mass-derived units (L/h, L) and therefore translate
    # directly to a mass-units encoding of the model.
    lcl     <- log(875)
    label("Apparent artemether elimination clearance, CL_ARM/F (L/h)")  # Tarning 2012 Table 2: CL_ARM/F = 875 (RSE 18.7%; 95% CI 625-1280)
    lvc     <- log(2160)
    label("Apparent artemether central volume of distribution, V_ARM/F (L)")  # Tarning 2012 Table 2: V_ARM/F = 2160 (RSE 17.4%; 95% CI 1620-3100)
    lcl_dha <- log(468)
    label("Apparent dihydroartemisinin elimination clearance, CL_DHA/F (L/h)")  # Tarning 2012 Table 2: CL_DHA/F = 468 (RSE 10.2%; 95% CI 387-588)
    lvc_dha <- log(57.1)
    label("Apparent dihydroartemisinin central volume of distribution, V_DHA/F (L)")  # Tarning 2012 Table 2: V_DHA/F = 57.1 (RSE 20.1%; 95% CI 41.7-88.8)

    # Absorption: mean transit time of the 6-compartment transit chain and
    # duration of the zero-order dissolution step preceding the chain.
    # Number of transit compartments NN = 6 was fixed by the paper
    # ("Six transit compartments were sufficient", Results, p.7). The
    # absorption rate constant from transit6 -> central was set equal to
    # the transit rate constant ktr ("In the final model, the absorption
    # rate constant was set to be identical to the rate constant between
    # transit compartments", Results, p.8). With (NN + 1) = 7 first-order
    # transitions at rate ktr from depot through transit1..transit6 into
    # central, the mean transit time is MTT = (NN + 1) / ktr, so
    # ktr = 7 / MTT (Savic & Karlsson 2007 convention; same lab-idiom as
    # Kloprogge 2013 lumefantrine and Birgersson 2019 artesunate).
    lmtt    <- log(0.274)
    label("Mean transit time of the 6-compartment transit chain, MTT (h)")  # Tarning 2012 Table 2: MTT = 0.274 (RSE 19.4%; 95% CI 0.174-0.378)
    ldur    <- log(0.687)
    label("Duration of zero-order dissolution into the depot, DUR (h)")  # Tarning 2012 Table 2: DUR = 0.687 (RSE 25.5%; 95% CI 0.380-1.14)

    # Relative bioavailability anchored at 1 (structural, fixed) with
    # log-normal IIV. The paper held F = 1 as a structural reference (no
    # absolute oral bioavailability data available) but estimated IIV on
    # F to capture between-subject variability in absorption.
    lfdepot <- fixed(log(1))
    label("Relative bioavailability F (unitless); fixed at 1")  # Tarning 2012 Table 2: F = "1 (fixed)" (no estimation, no RSE)

    # Inter-individual variability. Tarning 2012 Table 2 reports IIV as
    # %CV with footnote a: "IIV is presented as 100 * sqrt(exp(omega^2) -
    # 1)". The internal variances are recovered by
    # omega^2 = log((CV/100)^2 + 1):
    #   CL_ARM CV 28.0% -> log(0.280^2 + 1) = 0.07549
    #   CL_DHA CV 90.4% -> log(0.904^2 + 1) = 0.59716
    #   MTT    CV 75.2% -> log(0.752^2 + 1) = 0.44819
    #   DUR    CV 151%  -> log(1.510^2 + 1) = 1.18800
    #   F      CV 85.5% -> log(0.855^2 + 1) = 0.54864
    # IIV on V_ARM and V_DHA was fixed to zero in the final model
    # ("IIV for the distribution volume of artemether and dihydroartemisinin
    # were fixed to zero because of poor precision (RSE > 50%)", Results,
    # p.8) so etalvc and etalvc_dha are omitted.
    etalcl     ~ 0.07549   # Tarning 2012 Table 2: BSV on CL_ARM/F = 28.0% CV (RSE 47.6%)
    etalcl_dha ~ 0.59716   # Tarning 2012 Table 2: BSV on CL_DHA/F = 90.4% CV (RSE 39.0%)
    etalmtt    ~ 0.44819   # Tarning 2012 Table 2: BSV on MTT       = 75.2% CV (RSE 39.6%)
    etaldur    ~ 1.18800   # Tarning 2012 Table 2: BSV on DUR       = 151%  CV (RSE 24.1%)
    etalfdepot ~ 0.54864   # Tarning 2012 Table 2: BSV on F         = 85.5% CV (RSE 24.8%)

    # Residual error. The paper modelled the natural logarithm of the
    # molar plasma concentration with an additive error on the log scale
    # ("A combined additive error model for both the drug and the
    # metabolite was sufficient to describe the random residual
    # variability in the data", Results, p.8; Table 2 footnote c). NONMEM
    # "additive on log scale" maps to proportional residual error in
    # nlmixr2's linear-concentration space (see
    # references/parameter-names.md "Residual error"). Table 2 reports
    # the variance on the log scale as sigma = 0.166; the corresponding
    # SD is sqrt(0.166) = 0.4074. The source paper uses a single combined
    # sigma for both observations, but nlmixr2 requires a distinct
    # endpoint parameter per output, so propSd and propSd_dha are
    # declared here with the same starting value (the published sigma).
    # The 23.1% CV IIV on sigma reported in Table 2 is not encoded here
    # (see vignette Assumptions and deviations).
    propSd     <- sqrt(0.166)
    label("Proportional residual SD for artemether plasma concentration (SD on log scale; combined sigma shared with DHA per Table 2)")  # Tarning 2012 Table 2: sigma = 0.166 (variance, RSE 6.87%; 95% CI 0.130-0.221); combined additive-on-log residual for both species
    propSd_dha <- sqrt(0.166)
    label("Proportional residual SD for DHA plasma concentration (SD on log scale; combined sigma shared with artemether per Table 2)")  # Tarning 2012 Table 2: sigma = 0.166 (variance, RSE 6.87%; 95% CI 0.130-0.221); combined additive-on-log residual for both species
  })

  model({
    # Molecular weights (g/mol; paper Methods, Non-compartmental analysis
    # section: MW_ARM = 298.4 g/mol, MW_DHA = 284.3 g/mol). Used at the
    # ARM -> DHA formation step to convert mass-rate of artemether
    # elimination to mass-rate of DHA formation under the source paper's
    # assumption of complete in-vivo 1:1 molar conversion. With doses in
    # mg the state variables central and central_dha are mass amounts (mg);
    # 1 mol ARM -> 1 mol DHA implies a mass-rate factor of MW_DHA / MW_ARM.
    mw_arm <- 298.4
    mw_dha <- 284.3

    # Individual PK parameters. No significant covariates were retained in
    # the final model (Results, "There were no statistically significant
    # covariates in this study"), so no allometric / categorical / linear-
    # deviation scaling is applied. Volumes carry no IIV (fixed to zero in
    # the source paper); only CL_ARM, CL_DHA, MTT, DUR, and F have etas.
    cl     <- exp(lcl     + etalcl)
    vc     <- exp(lvc)
    cl_dha <- exp(lcl_dha + etalcl_dha)
    vc_dha <- exp(lvc_dha)

    # Absorption rate constants. NN = 6 transit compartments fixed by the
    # paper; with the absorption chain depot -> transit1 -> ... -> transit6
    # -> central (NN + 1 = 7 transitions at rate ktr), the mean transit
    # time is MTT = (NN + 1) / ktr so ktr = 7 / MTT. ka was set equal to
    # ktr in the final model (Results, p.8). Zero-order dissolution into
    # depot has duration DUR (delivered via dur(depot) on the dose record).
    mtt <- exp(lmtt + etalmtt)
    ktr <- 7 / mtt
    dur_in <- exp(ldur + etaldur)

    # Elimination micro-constants. Under the paper's assumption of complete
    # in-vivo 1:1 molar conversion of artemether to DHA, all artemether
    # clearance is metabolic conversion: kel = cl / vc carries artemether
    # out of central and into central_dha (with the MW factor). DHA is
    # then eliminated linearly at kel_dha = cl_dha / vc_dha.
    kel     <- cl     / vc
    kel_dha <- cl_dha / vc_dha

    # ODE system: zero-order dose delivery into depot (duration dur_in),
    # six-compartment first-order transit chain at rate ktr, and two
    # parallel one-compartment disposition models for artemether (central)
    # and dihydroartemisinin (central_dha). The dose record sets
    # dur(depot) <- dur_in to produce the zero-order release; depot itself
    # then empties first-order at rate ktr into transit1.
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(transit3)    <-  ktr * transit2 - ktr * transit3
    d/dt(transit4)    <-  ktr * transit3 - ktr * transit4
    d/dt(transit5)    <-  ktr * transit4 - ktr * transit5
    d/dt(transit6)    <-  ktr * transit5 - ktr * transit6
    d/dt(central)     <-  ktr * transit6 - kel * central
    d/dt(central_dha) <-  kel * central * (mw_dha / mw_arm) -
                          kel_dha * central_dha

    # Zero-order dissolution: dose into depot is delivered over duration
    # dur_in (rxode2 reads dur() at solve time and converts a bolus dose
    # to a constant-rate input of that duration).
    dur(depot) <- dur_in

    # Bioavailability applied to the depot compartment. The structural
    # typical value lfdepot is log(1) (F fixed at 1); etalfdepot carries
    # the log-normal IIV around that anchor.
    f(depot) <- exp(lfdepot + etalfdepot)

    # Plasma concentrations in ng/mL. With doses in mg and volumes in L,
    # central / vc has units mg/L = ug/mL; multiply by 1000 to obtain
    # ng/mL, which matches the units reported in Tarning 2012 Tables 2,
    # 4, and 5 (Cmax_ARM = 32.9 ng/mL, Cmax_DHA = 45.2 ng/mL, AUC values
    # in ng*h/mL).
    Cc     <- 1000 * central     / vc
    Cc_dha <- 1000 * central_dha / vc_dha

    # Combined additive-on-log-scale residual maps to proportional residual
    # in nlmixr2's linear-concentration space. The source paper used a
    # single shared sigma for both observations ("combined additive error
    # model"); nlmixr2 requires distinct endpoint parameters per output,
    # so propSd and propSd_dha are declared with the same initial value
    # (sqrt(0.166)) in ini() to reproduce the published combined-error
    # behaviour for forward simulation.
    Cc     ~ prop(propSd)
    Cc_dha ~ prop(propSd_dha)
  })
}
