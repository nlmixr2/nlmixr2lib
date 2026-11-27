Drobnjak_2018_pp13_rabbit <- function() {
  description <- paste(
    "Preclinical (rabbit). Two-compartment pharmacokinetic model of",
    "recombinant human placental protein 13 (PP13, galectin-13) in serum",
    "after a single intravenous bolus (1.3, 2.6 or 12.8 ng/kg) or",
    "subcutaneous injection (12.8 ng/kg) in nonpregnant female New Zealand",
    "White rabbits (Drobnjak 2018). Each animal was fitted individually and",
    "the parameters are the dose-group means of the individual estimates, so",
    "the model carries no between-animal variability. The disposition is",
    "dose dependent (clearance rises from 14 to 34 mL/h between 1.3 and",
    "12.8 ng/kg), so the dose group selects the parameter set via",
    "DOSE_PP13_NGKG. Subcutaneous dosing uses first-order absorption",
    "(ka = 1.78 1/h) and bioavailability 57.1 percent on the 12.8 ng/kg",
    "intravenous disposition. Clearance and volumes are per animal (mean",
    "3.9 kg), dose in pg of PP13 per animal, concentration in pg/mL."
  )
  reference <- paste(
    "Drobnjak T, Jonsdottir AM, Helgadottir H, Runolfsdottir MS, Meiri H,",
    "Sammar M, Osol G, Mandala M, Huppertz B, Gizurarson S (2018).",
    "Pharmacokinetics of placental protein 13 after intravenous and",
    "subcutaneous administration in rabbits. Drug Des Devel Ther",
    "12:1977-1983. doi:10.2147/DDDT.S167926"
  )
  vignette <- "Drobnjak_2018_pp13_rabbit"
  units <- list(time = "h", dosing = "pg", concentration = "pg/mL")

  covariateData <- list(
    DOSE_PP13_NGKG = list(
      description = paste(
        "PP13 dose level of the animal's treatment group (ng per kg body",
        "weight). Selects the dose-group parameter set: below 1.84 ng/kg the",
        "1.3 ng/kg set, 1.84 to below 5.77 ng/kg the 2.6 ng/kg set, and",
        "5.77 ng/kg or above the 12.8 ng/kg set (breakpoints are the",
        "geometric midpoints between the studied levels)."
      ),
      units = "ng/kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Studied levels 1.3, 2.6 and 12.8 ng/kg IV and 12.8 ng/kg SC",
        "(Methods, Experimental procedure). The paper fitted each group",
        "separately and gives no dose-exposure function, so doses between",
        "or outside the studied levels take the nearest group's parameters;",
        "this is a lookup, not an interpolation. The dose amount itself is",
        "carried by the event table in pg per animal."
      ),
      source_name = "dose group (Group 1-4)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "PP13", units = "pg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "PP13", units = "pg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "PP13", units = "pg", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "rabbit (New Zealand White)",
    n_subjects = 12L,
    n_studies = 1L,
    weight_range = "3.9 +/- 0.3 kg (mean +/- SD)",
    sex_female_pct = 100,
    disease_state = "healthy nonpregnant female rabbits",
    dose_range = paste(
      "Single IV bolus of 1.3, 2.6 or 12.8 ng/kg (1 mL of 5, 10 or 50 ng/mL,",
      "i.e. about 5, 10 or 50 ng per animal) into the marginal ear vein, or",
      "a single SC injection of 12.8 ng/kg at the neck."
    ),
    regions = "Iceland (University of Iceland, Reykjavik)",
    notes = paste(
      "Twelve rabbits rotated through four groups of six animals each with",
      "a one-week washout (Methods). Serum sampled at 0, 10, 30, 60 and 90",
      "min and 2, 4, 6 and 24 h; PP13 by ELISA (Cusabio CSB-E12733h,",
      "working range 2.5-1000 pg/mL). Two-compartment fits per animal in",
      "GraphPad Prism / Excel; Table 1 reports group means +/- SD."
    )
  )

  ini({
    # Table 1 (means of individual two-compartment fits, 6 rabbits per group).
    # Table 1 prints Cltot as 'mL/kg/h' and Vd as 'L', but both are per-animal
    # quantities: Dose/AUC with the per-animal dose of 5, 10 and 50 ng
    # (1 mL of 5, 10 or 50 ng/mL) gives 15.9, 20.2 and 36.8 mL/h against
    # Cltot 14.0, 20.3 and 33.5, and the Table 1 'AUC/dose (h/mL)' column
    # (0.07, 0.05, 0.03) is AUC over the per-animal dose. The central volume
    # is CL / ke from the same table; the printed 'Vd' (110, 132, 228) is the
    # same quantity averaged differently (vignette, Assumptions). The
    # micro-constants follow in model() from alpha * beta = kel * k21 and
    # alpha + beta = kel + k12 + k21.

    # Group 1, 1.3 ng/kg IV
    lcl_dose1p3 <- log(14.0); label("Clearance, 1.3 ng/kg group (mL/h per animal)") # Table 1 Group 1 Cltot = 14.0 +/- 1.9
    lvc_dose1p3 <- log(93.333); label("Central volume, 1.3 ng/kg group (mL per animal)") # derived: Table 1 Group 1 Cltot / ke = 14.0 / 0.15 = 93.333 mL (printed Vd 110)
    lalpha_dose1p3 <- log(0.36); label("Distribution-phase rate constant alpha, 1.3 ng/kg group (1/h)") # Table 1 Group 1 alpha = 0.36 +/- 0.14
    lbeta_dose1p3 <- log(0.09); label("Terminal-phase rate constant beta, 1.3 ng/kg group (1/h)") # Table 1 Group 1 beta = 0.09 +/- 0.03

    # Group 2, 2.6 ng/kg IV
    lcl_dose2p6 <- log(20.3); label("Clearance, 2.6 ng/kg group (mL/h per animal)") # Table 1 Group 2 Cltot = 20.3 +/- 3.5
    lvc_dose2p6 <- log(126.875); label("Central volume, 2.6 ng/kg group (mL per animal)") # derived: Table 1 Group 2 Cltot / ke = 20.3 / 0.16 = 126.875 mL (printed Vd 132)
    lalpha_dose2p6 <- log(0.50); label("Distribution-phase rate constant alpha, 2.6 ng/kg group (1/h)") # Table 1 Group 2 alpha = 0.50 +/- 0.13
    lbeta_dose2p6 <- log(0.13); label("Terminal-phase rate constant beta, 2.6 ng/kg group (1/h)") # Table 1 Group 2 beta = 0.13 +/- 0.04

    # Group 3, 12.8 ng/kg IV (also the disposition for the 12.8 ng/kg SC group)
    lcl_dose12p8 <- log(33.5); label("Clearance, 12.8 ng/kg group (mL/h per animal)") # Table 1 Group 3 Cltot = 33.5 +/- 2.6
    lvc_dose12p8 <- log(223.333); label("Central volume, 12.8 ng/kg group (mL per animal)") # derived: Table 1 Group 3 Cltot / ke = 33.5 / 0.15 = 223.333 mL (printed Vd 228)
    lalpha_dose12p8 <- log(0.60); label("Distribution-phase rate constant alpha, 12.8 ng/kg group (1/h)") # Table 1 Group 3 alpha = 0.60 +/- 0.23
    lbeta_dose12p8 <- log(0.04); label("Terminal-phase rate constant beta, 12.8 ng/kg group (1/h)") # Table 1 Group 3 beta = 0.04 +/- 0.03

    # Group 4, 12.8 ng/kg SC
    lka <- log(1.78); label("First-order SC absorption rate constant (1/h)") # Table 1 Group 4 kabs = 1.78 +/- 0.34
    lfdepot <- log(0.571); label("SC bioavailability (fraction)") # Table 1 Group 4 F% = 57.1 +/- 3.2; Results 'F ... 57.1% +/- 3.2%'

    # Individual fits only: no residual-error model was reported.
    propSd <- fixed(0); label("Proportional residual error (fraction; not reported)") # not reported; ELISA intra-/inter-assay CV 10%/11% (Methods) is assay precision, not model error
  })

  model({
    # Dose-group lookup (no dose-exposure function is published)
    grp1 <- DOSE_PP13_NGKG < 1.84
    grp3 <- DOSE_PP13_NGKG >= 5.77
    grp2 <- 1 - grp1 - grp3

    cl <- exp(grp1 * lcl_dose1p3 + grp2 * lcl_dose2p6 + grp3 * lcl_dose12p8)
    vc <- exp(grp1 * lvc_dose1p3 + grp2 * lvc_dose2p6 + grp3 * lvc_dose12p8)
    alpha <- exp(grp1 * lalpha_dose1p3 + grp2 * lalpha_dose2p6 + grp3 * lalpha_dose12p8)
    beta <- exp(grp1 * lbeta_dose1p3 + grp2 * lbeta_dose2p6 + grp3 * lbeta_dose12p8)
    ka <- exp(lka)
    fdepot <- exp(lfdepot)

    # Micro-constants of the two-compartment disposition
    kel <- cl / vc
    k21 <- alpha * beta / kel
    k12 <- alpha + beta - kel - k21
    q <- k12 * vc
    vp <- q / k21

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    f(depot) <- fdepot

    # pg / mL
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
