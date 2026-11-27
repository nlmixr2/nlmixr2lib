Soliman_2025_omeprazole <- function() {
  description <- paste(
    "Semi-mechanistic intragastric-pH model for omeprazole in adults and",
    "children (Soliman 2025), the PD layer of the Soliman 2025 omeprazole",
    "PBPK/PD analysis. Plasma omeprazole irreversibly inactivates the gastric",
    "H+/K+-ATPase (relative activity E/E0, turnover kdeg, second-order",
    "inactivation kd * Cp acting after a 0.396 h lag); the active pump fraction",
    "scales a zero-order intragastric H+ secretion rate that carries a narrow",
    "night-time circadian dip, and H+ is cleared first-order (indirect-response",
    "model of Liu et al. 2016 for esomeprazole, with the omeprazole kd of",
    "Katashima et al. 1998). Meals dilute the gastric H+: each meal volume in",
    "mL is added to a dimensionless gastric-dilution state that relaxes back to",
    "1 at the gastric-emptying rate kfe, and the observed H+ is the secreted H+",
    "divided by that state; pH is capped at 7.4. Only this PD layer is",
    "reproduced here: the PK layer is a GastroPlus whole-body PBPK model whose",
    "organ volumes, blood flows, partition coefficients and enzyme abundances",
    "are platform database outputs the paper does not print, so the omeprazole",
    "plasma concentration is supplied by the user as the time-varying",
    "covariate CP_OMEPRAZOLE_NGML and meals are supplied as dose records into",
    "the gastric_volume compartment. The baseline H+ concentration defaults to",
    "the adult value; the paper's pediatric baselines are listed in the",
    "population notes and are applied by overriding lrbase. No",
    "inter-individual variability or residual error is reported, so both are",
    "zero.",
    sep = " "
  )
  reference <- paste(
    "Soliman A, Rodriguez-Vera L, Alarcia-Lacalle A, et al. Leveraging",
    "Omeprazole PBPK/PD Modeling to Inform Drug-Drug Interactions and Specific",
    "Recommendations for Pediatric Labeling. Pharmaceutics. 2025;17(3):373.",
    "doi:10.3390/pharmaceutics17030373. PD model structure and the kdeg,",
    "kout, circadian and lag values from Liu D, Yang H, Jiang J, et al.",
    "Pharmacokinetic and Pharmacodynamic Modeling Analysis of Intravenous",
    "Esomeprazole in Healthy Volunteers. J Clin Pharmacol. 2016;56:816-826",
    "(doi:10.1002/jcph.733); omeprazole kd from Katashima M, et al. Eur J Drug",
    "Metab Pharmacokinet. 1998;23:19-26 (doi:10.1007/BF03189822), both as",
    "tabulated in Soliman 2025 supplementary Table S1.",
    sep = " "
  )
  vignette <- "Soliman_2025_omeprazole"
  units <- list(
    time = "h",
    dosing = "mL of ingested meal / milk volume, dosed into the gastric_volume compartment (omeprazole itself enters only through the CP_OMEPRAZOLE_NGML covariate)",
    concentration = "ng/mL (omeprazole plasma concentration via CP_OMEPRAZOLE_NGML); hplus is in mmol/L and gastric_ph in pH units"
  )
  dosing <- "gastric_volume"

  # Mechanistic PD states that have no canonical compartment name: the
  # relative proton-pump activity, the secreted intragastric H+ concentration
  # and the dimensionless meal-dilution factor.
  paper_specific_compartments <- c("pump", "hplus", "gastric_volume")

  compartmentData <- list(
    pump = list(
      analyte = "active gastric H+/K+-ATPase relative to baseline (E/E0)",
      units = "(unitless)",
      specimen = "not applicable",
      verified = TRUE
    ),
    hplus = list(
      analyte = "secreted intragastric H+ before meal dilution",
      units = "mmol/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    gastric_volume = list(
      analyte = "meal-dilution factor (fasting = 1; meals add their volume in mL)",
      units = "(unitless)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    CP_OMEPRAZOLE_NGML = list(
      description = "Instantaneous omeprazole (racemic, total) plasma concentration, supplied as the time-varying driver of H+/K+-ATPase inactivation",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters the pump equation as kd * (E/E0) * Cp with kd in 1/h per",
        "ng/mL (Soliman 2025 supplementary Equation S6 and Table S1). The",
        "0.396 h lag between plasma concentration and intragastric H+ is",
        "applied inside the model, so supply the undelayed plasma",
        "concentration. In the source study Cp is the profile predicted by the",
        "authors' GastroPlus 9.8.2 whole-body PBPK model of omeprazole and its",
        "hydroxy and sulphone metabolites (CYP2C19 / CYP3A4 / CYP2C9 / CYP3A7",
        "Michaelis-Menten metabolism with auto-inhibition; main-text Table 2",
        "and Table 3); that model is not reproduced in this file and cannot be",
        "reproduced from the publication. Users supply the column from their",
        "own omeprazole PK model or from observed concentrations; set it to 0",
        "for drug-free records. The metabolites do not inhibit acid secretion",
        "and are not inputs."
      ),
      source_name = "Cp (Soliman 2025 supplementary Equation S6)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 6L,
    age_range = "simulated reference individuals of 4 months, 8 months, 2 years, 6 years and 37 years; the PBPK layer was verified against clinical studies from 0.3 to 45 years",
    weight_range = "simulated population means 6.6 kg (4 months), 9.11 kg (8 months), 12.4 kg (2 years), 23.5 kg (6 years) and 85.27 kg (37 years) (supplementary Table S6)",
    sex_female_pct = NA_real_,
    disease_state = "healthy virtual subjects; the PD parameters come from esomeprazole and omeprazole studies in healthy adult volunteers",
    dose_range = "omeprazole 20 mg once daily (adults) and 0.7 mg/kg once daily (children) for 8 days in the PD simulations",
    regions = "not reported (clinical PK studies digitised from the published literature)",
    notes = paste(
      "The PD model was not fitted by Soliman 2025: its structure and",
      "parameters were taken from Liu et al. 2016 (esomeprazole IV in healthy",
      "volunteers) and Katashima et al. 1998 (omeprazole kd), and applied",
      "unchanged across ages ('the dose-exposure-response relationship was",
      "assumed to be consistent across ages', Methods 2.5) except for the",
      "baseline intragastric H+ concentration and the meal pattern. Baseline",
      "H+ (supplementary Table S1): 40.1 mmol/L (pH 1.4) adults, 31.62",
      "(pH 1.5) 6 years, 16.59 (pH 1.78) 2 years, 1.2 (pH 2.92) 8 months, 0.79",
      "(pH 3.1) 4 months. Meal volumes (mL, Table S1) for adults / 6 years /",
      "2 years / 8 months: breakfast at 1 h 540/410/290/330, snack at 3 h",
      "270/180/120/240, lunch at 5 h 600/440/360/330, dinner at 11 h",
      "630/380/330/510, with a second snack at 8 h whose volume Table S1 does",
      "not list; 8-month infants also get 150 mL milk at 13 h, and 4-month",
      "infants get 120 mL milk every 4 h from 1 h (six feeds). Meals and milk",
      "are taken over 30 min and snacks as boluses (supplementary Equations",
      "S3-S5); clock times are hours after the once-daily morning dose. The",
      "clinical PK studies behind the PBPK layer (main-text Table 1) enrolled",
      "79 subjects."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Soliman 2025 supplementary section 'Semi-mechanistic model to
    # predict intragastric pH', Equations S1-S8, parameter values from
    # supplementary Table S1. Every value was taken from the literature
    # (Liu 2016, Katashima 1998, GastroPlus volumes) rather than estimated
    # in Soliman 2025, so all are fixed.
    #
    #   kin(0) = kout * BASE * (1 + MA / ((MTmax/MW)^4 + 1))        (S1)
    #   kin(t) = kin(0) / (1 + MA / (((T - MTmax)/MW)^4 + 1))        (S2)
    #   dVEf/dt = meal inputs + kfe - kfe * VEf                      (S3-S5)
    #   d(E/E0)/dt = kdeg - kdeg * E/E0 - kd * E/E0 * Cp             (S6)
    #   dH/dt  = kin(t) * E/E0 - kout * H                            (S7)
    #   Hobs   = H / VEf                                             (S8)
    # ------------------------------------------------------------------

    lkdeg <- fixed(log(0.00827))
    label("First-order turnover rate constant of the H+/K+-ATPase kdeg (1/h)")
    # Table S1 'kdeg (1/h)' = 0.00827 [44] (Liu 2016); restated in the text below Equation S8

    lkinact <- fixed(log(0.00388))
    label("Second-order irreversible H+/K+-ATPase inactivation rate constant kd of omeprazole (1/h per ng/mL)")
    # Table S1 'kd (1/h/ng/mL)' = 0.00388 [67] (Katashima 1998); restated in the text below Equation S8

    lrbase <- fixed(log(40.1))
    label("Baseline intragastric H+ concentration BASE, adults (mmol/L)")
    # Table S1 'Base (mM)', 37 year-old adults = 40.1 [44] (pH 1.4); pediatric values in population$notes

    lkout <- fixed(log(3.83))
    label("First-order elimination rate constant of intragastric H+ kout (1/h)")
    # Table S1 'kout (1/h)' = 3.83 (Liu 2016); restated in the supplementary text after Equation S2

    lamp_night <- fixed(log(537))
    label("Amplitude MA of the night-time dip in H+ secretion (unitless)")
    # Table S1 'MA (unit)' = 537. The supplementary text instead quotes 398; Table S1 is used because only it reproduces Figure 5 (see vignette)

    lwidth_night <- fixed(log(0.0580))
    label("Width MW of the night-time dip in H+ secretion (h)")
    # Table S1 'MW (h)' = 0.0580. The supplementary text instead quotes 1.25 h; Table S1 is used because only it reproduces Figure 5 (see vignette)

    ltpeak_night <- fixed(log(20.1))
    label("Time MTmax of the night-time dip in H+ secretion (h after the daily dose)")
    # Table S1 'MTmax (h)' = 20.1. The supplementary text instead quotes 22.4 h; Table S1 is used because only it reproduces Figure 5 (see vignette)

    lkfe <- fixed(log(0.694))
    label("First-order gastric-emptying rate constant of the meal-dilution state kfe (1/h)")
    # Table S1 'kfe (1/h)' = 0.694; restated in the supplementary text before Equation S3

    ltlag <- fixed(log(0.396))
    label("Lag between plasma omeprazole and its effect on intragastric H+ tlag (h)")
    # Table S1 'tlag (h)' = 0.396 (Liu 2016); restated in the text below Equation S8

    addSd_gastric_ph <- fixed(0)
    label("Additive residual SD on intragastric pH (pH units; zero - no residual-error model is reported)")
  })

  model({
    kdeg <- exp(lkdeg)
    kinact <- exp(lkinact)
    rbase <- exp(lrbase)
    kout <- exp(lkout)
    amp_night <- exp(lamp_night)
    width_night <- exp(lwidth_night)
    tpeak_night <- exp(ltpeak_night)
    kfe <- exp(lkfe)
    tlag <- exp(ltlag)

    # Circadian H+ secretion (Equations S1-S2). T is the time of day in hours
    # after the once-daily dose, so model time 0 must be a dosing time.
    tday <- time - 24 * floor(time / 24)
    kin0 <- kout * rbase * (1 + amp_night / ((tpeak_night / width_night)^4 + 1))
    kin <- kin0 / (1 + amp_night / (((tday - tpeak_night) / width_night)^4 + 1))

    # Relative H+/K+-ATPase activity E/E0 (Equation S6). The pump equation
    # depends on time only through Cp, so evaluating it on the undelayed Cp
    # and delaying the state by tlag equals driving it with Cp(t - tlag).
    d/dt(pump) <- kdeg - kdeg * pump - kinact * pump * CP_OMEPRAZOLE_NGML
    pump(0) <- 1
    pump_lag <- delay(pump, tlag)

    # Secreted intragastric H+ in mmol/L (Equation S7); kin(0) = kout * BASE
    # puts the drug-free, fasting state at BASE at time 0.
    d/dt(hplus) <- kin * pump_lag - kout * hplus
    hplus(0) <- rbase

    # Meal dilution (Equations S3-S5). Meal and milk volumes in mL arrive as
    # doses (0.5 h infusions for meals and milk, boluses for snacks) and are
    # added to a dimensionless state that relaxes back to 1.
    d/dt(gastric_volume) <- kfe - kfe * gastric_volume
    gastric_volume(0) <- 1

    # Observed H+ (Equation S8) and pH, capped at the physiological 7.4.
    hplus_obs <- hplus / gastric_volume
    gastric_ph <- -log10(hplus_obs / 1000)
    if (gastric_ph > 7.4) {
      gastric_ph <- 7.4
    }

    gastric_ph ~ add(addSd_gastric_ph)
  })
}
