Lignet_2023_m8891_human <- function() {
  description <- "Predicted-human translational PK/PD model for M8891, a selective and reversible methionine aminopeptidase 2 (MetAP2) inhibitor. One-compartment oral PK whose disposition (CL, Vss) is the mean of three preclinical-to-human scaling methods and whose absorption rate constant ka comes from a GastroPlus PBPK model, coupled to the effect-compartment plus Met-EF1a turnover PD model estimated in Caki-1 xenograft-bearing mice. No human subjects were dosed in this analysis; the model is the forward projection that was used to select the M8891 dose for the Phase Ia study NCT03138538. PD parameters are carried over unchanged from the mouse model on the hypothesis that the same Met-EF1a modulation level is associated with efficacy in humans."
  reference <- paste(
    "Lignet F, Friese-Hamim M, Jaehrling F, El Bawab S, Rohdich F.",
    "Preclinical Pharmacokinetics and Translational",
    "Pharmacokinetic/Pharmacodynamic Modeling of M8891, a Potent and",
    "Reversible Inhibitor of Methionine Aminopeptidase 2.",
    "Pharm Res. 2023;40(12):3011-3023.",
    "doi:10.1007/s11095-023-03611-z.",
    "PD parameters are inherited unchanged from the mouse Caki-1 xenograft",
    "fit in the same paper; see modellib('Lignet_2023_m8891_mouse').",
    sep = " "
  )
  vignette <- "Lignet_2023_m8891"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # See Lignet_2023_m8891_mouse.R for the rationale; `effect` is the paper's
  # effect compartment Ce and `metef1a` is the tumour Met-EF1a turnover pool.
  paper_specific_compartments <- c("metef1a")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    depot   = list(analyte = "M8891", units = "mg", specimen = "administration site", verified = FALSE),
    central = list(analyte = "M8891", units = "mg", specimen = "plasma", verified = FALSE),
    effect  = list(analyte = "Met-EF1a modulation level", units = "mg", specimen = "not applicable", verified = FALSE),
    metef1a = list(analyte = "Met-EF1a", units = "mg", specimen = "not applicable", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed. Lignet 2023 reports the predicted human CL and Vss per kilogram of body weight (0.020 L/h/kg and 0.21 L/kg, Table III), so both scale linearly with WT: cl = 1.4 * (WT/70) L/h and vc = 14.7 * (WT/70) L. The 70 kg reference is the body weight the authors used for the human compartmental PK simulations (Lignet 2023 Supplementary Materials, 'System Parameterization' table: 'Body weight (compartmental PK model) 70 kg'; the accompanying PBPK population was 'HumanAmerican30YO_70 kg'). This is a unit-conversion scaling implied by the published per-kilogram parameterisation, not an allometric exponent estimated from human data.",
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human (predicted; no human subjects dosed in this analysis)",
    n_subjects     = NA_integer_,
    n_studies      = 0L,
    age_range      = NA_character_,
    weight_range   = "70 kg reference body weight",
    weight_median  = "70 kg",
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state  = "Intended population: patients with advanced solid tumours (renal cell carcinoma was the proof-of-concept indication). The projection supported the Phase Ia dose-escalation study NCT03138538.",
    dose_range     = "Oral once-daily dosing. Lignet 2023 Fig. 5 simulates 150 mg QD, the dose predicted to hold tumour Met-EF1a above the 125 ug/mg-protein efficacy target; the supplement additionally simulates 50-400 mg QD (Table S7).",
    regions        = NA_character_,
    notes          = "Human CL and Vss are the mean of three preclinical-to-human scaling methods (Lignet 2023 Table III, bold row: NASfub / TME / SA-RoE for CL and SAfub / Oie-Tozer / human-dog proportionality for Vss), derived from single-dose i.v. PK in NMRI mice, Wistar rats, beagle dogs, and cynomolgus monkeys corrected for plasma protein binding and liver-microsomal CLint. The absorption rate constant and oral bioavailability come from a GastroPlus 9.5 PBPK/ACAT model calibrated on rat and dog data. The PD component (ke0, kin, kout, Imax, IC50) is transferred unchanged from the mouse Caki-1 xenograft fit. Lignet 2023 Discussion notes that the observed Phase Ia terminal half-life was about 30 h, roughly fourfold longer than the 7.3 h projected here, so this model over-predicts clearance and hence the required dose; the recommended Phase II dose was 35 mg, not the 150 mg predicted here (Carducci et al. 2023). The model is retained as the published translational projection, not as a description of observed human PK."
  )

  ini({
    # --- PK. ---
    # Lignet 2023 Table III bold row: mean predicted human CL = 0.020 L/h/kg
    # and Vss = 0.21 L/kg; at the 70 kg reference weight used for the human
    # compartmental simulations these are 1.4 L/h and 14.7 L. They are
    # entered here as APPARENT oral parameters (CL/F and V/F), matching
    # Lignet 2023 Eq. (3), in which the dose is not multiplied by F and the
    # bioavailability appears only inside the apparent volume V/F. Doing so
    # reproduces the paper's headline output: a steady-state trough of
    # 1584 ng/mL at 150 mg QD against the 1500 ng/mL reported in the
    # Results and drawn in Fig. 5a, and 1540 ng/mL (Cmin) in Supplementary
    # Table S7. Applying the separately-predicted F (>= 60% in the main
    # text, 70.8% at 150 mg in Supplementary Table S7) on top of these
    # values would give roughly 950-1120 ng/mL and would NOT reproduce the
    # published trough; see the vignette's "Assumptions and deviations".
    lka <- fixed(log(0.35)); label("Absorption rate constant ka (1/h)")                           # Lignet 2023 Results, "PBPK Modeling": "the rate of absorption (ka) was estimated to be 0.35 h-1" (GastroPlus prediction, not estimated from human data)
    lcl <- log(1.4); label("Apparent oral clearance CL/F (L/h) at 70 kg")                         # Lignet 2023 Table III bold row: mean CL = 0.020 +/- 0.002 L/h/kg; 0.020 * 70 kg = 1.4 L/h
    lvc <- log(14.7); label("Apparent oral volume of distribution V/F (L) at 70 kg")              # Lignet 2023 Table III bold row: mean Vss = 0.21 +/- 0.017 L/kg; 0.21 * 70 kg = 14.7 L

    # --- PD, transferred unchanged from the mouse Caki-1 xenograft fit. ---
    # Lignet 2023 Human PK and PD Prediction: "The previously established
    # xenograft PK/PD model, with the PK part replaced by the predicted human
    # PK model, was used to simulate Met-EF1a modulation after M8891 dosing
    # in humans." All five values are Lignet 2023 Table IV (mouse).
    lke0 <- log(0.0566); label("Effect-compartment equilibration rate constant ke0 (1/h)")        # Lignet 2023 Table IV (mouse estimate, carried over): ke0 = 0.0566 1/h
    lkin <- log(29.1); label("Zero-order Met-EF1a synthesis rate kin (ug/mg protein/h)")          # Lignet 2023 Table IV (mouse estimate, carried over): Kin = 29.1 ug/mg/h
    lkout <- log(1.45); label("First-order Met-EF1a degradation rate constant kout (1/h)")        # Lignet 2023 Table IV (mouse estimate, carried over): Kout = 1.45 1/h
    limax <- log(0.91); label("Maximum fractional inhibition of Met-EF1a degradation Imax (unitless)") # Lignet 2023 Table IV (mouse estimate, carried over): Imax = 0.91
    lic50 <- log(0.340); label("Effect-site concentration producing half-maximal inhibition IC50 (mg/L)") # Lignet 2023 Table IV (mouse estimate, carried over): IC50 = 340 ng/mL = 0.340 mg/L

    # --- Residual error. ---
    # This model is a deterministic forward projection: no human data were
    # fitted, so no residual error or inter-individual variability exists.
    # The proportional form is retained from the mouse fit with a magnitude
    # fixed at zero; see the vignette's "Assumptions and deviations".
    propSd <- fixed(0); label("Proportional residual error on plasma concentration (fraction)")   # not applicable: forward projection, no human data fitted
    propSd_MetEF1a <- fixed(0); label("Proportional residual error on Met-EF1a (fraction)")       # not applicable: forward projection, no human data fitted
  })

  model({
    ka <- exp(lka)
    # CL and Vss are published per kilogram of body weight, so both scale
    # linearly with WT about the 70 kg reference used by the authors.
    cl <- exp(lcl) * (WT / 70)
    vc <- exp(lvc) * (WT / 70)
    ke0 <- exp(lke0)
    kin <- exp(lkin)
    kout <- exp(lkout)
    imax <- exp(limax)
    ic50 <- exp(lic50)

    kel <- cl / vc

    # One-compartment oral PK, Lignet 2023 Eqs. (3)-(4); doses in mg and vc
    # in L give Cc in mg/L (= ug/mL = 1000 ng/mL, the unit of Fig. 5a).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    Cc <- central / vc

    # Effect compartment, Lignet 2023 Eq. (5).
    d/dt(effect) <- ke0 * (Cc - effect)

    # Met-EF1a turnover with inhibition of degradation, Lignet 2023 Eq. (6),
    # starting from the undrugged steady state kin/kout = 20.1 ug/mg protein
    # (the t = 0 value in Lignet 2023 Fig. 5b).
    d/dt(metef1a) <- kin - kout * metef1a * (1 - imax * effect / (effect + ic50))
    metef1a(0) <- kin / kout

    MetEF1a <- metef1a

    Cc ~ prop(propSd)
    MetEF1a ~ prop(propSd_MetEF1a)
  })
}
