Jia_2015_unfractionatedHeparin <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order elimination for",
    "unfractionated heparin (UFH) administered as multiple intravenous bolus",
    "injections during cardiopulmonary bypass (CPB) in adult Chinese cardiac",
    "surgery patients (Jia 2015). Plasma UFH exposure was inferred from",
    "anti-FIIa chromogenic activity. No covariates were retained in the final",
    "model (age, body weight, and sex were tested via forward inclusion /",
    "backward elimination and none met the p < 0.001 retention threshold).",
    "Concentrations are reported in IU/mL of anti-FIIa activity; doses are in",
    "IU (1 mg UFH = 125 IU). The published model also describes instantaneous",
    "neutralization of central-compartment UFH at protamine sulfate dosing",
    "(see vignette for the simulation pattern); the structural ODEs here are",
    "the standard two-compartment IV bolus form."
  )
  reference <- paste(
    "Jia Z, Tian G, Ren Y, Sun Z, Lu W, Hou X.",
    "Pharmacokinetic model of unfractionated heparin during and after",
    "cardiopulmonary bypass in cardiac surgery.",
    "J Transl Med. 2015 Feb 1;13:45.",
    "doi:10.1186/s12967-015-0404-5.",
    sep = " "
  )
  vignette <- "Jia_2015_unfractionatedHeparin"
  units <- list(
    time          = "h",
    dosing        = "IU",
    concentration = "IU/mL"
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 32L,
    n_studies      = 1L,
    age_range      = "18-74 years",
    age_median     = "53.4 years (mean)",
    weight_range   = "41-82 kg",
    weight_median  = "66 kg (mean)",
    sex_female_pct = 59.4,
    race_ethnicity = "Chinese (single-center study at Beijing Anzhen Hospital, Capital Medical University)",
    disease_state  = paste(
      "Adults undergoing cardiopulmonary bypass for cardiac surgery.",
      "Subjects with dysfunction of the kidney, liver, or blood coagulation",
      "were excluded prior to enrollment."
    ),
    dose_range     = paste(
      "Initial intravenous bolus of UFH at 375 IU/kg (= 3 mg/kg) before CPB",
      "(mean first dose 24,805 IU; range 18,750-31,250). A second bolus of",
      "1 mg/kg (mean 8,516 IU; range 6,250-10,000) was added to the priming",
      "fluid. Additional irregular boluses were given during CPB as needed.",
      "Mean total UFH dose 34,023 IU (range 25,000-45,000)."
    ),
    cpb_time       = "2.04 h mean (range 0.95-3.29 h)",
    regions        = "China (Beijing)",
    notes          = paste(
      "Baseline demographics from Jia 2015 Table 1 (n = 41 enrolled; 32",
      "completed the protocol and contributed to the final model). The Table",
      "1 sex row lists 19 females and 13 males, which sums to 32 completers.",
      "Sampling: before protamine neutralization, plasma was drawn at 30-min",
      "intervals; after protamine neutralization, samples were drawn at 2, 4,",
      "8, 12, and 24 h. Plasma anti-FIIa activity was measured with the",
      "Heparin Chromogenic Activity Kit 820 on the ACL-TOP coagulation",
      "platform; pre-neutralization samples were diluted 1:29 in normal",
      "pooled platelet-poor plasma to bring activity into the assay range",
      "(0.0-0.6 IU/mL). NONMEM VII level 2.0 with FOCE, PsN 3.4.0. Estimation",
      "via the $PRED block with explicit Laplace-transformed solution.",
      "Validation by bootstrap (1000 resamples) and VPC (1000 simulations)."
    )
  )

  ini({
    # Structural parameters from Jia 2015 Table 2 (final-model estimates).
    # The paper parameterizes the model in clearance and volume rather than
    # rate constants (Results: "the model was parameterized in terms of
    # volume of distribution and clearance rather than rate constants").
    # Units: CL and Q in L/h, Vc and Vp in L. Doses are in IU; concentration
    # Cc is converted to IU/mL in model() via division by 1000.
    lcl <- log(1.18);   label("Clearance (L/h)")                                  # Jia 2015 Table 2: CL = 1.18 L/h (RSE 7.25%)
    lvc <- log(3.04);   label("Central volume of distribution (L)")               # Jia 2015 Table 2: V_UFH-C = 3.04 L (RSE 8.09%)
    lq  <- log(0.171);  label("Inter-compartmental clearance (L/h)")              # Jia 2015 Table 2: Q_UFH = 0.171 L/h (RSE 14.40%)
    lvp <- log(8.01);   label("Peripheral volume of distribution (L)")            # Jia 2015 Table 2: V_UFH-P = 8.01 L (RSE 31.8%)

    # Inter-individual variability. The paper uses an exponential IIV model
    # (Eq. 15: P_i = P_pop * exp(eta_i), eta_i ~ N(0, omega^2)). Because the
    # structural parameters above are already log-transformed, the variance
    # of each eta on the log scale equals the reported omega^2 directly.
    # omega^2(V_UFH-P) was fixed at 0 in the final model (Table 2), so no
    # eta is applied to lvp.
    etalcl ~ 0.122   # Jia 2015 Table 2: omega^2(CL) = 0.122 (RSE 77.7%)
    etalvc ~ 0.105   # Jia 2015 Table 2: omega^2(V_UFH-C) = 0.105 (RSE 48.7%)
    etalq  ~ 0.0978  # Jia 2015 Table 2: omega^2(Q_UFH) = 0.0978 (RSE 76.9%)

    # Residual error: the paper describes a hybrid additive + proportional
    # residual structure (Eq. 16). In the final model, the additive variance
    # was fixed at 0 (Table 2: sigma^2_add = 0 FIX), so only the proportional
    # term contributes. propSd is the SD on the linear scale of the
    # proportional component: propSd = sqrt(sigma^2_prop) = sqrt(0.139).
    propSd <- sqrt(0.139)   # Jia 2015 Table 2: sigma^2_prop(UFH) = 0.139 (RSE 14.2%); SD = sqrt(0.139) ~= 0.3728
    label("Proportional residual error SD (fraction)")
  })
  model({
    # Individual structural PK parameters. No covariates were retained in
    # the final model (Results: "None of the tested covariates significantly
    # decreased the objective function"), so the per-subject value depends
    # only on the typical value and the eta.
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq  + etalq)
    vp <- exp(lvp)

    # Micro-constants for the two-compartment IV bolus model (Jia 2015
    # Eq. 1-2). The amount in the central compartment is dosed directly
    # via an IV bolus event on the central compartment.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Observation. The state `central` carries the amount of UFH in IU and
    # vc has units of L, so `central / vc` is in IU/L. The paper reports
    # plasma anti-FIIa activity in IU/mL, so divide by 1000 to convert.
    # Worked sanity check: an initial bolus of 24,805 IU into vc = 3.04 L
    # gives Cc(0+) = 24,805 / 3.04 / 1000 = 8.16 IU/mL, consistent with
    # the median plateau plasma anti-FIIa activity range of 2-19 IU/mL
    # reported during CPB (Discussion).
    Cc <- central / vc / 1000
    Cc ~ prop(propSd)
  })
}
