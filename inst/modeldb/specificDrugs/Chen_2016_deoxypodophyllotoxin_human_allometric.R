Chen_2016_deoxypodophyllotoxin_human_allometric <- function() {
  description <- paste0(
    "Two-compartment human deoxypodophyllotoxin (DPT) model obtained by ",
    "interspecies allometric scaling across mouse, rat and monkey (Chen et ",
    "al. 2016, Front Pharmacol). This is the paper's SECOND, independent ",
    "human prediction: the plasma concentration-time data of all dose groups ",
    "in the three species were collapsed onto a Dedrick plot (physiological ",
    "time t' = t / W^(b1 - b2), normalised concentration C' = C / (D / ",
    "W^b1)), the pooled profile was fitted with a two-compartment model, and ",
    "the fit was reverse-transformed to a 70 kg human. Dog was excluded from ",
    "the scaling because its clearance and steady-state volume are lower in ",
    "absolute terms than the monkey's despite a larger body weight. No human ",
    "data were fitted; the model is a forward prediction whose only ",
    "validation is that it agrees with the paper's whole-body PBPK ",
    "prediction (Chen_2016_deoxypodophyllotoxin_human_pbpk). Deterministic: ",
    "the publication reports no inter-individual variance and no ",
    "residual-error magnitude, so the model is intended for typical-value ",
    "simulation."
  )
  reference <- paste0(
    "Chen Y, Zhao K, Liu F, Xie Q, Zhong Z, Miao M, Liu X, Liu L. ",
    "Prediction of Deoxypodophyllotoxin Disposition in Mouse, Rat, ",
    "Monkey, and Dog by Physiologically Based Pharmacokinetic Model ",
    "and the Extrapolation to Human. Front Pharmacol. 2016;7:488. ",
    "doi:10.3389/fphar.2016.00488"
  )
  vignette <- "Chen_2016_deoxypodophyllotoxin"
  units <- list(
    time = "min",
    dosing = "mg",
    concentration = "mg/L", # numerically identical to the paper's ug/mL
    weight = "kg"
  )

  # Issue #482: what each ODE state holds.
  compartmentData <- list(
    # The peripheral compartment is internal to the linCmt() analytic
    # solution and is not an addressable ODE state, so it has no entry here.
    central = list(analyte = "deoxypodophyllotoxin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 0L,
    n_studies = 1L,
    age_range = NULL,
    weight_range = "70 kg (the reference body weight of the reverse transformation)",
    sex_female_pct = NA_real_,
    disease_state = "None. No clinical data for deoxypodophyllotoxin existed when the paper was written; this is a pure forward prediction.",
    dose_range = "Simulated single intravenous bolus of 16 mg (not mg/kg) in a 70 kg adult.",
    regions = NULL,
    notes = paste0(
      "The underlying animal data are the mouse (12.5 and 25.0 mg/kg), rat ",
      "(1.0, 2.0 and 4.0 mg/kg) and monkey (0.5, 1.0 and 2.0 mg/kg) ",
      "intravenous dose groups; the 0.3 mg/kg beagle dog group was excluded. ",
      "The 16 mg simulation dose was chosen from the monkey maximum ",
      "tolerated dose of 4 mg/kg."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # The paper reports this model as a biexponential plasma profile for a
    # 16 mg intravenous bolus in a 70 kg human (Results, "Interspecies
    # Allometric Scaling"):
    #
    #   C = 0.628 * exp(-0.147 * t) + 0.066 * exp(-0.011 * t)
    #   (units: C in ug/mL, t in min)
    #
    # The four macro-constants are converted to the canonical clearance /
    # volume parameterisation by the standard identities, with A = 0.628,
    # alpha = 0.147, B = 0.066, beta = 0.011 and D = 16 mg:
    #
    #   vc  = D / (A + B)                      = 16 / 0.694      = 23.055 L
    #   k21 = (A*beta + B*alpha) / (A + B)                       = 0.023934 1/min
    #   k10 = alpha * beta / k21                                 = 0.067561 1/min
    #   k12 = alpha + beta - k21 - k10                           = 0.066505 1/min
    #   cl  = k10 * vc                                           = 1.5576 L/min
    #   q   = k12 * vc                                           = 1.5333 L/min
    #   vp  = q / k21                                            = 64.057 L
    #
    # The conversion reproduces the three derived quantities the paper
    # prints for this model: CL = 1.56 L/min (here 1.5576), Vss = vc + vp =
    # 87.11 L (paper 87.67), terminal t1/2 = ln(2)/beta = 63.0 min (paper
    # 63.08) and AUC0-inf = D / cl = 10.27 ug*min/mL (paper 10.25). The
    # residual differences are the rounding of the printed macro-constants
    # to three decimal places.
    #
    # The paper also prints the two allometric regressions the Dedrick plot
    # was built from, across mouse, rat and monkey (W in kg):
    #   CL  = 44.28 * W^0.832  (mL/min)  -> 1.52 L/min at 70 kg
    #   Vss =  1.99 * W^0.857  (L)       -> 75.65 L    at 70 kg
    # Those give a slightly different human clearance and volume than the
    # Dedrick reverse transformation encoded here; the paper reports both
    # and treats the Dedrick profile as the predicted concentration-time
    # curve. Everything is FIXED: nothing in this model was estimated on
    # human data.
    # ------------------------------------------------------------------
    lcl <- fixed(log(1.5576))
    label("Clearance (L/min)")
    lvc <- fixed(log(23.055))
    label("Central volume of distribution (L)")
    lq <- fixed(log(1.5333))
    label("Intercompartmental clearance (L/min)")
    lvp <- fixed(log(64.057))
    label("Peripheral volume of distribution (L)")

    # The paper reports no residual-error magnitude for this model (it is a
    # deterministic reverse transformation of a pooled animal fit, not a
    # population fit to human data). Encoded as zero rather than invented;
    # see the vignette Errata.
    propSd <- fixed(0)
    label("Proportional residual error (fraction; ZERO - magnitude not reported in source)")
  })

  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)

    Cc <- linCmt()
    Cc ~ prop(propSd)
  })
}
