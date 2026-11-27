Ji_2020_benapenem_mouse <- function() {
  description <- paste(
    "Preclinical (mouse). Two-compartment population PK model for",
    "benapenem (a carbapenem) given as an IV bolus to neutropenic",
    "thigh-infected ICR mice (Ji 2020). Clearance and central volume both",
    "rise with the mg/kg dose as CL = 0.913 * [1 + (DOSE / 58.4)^1.46]",
    "mL/h and V1 = 3.62 * [1 + (DOSE / 58.4)^1.06] mL. The unbound plasma",
    "concentration is derived from the total with the paper's two-site",
    "(saturable plus linear) mouse plasma protein-binding model (Bmax 164",
    "ug/mL, C50 0.232 ug/mL, slope 1.02), which the authors used to compute",
    "%fT>MIC for the mouse PK/PD analysis.",
    sep = " "
  )
  reference <- paste(
    "Ji XW, Xue F, Kang ZS, Zhong W, Kuan IHS, Yang XP, Zhu X, Li Y,",
    "Lv Y (2020). Model-Informed Drug Development, Pharmacokinetic/",
    "Pharmacodynamic Cutoff Value Determination, and Antibacterial",
    "Efficacy of Benapenem against Enterobacteriaceae. Antimicrobial",
    "Agents and Chemotherapy 64(3):e01751-19. doi:10.1128/AAC.01751-19.",
    sep = " "
  )
  vignette <- "Ji_2020_benapenem"
  units <- list(time = "h", dosing = "ug", concentration = "ug/mL")

  compartmentData <- list(
    central = list(analyte = "benapenem (total, bound + unbound)", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "benapenem", units = "ug", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    DOSE_BENAPENEM_MGKG = list(
      description = "Benapenem dose per administration normalised to body weight",
      units = "mg/kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters clearance as [1 + (DOSE / 58.4)^1.46] and central volume as",
        "[1 + (DOSE / 58.4)^1.06] (Equations 1 and 2); 58.4 mg/kg is the",
        "third of the four PK-study dose levels (1.9, 14.6, 58.4 and 233",
        "mg/kg). The dose records themselves carry the amount per mouse in",
        "ug (mg/kg x body weight in kg x 1000); the paper does not report",
        "mouse body weights."
      ),
      source_name = "dose"
    )
  )

  population <- list(
    species = "mouse (ICR, neutropenic, thigh-infected)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    age_range = "5-6 weeks",
    sex_female_pct = NA_real_,
    disease_state = paste(
      "Cyclophosphamide-induced neutropenic ICR mice (150 then 100 mg/kg",
      "intraperitoneally, 3 days apart) infected intramuscularly in the",
      "thigh with Escherichia coli ATCC 25922 (ESBL-negative) 2 h before",
      "dosing."
    ),
    dose_range = paste(
      "Single IV (tail vein) doses of 1.9, 14.6, 58.4 and 233 mg/kg.",
      "Plasma sampled at 0, 0.083, 0.25, 0.5, 1, 2, 4, 8, 12 and 24 h,",
      "six mice per time point (destructive sampling)."
    ),
    regions = "China (Peking University First Hospital).",
    notes = paste(
      "Males and females; numbers by sex and body weights are not",
      "reported. Six mice per time point across 10 time points and 4 dose",
      "levels implies about 240 mice if each mouse gave one sample; the",
      "paper does not state the subject count."
    )
  )

  ini({
    # Structural PK, Table 3 (NONMEM FOCE-I). Equations 1 and 2.
    lcl <- log(0.913); label("Clearance coefficient CL_pop in CL = CL_pop * [1 + (DOSE / 58.4)^e] (mL/h)") # Table 3: CL = 0.913 mL/h (RSE 6%)
    lvc <- log(3.62); label("Central volume coefficient V1_pop in V1 = V1_pop * [1 + (DOSE / 58.4)^e] (mL)") # Table 3: V1 = 3.62 mL (RSE 3%)
    lq <- log(0.09); label("Intercompartmental clearance (mL/h)") # Table 3: Q = 0.09 mL/h (RSE 36%)
    lvp <- log(1.65); label("Peripheral volume (mL)") # Table 3: V2 = 1.65 mL (RSE 34%)

    e_dose_cl <- 1.46; label("Exponent of dose / 58.4 mg/kg in the CL dose term (unitless)") # Table 3: Dose_CL = 1.46 (RSE 4%); Equation 1
    e_dose_vc <- 1.06; label("Exponent of dose / 58.4 mg/kg in the V1 dose term (unitless)") # Table 3: Dose_V1 = 1.06 (RSE 5%); Equation 2

    # Two-site plasma protein-binding model (Equation 9), mouse column of
    # Table 5. Fitted to in vitro binding data before the PK model was
    # developed and not re-estimated with it, so held fixed here.
    lbmax <- fixed(log(164)); label("Maximum binding capacity of the high-affinity site Bmax (ug/mL)") # Table 5 mouse: Bmax = 164 ug/mL (RSE 3%)
    lkd <- fixed(log(0.232)); label("Apparent dissociation constant of the high-affinity site C50 (ug/mL)") # Table 5 mouse: C50 = 0.232 ug/mL (RSE 6%)
    kns <- fixed(1.02); label("Linear low-affinity binding slope (unitless)") # Table 5 mouse: Slope = 1.02 (RSE 5%)

    # IIV: exponential (Equation 5), on CL only. Table 3 prints IIV_CL 36.5
    # without a definition; read as a CV percentage and converted with
    # omega^2 = log(1 + CV^2).
    etalcl ~ log(1 + 0.365^2) # Table 3: IIV_CL 36.5 (RSE 15), shrinkage 11

    # Residual error: Equation 6 with the additive component fixed to zero
    # (Table 3 'Add. error 0 (fix)'), leaving a proportional error only.
    propSd <- 0.301; label("Proportional residual error (fraction)") # Table 3: Prop. error = 30.1% (RSE 7%)
  })
  model({
    cl <- exp(lcl + etalcl) * (1 + (DOSE_BENAPENEM_MGKG / 58.4)^e_dose_cl)
    vc <- exp(lvc) * (1 + (DOSE_BENAPENEM_MGKG / 58.4)^e_dose_vc)
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Total plasma concentration (ug/mL = mg/L).
    Cc <- central / vc

    # Unbound concentration (Equations 9-12). Total = unbound + bound with
    # bound = bmax * Cu / (kd + Cu) + kns * Cu, so a * Cu^2 + b * Cu + c = 0
    # with a = 1 + kns, b = (1 + kns) * kd + bmax - Cc and c = -kd * Cc. The
    # positive root (Equation 12) is written in the rationalised form, which
    # stays finite and non-negative for every Cc >= 0.
    bmax <- exp(lbmax)
    kd <- exp(lkd)
    bterm <- (1 + kns) * kd + bmax - Cc
    Cunbound <- 2 * kd * Cc / (bterm + sqrt(bterm * bterm + 4 * (1 + kns) * kd * Cc))

    Cc ~ prop(propSd)
  })
}
