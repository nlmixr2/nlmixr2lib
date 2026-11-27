Ji_2020_benapenem <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenously infused",
    "benapenem (a carbapenem) in healthy adult volunteers (Ji 2020).",
    "Clearance and central volume both scale as power functions of the",
    "administered dose, (DOSE / 250 mg)^0.183 and (DOSE / 250 mg)^0.184,",
    "and central volume scales linearly with body weight (WT / 70 kg).",
    "The infusion duration was estimated (D1 = 0.404 h against a nominal",
    "0.5 h). The unbound plasma concentration is derived from the total",
    "with the paper's two-site (saturable plus linear) human plasma",
    "protein-binding model (Bmax 115 ug/mL, C50 8.63 ug/mL, slope 0.666),",
    "which the authors used to compute %fT>MIC for probability of target",
    "attainment.",
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
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "benapenem (total, bound + unbound)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "benapenem", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    DOSE_BENAPENEM_MG = list(
      description = "Benapenem dose per administration",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters clearance as (DOSE / 250)^0.183 and central volume as",
        "(DOSE / 250)^0.184 (Equations 3 and 4; 250 mg is the lowest",
        "phase I dose). The phase I study gave single 250, 500 or 1,000 mg",
        "doses as 30-min IV infusions; the authors' Monte Carlo",
        "simulations used the same three doses every 6, 8, 12 or 24 h, so",
        "the covariate is the per-administration dose, held constant over a",
        "regimen."
      ),
      source_name = "dose"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Linear (exponent 1) on central volume, (WT / 70) (Equation 4).",
        "Phase I cohort median 61.5 kg, range 51.8-78.3 kg (Table 1)."
      ),
      source_name = "WT"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 12L,
    n_studies = 1L,
    age_range = "21-35 years (median 28)",
    weight_range = "51.8-78.3 kg (median 61.5)",
    height_range = "153.5-181 cm (median 165.3)",
    sex_female_pct = 50,
    disease_state = "Healthy volunteers (benapenem phase I trial).",
    dose_range = paste(
      "Single 250, 500 or 1,000 mg IV doses infused over 30 min. Sampling",
      "at 0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5, 8, 12, 24, 36 and 48 h."
    ),
    regions = "China.",
    notes = paste(
      "Demographics from Table 1 (12 subjects, 6 female / 6 male). The paper",
      "does not say how the 12 subjects were split across the three dose",
      "groups. Data were supplied by XuanZhu Pharma Co., Ltd."
    )
  )

  ini({
    # Structural PK, Table 4 (NONMEM FOCE-I). Equations 3 and 4.
    lcl <- log(0.825); label("Clearance at 250 mg (L/h)") # Table 4: CL = 0.825 L/h (RSE 3%)
    lvc <- log(5.03); label("Central volume at 250 mg and 70 kg (L)") # Table 4: V1 = 5.03 L (RSE 2%)
    lq <- log(1.14); label("Intercompartmental clearance (L/h)") # Table 4: Q = 1.14 L/h (RSE 8%)
    lvp <- log(2.93); label("Peripheral volume (L)") # Table 4: V2 = 2.93 L (RSE 4%)
    ld1 <- log(0.404); label("Estimated infusion duration D1 (h)") # Table 4: D1 = 0.404 h (RSE 1%)

    e_dose_cl <- 0.183; label("Power exponent of dose / 250 mg on CL (unitless)") # Table 4: Dose_CL = 0.183 (RSE 8%); Equation 3
    e_dose_vc <- 0.184; label("Power exponent of dose / 250 mg on V1 (unitless)") # Table 4: Dose_V1 = 0.184 (RSE 7%); Equation 4

    # Two-site plasma protein-binding model (Equation 9), human column of
    # Table 5. Fitted to in vitro binding data before the PK model was
    # developed and not re-estimated with it, so held fixed here.
    lbmax <- fixed(log(115)); label("Maximum binding capacity of the high-affinity site Bmax (ug/mL)") # Table 5 human: Bmax = 115 ug/mL (RSE 15%)
    lkd <- fixed(log(8.63)); label("Apparent dissociation constant of the high-affinity site C50 (ug/mL)") # Table 5 human: C50 = 8.63 ug/mL (RSE 16%)
    kns <- fixed(0.666); label("Linear low-affinity binding slope (unitless)") # Table 5 human: Slope = 0.666 (RSE 15%)

    # IIV: exponential (Equation 5). Table 4 prints IIV_CL 10.5 and IIV_V1
    # 7.6 without a definition; read as CV percentages and converted with
    # omega^2 = log(1 + CV^2).
    etalcl ~ log(1 + 0.105^2) # Table 4: IIV_CL 10.5 (RSE 15)
    etalvc ~ log(1 + 0.076^2) # Table 4: IIV_V1 7.6 (RSE 16)

    # Residual error: Equation 6, Cobs = Cpred * (1 + eps1) + eps2.
    propSd <- 0.074; label("Proportional residual error (fraction)") # Table 4: Prop. error = 7.4% (RSE 4%)
    addSd <- 0.233; label("Additive residual error (mg/L)") # Table 4: Add. error = 0.233 (RSE 15%)
  })
  model({
    cl <- exp(lcl + etalcl) * (DOSE_BENAPENEM_MG / 250)^e_dose_cl
    vc <- exp(lvc + etalvc) * (DOSE_BENAPENEM_MG / 250)^e_dose_vc * (WT / 70)
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Infusion duration applies to dose records with rate = -2; records that
    # carry an explicit rate or duration use that instead.
    d1 <- exp(ld1)
    dur(central) <- d1

    # Total plasma concentration (mg/L = ug/mL).
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

    Cc ~ add(addSd) + prop(propSd)
  })
}
