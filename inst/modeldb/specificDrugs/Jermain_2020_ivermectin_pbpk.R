Jermain_2020_ivermectin_pbpk <- function() {
  description <- paste(
    "PBPK (minimal, four-compartment; Phoenix WinNonlin 8.2). Oral",
    "ivermectin disposition in adults with a transit-delayed first-order",
    "absorption, a plasma compartment, a lung compartment perfused by the",
    "full cardiac output and a lumped rest-of-body (Other) compartment",
    "perfused by a fraction of the cardiac output (Jermain et al. 2020, J",
    "Pharm Sci 109:3574). Tissue uptake is perfusion-rate-limited. Plasma",
    "volume, lung volume, cardiac output and the remaining volume are fixed",
    "physiological values; the lung partition coefficient was fixed from",
    "the lung-to-plasma AUC ratio in calves given subcutaneous ivermectin;",
    "ka, CL/F, the Other partition coefficient, the fraction of cardiac",
    "output and the transit rate were fitted to a simulated geometric-mean",
    "plasma profile after a 15 mg oral dose. Clearance and the partition",
    "coefficients are apparent (divided by oral bioavailability F). The",
    "model carries no between-subject variability and was built to",
    "predict total ivermectin concentrations in lung tissue after single",
    "oral doses for COVID-19 repurposing."
  )
  reference <- paste(
    "Jermain B, Hanafin PO, Cao Y, Lifschitz A, Lanusse C, Rao GG.",
    "Development of a Minimal Physiologically-Based Pharmacokinetic Model",
    "to Simulate Lung Exposure in Humans Following Oral Administration of",
    "Ivermectin for COVID-19 Drug Repurposing. J Pharm Sci.",
    "2020;109(12):3574-3578. doi:10.1016/j.xphs.2020.08.024.",
    sep = " "
  )
  vignette <- "Jermain_2020_ivermectin_pbpk"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "ivermectin", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "ivermectin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ivermectin", units = "mg", specimen = "plasma", verified = TRUE),
    lung = list(analyte = "ivermectin", units = "mg", specimen = "tissue", verified = TRUE),
    other = list(analyte = "ivermectin", units = "mg", specimen = "tissue", verified = TRUE)
  )

  # No covariates: every physiological volume and flow is a fixed
  # population-average adult value and the fitted parameters describe one
  # geometric-mean profile.
  covariateData <- list()

  population <- list(
    species = "human (lung partition coefficient from Holstein calves)",
    n_subjects = 18L,
    n_studies = 1L,
    age_range = "adult",
    weight_range = "not reported (0.2 mg/kg dosing; average dose 15 mg)",
    disease_state = paste(
      "Healthy volunteers of the phase 1 study (El-Tahtawy et al., ref 13",
      "of the paper) whose population PK model was used to simulate the",
      "geometric-mean plasma profile the minimal PBPK model was fitted to."
    ),
    dose_range = paste(
      "Fitted to a single 15 mg oral dose (the study-average of 0.2 mg/kg);",
      "simulated at single oral doses of 12, 30 and 120 mg."
    ),
    notes = paste(
      "The human fit used no individual data: the plasma profile was a",
      "typical-value simulation (no IIV, no residual error) from a published",
      "population PK model. The lung partition coefficient (2.68) is the",
      "lung-to-plasma AUC ratio from calves dosed subcutaneously at",
      "0.2 mg/kg (Lifschitz et al., ref 12 of the paper), assumed",
      "species-independent."
    )
  )

  ini({
    # Fitted parameters: Jermain 2020 Table 1 (estimate with CV%).
    lka <- log(0.14)
    label("First-order absorption rate constant ka (1/h)")                        # Table 1 'Ka' = 0.14 (2.52%) 1/h
    lcl <- log(10.85)
    label("Apparent plasma clearance CLp/F (L/h)")                                # Table 1 'CLp/F' = 10.85 (3.78%) L/h
    lkp_other <- log(17.32)
    label("Apparent rest-of-body partition coefficient Kp,other/F (unitless)")   # Table 1 'Kp,other/F' = 17.32 (5.6%)
    fq_other <- 0.08
    label("Fraction of cardiac output perfusing the rest of the body (fraction)") # Table 1 'Fraction' = 0.08 (1.18%)
    lktr <- log(0.36)
    label("Transit rate constant ktr (1/h)")                                      # Table 1 'Ktr' = 0.36 (3.09%) 1/h

    # Fixed parameters: Table 1 footnote a, 'Parameter not estimated'.
    lkp_lung <- fixed(log(2.68))
    label("Apparent lung partition coefficient Kp,lung/F (unitless)")            # Table 1 'Kp,lung/F' = 2.68 (calf lung/plasma AUC ratio)
    lvc <- fixed(log(3))
    label("Plasma volume Vp (L)")                                                  # Table 1 'Vp' = 3 L (ref 15)
    lv_lung <- fixed(log(1.3))
    label("Lung volume VL (L)")                                                    # Table 1 'VL' = 1.3 L (ref 14)
    lv_other <- fixed(log(65.7))
    label("Remaining (rest-of-body) volume VOther (L)")                            # Table 1 'VOther' = 65.7 L
    q_co <- fixed(282)
    label("Cardiac output Qco (L/h)")                                              # Table 1 'Qco' = 282 L/h (ref 15)

    # The model was fitted to a noise-free typical-value profile, so no
    # IIV was estimated. The paper's lung-exposure simulations used 'an
    # assumed 20% RUV' (Methods, Human Lung Exposure Simulation), applied
    # here to both the plasma and the lung output.
    propSd <- fixed(0.2)
    label("Proportional residual error, plasma (fraction)")                       # Methods 'assumed 20% RUV'
    propSd_Clung <- fixed(0.2)
    label("Proportional residual error, lung (fraction)")                         # Methods 'assumed 20% RUV'
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    kp_other <- exp(lkp_other)
    kp_lung <- exp(lkp_lung)
    ktr <- exp(lktr)
    vc <- exp(lvc)
    v_lung <- exp(lv_lung)
    v_other <- exp(lv_other)

    # Page 3576 equations. The printed absorption and plasma equations read
    # 'Ka*Dose'; the Absorption state is the only quantity that can be
    # meant (Fig. 1 draws Absorption --Ka--> Plasma, and a constant
    # Ka*Dose input would never stop), so Ka multiplies the Absorption
    # amount here. The dose enters the transit compartment a0 (depot).
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ka * transit1
    d/dt(central) <- ka * transit1 +
      lung / v_lung * q_co / kp_lung -
      central / vc * cl -
      central / vc * q_co +
      other / v_other * q_co * fq_other / kp_other -
      central / vc * q_co * fq_other
    d/dt(lung) <- central / vc * q_co - lung / v_lung * q_co / kp_lung
    d/dt(other) <- central / vc * q_co * fq_other -
      other / v_other * q_co * fq_other / kp_other

    # mg/L -> ng/mL
    Cc <- 1000 * central / vc
    Clung <- 1000 * lung / v_lung

    Cc ~ prop(propSd)
    Clung ~ prop(propSd_Clung)
  })
}
