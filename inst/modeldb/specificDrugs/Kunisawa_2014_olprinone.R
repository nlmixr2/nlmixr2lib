Kunisawa_2014_olprinone <- function() {
  description <- "Two-compartment intravenous population PK model for olprinone (a phosphodiesterase III inhibitor) in healthy adult Japanese male volunteers with body-weight normalization on CL, Vc, Q and Vp (Kunisawa 2014)"
  reference <- paste(
    "Kunisawa T, Kasai H, Suda M, Yoshimura M, Sugawara A, Izumi Y,",
    "Iida T, Kurosawa A, Iwasaki H.",
    "Population pharmacokinetics of olprinone in healthy male volunteers.",
    "Clin Pharmacol Adv Appl. 2014;6:43-50.",
    "doi:10.2147/CPAA.S50626."
  )
  vignette <- "Kunisawa_2014_olprinone"
  units <- list(time = "hr", dosing = "ug", concentration = "ng/mL") # Methods: doses in ug/kg, plasma concentrations measured in ng/mL by HPLC; time converted from minutes-as-reported to hours

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear body-weight normalization per Methods: 'Dose and all pharmacokinetic parameters such as total clearance (CL), distribution volume of the central compartment (V1), intercompartmental clearance (Q ...) and distribution volume of the peripheral compartment (V2 ...) were adjusted for body weight.' Implemented as (WT/70)^1 with reference weight 70 kg; cohort mean weight was 64.1 kg.",
      source_name        = "WT"
    )
  )

  population <- list(
    n_subjects     = 9L,                                  # Table 1: nine unique healthy male volunteers (subject numbers 1-9); 39 cumulative subject-stages, 30 evaluable after exclusions
    n_studies      = 2L,                                  # Methods: Study I (single 5-min infusion, stages I-VII) and Study II (3-h infusion, stages VIII-X)
    age_range      = "24-32 years",                       # Table 1
    age_median     = "26 years",                          # Table 1 (median)
    weight_range   = "56.5-75.0 kg",                      # Table 1
    weight_median  = "61.5 kg",                           # Table 1 (median); mean 64.1 kg, SD 6.9
    sex_female_pct = 0,                                   # Title and Methods: 'healthy male volunteers' only
    race_ethnicity = c(Asian = 100),                      # Asahikawa Medical University Hospital, Japan; all Japanese (Methods + Acknowledgments)
    disease_state  = "Healthy adult male volunteers; no clinical abnormality on medical interview, physical examination, ECG, blood pressure, blood test, or urine test (Methods).",
    dose_range     = "1.25-50 ug/kg single 5-min IV infusion (Study I); 15-min IV loading (1.33-2.67 ug/kg/min) followed by 165-min continuous IV infusion (0.25-0.75 ug/kg/min) (Study II).",
    regions        = "Japan (Asahikawa Medical University, Hokkaido).",
    sampling       = "500 evaluable plasma concentrations (rich): Study I sampling at 5, 7.5, 10, 15, 30, 45 min and 1, 1.5, 2, 3, 4, 6, 8, 24 h post-dose; Study II at 15, 30 min and 1, 2, 3, 3.04, 3.08, 3.16, 3.25, 3.5, 4, 4.5, 5, 6, 7, 9, 11, 24 h post-dose.",
    notes          = "All subjects healthy adult Japanese males. Age was the only covariate tested; not retained (P=0.363, correlation coefficient with post hoc CL = 0.240). Quantification limit 0.1 ng/mL; concentrations <LLOQ excluded."
  )

  ini({
    # Structural PK parameters - reference weight 70 kg (Kunisawa 2014 Table 3 reports mL/min/kg; converted here to L/h at 70 kg)
    lcl <- log(30.954); label("Total clearance CL at 70 kg (L/h)")                            # Table 3: CL = 7.37 mL/min/kg = 7.37 * 60 / 1000 * 70 = 30.954 L/h at 70 kg (95% CI 6.63-8.11 mL/min/kg)
    lvc <- log(9.380);  label("Central volume of distribution Vc at 70 kg (L)")               # Table 3: V1 = 134 mL/kg = 134 * 70 / 1000 = 9.380 L at 70 kg (95% CI 115-153 mL/kg)
    lq  <- log(32.550); label("Intercompartmental clearance Q at 70 kg (L/h)")                # Table 3: Q = 7.75 mL/min/kg = 7.75 * 60 / 1000 * 70 = 32.550 L/h at 70 kg (95% CI 6.70-8.80 mL/min/kg)
    lvp <- log(19.250); label("Peripheral volume of distribution Vp at 70 kg (L)")            # Table 3: V2 = 275 mL/kg = 275 * 70 / 1000 = 19.250 L at 70 kg (95% CI 255-295 mL/kg)

    # Body-weight scaling exponents - fixed at 1 because the source paper reports parameters per-kg (linear scaling), not estimated allometric exponents
    e_wt_cl <- fixed(1); label("Body-weight scaling exponent on CL (linear normalization)")   # Methods: 'all pharmacokinetic parameters ... were adjusted for body weight' (linear, per-kg reporting in Table 3)
    e_wt_vc <- fixed(1); label("Body-weight scaling exponent on Vc (linear normalization)")   # Methods: linear body-weight normalization (Table 3 reports V1 in mL/kg)
    e_wt_q  <- fixed(1); label("Body-weight scaling exponent on Q (linear normalization)")    # Methods: linear body-weight normalization (Table 3 reports Q in mL/min/kg)
    e_wt_vp <- fixed(1); label("Body-weight scaling exponent on Vp (linear normalization)")   # Methods: linear body-weight normalization (Table 3 reports V2 in mL/kg)

    # IIV - only on CL per Table 3 base/final model selection (single eta retained)
    # Reported as variance omega^2 = 0.0153 on the internal (log) scale; back-transformed CV = 12.4%
    etalcl ~ 0.0153 # Table 3: omega_CL^2 = 0.0153 (SE 0.00557); 95% CI 0.00438-0.0262; reported %CV = 12.4

    # Residual error - mixed exponential + additive (NONMEM 'exponential' EPS approximates a proportional residual on the linear scale)
    # Table 3: sigma^2_exponential = 0.0495 -> SD = sqrt(0.0495) = 0.2225 (paper's '22.2%')
    # Table 3: sigma^2_additive    = 0.0166 -> SD = sqrt(0.0166) = 0.1288 ng/mL (paper's '0.129 SD')
    propSd <- 0.2225; label("Proportional residual error (fraction)")                          # Table 3: sigma^2_exponential = 0.0495; %CV reported as 22.2; mapped from NONMEM exponential to nlmixr2 proportional
    addSd  <- 0.1288; label("Additive residual error (ng/mL)")                                 # Table 3: sigma^2_additive = 0.0166; SD reported as 0.129
  })
  model({
    # Individual structural parameters with linear body-weight scaling to a 70 kg reference
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc)          * (WT / 70)^e_wt_vc
    q  <- exp(lq)           * (WT / 70)^e_wt_q
    vp <- exp(lvp)          * (WT / 70)^e_wt_vp

    Cc <- linCmt()
    Cc ~ prop(propSd) + add(addSd)
  })
}
