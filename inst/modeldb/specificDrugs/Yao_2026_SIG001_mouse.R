Yao_2026_SIG001_mouse <- function() {
  description <- "Preclinical (mouse, CD-1). Two-compartment population PK model with linear elimination for SIG001, a humanized IgG monoclonal antibody against cancer-derived sialylated IgG (SIA-IgG), after single IV doses of 5, 10 and 20 mg/kg in non-tumour-bearing CD-1 mice (Yao 2026). Fitted to sparse (two samples per animal) data with the M3 method for BQL observations. Parameters are body-weight normalised (L/kg, L/h/kg): dose in mg/kg and simulate a nominal 1-kg subject. Used by the authors to compute the steady-state exposure of the efficacious xenograft regimens for the human pharmacologically active dose projection."
  reference <- "Yao Q, Liu C, Zhang S, Liu Y, Zhou Y, Qiu X. Preclinical pharmacokinetic modeling and clinical translation of SIG001, a novel anti-cancer therapeutic antibody targeting on cancer-derived sialylated IgG. Front Pharmacol. 2026;17:1904283. doi:10.3389/fphar.2026.1904283"
  vignette <- "Yao_2026_SIG001"
  units <- list(time = "h", dosing = "mg/kg", concentration = "ug/mL")

  # States hold a body-weight-normalised amount (mg per kg of body weight):
  # the model was fitted to mg/kg doses with mL/kg volumes (Yao 2026 Table 2).
  compartmentData <- list(
    central = list(analyte = "SIG001", units = "mg/kg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "SIG001", units = "mg/kg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "mouse (CD-1)",
    n_subjects = 135L,
    n_studies = 1L,
    sex_female_pct = 50,
    disease_state = "Healthy (non-tumour-bearing) CD-1 mice.",
    dose_range = "Single IV dose of 5, 10 or 20 mg/kg.",
    regions = "WuXi AppTec (Suzhou, China).",
    notes = "Yao 2026 Methods 'Preclinical PK analysis': crossover (sparse) sampling, each animal sampled at two of the 14 monkey-study time points (0.5 h to 28 days); 6 animals per time point per dose group (12 at 28 days), equal numbers of males and females. n_subjects = 135 is derived (13 x 6 + 12 = 90 samples per dose group / 2 samples per animal x 3 dose groups), not printed. 216 quantifiable plus 54 BQL observations (Results); ELISA LLOQ 100 ng/mL; BQL handled by M3. NONMEM 7.5.0 FOCE-I. No covariates in the final model."
  )

  ini({
    # Yao 2026 Table 2 (mice), reported in mL/h/kg and mL/kg; converted to
    # L/h/kg and L/kg here so that mg/kg / (L/kg) = mg/L = ug/mL.
    lcl <- log(0.00137)
    label("Clearance CL (L/h/kg)") # Table 2 CL = 1.37 mL/h/kg (RSE 3.7%)
    lvc <- log(0.0577)
    label("Central volume of distribution VC (L/kg)") # Table 2 VC = 57.7 mL/kg (RSE 3.1%)
    lq <- log(0.00348)
    label("Intercompartmental clearance Q (L/h/kg)") # Table 2 Q = 3.48 mL/h/kg (RSE 15.8%)
    lvp <- log(0.0329)
    label("Peripheral volume of distribution VP (L/kg)") # Table 2 VP = 32.9 mL/kg (RSE 5.8%)

    # Only IIV on CL was estimable from two samples per animal (Discussion).
    # Table 2 reports it as a percentage of an exponential random effect;
    # converted as omega^2 = log(CV^2 + 1).
    etalcl ~ 0.17257 # Table 2 IIV CL = 43.4%: log(0.434^2 + 1)

    propSd <- 0.206
    label("Proportional residual error (fraction)") # Table 2 proportional error = 20.6%
  })
  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
