Yao_2026_SIG001_cyno <- function() {
  description <- "Preclinical (cynomolgus monkey). Two-compartment population PK model with linear elimination for SIG001, a humanized IgG monoclonal antibody against cancer-derived sialylated IgG (SIA-IgG), after single IV doses of 5, 15 and 30 mg/kg in healthy cynomolgus monkeys (Yao 2026). Parameters are body-weight normalised (L/kg, L/h/kg): dose in mg/kg and simulate a nominal 1-kg subject. This fit is the source of the allometric human projection Yao_2026_SIG001_human."
  reference <- "Yao Q, Liu C, Zhang S, Liu Y, Zhou Y, Qiu X. Preclinical pharmacokinetic modeling and clinical translation of SIG001, a novel anti-cancer therapeutic antibody targeting on cancer-derived sialylated IgG. Front Pharmacol. 2026;17:1904283. doi:10.3389/fphar.2026.1904283"
  vignette <- "Yao_2026_SIG001"
  units <- list(time = "h", dosing = "mg/kg", concentration = "ug/mL")

  # States hold a body-weight-normalised amount (mg per kg of body weight):
  # the model was fitted to mg/kg doses with mL/kg volumes (Yao 2026 Table 1).
  compartmentData <- list(
    central = list(analyte = "SIG001", units = "mg/kg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "SIG001", units = "mg/kg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "cynomolgus monkey",
    n_subjects = 18L,
    n_studies = 1L,
    weight_range = "mean 2.24 kg (Methods, 'Allometric scaling of PK parameters in humans')",
    sex_female_pct = 50,
    disease_state = "Healthy (non-tumour-bearing) animals.",
    dose_range = "Single IV dose of 5, 15 or 30 mg/kg; 6 animals (3 male, 3 female) per dose group.",
    regions = "WuXi AppTec (Suzhou, China).",
    notes = "Yao 2026 Methods 'Preclinical PK analysis': serial sampling pre-dose and at 0.5, 1, 2, 4, 8, 24 h and 2, 3, 4, 7, 10, 14, 21, 28 days; ELISA with LLOQ 100 ng/mL. 251 quantifiable observations plus 1 BQL observation (Results). NONMEM 7.5.0 FOCE-I. No sex difference in PK. No covariates in the final model."
  )

  ini({
    # Yao 2026 Table 1 (cynomolgus monkey), reported in mL/h/kg and mL/kg;
    # converted to L/h/kg and L/kg here so that mg/kg / (L/kg) = mg/L = ug/mL.
    lcl <- log(0.000386)
    label("Clearance CL (L/h/kg)") # Table 1 CL = 0.386 mL/h/kg (RSE 4.7%)
    lvc <- log(0.0472)
    label("Central volume of distribution VC (L/kg)") # Table 1 VC = 47.2 mL/kg (RSE 3.3%)
    lq <- log(0.00105)
    label("Intercompartmental clearance Q (L/h/kg)") # Table 1 Q = 1.05 mL/h/kg (RSE 17.0%)
    lvp <- log(0.0314)
    label("Peripheral volume of distribution VP (L/kg)") # Table 1 VP = 31.4 mL/kg (RSE 7.9%)

    # Table 1 reports IIV as a percentage of an exponential (log-normal)
    # random effect; converted as omega^2 = log(CV^2 + 1).
    etalcl ~ 0.03474 # Table 1 IIV CL = 18.8%: log(0.188^2 + 1)
    etalvc ~ 0.02590 # Table 1 IIV VC = 16.2%: log(0.162^2 + 1)

    propSd <- 0.126
    label("Proportional residual error (fraction)") # Table 1 proportional error = 12.6%
  })
  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
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
