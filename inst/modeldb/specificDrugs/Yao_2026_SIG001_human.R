Yao_2026_SIG001_human <- function() {
  description <- "Human first-in-human PROJECTION (no human data were fitted). Two-compartment PK model with linear elimination for SIG001, a humanized IgG monoclonal antibody against cancer-derived sialylated IgG (SIA-IgG), obtained by single-species allometric scaling of the cynomolgus monkey fit (Yao_2026_SIG001_cyno) from a 2.24-kg monkey to a 70-kg human with exponents 0.8 (CL), 0.75 (Q), 1.0 (VC) and 0.95 (VP) (Yao 2026). Inter-individual variability (30% on CL) and residual error (10% proportional) are the authors' assumed simulation values, not estimates. Parameters are body-weight normalised (L/kg, L/h/kg) for a 70-kg adult: dose in mg/kg and simulate a nominal 1-kg subject. Target-mediated clearance in patients with SIA-IgG-expressing tumours is not represented."
  reference <- "Yao Q, Liu C, Zhang S, Liu Y, Zhou Y, Qiu X. Preclinical pharmacokinetic modeling and clinical translation of SIG001, a novel anti-cancer therapeutic antibody targeting on cancer-derived sialylated IgG. Front Pharmacol. 2026;17:1904283. doi:10.3389/fphar.2026.1904283"
  vignette <- "Yao_2026_SIG001"
  units <- list(time = "h", dosing = "mg/kg", concentration = "ug/mL")

  # States hold a body-weight-normalised amount (mg per kg of body weight).
  compartmentData <- list(
    central = list(analyte = "SIG001", units = "mg/kg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "SIG001", units = "mg/kg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = 0L,
    n_studies = 0L,
    weight_range = "70 kg (assumed human body weight in the allometric scaling, Equation 1)",
    disease_state = "Projected adult patients with epithelial-derived cancers (simulation only; no human subjects were studied).",
    dose_range = "Simulated 1.0 and 1.5 mg/kg IV every 2 weeks for the first-in-human safety projection (Figure 5); 0.15 mg/kg Q2W proposed starting dose; 1.8-7.0 mg/kg Q2W projected pharmacologically active dose range.",
    regions = "Not applicable (in silico).",
    notes = "IMPORTANT: this file contains no fitted human parameters. Yao 2026 Equation 1: P_human (absolute) = P_monkey (absolute) x (70 / 2.24)^alpha, with alpha = 0.8 (CL), 0.75 (Q), 1.0 (VC), 0.95 (VP); 2.24 kg is the mean body weight of the monkeys in the PK study. The per-kg human values printed in Results are used here. IIV on CL of 30% and 10% proportional residual error are assumed values the authors took from the monoclonal-antibody literature (Discussion, citing Dirks and Meibohm 2010). The authors caution that higher SIA-IgG expression in patients could add target-mediated clearance not captured by this linear model."
  )

  ini({
    # Yao 2026 Results 'Projection of FIH dose of SIG001' prints the scaled
    # human values in mL/h/kg and mL/kg; converted to L/h/kg and L/kg.
    # All values are allometric projections, so every parameter is fixed.
    # Equation 1 reproduces them from Table 1: 0.386 x (70/2.24)^-0.2 = 0.194;
    # 47.2 x (70/2.24)^0 = 47.2; 1.05 x (70/2.24)^-0.25 = 0.444;
    # 31.4 x (70/2.24)^-0.05 = 26.4.
    lcl <- fixed(log(0.000194))
    label("Clearance CL (L/h/kg)") # Results: CL = 0.194 mL/h/kg
    lvc <- fixed(log(0.0472))
    label("Central volume of distribution VC (L/kg)") # Results: VC = 47.2 mL/kg
    lq <- fixed(log(0.00044))
    label("Intercompartmental clearance Q (L/h/kg)") # Results: Q = 0.44 mL/h/kg
    lvp <- fixed(log(0.0264))
    label("Peripheral volume of distribution VP (L/kg)") # Results: VP = 26.4 mL/kg

    # Assumed (not estimated) simulation variability, Results and Discussion:
    # '30% IIV of CL and a 10% proportional residual variability'.
    # 30% treated as a log-normal CV: omega^2 = log(0.3^2 + 1).
    etalcl ~ fixed(0.08618)

    propSd <- fixed(0.10)
    label("Proportional residual error (fraction)") # Results: 10% proportional residual variability (assumed)
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
