Pereira_2024_M6495_human <- function() {
  description <- "Human allometric projection (no human data in the fit). Cynomolgus-monkey-to-human scaled PK/PD model for the anti-ADAMTS-5 NANOBODY M6495 (Pereira 2024): two-compartment PK with first-order SC absorption and parallel linear plus Michaelis-Menten elimination, linked to an indirect-response model in which M6495 inhibits the production of the serum aggrecan neo-epitope ARGS. V1 and V2 scale with body weight to the power 1, CL, Q and Vmax to 0.75 and Ka to -0.25; Km, F, kout, Imax and IC50 are the monkey values, and the ARGS baseline is the measured healthy-human value. Used to choose the doses of the first-in-human single ascending dose study."
  reference <- "Pereira JNS, Ottevaere I, Serruys B, Guehring H, Ladel C, Lindemann S. Translational pharmacokinetic and pharmacodynamic modelling of the anti-ADAMTS-5 NANOBODY(R) (M6495) using the neo-epitope ARGS as a biomarker. J Pharmacokinet Pharmacodyn. 2024;52(1):8. doi:10.1007/s10928-024-09958-z. PMCID: PMC11662058. Scaled values from Table 2 (70 kg human column) and Supplementary Table 3; interindividual variances from the human simulation R code in Supplementary Material 1."
  vignette <- "Pereira_2024_M6495"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # ARGS is the serum aggrecan neo-epitope released by ADAMTS-5 cleavage; it
  # is a paper-specific biomarker turnover state with no canonical name.
  paper_specific_compartments <- c("args")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling covariate with reference 70 kg (Table 2 'Parameters for a 70 Kg Human'). Exponents are fixed by the authors' scaling rule, Eq. 11 and Supplementary Table 3: 1 for V1 and V2, 0.75 for CL, Q and Vmax, -0.25 for Ka. The paper applied the rule to a single 70 kg human; the supplementary simulation code hardcodes the 70 kg values, so simulate with WT = 70 to reproduce Figures 4 and 5.",
      source_name = "BW_human (Eq. 11); not a fitted covariate"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "M6495", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "M6495", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "M6495", units = "mg", specimen = "serum", verified = TRUE),
    args = list(analyte = "aggrecan ARGS neo-epitope", units = "nM", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human (allometric projection from cynomolgus monkey; no human PK or PD data in the fit)",
    n_subjects = NA_integer_,
    n_studies = 0L,
    age_range = "Not applicable; no human subjects contributed to the model fit.",
    weight_range = "Not applicable; projection anchored at a 70 kg reference adult.",
    sex_female_pct = NA_real_,
    disease_state = "Healthy adults (projection target). The ARGS baseline and its variability come from serum samples of 48 healthy humans measured in-house.",
    dose_range = "Projected SC single doses 0.2 to 300 mg; the first-in-human single ascending dose study (NCT03224702) used 1, 5, 20, 75, 150 and 300 mg.",
    notes = "Source fit: Pereira_2024_M6495_cyno (70 cynomolgus monkeys, 0.01-150 mg/kg). Table 2 converts the per-kg monkey estimates to a 3.9 kg monkey before scaling to 70 kg. The authors overlaid the projection on the published Phase I single ascending dose data (Fig. 4 and 5)."
  )

  ini({
    lka <- log(0.0191); label("First-order SC absorption rate constant, Ka, for a 70 kg adult (1/h)")      # Table 2 human column: Ka 0.0191 /h (monkey 0.0394 x (70/3.9)^-0.25)
    lfdepot <- log(1); label("SC bioavailability, F (fraction)")                                              # Table 2 human column: F 1 (assumed equal to monkey)
    lcl <- log(0.0085); label("Linear clearance, CL, for a 70 kg adult (L/h)")                              # Table 2 human column: CL 0.0085 L/h
    lvc <- log(3.22); label("Central volume, V1, for a 70 kg adult (L)")                                    # Table 2 human column: V1 3.22 L
    lq <- log(0.0558); label("Intercompartmental clearance, Q, for a 70 kg adult (L/h)")                    # Table 2 human column: Q 0.0558 L/h
    lvp <- log(2.177); label("Peripheral volume, V2, for a 70 kg adult (L)")                                # Table 2 human column: V2 2.177 L
    lvmax <- log(0.0142); label("Maximum rate of non-linear elimination, Vmax, for a 70 kg adult (mg/h)")    # Table 2 human column: Vmax 0.0142 mg/h
    lkm <- log(0.265); label("Concentration at half-maximal non-linear elimination, Km (mg/L)")             # Table 2 human column: Km 0.265 mg/L (assumed equal to monkey)

    e_wt_ka <- fixed(-0.25); label("Allometric exponent on Ka for WT/70 (unitless)")                         # Table 2 and Supplementary Table 3: Ka exponent -0.25
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL for WT/70 (unitless)")                          # Table 2 and Supplementary Table 3: CL exponent 0.75
    e_wt_vc <- fixed(1); label("Allometric exponent on V1 for WT/70 (unitless)")                             # Table 2 and Supplementary Table 3: V1 exponent 1
    e_wt_q <- fixed(0.75); label("Allometric exponent on Q for WT/70 (unitless)")                            # Table 2 and Supplementary Table 3: Q exponent 0.75
    e_wt_vp <- fixed(1); label("Allometric exponent on V2 for WT/70 (unitless)")                             # Table 2 and Supplementary Table 3: V2 exponent 1
    e_wt_vmax <- fixed(0.75); label("Allometric exponent on Vmax for WT/70 (unitless)")                      # Table 2 and Supplementary Table 3: Vmax exponent 0.75

    lkout <- log(0.0401); label("First-order ARGS elimination rate constant, kout (1/h)")                   # Table 2 human column: Kout 0.0401 /h (assumed equal to monkey)
    limax <- log(0.725); label("Maximum fractional inhibition of ARGS production, Imax (fraction)")         # Table 2 human column: IMax 0.725 (assumed equal to monkey)
    lic50 <- log(1.17); label("M6495 concentration at half-maximal inhibition, IC50 (mg/L)")                # Table 2 human column: IC50 1.17 mg/L (assumed equal to monkey)
    lrbase <- log(0.14); label("Baseline serum ARGS concentration in healthy humans, BASE (nM)")            # Table 2 human column: BASE 0.14 nM (measured in 48 healthy humans)
    boxcox_limax <- -3.16; label("Box-Cox shape parameter for the interindividual distribution of Imax (unitless)") # Table 1 Boxcox parameter -3.16; Supplementary R code BXPAR = -3.16

    # Variances as coded in the supplementary R simulation, sigma = diag(c(0.1,
    # 0.1, 0.055, 0.21, 0.039)). The Methods give the CL and Ka variance of 0.1
    # and the 46% CV of the human ARGS baseline (0.46^2 = 0.21); Vmax and Imax
    # are stated to keep the monkey variability.
    etalka ~ 0.1 # Methods and Supplementary R code: Ka variance 0.1 (inflated for human PK variability)
    etalcl ~ 0.1 # Methods and Supplementary R code: CL variance 0.1 (inflated for human PK variability)
    etalvmax ~ 0.055 # Supplementary R code: Vmax variance 0.055
    etalimax ~ 0.039 # Supplementary R code: Imax variance 0.039 (Box-Cox transformed in model())
    etalrbase ~ 0.21 # Supplementary R code: BASE variance 0.21 (Methods: 46% CV in 48 healthy humans)

    # Not part of the human projection: the published simulations carry no
    # residual error. Carried from the monkey fit (Table 1 RE1 and RE2) so that
    # simulated observations can be generated.
    propSd <- 0.152; label("Proportional residual error for M6495 (fraction)")                         # Table 1: RE1 (PK) 15.2%, carried from the monkey fit
    propSd_ARGS <- 0.221; label("Proportional residual error for ARGS (fraction)")                     # Table 1: RE2 (PD) 22.1%, carried from the monkey fit
  })

  model({
    ka <- exp(lka + etalka) * (WT / 70)^e_wt_ka
    fdepot <- exp(lfdepot)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc) * (WT / 70)^e_wt_vc
    q <- exp(lq) * (WT / 70)^e_wt_q
    vp <- exp(lvp) * (WT / 70)^e_wt_vp
    vmax <- exp(lvmax + etalvmax) * (WT / 70)^e_wt_vmax
    km <- exp(lkm)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Imax: Box-Cox transform of its eta (Eq. 8-10) keeps Imax below 1.
    phi_limax <- (exp(etalimax)^boxcox_limax - 1) / boxcox_limax
    imax <- exp(limax + phi_limax)
    ic50 <- exp(lic50)
    kout <- exp(lkout)
    rbase <- exp(lrbase + etalrbase)
    kin <- rbase * kout # Eq. 6, steady state at time 0

    Cc <- central / vc

    # Eq. 2-5, as coded in the supplementary R simulation (bioavailability in
    # the absorption flux; identical to scaling the SC dose by F).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - k12 * central + k21 * peripheral1 - kel * central - vmax * Cc / (km + Cc)
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(args) <- kin * (1 - imax * Cc / (ic50 + Cc)) - kout * args
    args(0) <- rbase # Eq. 7
    f(depot) <- fdepot

    ARGS <- args

    Cc ~ prop(propSd)
    ARGS ~ prop(propSd_ARGS)
  })
}
