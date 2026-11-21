Carrot_2025_ovarianCancer_ca125 <- function() {
  description <- "K-PD model of CA-125 kinetics (KELIM-B model) during bevacizumab + olaparib + durvalumab in recurrent platinum-resistant or platinum-sensitive high-grade ovarian cancer (BOLD phase II GINECO trial). Arbitrary unit dose into a treatment compartment and a transit compartment; the transit amount inhibits zero-order CA-125 production through 1 - C2/(EC50 + C2); first-order CA-125 elimination KELIM. Baseline CA-125 and KELIM are estimated separately for platinum-resistant and platinum-sensitive relapse."
  reference <- "Carrot A, Corbaux P, Colomban O, Lebreton C, Gladieff L, Tredan O, Selle F, Abdeddaim C, Leary A, D'Hondt V, Dubot C, Reverdy T, You B, Freyer G. Prognostic Value of a Joint K-PD Model With Tumor Size Dynamics and CA-125 Kinetics in Recurrent Ovarian Cancer Patients: BOLD Phase II GINECO Study. CPT Pharmacometrics Syst Pharmacol. 2025;14(8):1310-1321. doi:10.1002/psp4.70025"
  vignette <- "Carrot_2025_ovarianCancer"
  paper_specific_compartments <- c("ca125")
  paper_specific_residual_sds <- c("expSd_ca125")
  units <- list(
    time = "day",
    dosing = "a.u. (arbitrary K-PD dose, fixed to 1 per treatment administration)",
    concentration = "kU/L (= U/mL) for CA-125"
  )

  compartmentData <- list(
    depot = list(
      analyte = "bevacizumab + olaparib + durvalumab (K-PD virtual treatment amount)",
      units = "a.u.",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit1 = list(
      analyte = "bevacizumab + olaparib + durvalumab (K-PD virtual treatment amount)",
      units = "a.u.",
      specimen = "not applicable",
      verified = TRUE
    ),
    ca125 = list(analyte = "CA-125", units = "kU/L", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    PLATINUM_RESISTANT = list(
      description = "Platinum-resistant relapse indicator: 1 = platinum-resistant relapse (PRR), 0 = platinum-sensitive relapse (PSR)",
      units = "(binary)",
      type = "binary",
      reference_category = "No reference category: the paper estimates a separate typical baseline CA-125 and KELIM in each stratum",
      notes = "Carrot 2025 Results 3.2: the platinum sensitivity status (PRR vs PSR) 'was used as a baseline covariate on CA-125_0 (CA-125 initial concentration) and on KELIM estimations in the CA-125 kinetic model'. Supplementary Material S1 selects KELIM = KELIM_PRR if Cohort = PRR and KELIM_PSR if Cohort = PSR, and likewise for C3(t=0). The cohort had 28 PRR and 22 PSR patients (Table 1). PRR / PSR follow the BOLD trial definitions (progression < 6 months vs >= 6 months after the last platinum dose).",
      source_name = "Cohort (PRR / PSR)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 50,
    n_studies = 1,
    sex_female_pct = 100,
    disease_state = "Recurrent high-grade ovarian cancer previously treated with at least one line of platinum-taxane chemotherapy; 28 platinum-resistant relapse (PRR) and 22 platinum-sensitive relapse (PSR) patients",
    dose_range = "Bevacizumab + olaparib + durvalumab (BOLD regimen); no drug concentrations available, so the K-PD dose is fixed to 1 a.u. per treatment administration",
    regions = "France (GINECO group, BOLD phase II trial NCT04015739)",
    n_observations = "198 CA-125 values (median 4 per patient) over the first 100 treatment days",
    notes = "50 of the 74 BOLD patients had >= 3 CA-125 values and >= 2 tumour-size measurements in the first 100 days (Figure 2). Table 1: FIGO stage III in 50.0% (PRR) / 72.7% (PSR); BRCA1/2-mutated tumours in 3.6% (PRR) / 36.4% (PSR); >= 3 previous chemotherapy lines in 53.5% / 50.0%. Age and body size are not reported in this paper."
  )

  ini({
    # K-PD treatment compartment and transit compartment share one rate K.
    lkel <- log(1.16); label("K-PD treatment kinetic rate constant K (1/day)") # Table 2, CA-125 kinetics model, K = 1.16 (RSE 48.19%)
    lec50 <- log(0.29); label("K-PD transit amount giving 50% inhibition of CA-125 production, EC50 (a.u.)") # Table 2, CA-125 kinetics model, EC50 = 0.29 AU (RSE 73.33%)
    lkin <- log(1.06); label("Zero-order CA-125 production rate KPROD (kU/L/day)") # Table 2, CA-125 kinetics model, KPROD = 1.06 (RSE 60.75%)

    # Separate typical values per platinum-sensitivity stratum.
    lrbase_ca125_resistant <- log(559); label("Baseline CA-125, platinum-resistant relapse (kU/L)") # Table 2, CA-125 kinetics model, CA-125 PRR,0 = 559.00 (RSE 30.95%)
    lrbase_ca125_sensitive <- log(141); label("Baseline CA-125, platinum-sensitive relapse (kU/L)") # Table 2, CA-125 kinetics model, CA-125 PSR,0 = 141.00 (RSE 35.67%)
    lkout_resistant <- log(0.01); label("CA-125 elimination rate constant KELIM, platinum-resistant relapse (1/day)") # Table 2, CA-125 kinetics model, KELIM PRR = 0.01 (RSE 45.92%)
    lkout_sensitive <- log(0.02); label("CA-125 elimination rate constant KELIM, platinum-sensitive relapse (1/day)") # Table 2, CA-125 kinetics model, KELIM PSR = 0.02 (RSE 50.05%)

    # IIV: Table 2 'IIV' column is the omega^2 variance of the log-normal
    # random effect (Methods Eq. 2). One eta per parameter is shared by both
    # strata (Table 2 prints the same IIV on the PRR and PSR rows).
    etalkel ~ 1.26 # Table 2, CA-125 kinetics model, K IIV = 1.26 (shrinkage 91.02%)
    etalkin ~ 2.92 # Table 2, CA-125 kinetics model, KPROD IIV = 2.92 (shrinkage 39.11%)
    etalec50 ~ 0.31 # Table 2, CA-125 kinetics model, EC50 IIV = 0.31 (shrinkage 93.34%)
    etalrbase_ca125 ~ 2.75 # Table 2, CA-125 kinetics model, CA-125 0 IIV = 2.75 (shrinkage 1.20%)
    etalkout ~ 1.67 # Table 2, CA-125 kinetics model, KELIM IIV = 1.67 (shrinkage 22.53%)

    # CA-125 was log-transformed and W is an additive error on the log scale.
    expSd_ca125 <- 0.29; label("CA-125 additive residual SD on the log scale, W (log kU/L)") # Table 2, CA-125 kinetics model, W = 0.29 (RSE 13.76%)
  })

  model({
    kel <- exp(lkel + etalkel)
    ec50 <- exp(lec50 + etalec50)
    kin <- exp(lkin + etalkin)

    # PLATINUM_RESISTANT selects the stratum-specific typical value.
    rbase_ca125 <- exp(
      lrbase_ca125_resistant * PLATINUM_RESISTANT +
        lrbase_ca125_sensitive * (1 - PLATINUM_RESISTANT) +
        etalrbase_ca125
    )
    kout <- exp(
      lkout_resistant * PLATINUM_RESISTANT +
        lkout_sensitive * (1 - PLATINUM_RESISTANT) +
        etalkout
    )

    # Supplementary Material S1: C1 = depot, C2 = transit1, C3 = ca125.
    d/dt(depot) <- -kel * depot
    d/dt(transit1) <- kel * depot - kel * transit1
    effect <- 1 - transit1 / (ec50 + transit1)
    ca125(0) <- rbase_ca125
    d/dt(ca125) <- kin * effect - kout * ca125

    ca125 ~ lnorm(expSd_ca125)
  })
}
