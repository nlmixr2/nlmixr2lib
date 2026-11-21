Carrot_2025_ovarianCancer_ca125_tgi <- function() {
  description <- "Joint K-PD model of tumour size and CA-125 kinetics (KELIM-joint model, Wilbaux 2014 form) during bevacizumab + olaparib + durvalumab in recurrent platinum-resistant or platinum-sensitive high-grade ovarian cancer (BOLD phase II GINECO trial). Arbitrary unit dose into a treatment compartment and a transit compartment; the transit amount inhibits zero-order tumour proliferation KPROL through 1 - C2/(EC50 + C2) with first-order tumour reduction KREDUC; CA-125 is produced at a basal rate KPROD1 plus a tumour-driven rate KPROD2 * exp(K2 * dTS/dt) and eliminated at KELIM. Baseline CA-125 and KELIM are estimated separately for platinum-resistant and platinum-sensitive relapse."
  reference <- "Carrot A, Corbaux P, Colomban O, Lebreton C, Gladieff L, Tredan O, Selle F, Abdeddaim C, Leary A, D'Hondt V, Dubot C, Reverdy T, You B, Freyer G. Prognostic Value of a Joint K-PD Model With Tumor Size Dynamics and CA-125 Kinetics in Recurrent Ovarian Cancer Patients: BOLD Phase II GINECO Study. CPT Pharmacometrics Syst Pharmacol. 2025;14(8):1310-1321. doi:10.1002/psp4.70025"
  vignette <- "Carrot_2025_ovarianCancer"
  paper_specific_compartments <- c("ca125")
  paper_specific_etas <- c("etalkprol", "etalkreduc", "etalkprod1", "etalkprod2", "etalk2")
  paper_specific_residual_sds <- c("expSd_tumor_size", "expSd_ca125")
  units <- list(
    time = "day",
    dosing = "a.u. (arbitrary K-PD dose, fixed to 1 per treatment administration)",
    concentration = "kU/L (= U/mL) for CA-125; mm for tumour size"
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
    tumor_size = list(
      analyte = "sum of longest diameters of target lesions",
      units = "mm",
      specimen = "tumor",
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
      notes = "Carrot 2025 Results 3.2: the platinum sensitivity status (PRR vs PSR) was used as a baseline covariate on CA-125_0 and on KELIM 'in the CA-125 kinetic model, and in the joint model'. Supplementary Material S1 (joint model) selects KELIM_CA125 = KELIM_PRR if Cohort = PRR and KELIM_PSR if Cohort = PSR, and likewise for C4(t=0). It does not act on the tumour-size parameters. The cohort had 28 PRR and 22 PSR patients (Table 1).",
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
    n_observations = "198 CA-125 values and 120 tumour sizes (medians 4 and 2 per patient) over the first 100 treatment days; diameters below 10 mm (3.92% of values) set to LOQ/2",
    notes = "50 of the 74 BOLD patients had >= 3 CA-125 values and >= 2 tumour-size measurements in the first 100 days (Figure 2). Table 1: FIGO stage III in 50.0% (PRR) / 72.7% (PSR); BRCA1/2-mutated tumours in 3.6% (PRR) / 36.4% (PSR). Age and body size are not reported in this paper."
  )

  ini({
    lkel <- log(0.30); label("K-PD treatment kinetic rate constant K (1/day)") # Table 2, Joint model, K = 0.30 (RSE 41.69%)
    lec50 <- log(0.40); label("K-PD transit amount giving 50% inhibition of tumour proliferation, EC50 (a.u.)") # Table 2, Joint model, EC50 = 0.40 AU (RSE 80.15%)
    lrbase_tumor_size <- log(57.8); label("Baseline tumour size TS0 (mm)") # Table 2, Joint model, TS0 = 57.80 (RSE 11.26%)
    lkreduc <- log(0.02); label("First-order tumour reduction rate constant KREDUC (1/day)") # Table 2, Joint model, KREDUC = 0.02 (RSE 51.23%)
    lkprol <- log(3.35); label("Zero-order tumour proliferation rate KPROL (mm/day)") # Table 2, Joint model, KPROL = 3.35 (RSE 23.85%)
    lkprod1 <- log(0.19); label("Basal (normal-tissue) CA-125 production rate KPROD1 (kU/L/day)") # Table 2, Joint model, KPROD 1 = 0.19 (RSE 128.65%)
    lkprod2 <- log(1.43); label("Tumour-driven CA-125 production rate KPROD2 (kU/L/day)") # Table 2, Joint model, KPROD 2 = 1.43 (RSE 61.61%)
    lk2 <- log(0.18); label("Factor linking KPROD2 to the tumour-size rate of change K2 (day/mm)") # Table 2, Joint model, K 2 = 0.18 (RSE 51.64%)

    # Separate typical values per platinum-sensitivity stratum.
    lrbase_ca125_resistant <- log(567); label("Baseline CA-125, platinum-resistant relapse (kU/L)") # Table 2, Joint model, CA-125 0,PRR = 567.00 (RSE 31.04%)
    lrbase_ca125_sensitive <- log(144); label("Baseline CA-125, platinum-sensitive relapse (kU/L)") # Table 2, Joint model, CA-125 0,PSR = 144.00 (RSE 36.18%)
    lkout_resistant <- log(0.01); label("CA-125 elimination rate constant KELIM, platinum-resistant relapse (1/day)") # Table 2, Joint model, KELIM PRR = 0.01 (RSE 46.95%)
    lkout_sensitive <- log(0.03); label("CA-125 elimination rate constant KELIM, platinum-sensitive relapse (1/day)") # Table 2, Joint model, KELIM PSR = 0.03 (RSE 46.13%)

    # IIV: Table 2 'IIV' column is the omega^2 variance (Methods Eq. 2). One
    # eta per parameter is shared by both platinum strata.
    etalkel ~ 0.09 # Table 2, Joint model, K IIV = 0.09 (shrinkage 72.70%)
    etalec50 ~ 2.79 # Table 2, Joint model, EC50 IIV = 2.79 (shrinkage 34.49%)
    etalrbase_tumor_size ~ 0.61 # Table 2, Joint model, TS0 IIV = 0.61 (shrinkage 0.90%)
    etalkreduc ~ 0.49 # Table 2, Joint model, KREDUC IIV = 0.49 (shrinkage 73.02%)
    etalkprol ~ 0.81 # Table 2, Joint model, KPROL IIV = 0.81 (shrinkage 83.04%)
    etalkprod1 ~ 1.51 # Table 2, Joint model, KPROD 1 IIV = 1.51 (shrinkage 82.84%)
    etalkprod2 ~ 1.71 # Table 2, Joint model, KPROD 2 IIV = 1.71 (shrinkage 38.09%)
    etalrbase_ca125 ~ 2.72 # Table 2, Joint model, CA-125 0 IIV = 2.72 (shrinkage 0.85%)
    etalkout ~ 1.93 # Table 2, Joint model, KELIM IIV = 1.93 (shrinkage 18.33%)
    etalk2 ~ 0.13 # Table 2, Joint model, K 2 IIV = 0.13 (shrinkage 90.57%)

    # Both outputs were log-transformed; W is an additive error on the log scale.
    expSd_tumor_size <- 0.084; label("Tumour-size additive residual SD on the log scale, W TS (log mm)") # Table 2, Joint model, W TS = 0.084 (RSE 11.81%)
    expSd_ca125 <- 0.29; label("CA-125 additive residual SD on the log scale, W CA (log kU/L)") # Table 2, Joint model, W CA = 0.29 (RSE 15.21%)
  })

  model({
    kel <- exp(lkel + etalkel)
    ec50 <- exp(lec50 + etalec50)
    kreduc <- exp(lkreduc + etalkreduc)
    kprol <- exp(lkprol + etalkprol)
    kprod1 <- exp(lkprod1 + etalkprod1)
    kprod2 <- exp(lkprod2 + etalkprod2)
    k2 <- exp(lk2 + etalk2)

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

    # Supplementary Material S1: C1 = depot, C2 = transit1, C3 = tumor_size,
    # C4 = ca125. VARTS is the instantaneous tumour-size rate of change.
    d/dt(depot) <- -kel * depot
    d/dt(transit1) <- kel * depot - kel * transit1
    effect <- 1 - transit1 / (ec50 + transit1)
    tumor_size(0) <- exp(lrbase_tumor_size + etalrbase_tumor_size)
    varts <- kprol * effect - kreduc * tumor_size
    d/dt(tumor_size) <- varts
    ca125(0) <- rbase_ca125
    d/dt(ca125) <- kprod1 + kprod2 * exp(k2 * varts) - kout * ca125

    tumor_size ~ lnorm(expSd_tumor_size)
    ca125 ~ lnorm(expSd_ca125)
  })
}
