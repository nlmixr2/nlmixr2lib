Carrot_2025_ovarianCancer_tgi <- function() {
  description <- "K-PD tumour growth inhibition model of the sum of longest target-lesion diameters during bevacizumab + olaparib + durvalumab in recurrent high-grade ovarian cancer (BOLD phase II GINECO trial). Arbitrary unit dose into a treatment compartment and a transit compartment; the transit amount inhibits zero-order tumour proliferation KPROL through 1 - C2/(EC50 + C2); first-order tumour reduction KREDUC. Fit on its own to check the feasibility of the joint CA-125 + tumour-size model; platinum sensitivity status was tested and not retained."
  reference <- "Carrot A, Corbaux P, Colomban O, Lebreton C, Gladieff L, Tredan O, Selle F, Abdeddaim C, Leary A, D'Hondt V, Dubot C, Reverdy T, You B, Freyer G. Prognostic Value of a Joint K-PD Model With Tumor Size Dynamics and CA-125 Kinetics in Recurrent Ovarian Cancer Patients: BOLD Phase II GINECO Study. CPT Pharmacometrics Syst Pharmacol. 2025;14(8):1310-1321. doi:10.1002/psp4.70025"
  vignette <- "Carrot_2025_ovarianCancer"
  paper_specific_etas <- c("etalkprol", "etalkreduc")
  paper_specific_residual_sds <- c("expSd_tumor_size")
  units <- list(
    time = "day",
    dosing = "a.u. (arbitrary K-PD dose, fixed to 1 per treatment administration)",
    concentration = "mm for tumour size (sum of longest diameters, RECIST 1.1)"
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
    )
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    PLATINUM_RESISTANT = list(
      description = "Platinum-resistant relapse indicator: 1 = platinum-resistant relapse (PRR), 0 = platinum-sensitive relapse (PSR)",
      units = "(binary)",
      type = "binary",
      notes = "Carrot 2025 Results 3.2: 'The integration of this covariate did not improve the TGI model and was not kept in this model.' It is retained in the CA-125 and joint models (Carrot_2025_ovarianCancer_ca125, Carrot_2025_ovarianCancer_ca125_tgi).",
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
    n_observations = "120 tumour sizes (median 2 per patient) over the first 100 treatment days; diameters below 10 mm (3.92% of values) set to LOQ/2",
    notes = "50 of the 74 BOLD patients had >= 3 CA-125 values and >= 2 tumour-size measurements in the first 100 days (Figure 2). Table 1: FIGO stage III in 50.0% (PRR) / 72.7% (PSR); BRCA1/2-mutated tumours in 3.6% (PRR) / 36.4% (PSR). Age and body size are not reported in this paper."
  )

  ini({
    lkel <- log(0.33); label("K-PD treatment kinetic rate constant K (1/day)") # Table 2, TGI model, K = 0.33 (RSE 31.83%)
    lec50 <- log(0.29); label("K-PD transit amount giving 50% inhibition of tumour proliferation, EC50 (a.u.)") # Table 2, TGI model, EC50 = 0.29 AU (RSE 81.12%)
    lrbase_tumor_size <- log(58.1); label("Baseline tumour size TS0 (mm)") # Table 2, TGI model, TS0 = 58.10 (RSE 11.24%)
    lkreduc <- log(0.02); label("First-order tumour reduction rate constant KREDUC (1/day)") # Table 2, TGI model, KREDUC = 0.02 (RSE 102.25%)
    lkprol <- log(3.67); label("Zero-order tumour proliferation rate KPROL (mm/day)") # Table 2, TGI model, KPROL = 3.67 (RSE 27.25%)

    # IIV: Table 2 'IIV' column is the omega^2 variance (Methods Eq. 2).
    etalkel ~ 0.83 # Table 2, TGI model, K IIV = 0.83 (shrinkage 43.10%)
    etalec50 ~ 0.55 # Table 2, TGI model, EC50 IIV = 0.55 (shrinkage 64.69%)
    etalrbase_tumor_size ~ 0.61 # Table 2, TGI model, TS0 IIV = 0.61 (shrinkage 1.25%)
    etalkreduc ~ 0.82 # Table 2, TGI model, KREDUC IIV = 0.82 (shrinkage 69.16%)
    etalkprol ~ 1.73 # Table 2, TGI model, KPROL IIV = 1.73 (shrinkage 81.27%)

    # Tumour size was log-transformed and W is an additive error on the log scale.
    expSd_tumor_size <- 0.084; label("Tumour-size additive residual SD on the log scale, W (log mm)") # Table 2, TGI model, W = 0.084 (RSE 19.65%)
  })

  model({
    kel <- exp(lkel + etalkel)
    ec50 <- exp(lec50 + etalec50)
    kreduc <- exp(lkreduc + etalkreduc)
    kprol <- exp(lkprol + etalkprol)

    # Figure 1 / Discussion: C1 = depot, C2 = transit1, C3 = tumor_size.
    d/dt(depot) <- -kel * depot
    d/dt(transit1) <- kel * depot - kel * transit1
    effect <- 1 - transit1 / (ec50 + transit1)
    tumor_size(0) <- exp(lrbase_tumor_size + etalrbase_tumor_size)
    d/dt(tumor_size) <- kprol * effect - kreduc * tumor_size

    tumor_size ~ lnorm(expSd_tumor_size)
  })
}
