Saravanakumar_2026_tak_500_human <- function() {
  description <- "Human PK/receptor-occupancy projection for TAK-500, an anti-CCR2 immune-cell-directed antibody-drug conjugate carrying the STING agonist payload dazostinag (TAK-676), used to plan the first-in-human study (Saravanakumar 2026). Total-antibody disposition is taken from the human population PK model of the parental antibody plozalizumab (TAK-202; two compartments, linear kel plus Michaelis-Menten CCR2-mediated elimination, Supplementary Table S5). The ADC structure of the monkey TAK-500 model is kept: the conjugate (central/peripheral1) deconjugates at the monkey kdec to deconjugated antibody (central_dar0/peripheral1_dar0), both species share the disposition and the Michaelis-Menten saturation term. Outputs are conjugated payload (DAR 4.6 x ADC, Cc) and total antibody (Cc_tab) in nM, and unbound CCR2 on monocytes (freeCCR2, MESF) and receptor occupancy (RO, %) driven by total antibody through the plozalizumab human Imax model. All values are carried over and fixed; parameters are absolute for a 70 kg adult. Dose in nmol of antibody (1 mg = 6.725 nmol, MW 148,700 g/mol)."
  reference <- "Saravanakumar A, Deng R, Dong L, Haridas S, Xia CQ, Appleman VA, Abu-Yousif AO, Piatkov K, Wang H. Population-Based Modeling to Predict Human PK/PD of TAK-500, an Anti-CCR2 Antibody-Drug Conjugate for First-in-Human Study in Cancer Patients. Clin Transl Sci. 2026;19(4):e70522. doi:10.1111/cts.70522"
  vignette <- "Saravanakumar_2026_tak_500"
  units <- list(
    time = "h",
    dosing = "nmol",
    concentration = "nM (Cc conjugated TAK-676, Cc_tab total antibody); freeCCR2 in MESF; RO in %"
  )

  # Saravanakumar 2026 Methods 'TAK-500 Human PK Projection' applies the
  # monkey TAK-500 structure (Equations 1-6) with human TAK-202 disposition.
  compartmentData <- list(
    central = list(analyte = "TAK-500 (conjugated antibody)", units = "nmol", specimen = "serum", verified = TRUE),
    peripheral1 = list(
      analyte = "TAK-500 (conjugated antibody)",
      units = "nmol",
      specimen = "not applicable",
      verified = TRUE
    ),
    central_dar0 = list(
      analyte = "deconjugated anti-CCR2 antibody",
      units = "nmol",
      specimen = "serum",
      verified = TRUE
    ),
    peripheral1_dar0 = list(
      analyte = "deconjugated anti-CCR2 antibody",
      units = "nmol",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = 36,
    n_studies = 1,
    disease_state = "Projection for adults with cancer entering the first-in-human study of TAK-500; the disposition parameters come from healthy volunteers given plozalizumab.",
    dose_range = "Projected single 1-h IV infusions of 50, 150 and 500 ug/kg in a 70 kg adult (Figure 4, Table 3).",
    notes = "No TAK-500 human data were fitted. Saravanakumar 2026 projection assumptions: (1) TAK-500 total antibody has the human PK of plozalizumab (36 healthy volunteers, Supplementary Table S5); (2) kdec in humans equals the monkey estimate (Table 2); (3) clinical TAK-500 has a mean DAR of 4.6. n_subjects is the size of the plozalizumab dataset behind the disposition parameters."
  )

  ini({
    # Human plozalizumab disposition -- Supplementary Table S5 point estimates.
    lvc <- fixed(log(3.184)); label("Central volume of distribution (L)") # Table S5 'TV (Vc), L' 3.184
    lkel <- fixed(log(0.002)); label("Linear elimination rate constant from the central compartment (1/h)") # Table S5 'TV (kel), 1/h' 0.002
    lvmax <- fixed(log(1.72)); label("Maximum Michaelis-Menten (CCR2-mediated) elimination rate (nmol/h)") # Table S5 'TV (Vmax), nmole/h' 1.72
    lkm <- fixed(log(46.6)); label("Antibody concentration at half-maximal Michaelis-Menten elimination (nM)") # Table S5 'TV (Km), nM' 46.6
    lk12 <- fixed(log(0.0091)); label("Distribution rate constant from central to peripheral compartment, kpt (1/h)") # Table S5 'TV (kpt), 1/h' 0.0091
    lk21 <- fixed(log(0.0128)); label("Distribution rate constant from peripheral to central compartment, ktp (1/h)") # Table S5 'TV (ktp),1/h' 0.0128

    # ADC-specific values carried from the monkey TAK-500 model.
    lkdec <- fixed(log(0.0324)); label("First-order deconjugation rate constant of the ADC to deconjugated antibody (1/h)") # Table 2 'TV (kdec), 1/h' 0.0324; Methods assumption 2 'kdec in human ... same as that in monkey'
    dar <- fixed(4.6); label("Mean drug-to-antibody ratio of clinical TAK-500 (mol TAK-676 per mol antibody)") # Methods assumption 3 'DAR of 4.6 (average)'; Figure 4 caption says DAR 4 (see vignette)

    # Human CCR2 occupancy -- Supplementary Table S5 'RO parameters'.
    lrbase <- fixed(log(21509)); label("Baseline unbound CCR2 on monocytes, RO0 (MESF)") # Table S5 'TV (ROb), MESF' 21509
    limax <- fixed(log(0.522)); label("Maximal fractional reduction of unbound CCR2, Emax (fraction)") # Table S5 'TV (Emax), 1/h' 0.522 (a fraction; unit misprinted)
    lic50 <- fixed(log(5.71)); label("Total-antibody concentration giving half of Emax, IC50 (nM)") # Table S5 'TV (IC50), nM' 5.71

    # IIV carried from Table S5 (eta column = Phoenix Omega variance x 100).
    etalvc ~ fixed(0.02) # Table S5 'TV (Vc)' eta 2
    etalkel + etalkm ~ fixed(c(0.25, 0.21, 0.60)) # Table S5 'TV (kel)' eta 25, 'omega Km,Kel' 21, 'TV (Km)' eta 60
    etalrbase + etalimax ~ fixed(c(0.02, 0.02, 0.04)) # Table S5 'TV (ROb)' eta 2, 'omega ROb,Emax' 2, 'TV (Emax)' eta 4

    # Residual error carried from the plozalizumab human model (Table S5);
    # the conjugated-payload error is assumed equal to the antibody error.
    propSd <- fixed(0.16); label("Proportional residual error, conjugated payload (fraction)") # Table S5 PK 'Proportional error (%)' 16, applied to conjugated payload by assumption
    propSd_tab <- fixed(0.16); label("Proportional residual error, total antibody (fraction)") # Table S5 PK 'Proportional error (%)' 16
    propSd_freeCCR2 <- fixed(0.23); label("Proportional residual error, unbound CCR2 (fraction)") # Table S5 RO 'Proportional error (%)' 23
  })

  model({
    vc <- exp(lvc + etalvc)
    kel <- exp(lkel + etalkel)
    vmax <- exp(lvmax)
    km <- exp(lkm + etalkm)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    kdec <- exp(lkdec)
    rbase <- exp(lrbase + etalrbase)
    imax <- exp(limax + etalimax)
    ic50 <- exp(lic50)

    # Equations 1-4 with human disposition parameters.
    mm_denom <- km * vc + central + central_dar0
    d/dt(central) <- -kel * central - vmax * central / mm_denom - k12 * central + k21 * peripheral1 - kdec * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(central_dar0) <- -kel * central_dar0 - vmax * central_dar0 / mm_denom - k12 * central_dar0 + k21 * peripheral1_dar0 + kdec * central
    d/dt(peripheral1_dar0) <- k12 * central_dar0 - k21 * peripheral1_dar0

    # Equations 5-6
    Cc_tab <- (central + central_dar0) / vc
    Cc <- dar * central / vc

    # Equation 22 driven by total antibody; RO (%) as plotted in Figure 4B.
    freeCCR2 <- rbase * (1 - imax * Cc_tab / (ic50 + Cc_tab))
    RO <- 100 * (1 - freeCCR2 / rbase)

    Cc ~ prop(propSd)
    Cc_tab ~ prop(propSd_tab)
    freeCCR2 ~ prop(propSd_freeCCR2)
  })
}
