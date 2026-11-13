Saravanakumar_2026_plozalizumab_human <- function() {
  description <- "Two-compartment population PK/receptor-occupancy model of plozalizumab (TAK-202, MLN1202), a humanized anti-CCR2 monoclonal antibody and the parental antibody of the antibody-drug conjugate TAK-500, in healthy adults given single 1-h IV infusions of 0.1-10 mg/kg (Saravanakumar 2026, Supplementary Table S5). Parallel linear (kel) and Michaelis-Menten (CCR2 target-mediated) elimination from the central compartment; correlated IIV on kel and Km. Serum concentration drives an inhibitory Emax reduction of unbound CCR2 on monocytes, freeCCR2 = RO0 * (1 - Emax * Cc / (IC50 + Cc)) in MESF, with RO = 100 * (1 - freeCCR2 / RO0) the derived receptor occupancy (%). Parameters are absolute (no body-weight covariate); the authors' simulations use a 70 kg adult. Dose in nmol (1 mg = 6.725 nmol, MW 148,700 g/mol)."
  reference <- "Saravanakumar A, Deng R, Dong L, Haridas S, Xia CQ, Appleman VA, Abu-Yousif AO, Piatkov K, Wang H. Population-Based Modeling to Predict Human PK/PD of TAK-500, an Anti-CCR2 Antibody-Drug Conjugate for First-in-Human Study in Cancer Patients. Clin Transl Sci. 2026;19(4):e70522. doi:10.1111/cts.70522"
  vignette <- "Saravanakumar_2026_tak_500"
  units <- list(time = "h", dosing = "nmol", concentration = "nM (Cc); freeCCR2 in MESF; RO in %")

  # Saravanakumar 2026 Equations 19-21 and Figure 3B: Ac/At = plozalizumab
  # amount in the central / peripheral compartment.
  compartmentData <- list(
    central = list(analyte = "plozalizumab", units = "nmol", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "plozalizumab", units = "nmol", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = 36,
    n_studies = 1,
    disease_state = "Healthy adult volunteers (phase I single-ascending-dose study of MLN1202).",
    dose_range = "Single 1-h IV infusion of 0.1, 0.3, 1, 3, 6 or 10 mg/kg (6 active subjects per dose).",
    notes = "Saravanakumar 2026 Methods 'TAK-202 From Phase I Clinical Study': 48 subjects enrolled (6 active + 2 placebo per cohort), data available from 36. Dense PK sampling; unbound CCR2 on CD14+ monocytes measured with a FITC-labelled competing antibody (mu1D9-FITC) and reported in MESF. Demographics are not reported. Fitted in Phoenix NLME (FOCE-ELS)."
  )

  ini({
    # PK parameters -- Saravanakumar 2026 Supplementary Table S5 'Point estimate'.
    lvc <- log(3.184); label("Central volume of distribution (L)") # Table S5 'TV (Vc), L' 3.184 (%RSE 4)
    lkel <- log(0.002); label("Linear elimination rate constant from the central compartment (1/h)") # Table S5 'TV (kel), 1/h' 0.002 (%RSE 23)
    lvmax <- log(1.72); label("Maximum Michaelis-Menten (CCR2-mediated) elimination rate (nmol/h)") # Table S5 'TV (Vmax), nmole/h' 1.72 (%RSE 29)
    lkm <- log(46.6); label("Antibody concentration at half-maximal Michaelis-Menten elimination (nM)") # Table S5 'TV (Km), nM' 46.6 (%RSE 37)
    lk12 <- log(0.0091); label("Distribution rate constant from central to peripheral compartment, kpt (1/h)") # Table S5 'TV (kpt), 1/h' 0.0091 (%RSE 22)
    lk21 <- log(0.0128); label("Distribution rate constant from peripheral to central compartment, ktp (1/h)") # Table S5 'TV (ktp),1/h' 0.0128 (%RSE 33)

    # RO parameters -- Table S5 'RO parameters'.
    lrbase <- log(21509); label("Baseline unbound CCR2 on monocytes, RO0 (MESF)") # Table S5 'TV (ROb), MESF' 21509 (%RSE 3)
    limax <- log(0.522); label("Maximal fractional reduction of unbound CCR2, Emax (fraction)") # Table S5 'TV (Emax), 1/h' 0.522 (%RSE 5); unit printed as 1/h, the quantity is a fraction (Equation 22)
    lic50 <- log(5.71); label("Plozalizumab concentration giving half of Emax, IC50 (nM)") # Table S5 'TV (IC50), nM' 5.71 (%RSE 58)

    # IIV -- Table S5 'Inter-subject variability, eta' and 'omega' rows, read
    # as the Phoenix Omega variance x 100 (see the vignette for the evidence).
    etalvc ~ 0.02 # Table S5 'TV (Vc)' eta 2 (shrinkage 5%)
    etalkel + etalkm ~ c(0.25, 0.21, 0.60) # Table S5 'TV (kel)' eta 25, 'omega Km,Kel' 21, 'TV (Km)' eta 60
    etalrbase + etalimax ~ c(0.02, 0.02, 0.04) # Table S5 'TV (ROb)' eta 2, 'omega ROb,Emax' 2, 'TV (Emax)' eta 4

    # Residual error -- Equations 23-24, proportional.
    propSd <- 0.16; label("Proportional residual error, plozalizumab concentration (fraction)") # Table S5 PK 'Proportional error (%)' 16
    propSd_freeCCR2 <- 0.23; label("Proportional residual error, unbound CCR2 (fraction)") # Table S5 RO 'Proportional error (%)' 23
  })

  model({
    vc <- exp(lvc + etalvc)
    kel <- exp(lkel + etalkel)
    vmax <- exp(lvmax)
    km <- exp(lkm + etalkm)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    rbase <- exp(lrbase + etalrbase)
    imax <- exp(limax + etalimax)
    ic50 <- exp(lic50)

    # Equations 19-20
    d/dt(central) <- k21 * peripheral1 - kel * central - vmax * central / (km * vc + central) - k12 * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Equations 21-22: concentration and unbound CCR2 (MESF); RO (%) is the
    # fractional loss of unbound receptor (Figure 4B).
    Cc <- central / vc
    freeCCR2 <- rbase * (1 - imax * Cc / (ic50 + Cc))
    RO <- 100 * (1 - freeCCR2 / rbase)

    # Equations 23-24
    Cc ~ prop(propSd)
    freeCCR2 ~ prop(propSd_freeCCR2)
  })
}
