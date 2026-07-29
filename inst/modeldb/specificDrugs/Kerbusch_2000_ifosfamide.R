Kerbusch_2000_ifosfamide <- function() {
  description <- "One-compartment population PK model for ifosfamide with autoinduction of CYP3A4-mediated metabolism implemented as ifosfamide-driven inhibition of enzyme-pool degradation (no lag time); estimated in 15 adults with soft tissue sarcoma receiving 9 or 12 g/m^2 as a 72-h continuous IV infusion."
  reference <- paste(
    "Kerbusch T., Huitema A. D. R., Ouwerkerk J., Keizer H. J.,",
    "Mathot R. A. A., Schellens J. H. M., Beijnen J. H. (2000).",
    "Evaluation of the autoinduction of ifosfamide metabolism by a",
    "population pharmacokinetic approach using NONMEM.",
    "British Journal of Clinical Pharmacology 49(6):555-561.",
    "doi:10.1046/j.1365-2125.2000.00217.x.",
    sep = " "
  )
  vignette <- "Kerbusch_2000_ifosfamide"
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 15L,
    n_studies      = 1L,
    age_range      = "23 to 72 years",
    age_median     = "49 years (mean)",
    weight_range   = "49 to 82 kg",
    weight_median  = "59 kg (mean)",
    sex_female_pct = 40,
    race_ethnicity = "not reported",
    disease_state  = "Soft tissue sarcoma (rhabdomyosarcoma n=3, neurofibrosarcoma n=1, osteosarcoma n=3, synoviumsarcoma n=2, leiomyosarcoma n=4, endometriumsarcoma n=2).",
    dose_range     = "9 or 12 g/m^2 ifosfamide as a 72-h continuous IV infusion, once every 4 weeks (one patient stopped at 48 h due to severe neurotoxicity, receiving 66% of the planned dose).",
    regions        = "Netherlands (Netherlands Cancer Institute, Amsterdam; Leiden University Medical Center)",
    notes          = "Open non-randomised phase II trial. Baseline demographics from Kerbusch 2000 Results paragraph 1. Comedication: mesna and bicarbonate (supportive care), anti-emetics, methylene blue (neurotoxicity antidote), plus 21 other drugs (mean 7 per patient) including anticoagulants, H2 antagonists, glucocorticosteroids, tricyclic antidepressants; none reported as CYP3A4 inhibitors or inducers."
  )

  ini({
    # All parameters from Kerbusch 2000 Table 1. Estimated in NONMEM V (FOCE)
    # with ADVAN6 TOL=5 and a proportional (exponential) IIV model on CL, V,
    # and K_enz,out; combined additive + proportional residual error.
    lcl    <- log(2.94)   ; label("Initial (preinduced) ifosfamide clearance CL (L/h)")                                                     # Kerbusch 2000 Table 1: CL = 2.94 L/h (SE 0.27)
    lvc    <- log(43.5)   ; label("Volume of distribution V (L)")                                                                            # Kerbusch 2000 Table 1: V  = 43.5 L (SE 2.9)
    lkdeg  <- log(0.0546) ; label("First-order enzyme degradation rate K_enz,out (1/h); equals K_enz,in by the A2(t=0)=1 steady-state constraint") # Kerbusch 2000 Table 1: K_enz,out = 0.0546 /h (SE 0.0078); induction half-life ln(2)/K_enz,out = 12.7 h
    lic50  <- log(30.7)   ; label("Ifosfamide concentration at 50% of the maximum inhibition of enzyme degradation IC50 (umol/L)")           # Kerbusch 2000 Table 1: IC50 = 30.7 umol/L (SE 4.8)

    # IIV: proportional / exponential on log-scale. omega^2 = log(1 + CV^2).
    # No IIV on IC50 ("the inclusion of interindividual variability did not
    # improve the fit"; see Kerbusch 2000 Discussion paragraph 5).
    etalcl   ~ 0.0583  # Kerbusch 2000 Table 1: IIV CL  = 24.5% CV; omega^2 = log(1 + 0.245^2) = 0.0583
    etalvc   ~ 0.0533  # Kerbusch 2000 Table 1: IIV V   = 23.4% CV; omega^2 = log(1 + 0.234^2) = 0.0533
    etalkdeg ~ 0.0502  # Kerbusch 2000 Table 1: IIV K_enz,out = 22.7% CV; omega^2 = log(1 + 0.227^2) = 0.0502

    # Combined additive + proportional residual error (Kerbusch 2000
    # Table 1 'Residual variability'). Proportional 13.6% of the
    # observed concentration; additive 0.0763 umol/L.
    propSd <- 0.136   ; label("Proportional residual error (fraction)")          # Kerbusch 2000 Table 1: P.E. = 13.6%
    addSd  <- 0.0763  ; label("Additive residual error (umol/L)")                # Kerbusch 2000 Table 1: A.E. = 0.0763 umol/L
  })

  model({
    # --- Individual parameters (log-normal IIV).
    cl    <- exp(lcl   + etalcl)
    vc    <- exp(lvc   + etalvc)
    kdeg  <- exp(lkdeg + etalkdeg)
    ic50  <- exp(lic50)

    # --- Observed plasma concentration. The model carries ifosfamide as an
    #     amount (umol) in central; dividing by V gives the plasma
    #     concentration in umol/L = uM, matching the IC50 unit.
    Cc <- central / vc

    # --- PK ODE system (Kerbusch 2000 Equations 1-2 and 4; one-compartment
    #     with autoinduction driven by the relative enzyme amount A2 in
    #     a hypothetical compartment 2).
    #
    #     Equation 1 (drug):
    #       dA1/dt = R - (CLapp / V) * A1 = R - CLapp * Cc
    #     Equation 2:
    #       CLapp = CL * A2
    #     Equation 4 (enzyme; IC50-mediated inhibition of degradation, with
    #     K_enz,in = K_enz,out enforced by the A2(0)=1 steady-state constraint
    #     stated below Eq. 3 of the paper):
    #       dA2/dt = K_enz,out * (1 - A2 * IC50 / (IC50 + Cp))
    #
    #     R (umol/h) enters as the infusion rate set on the dosing record.
    d/dt(central)  <- -cl * enz_pool * Cc
    d/dt(enz_pool) <-  kdeg * (1 - enz_pool * ic50 / (ic50 + Cc))

    # Relative enzyme amount A2 is normalised to unity at t = 0 (Kerbusch
    # 2000 Methods, last paragraph of 'Model development': 'Before
    # treatment with ifosfamide (t = 0) steady-state levels of the enzyme
    # were assumed... By defining A2,t=0 at a value of 1, K_enz,in could
    # be substituted for K_enz,out').
    enz_pool(0) <- 1.0

    # --- Combined residual error on plasma ifosfamide concentration.
    Cc ~ add(addSd) + prop(propSd)
  })
}
