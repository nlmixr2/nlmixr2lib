Yates_2023_tgi_exponential <- function() {
  description <- paste(
    "Theoretical (illustrative; no data fitted). Exponential tumor-growth",
    "law with a sigmoidal Emax drug effect on the growth rate, driven by a",
    "one-compartment i.v. bolus PK model. Case study 1 of Yates and Mistry",
    "(2023), used to show that the location (ED50) and steepness of a",
    "dose-response curve are set jointly by potency (EC50), efficacy",
    "(Emax) and the drug's elimination rate constant a = CL/Vd.",
    sep = " "
  )
  reference <- paste(
    "Yates JWT, Mistry HB. Skipping a pillar does not make for strong",
    "foundations: Pharmacokinetic-pharmacodynamic reasoning behind the",
    "shape of dose-response relationships in oncology.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12:1591-1601.",
    "doi:10.1002/psp4.13020.",
    sep = " "
  )
  vignette <- "Yates_2023_dose_response_shape"
  units <- list(
    time = "day",
    dosing = "arbitrary dose unit",
    concentration = "arbitrary dose unit/L"
  )

  compartmentData <- list(
    central   = list(analyte = "Drug (generic)", units = "arbitrary dose unit", specimen = "plasma", verified = FALSE),
    tumor_vol = list(analyte = "tumour_size", units = "L", specimen = "tumor", verified = FALSE)
  )

  population <- list(
    species = "not applicable (theoretical illustration; no subjects, no data fitted)",
    n_subjects = 0,
    n_studies = 0,
    disease_state = "solid tumor (generic; exponentially growing tumor burden)",
    dose_range = paste(
      "illustrative single and repeat i.v. bolus doses; the published figures",
      "sweep dose over several orders of magnitude and sweep CL from 0.01 to",
      "10 L/day at Vd = 1 L to vary the PK half-life"
    ),
    notes = paste(
      "Yates and Mistry (2023) is a theoretical / tutorial paper. No",
      "clinical or preclinical data were fitted, so there is no study",
      "population, no inter-individual variability and no residual-error",
      "model. Every parameter below is an illustrative constant chosen by",
      "the authors for the published figures, and is therefore encoded with",
      "fixed(). Full derivations are in Appendix S1 of the source."
    )
  )

  ini({
    # ---- One-compartment i.v. bolus PK ----------------------------------
    # The paper writes the PK analytically as C(t) = (D / Vd) * exp(-a * t)
    # with a = CL / Vd (p. 1595, "Case study: Exponential growth";
    # Appendix S1 Eq. between (1) and (2)). It is encoded here as the
    # equivalent one-compartment ODE so that repeat dosing and the
    # accumulation results of the paper follow directly.
    lcl <- fixed(log(0.01)); label("Clearance CL (L/day)")                     # Figure 2 caption (CL = 0.01; Figures 1/3/4 sweep CL = 0.01-10)
    lvc <- fixed(log(1))   ; label("Volume of distribution Vd (L)")            # Figure 2 caption (Vd = 1)

    # ---- Tumor growth ----------------------------------------------------
    lrbase <- fixed(log(1))    ; label("Initial tumor volume V0 (L)")          # Figure 2 caption (V0 = 1); initial condition V(0) = V0
    lp     <- fixed(log(0.005)); label("Exponential tumor growth rate k (1/day)")  # Figure 2 caption (k = 0.005)

    # ---- Sigmoidal Emax drug effect on the growth rate -------------------
    # dV/dt = V * (k - Emax * C^n / (EC50^n + C^n)).
    # Figure 2 uses Emax = 0.01 > k = 0.005, the regime in which the drug
    # can shrink the tumor. Figures 1, 3 and 4 instead use Emax = 0.0025 < k,
    # where the drug can only slow growth.
    lemax <- fixed(log(0.01)); label("Maximum fractional growth-inhibition rate Emax (1/day)")  # Figure 2 caption (Emax = 0.01)
    lec50 <- fixed(log(1))   ; label("Concentration for 50% of maximal effect EC50 (dose unit/L)")  # Figure 2 caption (EC50 = 1)
    # The base case of every published figure is the plain Emax model
    # (n = 1). The paper generalises the same case study to a steep PK/PD
    # relationship with Hill coefficient n (p. 1596; Appendix S1); set hill
    # above 1 to reproduce that variant.
    lhill <- fixed(log(1)); label("Hill coefficient n on the PK/PD relationship (unitless)")  # p. 1596, sigmoidal PK/PD generalisation (base case n = 1)
  })

  model({
    cl    <- exp(lcl)
    vc    <- exp(lvc)
    rbase <- exp(lrbase)
    p     <- exp(lp)
    emax  <- exp(lemax)
    ec50  <- exp(lec50)
    hill  <- exp(lhill)

    kel <- cl / vc

    tumor_vol(0) <- rbase

    d/dt(central) <- -kel * central
    Cc <- central / vc

    # Sigmoidal Emax inhibition of the fractional growth rate.
    drugEffect <- emax * Cc^hill / (ec50^hill + Cc^hill)

    # Paper Eq. dV/dt = V * (k - Emax * C(t) / (EC50 + C(t))), p. 1595.
    d/dt(tumor_vol) <- tumor_vol * (p - drugEffect)
  })
}
