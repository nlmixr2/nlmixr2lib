Yates_2023_tgi_mayneord <- function() {
  description <- paste(
    "Theoretical (illustrative; no data fitted). Mayneord linear-radial",
    "tumor-growth law (dV/dt = k * V^(2/3)) with an Emax drug effect on the",
    "radial growth rate, driven by a one-compartment i.v. bolus PK model.",
    "Case study 2 of Yates and Mistry (2023): a sub-exponential growth law",
    "yields a log-linear rather than sigmoidal dose-response shape.",
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
    disease_state = "solid tumor (generic; sub-exponential linear-radial growth)",
    dose_range = paste(
      "illustrative single and repeat i.v. bolus doses; the published figures",
      "sweep CL from 0.01 to 10 L/day at Vd = 1 L to vary the PK half-life"
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
    # C(t) = (D / Vd) * exp(-a * t) with a = CL / Vd, as for case study 1
    # (p. 1595). Encoded as the equivalent one-compartment ODE.
    lcl <- fixed(log(0.01)); label("Clearance CL (L/day)")           # Figure 5 caption (CL = 0.01; the figure sweeps CL = 0.01-10)
    lvc <- fixed(log(1))   ; label("Volume of distribution Vd (L)")  # Figure 5 caption (Vd = 1)

    # ---- Mayneord linear-radial tumor growth ------------------------------
    # dV/dt = k * V^(2/3); the cube-root of volume (the tumor radius) grows
    # linearly at k/3 (p. 1598; Appendix S1 "Analytical solution to
    # Mayneord growth law").
    lrbase <- fixed(log(1))    ; label("Initial tumor volume V0 (L)")                       # Figure 5 caption (V0 = 1)
    lp     <- fixed(log(0.005)); label("Mayneord radial growth constant k (L^(1/3)/day)")   # Figure 5 caption (k = 0.005)

    # ---- Emax drug effect on the radial growth rate ----------------------
    # dV/dt = V^(2/3) * (k - Emax * C / (EC50 + C)).
    # The paper notes [Emax] = L.T^-1 in this model (i.e. Emax carries the
    # same L^(1/3)/day dimensions as k, acting on the radius), and that
    # [Emax/a] = L so the drug effect is a reduction of tumor radius
    # (p. 1598). Figure 5 uses Emax = 0.01; Figure 6 uses Emax = 0.0025.
    lemax <- fixed(log(0.01)); label("Maximum radial growth-inhibition rate Emax (L^(1/3)/day)")  # Figure 5 caption (Emax = 0.01)
    lec50 <- fixed(log(1))   ; label("Concentration for 50% of maximal effect EC50 (dose unit/L)")  # Figure 5 caption (EC50 = 1)
  })

  model({
    cl    <- exp(lcl)
    vc    <- exp(lvc)
    rbase <- exp(lrbase)
    p     <- exp(lp)
    emax  <- exp(lemax)
    ec50  <- exp(lec50)

    kel <- cl / vc

    tumor_vol(0) <- rbase

    d/dt(central) <- -kel * central
    Cc <- central / vc

    drugEffect <- emax * Cc / (ec50 + Cc)

    # Paper Eq. dV/dt = V^(2/3) * (k - Emax * C(t) / (EC50 + C(t))), p. 1598.
    d/dt(tumor_vol) <- tumor_vol^(2 / 3) * (p - drugEffect)
  })
}
