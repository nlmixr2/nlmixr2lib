Hurtado_2014_levofloxacin_rat <- function() {
  description <- paste(
    "Preclinical (rat). Three-compartment population PK model for unbound",
    "levofloxacin in plasma and prostate interstitial fluid in male Wistar",
    "rats after a single 7 mg/kg IV bolus, with simultaneous fitting of",
    "total plasma concentrations (central, Vc) and free prostate ISF",
    "concentrations measured by microdialysis (effect compartment, apparent",
    "volume V3* = V_prostate / fu_prostate). Prostate kinetics are",
    "asymmetric: uptake from central is first-order (k13), efflux back to",
    "central combines a linear first-order term (k31) with a saturable",
    "Michaelis-Menten efflux (Vmax, kM) consistent with active transporter",
    "involvement. The standard central <-> peripheral1 disposition uses",
    "macro-constants CL, Q, Vc, Vp (Hurtado 2014)."
  )
  reference <- paste(
    "Hurtado FK, Weber B, Derendorf H, Hochhaus G, Dalla Costa T. (2014).",
    "Population pharmacokinetic modeling of the unbound levofloxacin",
    "concentrations in rat plasma and prostate tissue measured by",
    "microdialysis. Antimicrob Agents Chemother 58(2):678-685.",
    "doi:10.1128/AAC.01884-13"
  )
  vignette <- "Hurtado_2014_levofloxacin_rat"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list()

  population <- list(
    species        = "rat (Wistar, male)",
    n_subjects     = 7L,
    n_studies      = 1L,
    age_range      = "Adult; specific age not reported",
    weight_range   = "0.25-0.35 kg",
    sex_female_pct = 0,
    race_ethnicity = NA,
    disease_state  = paste(
      "Healthy male Wistar rats anaesthetised with urethane (1.25 g/kg IP),",
      "carotid artery cannulated for blood sampling and microdialysis probes",
      "implanted in the prostate. Single 7 mg/kg IV bolus of levofloxacin",
      "administered via the femoral vein."
    ),
    dose_range     = "7 mg/kg single IV bolus (~1.75-2.45 mg per rat); typical individual dose 2.198 mg",
    regions        = "Porto Alegre, Brazil (UFRGS)",
    notes          = paste(
      "Demographics from Hurtado 2014 Materials and Methods (Animal",
      "experiments). Plasma samples collected at 0.083, 0.25, 0.5, 0.75,",
      "1, 1.5, 2, 4, 6, 8, 12 h post-dose; prostate microdialysate",
      "collected at 20-min intervals up to 12 h. Plasma protein binding",
      "fb_plasma = 45.5 +/- 9.4% (fu_plasma = 0.545). Free prostate",
      "concentrations corrected by in vivo no-net-flux probe recovery",
      "RR_NNF = 17.5%."
    )
  )

  ini({
    # --------------------------------------------------------------------
    # Structural parameters. Final NONMEM 6 (ADVAN6 TRANS1) point estimates
    # from Hurtado 2014 Table 2. The three-compartment model has:
    #   - central (X1) total plasma
    #   - peripheral1 (X2) other tissues
    #   - effect (X3) prostate ISF (free drug), apparent volume V3* =
    #     V_prostate / fu_prostate (un-deconvolvable because microdialysis
    #     measures unbound concentrations only).
    # Central <-> peripheral1 transport is symmetric (single Q derivable
    # from V1, k12, k21), so the macro-constants CL / Vc / Q / Vp are used.
    # Central <-> effect (prostate) is asymmetric: the paper reports
    # separate k13 (uptake) and k31 (linear efflux), and there is an
    # additional saturable (Vmax, kM) efflux from prostate to central
    # consistent with active transporter involvement. The asymmetric arm
    # cannot be re-expressed under a single Q, so lk13 and lk31 are kept
    # as primary log-transformed parameters (precedent: Clinckers 2008
    # MHD rat lk23/lk32). The deviation is listed in the vignette
    # Assumptions and deviations.
    # --------------------------------------------------------------------
    lcl <- log(0.22)
    label("Total plasma clearance CL (L/h)")                                   # Hurtado 2014 Table 2: CL = 0.22 L/h

    lvc <- log(0.38)
    label("Central volume of distribution V1 (L)")                             # Hurtado 2014 Table 2: V1 = 0.38 L

    lq <- log(2.27 * 0.38)
    label("Inter-compartmental clearance Q (L/h) [Q = k12 * V1]")               # Hurtado 2014 Table 2: k12 = 2.27 /h, V1 = 0.38 L -> Q = 0.8626 L/h

    lvp <- log(2.27 * 0.38 / 1.44)
    label("Peripheral volume of distribution Vp (L) [Vp = k12 * V1 / k21]")     # Hurtado 2014 Table 2: k12 = 2.27 /h, V1 = 0.38 L, k21 = 1.44 /h -> Vp = 0.5990 L

    lk13 <- log(0.69)
    label("First-order uptake rate central -> prostate k13 (1/h)")             # Hurtado 2014 Table 2: k13 = 0.69 /h

    lk31 <- log(3.67)
    label("First-order linear efflux rate prostate -> central k31 (1/h)")      # Hurtado 2014 Table 2: k31 = 3.67 /h

    lveff <- log(0.05)
    label("Apparent prostate volume V3* = V_prostate / fu_prostate (L)")       # Hurtado 2014 Table 2: V3/fu_prostate = 0.05 L

    lvmax <- log(7.19e-3)
    label("Maximum saturable efflux velocity Vmax (mg/h)")                     # Hurtado 2014 Table 2: Vmax = 7.19 ug/h = 7.19e-3 mg/h (converted to internal mg units)

    lkm <- log(0.35)
    label("Michaelis-Menten constant kM (ug/mL = mg/L)")                       # Hurtado 2014 Table 2: kM = 0.35 ug/mL (numerically equal to 0.35 mg/L)

    # --------------------------------------------------------------------
    # Inter-individual variability. Hurtado 2014 Methods reports lognormal
    # IIV (P_i = theta * exp(eta_i)). Table 2 reports IIV in %CV for the
    # four parameters where IIV was retained (CL, V1, Vmax, kM); IIV on
    # k12, k21, k13, k31, V3*, Q, Vp was estimated as negligible and
    # dropped. Convert %CV to variance via omega^2 = log(1 + CV^2):
    #   V1   21.0 %CV -> omega^2 = log(1 + 0.210^2) = 0.04316
    #   CL   36.7 %CV -> omega^2 = log(1 + 0.367^2) = 0.12640
    #   Vmax 41.6 %CV -> omega^2 = log(1 + 0.416^2) = 0.15967
    #   kM   76.0 %CV -> omega^2 = log(1 + 0.760^2) = 0.45593
    # --------------------------------------------------------------------
    etalvc   ~ 0.04316                                                          # Hurtado 2014 Table 2: omega^2(V1)   = 21.0 %CV -> log(1 + 0.210^2)
    etalcl   ~ 0.12640                                                          # Hurtado 2014 Table 2: omega^2(CL)   = 36.7 %CV -> log(1 + 0.367^2)
    etalvmax ~ 0.15967                                                          # Hurtado 2014 Table 2: omega^2(Vmax) = 41.6 %CV -> log(1 + 0.416^2)
    etalkm   ~ 0.45593                                                          # Hurtado 2014 Table 2: omega^2(kM)   = 76.0 %CV -> log(1 + 0.760^2)

    # --------------------------------------------------------------------
    # Residual error. Hurtado 2014 Methods uses separate combined
    # proportional + additive error models for plasma (total) and free
    # prostate concentrations. Table 2 reports the proportional terms in
    # %CV and the additive terms in ug/mL (concentration units), i.e. the
    # values are residual SDs on the linear scale, not variances.
    #   sigma1 plasma prop   10.2 %CV  -> propSd          = 0.102
    #   sigma2 plasma add    0.085 ug/mL -> addSd         = 0.085
    #   sigma3 tissue prop   15.2 %CV  -> propSd_Cprostate = 0.152
    #   sigma4 tissue add    0.015 ug/mL -> addSd_Cprostate = 0.015
    # --------------------------------------------------------------------
    propSd <- 0.102
    label("Plasma proportional residual SD (fraction)")                        # Hurtado 2014 Table 2: sigma1 plasma proportional = 10.2 %CV

    addSd <- 0.085
    label("Plasma additive residual SD (ug/mL)")                               # Hurtado 2014 Table 2: sigma2 plasma additive = 0.085 ug/mL

    propSd_Cprostate <- 0.152
    label("Prostate ISF proportional residual SD (fraction)")                  # Hurtado 2014 Table 2: sigma3 tissue proportional = 15.2 %CV

    addSd_Cprostate <- 0.015
    label("Prostate ISF additive residual SD (ug/mL)")                         # Hurtado 2014 Table 2: sigma4 tissue additive = 0.015 ug/mL
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual parameters (lognormal IIV on CL, V1, Vmax, kM; no IIV
    #    on Q, Vp, k13, k31, V3* per Table 2 footnote).
    # ------------------------------------------------------------------
    cl   <- exp(lcl + etalcl)
    vc   <- exp(lvc + etalvc)
    q    <- exp(lq)
    vp   <- exp(lvp)
    k13  <- exp(lk13)
    k31  <- exp(lk31)
    veff <- exp(lveff)
    vmax <- exp(lvmax + etalvmax)
    km   <- exp(lkm + etalkm)

    # ------------------------------------------------------------------
    # 2. Micro-constants for the central <-> peripheral1 disposition.
    # ------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ------------------------------------------------------------------
    # 3. Prostate ISF free concentration (amount in mg, veff in L ->
    #    mg/L = ug/mL, matching the paper's reported units).
    # ------------------------------------------------------------------
    Cprostate <- effect / veff

    # ------------------------------------------------------------------
    # 4. Saturable (Michaelis-Menten) efflux from prostate to central.
    #    vmax is in mg/h (converted from the paper's 7.19 ug/h);
    #    Cprostate and km are in ug/mL = mg/L; the product is mg/h.
    # ------------------------------------------------------------------
    mm_efflux <- vmax * Cprostate / (km + Cprostate)

    # ------------------------------------------------------------------
    # 5. ODE system (Hurtado 2014 equations 3a-3c). The 7 mg/kg IV bolus
    #    goes directly into central (no depot).
    # ------------------------------------------------------------------
    d/dt(central)     <- -kel * central -
                          k12 * central + k21 * peripheral1 -
                          k13 * central + k31 * effect      +
                          mm_efflux
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(effect)      <-  k13 * central - k31 * effect      -
                          mm_efflux

    # ------------------------------------------------------------------
    # 6. Observations (plasma total + free prostate ISF) and combined
    #    proportional + additive residual error per output.
    # ------------------------------------------------------------------
    Cc <- central / vc

    Cc        ~ prop(propSd)           + add(addSd)
    Cprostate ~ prop(propSd_Cprostate) + add(addSd_Cprostate)
  })
}
