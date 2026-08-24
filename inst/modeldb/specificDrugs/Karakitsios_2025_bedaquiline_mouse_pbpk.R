Karakitsios_2025_bedaquiline_mouse_pbpk <- function() {
  description <- paste(
    "PBPK (middle-out, permeability-limited lung; Monolix Suite 2023R1).",
    "Preclinical (mouse). Bedaquiline distribution into healthy/uninvolved lung and into",
    "tuberculosis-infected lung tissue in mice. An empirical three-compartment oral plasma",
    "model drives a permeability-limited lung model split into pulmonary blood, extracellular",
    "water, intracellular water and lysosomes, with pH-dependent ion trapping in the lysosomes.",
    "The infected-lung extension adds a cellular-lesion (macrophage) compartment and a",
    "caseous-granuloma arm in which an outer caseum of foamy macrophages is in instantaneous",
    "equilibrium with a six-compartment catenary chain representing 300-micrometre rings of the",
    "necrotic core, with anomalous (fractal) diffusion described by a time-dependent rate",
    "coefficient. The lung states hold CONCENTRATIONS, not amounts. This is the species in which",
    "every lung, lesion and caseum parameter was estimated.",
    sep = " "
  )
  reference <- paste(
    "Karakitsios E, Della Pasqua O, Dokoumetzidis A.",
    "Extrapolation of lung pharmacokinetics of bedaquiline across species using",
    "physiologically-based pharmacokinetic modelling.",
    "Br J Clin Pharmacol. 2025;91(11):3167-3178.",
    "doi:10.1002/bcp.70163.",
    sep = " "
  )
  vignette <- "Karakitsios_2025_bedaquiline_lung_pbpk"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # caseum1-caseum6 (the paper's cas1-cas6) are its spatial binning of the
  # necrotic caseous core into 300-micrometre rings measured outward-to-inward
  # from the outer caseum edge (Supporting Information, "PBPK models for the
  # infected lung"; Walter 2023 laser-capture ring dissection). `caseum<n>` is
  # a canonical chain prefix -- see inst/references/compartment-names.md.

  compartmentData <- list(
    depot       = list(analyte = "bedaquiline", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "bedaquiline", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "bedaquiline", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "bedaquiline", units = "mg", specimen = "plasma", verified = TRUE),
    lung_ew     = list(analyte = "bedaquiline", units = "ug/mL", specimen = "tissue", verified = TRUE),
    lung_iw     = list(analyte = "bedaquiline", units = "ug/mL", specimen = "tissue", verified = TRUE),
    lesion      = list(analyte = "bedaquiline", units = "ug/mL", specimen = "tissue", verified = TRUE),
    caseum1     = list(analyte = "bedaquiline", units = "ug/mL", specimen = "tissue", verified = TRUE),
    caseum2     = list(analyte = "bedaquiline", units = "ug/mL", specimen = "tissue", verified = TRUE),
    caseum3     = list(analyte = "bedaquiline", units = "ug/mL", specimen = "tissue", verified = TRUE),
    caseum4     = list(analyte = "bedaquiline", units = "ug/mL", specimen = "tissue", verified = TRUE),
    caseum5     = list(analyte = "bedaquiline", units = "ug/mL", specimen = "tissue", verified = TRUE),
    caseum6     = list(analyte = "bedaquiline", units = "ug/mL", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Not a model covariate. Every plasma and lung parameter below is the ABSOLUTE value for",
        "a typical 0.025 kg mouse: Table S5 reports the plasma model in absolute L and L/h",
        "alongside a weight-normalised column, and Table S3 tabulates lung blood flow, the four",
        "lung sub-volumes and the lung-cell surface area as fixed per-species quantities that are",
        "not weight-scalable within a species. Dose in mg for a 0.025 kg animal (e.g. 25 mg/kg =",
        "0.625 mg). Table S3 footnote a gives the Simcyp reference weight used for mice (0.025 kg).",
        sep = " "
      )
    )
  )

  population <- list(
    species        = "mouse (BALB/c, C3HeB/FeJ 'Kramnik', and Swiss SPF (CD1))",
    n_subjects     = NA_integer_,
    n_studies      = 3L,
    sex_female_pct = NA_real_,
    age_range      = "8-10 weeks for the BALB/c and C3HeB/FeJ infection protocols",
    weight_median  = "0.025 kg (Simcyp reference mouse; Table S3 footnote a)",
    disease_state  = paste(
      "Healthy/uninvolved lung was characterised in uninfected and infected animals across three",
      "strains. Infected-lung data come from Mycobacterium tuberculosis-infected mice: BALB/c",
      "develop type III cellular non-necrotising lesions only, whereas C3HeB/FeJ (Kramnik) mice",
      "develop type I caseous necrotic granulomas with a foamy-macrophage cellular rim around a",
      "necrotic caseous core.",
      sep = " "
    ),
    dose_range     = paste(
      "Oral gavage. 25 mg/kg single dose (BALB/c and Kramnik, Irwin 2016); 25 mg/kg/day single and",
      "multiple doses (Kramnik, Walter 2023; up to 17 daily doses); 6.25 and 25 mg/kg single and",
      "multiple doses (Swiss SPF (CD1), Janssen R&D). Table S1.",
      sep = " "
    ),
    notes          = paste(
      "Mean concentrations were modelled rather than individual data because most protocols used",
      "only about 3 animals per sampling time (Methods 2.1). Plasma data were Janssen R&D rich",
      "sampling; lung data pooled Kramnik, BALB/c and Swiss SPF strains with a random effect on",
      "cellular permeability (Methods 2.2). Strain-specific empirical Bayes permeabilities are",
      "reported in Methods 2.2 as 29.6 cm/h (Kramnik), 43.5-68.2 cm/h (Swiss SPF) and 79.1 cm/h",
      "(BALB/c) around the 55.7 cm/h typical value of Table S8. IMPORTANT: the infected-lung",
      "parameters of Tables S11 and S13 were estimated against the BALB/c and Kramnik plasma",
      "models of Irwin 2016, which are not tabulated in this paper; this file uses the Table S5",
      "mouse plasma model, which yields a lower plasma exposure. See the vignette Errata.",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # 1. Empirical plasma model (Table S5). Absolute values for a typical
    #    0.025 kg mouse. Compartment 2 is the deep peripheral space and
    #    compartment 3 the shallow one; the paper's V2/Q2 map to vp/q and
    #    its V3/Q3 to vp2/q2, preserving the published numbering.
    # =====================================================================
    lka     <- log(0.13)     ; label("Absorption rate constant ka (1/h)")                        # Table S5 ka = 0.13 (RSE 12.1%)
    lcl     <- log(0.0089)   ; label("Clearance CL (L/h, typical 0.025 kg mouse)")               # Table S5 Cl = 0.0089 L/h (RSE 14.9%); 0.356 L/h/kg
    lvc     <- log(0.00033)  ; label("Central volume V1 (L, typical 0.025 kg mouse)")            # Table S5 V1 = 0.00033 L (RSE 31.7%); 0.0132 L/kg
    lq      <- log(0.0024)   ; label("Inter-compartmental clearance Q2 (L/h, typical mouse)")    # Table S5 Q2 = 0.0024 L/h (RSE 22.7%)
    lvp     <- log(0.14)     ; label("First peripheral volume V2 (L, typical mouse)")            # Table S5 V2 = 0.14 L (RSE 16.4%)
    lq2     <- log(0.006)    ; label("Inter-compartmental clearance Q3 (L/h, typical mouse)")    # Table S5 Q3 = 0.006 L/h (RSE 43.0%)
    lvp2    <- log(0.0094)   ; label("Second peripheral volume V3 (L, typical mouse)")           # Table S5 V3 = 0.0094 L (RSE 35.2%)
    lfdepot <- log(0.2)      ; label("Oral bioavailability F (unitless)")                        # Table S5 F = 0.2 (RSE 22.4%)

    # =====================================================================
    # 2. Drug-related parameters, independent of species (Table S2).
    #    Bedaquiline is a diprotic BASE; both pKa values are basic and are
    #    labelled as such on the structure in Figure S8.
    # =====================================================================
    pka1 <- fixed(8.5) ; label("First (strong) basic pKa (unitless)")   # Table S2 pKa1 = 8.5 (base), Janssen R&D data
    pka2 <- fixed(1.1) ; label("Second (weak) basic pKa (unitless)")    # Table S2 pKa2 = 1.1 (base), ACD predicted

    # =====================================================================
    # 3. Species-related lung physiology for the mouse (Table S3).
    #    Flows in mL/h, volumes in mL, surface area in cm^2 -- the units the
    #    permeability-limited lung equations (S1) and (S4) are written in.
    # =====================================================================
    qlung        <- fixed(840.6)  ; label("Lung blood flow rate Qlung (mL/h)")                     # Table S3 mouse blood flow rate = 840.6 mL/h
    sa_cells     <- fixed(500)    ; label("Surface area of lung cells SAcells (cm^2)")             # Table S3 mouse surface area lung cells = 500 cm^2
    v_pb         <- fixed(0.0314) ; label("Volume of pulmonary blood VPB (mL)")                    # Table S3 mouse volume of pulmonary blood = 0.0314 mL
    v_ew         <- fixed(0.0192) ; label("Volume of lung extracellular water VEW (mL)")           # Table S3 mouse volume of EW = 0.0192 mL
    v_lysosome   <- fixed(0.0015) ; label("Volume of lung lysosomes VLYS (mL)")                    # Table S3 mouse volume of lysosomes = 0.0015 mL (1% of total lung volume)
    v_iw         <- fixed(0.098)  ; label("Volume of lung intracellular water VIW (mL)")           # Table S3 mouse volume of IW = 0.098 mL
    ph_plasma    <- fixed(7.39)   ; label("Plasma pH (unitless)")                                  # Table S3 mouse pH_plasma = 7.39 (Iversen 2012)
    ph_ew        <- fixed(7.4)    ; label("Lung extracellular-water pH (unitless)")                # Table S3 mouse pH_EW = 7.4 (Simcyp 2023)
    ph_iw        <- fixed(7.36)   ; label("Lung intracellular-water pH (unitless)")                # Table S3 mouse pH_IW = 7.36 (Simcyp 2023)
    ph_lysosome  <- fixed(5.0)    ; label("Lysosomal pH (unitless)")                               # Table S3 mouse pH_Lysosomes = 5.0; Methods 2.2 "lysosomes had a constant pH of 5.0 across species"
    # Methods 2.3 reports pH 5.84 for foamy macrophages (Kempker 2017) and adds
    # that mice and humans "share the same pH in cellular lesions as well as in
    # the cellular rim of caseous granulomas". No separate pH is given for the
    # macrophages of cellular lesions, so the same 5.84 is applied to both; this
    # affects only the derived UNBOUND observables, never the ODE states or any
    # total concentration. See vignette Errata.
    ph_lesion_iw <- fixed(5.84)   ; label("Cytosolic pH of macrophages and foamy macrophages (unitless)")  # Methods 2.3: "A pH of 5.84 was used for foamy macrophages" (Kempker 2017)

    # =====================================================================
    # 4. Drug- AND species-related parameters (Table S4).
    # =====================================================================
    bp        <- fixed(1.0)     ; label("Blood-to-plasma concentration ratio B/P (unitless)")      # Table S4 mouse blood to plasma ratio = 1.0 (Janssen R&D data)
    fu_plasma <- fixed(7.2e-4)  ; label("Fraction unbound in plasma fupls (unitless)")             # Table S4 mouse fupls = 7.2E-4, derived from albumin binding via equations (S5)-(S6)
    fu_ew     <- fixed(1.1e-3)  ; label("Fraction unbound in lung extracellular water fuEW")       # Table S4 mouse fuEW = 1.1E-3, from equation (S7)
    fu_iw     <- fixed(1.0)     ; label("Fraction unbound in lung intracellular water fuIW")       # Supporting Information: "a fixed value of 1 was used for fuIW"; sensitivity analysis in Table S9

    # =====================================================================
    # 5. Lung parameters estimated on mouse lung PK profiles (Table S8).
    # =====================================================================
    lperm_cells <- log(55.7)  ; label("Permeability of lung cells PermCells (cm/h, log scale)")    # Table S8 Permcells = 5.57E+1 cm/h (RSE 18.9%); Table S10 control row = 55.7
    # KpLYS:IW/fuIW was fitted in mice as 3.48E+6 (Table S8). Equation (S8)
    # converts it to the species-invariant product fuLYS*fniLYS, which the
    # paper holds constant across species; equation (S9) then recovers
    # KpLYS:IW/fuIW for any species from its own pH_IW. Stored as that product
    # so the same constant appears in every species file:
    #   fuLYS*fniLYS = fniIW(mouse) / 3.48e6
    #                = (1/((1+10^(8.5-7.36))*(1+10^(1.1-7.36)))) / 3.48e6
    #                = 0.0675502 / 3.48e6 = 1.9411e-8
    # In this mouse file model() therefore recovers KpLYS:IW/fuIW = 3.48e6 exactly.
    lfulys_fnilys <- fixed(log(1.9411e-8)) ; label("Product fuLYS*fniLYS in lysosomes (unitless, log scale)")  # derived from Table S8 KpLYS:IW/fuIW = 3.48E+6 via equations (S8)-(S9)

    # =====================================================================
    # 6. Infected lung -- cellular lesion (macrophages), fitted on BALB/c
    #    mouse data (Table S11). fuIW(mphi) and fniIW(mphi) were both fixed
    #    to 1, so the fitted ratio Kpmphi/(fuIW*fniIW) IS Kpmphi.
    # =====================================================================
    lrate_macrophage <- log(143e3)  ; label("Macrophage permeation rate constant rate_mphi (1/h, log scale)")  # Table S11 rate_mphi = 143E+3 1/h (RSE 7.84%)
    lkp_macrophage   <- log(568e4)  ; label("Lysosome:cytosol partition in macrophages Kpmphi (unitless, log)") # Table S11 Kpmphi/(fuIW(mphi)*fniIW(mphi)) = 568E+4 (RSE 3.58%), with fuIW = fniIW = 1
    fu_iw_macrophage  <- fixed(1.0) ; label("Unbound fraction in macrophage cytosol fuIW(mphi)")                # Supporting Information: "a fixed value of 1 was used for fuIW(mphi), fniIW(mphi), fuIW(fmphi) and fniIW(fmphi)"
    fni_iw_macrophage <- fixed(1.0) ; label("Non-ionised fraction in macrophage cytosol fniIW(mphi)")           # Supporting Information, same sentence

    # =====================================================================
    # 7. Infected lung -- caseous granuloma, fitted on Kramnik mouse data
    #    (Table S13).
    # =====================================================================
    lrate_foamy <- log(5.78e3) ; label("Foamy-macrophage permeation rate constant rate_fmphi (1/h, log)")  # Table S13 rate_fmphi = 5.78E+3 1/h (RSE 25.6%)
    lkp_foamy   <- log(1.44e6) ; label("Lysosome:cytosol partition in foamy macrophages Kpfmphi (unitless, log)") # Table S13 Kpfmphi/(fuIW(fmphi)*fniIW(fmphi)) = 1.44E+6 (RSE 27.4%), with fuIW = fniIW = 1
    lkp_caseum  <- log(0.71)   ; label("Caseum-1 : outer-caseum partition KpCaseum (unitless, log scale)")  # Table S13 KpCaseum = 0.71 (RSE 29.5%)
    lkcas       <- log(0.25)   ; label("Anomalous-diffusion rate constant kcas (h^(h-1), log scale)")       # Table S13 kcas = 0.25 (RSE 18.9%)
    h_caseum    <- 0.75        ; label("Heterogeneity exponent h of the fractal caseum diffusion (unitless)") # Table S13 h = 0.75 (RSE 6.69%); equation (S15) requires 0 <= h <= 1

    # Equation (S15) writes kcaseum(t) = kcas / (t + tau)^h, where tau is
    # "fixed to a very small positive constant to avoid infinities for t = 0".
    # The paper does not report its value; 1e-3 h is used here. Because h < 1
    # the time integral of kcaseum converges, so the choice is weakly
    # influential: over 0-24 h, tau = 1e-3 and tau = 1e-6 differ by ~7% in
    # the cumulative diffusion, and less thereafter. See vignette Errata.
    tau_caseum <- fixed(1e-3) ; label("Onset-time regulariser tau of equation (S15) (h)")  # equation (S15): "Critical time tau was fixed to a very small positive constant"; value not reported

    # Fixed geometry of the infected lung.
    flys_cell   <- fixed(0.078) ; label("Lysosomal contribution to cellular volume fLys(AM) (unitless)")  # Supporting Information: 7.8% average in rat alveolar macrophages (Ufuk 2017), assumed equal across species
    ratio_cas1  <- fixed(1.67)  ; label("Volume ratio outer caseum : caseum-1 (ratio1, unitless)")        # Supporting Information: "A value of 1.67 was used for ratio1", from a representative Kramnik LCM image (Walter 2023); same value used in humans
    ratio_inner <- fixed(1.0)   ; label("Volume ratio caseum(n-1) : caseum(n) for n = 2..6 (unitless)")   # Supporting Information: "a value equal to 1 was used for all other ratios, namely ratio2 to ratio6"

    # =====================================================================
    # 8. Random effects.
    # =====================================================================
    etalfdepot   ~ 0.3136   # Table S5 omega_F = 0.56 (Monolix SD of the log-scale random effect); variance = 0.56^2
    etalka       ~ 0.0529   # Table S5 omega_ka = 0.23; variance = 0.23^2
    etalperm_cells ~ 0.1521 # Table S8 omega_Permcells = 0.39; INTER-STRAIN (not inter-individual) variability, see Table S8 footnote

    # =====================================================================
    # 9. Residual error. Each layer was a separate Monolix estimation, so the
    #    four proportional error terms come from four different tables.
    # =====================================================================
    propSd         <- 0.2   ; label("Plasma proportional residual error (fraction)")           # Table S5 b = 0.2 (RSE 9.53%)
    propSd_c_lung    <- 0.24  ; label("Total lung proportional residual error (fraction)")       # Table S8 b = 0.24 (RSE 10.3%)
    propSd_c_lesion  <- 0.096 ; label("Cellular-lesion proportional residual error (fraction)")  # Table S11 b = 0.096 (RSE 25.0%)
    propSd_c_outer_caseum  <- 0.48  ; label("Caseum proportional residual error (fraction)")           # Table S13 b = 0.48 (RSE 14.4%)
  })

  model({
    # -------------------------------------------------------------------
    # 1. Individual plasma parameters
    # -------------------------------------------------------------------
    ka  <- exp(lka + etalka)
    cl  <- exp(lcl)
    vc  <- exp(lvc)
    q   <- exp(lq)
    vp  <- exp(lvp)
    q2  <- exp(lq2)
    vp2 <- exp(lvp2)

    perm_cells     <- exp(lperm_cells + etalperm_cells)
    rate_macrophage <- exp(lrate_macrophage)
    kp_macrophage  <- exp(lkp_macrophage)
    rate_foamy     <- exp(lrate_foamy)
    kp_foamy       <- exp(lkp_foamy)
    kp_caseum      <- exp(lkp_caseum)
    kcas           <- exp(lkcas)

    # -------------------------------------------------------------------
    # 2. Ionisation. Supporting Information gives the general four-term
    #    Henderson-Hasselbalch factor
    #      W = 1 + 10^(pKa1-pH) + 10^(pH-pKa2) + 10^(pKa1+pKa2-2pH)
    #    whose four terms are exactly the expansion of the DIPROTIC BASE
    #    product (1+10^(pKa1-pH)) * (1+10^(pKa2-pH)); the printed third term
    #    is the monoprotic-ACID term and carries a reversed exponent. Table S2
    #    and Figure S8 both declare bedaquiline's two pKa values to be BASIC,
    #    so the product form is used here. Taking the third term literally
    #    inverts the lysosomal pH gradient that drives the whole model and
    #    reproduces the paper's own rat prediction (Figure S4) 80-fold low;
    #    the product form reproduces it within 5%. See vignette Errata.
    x_plasma <- (1 + 10^(pka1 - ph_plasma)) * (1 + 10^(pka2 - ph_plasma))
    x_ew     <- (1 + 10^(pka1 - ph_ew))     * (1 + 10^(pka2 - ph_ew))
    x_iw     <- (1 + 10^(pka1 - ph_iw))     * (1 + 10^(pka2 - ph_iw))

    fi_plasma <- 1 / x_plasma      # non-ionised fraction in plasma
    fi_ew     <- 1 / x_ew          # non-ionised fraction in lung extracellular water
    fi_iw     <- 1 / x_iw          # non-ionised fraction in lung intracellular water

    # Equation (S9): KpLYS:IW / fuIW = fniIW / (fuLYS * fniLYS), with the
    # denominator fitted in mice and held constant across species.
    kp_lysosome <- fi_iw / exp(lfulys_fnilys)

    # -------------------------------------------------------------------
    # 3. Derived lung constants
    # -------------------------------------------------------------------
    # KEW:B follows from the stated instantaneous PB/EW equilibrium of the
    # unbound-unionised drug (Supporting Information: "an instantaneous
    # equilibrium was assumed between PB and EW, with unbound-unionized drug
    # concentrations in these two compartments being equal"), i.e.
    # fupls*CPB*fi_plasma = fuEW*CEW*fi_ew, so CEW/CPB = (fupls/fuEW)*(x_ew/x_plasma).
    kewb <- (fu_plasma / fu_ew) * (x_ew / x_plasma)

    ps      <- perm_cells * sa_cells           # cm^3/h = mL/h
    v_ew_eff <- (v_ew + v_pb / kewb) / fu_ew   # left-hand-side volume term of equation (S1)
    v_iw_eff <- v_iw + kp_lysosome * v_lysosome # left-hand-side volume term of equation (S4)
    f_iw    <- v_iw / (v_pb + v_ew + v_iw + v_lysosome)
    r_lys   <- flys_cell / f_iw                # the 0.078/fIW term of equations (S11)-(S12)

    cblood <- (central / vc) * bp
    jperm  <- ps * (lung_ew * fi_ew - fu_iw * lung_iw * fi_iw)

    # -------------------------------------------------------------------
    # 4. Plasma ODEs
    # -------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - cl * central / vc -
                          q * central / vc + q * peripheral1 / vp -
                          q2 * central / vc + q2 * peripheral2 / vp2
    d/dt(peripheral1) <-  q * central / vc - q * peripheral1 / vp
    d/dt(peripheral2) <-  q2 * central / vc - q2 * peripheral2 / vp2

    f(depot) <- exp(lfdepot + etalfdepot)

    # -------------------------------------------------------------------
    # 5. Healthy/uninvolved lung -- equations (S1) and (S4). lung_ew holds
    #    the UNBOUND extracellular-water concentration; lung_iw holds the
    #    TOTAL intracellular-water concentration.
    # -------------------------------------------------------------------
    d/dt(lung_ew) <- (qlung * cblood - (qlung / fu_ew) * (lung_ew / kewb) - jperm) / v_ew_eff
    d/dt(lung_iw) <- jperm / v_iw_eff

    # Unbound-unionised extracellular driver shared by both infected-lung arms.
    cu_ew_ni <- lung_ew * fi_ew

    # -------------------------------------------------------------------
    # 6. Infected lung -- cellular lesion, equation (S11).
    # -------------------------------------------------------------------
    denom_macrophage <- 1 / (r_lys + 1) + kp_macrophage * r_lys / (r_lys + 1)
    d/dt(lesion) <- rate_macrophage *
      (cu_ew_ni - fu_iw_macrophage * fni_iw_macrophage / denom_macrophage * lesion)

    # -------------------------------------------------------------------
    # 7. Infected lung -- caseous granuloma, equations (S12)-(S15).
    #    The outer caseum (foamy macrophages) is in instantaneous equilibrium
    #    with caseum1, so it is not a separate state: its concentration is
    #    caseum1 / KpCaseum, and the effective volume factor (1 + ratio1/KpCaseum)
    #    divides the caseum1 derivative.
    # -------------------------------------------------------------------
    kcaseum   <- kcas / (t + tau_caseum)^h_caseum
    x_foamy   <- (r_lys + 1) / (1 + kp_foamy * r_lys)
    vfac_cas1 <- 1 + ratio_cas1 / kp_caseum

    d/dt(caseum1) <- (rate_foamy * ratio_cas1 *
                        (cu_ew_ni - x_foamy * caseum1 / kp_caseum) -
                      kcaseum * caseum1 + kcaseum * caseum2 / ratio_inner) / vfac_cas1
    d/dt(caseum2) <- kcaseum * caseum1 * ratio_inner - 2 * kcaseum * caseum2 +
                     kcaseum * caseum3 / ratio_inner
    d/dt(caseum3) <- kcaseum * caseum2 * ratio_inner - 2 * kcaseum * caseum3 +
                     kcaseum * caseum4 / ratio_inner
    d/dt(caseum4) <- kcaseum * caseum3 * ratio_inner - 2 * kcaseum * caseum4 +
                     kcaseum * caseum5 / ratio_inner
    d/dt(caseum5) <- kcaseum * caseum4 * ratio_inner - 2 * kcaseum * caseum5 +
                     kcaseum * caseum6 / ratio_inner
    d/dt(caseum6) <- kcaseum * caseum5 * ratio_inner - kcaseum * caseum6

    # -------------------------------------------------------------------
    # 8. Observations.
    # -------------------------------------------------------------------
    Cc <- central / vc                                   # plasma bedaquiline

    c_lung_ew  <- lung_ew / fu_ew                        # total EW concentration
    c_lung_pb  <- c_lung_ew / kewb                       # total pulmonary-blood concentration
    c_lysosome <- kp_lysosome * lung_iw                  # total lysosomal concentration
    # Total lung homogenate = amount in the whole tissue / total lung volume.
    c_lung <- (c_lung_pb * v_pb + c_lung_ew * v_ew +
               lung_iw * v_iw + c_lysosome * v_lysosome) /
              (v_pb + v_ew + v_iw + v_lysosome)

    c_lesion <- lesion                                   # total macrophage (cellular lesion)
    c_outer_caseum <- caseum1 / kp_caseum                # total foamy-macrophage (outer caseum)

    # Equations (S16)-(S19). (S17) gives the TOTAL lysosomal concentration from
    # the whole-cell concentration, (S18) the unbound cytosolic concentration
    # and (S19) the unbound lysosomal concentration.
    #
    # Tables S11 and S13 report the RATIO Kp/(fuIW*fniIW), estimated with fuIW
    # and fniIW fixed to 1, so `kp_macrophage` and `kp_foamy` above are that
    # ratio and not Kp itself. Equations (S16)-(S19) are statements about the
    # true Kp = CLYS/CIW, so the ratio is un-normalised here with the
    # prediction-time fniIW evaluated at the cellular-lesion pH of 5.84
    # (Methods 2.3). Skipping this step leaves the unbound lysosomal
    # concentrations about 460-fold low and contradicts the paper's own
    # conclusion (Results 3.3, Figure 4) that they exceed the MIC50 and MIC90;
    # see vignette Errata. The ODE states are unaffected -- equation (S11) is
    # written and was fitted in terms of the ratio, and is left as published.
    x_lesion_iw <- (1 + 10^(pka1 - ph_lesion_iw)) * (1 + 10^(pka2 - ph_lesion_iw))
    x_lysosome  <- (1 + 10^(pka1 - ph_lysosome))  * (1 + 10^(pka2 - ph_lysosome))
    fni_lesion_iw <- 1 / x_lesion_iw
    fni_lysosome  <- 1 / x_lysosome

    kp_macrophage_true <- kp_macrophage * fu_iw_macrophage * fni_lesion_iw
    kp_foamy_true      <- kp_foamy      * fu_iw_macrophage * fni_lesion_iw

    c_lesion_lys  <- c_lesion * (r_lys + 1) / (1 / kp_macrophage_true + r_lys)   # (S17), macrophages
    cu_lesion_iw  <- fu_iw_macrophage * c_lesion_lys / kp_macrophage_true        # (S18)
    cu_lesion_lys <- cu_lesion_iw * (fni_lesion_iw / fni_lysosome)               # (S19)

    c_caseum_lys  <- c_outer_caseum * (r_lys + 1) / (1 / kp_foamy_true + r_lys)  # (S17), foamy macrophages
    cu_caseum_iw  <- fu_iw_macrophage * c_caseum_lys / kp_foamy_true             # (S18)
    cu_caseum_lys <- cu_caseum_iw * (fni_lesion_iw / fni_lysosome)               # (S19)

    cu_plasma <- Cc * fu_plasma                                            # unbound plasma

    Cc             ~ prop(propSd)
    c_lung         ~ prop(propSd_c_lung)
    c_lesion       ~ prop(propSd_c_lesion)
    c_outer_caseum ~ prop(propSd_c_outer_caseum)
  })
}
