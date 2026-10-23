Byun_2024_efineptakin_alfa_ptmdd <- function() {
  description <- paste(
    "First-order (Pade) approximation of the total-QSSA target-mediated drug disposition model",
    "(pTMDD) for efineptakin alfa (rhIL-7-hyFc, GX-I7) given by intramuscular or subcutaneous",
    "injection to healthy adults, with FcRn-mediated recycling; case 2 'pTMDD' model of",
    "Byun 2024, the paper's novel approximation. Same structure as",
    "modellib('Byun_2024_efineptakin_alfa_qtmdd') but the bound drug is the first-order",
    "approximation of the QSS quadratic with no square root: FcRn_tot * C / (KSS1 + C +",
    "FcRn_tot) in the endosomal space and Rtot * C / (KSS2 + C + Rtot) in central, free drug",
    "being the total minus the bound. Doses are pmol of GX-I7 with bioavailability 2 (IL-7",
    "basis); the observation is the free central concentration plus the observed baseline",
    "IL-7 (BL_IL7)."
  )
  reference <- paste(
    "Byun JH, Jeon HS, Yun H-y, Kim JK. Validity conditions of approximations for a",
    "target-mediated drug disposition model: A novel first-order approximation and its",
    "comparison to other approximations. PLoS Comput Biol. 2024;20(4):e1012066.",
    "doi:10.1371/journal.pcbi.1012066. Case 2 pTMDD model: Table 4 (estimates) and",
    "S1 Text (h) (NONMEM control stream). The clinical study is Lee SW et al. hIL-7-hyFc, A",
    "Long-Acting IL-7, Increased Absolute Lymphocyte Count in Healthy Subjects. Clin Transl",
    "Sci. 2020;13(6):1161-1169. doi:10.1111/cts.12800."
  )
  vignette <- "Byun_2024_tmdd_approximations"
  units <- list(time = "h", dosing = "pmol", concentration = "pmol/L")

  compartmentData <- list(
    depot = list(
      analyte = "efineptakin alfa (IL-7 basis)",
      units = "pmol",
      specimen = "administration site",
      verified = TRUE
    ),
    depot2 = list(
      analyte = "efineptakin alfa (IL-7 basis)",
      units = "pmol",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "efineptakin alfa (total: free + IL-7R-bound, IL-7 basis)",
      units = "pmol",
      specimen = "serum",
      verified = TRUE
    ),
    isf = list(analyte = "efineptakin alfa (IL-7 basis)", units = "pmol", specimen = "tissue", verified = TRUE),
    endosome = list(
      analyte = "efineptakin alfa (total: free + FcRn-bound, IL-7 basis)",
      units = "pmol",
      specimen = "endosome",
      verified = TRUE
    )
  )

  covariateData <- list(
    ROUTE_SC = list(
      description = "Route of the subject's single dose: 1 = subcutaneous (cohorts SC 20 and SC 60 ug/kg), 0 = intramuscular (cohort IM 60 ug/kg)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (IM)",
      notes = "Per-subject cohort indicator. Selects the volume of the interstitial/vascular absorption space: S1 Text (h) $PK 'VI = VISC; IF(COHORT.EQ.3) VI = VIIM'. SC doses go to depot2 and IM doses to depot.",
      source_name = "COHORT (1 = SC20, 2 = SC60, 3 = IM60); ROUTE_SC = 1 when COHORT is 1 or 2"
    ),
    BL_IL7 = list(
      description = "Observed baseline endogenous serum IL-7 concentration of the subject",
      units = "pmol/L",
      type = "continuous",
      reference_category = "n/a",
      notes = "Added to the predicted drug concentration: S1 Text (h) $ERROR 'IPRED = BASE + Cf ; add baseline IL7 conc'. Range 1.51-4.98 pmol/L in the S2 Table dataset.",
      source_name = "BASE"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 24,
    n_studies = 1,
    disease_state = "Healthy adult subjects",
    dose_range = "Single dose of efineptakin alfa: 20 ug/kg SC (n = 8), 60 ug/kg SC (n = 8), or 60 ug/kg IM (n = 8)",
    regions = "Republic of Korea",
    n_observations = 410,
    notes = paste(
      "Byun 2024 case 2 refits the healthy-subject study of Lee 2020 (Clin Transl Sci",
      "13:1161). Byun 2024 states 456 records; the S2 Table dataset holds 434 rows (24 dose",
      "records and 410 observations over 672 h). Byun 2024 does not report demographics. The",
      "dataset AMT is in pmol of GX-I7 (MW 104 kDa; the IL-7 part is 20.5 kDa) and the model",
      "converts to IL-7 pmol with F = 2 on both depots."
    )
  )

  ini({
    # Final estimates are Byun 2024 Table 4 (case 2, pTMDD column). VIIM and VISC
    # are not in Table 4 and are taken from the S1 Text (h) $THETA, which carries
    # the final estimates re-used as starting values.
    lcl <- log(5.27); label("Clearance of free drug from central (L/h)") # Table 4 case 2 pTMDD CL = 5.27; S1 Text (h) THETA(1)
    lvc <- log(2.07); label("Central volume of distribution (L)") # Table 4 case 2 pTMDD Vc = 2.07; S1 Text (h) THETA(2) 2.08
    lka_im <- log(1.24); label("Absorption rate from the intramuscular depot (1/h)") # Table 4 case 2 pTMDD Ka = 1.24 (im); S1 Text (h) THETA(3) KAIM
    lka_sc <- log(0.679); label("Absorption rate from the subcutaneous depot (1/h)") # Table 4 case 2 pTMDD Ka = 0.679 (sc); S1 Text (h) THETA(4) KASC
    lq <- log(1.3); label("Inter-compartmental clearance between central and the absorption space (L/h)") # Table 4 case 2 pTMDD Q = 1.3; S1 Text (h) THETA(5)
    lv_endosome <- log(35.9); label("Apparent volume of the endosomal FcRn space (L)") # Table 4 case 2 pTMDD Vd = 35.9; S1 Text (h) THETA(7) VD
    lkrec <- log(0.000619); label("Recycling rate of the drug-FcRn complex (1/h)") # Table 4 case 2 pTMDD Krecycle = 0.000619; S1 Text (h) THETA(8) 0.000621
    lv_isf_im <- log(358); label("Volume of the absorption space, intramuscular cohort (L)") # S1 Text (h) THETA(13) VIIM = 358 (not in Table 4)
    lv_isf_sc <- log(1020); label("Volume of the absorption space, subcutaneous cohorts (L)") # S1 Text (h) THETA(14) VISC = 1020 (not in Table 4)

    # Held constant at literature values (Table 4 asterisk, or hard-coded in the
    # S1 Text (h) $PK block).
    lkdeg <- fixed(log(0.642)); label("Degradation rate of free drug in the endosomal space (1/h)") # Table 4 case 2 Kdeg = 0.642 (literature); S1 Text (h) THETA(6) FIX
    lkint <- fixed(log(0.642)); label("Internalisation rate of the IL-7R-drug complex (1/h)") # Table 4 case 2 Kint = 0.642 (literature); S1 Text (h) THETA(9) FIX
    lc_il7r_t <- fixed(log(1060)); label("Total IL-7R concentration in central (pmol/L)") # Table 4 case 2 Rtot = 1060 (literature); S1 Text (h) THETA(10) FIX
    lkup <- fixed(log(0.00952)); label("Uptake rate from the absorption space into the endosomal space (1/h)") # Table 4 case 2 Kuptake = 0.00952; S1 Text (h) $PK Kup = 0.00952
    lkon1 <- fixed(log(4.068e-4)); label("Drug-FcRn association rate constant (L/pmol/h)") # S1 Text (h) $PK KON1 = 4.068/10000
    lkd1 <- fixed(log(29823)); label("Drug-FcRn equilibrium dissociation constant (pmol/L)") # S1 Text (h) $PK KD = 29823
    lc_fcrn_t <- fixed(log(4e7)); label("Total FcRn concentration in the endosomal space (pmol/L)") # S1 Text (h) $PK FcRn = 40000000
    fr <- fixed(0.715); label("Fraction of recycled drug-FcRn complex returned to central (fraction)") # S1 Text (h) $PK FR = 0.715
    lkon2 <- fixed(log(2.0196e-3)); label("Drug-IL-7R association rate constant (L/pmol/h)") # S1 Text (h) $PK KON2 = 20.196/10000
    lkoff2 <- fixed(log(59.76)); label("Drug-IL-7R dissociation rate constant (1/h)") # S1 Text (h) $PK KOFF2 = 59.76

    etalq ~ 0.303 # Table 4 case 2 pTMDD IIV Q = 0.303; S1 Text (h) OMEGA 1 = 0.302

    # Residual error: W = SQRT(THETA(11)^2 * IPRED^2 + THETA(12)^2), SIGMA 1 FIX.
    propSd <- 0.332; label("Proportional residual error (fraction)") # Table 4 case 2 pTMDD proportional = 0.332; S1 Text (h) THETA(11)
    addSd <- 0.00486; label("Additive residual error (pmol/L)") # Table 4 case 2 pTMDD additive = 0.00486; S1 Text (h) THETA(12)
  })

  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    ka_im <- exp(lka_im)
    ka_sc <- exp(lka_sc)
    q <- exp(lq + etalq)
    v_endosome <- exp(lv_endosome)
    krec <- exp(lkrec)
    kdeg <- exp(lkdeg)
    kint <- exp(lkint)
    c_il7r_t <- exp(lc_il7r_t)
    c_fcrn_t <- exp(lc_fcrn_t)
    kon1 <- exp(lkon1)
    kd1 <- exp(lkd1)
    kon2 <- exp(lkon2)
    koff2 <- exp(lkoff2)
    # S1 Text (h) also sets 'IF(TIME.EQ.0) Kup=0', which acts only on the empty
    # interval up to time 0 under NONMEM's end-of-interval convention.
    kup <- exp(lkup)

    # S1 Text (h) $PK: KSS1 = KD + Krec/KON1; KSS2 = (Kint + KOFF2)/KON2.
    kss1 <- kd1 + krec / kon1
    kss2 <- (kint + koff2) / kon2

    v_isf <- ROUTE_SC * exp(lv_isf_sc) + (1 - ROUTE_SC) * exp(lv_isf_im)
    kel <- cl / vc
    k43 <- q / v_isf
    k34 <- q / vc

    # S1 Text (h) $DES, first-order approximation (Byun 2024 Table 1, pTMDD).
    # Endosomal space: Atot = A(5)/VD; Acpx = FcRn*Atot/(KSS1+Atot+FcRn); Afree = Atot-Acpx.
    c_endo_tot <- endosome / v_endosome
    c_fcrn_d <- c_fcrn_t * c_endo_tot / (kss1 + c_endo_tot + c_fcrn_t)
    c_endo_free <- c_endo_tot - c_fcrn_d
    # Central: Ctot = A(3)/VC; Ccpx = Rtot*Ctot/(KSS2+Ctot+Rtot); Cfree = Ctot-Ccpx.
    c_tot <- central / vc
    c_cpx <- c_il7r_t * c_tot / (kss2 + c_tot + c_il7r_t)
    c_free <- c_tot - c_cpx

    # States in the S1 Text (h) $MODEL order: IMDEPOT, SCDEPOT, CENTRAL,
    # DISTRIBUTION1, DISTRIBUTION2.
    d/dt(depot) <- -ka_im * depot
    d/dt(depot2) <- -ka_sc * depot2
    d/dt(central) <- k43 * isf - k34 * c_free * vc + krec * c_fcrn_d * fr * v_endosome -
      kel * c_free * vc - kint * c_cpx * vc
    d/dt(isf) <- ka_im * depot + ka_sc * depot2 - kup * isf + krec * c_fcrn_d * (1 - fr) * v_endosome -
      k43 * isf + k34 * c_free * vc
    d/dt(endosome) <- kup * isf - krec * c_fcrn_d * fr * v_endosome - krec * c_fcrn_d * (1 - fr) * v_endosome -
      kdeg * c_endo_free * v_endosome

    f(depot) <- 2
    f(depot2) <- 2

    # S1 Text (h) $ERROR: Cf = Ct - Ccp, the free concentration; IPRED = BASE + Cf.
    Cc <- BL_IL7 + c_free
    Cc ~ add(addSd) + prop(propSd)
  })
}
