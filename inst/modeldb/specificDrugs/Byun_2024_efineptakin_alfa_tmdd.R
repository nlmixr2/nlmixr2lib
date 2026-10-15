Byun_2024_efineptakin_alfa_tmdd <- function() {
  description <- paste(
    "Full target-mediated drug disposition (TMDD) population PK model for efineptakin alfa",
    "(rhIL-7-hyFc, GX-I7) given by intramuscular or subcutaneous injection to healthy adults,",
    "with FcRn-mediated recycling; case 2 'TMDD' model of Byun 2024. Each route has its own",
    "depot (IM: depot, SC: depot2) feeding an interstitial/vascular absorption space (isf) whose",
    "volume differs by route and which exchanges with central (Q). Drug is taken up from isf",
    "into a lumped endosomal space (Kup), where it is degraded (Kdeg) or binds FcRn by explicit",
    "association / dissociation (kon1, koff1); the drug-FcRn complex recycles (Krec) to central",
    "(fraction FR) and to isf (1 - FR). In central the drug binds IL-7R by explicit kinetics",
    "(kon2, koff2) into a complex internalised at Kint, and is eliminated linearly (CL). Doses",
    "are pmol of GX-I7 and enter with bioavailability 2 (two IL-7 moieties per molecule); the",
    "model works in IL-7 pmol and pmol/L. The observation is the free drug concentration plus",
    "the subject's observed baseline endogenous IL-7 (BL_IL7)."
  )
  reference <- paste(
    "Byun JH, Jeon HS, Yun H-y, Kim JK. Validity conditions of approximations for a",
    "target-mediated drug disposition model: A novel first-order approximation and its",
    "comparison to other approximations. PLoS Comput Biol. 2024;20(4):e1012066.",
    "doi:10.1371/journal.pcbi.1012066. Case 2 TMDD model: Table 4 (estimates) and",
    "S1 Text (e) (NONMEM control stream). The clinical study is Lee SW et al. hIL-7-hyFc, A",
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
      analyte = "efineptakin alfa (free, IL-7 basis)",
      units = "pmol",
      specimen = "serum",
      verified = TRUE
    ),
    isf = list(analyte = "efineptakin alfa (IL-7 basis)", units = "pmol", specimen = "tissue", verified = TRUE),
    endosome = list(
      analyte = "efineptakin alfa (free, IL-7 basis)",
      units = "pmol",
      specimen = "endosome",
      verified = TRUE
    ),
    complex = list(analyte = "efineptakin alfa-IL-7R complex", units = "pmol", specimen = "serum", verified = TRUE),
    complex_fcrn = list(
      analyte = "efineptakin alfa-FcRn complex",
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
      notes = "Per-subject cohort indicator. Selects the volume of the interstitial/vascular absorption space: S1 Text (e) $PK 'VI = VISC; IF(COHORT.EQ.3) VI = VIIM'. SC doses go to depot2 and IM doses to depot.",
      source_name = "COHORT (1 = SC20, 2 = SC60, 3 = IM60); ROUTE_SC = 1 when COHORT is 1 or 2"
    ),
    BL_IL7 = list(
      description = "Observed baseline endogenous serum IL-7 concentration of the subject",
      units = "pmol/L",
      type = "continuous",
      reference_category = "n/a",
      notes = "Added to the predicted drug concentration: S1 Text (e) $ERROR 'IPRED = BASE + Cf ; add baseline IL7 conc'. Range 1.51-4.98 pmol/L in the S2 Table dataset.",
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
    # Final estimates are Byun 2024 Table 4 (case 2, TMDD column). VIIM and VISC
    # are not in Table 4 and are taken from the S1 Text (e) $THETA, which carries
    # the final estimates re-used as starting values.
    lcl <- log(5.26); label("Clearance of free drug from central (L/h)") # Table 4 case 2 TMDD CL = 5.26; S1 Text (e) THETA(1)
    lvc <- log(1.96); label("Central volume of distribution (L)") # Table 4 case 2 TMDD Vc = 1.96; S1 Text (e) THETA(2)
    lka_im <- log(1.22); label("Absorption rate from the intramuscular depot (1/h)") # Table 4 case 2 TMDD Ka = 1.22 (im); S1 Text (e) THETA(3) KAIM
    lka_sc <- log(0.665); label("Absorption rate from the subcutaneous depot (1/h)") # Table 4 case 2 TMDD Ka = 0.665 (sc); S1 Text (e) THETA(4) KASC 0.667
    lq <- log(1.21); label("Inter-compartmental clearance between central and the absorption space (L/h)") # Table 4 case 2 TMDD Q = 1.21; S1 Text (e) THETA(5)
    lv_endosome <- log(31.4); label("Apparent volume of the endosomal FcRn space (L)") # Table 4 case 2 TMDD Vd = 31.4; S1 Text (e) THETA(7) VD
    lkrec <- log(0.000607); label("Recycling rate of the drug-FcRn complex (1/h)") # Table 4 case 2 TMDD Krecycle = 0.000607; S1 Text (e) THETA(8) 0.000608
    lv_isf_im <- log(338); label("Volume of the absorption space, intramuscular cohort (L)") # S1 Text (e) THETA(13) VIIM = 338 (not in Table 4)
    lv_isf_sc <- log(956); label("Volume of the absorption space, subcutaneous cohorts (L)") # S1 Text (e) THETA(14) VISC = 956 (not in Table 4)

    # Held constant at literature values (Table 4 asterisk, or hard-coded in the
    # S1 Text (e) $PK block).
    lkdeg <- fixed(log(0.642)); label("Degradation rate of free drug in the endosomal space (1/h)") # Table 4 case 2 Kdeg = 0.642 (literature); S1 Text (e) THETA(6) FIX
    lkint <- fixed(log(0.642)); label("Internalisation rate of the IL-7R-drug complex (1/h)") # Table 4 case 2 Kint = 0.642 (literature); S1 Text (e) THETA(9) FIX
    lc_il7r_t <- fixed(log(1060)); label("Total IL-7R concentration in central (pmol/L)") # Table 4 case 2 Rtot = 1060 (literature); S1 Text (e) THETA(10) FIX
    lkup <- fixed(log(0.00952)); label("Uptake rate from the absorption space into the endosomal space (1/h)") # Table 4 case 2 Kuptake = 0.00952; S1 Text (e) $PK Kup = 0.00952
    lkon1 <- fixed(log(4.068e-4)); label("Drug-FcRn association rate constant (L/pmol/h)") # S1 Text (e) $PK KON1 = 4.068/10000
    lkoff1 <- fixed(log(12.132)); label("Drug-FcRn dissociation rate constant (1/h)") # S1 Text (e) $PK KOFF1 = 12.132
    lc_fcrn_t <- fixed(log(4e7)); label("Total FcRn concentration in the endosomal space (pmol/L)") # S1 Text (e) $PK FcRn = 40000000
    fr <- fixed(0.715); label("Fraction of recycled drug-FcRn complex returned to central (fraction)") # S1 Text (e) $PK FR = 0.715
    lkon2 <- fixed(log(2.0196e-3)); label("Drug-IL-7R association rate constant (L/pmol/h)") # S1 Text (e) $PK KON2 = 20.196/10000
    lkoff2 <- fixed(log(59.76)); label("Drug-IL-7R dissociation rate constant (1/h)") # S1 Text (e) $PK KOFF2 = 59.76

    etalq ~ 0.292 # Table 4 case 2 TMDD IIV Q = 0.292; S1 Text (e) OMEGA 1

    # Residual error: W = SQRT(THETA(11)^2 * IPRED^2 + THETA(12)^2), SIGMA 1 FIX.
    propSd <- 0.332; label("Proportional residual error (fraction)") # Table 4 case 2 TMDD proportional = 0.332; S1 Text (e) THETA(11)
    addSd <- 0.00303; label("Additive residual error (pmol/L)") # Table 4 case 2 TMDD additive = 0.00303; S1 Text (e) THETA(12)
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
    koff1 <- exp(lkoff1)
    kon2 <- exp(lkon2)
    koff2 <- exp(lkoff2)
    # S1 Text (e) also sets 'IF(TIME.EQ.0) Kup=0'. NONMEM applies the parameters
    # of a record over the interval ENDING at that record, so this touches only
    # the empty interval up to time 0 and is not reproduced.
    kup <- exp(lkup)

    # S1 Text (e) $PK: VI = VISC; IF(COHORT.EQ.3) VI = VIIM.
    v_isf <- ROUTE_SC * exp(lv_isf_sc) + (1 - ROUTE_SC) * exp(lv_isf_im)
    kel <- cl / vc
    k43 <- q / v_isf
    k34 <- q / vc

    # States in the S1 Text (e) $MODEL order: IMDEPOT, SCDEPOT, CENTRAL,
    # DISTRIBUTION1, DISTRIBUTION2, CPX1 (IL-7R complex), CPX2 (FcRn complex).
    # central and endosome hold FREE drug; binding fluxes are written as in the
    # control stream, kon * A * R_total - kon * C * complex.
    d/dt(depot) <- -ka_im * depot
    d/dt(depot2) <- -ka_sc * depot2
    d/dt(central) <- k43 * isf - k34 * central + krec * complex_fcrn * fr - kel * central -
      kon2 * central * c_il7r_t + (kon2 * (central / vc) + koff2) * complex
    d/dt(isf) <- ka_im * depot + ka_sc * depot2 - kup * isf + krec * complex_fcrn * (1 - fr) -
      k43 * isf + k34 * central
    d/dt(endosome) <- kup * isf - kdeg * endosome - kon1 * endosome * c_fcrn_t +
      (kon1 * (endosome / v_endosome) + koff1) * complex_fcrn
    d/dt(complex) <- kon2 * central * c_il7r_t - (kon2 * (central / vc) + koff2 + kint) * complex
    d/dt(complex_fcrn) <- kon1 * endosome * c_fcrn_t - (kon1 * (endosome / v_endosome) + koff1 + krec) * complex_fcrn

    # S1 Text (e) $PK: F1 = 2, F2 = 2 -- dataset AMT is GX-I7 pmol, the model
    # works in IL-7 pmol (two IL-7 moieties per GX-I7 molecule).
    f(depot) <- 2
    f(depot2) <- 2

    # S1 Text (e) $ERROR: Cf = A(3)/VC (free drug); IPRED = BASE + Cf.
    c_free <- central / vc
    Cc <- BL_IL7 + c_free
    Cc ~ add(addSd) + prop(propSd)
  })
}
