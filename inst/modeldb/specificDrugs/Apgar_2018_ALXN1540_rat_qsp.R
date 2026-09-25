Apgar_2018_ALXN1540_rat_qsp <- function() {
  description <- "QSP (preclinical, Gunn rat). Mechanistic model of ALXN1540 (hUGT1A1-modRNA in a lipid nanoparticle) for Crigler-Najjar syndrome type 1. Eleven states chain lipid-nanoparticle plasma disposition and hepatocyte attachment, endocytosis and endosomal escape, to cytoplasmic mRNA, UGT1A1 protein translation and turnover, and enzymatic glucuronidation of bilirubin to mono- and diglucuronide with their faster elimination. Amounts in nmol; second-order terms use amount/plasma-volume concentrations in nM."
  reference <- "Apgar JF, Tang JP, Singh P, Balasubramanian N, Burke JM, Hodges MR, Lasaro MA, Lin L, Millard BL, Moore K, Jun LS, Sobolov S, Wilkins AK, Gao X. Quantitative Systems Pharmacology Model of hUGT1A1-modRNA Encoding for the UGT1A1 Enzyme to Treat Crigler-Najjar Syndrome Type 1. CPT Pharmacometrics Syst Pharmacol. 2018;7(6):404-412. doi:10.1002/psp4.12301. Author names as corrected by the Corrigendum, CPT Pharmacometrics Syst Pharmacol. 2020;9:185, doi:10.1002/psp4.12484 (author-name corrections only; no parameter, equation or unit is revised). Rate constants and plasma volume from Table 1 (Gunn rat column); reaction network from Supplementary Model S1 (Gunn Rat KroneckerBio Model File)."
  vignette <- "Apgar_2018_ALXN1540"
  # Mechanistic states of the published reaction network (Supplementary Model
  # S1). These are genuinely paper-specific: an LNP delivery chain feeding
  # hepatic mRNA translation and a bilirubin glucuronidation cascade does not
  # map onto the canonical PK compartment vocabulary.
  paper_specific_compartments <- c(
    "lnp_central",
    "lnp_peripheral",
    "lnp_attached",
    "lnp_endosome",
    "mrna_cyto",
    "ugt_cyto",
    "bil",
    "bil_ugt",
    "mgt",
    "dgt",
    "bil_surge"
  )
  units <- list(time = "h", dosing = "mg/kg", concentration = "nM")

  compartmentData <- list(
    lnp_central = list(
      analyte = "intact ALXN1540 lipid nanoparticle (hUGT1A1-modRNA)",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    lnp_peripheral = list(
      analyte = "intact ALXN1540 lipid nanoparticle (hUGT1A1-modRNA)",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    lnp_attached = list(
      analyte = "ALXN1540 lipid nanoparticle bound to the hepatocyte surface",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    lnp_endosome = list(
      analyte = "ALXN1540 lipid nanoparticle inside the hepatocyte endosome",
      units = "nmol",
      specimen = "endosome",
      verified = TRUE
    ),
    mrna_cyto = list(
      analyte = "hUGT1A1-modRNA in the hepatocyte cytoplasm",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    ugt_cyto = list(
      analyte = "free UGT1A1 enzyme in the hepatocyte cytoplasm",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    bil = list(
      analyte = "unconjugated bilirubin",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    bil_ugt = list(
      analyte = "bilirubin-UGT1A1 enzyme-substrate complex",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    mgt = list(
      analyte = "bilirubin monoglucuronide",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    dgt = list(
      analyte = "bilirubin diglucuronide",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    bil_surge = list(
      analyte = "not applicable",
      units = "nmol",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units = "unitless",
      type = "binary",
      notes = "The repeat-dose PD study enrolled male and female animals, but the paper reports 'no clear gender or age differences were noted in the PK data' and pooled all four groups for calibration; sex is therefore not a model covariate."
    ),
    AGE = list(
      description = "Age group (juvenile vs adult Gunn rat)",
      units = "unitless",
      type = "categorical",
      notes = "Adolescent/juvenile and adult animals were pooled for calibration. Per the Figure 2 caption the ONLY quantity that differed between the adult and juvenile fits was the initial total bilirubin level, which is carried by the bl_bilirubin parameter rather than by an age covariate."
    )
  )

  population <- list(
    species = "rat (Gunn, Gunn-UGT1a1j/BluHsdRrrc)",
    n_subjects = NA_integer_,
    n_studies = 3L,
    disease_state = "UGT1A1-deficient Gunn rat, the animal model of Crigler-Najjar syndrome type 1 (unconjugated hyperbilirubinemia)",
    dose_range = "0.3 mg/kg single i.v. bolus (single-dose PK/PD study); 0.1, 0.2 and 0.5 mg/kg i.v. Q2W and 0.5 mg/kg i.v. Q4W for five doses (repeat-dose PD study); 0.1, 0.2 and 0.5 mg/kg single i.v. bolus (escalating single-dose prediction study)",
    regions = "Preclinical (Alexion Pharmaceuticals / Applied BioMath)",
    notes = "Three studies: a single-dose 0.3 mg/kg PK/PD study in adult (N = 4-5) and juvenile (N = 5) rats (Figure 2); a repeat-dose PD study in adult rats with a 0.5 mg/kg Q2W luciferase-mRNA vehicle control (N = 4-5 per time point, Figure 3); and an escalating single-dose PD study used as a naive prediction check (N = 4-6 per time point, Figure 4). Body weights recorded in the deposited single-dose dataset were 155-262 g (adults) and 70-116 g (juveniles). All doses were given i.v. via the tail vein. The model is deterministic: the paper reports no between-animal random effects and no residual-error estimates, and states only that a 20% variation around the parameter values captured the observed data variability."
  )

  ini({
    # =====================================================================
    # All rate constants are published in 1/s (Table 1, Gunn rat column) and
    # are multiplied by 3600 here so the model runs on a time base of hours.
    # The published 1/s value is kept visible inside the log() so the
    # source-trace is a direct read of Table 1.
    # =====================================================================

    # ---- Lipid-nanoparticle plasma disposition -------------------------
    lkw <- log(2.41e-05 * 3600); label("First-order elimination of LNP from plasma, kw (1/h)") # Table 1, kw = 2.41E-05 1/s (Estimated); plasma half-life ~8 h
    lk12 <- log(4.79e-05 * 3600); label("LNP distribution central to peripheral, k12 (1/h)") # Table 1, k12 = 4.79E-05 1/s (Estimated)
    lk21 <- log(2.65e-07 * 3600); label("LNP distribution peripheral to central, k21 (1/h)") # Table 1, k21 = 2.65E-07 1/s (Estimated)

    # ---- Hepatocyte uptake and endosomal escape ------------------------
    lka <- log(1.17e-05 * 3600); label("LNP attachment to hepatocyte, ka (1/h)") # Table 1, ka = 1.17E-05 1/s (Estimated)
    lke <- log(7.70e-05 * 3600); label("LNP endocytosis, ke (1/h)") # Table 1, ke = 7.70E-05 1/s (Estimated)
    lde <- log(9.32e-05 * 3600); label("Endosomal degradation, de (1/h)") # Table 1, de = 9.32E-05 1/s (Estimated)
    lkl <- log(1.93e-05 * 3600); label("Escape from endosome to cytoplasm, kl (1/h)") # Table 1, kl = 1.93E-05 1/s (Estimated)

    # ---- mRNA and UGT1A1 protein turnover ------------------------------
    ldmrna <- log(1.07e-05 * 3600); label("Cytoplasmic mRNA degradation, dmRNA (1/h)") # Table 1, dmRNA = 1.07E-05 1/s (Estimated); ~18 h half-life
    lkt <- log(17.73 * 3600); label("Translation rate of UGT1A1 from mRNA, kt (1/h)") # Table 1, kt = 17.73 1/s (Estimated)
    ldugt <- log(6.76e-06 * 3600); label("Cytoplasmic UGT1A1 protein degradation, dUGTc (1/h)") # Table 1, dUGTc = 6.76E-06 1/s (Estimated)

    # Endogenous (drug-independent) UGT1A1 translation. Reaction (8) of
    # Supplementary Model S1 carries a zero-order background translation
    # constant ktbg, which is NOT reported in Table 1 or anywhere else in the
    # paper. It is fixed to zero here: Gunn rats are the UGT1a1-deficient
    # animal model of CN1, the paper's premise is the absence of UGT1A1
    # activity, and a non-zero ktbg would conjugate bilirubin before dosing
    # and so contradict Table 1's own rule that kprod is set to reproduce a
    # steady-state total bilirubin governed by kclearBil alone.
    ktbg <- fixed(0); label("Endogenous background UGT1A1 translation (nmol/L/h)") # Supplementary Model S1 reaction (8); value not reported - see vignette Errata

    # ---- Bilirubin conjugation kinetics --------------------------------
    lkon <- log(0.001 * 3600); label("Bilirubin-UGT1A1 association rate, kon (1/nM/h)") # Table 1, kon = 0.001 1/nM/s (Fixed); typical protein-protein on-rate
    lkoff <- log(0.2589 * 3600); label("Bilirubin-UGT1A1 dissociation rate, koff (1/h)") # Table 1, koff = 0.2589 1/s (Fixed)
    lkcat <- log(0.0011 * 3600); label("UGT1A1 glucuronidation catalytic rate, kcat (1/h)") # Table 1, kcat = 0.0011 1/s (Fixed); measured in liver microsomes

    # ---- Bilirubin and conjugate elimination ---------------------------
    lkclearbil <- fixed(log(3.5e-06 * 3600)); label("Elimination of unconjugated bilirubin, kclearBil (1/h)") # Table 1, kclearBil = 3.5E-06 1/s (Fixed); Gunn rat half-life 55 h
    lkclearmgt <- fixed(log(3.5e-05 * 3600)); label("Elimination of bilirubin monoglucuronide, kclearMGT (1/h)") # Table 1, kclearMGT = 3.5E-05 1/s (Fixed); assumed equal to DGT
    lkcleardgt <- fixed(log(3.5e-05 * 3600)); label("Elimination of bilirubin diglucuronide, kclearDGT (1/h)") # Table 1, kclearDGT = 3.5E-05 1/s (Fixed); assumed 10x the unmodified rate

    # ---- Bilirubin production ------------------------------------------
    # Table 1 gives kprod as 'Variable ... Set to match the steady-state
    # bilirubin levels in the modeled individual or group'. It is therefore
    # parameterised here by the pre-dose steady-state total bilirubin
    # (bl_bilirubin) and derived in model() as kprod = kclearBil * bl_bilirubin.
    bl_bilirubin <- 455; label("Pre-dose steady-state total bilirubin (nmol)") # Figure 2A fitted curve at t = 0 (adult Gunn rat), digitized; 568 nmol for the juvenile fit in Figure 2B

    # Juvenile natural-history surge in bilirubin production. Reactions (15)
    # and (16) of Supplementary Model S1 add a decaying extra production term
    # driven by the latent state sBil, used to describe the drop in plasma
    # bilirubin from day 0 to day 21 seen in the luciferase control group.
    # Neither its production constant (ksynhigh) nor its decay constant
    # (kelSbil) is reported anywhere in the paper, so both are fixed to zero,
    # which leaves the constant production rate the paper itself assumed for
    # the adult and escalating-dose studies. Override them to explore the
    # juvenile time-varying production.
    ksynhigh <- fixed(0); label("Surge bilirubin production rate, ksynhigh (1/h)") # Supplementary Model S1 reaction (15); value not reported - see vignette Errata
    kelsurge <- fixed(0); label("Decay of the bilirubin production surge, kelSbil (1/h)") # Supplementary Model S1 reaction (16); value not reported - see vignette Errata

    # ---- Volumes and unit conversions ----------------------------------
    lvc <- fixed(log(0.0078)); label("Plasma volume of the central compartment, Vc (L)") # Supplementary Model S1, % Compartments: Plasma 0.0078 L
    mw_bil <- fixed(584.66); label("Bilirubin molar mass (g/mol)") # Physical constant for bilirubin (C33H36N4O6); used only to report total bilirubin in mg/dL

    # Amount of LNP-encapsulated hUGT1A1-modRNA delivered per mg/kg of dose.
    # The paper does not state the modRNA molar mass or the animal body
    # weight used in the simulations, so this conversion is digitized from
    # the fitted plasma-mRNA curve of Figure 2A: a 0.3 mg/kg i.v. bolus
    # starts at 9.3E-02 nmol, i.e. 0.31 nmol per mg/kg. The digitization was
    # cross-checked against the initial log-linear slope of the same curve,
    # which recovers kw + k12 + ka to within 1%.
    dose_scale <- fixed(0.31); label("LNP delivered per unit dose (nmol per mg/kg)") # Digitized from Figure 2A (plasma mRNA at t = 0 for 0.3 mg/kg) - see vignette Errata
  })

  model({
    # =====================================================================
    # 1. Parameters on the natural scale
    # =====================================================================
    kw <- exp(lkw)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    ka <- exp(lka)
    ke <- exp(lke)
    de <- exp(lde)
    kl <- exp(lkl)
    dmrna <- exp(ldmrna)
    kt <- exp(lkt)
    dugt <- exp(ldugt)
    kon <- exp(lkon)
    koff <- exp(lkoff)
    kcat <- exp(lkcat)
    kclearbil <- exp(lkclearbil)
    kclearmgt <- exp(lkclearmgt)
    kcleardgt <- exp(lkcleardgt)
    vc <- exp(lvc)

    # =====================================================================
    # 2. Bilirubin production derived from the pre-dose steady state
    # =====================================================================
    # Before dosing there is no UGT1A1 (ktbg = 0), so the only bilirubin
    # species present is the unconjugated form and d(bil)/dt = 0 gives
    # kprod = kclearBil * bl_bilirubin. This is exactly Table 1's rule for
    # kprod, expressed so the input is the observable baseline.
    kprod <- kclearbil * bl_bilirubin

    # =====================================================================
    # 3. Reaction rates (nmol/h). Second-order reactions use concentrations
    #    (nmol/L = nM), hence the division by vc.
    #    Numbering follows Supplementary Model S1, % Reactions Plasma.
    # =====================================================================
    v01_elim_c <- kw * lnp_central # (1) first-order elimination of plasma LNP
    v01_c2p <- k12 * lnp_central # (1) distribution, central to peripheral
    v01_p2c <- k21 * lnp_peripheral # (1) distribution, peripheral to central
    v01_elim_p <- kw * lnp_peripheral # (1) first-order elimination of peripheral LNP
    v02_attach <- ka * lnp_central # (2) attachment to hepatocyte
    v03_endo <- ke * lnp_attached # (3) endocytosis
    v04_edeg <- de * lnp_endosome # (4) fails to escape from endosome
    v05_escape <- kl * lnp_endosome # (5) escape from endosome
    v07_mdeg <- dmrna * mrna_cyto # (7) mRNA degrades
    v08_bgtrans <- ktbg * vc # (8) endogenous translation (zero order)
    v09_trans <- kt * mrna_cyto # (9) translation
    v10_udeg <- dugt * ugt_cyto # (10) cytoplasmic protein degrades
    v10_cdeg <- dugt * bil_ugt # (10) complexed protein degrades, releasing bilirubin
    v11_bind <- kon * ugt_cyto * bil / vc # (11) glucuronidation, binding
    v11_unbind <- koff * bil_ugt # (11) glucuronidation, unbinding
    v12_cat <- kcat * bil_ugt # (12) glucuronidation, bilirubin to monoglucuronide
    # (13) In the Gunn rat file the monoglucuronide step is a single
    # second-order reaction UGTc + MGT -> UGTc + DGT that reuses kcat as its
    # rate constant. The human file instead resolves this step into an
    # explicit UGTc:MGT complex. Both are reproduced as published; see the
    # vignette Errata for the dimensional inconsistency this introduces here.
    v13_cat2 <- kcat * ugt_cyto * mgt / vc # (13) glucuronidation, monoglucuronide to diglucuronide
    v15_surge <- ksynhigh * bil_surge # (15) high production of bilirubin
    v16_sdecay <- kelsurge * bil_surge # (16) decay of the production surge
    v17_clearbil <- kclearbil * bil # (17) elimination of bilirubin
    v18_clearmgt <- kclearmgt * mgt # (18) elimination of monoglucuronide
    v19_cleardgt <- kcleardgt * dgt # (19) elimination of diglucuronide

    # =====================================================================
    # 4. ODE system
    # =====================================================================
    d/dt(lnp_central) <- -v01_elim_c - v01_c2p + v01_p2c - v02_attach
    d/dt(lnp_peripheral) <- v01_c2p - v01_p2c - v01_elim_p
    d/dt(lnp_attached) <- v02_attach - v03_endo
    d/dt(lnp_endosome) <- v03_endo - v04_edeg - v05_escape
    d/dt(mrna_cyto) <- v05_escape - v07_mdeg
    d/dt(ugt_cyto) <- v08_bgtrans + v09_trans - v10_udeg - v11_bind + v11_unbind + v12_cat
    d/dt(bil_ugt) <- v11_bind - v11_unbind - v12_cat - v10_cdeg
    d/dt(bil) <- kprod + v15_surge - v17_clearbil - v11_bind + v11_unbind + v10_cdeg
    d/dt(mgt) <- v12_cat - v13_cat2 - v18_clearmgt
    d/dt(dgt) <- v13_cat2 - v19_cleardgt
    d/dt(bil_surge) <- -v16_sdecay

    # =====================================================================
    # 5. Initial conditions
    # =====================================================================
    # All bilirubin is unconjugated before dosing, and sits at the
    # steady state implied by kprod above.
    bil(0) <- bl_bilirubin
    # Supplementary Model S1, % Seeds: every seed is 0 except sBil_0 = 1, and
    # each state is initialised as seed x compartment volume.
    bil_surge(0) <- vc

    # =====================================================================
    # 6. Dose conversion
    # =====================================================================
    # Doses are given in mg/kg; dose_scale converts to the nmol of
    # LNP-encapsulated modRNA entering the plasma compartment.
    f(lnp_central) <- dose_scale

    # =====================================================================
    # 7. Observations
    # =====================================================================
    # Cc is the intact plasma LNP concentration (the model's PK observable).
    Cc <- lnp_central / vc # nM
    # Supplementary Model S1, % Outputs.
    plasma_mrna <- lnp_central # PlasmaDrug = LNP (nmol; Figure 2 left panels)
    liver_mrna <- mrna_cyto + lnp_attached + lnp_endosome # mRNA = mRNAc LNPa LNPe (nmol; Figure 2 middle panels)
    cyto_mrna <- mrna_cyto # cytomRNA = mRNAc (nmol)
    liver_ugt <- ugt_cyto + bil_ugt # Enzyme = UGTc Bil:UGTc (nmol; Figure 4 protein panels)
    tbili_nmol <- bil + mgt + dgt # TotalBilirubin = Bil MGT DGT (nmol; Figure 2 and 3 right panels)
    # Total bilirubin expressed in the clinical unit. 1 nmol in vc litres is
    # mw_bil * 1e-6 mg, and dividing by 10 converts mg/L to mg/dL.
    tbili_mgdl <- tbili_nmol * mw_bil * 1e-7 / vc

    # The paper reports no between-animal variability and no residual-error
    # estimates (only that a 20% parameter variation spanned the observed
    # data spread), so no error model is attached. The model is intended for
    # simulation rather than estimation.
  })
}
