Apgar_2018_ALXN1540_human_qsp <- function() {
  description <- "QSP. Humanized (Crigler-Najjar syndrome type 1) translation of the ALXN1540 (hUGT1A1-modRNA lipid nanoparticle) mechanistic model used for first-in-human dose projection. Twelve states chain lipid-nanoparticle plasma disposition and hepatocyte attachment, endocytosis and endosomal escape, to cytoplasmic mRNA, UGT1A1 protein translation and turnover, and two-step enzymatic glucuronidation of bilirubin via explicit enzyme-substrate complexes. Biophysical and cellular rate constants are carried over from the Gunn rat fit; bilirubin production and clearance are human CN1 values. Amounts in nmol; second-order terms use amount/plasma-volume concentrations in nM."
  reference <- "Apgar JF, Tang JP, Singh P, Balasubramanian N, Burke JM, Hodges MR, Lasaro MA, Lin L, Millard BL, Moore K, Jun LS, Sobolov S, Wilkins AK, Gao X. Quantitative Systems Pharmacology Model of hUGT1A1-modRNA Encoding for the UGT1A1 Enzyme to Treat Crigler-Najjar Syndrome Type 1. CPT Pharmacometrics Syst Pharmacol. 2018;7(6):404-412. doi:10.1002/psp4.12301. Author names as corrected by the Corrigendum, CPT Pharmacometrics Syst Pharmacol. 2020;9:185, doi:10.1002/psp4.12484 (author-name corrections only; no parameter, equation or unit is revised). Rate constants from Table 1 (Human column); reaction network and plasma volume from Supplemental Model 2 (CN1 Human KroneckerBio Model File); bilirubin half-life and production from Supplementary Table S1."
  vignette <- "Apgar_2018_ALXN1540"
  # Mechanistic states of the published reaction network (Supplemental Model
  # 2). These are genuinely paper-specific: an LNP delivery chain feeding
  # hepatic mRNA translation and a bilirubin glucuronidation cascade does not
  # map onto the canonical PK compartment vocabulary. ugt_mgt is the extra
  # enzyme-substrate complex the human file carries and the rat file does not.
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
    "ugt_mgt",
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
    ugt_mgt = list(
      analyte = "monoglucuronide-UGT1A1 enzyme-substrate complex",
      units = "nmol",
      specimen = "tissue",
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

  covariatesDataExcluded <- list()

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 0L,
    disease_state = "Crigler-Najjar syndrome type 1 (autosomal recessive absence of UGT1A1 activity, chronic unconjugated hyperbilirubinemia)",
    dose_range = "Simulated 1.0E-05 to 1.0 mg/kg i.v., as a single dose or Q1W, Q2W, Q3W or Q4W (Figure 5); 0.5 mg/kg Q4W identified as the likely efficacious regimen",
    regions = "Not applicable (model-based first-in-human projection; no clinical data)",
    notes = "This is a forward projection, not a fit: no CN1 patients were dosed. The humanization is a hybrid described in the paper's 'Translation of the QSP model to predict response in the CN1 subjects' section. Biophysical and cellular parameters (LNP disposition, endocytosis, endosomal escape, mRNA and UGT1A1 turnover, binding and catalysis) are held at the Gunn rat estimates because they represent processes expected to be preserved between species; bilirubin clearance is set to the 156 h half-life measured by Schmid and Hammaker (1963) in a 4.5-year-old subject with congenital nonhemolytic jaundice; plasma volume is the physiological 3 L. Reported bilirubin production rates are similar across species (3.0 mg/kg/day in the CN1 case study, 3.25 in the Gunn rat, 3.8 in healthy adults; Supplementary Table S1)."
  )

  ini({
    # =====================================================================
    # All rate constants are published in 1/s (Table 1, Human column) and are
    # multiplied by 3600 here so the model runs on a time base of hours. The
    # published 1/s value is kept visible inside the log() so the
    # source-trace is a direct read of Table 1.
    #
    # Every parameter except the three bilirubin/conjugate elimination
    # constants and the plasma volume is identical to the Gunn rat column:
    # the paper holds the biophysical and cellular processes fixed across
    # species and changes only bilirubin production and clearance.
    # =====================================================================

    # ---- Lipid-nanoparticle plasma disposition -------------------------
    lkw <- log(2.41e-05 * 3600); label("First-order elimination of LNP from plasma, kw (1/h)") # Table 1, kw = 2.41E-05 1/s (Estimated in the rat, assumed identical in human)
    lk12 <- log(4.79e-05 * 3600); label("LNP distribution central to peripheral, k12 (1/h)") # Table 1, k12 = 4.79E-05 1/s
    lk21 <- log(2.65e-07 * 3600); label("LNP distribution peripheral to central, k21 (1/h)") # Table 1, k21 = 2.65E-07 1/s

    # ---- Hepatocyte uptake and endosomal escape ------------------------
    lka <- log(1.17e-05 * 3600); label("LNP attachment to hepatocyte, ka (1/h)") # Table 1, ka = 1.17E-05 1/s
    lke <- log(7.70e-05 * 3600); label("LNP endocytosis, ke (1/h)") # Table 1, ke = 7.70E-05 1/s
    lde <- log(9.32e-05 * 3600); label("Endosomal degradation, de (1/h)") # Table 1, de = 9.32E-05 1/s
    lkl <- log(1.93e-05 * 3600); label("Escape from endosome to cytoplasm, kl (1/h)") # Table 1, kl = 1.93E-05 1/s

    # ---- mRNA and UGT1A1 protein turnover ------------------------------
    ldmrna <- log(1.07e-05 * 3600); label("Cytoplasmic mRNA degradation, dmRNA (1/h)") # Table 1, dmRNA = 1.07E-05 1/s
    lkt <- log(17.73 * 3600); label("Translation rate of UGT1A1 from mRNA, kt (1/h)") # Table 1, kt = 17.73 1/s
    ldugt <- log(6.76e-06 * 3600); label("Cytoplasmic UGT1A1 protein degradation, dUGTc (1/h)") # Table 1, dUGTc = 6.76E-06 1/s

    # Endogenous (drug-independent) UGT1A1 translation. Reaction (8) of
    # Supplemental Model 2 carries a zero-order background translation
    # constant ktbg, which is NOT reported in Table 1 or anywhere else in the
    # paper. It is fixed to zero here because CN1 is defined by the absence
    # of UGT1A1 activity, and a non-zero ktbg would conjugate bilirubin
    # before dosing and so contradict Table 1's own rule that kprod is set to
    # reproduce a steady-state total bilirubin governed by kclearBil alone.
    ktbg <- fixed(0); label("Endogenous background UGT1A1 translation (nmol/L/h)") # Supplemental Model 2 reaction (8); value not reported - see vignette Errata

    # ---- Bilirubin conjugation kinetics --------------------------------
    lkon <- log(0.001 * 3600); label("Enzyme-substrate association rate, kon (1/nM/h)") # Table 1, kon = 0.001 1/nM/s (Fixed); typical protein-protein on-rate
    lkoff <- log(0.2589 * 3600); label("Enzyme-substrate dissociation rate, koff (1/h)") # Table 1, koff = 0.2589 1/s (Fixed)
    lkcat <- log(0.0011 * 3600); label("UGT1A1 glucuronidation catalytic rate, kcat (1/h)") # Table 1, kcat = 0.0011 1/s (Fixed); measured in liver microsomes

    # ---- Bilirubin and conjugate elimination ---------------------------
    lkclearbil <- fixed(log(1.23e-06 * 3600)); label("Elimination of unconjugated bilirubin, kclearBil (1/h)") # Table 1, kclearBil = 1.23E-6 1/s (Fixed); CN1 half-life 156 h (Schmid 1963)
    lkclearmgt <- fixed(log(1.23e-05 * 3600)); label("Elimination of bilirubin monoglucuronide, kclearMGT (1/h)") # Table 1, kclearMGT = 1.23E-5 1/s (Fixed); assumed equal to DGT
    lkcleardgt <- fixed(log(1.23e-05 * 3600)); label("Elimination of bilirubin diglucuronide, kclearDGT (1/h)") # Table 1, kclearDGT = 1.23E-5 1/s (Fixed); assumed 10x the unmodified rate

    # ---- Bilirubin production ------------------------------------------
    # Table 1 gives kprod as 'Variable ... Set to match the steady-state
    # bilirubin levels in the modeled individual or group'. It is therefore
    # parameterised here by the pre-dose steady-state total bilirubin
    # (bl_bilirubin) and derived in model() as kprod = kclearBil * bl_bilirubin.
    # Every panel of Figure 5 starts from 20 mg/dL, which at vc = 3.0 L and a
    # bilirubin molar mass of 584.66 g/mol is 20 * 3.0 / 5.8466E-05 = 1.026E+06 nmol.
    bl_bilirubin <- 1.0262e6; label("Pre-dose steady-state total bilirubin (nmol)") # Figure 5a-e, all curves start at 20 mg/dL; converted at vc = 3.0 L

    # Latent surge in bilirubin production. Reaction (15) of Supplemental
    # Model 2 adds a decaying extra production term through the latent state
    # sBil. Neither its production input (ksynhigh) nor its transfer constant
    # (kelSbil) is reported anywhere in the paper, and the deposited human
    # file sets the ksynhigh input to a placeholder with sBil_0 = 0, so both
    # are fixed to zero. This leaves the constant bilirubin production rate
    # the paper used for all first-in-human projections.
    ksynhigh <- fixed(0); label("Surge bilirubin production input, ksynhigh (nmol/h)") # Supplemental Model 2 reaction (15); value not reported - see vignette Errata
    kelsurge <- fixed(0); label("Transfer of the production surge into bilirubin, kelSbil (1/h)") # Supplemental Model 2 reaction (15); value not reported - see vignette Errata

    # ---- Volumes and unit conversions ----------------------------------
    lvc <- fixed(log(3.0)); label("Plasma volume of the central compartment, Vc (L)") # Supplemental Model 2, % Compartments: Plasma 3.0 L; 'Physiological value of 3L was used for the plasma volume in the human model'
    mw_bil <- fixed(584.66); label("Bilirubin molar mass (g/mol)") # Physical constant for bilirubin (C33H36N4O6); used only to report total bilirubin in mg/dL

    # Amount of LNP-encapsulated hUGT1A1-modRNA delivered per mg/kg of dose.
    # The paper states neither the modRNA molar mass nor the body weight used
    # in the human simulations. The deposited KroneckerBio files initialise
    # every state as seed x compartment volume, so a given mg/kg dose enters
    # as the same seed CONCENTRATION in both species and the delivered amount
    # scales with plasma volume. The rat anchor digitized from Figure 2A
    # (0.31 nmol per mg/kg at vc = 0.0078 L) therefore becomes
    # 0.31 * 3.0 / 0.0078 = 119.2 nmol per mg/kg in the human model.
    dose_scale <- fixed(119.2); label("LNP delivered per unit dose (nmol per mg/kg)") # Derived from the Figure 2A rat anchor by the deposited seed x volume convention - see vignette Errata
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
    #    Numbering follows Supplemental Model 2, % Reactions Plasma.
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
    v11_bind <- kon * ugt_cyto * bil / vc # (11) glucuronidation, bilirubin binding
    v11_unbind <- koff * bil_ugt # (11) glucuronidation, bilirubin unbinding
    v12_cat <- kcat * bil_ugt # (12) glucuronidation, bilirubin to monoglucuronide
    # (13) The human file resolves the second glucuronidation step into an
    # explicit UGTc:MGT complex with the same kon/koff/kcat constants. The
    # Gunn rat file instead collapses it into a single second-order step; see
    # the vignette Errata.
    v13_bind <- kon * ugt_cyto * mgt / vc # (13) glucuronidation, monoglucuronide binding
    v13_unbind <- koff * ugt_mgt # (13) glucuronidation, monoglucuronide unbinding
    v13_cat <- kcat * ugt_mgt # (13) glucuronidation, monoglucuronide to diglucuronide
    v15_surge <- kelsurge * bil_surge # (15) transfer of the production surge into bilirubin
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
    d/dt(ugt_cyto) <-
      v08_bgtrans + v09_trans - v10_udeg -
      v11_bind + v11_unbind + v12_cat -
      v13_bind + v13_unbind + v13_cat
    d/dt(bil_ugt) <- v11_bind - v11_unbind - v12_cat - v10_cdeg
    d/dt(bil) <- kprod + v15_surge - v17_clearbil - v11_bind + v11_unbind + v10_cdeg
    d/dt(mgt) <- v12_cat - v13_bind + v13_unbind - v18_clearmgt
    d/dt(ugt_mgt) <- v13_bind - v13_unbind - v13_cat
    d/dt(dgt) <- v13_cat - v19_cleardgt
    d/dt(bil_surge) <- ksynhigh - v15_surge

    # =====================================================================
    # 5. Initial conditions
    # =====================================================================
    # All bilirubin is unconjugated before dosing, and sits at the
    # steady state implied by kprod above. Supplemental Model 2, % Seeds:
    # every seed including sBil_0 is 0.
    bil(0) <- bl_bilirubin

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
    # Supplemental Model 2, % Outputs.
    plasma_mrna <- lnp_central # PlasmaDrug = LNP (nmol)
    liver_mrna <- mrna_cyto + lnp_attached + lnp_endosome # mRNA = mRNAc LNPa LNPe (nmol)
    cyto_mrna <- mrna_cyto # cytomRNA = mRNAc (nmol)
    liver_ugt <- ugt_cyto + bil_ugt + ugt_mgt # Enzyme = UGTc Bil:UGTc, extended with the UGTc:MGT complex (nmol)
    tbili_nmol <- bil + mgt + dgt # TotalBilirubin = Bil MGT DGT (nmol)
    # Total bilirubin expressed in the clinical unit used throughout Figure 5.
    # 1 nmol in vc litres is mw_bil * 1e-6 mg, and dividing by 10 converts
    # mg/L to mg/dL.
    tbili_mgdl <- tbili_nmol * mw_bil * 1e-7 / vc

    # The paper reports no between-subject variability and no residual-error
    # estimates for the human projection - it is a deterministic forward
    # simulation - so no error model is attached.
  })
}
