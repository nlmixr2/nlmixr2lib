Marcantonio_2022_panitumumab <- function() {
  description <- paste(
    "QSP. Two-compartment monospecific anti-receptor mechanistic PKPD model of",
    "panitumumab-EGFR binding in adults with metastatic colorectal cancer",
    "(Marcantonio 2022 Early Feasibility Assessment, Case Study 2). Drug",
    "administered IV as a bivalent (valency = 2) antibody that binds",
    "membrane EGFR (in central and peripheral compartments). Soluble EGFR",
    "is not modelled (JSON shed_css = 0). All species eliminate first-order;",
    "drug bound to membrane receptor eliminates at the receptor's rate.",
    "Structure identical to the Marcantonio 2022 trastuzumab anti-receptor",
    "model; differs only in the drug-specific parameters. Parameters FIXED",
    "from paper Tables 3 (target) and 4 (drug) and the Assess run file JSON.",
    sep = " "
  )
  reference <- paste(
    "Marcantonio DH et al. (2022). Front Pharmacol 13:864768.",
    "doi:10.3389/fphar.2022.864768. Case Study 2 (panitumumab EGFR,",
    "colorectal cancer). Drug-specific parameters from paper Tables 3 and 4",
    "and the Assess run file JSON (Data Sheet 2).",
    sep = " "
  )
  vignette <- "Marcantonio_2022_efa"
  units <- list(
    time          = "day",
    dosing        = "Panitumumab dose amount into Ab_00_c (IV bolus) in nmol; MW = 150000 Da so 420 mg = 2800 nmol.",
    concentration = "Free panitumumab plasma concentration Cc = Ab_00_c / Vc in nM; central volume Vc = 3 L, peripheral Vp = 13 L. Target engagement = drug-bound R1 / (drug-bound R1 + free R1)."
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = NA_integer_,
    disease_state  = "Adults with metastatic KRAS wild-type colorectal cancer.",
    dose_range     = "Clinically approved dose: 6 mg/kg IV every 2 weeks (equal to 420 mg for a 70 kg patient); Vectibix USPI. Marcantonio 2022 predicts an effective dose (98% TE peripheral) of 162 mg Q2W IV.",
    regions        = NA_character_,
    notes          = "See Marcantonio 2022 trastuzumab sibling for shared 2-cpt anti-receptor structure and Case Study 2 methodology."
  )

  ini({
    lka       <- fixed(log(log(2) / 2.5));         label("First-order SC absorption placeholder rate constant (1/day; unused for IV dosing)")     # Assess default
    lkclearAb <- fixed(log(log(2) / 16));          label("First-order panitumumab elimination rate constant (1/day; from t1/2 = 16 days)")        # Marcantonio 2022 Table 4 Panitumumab Half-Life (Yang 2001, Ma 2009)
    kd_drug   <- fixed(0.05);                      label("Panitumumab-EGFR equilibrium dissociation constant (nM)")                                # Marcantonio 2022 Table 4 Panitumumab KD for EGFR (Yang 2001, Ma 2009)
    valency   <- fixed(2);                         label("Effective drug valency for EGFR binding (dimensionless)")                                 # Marcantonio 2022 Table 4 Panitumumab Valency (Yang 2001, Ma 2009)

    R1_conc_c    <- fixed(0.0152);                 label("EGFR membrane receptor concentration in central compartment (nM)")                        # Marcantonio 2022 Table 3 EGFR expression central = 4.57e-2 nmoles; 0.0457/3 = 0.01523 nM
    R1_conc_p    <- fixed(1.13);                   label("EGFR membrane receptor concentration in peripheral compartment (nM)")                     # Marcantonio 2022 Table 3 EGFR expression peripheral = 1.47e+01 nmoles; 14.7/13 = 1.131 nM
    S1_conc_c    <- fixed(0);                      label("Soluble shed EGFR concentration in central compartment (nM; not modelled)")               # Marcantonio 2022 Assess run file shed_css_1_central = 0 (soluble EGFR not modelled)
    S1_conc_p    <- fixed(0);                      label("Soluble shed EGFR concentration in peripheral compartment (nM; not modelled)")            # Assess run file shed_css_1_peripheral = 0

    lkclearR1    <- fixed(log(log(2) / (5 / 24)));     label("First-order EGFR membrane elimination rate constant (1/day; from t1/2 = 5 hr)")       # Marcantonio 2022 Table 3 EGFR receptor half-life (Sigismund 2008)
    lkclearS1    <- fixed(log(log(2) / (1000 / 24)));  label("First-order soluble EGFR elimination rate constant (placeholder; unused)")             # Marcantonio 2022 Assess run file shed_half_1 = 1000 hr (essentially disabled; unused when S1_conc = 0)

    Vc          <- fixed(3);                       label("Central compartment volume (L)")                                                          # Marcantonio 2022 Table 3 Central compartment volume (plasma volume)
    Vp          <- fixed(13);                      label("Peripheral compartment volume (L)")                                                       # Marcantonio 2022 Table 3 Peripheral compartment volume (interstitial)
    Tdist_Ab_hr <- fixed(35);                      label("Drug distribution half-life to peripheral (hr)")                                          # Marcantonio 2022 Table 4 Tdist12 = 35 hr (Betts 2018)
    Pdist_Ab    <- fixed(0.190625);                label("Drug partition coefficient (peripheral vs central; unitless)")                            # Marcantonio 2022 Table 4 Pdist12 = 0.19 (Betts 2018); Assess uses 0.190625
    Tdist_S1_hr <- fixed(30);                      label("Soluble receptor distribution half-life to peripheral (hr; unused when S1_conc = 0)")     # Assess default

    kon         <- fixed(0.001 * 86400);           label("Bimolecular association rate constant (L/nmol/day)")                                      # Assess default

    propSd    <- fixed(0.01);                      label("Placeholder proportional residual error on Cc (not paper-derived)")                      # Placeholder
  })

  model({
    ka         <- exp(lka)
    kclearAb   <- exp(lkclearAb)
    kclearR1   <- exp(lkclearR1)
    kclearS1   <- exp(lkclearS1)
    kon1Ab     <- kon
    kon2Ab     <- floor(valency / 2) * kon
    koffAb     <- kon * kd_drug

    # Drug distribution rate constants.
    Tdist_Ab_d <- Tdist_Ab_hr / 24
    kout_Ab    <- (log(2) / Tdist_Ab_d) * Pdist_Ab / (Pdist_Ab + Vc / Vp)
    kin_Ab     <- (log(2) / Tdist_Ab_d) / (1 + Pdist_Ab * Vp / Vc)

    # Soluble-receptor distribution (unused for panitumumab since S1_conc = 0
    # but retained to keep structural symmetry with the sibling trastuzumab
    # anti-receptor model).
    Tdist_S1_d <- Tdist_S1_hr / 24
    QS1        <- log(2) / Tdist_S1_d
    Pdist_S1   <- (S1_conc_p + 1e-16) / (S1_conc_c + 1e-16)
    kout_S1    <- QS1 * Pdist_S1 / (Pdist_S1 + Vc / Vp)
    kin_S1     <- QS1 / (1 + Pdist_S1 * Vp / Vc)

    total_R1_c <- R1_conc_c * Vc
    total_R1_p <- R1_conc_p * Vp
    S1_0_c     <- S1_conc_c * Vc
    S1_0_p     <- S1_conc_p * Vp

    # Zero-shed regime (S1_conc = 0) -> avoid divide-by-zero in kshed.
    kshed_R1_c <- kclearS1 * S1_0_c / (total_R1_c + 1e-16)
    kshed_R1_p <- kclearS1 * S1_0_p / (total_R1_p + 1e-16)
    ksynth_R1_c <- (kclearR1 + kshed_R1_c) * total_R1_c
    ksynth_R1_p <- (kclearR1 + kshed_R1_p) * total_R1_p

    d/dt(R1_c) <-  ksynth_R1_c - kclearR1 * R1_c - kshed_R1_c * R1_c -
                    kon1Ab * R1_c * Ab_00_c / Vc + koffAb * Ab_R0_c -
                    kon2Ab * R1_c * Ab_00_c / Vc + koffAb * Ab_0R_c -
                    kon1Ab * R1_c * Ab_0R_c / Vc + koffAb * Ab_RR_c -
                    kon2Ab * R1_c * Ab_R0_c / Vc + koffAb * Ab_RR_c -
                    kon1Ab * R1_c * Ab_0S_c / Vc + koffAb * Ab_RS_c -
                    kon2Ab * R1_c * Ab_S0_c / Vc + koffAb * Ab_SR_c +
                    2 * kclearR1 * Ab_RR_c

    d/dt(S1_c) <-  kshed_R1_c * R1_c - kclearS1 * S1_c -
                    kout_S1 * S1_c + kin_S1 * S1_p -
                    kon1Ab * S1_c * Ab_00_c / Vc + koffAb * Ab_S0_c -
                    kon2Ab * S1_c * Ab_00_c / Vc + koffAb * Ab_0S_c -
                    kon1Ab * S1_c * Ab_0R_c / Vc + koffAb * Ab_SR_c -
                    kon2Ab * S1_c * Ab_R0_c / Vc + koffAb * Ab_RS_c -
                    kon1Ab * S1_c * Ab_0S_c / Vc + koffAb * Ab_SS_c -
                    kon2Ab * S1_c * Ab_S0_c / Vc + koffAb * Ab_SS_c

    d/dt(Ab_00_c) <-  ka * depot -
                       kclearAb * Ab_00_c -
                       kout_Ab * Ab_00_c + kin_Ab * Ab_00_p -
                       kon1Ab * R1_c * Ab_00_c / Vc + koffAb * Ab_R0_c -
                       kon2Ab * R1_c * Ab_00_c / Vc + koffAb * Ab_0R_c -
                       kon1Ab * S1_c * Ab_00_c / Vc + koffAb * Ab_S0_c -
                       kon2Ab * S1_c * Ab_00_c / Vc + koffAb * Ab_0S_c

    d/dt(Ab_0R_c) <-  kon2Ab * R1_c * Ab_00_c / Vc - koffAb * Ab_0R_c -
                       kclearR1 * Ab_0R_c -
                       kon1Ab * R1_c * Ab_0R_c / Vc + koffAb * Ab_RR_c -
                       kon1Ab * S1_c * Ab_0R_c / Vc + koffAb * Ab_SR_c

    d/dt(Ab_R0_c) <-  kon1Ab * R1_c * Ab_00_c / Vc - koffAb * Ab_R0_c -
                       kclearR1 * Ab_R0_c -
                       kon2Ab * R1_c * Ab_R0_c / Vc + koffAb * Ab_RR_c -
                       kon2Ab * S1_c * Ab_R0_c / Vc + koffAb * Ab_RS_c

    d/dt(Ab_RR_c) <-  kon1Ab * R1_c * Ab_0R_c / Vc - koffAb * Ab_RR_c +
                       kon2Ab * R1_c * Ab_R0_c / Vc - koffAb * Ab_RR_c -
                       kclearR1 * Ab_RR_c - kclearR1 * Ab_RR_c

    d/dt(Ab_0S_c) <-  kon2Ab * S1_c * Ab_00_c / Vc - koffAb * Ab_0S_c -
                       kclearAb * Ab_0S_c -
                       kout_Ab * Ab_0S_c + kin_Ab * Ab_0S_p -
                       kon1Ab * R1_c * Ab_0S_c / Vc + koffAb * Ab_RS_c -
                       kon1Ab * S1_c * Ab_0S_c / Vc + koffAb * Ab_SS_c

    d/dt(Ab_S0_c) <-  kon1Ab * S1_c * Ab_00_c / Vc - koffAb * Ab_S0_c -
                       kclearAb * Ab_S0_c -
                       kout_Ab * Ab_S0_c + kin_Ab * Ab_S0_p -
                       kon2Ab * R1_c * Ab_S0_c / Vc + koffAb * Ab_SR_c -
                       kon2Ab * S1_c * Ab_S0_c / Vc + koffAb * Ab_SS_c

    d/dt(Ab_RS_c) <-  kon2Ab * S1_c * Ab_R0_c / Vc - koffAb * Ab_RS_c +
                       kon1Ab * R1_c * Ab_0S_c / Vc - koffAb * Ab_RS_c -
                       kclearR1 * Ab_RS_c

    d/dt(Ab_SR_c) <-  kon2Ab * R1_c * Ab_S0_c / Vc - koffAb * Ab_SR_c +
                       kon1Ab * S1_c * Ab_0R_c / Vc - koffAb * Ab_SR_c -
                       kclearR1 * Ab_SR_c

    d/dt(Ab_SS_c) <-  kon1Ab * S1_c * Ab_0S_c / Vc - koffAb * Ab_SS_c +
                       kon2Ab * S1_c * Ab_S0_c / Vc - koffAb * Ab_SS_c -
                       kclearAb * Ab_SS_c -
                       kout_Ab * Ab_SS_c + kin_Ab * Ab_SS_p

    d/dt(R1_p) <-  ksynth_R1_p - kclearR1 * R1_p - kshed_R1_p * R1_p -
                    kon1Ab * R1_p * Ab_00_p / Vp + koffAb * Ab_R0_p -
                    kon2Ab * R1_p * Ab_00_p / Vp + koffAb * Ab_0R_p -
                    kon1Ab * R1_p * Ab_0R_p / Vp + koffAb * Ab_RR_p -
                    kon2Ab * R1_p * Ab_R0_p / Vp + koffAb * Ab_RR_p -
                    kon1Ab * R1_p * Ab_0S_p / Vp + koffAb * Ab_RS_p -
                    kon2Ab * R1_p * Ab_S0_p / Vp + koffAb * Ab_SR_p +
                    2 * kclearR1 * Ab_RR_p

    d/dt(S1_p) <-  kshed_R1_p * R1_p - kclearS1 * S1_p +
                    kout_S1 * S1_c - kin_S1 * S1_p -
                    kon1Ab * S1_p * Ab_00_p / Vp + koffAb * Ab_S0_p -
                    kon2Ab * S1_p * Ab_00_p / Vp + koffAb * Ab_0S_p -
                    kon1Ab * S1_p * Ab_0R_p / Vp + koffAb * Ab_SR_p -
                    kon2Ab * S1_p * Ab_R0_p / Vp + koffAb * Ab_RS_p -
                    kon1Ab * S1_p * Ab_0S_p / Vp + koffAb * Ab_SS_p -
                    kon2Ab * S1_p * Ab_S0_p / Vp + koffAb * Ab_SS_p

    d/dt(Ab_00_p) <-  -kclearAb * Ab_00_p +
                       kout_Ab * Ab_00_c - kin_Ab * Ab_00_p -
                       kon1Ab * R1_p * Ab_00_p / Vp + koffAb * Ab_R0_p -
                       kon2Ab * R1_p * Ab_00_p / Vp + koffAb * Ab_0R_p -
                       kon1Ab * S1_p * Ab_00_p / Vp + koffAb * Ab_S0_p -
                       kon2Ab * S1_p * Ab_00_p / Vp + koffAb * Ab_0S_p

    d/dt(Ab_0R_p) <-  kon2Ab * R1_p * Ab_00_p / Vp - koffAb * Ab_0R_p -
                       kclearR1 * Ab_0R_p -
                       kon1Ab * R1_p * Ab_0R_p / Vp + koffAb * Ab_RR_p -
                       kon1Ab * S1_p * Ab_0R_p / Vp + koffAb * Ab_SR_p

    d/dt(Ab_R0_p) <-  kon1Ab * R1_p * Ab_00_p / Vp - koffAb * Ab_R0_p -
                       kclearR1 * Ab_R0_p -
                       kon2Ab * R1_p * Ab_R0_p / Vp + koffAb * Ab_RR_p -
                       kon2Ab * S1_p * Ab_R0_p / Vp + koffAb * Ab_RS_p

    d/dt(Ab_RR_p) <-  kon1Ab * R1_p * Ab_0R_p / Vp - koffAb * Ab_RR_p +
                       kon2Ab * R1_p * Ab_R0_p / Vp - koffAb * Ab_RR_p -
                       kclearR1 * Ab_RR_p - kclearR1 * Ab_RR_p

    d/dt(Ab_0S_p) <-  kon2Ab * S1_p * Ab_00_p / Vp - koffAb * Ab_0S_p -
                       kclearAb * Ab_0S_p +
                       kout_Ab * Ab_0S_c - kin_Ab * Ab_0S_p -
                       kon1Ab * R1_p * Ab_0S_p / Vp + koffAb * Ab_RS_p -
                       kon1Ab * S1_p * Ab_0S_p / Vp + koffAb * Ab_SS_p

    d/dt(Ab_S0_p) <-  kon1Ab * S1_p * Ab_00_p / Vp - koffAb * Ab_S0_p -
                       kclearAb * Ab_S0_p +
                       kout_Ab * Ab_S0_c - kin_Ab * Ab_S0_p -
                       kon2Ab * R1_p * Ab_S0_p / Vp + koffAb * Ab_SR_p -
                       kon2Ab * S1_p * Ab_S0_p / Vp + koffAb * Ab_SS_p

    d/dt(Ab_RS_p) <-  kon2Ab * S1_p * Ab_R0_p / Vp - koffAb * Ab_RS_p +
                       kon1Ab * R1_p * Ab_0S_p / Vp - koffAb * Ab_RS_p -
                       kclearR1 * Ab_RS_p

    d/dt(Ab_SR_p) <-  kon2Ab * R1_p * Ab_S0_p / Vp - koffAb * Ab_SR_p +
                       kon1Ab * S1_p * Ab_0R_p / Vp - koffAb * Ab_SR_p -
                       kclearR1 * Ab_SR_p

    d/dt(Ab_SS_p) <-  kon1Ab * S1_p * Ab_0S_p / Vp - koffAb * Ab_SS_p +
                       kon2Ab * S1_p * Ab_S0_p / Vp - koffAb * Ab_SS_p -
                       kclearAb * Ab_SS_p +
                       kout_Ab * Ab_SS_c - kin_Ab * Ab_SS_p

    d/dt(depot) <- -ka * depot

    R1_c(0)  <- total_R1_c
    R1_p(0)  <- total_R1_p
    S1_c(0)  <- S1_0_c
    S1_p(0)  <- S1_0_p

    Cc <- Ab_00_c / Vc
    Cc ~ prop(propSd)
  })
}
