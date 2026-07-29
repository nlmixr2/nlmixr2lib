Marcantonio_2022_trastuzumab <- function() {
  description <- paste(
    "QSP. Two-compartment monospecific anti-receptor mechanistic PKPD model of",
    "trastuzumab-HER2 binding in adults with HER2-overexpressing breast cancer",
    "(Marcantonio 2022 Early Feasibility Assessment, Case Study 2/3). Drug",
    "administered IV as a bivalent (valency = 2) antibody that binds",
    "membrane HER2 (in central and peripheral compartments) and soluble HER2",
    "(the shed extracellular domain) via independent binding events with",
    "identical Kd. All species eliminate first-order; drug bound to membrane",
    "receptor eliminates at the receptor's rate; drug bound to soluble",
    "receptor eliminates at the drug's rate. Only soluble species (Ab_00,",
    "Ab_0S, Ab_S0, Ab_SS, S1) distribute between central and peripheral;",
    "membrane-bound drug forms (Ab_0R, Ab_R0, Ab_RR, Ab_RS, Ab_SR) stay in",
    "their compartment. The endogenous HER2 ligand is not resolved (paper",
    "sets its concentration to effectively zero) so L1 and L1R1 species are",
    "omitted here. Parameters FIXED from paper Table S7 and JSON run file;",
    "disease and toxicity compartments are omitted per paper Case Study 2",
    "text.",
    sep = " "
  )
  reference <- paste(
    "Marcantonio DH et al. (2022). Front Pharmacol 13:864768.",
    "doi:10.3389/fphar.2022.864768. Case Study 3 (trastuzumab HER2, breast",
    "cancer). Drug-specific parameters from paper Supplementary Table S7",
    "and the Assess run file JSON (Data Sheet 2).",
    sep = " "
  )
  vignette <- "Marcantonio_2022_efa"
  units <- list(
    time          = "day",
    dosing        = "Trastuzumab dose amount into Ab_00_c (IV bolus) in nmol; MW = 145531.5 Da so 140 mg = 962.0 nmol.",
    concentration = "Free trastuzumab plasma concentration Cc = Ab_00_c / Vc in nM; central volume Vc = 3 L, peripheral Vp = 13 L. Target engagement = drug-bound R1 / (drug-bound R1 + free R1)."
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = NA_integer_,
    disease_state  = "Adults with metastatic HER2-overexpressing breast cancer.",
    dose_range     = "For metastatic HER2-overexpressing breast cancer: initial 4 mg/kg IV loading, then 2 mg/kg IV weekly (Herceptin USPI). Model prediction targets the weekly maintenance dose in a 70 kg patient.",
    regions        = NA_character_,
    notes          = "See Marcantonio 2022 anti-ligand siblings for shared bottom-up target-expression methodology. Soluble HER2 (shed) concentration is 7 nM per Baselga 2001 (patient plasma levels)."
  )

  ini({
    lka       <- fixed(log(log(2) / 2.5));         label("First-order SC absorption placeholder rate constant (1/day; unused for IV dosing)")   # Assess default
    lkclearAb <- fixed(log(log(2) / 25));          label("First-order trastuzumab elimination rate constant (1/day; from t1/2 = 25 days)")     # Marcantonio 2022 Table S7 Drug Half Life (Quartino 2019 popPK)
    kd_drug   <- fixed(0.1);                       label("Trastuzumab-HER2 equilibrium dissociation constant (nM)")                             # Marcantonio 2022 Table S7 Trastuzumab KD (Li 2013, Troise 2008)
    valency   <- fixed(2);                         label("Effective drug valency for HER2 binding (dimensionless)")                              # Marcantonio 2022 Table S7 Drug Valency (assumed mAb)

    R1_conc_c    <- fixed(3.60e-4);                label("HER2 membrane receptor concentration in central compartment (nM)")                    # Marcantonio 2022 Table S7 HER2 expression central (You 2008, Suzuki 2015; bottom-up)
    R1_conc_p    <- fixed(5.89e-2);                label("HER2 membrane receptor concentration in peripheral compartment (nM)")                 # Marcantonio 2022 Table S7 HER2 expression peripheral (Li 2019, Onsum 2013, Ruiz 2018; bottom-up)
    S1_conc_c    <- fixed(7);                      label("Soluble shed HER2 concentration in central compartment (nM)")                          # Marcantonio 2022 Table S7 soluble HER2 receptor expression central (Baselga 2001)
    S1_conc_p    <- fixed(7);                      label("Soluble shed HER2 concentration in peripheral compartment (nM)")                       # Marcantonio 2022 Table S7 soluble HER2 receptor expression peripheral (assumed equal to central)

    lkclearR1    <- fixed(log(log(2) / (24 / 24)));    label("First-order HER2 membrane elimination rate constant (1/day; from t1/2 = 24 hr)")   # Marcantonio 2022 Table S7 HER2 receptor half-life (Pereira 2018 in vitro)
    lkclearS1    <- fixed(log(log(2) / (1 / 24)));     label("First-order soluble HER2 elimination rate constant (1/day; from t1/2 = 1 hr)")     # Marcantonio 2022 Table S7 soluble HER2 half-life (Betts 2020, Li 2014 PBPK estimates)

    Vc          <- fixed(3);                       label("Central compartment volume (L)")                                                        # Marcantonio 2022 Table S7 Central compartment volume (plasma volume; Shah 2012)
    Vp          <- fixed(13);                      label("Peripheral compartment volume (L)")                                                     # Marcantonio 2022 Table S7 Peripheral compartment volume (interstitial peripheral tissues; Shah 2012)
    Tdist_Ab_hr <- fixed(35);                      label("Drug distribution half-life to peripheral (hr)")                                        # Marcantonio 2022 JSON Tdist_Ab_hr_peripheral = 35 hr (typical mAb; Betts 2018)
    Pdist_Ab    <- fixed(0.190625);                label("Drug partition coefficient (peripheral vs central; unitless)")                          # Marcantonio 2022 JSON Pdist_Ab_peripheral = 0.190625 (typical mAb; Betts 2018)
    Tdist_S1_hr <- fixed(30);                      label("Soluble receptor distribution half-life to peripheral (hr)")                            # Marcantonio 2022 JSON Tdist_S1_hr_peripheral = 30 hr

    kon         <- fixed(0.001 * 86400);           label("Bimolecular association rate constant (L/nmol/day)")                                    # Assess default (0.001 nM-1 s-1)

    propSd    <- fixed(0.01);                      label("Placeholder proportional residual error on Cc (not paper-derived)")                    # Placeholder
  })

  model({
    ka         <- exp(lka)
    kclearAb   <- exp(lkclearAb)
    kclearR1   <- exp(lkclearR1)
    kclearS1   <- exp(lkclearS1)
    kon1Ab     <- kon
    kon2Ab     <- floor(valency / 2) * kon
    koffAb     <- kon * kd_drug

    # Drug distribution rate constants (Assess run file relationships block).
    # Tdist in HOURS -> convert to DAYS.
    Tdist_Ab_d <- Tdist_Ab_hr / 24
    kout_Ab    <- (log(2) / Tdist_Ab_d) * Pdist_Ab / (Pdist_Ab + Vc / Vp)
    kin_Ab     <- (log(2) / Tdist_Ab_d) / (1 + Pdist_Ab * Vp / Vc)

    # Soluble-receptor distribution assumes symmetric partition (Pdist_S1 = 1).
    Tdist_S1_d <- Tdist_S1_hr / 24
    QS1        <- log(2) / Tdist_S1_d
    Pdist_S1   <- S1_conc_p / (S1_conc_c + 1e-16)
    kout_S1    <- QS1 * Pdist_S1 / (Pdist_S1 + Vc / Vp)
    kin_S1     <- QS1 / (1 + Pdist_S1 * Vp / Vc)

    # Baseline (steady-state) values -- no ligand-receptor complex since L1
    # concentration is effectively zero for trastuzumab HER2 (no explicit
    # endogenous ligand modeled).
    total_R1_c <- R1_conc_c * Vc
    total_R1_p <- R1_conc_p * Vp
    S1_0_c     <- S1_conc_c * Vc
    S1_0_p     <- S1_conc_p * Vp

    # Shedding rate constants: at SS, R1 loss to S1 balances S1 clearance.
    kshed_R1_c <- kclearS1 * S1_0_c / total_R1_c
    kshed_R1_p <- kclearS1 * S1_0_p / total_R1_p

    # Synthesis rates: at SS with zero drug, R1 balance requires
    # ksynth = (kclear + kshed) * R1_0 (no ligand consumption).
    ksynth_R1_c <- (kclearR1 + kshed_R1_c) * total_R1_c
    ksynth_R1_p <- (kclearR1 + kshed_R1_p) * total_R1_p

    # ---------- CENTRAL COMPARTMENT ----------
    # NOTE: Assess source lists two "Ab_RR -> R1, kclearR1" reactions, so R1
    # gains 2*kclearR1*Ab_RR from Ab_RR clearance (each arm of Ab_RR releases
    # one R1 back to the free pool when internalized).
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

    # ---------- PERIPHERAL COMPARTMENT ----------
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

    # Initial conditions
    R1_c(0)  <- total_R1_c
    R1_p(0)  <- total_R1_p
    S1_c(0)  <- S1_0_c
    S1_p(0)  <- S1_0_p

    # Observations. Cc = free trastuzumab concentration in central (nM).
    Cc <- Ab_00_c / Vc
    Cc ~ prop(propSd)
  })
}
