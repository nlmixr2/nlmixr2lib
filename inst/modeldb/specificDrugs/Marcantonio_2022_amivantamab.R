Marcantonio_2022_amivantamab <- function() {
  description <- paste(
    "QSP. Two-compartment bispecific anti-receptor x anti-receptor",
    "mechanistic PKPD model of amivantamab-(EGFR, c-Met) binding in adults",
    "with EGFR-exon-20-mutant non-small-cell lung cancer (Marcantonio 2022",
    "Early Feasibility Assessment, Case Study 2). Amivantamab has one arm",
    "for EGFR (Kd 1.4 nM) and one arm for c-Met (Kd 0.04 nM) with valency 1",
    "per target (per Jarantow 2015). Model tracks membrane EGFR (R1),",
    "membrane c-Met (R2), and soluble c-Met (S2) at 5.9 nM in both central",
    "and peripheral compartments. Soluble EGFR is NOT modelled",
    "(shed_css_1 = 0 in JSON). Drug states: free drug (Ab_00) plus five",
    "populated bound forms Ab_R1, Ab_R2, Ab_S2, Ab_R1R2, Ab_R1S2 per",
    "compartment. Since each arm is monovalent, the states with double-",
    "target arm binding (Ab_R1R1, Ab_R2R2, Ab_S2S2) do not populate and are",
    "omitted. Parameters FIXED from paper Tables 3 (target-side EGFR/c-Met)",
    "and 4 (drug-side) plus the Assess run file JSON; disease and toxicity",
    "compartments omitted per paper Case Study 2 text.",
    sep = " "
  )
  reference <- paste(
    "Marcantonio DH et al. (2022). Front Pharmacol 13:864768.",
    "doi:10.3389/fphar.2022.864768. Case Study 2 (amivantamab EGFR + c-Met",
    "bispecific, NSCLC). Drug-specific parameters from paper Tables 3 and",
    "4 and the Assess run file JSON (Data Sheet 2).",
    sep = " "
  )
  vignette <- "Marcantonio_2022_efa"
  units <- list(
    time          = "day",
    dosing        = "Amivantamab dose amount into Ab_00_c (IV bolus) in nmol; MW = 150000 Da so 1050 mg = 7000 nmol.",
    concentration = "Free amivantamab plasma concentration Cc = Ab_00_c / Vc in nM; central Vc = 3 L, peripheral Vp = 13 L."
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = NA_integer_,
    disease_state  = "Adults with locally advanced or metastatic non-small-cell lung cancer with EGFR exon-20 insertion mutations.",
    dose_range     = "Clinically approved dose: 1050 mg IV for < 80 kg patients (1400 mg if >= 80 kg) weekly for the first 4 weeks then every 2 weeks (Rybrevant USPI). Model prediction targets 326 mg Q1W or 740 mg Q2W for 98% target engagement of both targets.",
    regions        = NA_character_,
    notes          = "See Marcantonio 2022 trastuzumab and panitumumab siblings for shared 2-cpt anti-receptor structure."
  )

  ini({
    lka       <- fixed(log(log(2) / 2.5));         label("First-order SC absorption placeholder rate constant (1/day; unused for IV dosing)")     # Assess default
    lkclearAb <- fixed(log(log(2) / 11));          label("First-order amivantamab elimination rate constant (1/day; from t1/2 = 11 days)")        # Marcantonio 2022 Table 4 Amivantamab Half-Life (Rybrevant USPI)
    kd_T1     <- fixed(1.4);                       label("Amivantamab-EGFR equilibrium dissociation constant (nM)")                                # Marcantonio 2022 Table 4 Amivantamab KD for EGFR (Jarantow 2015)
    kd_T2     <- fixed(0.04);                      label("Amivantamab-cMet equilibrium dissociation constant (nM)")                                # Marcantonio 2022 Table 4 Amivantamab KD for c-Met (Jarantow 2015)

    R1_conc_c <- fixed(0.0152);                    label("EGFR receptor concentration in central compartment (nM)")                                # Marcantonio 2022 Table 3 EGFR expression central (bottom-up)
    R1_conc_p <- fixed(1.13);                      label("EGFR receptor concentration in peripheral compartment (nM)")                             # Marcantonio 2022 Table 3 EGFR expression peripheral (bottom-up)
    R2_conc_c <- fixed(0.011);                     label("c-Met receptor concentration in central compartment (nM)")                               # Marcantonio 2022 Table 3 Met expression central (bottom-up)
    R2_conc_p <- fixed(0.45);                      label("c-Met receptor concentration in peripheral compartment (nM)")                            # Marcantonio 2022 Table 3 Met expression peripheral (bottom-up)
    S2_conc_c <- fixed(5.9);                       label("Soluble c-Met concentration in central compartment (nM)")                                # Marcantonio 2022 Table 3 soluble Met concentration (Rosen 2017, Gao 2016)
    S2_conc_p <- fixed(5.9);                       label("Soluble c-Met concentration in peripheral compartment (nM)")                             # Marcantonio 2022 Table 3 (assumed equal to central)

    lkclearR1 <- fixed(log(log(2) / (5 / 24)));    label("First-order EGFR membrane elimination rate constant (1/day; from t1/2 = 5 hr)")          # Marcantonio 2022 Table 3 EGFR receptor half-life (Sigismund 2008)
    lkclearR2 <- fixed(log(log(2) / (4 / 24)));    label("First-order c-Met membrane elimination rate constant (1/day; from t1/2 = 4 hr)")        # Marcantonio 2022 Table 3 Met receptor half-life (Li 2008, DaSilva 2020)
    lkclearS2 <- fixed(log(log(2) / (48 / 24)));   label("First-order soluble c-Met elimination rate constant (1/day; from t1/2 = 48 hr)")         # Marcantonio 2022 Table 3 soluble Met half-life (Li 2017 protein MW estimate)

    Vc          <- fixed(3);                       label("Central compartment volume (L)")                                                          # Marcantonio 2022 Table 3
    Vp          <- fixed(13);                      label("Peripheral compartment volume (L)")                                                       # Marcantonio 2022 Table 3
    Tdist_Ab_hr <- fixed(35);                      label("Drug distribution half-life to peripheral (hr)")                                          # Marcantonio 2022 Table 4 (Betts 2018)
    Pdist_Ab    <- fixed(0.190625);                label("Drug partition coefficient (unitless)")                                                    # Marcantonio 2022 Table 4 (Betts 2018)
    Tdist_S2_hr <- fixed(30);                      label("Soluble c-Met distribution half-life to peripheral (hr)")                                  # Marcantonio 2022 JSON default

    kon         <- fixed(0.001 * 86400);           label("Bimolecular association rate constant (L/nmol/day)")                                      # Assess default

    propSd    <- fixed(0.01);                      label("Placeholder proportional residual error on Cc (not paper-derived)")                      # Placeholder
  })

  model({
    ka         <- exp(lka)
    kclearAb   <- exp(lkclearAb)
    kclearR1   <- exp(lkclearR1)
    kclearR2   <- exp(lkclearR2)
    kclearS2   <- exp(lkclearS2)

    koff_T1    <- kon * kd_T1
    koff_T2    <- kon * kd_T2

    # Drug distribution
    Tdist_Ab_d <- Tdist_Ab_hr / 24
    kout_Ab    <- (log(2) / Tdist_Ab_d) * Pdist_Ab / (Pdist_Ab + Vc / Vp)
    kin_Ab     <- (log(2) / Tdist_Ab_d) / (1 + Pdist_Ab * Vp / Vc)

    # Soluble c-Met distribution (Pdist_S2 = 1 since S2 is at 5.9 nM in both compartments)
    Tdist_S2_d <- Tdist_S2_hr / 24
    QS2        <- log(2) / Tdist_S2_d
    Pdist_S2   <- S2_conc_p / (S2_conc_c + 1e-16)
    kout_S2    <- QS2 * Pdist_S2 / (Pdist_S2 + Vc / Vp)
    kin_S2     <- QS2 / (1 + Pdist_S2 * Vp / Vc)

    # Baseline steady-state amounts
    total_R1_c <- R1_conc_c * Vc
    total_R1_p <- R1_conc_p * Vp
    total_R2_c <- R2_conc_c * Vc
    total_R2_p <- R2_conc_p * Vp
    S2_0_c     <- S2_conc_c * Vc
    S2_0_p     <- S2_conc_p * Vp

    # Shedding rates (soluble c-Met arises from R2 shedding at SS)
    kshed_R2_c <- kclearS2 * S2_0_c / total_R2_c
    kshed_R2_p <- kclearS2 * S2_0_p / total_R2_p

    # Synthesis rates to maintain baseline
    ksynth_R1_c <- kclearR1 * total_R1_c
    ksynth_R1_p <- kclearR1 * total_R1_p
    ksynth_R2_c <- (kclearR2 + kshed_R2_c) * total_R2_c
    ksynth_R2_p <- (kclearR2 + kshed_R2_p) * total_R2_p

    # ---------- CENTRAL COMPARTMENT ----------
    d/dt(R1_c) <-  ksynth_R1_c - kclearR1 * R1_c -
                    kon * R1_c * Ab_00_c / Vc + koff_T1 * Ab_R1_c -
                    kon * R1_c * Ab_R2_c / Vc + koff_T1 * Ab_R1R2_c -
                    kon * R1_c * Ab_S2_c / Vc + koff_T1 * Ab_R1S2_c +
                    kclearR2 * Ab_R1R2_c

    d/dt(R2_c) <-  ksynth_R2_c - kclearR2 * R2_c - kshed_R2_c * R2_c -
                    kon * R2_c * Ab_00_c / Vc + koff_T2 * Ab_R2_c -
                    kon * R2_c * Ab_R1_c / Vc + koff_T2 * Ab_R1R2_c +
                    kclearR1 * Ab_R1R2_c

    d/dt(S2_c) <-  kshed_R2_c * R2_c - kclearS2 * S2_c -
                    kout_S2 * S2_c + kin_S2 * S2_p -
                    kon * S2_c * Ab_00_c / Vc + koff_T2 * Ab_S2_c -
                    kon * S2_c * Ab_R1_c / Vc + koff_T2 * Ab_R1S2_c

    d/dt(Ab_00_c) <-  ka * depot -
                       kclearAb * Ab_00_c -
                       kout_Ab * Ab_00_c + kin_Ab * Ab_00_p -
                       kon * R1_c * Ab_00_c / Vc + koff_T1 * Ab_R1_c -
                       kon * R2_c * Ab_00_c / Vc + koff_T2 * Ab_R2_c -
                       kon * S2_c * Ab_00_c / Vc + koff_T2 * Ab_S2_c

    d/dt(Ab_R1_c) <-  kon * R1_c * Ab_00_c / Vc - koff_T1 * Ab_R1_c -
                       kclearR1 * Ab_R1_c -
                       kon * R2_c * Ab_R1_c / Vc + koff_T2 * Ab_R1R2_c -
                       kon * S2_c * Ab_R1_c / Vc + koff_T2 * Ab_R1S2_c

    d/dt(Ab_R2_c) <-  kon * R2_c * Ab_00_c / Vc - koff_T2 * Ab_R2_c -
                       kclearR2 * Ab_R2_c -
                       kon * R1_c * Ab_R2_c / Vc + koff_T1 * Ab_R1R2_c

    d/dt(Ab_S2_c) <-  kon * S2_c * Ab_00_c / Vc - koff_T2 * Ab_S2_c -
                       kclearAb * Ab_S2_c -
                       kout_Ab * Ab_S2_c + kin_Ab * Ab_S2_p -
                       kon * R1_c * Ab_S2_c / Vc + koff_T1 * Ab_R1S2_c

    d/dt(Ab_R1R2_c) <- kon * R2_c * Ab_R1_c / Vc - koff_T2 * Ab_R1R2_c +
                        kon * R1_c * Ab_R2_c / Vc - koff_T1 * Ab_R1R2_c -
                        (kclearR1 + kclearR2) * Ab_R1R2_c

    d/dt(Ab_R1S2_c) <- kon * S2_c * Ab_R1_c / Vc - koff_T2 * Ab_R1S2_c +
                        kon * R1_c * Ab_S2_c / Vc - koff_T1 * Ab_R1S2_c -
                        kclearR1 * Ab_R1S2_c

    # ---------- PERIPHERAL COMPARTMENT ----------
    d/dt(R1_p) <-  ksynth_R1_p - kclearR1 * R1_p -
                    kon * R1_p * Ab_00_p / Vp + koff_T1 * Ab_R1_p -
                    kon * R1_p * Ab_R2_p / Vp + koff_T1 * Ab_R1R2_p -
                    kon * R1_p * Ab_S2_p / Vp + koff_T1 * Ab_R1S2_p +
                    kclearR2 * Ab_R1R2_p

    d/dt(R2_p) <-  ksynth_R2_p - kclearR2 * R2_p - kshed_R2_p * R2_p -
                    kon * R2_p * Ab_00_p / Vp + koff_T2 * Ab_R2_p -
                    kon * R2_p * Ab_R1_p / Vp + koff_T2 * Ab_R1R2_p +
                    kclearR1 * Ab_R1R2_p

    d/dt(S2_p) <-  kshed_R2_p * R2_p - kclearS2 * S2_p +
                    kout_S2 * S2_c - kin_S2 * S2_p -
                    kon * S2_p * Ab_00_p / Vp + koff_T2 * Ab_S2_p -
                    kon * S2_p * Ab_R1_p / Vp + koff_T2 * Ab_R1S2_p

    d/dt(Ab_00_p) <-  -kclearAb * Ab_00_p +
                       kout_Ab * Ab_00_c - kin_Ab * Ab_00_p -
                       kon * R1_p * Ab_00_p / Vp + koff_T1 * Ab_R1_p -
                       kon * R2_p * Ab_00_p / Vp + koff_T2 * Ab_R2_p -
                       kon * S2_p * Ab_00_p / Vp + koff_T2 * Ab_S2_p

    d/dt(Ab_R1_p) <-  kon * R1_p * Ab_00_p / Vp - koff_T1 * Ab_R1_p -
                       kclearR1 * Ab_R1_p -
                       kon * R2_p * Ab_R1_p / Vp + koff_T2 * Ab_R1R2_p -
                       kon * S2_p * Ab_R1_p / Vp + koff_T2 * Ab_R1S2_p

    d/dt(Ab_R2_p) <-  kon * R2_p * Ab_00_p / Vp - koff_T2 * Ab_R2_p -
                       kclearR2 * Ab_R2_p -
                       kon * R1_p * Ab_R2_p / Vp + koff_T1 * Ab_R1R2_p

    d/dt(Ab_S2_p) <-  kon * S2_p * Ab_00_p / Vp - koff_T2 * Ab_S2_p -
                       kclearAb * Ab_S2_p +
                       kout_Ab * Ab_S2_c - kin_Ab * Ab_S2_p -
                       kon * R1_p * Ab_S2_p / Vp + koff_T1 * Ab_R1S2_p

    d/dt(Ab_R1R2_p) <- kon * R2_p * Ab_R1_p / Vp - koff_T2 * Ab_R1R2_p +
                        kon * R1_p * Ab_R2_p / Vp - koff_T1 * Ab_R1R2_p -
                        (kclearR1 + kclearR2) * Ab_R1R2_p

    d/dt(Ab_R1S2_p) <- kon * S2_p * Ab_R1_p / Vp - koff_T2 * Ab_R1S2_p +
                        kon * R1_p * Ab_S2_p / Vp - koff_T1 * Ab_R1S2_p -
                        kclearR1 * Ab_R1S2_p

    d/dt(depot) <- -ka * depot

    R1_c(0) <- total_R1_c
    R1_p(0) <- total_R1_p
    R2_c(0) <- total_R2_c
    R2_p(0) <- total_R2_p
    S2_c(0) <- S2_0_c
    S2_p(0) <- S2_0_p

    Cc <- Ab_00_c / Vc
    Cc ~ prop(propSd)
  })
}
