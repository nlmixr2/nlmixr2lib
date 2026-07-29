Harrold_2020_filgrastim <- function() {
  description <- "Semi-mechanistic population PK / absolute-neutrophil-count / overall-survival model for subcutaneous filgrastim treatment of hematopoietic syndrome of acute radiation syndrome (HS-ARS) in adult and pediatric humans. PK is one-compartment (subcutaneous depot -> central drug amount) with target-mediated disposition through quadratic-equilibrium free / bound filgrastim partitioning against the time-varying G-CSF receptor pool. PD is a 5-stage granulopoiesis cascade (progenitor stem -> mitotic stem -> two precursor stages -> circulating neutrophils); bound drug stimulates receptor production (ST1) and transit between bone-marrow stages (ST2). Acute radiation effect is a kinetic-pharmacodynamic depot (depot_kpd) seeded by the radiation dose in Gy that decays first-order at rate kpde and kills the mitotic-stem stage at rate kpdkill * kpd ^ gamma; gamma depends on the radiation dose rate via a Hill-type function gamma = tgamma * DR / (DR + dr50). Overall survival is integrated as a Cox cumulative hazard (cumhaz_os) on a Box-Cox transformation of an effect-compartment ANC. All structural and IIV parameters fixed at the values from Harrold 2020 Table 2 (granulopoiesis values from Melhem 2018 popPK / ANC in healthy adults and chemotherapy-induced neutropenia; radiation / OS values scaled from rhesus-macaque NHP study with kpde and kpdkill multiplied by 0.72 to match the human LD50 and the kpdkill IIV omega halved per Harrold 2020 Methods 1.3 to address NHP-data sparsity)."
  reference <- paste(
    "Harrold JM, Olsson Gisleskog P, Perez-Ruixo JJ, Delor I, Chow A,",
    "Jacqmin P, Melhem M. (2020). Prediction of Survival Benefit of",
    "Filgrastim in Adult and Pediatric Patients With Acute Radiation",
    "Syndrome. Clin Transl Sci 13(4):807-817.",
    "doi:10.1111/cts.12777.",
    "Granulopoiesis sub-model parameters from Melhem M, et al. (2018)",
    "Br J Clin Pharmacol 84:911-925 (doi:10.1111/bcp.13504);",
    "radiation / overall-survival sub-model parameters scaled from the",
    "rhesus-macaque NHP analysis of Harrold J et al. (2015)",
    "J Pharmacokinet Pharmacodyn 42:S47-S48.",
    sep = " "
  )
  vignette <- "Harrold_2020_filgrastim_ars"
  units <- list(
    time          = "hour",
    dosing        = "nmol",
    concentration = "nmol/L",
    notes         = "Filgrastim dose in nmol to compartment 'depot' (1 ug = 1e-6/18800 mol = 53.2 pmol for filgrastim; molecular weight 18.8 kDa). Radiation dose in Gy delivered as a bolus to compartment 'depot_kpd'. Filgrastim concentration Cc reported in nM. ANC reported in cells/uL (= 10^9 cells/L when divided by 1000)."
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric size descriptor scaling filgrastim volume of distribution (exponent 0.943) and clearance (exponent 0.641) with the reference weight 70 kg. For pediatric simulations Harrold 2020 Methods 1.3 (Eq. 3) derives WT from age via Luscombe 2011 (APLS / 'Weight = 3*age + 7') for ages 1-16 years; for adult simulations the paper uses a Normal(mean = 70, SD = 15) distribution truncated to 45-125 kg (Methods 1.3).",
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1000L,
    n_studies      = 1L,
    age_range      = "1-16 years (pediatric subgroups 1-<6, 6-<12, 12-<16) plus adults; adult cohort defined by weight only (45-125 kg)",
    weight_range   = "45-125 kg (adults); 10-25, 25-43, and 43-55 kg for the three pediatric subgroups per Table S1",
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state  = "Adult and pediatric humans at risk of hematopoietic syndrome of acute radiation syndrome (HS-ARS) after acute whole-body ionising-radiation exposure. The model is a simulation framework: human ANC data after radiation are unavailable, so radiation / overall-survival parameters are scaled from a published rhesus-macaque NHP study (Harrold 2015) and calibrated against the historical human LD50 mortality curve of Scott & Dillehay 1990.",
    dose_range     = "Filgrastim subcutaneous 5, 7.5, 10, or 15 ug/kg once daily for 1-5 weeks, starting 1-21 days after radiation. Base scenario: 5 ug/kg q.d. for 28 days starting 1 day after radiation. Radiation dose 3-10.2 Gy delivered as a single bolus at dose rates 0.01-1000 Gy/h (base scenario 3.07 Gy at 1 Gy/h, the human LD50).",
    regions        = NA_character_,
    notes          = "1000 virtual subjects per arm in the base scenario (Methods 1.3). The granulopoiesis sub-model (Melhem 2018) was informed by healthy adult volunteers (75-750 ug or 5 ug/kg), adult chemotherapy patients (5 ug/kg), and pediatric chemotherapy patients (5, 10, or 15 ug/kg); the radiation / OS sub-model (Harrold 2015) was informed by rhesus-macaque NHPs receiving 10 ug/kg filgrastim after a lethal acute irradiation. Demographics (sex, race) are not reported in Harrold 2020 because the model is a simulation framework rather than a re-analysis of patient-level data."
  )

  ini({
    # ------------------------------------------------------------------
    # PK structural parameters -- Harrold 2020 Table 2 (granulopoiesis
    # sub-model values inherited from Melhem 2018; the same values are
    # restated in supplement Table S2 under 'Replicate Level / Variability').
    # All parameters are simulation inputs and held fixed in this library
    # implementation (Harrold 2020 Methods 1.1 describes the work as a
    # simulation built on the previously developed and validated
    # Melhem 2018 and Harrold 2015 components).
    # ------------------------------------------------------------------
    lfsc       <- fixed(log(1))     ; label("Relative bioavailability F_SC after s.c. filgrastim (fixed at 1)")   # Harrold 2020 Table 2 (FSC FIL = 1)
    lksc       <- fixed(log(0.123)) ; label("Subcutaneous absorption rate constant K_SC (1/h)")                    # Harrold 2020 Table 2 (KSC FIL = 0.123)
    lvd        <- fixed(log(3.12))  ; label("Filgrastim volume of distribution VD (L) at WT = 70 kg")              # Harrold 2020 Table 2 (VD FIL = 3.12 L)
    e_wt_vd    <- fixed(0.943)      ; label("Allometric exponent on VD with reference 70 kg")                      # Harrold 2020 Table 2 (beta_VD(WT/70) = 0.943)
    lcld       <- fixed(log(0.833)) ; label("Filgrastim non-receptor (renal + linear) clearance CLD (L/h) at WT = 70 kg")  # Harrold 2020 Table 2 (CLD FIL = 0.833 L/h)
    e_wt_cld   <- fixed(0.641)      ; label("Allometric exponent on CLD with reference 70 kg")                     # Harrold 2020 Table 2 (beta_CLD(WT/70) = 0.641)
    lkd        <- fixed(log(0.0237)); label("Filgrastim / G-CSFR dissociation constant K_D (nM)")                  # Harrold 2020 Table 2 (KD FIL = 0.0237 nM)
    lkint      <- fixed(log(0.113)) ; label("G-CSFR / drug complex internalisation rate K_INT (1/h)")              # Harrold 2020 Table 2 (KINT PT = 0.113)
    lbsld      <- fixed(log(0.00299)); label("Baseline endogenous G-CSF concentration BSLD (nM)")                   # Harrold 2020 Table 2 (BSLD = 0.00299 nM)

    # ------------------------------------------------------------------
    # Granulopoiesis ODE rate constants -- Harrold 2020 Table 2.
    # ------------------------------------------------------------------
    lkp        <- fixed(log(0.0276)); label("G-CSFR production rate K_P (nM/h)")                                    # Harrold 2020 Table 2 (K_P = 0.0276 nM/h)
    lktr       <- fixed(log(0.0330)); label("Bone-marrow transit rate K_TR (1/h) between receptor compartments")    # Harrold 2020 Table 2 (K_TR = 0.0330)
    lkc        <- fixed(log(0.120)) ; label("Neutrophil elimination rate K_C (1/h) from blood into tissues")        # Harrold 2020 Table 2 (K_C = 0.120)
    lsr        <- fixed(log(0.0590)); label("Receptor / ANC scaling factor S_R (nM per 1000 cells/uL units)")       # Harrold 2020 Table 2 (S_R = 0.0590)
    stm1       <- fixed(7.53)       ; label("Stimulation coefficient STM1 on receptor production (ST1 multiplier)")  # Harrold 2020 Table 2 (STM1 = 7.53)
    stm2       <- fixed(3.89)       ; label("Stimulation coefficient STM2 on transit rate (ST2 multiplier)")         # Harrold 2020 Table 2 (STM2 PT = 3.89)

    # ------------------------------------------------------------------
    # Acute-radiation effect parameters -- Harrold 2020 Table 2 and
    # Methods 1.2 / 1.3. K_PDE and K_PDKILL are reported as the
    # NHP-derived values 0.0141 and 0.425 in Table 2 footnotes a;
    # supplement Table S2 makes the *0.72 NHP-to-human scaling explicit
    # by setting TKPDE = 0.0141*0.72 and TKPDKILL = 0.425*0.72. The
    # values used in this model file are the scaled human-calibrated
    # values per Methods 1.2.
    # ------------------------------------------------------------------
    lkpde      <- fixed(log(0.0141 * 0.72)); label("Radiation-effect elimination rate K_PDE (1/h) -- NHP value 0.0141 scaled by 0.72")     # Harrold 2020 Table 2 (K_PDE = 0.0141) and Methods 1.2 (scaling factor 0.72); supplement Table S2 TKPDE
    lkpdkill   <- fixed(log(0.425  * 0.72)); label("Radiation cell-kill rate K_PDKILL (1/h KPD^-gamma) -- NHP value 0.425 scaled by 0.72") # Harrold 2020 Table 2 (K_PDKILL = 0.425) and Methods 1.2 (scaling factor 0.72); supplement Table S2 TKPDKILL
    ltgamma    <- fixed(log(2.20))         ; label("Maximum radiation-sensitivity exponent TGAMMA (unitless)")                            # Harrold 2020 Table 2 footnote c and Eq. 2 (TGAMMA = 2.20)
    ldr50      <- fixed(log(0.028))        ; label("Dose rate at 50% of TGAMMA (Gy/h)")                                                   # Harrold 2020 Table 2 footnote c and Eq. 2 (dose rate at half-max = 0.028 Gy/h)
    ldr        <- fixed(log(1))            ; label("Radiation dose rate DR (Gy/h); base scenario = 1, override per simulation")            # Harrold 2020 Methods 1.3 base scenario (DR = 1 Gy/h); see Table 1 for other DR-scenario values

    # ------------------------------------------------------------------
    # Overall-survival sub-model parameters -- Harrold 2020 Table 2
    # (carried from Harrold 2015 NHP analysis).
    # lambda_anc is the linear slope on the Box-Cox transformed ANCe
    # (negative because higher ANC -> lower hazard); lambda_bc is the
    # Box-Cox power parameter; both are bare (not log-transformed)
    # because the published values are negative.
    # ------------------------------------------------------------------
    lke0       <- fixed(log(0.0278)); label("Effect-compartment equilibration rate k_e0 for ANCe (1/h)")  # Harrold 2020 Table 2 (k_e0 = 0.0278)
    lambda_anc <- fixed(-2.14)      ; label("Hazard slope on Box-Cox(ANCe) -- negative")                  # Harrold 2020 Table 2 (lambda_ANC = -2.14)
    lambda_bc  <- fixed(-0.347)     ; label("Box-Cox power parameter on ANCe (unitless)")                 # Harrold 2020 Table 2 (lambda_BC = -0.347)

    # ------------------------------------------------------------------
    # Inter-individual variability -- Harrold 2020 Table 2 (Random
    # effects). All values are SDs on the log scale; variance = SD^2.
    # The kpdkill omega is divided by 2 per Methods 1.3 ("we reduced
    # variability by 50% (208%), which is closer to the IIV in CIN").
    # The correlation between K_PDE and K_PDKILL is 0.91 (Table 2);
    # off-diagonal covariance = 0.91 * sd_kpde * sd_kpdkill.
    # All IIV is fixed at the published values (simulation use).
    # ------------------------------------------------------------------
    etalfsc    ~ fixed(0.440^2)     # Harrold 2020 Table 2 (Omega FSC SD = 0.440)
    etalksc    ~ fixed(0.225^2)     # Harrold 2020 Table 2 (Omega KSC SD = 0.225)
    etalvd     ~ fixed(0.282^2)     # Harrold 2020 Table 2 (Omega VD SD = 0.282)
    etalcld    ~ fixed(0.370^2)     # Harrold 2020 Table 2 (Omega CLD SD = 0.370)
    etalkp     ~ fixed(0.265^2)     # Harrold 2020 Table 2 (Omega KP SD = 0.265)
    etalkd     ~ fixed(0.726^2)     # Harrold 2020 Table 2 (Omega KD SD = 0.726)
    etastm1    ~ fixed(0.315^2)     # Harrold 2020 Table 2 (Omega STM1 SD = 0.315) -- IIV on STM1 (receptor-production stimulation)
    etastm2    ~ fixed(0.273^2)     # Harrold 2020 Table 2 (Omega STM2 SD = 0.273) -- IIV on STM2 (transit stimulation)
    etalkint   ~ fixed(0.570^2)     # Harrold 2020 Table 2 (Omega KINT SD = 0.570)
    etalbsld   ~ fixed(0.260^2)     # Harrold 2020 Table 2 (Omega BSLD SD = 0.260)

    # Correlated IIV on K_PDE and K_PDKILL.
    # SD_kpde    = 0.314 (Harrold 2020 Table 2 Omega K_PDE)
    # SD_kpdkill = 4.16 / 2 = 2.08 (Harrold 2020 Table 2 Omega K_PDKILL with footnote b applied)
    # cov(kpde, kpdkill) = 0.91 * 0.314 * 2.08 = 0.5942
    etalkpde + etalkpdkill ~ fixed(c(0.314^2,
                                     0.91 * 0.314 * 2.08, 2.08^2))   # Harrold 2020 Table 2 (Omega K_PDE = 0.314, Omega K_PDKILL = 4.16 halved per footnote b, correlation 0.910)

    # ------------------------------------------------------------------
    # Residual error -- Harrold 2020 Table 2 ("Exponential residual
    # error model"; the source paper labels these as 'log domain'
    # residuals which in the linear-space simulation enter as
    # multiplicative noise on the predicted CONC and ANC). The
    # exploratory simulations in Methods 1.3 set the residuals to zero
    # ("No residual variability was applied to minimise random noise");
    # the published estimates are recorded here as fixed values and the
    # vignette demonstrates how to use rxode2::zeroRe() to suppress
    # residual variability for typical-value simulations.
    # ------------------------------------------------------------------
    expSd      <- fixed(0.537)      ; label("Residual SD on log-domain CONC (a1, PK)")    # Harrold 2020 Table 2 (a1 PK = 0.537)
    expSd_ANC  <- fixed(0.298)      ; label("Residual SD on log-domain ANC (a2, PD)")     # Harrold 2020 Table 2 (a2 PD = 0.298)
  })

  model({
    # ------------------------------------------------------------------
    # Individual parameters (apply log-normal IIV around the typical
    # values from Table 2). All structural and IIV parameters are
    # fixed; eta draws still vary subject to subject.
    # ------------------------------------------------------------------
    fsc      <- exp(lfsc      + etalfsc)
    ksc      <- exp(lksc      + etalksc)
    vd       <- exp(lvd       + etalvd)       * (WT / 70) ^ e_wt_vd      # Harrold 2020 Table 2 VD allometric on WT (ref 70 kg)
    cld      <- exp(lcld      + etalcld)      * (WT / 70) ^ e_wt_cld     # Harrold 2020 Table 2 CLD allometric on WT (ref 70 kg)
    kd       <- exp(lkd       + etalkd)
    kint     <- exp(lkint     + etalkint)
    bsld     <- exp(lbsld     + etalbsld)
    kp       <- exp(lkp       + etalkp)
    ktr      <- exp(lktr)                                                # KTR fixed (no IIV in Table 2)
    kc       <- exp(lkc)                                                 # KC fixed (no IIV in Table 2)
    sr       <- exp(lsr)                                                 # S_R fixed (no IIV in Table 2)
    stm1_i   <- stm1 * exp(etastm1)                                      # IIV applied multiplicatively on the linear-space STM1
    stm2_i   <- stm2 * exp(etastm2)                                      # IIV applied multiplicatively on the linear-space STM2
    kpde     <- exp(lkpde     + etalkpde)
    kpdkill  <- exp(lkpdkill  + etalkpdkill)
    tgamma   <- exp(ltgamma)
    dr50     <- exp(ldr50)
    dr       <- exp(ldr)                                                  # radiation dose rate (Gy/h); override per simulation
    gamma    <- tgamma * dr / (dr + dr50)                                 # Harrold 2020 Eq. 2 (Hill on dose rate)
    ke0      <- exp(lke0)

    # ------------------------------------------------------------------
    # Drug clearance from central. CLD is the non-receptor (renal +
    # linear) clearance term that operates on the FREE drug
    # concentration only (Harrold 2020 Figure 1: CL_D acts on the free
    # filgrastim arm; receptor-bound drug is cleared via internalisation
    # KINT on the complex). The effective first-order rate is
    # ke_drug = cld / vd applied to the free amount (FDC * VD), and the
    # receptor-internalisation term is kint applied to the bound amount
    # (RDC * VD).
    # ------------------------------------------------------------------
    ke_drug  <- cld / vd

    # ------------------------------------------------------------------
    # Free / bound filgrastim from quadratic binding equilibrium
    # (Harrold 2020 supplement Table S2, Structural Equations). At
    # equilibrium between total drug concentration tdc = central / vd
    # and total receptor concentration trc = circ (receptors track the
    # circulating neutrophil pool 1:1 in the model frame), the free
    # filgrastim concentration is the positive root of the binding
    # quadratic:
    #   fdc = 0.5 * (tdc - trc - kd + sqrt((tdc - trc - kd)^2 + 4*kd*tdc))
    # and the bound drug concentration is rdc = tdc - fdc.
    # ------------------------------------------------------------------
    tdc      <- central / vd
    trc      <- circ
    qbind    <- tdc - trc - kd
    fdc      <- 0.5 * (qbind + sqrt(qbind * qbind + 4 * kd * tdc))
    rdc      <- tdc - fdc

    # ------------------------------------------------------------------
    # Receptor-mediated stimulation of granulopoiesis (Harrold 2020
    # supplement Table S2). ST1 scales receptor production K_P; ST2
    # scales transit between bone-marrow compartments.
    # ------------------------------------------------------------------
    st1      <- 1 + stm1_i * (rdc / trc)
    st2      <- 1 + stm2_i * (rdc / trc)

    # ------------------------------------------------------------------
    # Radiation cell-kill term (Harrold 2020 supplement Table S2).
    # STM = kpdkill * kpd ^ gamma drives the additional loss of mitotic
    # stem cells (precursor2). The depot_kpd compartment carries the
    # decaying radiation effect; the radiation dose in Gy is delivered
    # as a bolus to depot_kpd at the radiation exposure event.
    # ------------------------------------------------------------------
    stm_rad  <- kpdkill * depot_kpd ^ gamma

    # ------------------------------------------------------------------
    # ODE system.
    #   depot      : s.c. depot for filgrastim (amount, nmol)
    #   central    : central filgrastim amount (nmol; total drug)
    #   precursor1 : SM, progenitor stem (concentration, nM)
    #   precursor2 : MT, mitotic stem (concentration, nM) -- radiation acts here
    #   precursor3 : PM1, precursor cells (concentration, nM)
    #   precursor4 : PM2, precursor cells (concentration, nM)
    #   circ       : RB, blood neutrophils (concentration, nM)
    #   depot_kpd  : K-PD radiation effect (Gy)
    #   effect     : ANCe effect compartment (cells/uL)
    #   cumhaz_os  : cumulative hazard for overall survival (unitless)
    # ------------------------------------------------------------------
    d/dt(depot)      <- -ksc * depot
    d/dt(central)    <-  ksc * depot - ke_drug * (fdc * vd) - kint * (rdc * vd)
    d/dt(precursor1) <-  kp * st1               - ktr * st2 * precursor1
    d/dt(precursor2) <-  ktr * st2 * precursor1 - (ktr * st2 + stm_rad) * precursor2
    d/dt(precursor3) <-  ktr * st2 * precursor2 - ktr * st2 * precursor3
    d/dt(precursor4) <-  ktr * st2 * precursor3 - ktr * st2 * precursor4
    d/dt(circ)       <-  ktr * st2 * precursor4 - kc * circ
    d/dt(depot_kpd)  <- -kpde * depot_kpd
    d/dt(effect)     <-  ke0 * (1000 * circ / sr - effect)
    d/dt(cumhaz_os)  <-  exp(lambda_anc * ((effect ^ lambda_bc - 1) / lambda_bc)) / 24

    # ------------------------------------------------------------------
    # Steady-state initial conditions (Harrold 2020 supplement Table S2
    # 'Model Parameters' block: SM = MT = PM1 = PM2 = KP/KTR; RB = KP/KC;
    # ANCE = 1000*(KP/KC)/SR). Per-individual ICs computed from the
    # individual kp / ktr / kc / sr values.
    # ------------------------------------------------------------------
    precursor1(0) <- kp / ktr
    precursor2(0) <- kp / ktr
    precursor3(0) <- kp / ktr
    precursor4(0) <- kp / ktr
    circ(0)       <- kp / kc
    effect(0)     <- 1000 * (kp / kc) / sr
    cumhaz_os(0)  <- 0

    # ------------------------------------------------------------------
    # Bioavailability of the s.c. depot (FSC = 1 typical; IIV via fsc).
    # ------------------------------------------------------------------
    f(depot)      <- fsc

    # ------------------------------------------------------------------
    # Observations.
    #   Cc          : serum filgrastim concentration in nM = FDC + BSLD
    #                 (free drug plus the endogenous G-CSF baseline,
    #                 Harrold 2020 supplement Table S2: CONC = FDC + BSLD).
    #   ANC         : circulating absolute neutrophil count in cells/uL
    #                 = 1000 * circ / sr (supplement Table S2: ANC = RB/SR).
    #   survival_os : Kaplan-Meier-equivalent survival probability
    #                 S(t) = exp(-cumhaz_os) per the Cox hazard model.
    # ------------------------------------------------------------------
    Cc          <- fdc + bsld
    Cc          ~ lnorm(expSd)
    ANC         <- 1000 * circ / sr
    ANC         ~ lnorm(expSd_ANC)
    survival_os <- exp(-cumhaz_os)
  })
}
