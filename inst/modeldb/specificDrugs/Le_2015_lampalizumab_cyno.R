Le_2015_lampalizumab_cyno <- function() {
  description <- "Preclinical (cynomolgus monkey). Semi-mechanistic and integrated ocular-systemic target-mediated drug-disposition (TMDD) PK/PD model of lampalizumab (anti-complement factor D antigen-binding fragment) and complement factor D (CFD) after intravenous or bilateral intravitreal (ITV) administration. Uses the quasi-steady-state (QSS) approximation of Gibiansky 2008. States: vitreous total lampalizumab amount (depot), serum total lampalizumab amount (central), a lampalizumab peripheral distribution compartment (peripheral1), and total CFD concentrations in vitreous (total_target), serum (total_target_central) and a target peripheral distribution compartment (total_target_peripheral1). Total lampalizumab concentrations in aqueous humor and retina are algebraic partitions of vitreous total drug. The model represents one eye; vitreous-to-serum flows in the serum ODEs carry a factor of 2 to reflect bilateral ITV dosing. Naive-pool fit (no IIV); residual error was reported qualitatively (proportional) but no numeric variance is given in the paper, so residual SDs are encoded as fixed(0) so deterministic simulation reproduces the paper exactly (see vignette Errata)."
  reference <- "Le KN, Gibiansky L, Good J, Davancaze T, van Lookeren Campagne M, Loyet KM, Morimoto A, Jin JY, Damico-Beyer LA, Hanley WD. A Mechanistic Pharmacokinetic/Pharmacodynamic Model of Factor D Inhibition in Cynomolgus Monkeys by Lampalizumab for the Treatment of Geographic Atrophy. J Pharmacol Exp Ther. 2015;355(2):288-296. doi:10.1124/jpet.115.227223. PMID: 26359312."
  vignette <- "Le_2015_lampalizumab_cyno"
  units <- list(time = "day", dosing = "nmol", concentration = "nmol/mL (= uM); CFD concentrations in nmol/mL")

  # Paper-specific compartments -- serum and peripheral total-CFD concentration
  # states extend the canonical `total_target` (vitreous) with two additional
  # compartments. Not covered by targetLocationRegex (which allows only
  # csf / isf / peripheralN suffixes on `target` / `complex`).
  paper_specific_compartments <- c(
    "total_target_central",
    "total_target_peripheral1"
  )

  covariateData <- list()

  population <- list(
    species        = "cynomolgus monkey (Macaca fascicularis)",
    n_subjects     = 116L,
    n_studies      = 3L,
    weight_range   = "2-4 kg (paper Results: serum V_c 127 mL is ~32-64 mL/kg for a 2-4 kg animal)",
    sex_female_pct = NULL,
    disease_state  = "Healthy cynomolgus monkeys (no disease model).",
    dose_range     = "Single i.v. 0.2, 2, or 20 mg per animal (study 08-1021, n=9). Single bilateral ITV 1, 2, 3, 5, 10, or 20 mg per eye across studies 08-1021 (n=20), 08-0782 (n=8), and 08-0783 (n=41). Studies conducted at Covance (Madison, WI).",
    regions        = "Preclinical, US site (Covance, Madison, WI).",
    notes          = "Naive-pool analysis (no IIV). Data from three GLP studies (08-1021, 08-0782, 08-0783). Sample matrices: serum, vitreous humor, aqueous humor, and retinal tissue (retinal tissue results normalised by tissue weight). Both males and females included; per-study sex distributions in Le 2015 Table 1. See tracking/operator_followups.md if a sex-stratified breakdown is later required.",
    trial_ids      = "Genentech internal: 08-1021, 08-0782, 08-0783."
  )

  ini({
    # ------------------------------------------------------------------
    # Ocular structural parameters. Le 2015 Table 2 ("Ocular parameters"
    # block). All time constants in day^-1; volumes in mL; concentrations
    # in nmol/mL (which equals uM).
    # ------------------------------------------------------------------
    lVVITR      <- log(2.2);          label("Vitreous volume V_VITR (mL)")                                                        # Le 2015 Table 2 row 'V_VITR' (2.2 mL, RSE 7.5%)
    lkout       <- log(0.24);         label("Ocular elimination rate constant of unbound lampalizumab k_out (1/day)")             # Le 2015 Table 2 row 'k_out' (0.24 /day, RSE 2.1%)
    lkoutC      <- log(0.75);         label("Ocular elimination rate constant of lampalizumab-CFD complex k_outC (1/day)")        # Le 2015 Table 2 row 'k_outC' (0.75 /day, RSE 17%)
    lkoutT      <- fixed(log(0.30));  label("Ocular elimination rate constant of free CFD k_outT (1/day; FIXED)")                  # Le 2015 Table 2 row 'k_outT' (0.30 /day FIXED to in-house vitreous CFD t1/2 data)
    kinC        <- fixed(0);          label("Ocular influx rate constant of complex k_inC (1/day; FIXED to 0)")                    # Le 2015 Table 2 row 'K_inC' (0 FIXED; complex influx assumed negligible)
    llambda_VA  <- log(4.4);          label("Vitreous-to-aqueous partition coefficient for lampalizumab lambda_VA (unitless)")     # Le 2015 Table 2 row 'lambda_VA' (4.4, RSE 12%)
    llambda_TVA <- log(0.47);         label("Vitreous-to-aqueous partition coefficient for CFD lambda_TVA (unitless)")             # Le 2015 Table 2 row 'lambda_TVA' (0.47, RSE 24%)
    llambda_VR  <- log(7.3);          label("Vitreous-to-retina partition coefficient for lampalizumab lambda_VR (unitless)")      # Le 2015 Table 2 row 'lambda_VR' (7.3, RSE 15%)
    lF1         <- log(0.84);         label("Intravitreal dose fraction reaching the vitreous F1 (unitless; the remaining 1-F1 leaks to serum by fast absorption)")  # Le 2015 Table 2 row 'F1' (0.84, RSE 1.2%)

    # ------------------------------------------------------------------
    # Systemic structural parameters. Le 2015 Table 2 ("Systemic parameters"
    # block).
    # ------------------------------------------------------------------
    lVser <- log(130);        label("Serum (central) volume V_ser (mL)")                                                             # Le 2015 Table 2 row 'V_ser' (130 mL, RSE 5.8%)
    lk    <- log(22);         label("Systemic first-order elimination rate constant k of unbound lampalizumab (1/day)")               # Le 2015 Table 2 row 'k' (22 /day, RSE 6.6%; paper text gives 21.3 /day for the derived 0.8-hour half-life)
    lkdeg <- log(96);         label("Systemic elimination rate constant of free CFD k_deg (1/day)")                                   # Le 2015 Table 2 row 'k_deg' (96 /day, RSE 10%)
    lkint <- log(4.2);        label("Systemic elimination rate constant of lampalizumab-CFD complex k_int (1/day)")                   # Le 2015 Table 2 row 'k_int' (4.2 /day, RSE 5.6%; paper text uses symbol k_inT interchangeably)
    lksyn <- log(2.6);        label("Zero-order synthesis rate of CFD in serum k_syn (nmol/mL/day)")                                   # Le 2015 Table 2 row 'k_syn' (2.6 nmol/mL/day, RSE 12%)
    lkSV  <- log(0.00032);    label("Serum-to-vitreous transfer rate constant of CFD k_SV (1/day)")                                    # Le 2015 Table 2 row 'k_sv' (0.00032 /day, RSE 22%; also written 3e-4 /day in the paper text)
    lk12  <- log(3.4);        label("Central-to-peripheral distribution rate of lampalizumab k_12 (1/day)")                            # Le 2015 Table 2 row 'k_12' (3.4 /day, RSE 12%)
    lk21  <- log(1.8);        label("Peripheral-to-central distribution rate of lampalizumab k_21 (1/day)")                            # Le 2015 Table 2 row 'k_21' (1.8 /day, RSE 8.2%)
    lkt12 <- log(24);         label("Symmetric distribution rate of CFD between serum and target-peripheral k_t12 = k_t21 (1/day)")    # Le 2015 Table 2 row 'k_t12 = k_t21' (24 /day, RSE 54%; constrained symmetric by the paper -- see Methods paragraph after Eq. 11)

    # ------------------------------------------------------------------
    # QSS binding constant. Fixed by the paper to the in-vitro K_D of
    # 11.7 pM (Loyet et al. 2014). Convert 11.7 pM to nmol/mL:
    #   11.7 pM = 11.7e-12 mol/L * 1e9 nmol/mol / 1000 mL/L = 1.17e-5 nmol/mL.
    # Le 2015 Discussion (page 293): "the quasi-steady state equilibrium
    # constant K_ss was fixed to the in vitro measured value of the
    # dissociation equilibrium constant (K_D) of 11.7 pM (Loyet et al.
    # 2014)". log() inside fixed() per parameter-names.md.
    # ------------------------------------------------------------------
    lKss <- fixed(log(1.17e-5));  label("Quasi-steady-state binding constant K_ss (nmol/mL; FIXED to in-vitro K_D 11.7 pM = 1.17e-5 nmol/mL from Loyet 2014)")  # Le 2015 text page 290, Loyet 2014 K_D

    # ------------------------------------------------------------------
    # Residual error. The paper (Methods, "Target-Mediated Drug
    # Disposition Model") states: "The residual error model for the
    # observations was best described by a proportional error model,
    # and the proportional residual errors were assumed to be independent
    # and normally distributed with zero means." No numeric variance is
    # reported for any of the seven observation streams. Per the skill
    # policy for unreported RUV with structural values present, encode
    # each proportional residual SD as fixed(0). Documented in vignette
    # Errata.
    # ------------------------------------------------------------------
    propSd_CVITR <- fixed(0); label("Proportional residual SD for total vitreous lampalizumab (fraction; FIXED to 0 -- not reported)")   # Le 2015 Methods (proportional error, no numeric value)
    propSd_CAQ   <- fixed(0); label("Proportional residual SD for total aqueous lampalizumab (fraction; FIXED to 0 -- not reported)")    # Le 2015 Methods
    propSd_CRET  <- fixed(0); label("Proportional residual SD for total retinal lampalizumab (fraction; FIXED to 0 -- not reported)")    # Le 2015 Methods
    propSd       <- fixed(0); label("Proportional residual SD for total serum lampalizumab (fraction; FIXED to 0 -- not reported)")      # Le 2015 Methods
    propSd_RVITR <- fixed(0); label("Proportional residual SD for total vitreous CFD (fraction; FIXED to 0 -- not reported)")            # Le 2015 Methods
    propSd_RAQ   <- fixed(0); label("Proportional residual SD for total aqueous CFD (fraction; FIXED to 0 -- not reported)")             # Le 2015 Methods
    propSd_RSER  <- fixed(0); label("Proportional residual SD for total serum CFD (fraction; FIXED to 0 -- not reported)")               # Le 2015 Methods
  })

  model({
    # ---- Individual (typical-value) parameters ---------------------------
    VVITR      <- exp(lVVITR)
    kout       <- exp(lkout)
    koutC      <- exp(lkoutC)
    koutT      <- exp(lkoutT)
    lambda_VA  <- exp(llambda_VA)
    lambda_TVA <- exp(llambda_TVA)
    lambda_VR  <- exp(llambda_VR)
    F1         <- exp(lF1)

    Vser  <- exp(lVser)
    k     <- exp(lk)
    kdeg  <- exp(lkdeg)
    kint  <- exp(lkint)
    ksyn  <- exp(lksyn)
    kSV   <- exp(lkSV)
    k12   <- exp(lk12)
    k21   <- exp(lk21)
    kt12  <- exp(lkt12)
    kt21  <- kt12                                             # Le 2015 Methods (page 290, paragraph after Eq. 11): "k12T was assumed to be the same as k21T"

    Kss   <- exp(lKss)

    # ---- Total-concentration observables (Le 2015 Eqs. 12-13) ------------
    CVITR <- depot   / VVITR                                 # Total vitreous lampalizumab concentration (nmol/mL) (Le 2015 Eq. 12)
    CSER  <- central / Vser                                  # Total serum lampalizumab concentration (nmol/mL) (Le 2015 Eq. 13)
    RVITR <- total_target                                    # Total vitreous CFD concentration (nmol/mL); state IS a concentration
    RSER  <- total_target_central                            # Total serum CFD concentration (nmol/mL); state IS a concentration
    RPER  <- total_target_peripheral1                        # Total CFD in target peripheral compartment (nmol/mL); Le 2015 Eq. 6

    # ---- QSS approximation for unbound drug concentrations (Le 2015 Eqs. 7-8) ----
    # Rearranged from the equilibrium C_free + drug-target complex = C_total.
    # C_free = 0.5 * [(C - R - Kss) + sqrt((C - R - Kss)^2 + 4*Kss*C)].
    Cvu <- 0.5 * ((CVITR - RVITR - Kss) + sqrt((CVITR - RVITR - Kss)^2 + 4 * Kss * CVITR))   # Le 2015 Eq. 7
    Csu <- 0.5 * ((CSER  - RSER  - Kss) + sqrt((CSER  - RSER  - Kss)^2 + 4 * Kss * CSER))    # Le 2015 Eq. 8

    # ---- Baseline (pre-dose) steady states for CFD compartments ----------
    # Serum CFD steady state (drug = 0):
    #   0 = ksyn - kdeg * R_ser  =>  R_ser_ss = ksyn / kdeg.
    # Vitreous CFD steady state (drug = 0, C_vu = 0):
    #   0 = kSV * R_ser * V_ser / V_VITR - koutT * R_vitr
    #   =>  R_vitr_ss = kSV * R_ser_ss * V_ser / (koutT * V_VITR).
    # Peripheral CFD steady state (symmetric kt12 = kt21, drug = 0):
    #   0 = kt12 * R_ser - kt21 * R_per  =>  R_per_ss = R_ser_ss.
    RSER_ss  <- ksyn / kdeg
    RVITR_ss <- kSV * RSER_ss * Vser / (koutT * VVITR)
    RPER_ss  <- RSER_ss
    total_target(0)             <- RVITR_ss
    total_target_central(0)     <- RSER_ss
    total_target_peripheral1(0) <- RPER_ss

    # ---- ODE system (Le 2015 Eqs. 1-6) -----------------------------------

    # Eq. 1: total vitreous lampalizumab amount (nmol/day).
    d/dt(depot) <-
      -kout  * Cvu * VVITR -
       koutC * (RVITR * Cvu * VVITR) / (Kss + Cvu) +
       kinC  * (RSER  * Csu * Vser)  / (Kss + Csu)

    # Eq. 2: total vitreous CFD concentration (nmol/mL/day). The k_outC and
    # k_outT complex-loss terms combine RVITR * K_ss / (K_ss + C_vu) and
    # RVITR * C_vu / (K_ss + C_vu); the V_VITR factors cancel out in the
    # koutC term as printed in the paper.
    d/dt(total_target) <-
       kSV   * (RSER  * Kss * Vser) / ((Kss + Cvu) * VVITR) -
       koutT * (RVITR * Kss)        / (Kss + Cvu) -
       koutC * (RVITR * Cvu)        / (Kss + Cvu) +
       kinC  * (RSER  * Csu * Vser) / ((Kss + Csu) * VVITR)

    # Eq. 3: total serum lampalizumab amount (nmol/day). The factor of 2 on
    # the vitreous-to-serum flow terms represents bilateral ITV dosing --
    # both eyes contribute to serum via the same rate constants.
    d/dt(central) <-
       2 * ( kout  * Cvu * VVITR +
             koutC * (RVITR * Cvu * VVITR) / (Kss + Cvu) ) -
       k    * Csu * Vser -
       kint * (RSER * Csu / (Kss + Csu)) * Vser -
       2 * kinC * (RSER * Csu * Vser) / (Kss + Csu) -
       k12  * central + k21 * peripheral1

    # Eq. 4: total serum CFD concentration (nmol/mL/day). Same bilateral
    # factor of 2 on the vitreous-derived source terms.
    d/dt(total_target_central) <-
       ksyn -
       kdeg * (RSER * Kss / (Kss + Csu)) -
       kint * (RSER * Csu / (Kss + Csu)) -
       2 * kSV * (RSER * Kss) / (Kss + Cvu) +
       2 * koutC * (RVITR * Cvu * VVITR) / ((Kss + Cvu) * Vser) +
       2 * koutT * (RVITR * Kss * VVITR) / ((Kss + Cvu) * Vser) -
       2 * kinC * (RSER * Csu / (Kss + Csu)) -
       kt12 * RSER + kt21 * RPER

    # Eq. 5: peripheral lampalizumab (nmol/day).
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Eq. 6: peripheral CFD (nmol/mL/day).
    d/dt(total_target_peripheral1) <- kt12 * RSER - kt21 * RPER

    # ---- Bioavailability of intravitreal dosing (Le 2015 Fig. 1) ---------
    # ITV dosing splits into two paths: F1 of the per-eye dose enters the
    # vitreous humor directly (encoded as bioavailability on `depot`), and
    # 2*(1-F1) of the per-eye dose enters serum via fast microcirculation
    # absorption (the factor 2 accounts for the bilateral dose). The serum-
    # leak fraction must be delivered by a second event-table row targeting
    # `central`; users administering IV doses target `central` directly with
    # F = 1 (see vignette for event-table helpers).
    f(depot) <- F1

    # ---- Aqueous / retinal derived observations (Le 2015 Eqs. 9-11) ------
    CAQ  <- CVITR / lambda_VA                                 # Aqueous total lampalizumab (Le 2015 Eq. 9): C_AQ_tot = C_VITR / lambda_VA
    RAQ  <- RVITR / lambda_TVA                                # Aqueous total CFD (Le 2015 Eq. 10): R_AQ = R_VITR / lambda_TVA
    CRET <- CVITR / lambda_VR                                 # Retinal total lampalizumab (Le 2015 Eq. 11): C_RET_tot = C_VITR / lambda_VR

    # ---- Observation-and-error model (proportional per output) -----------
    # The paper reports a single proportional error structure ("The
    # residual error model for the observations was best described by a
    # proportional error model"). One residual per observed stream; all
    # SDs FIXED to zero (values not reported -- vignette Errata).
    CVITR ~ prop(propSd_CVITR)   # Total lampalizumab in vitreous humor
    CAQ   ~ prop(propSd_CAQ)     # Total lampalizumab in aqueous humor
    CRET  ~ prop(propSd_CRET)    # Total lampalizumab in retina
    CSER  ~ prop(propSd)         # Total lampalizumab in serum (canonical bare-suffix propSd)
    RVITR ~ prop(propSd_RVITR)   # Total CFD in vitreous humor
    RAQ   ~ prop(propSd_RAQ)     # Total CFD in aqueous humor
    RSER  ~ prop(propSd_RSER)    # Total CFD in serum
  })
}
