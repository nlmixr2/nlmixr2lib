AitOudhia_2012_canakinumab <- function() {
  description <- "Integrated population PK/PD model of canakinumab (anti-IL-1beta IgG1/k mAb) in adults with rheumatoid arthritis (Ait-Oudhia 2012). Two-compartment popPK for total canakinumab is coupled to a quasi-equilibrium target-binding model with endogenous IL-1beta (zero-order production ksyn, linear clearance CLL). Predicted free IL-1beta drives downstream PD: (1) a three-compartment CRP transduction chain with a power-law stimulation (beta) on free-IL-1beta ratio and an empirical amplification (gamma) on the input to the third compartment, and (2) a single-compartment ACR latent variable (ACRL) regulated by a sigmoid Emax on the drop in free IL-1beta below baseline plus a first-order placebo build-up; the latent is mapped to ACR20/50/70 response probabilities via a logit transform with a between-subject random effect. Only body weight was a significant covariate (allometric on CL, CLL, CLDL, Vc, Vp with reference 70 kg)."
  reference <- "Ait-Oudhia S, Lowe PJ, Mager DE. Bridging Clinical Outcomes of Canakinumab Treatment in Patients With Rheumatoid Arthritis With a Population Model of IL-1beta Kinetics. CPT Pharmacometrics Syst Pharmacol. 2012;1:e5. doi:10.1038/psp.2012.6"
  vignette <- "AitOudhia_2012_canakinumab"
  units <- list(
    time          = "day",
    dosing        = "mg",
    concentration = "ug/mL (total canakinumab); pg/mL (total / free IL-1beta); mg/L (CRP)"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on CL (and CLDL), CLL with exponent 3/4 and on Vc, Vp with exponent 1, centred at 70 kg (Ait-Oudhia 2012 Methods, Data analysis paragraph, page 9: 'CL = thetaCL * exp(etaCL) * (BWT/70)^(3/4) and V = thetav * exp(etav) * (BWT/70)').",
      source_name        = "BWT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 472L,
    n_active       = 349L,
    n_placebo      = 123L,
    n_studies      = 4L,
    age_range      = "18-87 years",
    age_median     = "57 years",
    weight_range   = "40-111 kg",
    weight_median  = "74 kg",
    sex_female_pct = 80,
    race_ethnicity = NA,
    disease_state  = "Active rheumatoid arthritis",
    dose_range     = "0.1 mg/kg to 900 mg as a 2-hour IV infusion or SC injection Q2W or Q4W; alone or with methotrexate",
    regions        = NA,
    notes          = "Demographics from Ait-Oudhia 2012 Results, 'Data were obtained...' paragraph (page 2). Body-weight, age, gender, and methotrexate were tested as covariates; only body weight was retained (page 2-3). Per-study breakdown is described in Supplementary Table S1 (not on disk for this extraction). Other covariate-effect coefficients (age, gender, methotrexate) were tested and dropped during model building and are not encoded here."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters -- typical 70 kg RA patient. All values from
    # Ait-Oudhia 2012 Table 1 unless otherwise noted.
    # ------------------------------------------------------------------
    lka     <- log(0.266);  label("First-order SC absorption rate ka (1/day)")                                 # Ait-Oudhia 2012 Table 1
    lvc     <- log(3.71);   label("Central volume of distribution Vc (L) for 70 kg patient")                   # Ait-Oudhia 2012 Table 1
    lvp     <- log(2.24);   label("Peripheral volume of distribution Vp (L) for 70 kg patient")                # Ait-Oudhia 2012 Table 1
    lcl     <- log(0.104);  label("Total canakinumab clearance CL = CLDL (L/day) for 70 kg patient")           # Ait-Oudhia 2012 Table 1
    lcll    <- log(13.7);   label("Free IL-1beta clearance CLL (L/day)")                                       # Ait-Oudhia 2012 Table 1
    lq      <- log(0.165);  label("Intercompartmental clearance Q (L/day) for 70 kg patient")                  # Ait-Oudhia 2012 Table 1
    lkd     <- log(0.38);   label("Equilibrium dissociation constant Kd (nmol/L)")                             # Ait-Oudhia 2012 Table 1
    lfdepot <- log(0.667);  label("Subcutaneous bioavailability F (fraction)")                                 # Ait-Oudhia 2012 Table 1 (paper used logit transform; typical F = 0.667)
    lksyn   <- log(7.4);    label("IL-1beta zero-order production rate ksyn (ng/day)")                         # Ait-Oudhia 2012 Table 1 footnote (secondary parameter)

    # ------------------------------------------------------------------
    # Allometric exponents on body weight, fixed at the canonical values
    # (Ait-Oudhia 2012 Methods, Data analysis paragraph, page 9).
    # ------------------------------------------------------------------
    e_wt_cl_q   <- fixed(0.75); label("Fixed allometric exponent on CL, CLDL, CLL, Q (unitless)")              # Ait-Oudhia 2012 Methods page 9 (fixed at 3/4)
    e_wt_vc_vp  <- fixed(1);    label("Fixed allometric exponent on Vc, Vp (unitless)")                        # Ait-Oudhia 2012 Methods page 9 (fixed at 1)

    # ------------------------------------------------------------------
    # CRP transduction parameters (Ait-Oudhia 2012 Table 2).
    # ------------------------------------------------------------------
    lcrp0   <- log(8.44);   label("Baseline CRP concentration (mg/L)")                                         # Ait-Oudhia 2012 Table 2
    lkout   <- log(1.06);   label("First-order transit rate constant for CRP chain (1/day)")                   # Ait-Oudhia 2012 Table 2
    lbeta   <- log(0.25);   label("Power coefficient beta for the free-IL-1beta stimulatory function (unitless)")    # Ait-Oudhia 2012 Table 2
    lgamma  <- log(1.92);   label("Power coefficient gamma applied to the input of the third CRP transit (unitless)") # Ait-Oudhia 2012 Table 2

    # ------------------------------------------------------------------
    # ACRx parameters (Ait-Oudhia 2012 Table 2). Values for the ACR
    # 'no-BSV reported' parameters are entered as point estimates.
    # ------------------------------------------------------------------
    lemax     <- log(0.741);          label("Maximum drug-driven effect on ACR latent variable (unitless)")    # Ait-Oudhia 2012 Table 2
    lec50     <- log(0.204);          label("Free-IL-1beta concentration inducing 50% of Emax on ACR (pg/mL)") # Ait-Oudhia 2012 Table 2
    ltau      <- log(55.9);           label("Mean transit time tau for ACR latent variable (day)")             # Ait-Oudhia 2012 Table 2
    lkplb     <- log(0.0524);         label("First-order rate constant for placebo onset on ACR (1/day)")      # Ait-Oudhia 2012 Table 2 (0.524 x 10^-1)
    plbmax    <- 0.259;               label("Maximum placebo effect on ACR latent variable (unitless)")        # Ait-Oudhia 2012 Table 2 (no RSE reported)

    # Population-level logit-scale shift for the ACR response. Fixed at
    # zero; the corresponding random effect etaacrshift carries the
    # between-subject variability reported in Table 2 (54.3% CV).
    acrshift  <- fixed(0);            label("Population-mean logit shift on ACR (logit units)")                # Structural anchor (etaacrshift carries the BSV)

    # ------------------------------------------------------------------
    # IIV -- diagonal (no inter-parameter correlations reported).
    # Apparent CV values from the 'Variability (%RSE)' column of Table 1
    # / Table 2 are converted to log-scale variance via
    # omega^2 = log(1 + CV^2).
    # ------------------------------------------------------------------
    etalka     ~ 0.002497  # log(1 + 0.05^2);  Ait-Oudhia 2012 Table 1 ka  CV 5%
    etalvc     ~ 0.128356  # log(1 + 0.37^2);  Ait-Oudhia 2012 Table 1 Vc  CV 37%
    etalvp     ~ 0.147654  # log(1 + 0.399^2); Ait-Oudhia 2012 Table 1 Vp  CV 39.9%
    etalcl     ~ 0.013832  # log(1 + 0.118^2); Ait-Oudhia 2012 Table 1 CL  CV 11.8%
    etalcll    ~ 0.099887  # log(1 + 0.324^2); Ait-Oudhia 2012 Table 1 CLL CV 32.4%
    etalq      ~ 0.017298  # log(1 + 0.132^2); Ait-Oudhia 2012 Table 1 Q   CV 13.2%
    etalkd     ~ 0.109419  # log(1 + 0.34^2);  Ait-Oudhia 2012 Table 1 Kd  CV 34%
    etalfdepot ~ 0.001353  # log(1 + 0.0368^2); Ait-Oudhia 2012 Table 1 F  CV 3.68%
    etalcrp0   ~ 0.388700  # log(1 + 0.669^2); Ait-Oudhia 2012 Table 2 CRP0  CV 66.9%
    etalkout   ~ 0.099887  # log(1 + 0.324^2); Ait-Oudhia 2012 Table 2 kout  CV 32.4%
    etalbeta   ~ 0.435942  # log(1 + 0.753^2); Ait-Oudhia 2012 Table 2 beta  CV 75.3%
    etalgamma  ~ 0.330999  # log(1 + 0.636^2); Ait-Oudhia 2012 Table 2 gamma CV 63.6%
    etaacrshift ~ 0.258806  # log(1 + 0.543^2); Ait-Oudhia 2012 Table 2 ACR BSV CV 54.3% (added on the logit scale of the ACRx probability)

    # ------------------------------------------------------------------
    # Residual error -- combined additive + proportional for total
    # canakinumab and total IL-1beta (Ait-Oudhia 2012 Eq. 14); the
    # CRP residual is exponential (Eq. 15) and is implemented as
    # proportional (the small-sigma equivalent in nlmixr2 linear space).
    #
    # NOTE on IL-1beta additive units: Table 1 reports b2 = 0.317 nmol/L
    # for the additive component on total IL-1beta. Converting via
    # MW = 17 kDa gives 0.317 * 17000 = 5389 pg/mL, which is large
    # relative to the observed baseline (~0.5 pg/mL) and peak (~100 pg/mL)
    # IL-1beta concentrations shown in Figure 2. The value is recorded
    # here verbatim from the publication and may reflect either a
    # publication unit-labeling error or the residual variability
    # actually estimated on the model's internal molar scale; see the
    # vignette Assumptions and deviations section. Proportional error
    # (61.6% CV) dominates the predicted IL-1beta variability either way.
    # ------------------------------------------------------------------
    propSd            <- 0.111;     label("Proportional residual error on total canakinumab (fraction)")     # Ait-Oudhia 2012 Table 1 a1 = 11.1%
    addSd             <- 0.0326;    label("Additive residual error on total canakinumab (ug/mL)")            # Ait-Oudhia 2012 Table 1 b1 = 0.217 nmol/L -> 0.217 * 0.150 = 0.0326 ug/mL
    propSd_totalIL1b  <- 0.616;     label("Proportional residual error on total IL-1beta (fraction)")        # Ait-Oudhia 2012 Table 1 a2 = 61.6%
    addSd_totalIL1b   <- 5389;      label("Additive residual error on total IL-1beta (pg/mL); see Errata")   # Ait-Oudhia 2012 Table 1 b2 = 0.317 nmol/L -> 0.317 * 17000 = 5389 pg/mL (apparent unit anomaly noted)
    propSd_crp        <- 0.111;     label("Proportional residual error on CRP (fraction; from exponential)") # Ait-Oudhia 2012 Table 2 a = 11.1% (exponential error term)
  })

  model({
    # ------------------------------------------------------------------
    # Constants -- molecular weights and bookkeeping. The paper
    # explicitly cites canakinumab MW = 150 kDa (Methods page 7) and
    # IL-1beta MW = 17 kDa (Methods page 8). Dimensional convention:
    # 1 nmol/L * 1 kDa = 1 ng/mL.
    # ------------------------------------------------------------------
    MWX <- 150.0  # canakinumab molecular weight (kDa = ng/nmol)
    MWE <- 17.0   # IL-1beta molecular weight (kDa = ng/nmol)

    # ------------------------------------------------------------------
    # 1. Individual parameters with allometric scaling on body weight
    #    (Ait-Oudhia 2012 Data analysis paragraph, page 9).
    # ------------------------------------------------------------------
    ka     <- exp(lka  + etalka)
    cl     <- exp(lcl  + etalcl)  * (WT / 70)^e_wt_cl_q
    vc     <- exp(lvc  + etalvc)  * (WT / 70)^e_wt_vc_vp
    vp     <- exp(lvp  + etalvp)  * (WT / 70)^e_wt_vc_vp
    q      <- exp(lq   + etalq)   * (WT / 70)^e_wt_cl_q
    cll    <- exp(lcll + etalcll) * (WT / 70)^e_wt_cl_q
    kd     <- exp(lkd  + etalkd)
    fdepot <- exp(lfdepot + etalfdepot)
    cl_dl  <- cl                   # CLDL = CL (Ait-Oudhia 2012 Results, 'During the model building process...' page 2)

    # IL-1beta production in molar units: ksyn (ng/day) -> nmol/day
    # divides by MWE * 1000 = 17000 ng/nmol.
    ksyn   <- exp(lksyn) / (MWE * 1000)

    # CRP and ACR individual parameters
    crp0   <- exp(lcrp0  + etalcrp0)
    kout   <- exp(lkout  + etalkout)
    beta   <- exp(lbeta  + etalbeta)
    gamma  <- exp(lgamma + etalgamma)
    emax   <- exp(lemax)
    ec50   <- exp(lec50)
    tau    <- exp(ltau)
    kplb   <- exp(lkplb)

    # ------------------------------------------------------------------
    # 2. Quasi-equilibrium TMDD binding (Ait-Oudhia 2012 Eq. 5).
    #    Tc = central (nmol total drug), TL = central_il1b (nmol total
    #    IL-1beta), Kd in nmol/L, Vc in L. Free drug amount AC (nmol):
    #      AC = 0.5 * [(Tc - TL - Kd*Vc)
    #                  + sqrt((Tc - TL - Kd*Vc)^2 + 4 * Kd * Vc * Tc)]
    # ------------------------------------------------------------------
    diff_tdl <- central - central_il1b - kd * vc
    ac       <- 0.5 * (diff_tdl + sqrt(diff_tdl * diff_tdl + 4 * kd * vc * central))

    # Total / free IL-1beta concentrations in molar units (Ait-Oudhia
    # 2012 Eq. 5 continuation): drug-ligand complex CDL is the bound
    # fraction; free IL-1beta = total IL-1beta - CDL.
    ctl_nmol     <- central_il1b / vc
    free_drug    <- ac / vc                          # nmol/L free canakinumab
    cdl_nmol     <- ctl_nmol * free_drug / (kd + free_drug)
    free_il1b    <- ctl_nmol - cdl_nmol              # nmol/L free IL-1beta

    # Baseline free IL-1beta (steady-state with no drug): at t=0 there
    # is no complex, so free IL-1beta = total IL-1beta = ksyn / CLL.
    cfl_baseline <- ksyn / cll                       # nmol/L baseline free IL-1beta

    # Convert free IL-1beta concentrations to pg/mL for the PD driver
    # (EC50 is reported in pg/mL).
    free_il1b_pgmL    <- free_il1b   * MWE * 1000
    cfl_baseline_pgmL <- cfl_baseline * MWE * 1000

    # ------------------------------------------------------------------
    # 3. ODE system.
    #    States: depot (mg), central / peripheral1 (nmol total drug),
    #    central_il1b (nmol total IL-1beta), crp1..crp3 (mg/L; CRP
    #    transduction chain), acrl (unitless latent variable). The
    #    depot -> central transfer converts mg to nmol via 1000 / MWX.
    # ------------------------------------------------------------------
    d/dt(depot)        <- -ka * depot
    d/dt(central)      <-  ka * depot * (1000 / MWX) -
                           (cl + q) / vc * central + q / vp * peripheral1
    d/dt(peripheral1)  <-  q / vc * central - q / vp * peripheral1
    d/dt(central_il1b) <-  ksyn -
                           (cl_dl - cll) / vc * (central - ac) -
                           cll / vc * central_il1b

    # CRP zero-order production rate kin from the steady-state algebra:
    # at baseline, CRP1 = CRP2 = kin/kout and CRP3 = (kin/kout)^gamma =
    # CRP0; therefore kin = (CRP0 * kout^gamma)^(1/gamma). Stimulation
    # by free IL-1beta is the normalized power form (Methods page 8):
    # S(t) = (Cf,L / Cf,L(0))^beta.
    kin <- (crp0 * kout^gamma)^(1 / gamma)
    s_t <- (free_il1b / cfl_baseline)^beta

    d/dt(crp1) <- kin * s_t - kout * crp1
    d/dt(crp2) <- kout * (crp1 - crp2)
    d/dt(crp3) <- kout * (crp2^gamma - crp3)

    # ACR latent variable (Ait-Oudhia 2012 Eq. 9-11). PROD is the
    # production input shaped by the drug-driven drop in free IL-1beta
    # (Ef_L) and the first-order placebo build-up (Eplb).
    delta_il1b <- cfl_baseline_pgmL - free_il1b_pgmL
    ef_l       <- emax * delta_il1b / ((cfl_baseline_pgmL - ec50) + delta_il1b)
    eplb       <- plbmax * (1 - exp(-kplb * t))
    prod       <- 1 - ef_l - eplb
    d/dt(acrl) <- (prod - acrl) / tau

    # ------------------------------------------------------------------
    # 4. Initial conditions.
    #    Central / peripheral / depot start at zero (rxode2 default).
    #    central_il1b(0) = ksyn * vc / cll (Eq. 4 steady state, no drug).
    #    crp1(0) = crp2(0) = kin/kout; crp3(0) = (kin/kout)^gamma = CRP0.
    #    acrl(0) = 1 (Ait-Oudhia 2012 Eq. 9 text, page 8).
    # ------------------------------------------------------------------
    central_il1b(0) <- ksyn * vc / cll
    crp1(0)         <- kin / kout
    crp2(0)         <- kin / kout
    crp3(0)         <- crp0
    acrl(0)         <- 1

    # ------------------------------------------------------------------
    # 5. Bioavailability for the SC depot (Ait-Oudhia 2012 Methods,
    #    Eq. 1 initial condition Tsc(0) = F * Dose).
    # ------------------------------------------------------------------
    f(depot) <- fdepot

    # ------------------------------------------------------------------
    # 6. Observation outputs in assay units.
    #    Cc (ug/mL)        = (central / Vc) [nmol/L] * MWX / 1000
    #    totalIL1b (pg/mL) = (central_il1b / Vc) [nmol/L] * MWE * 1000
    #    freeIL1b  (pg/mL) = free_il1b_pgmL (derived above)
    #    crp (mg/L)        = crp3 directly
    # ------------------------------------------------------------------
    Cc        <- (central      / vc) * MWX / 1000
    totalIL1b <- ctl_nmol               * MWE * 1000
    freeIL1b  <- free_il1b_pgmL
    crp       <- crp3

    # Derived ACR response probabilities (logit transform, Eq. 12).
    # logit(1 - x/100): x = 20 -> log(4); x = 50 -> 0; x = 70 -> log(3/7).
    # logit(1 - ACRL) approaches -Inf as ACRL -> 1 (the t = 0 condition),
    # giving probability 0; clipped at small ACRL to avoid log(0) at the
    # other extreme.
    acrl_safe   <- 1 - acrl
    if (acrl_safe < 1e-12) acrl_safe <- 1e-12
    if (acrl_safe > 1 - 1e-12) acrl_safe <- 1 - 1e-12
    logit_acrl  <- log(acrl_safe / (1 - acrl_safe))

    logit_acr20 <- log(0.8 / 0.2) + logit_acrl + acrshift + etaacrshift
    logit_acr50 <- log(0.5 / 0.5) + logit_acrl + acrshift + etaacrshift
    logit_acr70 <- log(0.3 / 0.7) + logit_acrl + acrshift + etaacrshift

    prob_ACR20  <- 1 / (1 + exp(-logit_acr20))
    prob_ACR50  <- 1 / (1 + exp(-logit_acr50))
    prob_ACR70  <- 1 / (1 + exp(-logit_acr70))

    Cc        ~ prop(propSd)            + add(addSd)
    totalIL1b ~ prop(propSd_totalIL1b) + add(addSd_totalIL1b)
    crp       ~ prop(propSd_crp)
  })
}
