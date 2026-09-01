LledoGarcia_2022_rozanolixizumab <- function() {
  description <- "Two-compartment quasi-equilibrium TMDD PK model for the anti-FcRn monoclonal antibody rozanolixizumab with an indirect-response model in which free drug stimulates endogenous IgG catabolism; final model estimated from first-in-human intravenous data (Lledo-Garcia 2022)"
  reference <- "Lledo-Garcia R, Dixon K, Shock A, Oliver R. Pharmacokinetic-pharmacodynamic modelling of the anti-FcRn monoclonal antibody rozanolixizumab: Translation from preclinical stages to the clinic. CPT Pharmacometrics Syst Pharmacol. 2022;11(1):116-128. doi:10.1002/psp4.12739"
  vignette <- "LledoGarcia_2022_rozanolixizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    central     = list(analyte = "rozanolixizumab", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "rozanolixizumab", units = "mg", specimen = "tissue", verified = TRUE),
    total_igg   = list(analyte = "endogenous immunoglobulin G", units = "g/L", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species       = "human",
    n_subjects    = 18L,
    n_studies     = 1L,
    disease_state = "healthy adult volunteers",
    dose_range    = "1, 4 and 7 mg/kg single intravenous infusion (6 subjects per dose cohort)",
    notes         = paste(
      "First-in-human study NCT02220153 (Kiessling et al. 2017 Sci Transl Med 9:eaan1208; reference 11 of the",
      "source paper). That study randomised 49 subjects to rozanolixizumab (n = 36) or placebo (n = 13) across",
      "six cohorts; the first three cohorts were intravenous at 1, 4 and 7 mg/kg (n = 6 each) and the last three",
      "were subcutaneous. Lledo-Garcia 2022 Methods states that only the intravenous cohorts were analysed for",
      "this model, so n_subjects counts the 18 rozanolixizumab-treated intravenous subjects. The 7 intravenous",
      "placebo subjects contributed the IgG time course shown in the 'Placebo SD' panel of Figure 3. Dose",
      "escalation was halted at 7 mg/kg on tolerability grounds, so the model is informed by a narrower dose",
      "range than originally planned (Lledo-Garcia 2022 Methods, 'Update of the monkey-to-human translated",
      "PK/PD model with FIH data'). Baseline demographics are not tabulated in Lledo-Garcia 2022."
    )
  )

  ini({
    # ---- PK: two-compartment TMDD with the quasi-equilibrium approximation ----
    # All values are the FIH-estimated column of Table 2. Q, V2, KD and [FcRn]
    # were estimated with prior information (NWPRIOR) rather than held fixed --
    # each carries an RSE in Table 2 -- so none of them takes fixed().
    lvc   <- log(2.7);    label("Central volume of distribution V (L)")                                # Table 2 FIH: V = 2.7 L (RSE 6.3%)
    lq    <- log(0.271);  label("Intercompartmental flow of free drug Q (L/day)")                      # Table 2 FIH: Q = 0.271 L/day (RSE 16.6%); estimated with prior
    lvp   <- log(0.36);   label("Peripheral volume of distribution V2 (L)")                            # Table 2 FIH: V2 = 0.36 L (RSE 20.8%); estimated with prior
    lcl   <- log(0.968);  label("Clearance of free rozanolixizumab CLA (L/day)")                       # Table 2 FIH: CLA = 0.968 L/day (RSE 6.8%)
    lkd   <- log(1);      label("Rozanolixizumab-FcRn equilibrium binding constant KD (nmol/L)")       # Table 2 FIH: KD = 1 nM (RSE 14.6%); estimated with prior
    lrtot <- log(147);    label("Total FcRn concentration [FcRn] (nmol/L)")                            # Table 2 FIH: [FcRn] = 147 nM (RSE 14.4%); estimated with prior
    lkdeg <- log(0.88);   label("First-order degradation rate constant of FcRn and of the drug:FcRn complex Kdeg (1/day)") # Table 2 FIH: Kdeg = 0.88 /day (RSE 15.5%)

    # ---- PD: indirect response on total IgG ----
    lrbase <- log(9.88);   label("Baseline total IgG concentration IGSS (g/L)")                        # Table 2 FIH: IgG base = 9.88 mg/mL (RSE 1.8%)
    lkout  <- log(0.0364); label("First-order catabolic rate constant of IgG under normal conditions Kout (1/day)") # Table 2 FIH: Kout = 0.0364 /day (RSE 4.4%)
    lemax  <- log(4.24);   label("Maximum fractional stimulation of IgG catabolism Emax (unitless)")   # Table 2 FIH: Emax = 4.24 (RSE 5.5%)
    lec50  <- log(0.154);  label("Free rozanolixizumab concentration giving half-maximal stimulation of IgG catabolism EC50 (mg/L)") # Table 2 FIH: EC50 = 0.154 mg/L (RSE 15.3%)

    # ---- IIV (Table 2, FIH estimated column; omega^2 = log(CV^2 + 1)) ----
    #   CV 15.8% -> 0.0246575   (V)
    #   CV 30.3% -> 0.0878360   (Kdeg)
    #   CV 18.4% -> 0.0332955   (IgG base)
    etalvc    ~ 0.0246575   # Table 2 FIH: IIV V = 15.8% CV (RSE 21.2%)
    etalkdeg  ~ 0.0878360   # Table 2 FIH: IIV Kdeg = 30.3% CV (RSE 18.1%)
    etalrbase ~ 0.0332955   # Table 2 FIH: IIV IgG base = 18.4% CV (RSE 18.7%)

    # ---- Residual error (Table 2, FIH estimated column) ----
    propSd         <- 0.093; label("Proportional residual error on free rozanolixizumab concentration (fraction)") # Table 2 FIH: Prop RE = 9.3% CV (RSE 7.4%)
    propSd_IgG_obs <- 0.058; label("Proportional residual error on total IgG (fraction)")             # Table 2 FIH: Prop RE = 5.8% CV (RSE 2.6%)
  })

  model({
    # ---- Physical constants ----
    # Rozanolixizumab molecular weight. NOT reported anywhere in Lledo-Garcia
    # 2022, its supplement, or the cited characterisation paper (Smith 2018
    # mAbs 10:1111-1130); taken from the RYSTIGGO (rozanolixizumab-noli) US
    # prescribing information section 11 DESCRIPTION, 'approximate molecular
    # weight of 148 kDa' (DailyMed SPL setid 3c0eb8c2-c042-4954-b451-3baa77f5e6d1).
    # See the vignette Errata. The constant is load-bearing only through the
    # binding stoichiometry: the paper reports KD and [FcRn] in nmol/L but EC50
    # and the observed concentrations in mg/L, so a single conversion is needed
    # to put drug and target on one scale. Following the Said_2025_imatinib
    # precedent, a pure unit conversion is hardcoded here rather than declared
    # in ini().
    mwt <- 148000                 # g/mol
    nm2mgl <- mwt / 1e6           # multiply nmol/L by this to get mg/L (= ug/mL): 0.148

    # ---- Individual parameters ----
    vc   <- exp(lvc + etalvc)                 # L
    q    <- exp(lq)                           # L/day
    vp   <- exp(lvp)                          # L
    cl   <- exp(lcl)                          # L/day  (CLA, clearance of FREE drug)
    kdeg <- exp(lkdeg + etalkdeg)             # 1/day
    kd   <- exp(lkd)   * nm2mgl               # mg/L   (KD converted from nmol/L)
    rtot <- exp(lrtot) * nm2mgl               # mg/L   (total FcRn converted from nmol/L)

    # Supplementary Text S1: CLB = Kdeg*V is the clearance of free FcRn and of
    # the drug:FcRn complex. Because the two are equal the total FcRn pool is
    # constant in time, so the paper carries no differential equation for it.
    clb <- kdeg * vc                          # L/day

    rbase <- exp(lrbase + etalrbase)          # g/L
    kout  <- exp(lkout)                       # 1/day
    emax  <- exp(lemax)                       # unitless
    ec50  <- exp(lec50)                       # mg/L

    # ---- Quasi-equilibrium binding (Supplementary Text S1, third equation) ----
    # The paper writes the complex on the amount scale as
    #   CPLX = ((KD*V + TA + TB) - sqrt((KD*V + TA + TB)^2 - 4*TA*TB)) / 2
    # with TB = [FcRn]*V the (constant) total FcRn amount. Dividing through by V
    # gives the algebraically identical concentration-scale form used here. The
    # free-drug root is written directly, in the numerically stable arrangement
    # cfree = (disc + sqrt(disc^2 + 4*kd*ctot))/2, so that the deeply-bound
    # low-dose regime does not lose precision to cancellation.
    ctot  <- central / vc                     # mg/L, total (free + complexed) drug
    disc  <- ctot - rtot - kd
    cfree <- 0.5 * (disc + sqrt(disc * disc + 4 * kd * ctot))   # mg/L, FAC(t)
    cplx  <- ctot - cfree                     # mg/L, drug:FcRn complex

    # ---- PK ODEs (Supplementary Text S1, first two equations) ----
    # dTA/dt   = -FA*CLA/V - CPLX*CLB/V - Q*FA/V + Q*APERI/V2 + input
    # dAPERI/dt = Q*FA/V - Q*APERI/V2
    # Only free antibody distributes to the peripheral compartment.
    d/dt(central)     <- -cl * cfree - clb * cplx - q * cfree + (q / vp) * peripheral1
    d/dt(peripheral1) <-  q * cfree - (q / vp) * peripheral1

    # ---- PD ODE (Supplementary Text S2) ----
    # EFF      = EMAX*FAC(t) / (FAC(t) + EC50)
    # dIGG/dt  = IGSS*KOUT - KOUT*IGG(t)*(1 + EFF)
    eff <- emax * cfree / (cfree + ec50)
    total_igg(0) <- rbase
    d/dt(total_igg) <- rbase * kout - kout * total_igg * (1 + eff)

    # ---- Derived reporting quantities (not observations) ----
    # Target production rate, Table 2 footnote b: RB = [FcRn]*Kdeg*V (nmol/day).
    rb <- exp(lrtot) * kdeg * vc
    # Total drug concentration, for users whose assay measures total rather than free.
    Ctot <- ctot
    # Percentage change from baseline IgG -- the quantity the paper plots.
    pctChangeIgG <- 100 * (total_igg - rbase) / rbase

    # ---- Observations ----
    # The PK observation is the FREE rozanolixizumab concentration: Figure 3a
    # of the source paper labels its y axis 'Free drug concentration (ug/mL)',
    # and the cynomolgus assay captured drug on immobilised FcRn (Smith 2018),
    # i.e. it reports FcRn-binding-competent drug.
    Cc      <- cfree                          # ug/mL (= mg/L)
    IgG_obs <- total_igg                      # g/L

    Cc      ~ prop(propSd)
    IgG_obs ~ prop(propSd_IgG_obs)
  })
}
