LledoGarcia_2022_rozanolixizumab_translated <- function() {
  description <- "Cynomolgus-monkey-to-human translated two-compartment quasi-equilibrium TMDD PK model for the anti-FcRn monoclonal antibody rozanolixizumab with an indirect-response model on endogenous IgG; the allometrically scaled prediction used to design the first-in-human study, before any human data (Lledo-Garcia 2022)"
  reference <- "Lledo-Garcia R, Dixon K, Shock A, Oliver R. Pharmacokinetic-pharmacodynamic modelling of the anti-FcRn monoclonal antibody rozanolixizumab: Translation from preclinical stages to the clinic. CPT Pharmacometrics Syst Pharmacol. 2022;11(1):116-128. doi:10.1002/psp4.12739"
  vignette <- "LledoGarcia_2022_rozanolixizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    central     = list(analyte = "rozanolixizumab", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "rozanolixizumab", units = "mg", specimen = "tissue", verified = TRUE),
    total_igg   = list(analyte = "endogenous immunoglobulin G", units = "g/L", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species       = "human (predicted; no human data used)",
    n_subjects    = 0L,
    n_studies     = 0L,
    disease_state = "healthy adults (simulated 75 kg reference subject)",
    dose_range    = "1, 4 and 7 mg/kg single intravenous infusion simulated",
    notes         = paste(
      "This is a forward prediction, not a fit. The typical values are the cynomolgus monkey estimates of",
      "LledoGarcia_2022_rozanolixizumab_cyno scaled to a 75 kg human (Supplementary Text S4:",
      "Y = Ycyno*(BWhuman/BWcyno)^b, with b = 0.75 for clearance-related parameters and b = 1 for volumes),",
      "combined with in vitro MDCK-cell binding data for KD and EC50 and with literature values for the IgG",
      "turnover parameters (Table 1). The implied cynomolgus body weight is about 3.17 kg (V 0.108 -> 2.55 L and",
      "V2 0.0162 -> 0.383 L both give a 23.6-fold ratio at exponent 1). Interindividual variability was assumed",
      "from in-house monoclonal-antibody data and from the literature rather than estimated (Table 1,",
      "'IIV human (CV %)' column). n_subjects is 0 because no subjects were fitted; the model was used to run the",
      "stochastic Monte Carlo simulations behind Figures 2 and 3."
    )
  )

  ini({
    # ---- PK: two-compartment TMDD with the quasi-equilibrium approximation ----
    # Typical values are the 'Human translated mean value' column of Table 1.
    # Every one is an assumption carried into the simulation rather than an
    # estimate from data, so all are fixed().
    lvc   <- fixed(log(2.55));    label("Central volume of distribution V, allometrically scaled to 75 kg (L)")   # Table 1: V = 2.55 L (0.108 L cyno, exponent 1)
    lq    <- fixed(log(0.237));   label("Intercompartmental flow of free drug Q, allometrically scaled to 75 kg (L/day)") # Table 1: Q = 0.237 L/day (0.0221 L/day cyno, exponent 0.75); Table 2 prints 0.27 -- see vignette Errata
    lvp   <- fixed(log(0.383));   label("Peripheral volume of distribution V2, allometrically scaled to 75 kg (L)") # Table 1: V2 = 0.383 L (0.0162 L cyno, exponent 1)
    lcl   <- fixed(log(0.7605));  label("Clearance of free rozanolixizumab CLA, allometrically scaled to 75 kg (L/day)") # Table 1: CLA = 0.7605 L/day (0.0709 L/day cyno, exponent 0.75)
    lkd   <- fixed(log(0.359));   label("Rozanolixizumab-FcRn equilibrium binding constant KD, translated from in vitro MDCK-cell affinity differences (nmol/L)") # Table 1: KD = 0.359 nM (0.827 nM cyno)
    lrtot <- fixed(log(592.89));  label("Total FcRn concentration [FcRn], assumed the same in monkey and human (nmol/L)") # Table 1: [FcRn] = 592.89 nmol/L, same value assumed for both species
    lkdeg <- fixed(log(0.277));   label("First-order degradation rate constant of FcRn and of the drug:FcRn complex Kdeg (1/day)") # Table 1: Kdeg = 0.277 /day (85%-weighted mean of 0.3389 and 0.1536); Table 2 prints 0.227 -- see vignette Errata

    # ---- PD: indirect response on total IgG ----
    lrbase <- fixed(log(11.65));  label("Baseline total IgG concentration IGSS, from literature (g/L)")           # Table 1: IgG base = 11.65 mg/mL, from literature (references 47, 48)
    lkout  <- fixed(log(0.031));  label("First-order catabolic rate constant of IgG under normal conditions Kout (1/day)") # Table 1: Kout = 0.031 /day, average of the allometrically scaled monkey value and the literature endogenous-IgG value
    lemax  <- fixed(log(7.58));   label("Maximum fractional stimulation of IgG catabolism Emax, assumed the same as in cynomolgus monkey (unitless)") # Table 1: Emax = 7.58, assumed identical across species
    lec50  <- fixed(log(0.4959)); label("Free rozanolixizumab concentration giving half-maximal stimulation of IgG catabolism EC50 (mg/L)") # Table 1: EC50 = 0.4959 mg/L (1.14 mg/L cyno, corrected as for KD)

    # ---- IIV (Table 1, 'IIV human (CV %)' column; omega^2 = log(CV^2 + 1)) ----
    # Assumed from in-house monoclonal-antibody data and from the literature,
    # not estimated, so each variance is fixed(). Table 1 reports no IIV for
    # Q, KD or [FcRn], which therefore carry none.
    #   CV 16.35% -> 0.0263812   (V)
    #   CV 28.75% -> 0.0794175   (V2)
    #   CV 35.15% -> 0.1164953   (CLA)
    #   CV 18.00% -> 0.0318862   (Kdeg)
    #   CV 21.00% -> 0.0431553   (IgG base)
    #   CV 22.00% -> 0.0472652   (Kout)
    #   CV 25.00% -> 0.0606246   (Emax)
    #   CV 58.00% -> 0.2899794   (EC50)
    etalvc    ~ fixed(0.0263812)   # Table 1: IIV V = 16.35% CV
    etalvp    ~ fixed(0.0794175)   # Table 1: IIV V2 = 28.75% CV
    etalcl    ~ fixed(0.1164953)   # Table 1: IIV CLA = 35.15% CV
    etalkdeg  ~ fixed(0.0318862)   # Table 1: IIV Kdeg = 18% CV
    etalrbase ~ fixed(0.0431553)   # Table 1: IIV IgG base = 21% CV
    etalkout  ~ fixed(0.0472652)   # Table 1: IIV Kout = 22% CV
    etalemax  ~ fixed(0.0606246)   # Table 1: IIV Emax = 25% CV
    etalec50  ~ fixed(0.2899794)   # Table 1: IIV EC50 = 58% CV

    # ---- Residual error ----
    # Not reported: Table 2 prints '-' in the 'Human translated parameters'
    # column for both Prop RE rows, because this model was used for simulation
    # rather than fitted. Encoded as zero rather than invented; see the
    # vignette Errata.
    propSd         <- fixed(0); label("Proportional residual error on free rozanolixizumab concentration (fraction; ZERO - not reported in source)")  # Table 2: '-' for the human translated column
    propSd_IgG_obs <- fixed(0); label("Proportional residual error on total IgG (fraction; ZERO - not reported in source)")                            # Table 2: '-' for the human translated column
  })

  model({
    # ---- Physical constants ----
    # Rozanolixizumab molecular weight. NOT reported anywhere in Lledo-Garcia
    # 2022, its supplement, or the cited characterisation paper (Smith 2018
    # mAbs 10:1111-1130); taken from the RYSTIGGO (rozanolixizumab-noli) US
    # prescribing information section 11 DESCRIPTION, 'approximate molecular
    # weight of 148 kDa' (DailyMed SPL setid 3c0eb8c2-c042-4954-b451-3baa77f5e6d1).
    # See the vignette Errata.
    mwt <- 148000                 # g/mol
    nm2mgl <- mwt / 1e6           # multiply nmol/L by this to get mg/L (= ug/mL): 0.148

    # ---- Individual parameters ----
    vc   <- exp(lvc + etalvc)                 # L
    q    <- exp(lq)                           # L/day
    vp   <- exp(lvp + etalvp)                 # L
    cl   <- exp(lcl + etalcl)                 # L/day  (CLA, clearance of FREE drug)
    kdeg <- exp(lkdeg + etalkdeg)             # 1/day
    kd   <- exp(lkd)   * nm2mgl               # mg/L
    rtot <- exp(lrtot) * nm2mgl               # mg/L

    # Supplementary Text S1: CLB = Kdeg*V clears both free FcRn and the
    # drug:FcRn complex, so the total FcRn pool is constant and needs no ODE.
    clb <- kdeg * vc                          # L/day

    rbase <- exp(lrbase + etalrbase)          # g/L
    kout  <- exp(lkout + etalkout)            # 1/day
    emax  <- exp(lemax + etalemax)            # unitless
    ec50  <- exp(lec50 + etalec50)            # mg/L

    # ---- Quasi-equilibrium binding (Supplementary Text S1, third equation) ----
    # Amount-scale form in the paper, divided through by V to give the
    # algebraically identical concentration-scale form; the free-drug root is
    # written in the numerically stable arrangement.
    ctot  <- central / vc                     # mg/L, total (free + complexed) drug
    disc  <- ctot - rtot - kd
    cfree <- 0.5 * (disc + sqrt(disc * disc + 4 * kd * ctot))   # mg/L, FAC(t)
    cplx  <- ctot - cfree                     # mg/L

    # ---- PK ODEs (Supplementary Text S1, first two equations) ----
    d/dt(central)     <- -cl * cfree - clb * cplx - q * cfree + (q / vp) * peripheral1
    d/dt(peripheral1) <-  q * cfree - (q / vp) * peripheral1

    # ---- PD ODE (Supplementary Text S2) ----
    eff <- emax * cfree / (cfree + ec50)
    total_igg(0) <- rbase
    d/dt(total_igg) <- rbase * kout - kout * total_igg * (1 + eff)

    # ---- Derived reporting quantities (not observations) ----
    rb <- exp(lrtot) * kdeg * vc              # nmol/day; Table 1 footnote: RB = [FcRn]*Kdeg*V
    Ctot <- ctot
    pctChangeIgG <- 100 * (total_igg - rbase) / rbase

    # ---- Observations ----
    # Free rozanolixizumab concentration; Figure 3a y axis is
    # 'Free drug concentration (ug/mL)'.
    Cc      <- cfree                          # ug/mL (= mg/L)
    IgG_obs <- total_igg                      # g/L

    Cc      ~ prop(propSd)
    IgG_obs ~ prop(propSd_IgG_obs)
  })
}
