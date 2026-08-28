LledoGarcia_2022_rozanolixizumab_cyno <- function() {
  description <- "Preclinical (cynomolgus monkey). Two-compartment quasi-equilibrium TMDD PK model for the anti-FcRn monoclonal antibody rozanolixizumab with an indirect-response model in which free drug stimulates endogenous IgG catabolism; fitted to single-dose intravenous cynomolgus monkey data (Lledo-Garcia 2022)"
  reference <- "Lledo-Garcia R, Dixon K, Shock A, Oliver R. Pharmacokinetic-pharmacodynamic modelling of the anti-FcRn monoclonal antibody rozanolixizumab: Translation from preclinical stages to the clinic. CPT Pharmacometrics Syst Pharmacol. 2022;11(1):116-128. doi:10.1002/psp4.12739"
  vignette <- "LledoGarcia_2022_rozanolixizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    central     = list(analyte = "rozanolixizumab", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "rozanolixizumab", units = "mg", specimen = "tissue", verified = TRUE),
    total_igg   = list(analyte = "endogenous immunoglobulin G", units = "g/L", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species       = "cynomolgus monkey (Macaca fascicularis)",
    n_subjects    = 12L,
    n_studies     = 1L,
    sex_female_pct = 0,
    disease_state = "healthy young adult animals",
    dose_range    = "5, 10 and 30 mg/kg single intravenous infusion",
    weight_range  = "at least 2.3 kg",
    age_range     = "approximately 135-185 weeks at first dose",
    notes         = paste(
      "Only single-dose intravenous data were used for model development; multiple-dose data collected before",
      "antidrug antibodies were detected were used for external validation (Lledo-Garcia 2022 Methods, and",
      "Figure S2 of the supplement). The underlying animal study is Smith B, Kiessling A, Lledo-Garcia R, et al.",
      "mAbs 2018;10(7):1111-1130 (reference 9 of the source paper), which reports four healthy young adult male",
      "animals per dosing group, so the three single-dose intravenous groups give 12 animals. Antidrug antibodies",
      "were detectable in most animals during the dosing period but had no observed effect on PK or PD.",
      "The body weight implied by the paper's own allometric translation to a 75 kg human is about 3.17 kg."
    )
  )

  ini({
    # ---- PK: two-compartment TMDD with the quasi-equilibrium approximation ----
    # Typical values are the 'Cynomolgus monkey estimated typical mean value'
    # column of Table 1.
    lvc   <- log(0.108);   label("Central volume of distribution V (L)")                               # Table 1: V = 0.108 L
    lq    <- log(0.0221);  label("Intercompartmental flow of free drug Q (L/day)")                     # Table 1: Q = 0.0221 L/day
    lvp   <- log(0.0162);  label("Peripheral volume of distribution V2 (L)")                           # Table 1: V2 = 0.0162 L
    lcl   <- log(0.0709);  label("Clearance of free rozanolixizumab CLA (L/day)")                      # Table 1: CLA = 0.0709 L/day
    lkd   <- log(0.827);   label("Rozanolixizumab-FcRn equilibrium binding constant KD (nmol/L)")      # Table 1: KD = 0.827 nM (MDCK cells in vitro ~1 nM)
    lrtot <- log(592.89);  label("Total FcRn concentration [FcRn] (nmol/L)")                           # Table 1: [FcRn] = 592.89 nmol/L (0.6 uM)
    lkdeg <- log(0.3389);  label("First-order degradation rate constant of FcRn and of the drug:FcRn complex Kdeg (1/day)") # Table 1: Kdeg = 0.3389 /day (CLB = 0.0366 L/day = Kdeg*V)

    # ---- PD: indirect response on total IgG ----
    lrbase <- log(10.5);   label("Baseline total IgG concentration IGSS (g/L)")                        # Table 1: IgG base = 10.5 mg/mL
    lkout  <- log(0.0431); label("First-order catabolic rate constant of IgG under normal conditions Kout (1/day)") # Table 1: Kout = 0.0431 /day
    lemax  <- log(7.58);   label("Maximum fractional stimulation of IgG catabolism Emax (unitless)")   # Table 1: Emax = 7.58
    lec50  <- log(1.14);   label("Free rozanolixizumab concentration giving half-maximal stimulation of IgG catabolism EC50 (mg/L)") # Table 1: EC50 = 1.14 mg/L

    # ---- Residual error ----
    # Not reported. Table 1 gives only typical values for the cynomolgus model;
    # its 'IIV human (CV %)' and 'Uncertainty (RSE %)' columns describe the
    # HUMAN translation, not the monkey fit. Supplementary Text S3 states that a
    # proportional residual error model was used for the final model but gives no
    # magnitude for the monkey fit. Encoded as zero rather than invented; see the
    # vignette Errata.
    propSd         <- fixed(0); label("Proportional residual error on free rozanolixizumab concentration (fraction; ZERO - not reported in source)") # magnitude not reported for the cynomolgus fit
    propSd_IgG_obs <- fixed(0); label("Proportional residual error on total IgG (fraction; ZERO - not reported in source)")                           # magnitude not reported for the cynomolgus fit
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
    # No interindividual variability: none is reported for the cynomolgus fit.
    vc   <- exp(lvc)                          # L
    q    <- exp(lq)                           # L/day
    vp   <- exp(lvp)                          # L
    cl   <- exp(lcl)                          # L/day  (CLA, clearance of FREE drug)
    kdeg <- exp(lkdeg)                        # 1/day
    kd   <- exp(lkd)   * nm2mgl               # mg/L
    rtot <- exp(lrtot) * nm2mgl               # mg/L

    # Supplementary Text S1: CLB = Kdeg*V clears both free FcRn and the
    # drug:FcRn complex, so the total FcRn pool is constant and needs no ODE.
    clb <- kdeg * vc                          # L/day

    rbase <- exp(lrbase)                      # g/L
    kout  <- exp(lkout)                       # 1/day
    emax  <- exp(lemax)                       # unitless
    ec50  <- exp(lec50)                       # mg/L

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
    # Free rozanolixizumab concentration; the cynomolgus assay captured drug on
    # immobilised FcRn (Smith 2018 Methods), i.e. it reports FcRn-binding-competent
    # drug, and Figure 3a of the source paper labels its human counterpart
    # 'Free drug concentration (ug/mL)'.
    Cc      <- cfree                          # ug/mL (= mg/L)
    IgG_obs <- total_igg                      # g/L

    Cc      ~ prop(propSd)
    IgG_obs ~ prop(propSd_IgG_obs)
  })
}
