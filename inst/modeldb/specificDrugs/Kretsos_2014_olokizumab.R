Kretsos_2014_olokizumab <- function() {
  description <- "Two-compartment population PK with linear elimination and SC first-order absorption (depot, central, peripheral1) plus effect-compartment fractional sigmoid Imax PD model for C-reactive protein (CRP) suppression in mild-to-moderate rheumatoid arthritis patients receiving single-dose IV or SC olokizumab (anti-IL-6 monoclonal antibody, IgG4, CDP6038). Final-analysis estimates from Kretsos et al. 2014 Table 1 (Final column), pooling first-in-human (healthy volunteers, Hickling 2011) and first-in-patient (Cohorts 1+2, n=27 active-treatment subjects) data. The PK observation model adds a per-subject endogenous anti-IL-6 baseline ('endo') as an additive offset on the observed OKZ concentration. Body weight was reported as a significant covariate on CL and central volume (paper Discussion) but its functional form / exponents were not reported in main text or supplement; the body-weight covariate effect is omitted here -- see vignette Assumptions."
  reference <- "Kretsos K, Jullion A, Zamacona M, Harari O, Shaw S, Boulanger B, Oliver R. Model-Based Optimal Design and Execution of the First-Inpatient Trial of the Anti-IL-6, Olokizumab. CPT Pharmacometrics Syst Pharmacol. 2014;3:e119. doi:10.1038/psp.2014.17"
  vignette <- "Kretsos_2014_olokizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")
  # Notes on units: Table 1 reports CL in L/day, V1/V2 in L, ka in 1/day; doses in
  # the trial were mg/kg but the Final NONMEM fit uses absolute mg doses (Interim 1
  # was on a per-kg basis; Interim 2 and Final are absolute). Plasma OKZ and CRP
  # are both in ug/mL (1 mg/L = 1 ug/mL).

  covariateData <- list(
    # The Final model in Kretsos 2014 Table 1 has no estimated covariate effects:
    # body weight was reported as a significant covariate on CL and V1 (paper
    # Discussion p.5) but the functional form and exponents are not given in the
    # main text or in the supplement (which contains only the PD NONMEM control
    # stream, not the PK control stream). The model therefore exposes no
    # covariate columns. Downstream users who wish to apply weight scaling should
    # add it externally using their own convention (e.g., a priori allometric
    # (WT/70)^0.75 on CL/Q and (WT/70)^1 on V1/V2) and document that choice.
  )

  population <- list(
    species        = "human",
    n_subjects     = 27L,                                          # Kretsos 2014 Results: 27 patients on active treatment across Cohort 1 (2 randomization blocks, 18 active) and Cohort 2 (1 randomization block, 9 active); placebo subjects (9 across both cohorts) contribute to the PD fit. First-in-human PK data (Hickling et al. 2011, healthy male volunteers) were pooled with the first-in-patient PK for the Final model but FIH subject count is not stated in this paper.
    n_studies      = 2L,                                           # First-in-human (Hickling 2011, NCT not stated) + first-in-patient (NCT01009242)
    age_range      = "adult (specific range not reported)",       # Kretsos 2014: study population is adults with mild-to-moderate RA; demographic table not in main text or accessible supplement
    weight_range   = "not reported",                              # Kretsos 2014 main text and supplement: no demographic table accessible
    sex_female_pct = NA_real_,                                     # Not reported
    race_ethnicity = NULL,                                         # Not reported in main text or accessible supplement
    disease_state  = "Mild-to-moderate rheumatoid arthritis. Baseline CRP eligibility windowed during the trial (median 3.37 ug/mL, range 0.6-27.2 ug/mL per paper Discussion); used as randomization-stratification variable. First-in-human cohort that contributed PK data only was healthy male volunteers.",
    dose_range     = "Single-dose IV (0.1, 1, 3 mg/kg) or SC (1, 3 mg/kg); placebo IV or SC. First-in-human study contributed additional single-dose IV / SC PK data; doses not enumerated in Kretsos 2014.",
    regions        = "Clinical sites in Europe (UK, Germany) and USA (Kretsos 2014 Acknowledgments; investigators in TX, AZ, KS).",
    notes          = "Trial NCT01009242. Sponsor UCB Pharma. Cohort 1 (n=27 randomized = 18 active + 9 placebo across two i.v. doses + one s.c. dose plus three additional subjects); Cohort 2 (n=12 = 9 active 3 mg/kg s.c. + 3 placebo s.c.). 18 Cohort-1 patients with at least 4 weeks postdose data fed the first interim analysis. All available Cohort-1 and Cohort-2 data fed the second interim and final analyses. Drug development codename: CDP6038."
  )

  ini({
    # === Structural PK parameters (Kretsos 2014 Table 1, Final column) ===
    lcl      <- log(0.168);    label("Linear clearance CL (L/day)")                                  # Kretsos 2014 Table 1 Final: CL = 0.168 L/day (%RSE 3.85)
    lvc      <- log(4.2);      label("Central volume of distribution V1 (L)")                        # Kretsos 2014 Table 1 Final: V1 = 4.2 L (%RSE 2.69)
    lq       <- log(0.356);    label("Inter-compartmental clearance Q (L/day)")                      # Kretsos 2014 Table 1 Final: Q  = 0.356 L/day (%RSE 16.9); Table 1 reports no IIV on Q
    lvp      <- log(2.09);     label("Peripheral volume of distribution V2 (L)")                     # Kretsos 2014 Table 1 Final: V2 = 2.09 L (%RSE 8.37)
    lka      <- log(0.162);    label("SC first-order absorption rate ka (1/day)")                    # Kretsos 2014 Table 1 Final: ka = 0.162 day^-1 (%RSE 16.5)
    lfdepot  <- log(0.756);    label("SC bioavailability F (fraction)")                              # Kretsos 2014 Table 1 Final: F  = 75.6% (%RSE 4.81); held at 0.756 in the supplement PD NONMEM stream (F1=0.756) but estimated in the upstream PK fit
    lendo    <- log(0.0741);   label("Endogenous anti-IL-6 baseline offset on observed OKZ (ug/mL)") # Kretsos 2014 Table 1 Final: endogenous anti-IL-6 = 0.0741 ug/mL (%RSE 6.94); paper Methods (p.6): "The effect of endogenous anti-IL-6 levels was included as an additive term in the error model"

    # === Structural PD parameters (Kretsos 2014 Table 1, Final column) ===
    le0      <- log(3.32);     label("Baseline CRP E0 (ug/mL)")                                      # Kretsos 2014 Table 1 Final: baseline CRP = 3.32 ug/mL (%RSE 15.3)
    emax     <- fixed(0.99);   label("Maximum fractional CRP suppression Imax (unitless)")           # Kretsos 2014 Table 1 (all columns): Emax = 0.99 (fixed); supplement PD NONMEM stream hard-codes 0.99 inline (EFF = BASE*(1 - 0.99*CE^GAM/(EC50^GAM + CE^GAM)))
    lec50    <- log(0.415);    label("Effect-site OKZ for half-maximal CRP suppression EC50 (ug/mL)") # Kretsos 2014 Table 1 Final: EC50 = 0.415 ug/mL (%RSE 23.4)
    lhill   <- log(0.869);    label("Sigmoidicity factor hill (unitless)")                         # Kretsos 2014 Table 1 Final: hill = 0.869 (%RSE 17.4)
    lke0     <- log(0.0858);   label("Effect-compartment equilibration rate ke0 (1/day)")            # Kretsos 2014 Table 1 Final: ke0 = 0.0858 day^-1 (%RSE 27)

    # === IIV (Kretsos 2014 Table 1, Final column; diagonal $OMEGA per supplement PD NONMEM stream) ===
    # Paper footnote (Table 1): "The %CV for both intersubject and residual variability is an
    # approximation taken as the square root of the variance x 100." This is the simple
    # variance-on-log-scale approximation, so omega^2 = (CV/100)^2 (NOT the log-normal form
    # omega^2 = log(1 + CV^2)). The Final residual / IIV are encoded on that scale to match
    # the paper's reporting convention.
    etalcl    ~ 0.0590    # CV 24.3% on log CL    -> (0.243)^2 = 0.0590; Kretsos 2014 Table 1 Final IIV on CL
    etalvc    ~ 0.0161    # CV 12.7% on log V1    -> (0.127)^2 = 0.0161; Kretsos 2014 Table 1 Final IIV on V1
    etalvp    ~ 0.1832    # CV 42.8% on log V2    -> (0.428)^2 = 0.1832; Kretsos 2014 Table 1 Final IIV on V2
    etalka    ~ 0.6823    # CV 82.6% on log ka    -> (0.826)^2 = 0.6823; Kretsos 2014 Table 1 Final IIV on ka
    etalendo  ~ 0.1069    # CV 32.7% on log endo  -> (0.327)^2 = 0.1069; Kretsos 2014 Table 1 Final IIV on endogenous anti-IL-6
    etale0    ~ 1.0506    # CV 102.5% on log E0   -> (1.025)^2 = 1.0506; Kretsos 2014 Table 1 Final IIV on E0
    etalec50  ~ 0.6100    # CV 78.1% on log EC50  -> (0.781)^2 = 0.6100; Kretsos 2014 Table 1 Final IIV on EC50
    etalhill ~ 1.0302    # CV 101.5% on log hill-> (1.015)^2 = 1.0302; Kretsos 2014 Table 1 Final IIV on hill
    etalke0   ~ 0.7868    # CV 88.7% on log ke0   -> (0.887)^2 = 0.7868; Kretsos 2014 Table 1 Final IIV on ke0
    # NOTE: Kretsos 2014 Table 1 reports no IIV on Q (inter-compartmental clearance) or
    # F (SC bioavailability); both are population-only here.

    # === Residual error (Kretsos 2014 Table 1, Final column; proportional on each output) ===
    propSd     <- 0.240    # PK proportional residual SD on Cc (fraction); Kretsos 2014 Table 1 Final: 24.0% (%RSE 14.3). Supplement PD NONMEM stream shows Y=EFF+EFF*EPS(1) -> proportional. PK error model in main text Methods (Interim 1 used combined prop+add; Interim 2 / Final simplified to proportional only).
    propSd_crp <- 0.449    # PD proportional residual SD on CRP (fraction); Kretsos 2014 Table 1 Final: 44.9% (%RSE 22.8). Supplement PD NONMEM stream: Y=EFF+EFF*EPS(1).
  })

  model({
    # === Individual structural PK parameters ===
    cl      <- exp(lcl  + etalcl)
    vc      <- exp(lvc  + etalvc)
    q       <- exp(lq)                       # no IIV on Q per Table 1
    vp      <- exp(lvp  + etalvp)
    ka      <- exp(lka  + etalka)
    fdepot  <- exp(lfdepot)                  # no IIV on F per Table 1
    endo    <- exp(lendo + etalendo)         # per-subject endogenous anti-IL-6 baseline (additive offset on observed OKZ)

    # === Individual structural PD parameters ===
    e0      <- exp(le0    + etale0)
    ec50    <- exp(lec50  + etalec50)
    hill   <- exp(lhill + etalhill)
    ke0     <- exp(lke0   + etalke0)

    # === Micro-rate constants ===
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # === ODE system: 2-cmt PK with SC depot + effect compartment for CRP suppression ===
    # PK structure follows the supplement PD NONMEM stream:
    #   DADT(1) = -KA*A(1)
    #   DADT(2) =  KA*A(1) - CL/V2*A(2) - Q/V2*A(2) + Q/V3*A(3)     (NONMEM V2 = central -> nlmixr2 vc)
    #   DADT(3) =  Q/V2*A(2) - Q/V3*A(3)                            (NONMEM V3 = peripheral -> nlmixr2 vp)
    #   DADT(4) =  KE0 * (CP - A(4))    with CP = A(2)/V2
    # nlmixr2 form below.
    d/dt(depot)       <- -ka  * depot
    d/dt(central)     <-  ka  * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # SC bioavailability on depot (IV doses go directly to `central` with no `f()` -> F = 1)
    f(depot) <- fdepot

    # Model-predicted plasma OKZ concentration (drives the effect compartment;
    # the endogenous anti-IL-6 offset enters at the PK observation step, not here).
    ccent <- central / vc

    # Effect compartment driven by plasma OKZ (Sheiner effect-compartment formalism;
    # effect-compartment "concentration" carries the same unit as central concentration).
    d/dt(effect) <- ke0 * (ccent - effect)

    # === Observations and error models ===
    # PK: observed OKZ = model concentration + per-subject endogenous anti-IL-6 offset (additive in the error model)
    Cc  <- ccent + endo

    # PD: CRP follows a sigmoid Imax suppression of the baseline driven by the effect compartment
    crp <- e0 * (1 - emax * effect^hill / (ec50^hill + effect^hill))

    Cc  ~ prop(propSd)
    crp ~ prop(propSd_crp)
  })
}
