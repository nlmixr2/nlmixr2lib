Berges_2015_ozanezumab <- function() {
  description <- "Two-compartment IV population PK plus effect-compartment sigmoid Emax PKPD model for the proportion of skeletal-muscle membrane Nogo-A co-localized with ozanezumab in adults with amyotrophic lateral sclerosis (ALS), based on the GlaxoSmithKline first-in-human study NCT00875446 (Berges 2015, Table 2)"
  reference <- "Berges A, Bullman J, Bates S, Krull D, Williams N, Chen C. Ozanezumab Dose Selection for Amyotrophic Lateral Sclerosis by Pharmacokinetic-Pharmacodynamic Modelling of Immunohistochemistry Data from Patient Muscle Biopsies. PLoS ONE. 2015;10(2):e0117355. doi:10.1371/journal.pone.0117355"
  vignette <- "Berges_2015_ozanezumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 76L,
    n_studies      = 1L,
    age_range      = "Adults with ALS (NCT00875446); Berges 2015 reports n=76 enrolled (57 ozanezumab + 19 placebo).",
    weight_range   = "Not reported in Berges 2015; dosing is mg/kg so absolute mg dose depends on individual body weight.",
    sex_female_pct = NULL,
    race_ethnicity = NULL,
    disease_state  = "Amyotrophic lateral sclerosis (ALS); first-in-human single-ascending and repeat-dose study.",
    dose_range     = "Part 1: single IV doses 0.01-15 mg/kg (5 cohorts). Part 2: two IV doses 0.5-15 mg/kg given ~4 weeks apart (3 cohorts). Infusion duration 1 hour.",
    regions        = "Not reported in Berges 2015.",
    notes          = "Study reference: Meininger V et al. (NCT00875446). PK was sampled in all dosed subjects; the Nogo-A co-localization PD endpoint came from skeletal-muscle biopsy IHC in cohorts 7 and 8 (and pre-dose samples from cohorts 5-8). Muscle samples from cohort 3 were unsuitable for IHC due to freeze artefacts (Berges 2015 Methods; Table 1)."
  )

  ini({
    # Structural PK parameters - 2-compartment IV (Berges 2015 Table 2)
    # Paper reports CL=11.7 mL/h, Q=14.6 mL/h, Vc=3310 mL, Vp=3650 mL.
    # Converted to day-units for time and L-units for volumes:
    # CL: 11.7 mL/h * 24 h/day / 1000 mL/L = 0.2808 L/day
    # Q : 14.6 mL/h * 24 h/day / 1000 mL/L = 0.3504 L/day
    # Vc: 3310 mL / 1000 = 3.31 L
    # Vp: 3650 mL / 1000 = 3.65 L
    lcl <- log(0.2808); label("Elimination clearance (CL, L/day)")              # Berges 2015 Table 2: CL = 11.7 mL/h (RSE 4.2%)
    lvc <- log(3.31);   label("Central volume of distribution (Vc, L)")         # Berges 2015 Table 2: Vc = 3310 mL (RSE 4.3%)
    lq  <- log(0.3504); label("Inter-compartmental clearance (Q, L/day)")       # Berges 2015 Table 2: Q = 14.6 mL/h (RSE 8.6%)
    lvp <- log(3.65);   label("Peripheral volume of distribution (Vp, L)")      # Berges 2015 Table 2: Vp = 3650 mL (RSE 5.0%)

    # Effect-compartment delay (Berges 2015 Eq. 3; Table 2)
    # ke0 = 0.00359 1/h * 24 h/day = 0.08616 1/day
    lke0 <- log(0.08616); label("Effect-compartment equilibration rate (ke0, 1/day)") # Berges 2015 Table 2: ke0 = 0.00359 1/h (RSE 14%)

    # PD parameters - sigmoid Emax model on the proportion (percent) of membrane Nogo-A
    # co-localized with ozanezumab (Berges 2015 Eq. 4; Table 2). Asymptote is Emax (100%
    # FIXED); baseline at Ce = 0 is E0. The model interpolates E = E0 + (Emax - E0) * R,
    # where R = Ce^gamma / (EC50^gamma + Ce^gamma), so that E(0) = E0 and E(inf) = Emax.
    le0    <- log(6.04);       label("Baseline % membrane Nogo-A co-localized with drug at Ce=0 (E0, %)")   # Berges 2015 Table 2: E0 = 6.04% (RSE 26%)
    lemax  <- fixed(log(100)); label("Maximal % membrane Nogo-A co-localized with drug (Emax, %)")          # Berges 2015 Table 2: Emax = 100% (Fixed)
    lec50  <- log(24.9);       label("Effect-site concentration giving half-maximal response (EC50, ug/mL)") # Berges 2015 Table 2: EC50 = 24.9 ug/mL (RSE 10%)
    lgamma <- log(1.94);       label("Hill coefficient for sigmoid Emax (gamma, unitless)")                  # Berges 2015 Table 2: gamma = 1.94 (RSE 14%)

    # IIV - between-subject variance on log(CL) and log(Vc) only (Berges 2015 Table 2).
    # No IIV on Q, Vp, ke0, or any PD parameter (paper Methods: "between-subject variance
    # was applied only to the PK parameters CL and Vc"). No off-diagonal covariance reported.
    # omega^2 values from Table 2; CV% in parentheses confirms log-normal interpretation
    # (omega^2 = log(CV^2 + 1) gives 0.0606 for CV=25% and 0.0392 for CV=20%, matching).
    etalcl ~ 0.063        # Berges 2015 Table 2: variance of BSV on CL = 0.063 (CV 25%, RSE 21%)
    etalvc ~ 0.0401       # Berges 2015 Table 2: variance of BSV on Vc = 0.0401 (CV 20%, RSE 20%)

    # Residual error - proportional for PK, additive for PD (Berges 2015 Methods + Table 2).
    propSd          <- 0.254; label("Proportional residual error on Cc (fraction)")               # Berges 2015 Table 2: CV_PK = 25.4% (RSE 5.0%)
    addSd_nogoColoc <- 6.16;  label("Additive residual error on Nogo-A co-localization (%)")      # Berges 2015 Table 2: Sigma PD variance 38.0 (SD 6.16, RSE 56%). See vignette Assumptions for the additional Sigma PDr (replicate variance 89.7, SD 9.47) handled via NONMEM L2 in the original fit but not re-encoded here.
  })
  model({
    # Individual PK parameters
    cl  <- exp(lcl + etalcl)
    vc  <- exp(lvc + etalvc)
    q   <- exp(lq)
    vp  <- exp(lvp)
    ke0 <- exp(lke0)

    # Individual PD parameters
    e0    <- exp(le0)
    emax  <- exp(lemax)
    ec50  <- exp(lec50)
    gamma <- exp(lgamma)

    # Micro-rate constants for explicit 2-compartment ODE form
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 2-compartment IV disposition - dosing into central as IV infusion (Berges 2015 Eqs. 1-2)
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Effect compartment driven by plasma concentration; Ce equilibrates to Cc with rate ke0
    # (Berges 2015 Eq. 3). Effect-site units match Cc units (ug/mL = mg/L).
    Cc <- central / vc
    d/dt(effect) <- ke0 * (Cc - effect)

    # Sigmoid Emax PD on % membrane Nogo-A co-localized with ozanezumab (Berges 2015 Eq. 4).
    # Form is E = E0 + (Emax - E0) * Ce^gamma / (EC50^gamma + Ce^gamma) so that E(0) = E0
    # and the asymptote E(inf) = Emax (= 100%). This matches Table 2's definitions of E0
    # as "baseline level" and Emax as "maximal level" of % co-localization.
    nogoColoc <- e0 + (emax - e0) * effect^gamma / (ec50^gamma + effect^gamma)

    # Observation models
    Cc        ~ prop(propSd)
    nogoColoc ~ add(addSd_nogoColoc)
  })
}
