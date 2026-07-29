Rekic_2011_atazanavir <- function() {
  description <- "Population PK / PD model for atazanavir (boosted with ritonavir 100 mg QD) and its concentration-dependent effect on plasma bilirubin in adult antiretroviral-naive HIV-positive patients from the NORTHIV trial (Rekic 2011). Atazanavir disposition is described by a one-compartment model with first-order absorption and an absorption lag, fitted to log-transformed plasma atazanavir concentrations; ka and the lag time were fixed to the published values from the Colombo 2006 atazanavir popPK report (ref 27) because sparse absorption-phase sampling did not support their re-estimation, and CL/F and V/F were re-estimated under fixed allometric scaling on body weight centred at 70 kg (exponents 0.75 on CL/F and 1 on V/F, both fixed a priori). The bilirubin response is described by an indirect-response (turnover) model with concentration-dependent inhibition of the fractional turnover rate kout: dB/dt = kin - kout * (1 - Imax * Cc / (IC50 + Cc)) * B, with kin re-parameterised at steady state as kin = kout * Baseline. Inter-individual variability is supported only on V/F, CL/F (PK), and bilirubin baseline (PD); the paper notes that the data did not support IIV on the remaining PD parameters. PK residual variability is proportional; bilirubin residual variability is combined additive + proportional (the paper's 'slope-intercept' model)."
  reference <- "Rekic D, Clewe O, Roshammar D, Flamholc L, Sonnerborg A, Ormaasen V, Gisslen M, Abelo A, Ashton M. Bilirubin -- a potential marker of drug exposure in atazanavir-based antiretroviral therapy. The AAPS Journal. 2011;13(4):598-605. doi:10.1208/s12248-011-9299-0. Absorption-phase parameters (ka, lag time) fixed from Colombo S, Buclin T, Cavassini M, et al. Population pharmacokinetics of atazanavir in patients with human immunodeficiency virus infection. Antimicrob Agents Chemother. 2006;50(11):3801-3808 (cited as ref 27)."
  vignette <- "Rekic_2011_atazanavir"
  units <- list(time = "h", dosing = "mg", concentration = "umol/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed allometric covariate centred at the population median of 70 kg (Rekic 2011 Methods Eq. 1). Used to scale CL/F with a fixed exponent of 0.75 and V/F with a fixed exponent of 1.",
      source_name        = "BW"
    )
  )

  population <- list(
    species       = "human",
    n_subjects    = 82L,
    n_studies     = 1L,
    n_observations = "200 atazanavir steady-state plasma samples and 361 bilirubin observations (Rekic 2011 Methods 'The Data').",
    age_range     = "adults (specific age distribution not tabulated in the paper)",
    weight_range  = "median body weight 70 kg (Rekic 2011 Methods, reference body weight for the allometric scaling); individual weight range not tabulated in the paper",
    sex_female_pct = NA_real_,
    disease_state = "Antiretroviral-naive HIV-1-seropositive adults with normal baseline bilirubin (< 20 umol/L) starting atazanavir/ritonavir-based combination antiretroviral therapy in the NORTHIV trial.",
    dose_range    = "300 mg atazanavir + 100 mg ritonavir orally once daily; nucleoside reverse transcriptase inhibitor backbone allowed to vary per Swedish national guidelines, with a few protocol-permitted dose adjustments for clinical practice.",
    regions       = "Sweden and Norway (NORTHIV trial; multicentre, ICH-GCP, declaration of Helsinki).",
    notes         = "Atazanavir plasma sampling at trial weeks 4, 12, 48, 96, and 144; bilirubin measured at baseline and at the matching atazanavir time points. Baseline bilirubin mean (+/- SD) was 7.8 (+/- 3.3) umol/L; new steady-state bilirubin reached 34 (+/- 18.6) umol/L on therapy."
  )

  ini({
    # ====================================================================
    # Atazanavir population PK -- one-compartment with first-order
    # absorption and lag (Rekic 2011 Table I, 'PK model' rows).
    # ka and lag time were fixed from Colombo 2006 (ref 27) because too
    # few absorption-phase samples were available; CL/F and V/F were
    # re-estimated in NORTHIV (Rekic 2011 Methods 'Pharmacokinetic and
    # Pharmacodynamic Model Building').
    # ====================================================================
    ltlag <- fixed(log(0.96))
    label("Absorption lag time of atazanavir (h; fixed from Colombo 2006, ref 27)")  # Rekic 2011 Table I PK row: Lag time = 0.96 h, footnote a 'Fixed according to (27)'

    lka <- fixed(log(3.4))
    label("First-order absorption rate constant of atazanavir (1/h; fixed from Colombo 2006, ref 27)")  # Rekic 2011 Table I PK row: ka = 3.4 1/h, footnote a 'Fixed according to (27)'

    lvc <- log(93.6)
    label("Apparent central volume of distribution V/F of atazanavir at the reference 70 kg subject (L)")  # Rekic 2011 Table I PK row: V/F = 93.6 L (95% CI 62-125), RSE not reported

    lcl <- log(6.47)
    label("Apparent clearance CL/F of atazanavir at the reference 70 kg subject (L/h)")  # Rekic 2011 Table I PK row: CL/F = 6.47 L/h (95% CI 5.39-7.55)

    # Allometric exponents on body weight, fixed a priori to canonical
    # values per Rekic 2011 Methods Eq. 1 ('The scaling factor was a
    # priori set to 0.75 for clearance (CL/F) and to 1 for volume of
    # distribution (V/F)').
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent on CL/F (unitless; fixed a priori)")  # Rekic 2011 Methods 'Pharmacokinetic and Pharmacodynamic Model Building': "scaling factor was a priori set to 0.75 for clearance (CL/F)"
    e_wt_vc <- fixed(1)
    label("Allometric exponent on V/F (unitless; fixed a priori)")   # Rekic 2011 Methods 'Pharmacokinetic and Pharmacodynamic Model Building': "and to 1 for volume of distribution (V/F)"

    # IIV on the PK parameters. Rekic 2011 Table I reports between-subject
    # variability as CV%; convert to log-scale variance via omega^2 =
    # log(1 + CV^2) per the lognormal IIV convention.
    etalvc ~ log(1 + 0.531^2)  # Rekic 2011 Table I PK IIV row: omega_V/F = 53.1% CV (RSE 43.6%)
    etalcl ~ log(1 + 0.438^2)  # Rekic 2011 Table I PK IIV row: omega_CL/F = 43.8% CV (RSE 19.5%)

    # Atazanavir residual error -- proportional in linear (umol/L) space.
    # Rekic 2011 Table I PD row prints the PK sigma_prop (51.0%) ahead of
    # the bilirubin parameters in the same row; the in-text Methods
    # sentence "Proportional and slope intercept models were applied to
    # explain the residual variability" makes the assignment of the
    # proportional model to atazanavir (PK) and the slope-intercept model
    # to bilirubin (PD) unambiguous.
    propSd <- 0.510
    label("Proportional residual error on atazanavir concentration (fraction)")  # Rekic 2011 Table I 'PD model' row, sigma_prop column: 51.0% (95% CI 42.7-59.3%)

    # ====================================================================
    # Bilirubin population PD -- indirect-response model with
    # concentration-dependent inhibition of the fractional turnover rate
    # kout (Rekic 2011 Methods Eqs. 2 and 3, Fig. 1b; Table I PD rows).
    #
    #     E_inh = Imax * Cc / (IC50 + Cc)
    #     dB/dt = kin - kout * (1 - E_inh) * B
    #
    # At baseline (Cc = 0, E_inh = 0), dB/dt = 0 implies kin = kout *
    # Baseline; this steady-state relationship is enforced in model().
    # ====================================================================
    lrbase_bilirubin <- log(7.69)
    label("Predose typical bilirubin baseline (umol/L)")  # Rekic 2011 Table I PD row, Baseline: 7.69 umol/L (95% CI 6.99-8.39)

    lkout_bilirubin <- log(0.420)
    label("First-order bilirubin elimination rate constant kout (1/h)")  # Rekic 2011 Table I PD row, kout: 0.420 1/h (95% CI 0.36-0.48); uninhibited half-life ln(2)/0.420 = 1.65 h matches Rekic 2011 Results "Uninhibited, bilirubin displayed a half-life of 1.64 h"

    logitimax_bilirubin <- qlogis(0.910)
    label("Logit of maximum fractional inhibition Imax of bilirubin elimination by atazanavir (unitless; bounded to [0,1])")  # Rekic 2011 Table I PD row, Imax: 91.0% (95% CI 87-94%); logit transform enforces the paper-stated 0 <= Imax <= 1 constraint

    lic50_bilirubin <- log(0.300)
    label("Atazanavir concentration resulting in half-maximal Imax inhibition of bilirubin kout (umol/L)")  # Rekic 2011 Table I PD row, IC50: 0.300 umol/L (95% CI 0.24-0.37)

    # IIV on the PD parameters. Rekic 2011 Discussion: "The data
    # available could only support inter-individual variability on
    # bilirubin baseline in the pharmacodynamic model." kout, Imax, and
    # IC50 therefore carry no eta.
    etalrbase_bilirubin ~ log(1 + 0.326^2)  # Rekic 2011 Table I PD IIV row: omega_Baseline = 32.6% CV (RSE 20.2%)

    # Bilirubin residual error -- combined additive + proportional
    # ('slope-intercept' per Rekic 2011 Methods). Rekic 2011 Table I PD
    # 'Residual error' columns: sigma_prop 39.4% (95% CI 35.5-43.3%) and
    # sigma_add 2.39 umol/L (95% CI 1.96-2.82).
    propSd_bilirubin <- 0.394
    label("Proportional residual error on plasma bilirubin (fraction)")  # Rekic 2011 Table I PD row, sigma_prop (Residual error column): 39.4%
    addSd_bilirubin <- 2.39
    label("Additive residual error on plasma bilirubin (umol/L)")        # Rekic 2011 Table I PD row, sigma_add (Residual error column): 2.39 umol/L
  })

  model({
    # ====================================================================
    # Atazanavir molecular weight for the mg -> umol/L concentration
    # conversion. Rekic 2011 reports concentrations and IC50 in umol/L,
    # while the V/F estimate is in L and the user-facing dose is in mg.
    # Atazanavir base (the labelled mass on the 300 mg Reyataz capsule)
    # is C38H52N6O7 with MW = 704.86 g/mol = 704.86 mg/mmol; so
    # Cc [umol/L] = (central_amount [mg] / vc [L]) * (1000 umol / 704.86 mg).
    # ====================================================================
    mw_atz <- 704.86  # g/mol = mg/mmol; atazanavir base

    # ----- Individual PK parameters (allometric scaling on WT, fixed
    # canonical exponents per Rekic 2011 Methods Eq. 1) ---------------
    tlag <- exp(ltlag)
    ka   <- exp(lka)
    vc   <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    cl   <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    kel  <- cl / vc

    # ----- Individual PD parameters --------------------------------
    bilirubin_base <- exp(lrbase_bilirubin + etalrbase_bilirubin)
    kout_bilirubin <- exp(lkout_bilirubin)
    imax_bilirubin <- expit(logitimax_bilirubin)
    ic50_bilirubin <- exp(lic50_bilirubin)
    # Steady-state synthesis: at baseline (Cc = 0, no inhibition), the
    # turnover equation dB/dt = kin - kout * B = 0 implies kin = kout *
    # Baseline (Rekic 2011 Methods Eq. 2).
    kin_bilirubin <- kout_bilirubin * bilirubin_base

    # ----- PK ODE system: one-compartment, first-order absorption with
    # lag (Rekic 2011 Methods 'Pharmacokinetic and Pharmacodynamic Model
    # Building') --------------------------------------------------------
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central
    alag(depot)   <- tlag

    # Atazanavir plasma concentration in umol/L (paper units).
    Cc <- (central / vc) * 1000 / mw_atz

    # ----- Bilirubin indirect-response ODE (Rekic 2011 Methods Eqs.
    # 2 and 3; Fig. 1b) ----------------------------------------------
    # E_inh is the fractional inhibition of kout, bounded to [0, Imax].
    e_inh <- imax_bilirubin * Cc / (ic50_bilirubin + Cc)
    d/dt(effect) <- kin_bilirubin - kout_bilirubin * (1 - e_inh) * effect
    effect(0)    <- bilirubin_base

    # Bilirubin observation (umol/L) and residual error model.
    bilirubin <- effect

    Cc        ~ prop(propSd)
    bilirubin ~ add(addSd_bilirubin) + prop(propSd_bilirubin)
  })
}
