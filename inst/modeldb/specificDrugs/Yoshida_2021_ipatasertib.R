Yoshida_2021_ipatasertib <- function() {
  description <- paste(
    "Joint parent + metabolite population pharmacokinetic model for oral",
    "ipatasertib (AKT kinase inhibitor under development for breast and",
    "prostate cancer) and its primary active metabolite M1 (G-037720) in",
    "342 adult patients with cancer from five Phase 1 and 2 studies",
    "(Yoshida 2021). Each analyte is described by a 3-compartment",
    "disposition model with sequential zero-order then first-order",
    "absorption from its own depot. The two depots receive the oral",
    "parent dose simultaneously (the user supplies one event per depot",
    "with the same amount and time); both bioavailability anchors are",
    "fixed at F = 1 because absolute parent F and the fraction of parent",
    "metabolised to M1 are not separately identifiable from oral data",
    "alone, and the apparent M1 absorption parameters (kf, Dur, F)",
    "subsume formation, first-pass survival, and metabolite",
    "bioavailability per the source. Retained parent covariates: power",
    "effect of age on apparent CL/F, linear-additive effect of",
    "abiraterone coadministration on apparent CL/F, power effect of body",
    "weight on apparent F, and a +20.1% multiple-dose increment in",
    "apparent F representing CYP3A auto-inhibition by ipatasertib.",
    "Retained metabolite covariates: power effects of body weight on",
    "apparent V3 and Q3 of M1, a +33.1% multiple-dose increment in",
    "apparent F_M1, and an additional +61.5% abiraterone-by-multiple-",
    "dose effect on apparent F_M1. The paper fitted parent and",
    "metabolite in TWO SEPARATE NONMEM runs (Yoshida 2021 Discussion);",
    "this file collapses them into one rxode2 model with NO mechanistic",
    "fractional-conversion linkage, mirroring the paper's simulation",
    "strategy. See vignette Assumptions and deviations.")
  reference <- "Yoshida K, Wilkins J, Winkler J, Wade JR, Kotani N, Wang N, Sane R, Chanu P. Population pharmacokinetics of ipatasertib and its metabolite in cancer patients. J Clin Pharmacol. 2021;61(12):1579-1591. doi:10.1002/jcph.1942"
  vignette <- "Yoshida_2021_ipatasertib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    AGE = list(
      description        = "Subject baseline age in years.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling (AGE/64)^-0.382 on apparent ipatasertib CL/F; reference age 64 years per Yoshida 2021 Methods (median pooled-population age, Table 1). The age effect on M1 was tested in the full model but dropped from the final reduced model.",
      source_name        = "BAGE"
    ),
    WT = list(
      description        = "Subject baseline body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling (WT/75)^-0.617 on apparent ipatasertib F (Table 3 theta_FI,Weight); (WT/75)^0.870 on apparent M1 V3 (Table 4 theta_V3M1,Weight); (WT/75)^0.958 on apparent M1 Q3 (Table 4 theta_Q3M1,Weight). Reference weight 75 kg per Yoshida 2021 Methods (pooled-population median, Table 1). Body weight was tested as a covariate on all CL and V parameters in the full models for both parent and metabolite; only the listed effects were retained after backward elimination.",
      source_name        = "BWT"
    ),
    CONMED_ABI = list(
      description        = "Abiraterone coadministration indicator (1 = subject received abiraterone 1000 mg once daily alongside ipatasertib, 0 = otherwise).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant abiraterone)",
      notes              = "Linear-additive effect of -18.5% on apparent ipatasertib CL/F (Table 3 theta_CLI,Abi = -0.185); linear-additive effect of +61.5% on apparent M1 F applied only at the multiple-dose state (Table 4 theta_FM1,Abi = 0.615; the source NONMEM control stream enforced ABIRATER * SSFLAG as the indicator on this term). Of 342 subjects, 189 (55.3%) received abiraterone: all 183 from study GO27983, six from study JO29655, none from PAM4743g / GO29227 / PAM4983g.",
      source_name        = "ABIRATER"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 342L,
    n_studies       = 5L,
    age_range       = "26 to 88 years; median 64 years (geometric mean 62.6) per paper Table 1.",
    weight_range    = "41.5 to 160 kg; median 75 kg (geometric mean 76.1) per paper Table 1 (2 missing values).",
    sex_female_pct  = 33.3,
    race_ethnicity  = c(White = 76.6, `Black or African American` = 3.22, Asian = 14.9, Multiple = 0.292, Other = 0.877, Unknown = 4.09),
    disease_state   = "Adult patients with locally advanced or metastatic solid tumours (prostate 56.4%, breast 26.6%, other 17%) from five Phase 1 and 2 ipatasertib clinical studies.",
    dose_range      = "Oral ipatasertib 25 to 800 mg once daily on days 1-21 of 28-day cycles (or continuously in the abiraterone-combination studies); 400 mg is the dose under Phase 3 investigation. The 200 mg cohort (n = 96) came predominantly from GO27983, the 400 mg cohort (n = 188) from GO27983 and GO29227.",
    regions         = "Multinational (United States, Europe, and Japan). Studies: PAM4743g (NCT01090960), JO29655 (Japanese cohort), PAM4983g (NCT01362374 arm C), GO27983 (NCT01485861 A.MARTIN), GO29227 (NCT02162719 LOTUS).",
    n_observations  = "3050 ipatasertib observations and 2050 M1 observations across the five studies.",
    hepatic_impairment = c(Normal = 78.7, Mild = 20.2, Moderate = 0.585, Missing = 0.585),
    egfr_range      = "39.4 to 212 mL/min/1.73 m^2; median 92 (geometric mean 94.4).",
    notes           = "Mild and moderate renal impairment, mild hepatic impairment, and race were tested but not identified as significant covariates in the final reduced models per Yoshida 2021 Discussion. Ipatasertib was fitted in NONMEM 7.4.3 with FOCE-I (run230 reported MINIMIZATION SUCCESSFUL; the standard 'however, problems occurred with the minimization' caveat appeared but the covariance step ran successfully). The M1 model was fitted in a separate NONMEM run; the on-disk control stream is in Supplementary Text S2 and the final estimates come from paper Table 4."
  )

  ini({
    # ------------------------------------------------------------------
    # PARENT (ipatasertib) STRUCTURAL PARAMETERS
    # ------------------------------------------------------------------
    # All parent values come from the run230.lst FINAL PARAMETER
    # ESTIMATE block; they round to the paper Table 3 values to within
    # the precision the paper reports. The NONMEM listing was the
    # operator-provided authoritative source.

    lcl <- log(162)
    label("Ipatasertib apparent clearance CL/F (L/h)")                # run230.lst FINAL PARAMETER ESTIMATE TH1 = 1.62E+02; Table 3

    lvc <- log(1230)
    label("Ipatasertib apparent central volume V2/F (L)")             # run230.lst FINAL PARAMETER ESTIMATE TH2 = 1.23E+03; Table 3

    lvp <- log(2590)
    label("Ipatasertib apparent first peripheral volume V3/F (L)")    # run230.lst FINAL PARAMETER ESTIMATE TH3 = 2.59E+03; Table 3

    lvp2 <- log(4340)
    label("Ipatasertib apparent second peripheral volume V4/F (L)")   # run230.lst FINAL PARAMETER ESTIMATE TH4 = 4.34E+03; Table 3

    lq <- log(76.6)
    label("Ipatasertib apparent inter-compartmental clearance Q3/F (L/h)")  # run230.lst FINAL PARAMETER ESTIMATE TH5 = 7.66E+01; Table 3

    lq2 <- log(3.26)
    label("Ipatasertib apparent inter-compartmental clearance Q4/F (L/h)")  # run230.lst FINAL PARAMETER ESTIMATE TH6 = 3.26E+00; Table 3

    lka <- log(1.84)
    label("Ipatasertib first-order absorption rate ka (1/h)")          # run230.lst FINAL PARAMETER ESTIMATE TH7 = 1.84E+00; Table 3

    ldurdepot <- log(0.348)
    label("Ipatasertib zero-order absorption duration Dur (h)")        # run230.lst FINAL PARAMETER ESTIMATE TH9 = 3.48E-01; Table 3

    lfdepot <- fixed(log(1))
    label("Ipatasertib bioavailability anchor F (FIXED at 1; absolute F not identifiable from oral data)")  # run230.lst $THETA TH8 = 1 FIX; Yoshida 2021 Methods

    # ------------------------------------------------------------------
    # PARENT (ipatasertib) COVARIATE EFFECTS
    # ------------------------------------------------------------------

    e_md_fdepot <- 0.201
    label("Linear coefficient for multiple-dose state on ipatasertib F (unitless)")  # run230.lst FINAL PARAMETER ESTIMATE TH10 = 2.01E-01; Table 3 theta_FI,MD

    e_age_cl <- -0.382
    label("Power exponent of (AGE/64) on ipatasertib CL (unitless)")   # run230.lst FINAL PARAMETER ESTIMATE TH13 = -3.82E-01; Table 3 theta_CLI,Age

    e_abi_cl <- -0.185
    label("Linear coefficient for CONMED_ABI on ipatasertib CL (unitless)")  # run230.lst FINAL PARAMETER ESTIMATE TH20 = -1.85E-01; Table 3 theta_CLI,Abi

    e_wt_fdepot <- -0.617
    label("Power exponent of (WT/75) on ipatasertib F (unitless)")     # run230.lst FINAL PARAMETER ESTIMATE TH34 = -6.17E-01; Table 3 theta_FI,Weight

    # ------------------------------------------------------------------
    # PARENT (ipatasertib) INTER-INDIVIDUAL VARIABILITY
    # ------------------------------------------------------------------
    # OMEGA values are the variances reported in the .lst OMEGA matrix
    # (also in Table 3 as omega^2). The percent-CV the paper prints is
    # sqrt(omega^2) (not the log-normal CV formula), so the variance
    # IS the value to assign to the eta block in nlmixr2.

    etalcl       ~ 0.0729                                              # run230.lst FINAL PARAMETER ESTIMATE OMEGA(1,1) = 7.29E-02; Table 3
    etalka       ~ 1.74                                                # run230.lst FINAL PARAMETER ESTIMATE OMEGA(7,7) = 1.74E+00; Table 3
    etalfdepot   ~ 0.139                                               # run230.lst FINAL PARAMETER ESTIMATE OMEGA(8,8) = 1.39E-01; Table 3
    etaldurdepot ~ 2.23                                                # run230.lst FINAL PARAMETER ESTIMATE OMEGA(9,9) = 2.23E+00; Table 3

    # ------------------------------------------------------------------
    # PARENT (ipatasertib) RESIDUAL ERROR
    # ------------------------------------------------------------------
    # The source NONMEM control stream fits log-transformed observations
    # with additive normal error of variance SIGMA(1,1) = 0.248 (Table 3
    # sigma^2_I). This is mathematically a log-normal residual on
    # linear-space concentration with SD on the log scale = sqrt(0.248).

    expSd <- sqrt(0.248)
    label("Ipatasertib log-normal residual SD on log-concentration (unitless)")  # run230.lst FINAL PARAMETER ESTIMATE SIGMA = 2.48E-01 variance; Table 3 sigma^2_I

    # ------------------------------------------------------------------
    # METABOLITE (M1 / G-037720) STRUCTURAL PARAMETERS
    # ------------------------------------------------------------------
    # All M1 values come from paper Table 4 (final reduced model for M1
    # PK). No NONMEM listing file is on disk for the M1 fit; the source
    # NONMEM control stream is reproduced in Supplementary Text S2 (the
    # initial-value THETAs there are run-development starting values
    # that differ slightly from the Table 4 final estimates, so Table 4
    # is the authoritative source).

    lcl_m1 <- log(314)
    label("M1 apparent clearance CL/F (L/h)")                          # Yoshida 2021 Table 4 CL_M1

    lvc_m1 <- log(456)
    label("M1 apparent central volume V2/F (L)")                       # Yoshida 2021 Table 4 V2_M1

    lvp_m1 <- log(7270)
    label("M1 apparent first peripheral volume V3/F (L)")              # Yoshida 2021 Table 4 V3_M1

    lvp2_m1 <- log(14800)
    label("M1 apparent second peripheral volume V4/F (L)")             # Yoshida 2021 Table 4 V4_M1

    lq_m1 <- log(219)
    label("M1 apparent inter-compartmental clearance Q3/F (L/h)")      # Yoshida 2021 Table 4 Q3_M1

    lq2_m1 <- log(9.04)
    label("M1 apparent inter-compartmental clearance Q4/F (L/h)")      # Yoshida 2021 Table 4 Q4_M1

    lka_m1 <- log(0.191)
    label("M1 apparent first-order formation rate kf (1/h, treated as apparent ka)")  # Yoshida 2021 Table 4 kf_M1

    ldurdepot_m1 <- log(0.730)
    label("M1 apparent zero-order formation duration Dur (h, treated as apparent zero-order absorption duration)")  # Yoshida 2021 Table 4 Dur_M1

    lfdepot_m1 <- fixed(log(1))
    label("M1 bioavailability anchor F (FIXED at 1; subsumes fraction-formed and first-pass)")  # Yoshida 2021 Methods (M1 absorption parameters are apparent)

    # ------------------------------------------------------------------
    # METABOLITE (M1) COVARIATE EFFECTS
    # ------------------------------------------------------------------

    e_md_fdepot_m1 <- 0.331
    label("Linear coefficient for multiple-dose state on M1 F (unitless)")   # Yoshida 2021 Table 4 theta_FM1,MD

    e_abi_md_fdepot_m1 <- 0.615
    label("Additional linear coefficient for CONMED_ABI x multiple-dose on M1 F (unitless; active only when MULTI_DOSE = 1)")  # Yoshida 2021 Table 4 theta_FM1,Abi

    e_wt_vp_m1 <- 0.870
    label("Power exponent of (WT/75) on M1 V3 (unitless)")             # Yoshida 2021 Table 4 theta_V3M1,Weight

    e_wt_q_m1 <- 0.958
    label("Power exponent of (WT/75) on M1 Q3 (unitless)")             # Yoshida 2021 Table 4 theta_Q3M1,Weight

    # ------------------------------------------------------------------
    # METABOLITE (M1) INTER-INDIVIDUAL VARIABILITY
    # ------------------------------------------------------------------

    etalcl_m1       ~ 0.154                                            # Yoshida 2021 Table 4 omega^2(CL_M1)
    etalvc_m1       ~ 0.657                                            # Yoshida 2021 Table 4 omega^2(V2_M1)
    etalfdepot_m1   ~ 0.297                                            # Yoshida 2021 Table 4 omega^2(F_M1)
    etaldurdepot_m1 ~ 1.44                                             # Yoshida 2021 Table 4 omega^2(Dur_M1)

    # ------------------------------------------------------------------
    # METABOLITE (M1) RESIDUAL ERROR
    # ------------------------------------------------------------------

    expSd_m1 <- sqrt(0.228)
    label("M1 log-normal residual SD on log-concentration (unitless)")  # Yoshida 2021 Table 4 sigma^2_M1
  })

  model({
    # Reference values for power / linear effect equations
    # (Yoshida 2021 Methods, paper Table 1 medians).
    ref_age <- 64    # years
    ref_wt  <- 75    # kg

    # Multiple-dose state flag. Paper Simulations section: 'this change
    # was assumed to occur 48 hours after the first dose.' We encode
    # this directly from simulation time t. The convention is that the
    # user's event table places the FIRST dose at t = 0; see vignette
    # Assumptions and deviations.
    mdflag <- (t > 48) * 1.0

    # ------------------------------------------------------------------
    # PARENT (ipatasertib) individual parameters
    # ------------------------------------------------------------------
    cl       <- exp(lcl + etalcl) *
                (AGE / ref_age)^e_age_cl *
                (1 + e_abi_cl * CONMED_ABI)
    vc       <- exp(lvc)
    vp       <- exp(lvp)
    vp2      <- exp(lvp2)
    q        <- exp(lq)
    q2       <- exp(lq2)
    ka       <- exp(lka + etalka)
    durdepot <- exp(ldurdepot + etaldurdepot)
    fdepot   <- exp(lfdepot + etalfdepot) *
                (WT / ref_wt)^e_wt_fdepot *
                (1 + e_md_fdepot * mdflag)

    # ------------------------------------------------------------------
    # METABOLITE (M1) individual parameters
    # ------------------------------------------------------------------
    cl_m1       <- exp(lcl_m1 + etalcl_m1)
    vc_m1       <- exp(lvc_m1 + etalvc_m1)
    vp_m1       <- exp(lvp_m1)  * (WT / ref_wt)^e_wt_vp_m1
    vp2_m1      <- exp(lvp2_m1)
    q_m1        <- exp(lq_m1)   * (WT / ref_wt)^e_wt_q_m1
    q2_m1       <- exp(lq2_m1)
    ka_m1       <- exp(lka_m1)
    durdepot_m1 <- exp(ldurdepot_m1 + etaldurdepot_m1)
    fdepot_m1   <- exp(lfdepot_m1 + etalfdepot_m1) *
                   (1 + e_md_fdepot_m1     * mdflag) *
                   (1 + e_abi_md_fdepot_m1 * mdflag * CONMED_ABI)

    # Micro-constants
    kel  <- cl  / vc
    k12  <- q   / vc
    k21  <- q   / vp
    k13  <- q2  / vc
    k31  <- q2  / vp2

    kel_m1 <- cl_m1 / vc_m1
    k12_m1 <- q_m1  / vc_m1
    k21_m1 <- q_m1  / vp_m1
    k13_m1 <- q2_m1 / vc_m1
    k31_m1 <- q2_m1 / vp2_m1

    # ------------------------------------------------------------------
    # ODE system. The user supplies TWO dose records per oral ipatasertib
    # administration: one to cmt = depot (parent absorption) and one to
    # cmt = depot_m1 (apparent M1 "absorption" that subsumes formation +
    # first-pass + metabolite bioavailability; the NONMEM source models
    # the metabolite kinetics with its own kf, Dur, and F driven from
    # the parent dose). Both bioavailability anchors are F = 1 with
    # multiplicative covariate effects; the absolute parent F and the
    # fraction of parent metabolised to M1 are not separately
    # identifiable from oral data alone, so a single mechanistically-
    # coupled formation pathway was not estimated by the source authors
    # (Yoshida 2021 Discussion).
    # ------------------------------------------------------------------
    d/dt(depot)        <- -ka * depot
    d/dt(central)      <-  ka * depot - kel * central -
                            k12 * central + k21 * peripheral1 -
                            k13 * central + k31 * peripheral2
    d/dt(peripheral1)  <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2)  <-  k13 * central - k31 * peripheral2

    d/dt(depot_m1)        <- -ka_m1 * depot_m1
    d/dt(central_m1)      <-  ka_m1 * depot_m1 - kel_m1 * central_m1 -
                                k12_m1 * central_m1 + k21_m1 * peripheral1_m1 -
                                k13_m1 * central_m1 + k31_m1 * peripheral2_m1
    d/dt(peripheral1_m1)  <-  k12_m1 * central_m1 - k21_m1 * peripheral1_m1
    d/dt(peripheral2_m1)  <-  k13_m1 * central_m1 - k31_m1 * peripheral2_m1

    # Absorption modifiers: each depot receives sequential zero-order
    # input over its own duration, then transfers to its central
    # compartment at its own first-order rate. F is applied to the depot.
    dur(depot)    <- durdepot
    f(depot)      <- fdepot

    dur(depot_m1) <- durdepot_m1
    f(depot_m1)   <- fdepot_m1

    # Observations. Volumes are reported in L and doses in mg; multiplying
    # by 1000 converts mg/L to ng/mL (the unit the paper reports). The
    # source NONMEM control stream uses S2 = V2 / 1000 to achieve the same
    # scaling for the parent.
    Cc    <- central    / (vc    / 1000)
    Cc_m1 <- central_m1 / (vc_m1 / 1000)

    Cc    ~ lnorm(expSd)
    Cc_m1 ~ lnorm(expSd_m1)
  })
}
