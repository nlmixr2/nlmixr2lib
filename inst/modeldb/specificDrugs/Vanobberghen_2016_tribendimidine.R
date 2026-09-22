Vanobberghen_2016_tribendimidine <- function() {
  description <- "Joint population PK model for the two tribendimidine metabolites dADT (deacetylated amidantel, anthelminthically active) and adADT (acetylated dADT, marginally active) in Opisthorchis viverrini-infected Lao adults given single oral doses of 25-600 mg (Vanobberghen 2016). A six-transit-compartment absorption chain feeds a one-compartment dADT disposition model; a fixed 65% of dADT elimination is routed to a one-compartment adADT model (the remaining 35% is assumed renal). Allometric body-weight scaling (fixed 0.75 on clearances, 1 on volumes, reference 51.5 kg), a linear age effect on both clearances, and a 200-mg-tablet formulation effect on mean transit time and on both volumes. Full 4x4 variance-covariance block across the two clearances and two volumes. Fitted on natural-log-transformed molar concentrations, so amounts are nmol and concentrations nmol/L."
  reference <- "Vanobberghen F, Penny MA, Duthaler U, Odermatt P, Sayasone S, Keiser J, Tarning J. Population pharmacokinetic modeling of tribendimidine metabolites in Opisthorchis viverrini-infected adults. Antimicrob Agents Chemother. 2016;60(10):5695-5704. doi:10.1128/AAC.00655-16"
  vignette <- "Vanobberghen_2016_tribendimidine"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against the $MODEL block comments of the final
  # NONMEM control stream (supplemental File S1), which names COMP(1) as the
  # dose compartment, COMP(2) as dADT central, COMP(3) as adADT central and
  # COMP(4)-COMP(9) as transit compartments 1-6.
  compartmentData <- list(
    depot = list(
      analyte = "tribendimidine (dosed prodrug; never measured)",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    transit1 = list(
      analyte = "tribendimidine / dADT in transit",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    transit2 = list(
      analyte = "tribendimidine / dADT in transit",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    transit3 = list(
      analyte = "tribendimidine / dADT in transit",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    transit4 = list(
      analyte = "tribendimidine / dADT in transit",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    transit5 = list(
      analyte = "tribendimidine / dADT in transit",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    transit6 = list(
      analyte = "tribendimidine / dADT in transit",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "dADT (deacetylated amidantel)",
      units = "nmol",
      specimen = "whole blood",
      verified = TRUE
    ),
    central_adadt = list(
      analyte = "adADT (acetylated dADT)",
      units = "nmol",
      specimen = "whole blood",
      verified = TRUE
    )
  )

  covariateData <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters both clearances as the linear form (1 + slope * (AGE - 52)). Two independent places in the source give the centering value as 52 years: Vanobberghen 2016 Table 1 footnote a ('Results shown are for a typical patient aged 52 years, weighing 51.5 kg, and receiving 50-mg tablets'), and the final NONMEM control stream (supplemental File S1, CMAGE and CPAGE definition blocks), which reads `(AGE - 52.00)`. Materials and Methods loosely states that continuous covariates were centered on the median and the Results give a median age of 42 years (IQR 32-47), but the table footnote defining what the printed estimates mean says 52, so 52 is the age the published CL/F of 16.7 L/h refers to. Because the form is linear rather than exponential, it is only valid over the ages actually studied; extrapolating adADT clearance beyond roughly 99 years would make it non-positive.",
      source_name = "AGE"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling on both clearances (exponent fixed at 0.75) and both volumes (exponent fixed at 1), normalized to 51.5 kg. NOT listed as a row in Vanobberghen 2016 Table 1 because the exponents were fixed rather than estimated, but confirmed by two other places in the source: Table 1 footnote a states the printed estimates are for 'a typical patient aged 52 years, weighing 51.5 kg', and the final NONMEM control stream (supplemental File S1) applies `(WEIGHT/51.50)**0.75` and `(WEIGHT/51.50)**1.00` to every clearance and volume. The cohort median weight reported in the Results is 52 kg (IQR 47-57).",
      source_name = "WEIGHT"
    ),
    FORM_TRI_TAB200 = list(
      description = "Tribendimidine tablet-strength formulation indicator (1 = 200-mg enteric-coated tablets, as used in study 1; 0 = 50-mg enteric-coated tablets, as used in study 2)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (50-mg enteric-coated tablets, study 2)",
      notes = "Both arms are enteric-coated tablets from the same manufacturer, differing only in tablet strength; the 200-mg tablets were observed in vitro to float and hence to delay absorption (Vanobberghen 2016 Discussion, citing reference 27). Multiplicative effects relative to the 50-mg reference: +40.1% on mean transit time, +113% on dADT central volume and +364% on adADT central volume (Table 1). The control stream encodes this as `IF(STUDY.EQ.1) ... (1 + THETA)` with study 2 as the reference level, so the indicator is equivalent to a study indicator for this data set; it is named for the formulation because that is the mechanism the authors attribute the effect to. The 25-mg dose was given as split 50-mg tablets, which destroyed the enteric coating; that sub-level is NOT captured by this covariate (the authors' interaction model for split tablets did not converge).",
      source_name = "STUDY"
    )
  )

  # Screened during covariate model building but not retained in the final
  # model, so they are documented rather than implemented.
  covariatesDataExcluded <- list(
    CRCL = list(
      description = "Creatinine clearance estimated with the CKD-EPI equation",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      notes = "Evaluated as a linear covariate centered on the median in the stepwise procedure (Vanobberghen 2016 Materials and Methods) but not retained in the final model. Cohort median 66 mL/min/1.73 m^2 (IQR 50-112)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Evaluated in the stepwise covariate procedure but not retained. 35 of 68 participants (51%) were female."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 68L,
    n_studies = 2L,
    age_range = "median 42 years (IQR 32-47)",
    weight_range = "median 52 kg (IQR 47-57)",
    sex_female_pct = 51.5,
    renal_function = "Median creatinine clearance 66 mL/min/1.73 m^2 (IQR 50-112), CKD-EPI equation",
    disease_state = "Adults with confirmed Opisthorchis viverrini infection; median baseline egg burden 897 eggs per gram of stool (IQR 437-1,817). Cure, defined as no eggs detected at 21 days post-treatment, ranged from 11% at 25 mg to 100% at 400 mg.",
    dose_range = "Single oral doses of 25, 50, 100, 200, 400 and 600 mg tribendimidine. Study 1 gave 200, 400 and 600 mg using 200-mg enteric-coated tablets (n = 13, 9, 9); study 2 gave 25, 50, 100 and 200 mg using 50-mg enteric-coated tablets (n = 9, 9, 9, 10). The 25-mg dose was given as split 50-mg tablets.",
    regions = "Champasack Province, Lao People's Democratic Republic (Champasack Provincial Hospital, Pakse)",
    studies = "Two phase IIa open-label randomized ascending-dose-finding trials conducted in November 2012 (study 1) and October 2013 (study 2); ISRCTN Registry no. ISRCTN96948551.",
    notes = "Baseline demographics from Vanobberghen 2016 Results, first paragraph. 1,307 samples were analyzed for dADT and 1,303 for adADT, comprising 300 venous whole-blood, 669 plasma and 338 capillary dried-blood-spot samples; whole blood was collected in study 1 only. The three matrices were modeled jointly because a proportional transformation factor between matrices did not improve the fit. Assay LLOQ was 1 ng/mL for whole blood and plasma and for study 2 dried blood spots, and 10 ng/mL for study 1 dried blood spots; observations below the limit were handled with the M3 likelihood method."
  )

  ini({
    # Absorption. The chain is a dose compartment followed by six transit
    # compartments, i.e. seven first-order transfers in series, so the transit
    # rate constant is KTR = (n + 1) / MTT with n = 6 fixed (supplemental File
    # S1 $PK: `NN = 6` then `KTR = (NN+1)/MT`, with K14, K45, ... K92 all set
    # to KTR). The number of compartments is structural and is expressed by the
    # length of the ODE chain below rather than by an ini() parameter.
    lmtt <- log(3.38); label("Mean transit time MTT for the 50-mg-tablet reference formulation (h)") # Table 1, 'Mean transit time (h)' = 3.38 (95% CI 2.51-4.63)

    # dADT (deacetylated amidantel) disposition. Apparent values; relative
    # bioavailability was fixed to unity for the population.
    lcl <- log(16.7); label("Apparent dADT clearance CL/F at 51.5 kg and 52 years (L/h)") # Table 1, dADT 'CL/F (liters/h)' = 16.7 (95% CI 14.6-18.1)
    lvc <- log(93.3); label("Apparent dADT central volume Vc/F at 51.5 kg, 50-mg-tablet reference (L)") # Table 1, dADT 'V_C/F (liters)' = 93.3 (95% CI 74.3-112)

    # adADT (acetylated dADT) disposition. These are apparent values that also
    # absorb the assumed 65% metabolic fraction, so they are CL/(F * fm) and
    # V/(F * fm); the paper notes the fixed renal fraction 'simply acts as a
    # scaling factor for the pharmacokinetic parameters' (Discussion).
    lcl_adadt <- log(41.8); label("Apparent adADT clearance CL/F at 51.5 kg and 52 years (L/h)") # Table 1, adADT 'CL/F (liters/h)' = 41.8 (95% CI 31.9-55.2)
    lvc_adadt <- log(11.5); label("Apparent adADT central volume Vc/F at 51.5 kg, 50-mg-tablet reference (L)") # Table 1, adADT 'V_C/F (liters)' = 11.5 (95% CI 6.78-20.3)

    # Fraction of dADT elimination routed to adADT. Assumed, not estimated.
    fm <- fixed(0.65); label("Fraction of dADT elimination converted to adADT (unitless)") # Materials and Methods: renal clearance of dADT fixed at 35%, 'the remaining 65% was assumed to be completely metabolized into adADT'; File S1 $PK K20 = CP*0.35/V2, K23 = CP*0.65/V2

    # Allometric exponents, fixed at the canonical values. Present in the final
    # control stream but absent from Table 1 because they were not estimated.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on both clearances (unitless)") # Materials and Methods, 'exponents of 0.75 and 1 for clearance and volume parameters'; File S1 $PK (WEIGHT/51.50)**0.75
    e_wt_vc <- fixed(1); label("Allometric exponent on both central volumes (unitless)") # Materials and Methods, as above; File S1 $PK (WEIGHT/51.50)**1.00

    # Age on clearance, linear and centered on 52 years. Table 1 reports the
    # effect per 10 years older, so the per-year slope is one tenth of it.
    e_age_cl <- -0.0127; label("Linear age effect on dADT CL/F, per year older (fraction)") # Table 1, 'Age on dADT CL/F, per 10 yr older' = -12.7% (95% CI -19.9 to -8.08) => -0.127 per decade = -0.0127 per year
    e_age_cl_adadt <- -0.0212; label("Linear age effect on adADT CL/F, per year older (fraction)") # Table 1, 'Age on adADT CL/F, per 10 yr older' = -21.2% (95% CI -42.1 to -3.72) => -0.212 per decade = -0.0212 per year

    # 200-mg-tablet formulation effects, multiplicative on the 50-mg reference.
    e_tab200_mtt <- 0.401; label("200-mg-tablet effect on mean transit time (fraction longer)") # Table 1, 'Formulation on mean transit time' = 40.1% (95% CI 2.96-96.9)
    e_tab200_vc <- 1.13; label("200-mg-tablet effect on dADT central volume (fraction larger)") # Table 1, 'Formulation on dADT V_C/F' = 113% (95% CI 53.0-196)
    e_tab200_vc_adadt <- 3.64; label("200-mg-tablet effect on adADT central volume (fraction larger)") # Table 1, 'Formulation on adADT V_C/F' = 364% (95% CI 141-499)

    # Interindividual variability. Table 1 footnote d gives the printed % CV as
    # sqrt(exp(omega^2) - 1) * 100, which inverts to omega^2 = log(1 + CV^2).
    # The 'Correlations (% CV)' column is misleadingly headed: Table 1 footnote
    # f defines it as the covariance divided by sqrt(omega^2_i * omega^2_j),
    # i.e. a correlation coefficient as a percentage, so each covariance is
    # r * sd_i * sd_j. The $OMEGA BLOCK(4) of File S1 independently reproduces
    # the same six signs and comparable magnitudes.
    etalcl + etalvc + etalcl_adadt + etalvc_adadt ~ c(
      0.059687,
      0.188064, 0.703147,
      -0.132722, -0.335470, 0.753113,
      0.037652, 0.440412, 0.299168, 1.028047
    )
    etalmtt ~ 0.578424 # Table 1, IIV 'Mean transit time' = 88.5% CV; log(1 + 0.885^2) = 0.578424

    # Residual unexplained variability. Vanobberghen 2016 modeled a separate
    # additive error on natural-log-transformed concentrations for each
    # metabolite, 'which is essentially equivalent to an exponential residual
    # error on the arithmetic scale' (Materials and Methods) -- that is
    # lnorm() in nlmixr2. Table 1 prints sigma as % CV under the same footnote
    # d formula as the IIV rows, so the log-scale SD is sqrt(log(1 + CV^2)).
    expSd <- 0.923332; label("Log-scale residual SD for dADT (unitless)") # Table 1, dADT 'sigma (% CV)' = 116 (95% CI 94.1-141); sqrt(log(1 + 1.16^2)) = 0.923332
    expSd_adadt <- 0.585148; label("Log-scale residual SD for adADT (unitless)") # Table 1, adADT 'sigma (% CV)' = 63.9 (95% CI 50.1-77.3); sqrt(log(1 + 0.639^2)) = 0.585148
  })

  model({
    # Multiplicative covariate factors. Age is linear and centered on 52 years;
    # the formulation indicator is 0 for the 50-mg-tablet reference.
    cl_age <- 1 + e_age_cl * (AGE - 52)
    cl_age_adadt <- 1 + e_age_cl_adadt * (AGE - 52)
    mtt_form <- 1 + e_tab200_mtt * FORM_TRI_TAB200
    vc_form <- 1 + e_tab200_vc * FORM_TRI_TAB200
    vc_form_adadt <- 1 + e_tab200_vc_adadt * FORM_TRI_TAB200

    # Individual parameters. Allometric scaling is on disposition only; mean
    # transit time carries no weight term in the control stream.
    mtt <- exp(lmtt + etalmtt) * mtt_form
    ktr <- 7 / mtt # n = 6 transit compartments fixed, so KTR = (n + 1) / MTT

    cl <- exp(lcl + etalcl) * (WT / 51.5)^e_wt_cl * cl_age
    vc <- exp(lvc + etalvc) * (WT / 51.5)^e_wt_vc * vc_form
    cl_adadt <- exp(lcl_adadt + etalcl_adadt) * (WT / 51.5)^e_wt_cl * cl_age_adadt
    vc_adadt <- exp(lvc_adadt + etalvc_adadt) * (WT / 51.5)^e_wt_vc * vc_form_adadt

    kel <- cl / vc
    kel_adadt <- cl_adadt / vc_adadt

    # ODE system. Amounts are in nmol and volumes in L, so concentrations come
    # out in nmol/L -- the molar scale the model was fitted on. Tribendimidine
    # degrades to dADT without enzymatic involvement, and the dose compartment
    # is never observed, so the absorption chain is written in dADT-equivalent
    # moles and no molecular-weight factor appears anywhere in the system.
    #
    # Absorption: dose compartment plus six transit compartments in series,
    # every transfer at rate ktr (File S1: K14 = K45 = ... = K92 = KTR).
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(transit4) <- ktr * transit3 - ktr * transit4
    d/dt(transit5) <- ktr * transit4 - ktr * transit5
    d/dt(transit6) <- ktr * transit5 - ktr * transit6

    # dADT: one-compartment disposition. Total elimination is kel; the split
    # into a renal share (1 - fm) and a metabolic share (fm) does not change
    # the dADT profile, only how much mass reaches adADT.
    d/dt(central) <- ktr * transit6 - kel * central

    # adADT: one-compartment disposition formed from dADT at 1:1 stoichiometry
    # on the molar scale. Formation-rate-limited, so the terminal slope of
    # adADT tracks the dADT half-life (Vanobberghen 2016 Results, 'Secondary
    # PK parameters and outcomes').
    d/dt(central_adadt) <- fm * kel * central - kel_adadt * central_adadt

    Cc <- central / vc
    Cc_adadt <- central_adadt / vc_adadt

    Cc ~ lnorm(expSd)
    Cc_adadt ~ lnorm(expSd_adadt)
  })
}
