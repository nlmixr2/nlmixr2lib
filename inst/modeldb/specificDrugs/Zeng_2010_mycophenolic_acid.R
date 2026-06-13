Zeng_2010_mycophenolic_acid <- function() {
  description <- "Two-compartment population PK model with first-order absorption for mycophenolic acid (MPA, the active moiety of mycophenolate mofetil MMF) in children and young people undergoing blood or marrow, kidney, and liver transplantation (Zeng 2010, Br J Clin Pharmacol). Both intravenous (2 h infusion to the central compartment) and oral (depot with first-order absorption rate ka and bioavailability F) routes were modelled jointly. Apparent clearance CL/F (combined IV/oral typical value) varies linearly with body weight via (1 + theta_WT * WT/27.9) where 27.9 kg is the cohort median, and additively with concomitant calcineurin-inhibitor type via (1 + theta_CYTA * CYTA): CYTA = 0 on ciclosporin (the reference; n = 23) and CYTA = 1 on tacrolimus (n = 15) (Zeng 2010 Table 2 model 4). Ciclosporin inhibits MRP2-mediated biliary efflux of MPAG and thus suppresses enterohepatic recirculation of MPA, so paediatric MPA CL/F is approximately 2.5x higher under ciclosporin (CYTA = 0; multiplier 1.0) than under tacrolimus (CYTA = 1; multiplier 1 + (-0.60) = 0.40). Inter-individual variability is diagonal on CL, ka, and F; residual error is exponential (log-normal) with SD 0.48 on the log scale. Inter-occasion variability of 5.8% CV on CL/F reported in Zeng 2010 Table 3 (with occasion defined as 7 days in patients on daily MMF) is NOT encoded structurally here because the model-library use case does not define an operational occasion column; downstream users who want to simulate IOV can add an OCC indicator and a per-occasion eta in rxode2."
  reference <- paste(
    "Zeng L, Blair EYL, Nath CE, Shaw PJ, Earl JW, Stephen K, Montgomery K,",
    "Coakley JC, Hodson E, Stormon M, McLachlan AJ.",
    "Population pharmacokinetics of mycophenolic acid in children and young",
    "people undergoing blood or marrow and solid organ transplantation.",
    "Br J Clin Pharmacol. 2010;70(4):567-579.",
    "doi:10.1111/j.1365-2125.2010.03734.x.",
    sep = " "
  )
  vignette <- "Zeng_2010_mycophenolic_acid"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline per the source paper's analysis. Enters CL/F as the linear scaling factor (1 + theta_WT * WT / 27.9) where 27.9 kg is the cohort median (Zeng 2010 Table 2 model 3 footnote: 'The population CL term was standardized to 27.9 kg, which represents the median value of weight in this study group.'). Cohort range 3.4-87.7 kg, median 27.9 kg (Table 1). Zeng 2010 Section 'Covariate analysis' reports that linear weight scaling outperformed allometric scaling with a fixed 0.75 exponent (delta-OFV -17.40 vs -4.92 and IIV reduction 54.3 -> 36.9 vs 54.3 -> 46.8) so the linear form is what the final model uses.",
      source_name        = "WT"
    ),
    CONMED_CSA = list(
      description        = "Concomitant ciclosporin (CsA) indicator: 1 if the patient was co-administered ciclosporin as the calcineurin inhibitor (CNI), 0 if the patient was on tacrolimus (the alternative CNI in the source cohort).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (concomitant tacrolimus; lower MPA CL/F)",
      notes              = "Time-fixed per patient in the Zeng 2010 cohort because subjects were on one CNI regimen across the sampling window. Source column CYTA in Zeng 2010 carries the INVERTED value convention: CYTA = 0 if patient is on ciclosporin, CYTA = 1 if patient is without ciclosporin (i.e. on tacrolimus); the canonical CONMED_CSA flips this so that 1 = ciclosporin, matching the deWinter 2009 mycophenolic-acid precedent and the broader CONMED_* register convention (1 = exposed to the named conmed). The source paper's covariate coefficient e_csa_cl = theta_8 = -0.60 acts on the absence-of-ciclosporin indicator, so the model() block applies it via (1 - CONMED_CSA) to recover the paper's CYTA encoding. Ciclosporin inhibits MRP2-mediated biliary efflux of MPAG and suppresses enterohepatic recirculation of MPA, producing approximately 2.5x higher MPA CL/F under ciclosporin than under tacrolimus in this paediatric cohort (Zeng 2010 Discussion: 'co-administered ciclosporin as opposed to tacrolimus resulted in a mean increase in MPA CL of 63%'; the 2.5x ratio reflects both the direct CNI effect on EHC and the WT * CYTA interaction at the median weight). Cohort distribution: 23 patients on ciclosporin, 15 on tacrolimus (Table 1).",
      source_name        = "CYTA (inverted: CYTA = 1 - CONMED_CSA)"
    )
  )

  population <- list(
    species            = "human",
    n_subjects         = 38L,
    n_studies          = 1L,
    n_observations     = 859L,
    age_range          = "0.4-19.9 years",
    age_median         = "8.4 years (Table 1)",
    weight_range       = "3.4-87.7 kg",
    weight_median      = "27.9 kg (Table 1)",
    sex_female_pct     = NA_real_,
    race_ethnicity     = "Not reported in source paper.",
    disease_state      = "Children and young people receiving mycophenolate mofetil (MMF) after blood or marrow (n = 23), kidney (n = 5), or liver (n = 10) transplantation at The Children's Hospital at Westmead, Sydney. All subjects received MMF as part of an immunosuppressive regimen alongside either ciclosporin (n = 23) or tacrolimus (n = 15) as the calcineurin inhibitor; 9 also received concomitant acyclovir.",
    dose_range         = "MMF 10-15 mg/kg IV (2 h infusion) or oral, twice or three times daily.",
    regions            = "Australia (single centre: The Children's Hospital at Westmead, Sydney, NSW).",
    transplant_mix     = "Blood or marrow transplant: 23; kidney transplant: 5; liver transplant: 10 (Table 1). 13 patients had IV dosing only, 18 oral only, 7 received both routes.",
    sampling_design    = "Age-specific sampling protocol over 8-12 h: weighed > 20 kg children had 13-14 intensive samples per dose interval; weighed <= 20 kg children had a sparse 5-6 sample design. The remaining cohort had randomly timed but accurately recorded blood samples. All samples were collected at steady state.",
    cni_distribution   = "Concomitant ciclosporin: 23 patients (60.5%); concomitant tacrolimus: 15 patients (39.5%).",
    iov_structure      = "Inter-occasion variability (IOV) on CL/F (5.8% CV) was identified in addition to the diagonal IIV (Zeng 2010 Table 3 final-model). This model file does NOT encode IOV structurally -- the source paper does not define an operational occasion column for the model-library use case ('each occasion was defined as 7 days in patients who were administered MMF daily'). Downstream users who want to simulate IOV can add an OCC indicator in their event dataset and a per-occasion eta in rxode2; see vignette Assumptions and deviations.",
    notes              = "Prospective single-centre observational study. Total MPA was measured by HPLC (LLOQ 0.07 mg/L). 859 MPA concentrations from 38 subjects, median 23 samples/patient (range 2-49). 13 patients had intensive sampling and 25 had sparse or random sampling. Enterohepatic recirculation was observed in 9 of the 13 intensively sampled patients but could not be modelled separately given the limited number of such patients; this likely inflates the residual error in the final model (Zeng 2010 Discussion). Patient characteristics from Table 1."
  )

  ini({
    # Final-model parameter estimates from Zeng 2010 Table 3 ('Population
    # pharmacokinetic parameter estimates derived from the base and the final
    # models'). All apparent clearances (CL/F, Q/F) are in L/h; apparent
    # volumes (V1/F, V2/F) in L; ka in 1/h; F as a dimensionless fraction.
    # Reference subject for the typical-value structural parameters: WT = 0 kg
    # (the linear-WT scaling form (1 + theta_WT * WT/27.9) means theta_1 is
    # the extrapolated CL/F intercept, NOT the CL/F at the cohort median; at
    # WT = 27.9 kg on ciclosporin the typical CL/F is 6.42 * (1 + 1.09) = 13.42
    # L/h matching the Discussion).
    lka     <- log(0.39);    label("First-order absorption rate constant ka (1/h)")                                            # Zeng 2010 Table 3 final ka = 0.39 1/h
    lcl     <- log(6.42);    label("Apparent oral / IV clearance CL/F intercept (L/h) for the linear-WT scaling form")         # Zeng 2010 Table 3 final theta_1 (CL) = 6.42 L/h
    lvc     <- log(7.24);    label("Apparent central volume of distribution V1/F (L)")                                         # Zeng 2010 Table 3 final V1 = 7.24 L
    lq      <- log(3.74);    label("Apparent inter-compartmental clearance Q/F (L/h)")                                         # Zeng 2010 Table 3 final Q = 3.74 L/h
    lvp     <- log(16.80);   label("Apparent peripheral volume of distribution V2/F (L)")                                      # Zeng 2010 Table 3 final V2 = 16.80 L
    lfdepot <- log(0.48);    label("Oral bioavailability F (dimensionless fraction)")                                          # Zeng 2010 Table 3 final F = 0.48

    # Covariate effects (Zeng 2010 Table 2 final covariate model 4):
    #   TVCL = theta_1 * (1 + theta_7 * WT/27.9) * (1 + theta_8 * CYTA)
    # where CYTA = 0 if on ciclosporin and CYTA = 1 if on tacrolimus.
    # Canonical CONMED_CSA inverts this convention (CONMED_CSA = 1 if on
    # ciclosporin) so the model() block applies e_csa_cl via the absence
    # indicator (1 - CONMED_CSA) to recover the paper's CYTA encoding.
    e_wt_cl  <- 1.09;        label("Linear scaling slope of CL/F on WT/27.9 (unitless)")                                       # Zeng 2010 Table 3 final theta_7 (Weight factor) = 1.09
    e_csa_cl <- -0.60;       label("Additive effect of (1 - CONMED_CSA) on CL/F (unitless; -0.60 means 60% lower CL on tacrolimus relative to ciclosporin)") # Zeng 2010 Table 3 final theta_8 (Concomitant ciclosporin factor) = -0.60

    # Inter-individual variability. Zeng 2010 used an exponential IIV model
    # (theta_i = theta * exp(eta_i)) and reported per-parameter omega values
    # as %CV. CV% relates to log-normal variance via omega^2 = log(1 + CV^2).
    # IIV was only identifiable on CL, ka, and F; the source paper attempted
    # but could not reliably estimate IIV on V1, V2, Q (Zeng 2010 Section
    # 'Base model building').
    etalcl     ~ 0.0953     # 0.316^2; Zeng 2010 Table 3 final omega_CL  = 31.6% CV -> log(1 + 0.316^2)
    etalka     ~ 0.2994     # 0.591^2; Zeng 2010 Table 3 final omega_ka  = 59.1% CV -> log(1 + 0.591^2)
    etalfdepot ~ 0.1131     # 0.346^2; Zeng 2010 Table 3 final omega_F   = 34.6% CV -> log(1 + 0.346^2)

    # Residual variability. Zeng 2010 used an exponential (log-normal)
    # residual error model: Y = Y_pred * exp(eps_1) with eps_1 ~ N(0, sigma^2),
    # which is the same as nlmixr2's ~ lnorm(expSd). Table 3 reports
    # sigma_1 = 0.48 directly on the log scale.
    expSd      <- 0.48;      label("Log-normal residual error SD (on log scale)")                                              # Zeng 2010 Table 3 final sigma_1 = 0.48 (exponential error model)
  })

  model({
    # Reference / centering values from Zeng 2010 Table 2 model 4.
    ref_wt <- 27.9

    # Individual parameters. The source paper's covariate model uses a linear
    # weight scaling (1 + theta_WT * WT / ref_wt) on CL/F and an additive
    # calcineurin-inhibitor effect (1 + theta_CYTA * CYTA) where CYTA = 0 on
    # ciclosporin (the reference) and CYTA = 1 on tacrolimus. Canonical
    # CONMED_CSA = 1 means ciclosporin, so (1 - CONMED_CSA) reproduces the
    # source paper's CYTA indicator without further transformation.
    ka      <- exp(lka + etalka)
    cl      <- exp(lcl + etalcl) * (1 + e_wt_cl * (WT / ref_wt)) *
               (1 + e_csa_cl * (1 - CONMED_CSA))
    vc      <- exp(lvc)
    q       <- exp(lq)
    vp      <- exp(lvp)
    fdepot  <- exp(lfdepot + etalfdepot)

    # Micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system. State variables and units:
    #   depot       (mg)  -- oral MMF dose, expressed as MPA-equivalent mass
    #   central     (mg)  -- MPA central
    #   peripheral1 (mg)  -- MPA peripheral
    # IV doses bypass depot and load central directly via cmt = 2 (or the
    # rxode2 'central' name); oral doses target depot with bioavailability
    # fdepot applied below.
    d/dt(depot)        <- -ka * depot
    d/dt(central)      <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1)  <-  k12 * central - k21 * peripheral1

    # Oral bioavailability (applies only to depot-targeted doses; IV doses
    # delivered to central are unaffected).
    f(depot) <- fdepot

    # Observation: total MPA plasma concentration in mg/L.
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
