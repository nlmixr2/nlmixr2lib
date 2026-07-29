Nicholas_2009_ponezumab <- function() {
  description <- "Two-compartment intravenous population PK model for PF-04360365 (ponezumab), a humanized anti-amyloid IgG2 delta-a monoclonal antibody, in adults with mild-to-moderate Alzheimer's disease; allometric body-weight scaling with estimated exponents on every disposition parameter and a full 4x4 inter-individual block on (CL, V1, V2, Q) (Nicholas 2009 preliminary popPK)"
  reference <- "Nicholas T, Knebel W, Gastonguay MR, Bednar MM, Billing B, Landen JW, Kupiec JW, Corrigan B, Laurencot R, Zhao Q. Preliminary population pharmacokinetic modeling of PF-04360365, a humanized anti-amyloid monoclonal antibody, in patients with mild-to-moderate alzheimer's disease. Alzheimers Dement. 2009;5(4 Suppl):P425. doi:10.1016/j.jalz.2009.04.270"
  vignette <- "Nicholas_2009_ponezumab"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline body weight. Allometric power scaling on CL, V1, Q, and V2 with reference weight 70 kg (Nicholas 2009 Methods: 'Allometric scaling was implemented using a reference weight of 70 kg.'). The allometric exponents (Theta_5, Theta_6, Theta_7, Theta_8) are estimated rather than fixed at 0.75 / 1; CL and V1 are precisely estimated while Q and V2 exponents are imprecise (Theta_7 95% CI crosses zero, Theta_8 95% CI is wide). Covariate analysis in this preliminary model was limited to body weight because of the small sample size (n=26) and narrow age range (60-80 years).",
      source_name        = "WT"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 26L,
    n_studies        = 1L,
    age_range        = "60-80 years",
    weight_range     = "Not tabulated in the poster.",
    sex_female_pct   = NA_real_,
    race_ethnicity   = NA_character_,
    disease_state    = "Mild-to-moderate Alzheimer's disease (Mini Mental State Examination score 16-26).",
    dose_range       = "Single intravenous dose of PF-04360365 (ponezumab) 0.1-10 mg/kg dose-escalation (n=26 active; n=11 placebo not included in the analysis).",
    regions          = "Not tabulated in the poster.",
    notes            = "Plasma concentrations were measured by ELISA over an analytical range of 156-10,000 ng/mL with intra- and inter-assay precisions within 10% and accuracy within 16.0%. BLQ and missing concentrations were excluded; subjects with no concentration data were not included. Modeling was performed in NONMEM VI Level 2.0 with FOCE-I using ADVAN3 TRANS4 (Nicholas 2009 Methods)."
  )

  ini({
    # Structural disposition parameters at the reference 70 kg subject
    # (Nicholas 2009 Table 1, fixed-effect parameters with % SE).
    lcl <- log(0.00684); label("Clearance at WT=70 kg (CL, L/h)")                             # Nicholas 2009 Table 1 Theta_1: 0.00684 L/h (6% SE; 95% CI 0.00597-0.00771)
    lvc <- log(3.16);    label("Central volume of distribution at WT=70 kg (V1, L)")          # Nicholas 2009 Table 1 Theta_2: 3.16 L (3% SE; 95% CI 2.95-3.37)
    lq  <- log(0.0210);  label("Intercompartmental clearance at WT=70 kg (Q, L/h)")           # Nicholas 2009 Table 1 Theta_3: 0.0210 L/h (10% SE; 95% CI 0.0170-0.0250)
    lvp <- log(5.34);    label("Peripheral volume of distribution at WT=70 kg (V2, L)")       # Nicholas 2009 Table 1 Theta_4: 5.34 L (8% SE; 95% CI 4.49-6.18)

    # Allometric power exponents on body weight, estimated (not fixed)
    # (Nicholas 2009 Table 1, second row of each Theta pair).
    # Form per parameter: theta_i = theta_pop * (WT / 70) ^ exponent_i.
    e_wt_cl <- 0.911; label("Allometric exponent of WT on CL (unitless; reference 70 kg)")    # Nicholas 2009 Table 1 Theta_5: 0.911 (37% SE; 95% CI 0.370-1.45)
    e_wt_vc <- 0.573; label("Allometric exponent of WT on V1 (unitless; reference 70 kg)")    # Nicholas 2009 Table 1 Theta_6: 0.573 (34% SE; 95% CI 0.194-0.951)
    e_wt_q  <- 0.236; label("Allometric exponent of WT on Q (unitless; reference 70 kg)")     # Nicholas 2009 Table 1 Theta_7: 0.236 (126% SE; 95% CI -0.346-0.817)
    e_wt_vp <- 0.590; label("Allometric exponent of WT on V2 (unitless; reference 70 kg)")    # Nicholas 2009 Table 1 Theta_8: 0.590 (54% SE; 95% CI -0.0288-1.20)

    # Inter-individual variability - full 4x4 OMEGA block on (CL, V1, V2, Q).
    # Nicholas 2009 Methods: "Inter-individual random effects were modeled with
    # exponential variance models. Covariance was described with a full block
    # omega matrix." This is the standard log-normal interpretation, i.e.
    # P_i = P_typ * exp(eta_i) with eta ~ N(0, OMEGA), so the reported OMEGA
    # entries are variances/covariances on the log scale and translate
    # directly to the canonical etalcl + etalvc + etalvp + etalq variances.
    #
    # The OMEGA index ordering in Nicholas 2009 Table 1 is (CL, V1, V2, Q):
    # the table reports Omega_1.1 CL, Omega_1.2 CL-V1, Omega_2.2 V1,
    # Omega_1.3 CL-V2, Omega_2.3 V1-V2, Omega_3.3 V2, Omega_1.4 CL-Q,
    # Omega_2.4 V1-Q, Omega_3.4 V2-Q, Omega_4.4 Q. The canonical
    # etalcl + etalvc + etalvp + etalq ordering is (CL, V1=Vc, V2=Vp, Q),
    # which matches the table indices 1:4 exactly. The lower-triangular row
    # order below is therefore:
    #   var(CL)
    #   cov(CL,V1),  var(V1)
    #   cov(CL,V2),  cov(V1,V2),  var(V2)
    #   cov(CL,Q),   cov(V1,Q),   cov(V2,Q),   var(Q)
    etalcl + etalvc + etalvp + etalq ~ c(
      0.0714,
      0.0268, 0.0312,
      0.0756, 0.0465, 0.184,
      0.0421, 0.0424, 0.107, 0.0895)                                                   # Nicholas 2009 Table 1 inter-individual variance entries with reported % SE: var(CL)=0.0714 (44%); cov(CL,V1)=0.0268 (69%), var(V1)=0.0312 (36%); cov(CL,V2)=0.0756 (42%), cov(V1,V2)=0.0465 (57%), var(V2)=0.184 (53%); cov(CL,Q)=0.0421 (76%), cov(V1,Q)=0.0424 (54%), cov(V2,Q)=0.107 (61%), var(Q)=0.0895 (55%)

    # Residual error - proportional only (Nicholas 2009 Methods: "Additive and
    # proportional error structures were examined and a proportional error
    # model was utilized for the residual error model.").
    # Table 1 reports sigma^2 prop = 0.00998 -> SD = sqrt(0.00998) = 0.0999.
    propSd <- 0.0999; label("Proportional residual error (fraction)")                          # Nicholas 2009 Table 1 sigma^2 prop: 0.00998 (11% SE)
  })
  model({
    # Individual disposition parameters with log-normal IIV and
    # estimated allometric scaling to a 70 kg reference (Nicholas 2009 Methods
    # and Table 1).
    cl <- exp(lcl + etalcl) * (WT / 70) ^ e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70) ^ e_wt_vc
    vp <- exp(lvp + etalvp) * (WT / 70) ^ e_wt_vp
    q  <- exp(lq  + etalq ) * (WT / 70) ^ e_wt_q

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment IV disposition with linear elimination from the central
    # compartment - dose enters central directly via IV (Nicholas 2009 Methods:
    # "implemented using ADVAN 3 TRANS4"). Dose in mg, volumes in L -> central
    # concentration Cc = central / vc in mg/L = ug/mL.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
