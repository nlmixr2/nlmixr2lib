Prytula_2016_tacrolimus <- function() {
  description <- "Two-compartment population PK model with first-order absorption and a fixed absorption lag time for twice-daily oral tacrolimus (Prograft) in stable paediatric renal transplant recipients at least one year after kidney transplantation (Prytula 2016). All apparent-PK parameters (CL/F, Q/F, V1/F, V2/F, ka) scale allometrically with body weight at fixed exponents (0.75 on CL/F and Q/F, 1 on V1/F and V2/F, -0.25 on ka) referenced to a 70 kg adult; V2/F is fixed at 1090 L/70 kg during covariate analysis; CL/F additionally varies with CYP3A5*1 carrier status (1+0.45-fold higher in carriers vs *3/*3 nonexpressers), gamma-glutamyltransferase (power -0.21, centred at 13 U/L), and haematocrit (power -0.59, centred at 0.34); eta_Q is perfectly correlated with eta_CL and is constructed as iiv_q_scale * etalcl (iiv_q_scale = 2.0; the 'IIV-CL-Q' parameter in Table 2); inter-individual variability is a 3x3 correlated block on (ka, CL/F, V1/F); proportional residual error."
  reference <- "Prytula AA, Cransberg K, Bouts AHM, van Schaik RHN, de Jong H, de Wildt SN, Mathot RAA. The Effect of Weight and CYP3A5 Genotype on the Population Pharmacokinetics of Tacrolimus in Stable Paediatric Renal Transplant Recipients. Clin Pharmacokinet. 2016;55(9):1129-1143. doi:10.1007/s40262-016-0390-7"
  vignette <- "Prytula_2016_tacrolimus"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline in Prytula 2016 model-building dataset (Table 1). Allometric power scaling with reference weight 70 kg and theory-based fixed exponents: 0.75 on CL/F and Q/F, 1 on V1/F and V2/F, and -0.25 on ka (Prytula 2016 Section 3.2.1). Study median 38.6 kg, range 15-86 kg.",
      source_name        = "weight"
    ),
    GGT = list(
      description        = "Gamma-glutamyltransferase activity",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Enters CL/F via the power form (GGT / 13)^e_ggt_cl with reference 13 U/L (study median). A gGT increase from 13 to 30 U/L was reported to be associated with a 16% decrease in CL/F (Prytula 2016 Section 3.2.2 / Fig. 1d). Study median 13 U/L, range 4-118 U/L (Table 1).",
      source_name        = "gGT"
    ),
    HCT = list(
      description        = "Haematocrit, expressed as a fraction of total blood volume (0-1)",
      units              = "fraction",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Enters CL/F via the power form (HCT / 0.34)^e_hct_cl with reference 0.34 (study median). A haematocrit increase from 0.34 to 0.40 was reported to be associated with a 10% decrease in CL/F (Prytula 2016 Section 3.2.2 / Fig. 1c). Prytula 2016 reports HCT as a fraction (0-1), not as percent (0-100); the canonical-register HCT entry's units (%) are explicitly overridden here so the centring value 0.34 and the power exponent -0.59 reproduce the paper's equation directly. To use a dataset that records HCT in percent, multiply the column by 0.01 before passing it to this model. Study median 0.34, range 0.21-0.44 (Table 1).",
      source_name        = "Ht"
    ),
    CYP3A5_EXPR = list(
      description        = "CYP3A5 expresser indicator: 1 if the patient carries at least one functional CYP3A5*1 allele (genotype *1/*1 or *1/*3), 0 if homozygous *3/*3.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A5 *3/*3 nonexpresser)",
      notes              = "Time-fixed germline genotype determined from rs776746 (CYP3A5 6986A>G). In the Prytula 2016 model-building cohort, DNA was available for 49 of 54 patients; among the 49 genotyped subjects the distribution was *1/*1 = 1 (2%), *1/*3 = 12 (24.5%), *3/*3 = 36 (73.5%) (Table 1). The source paper estimated a separate CL/F shift for the 5 ungenotyped subjects (h(CYP3A5 missing) = 0.41 with 61% SE, bootstrap 0.52 95% CI -0.08 to 1.52; Table 2), but that arm is not encoded here -- the wide CI makes the estimate statistically indistinguishable from both the carrier effect and from zero, and downstream simulation use cases have known genotype data. Users with subjects of unknown genotype should impute the most prevalent reference category (CYP3A5_EXPR = 0) and inspect the impact qualitatively; see vignette Assumptions and deviations. CL/F enters as (1 + 0.45 * CYP3A5_EXPR), so expressers have 45% higher apparent oral clearance than nonexpressers (Prytula 2016 Section 3.2.2 / Table 2 / equations on p. 1136).",
      source_name        = "CYP3A5"
    )
  )

  population <- list(
    species               = "human",
    n_subjects            = 54L,
    n_studies             = 2L,
    n_profiles            = 120L,
    age_range             = "3.8-18.4 years",
    age_median            = "11.1 years",
    weight_range          = "15-86 kg",
    weight_median         = "38.6 kg",
    sex_female_pct        = 48,
    race_ethnicity        = c(Caucasian = 70, Other = 30),
    disease_state         = "Stable paediatric renal transplant recipients at least one year after kidney transplantation (median 16.2 months, range 11.4-124 months posttransplant). All patients received twice-daily oral tacrolimus (Prograft) for at least 6 weeks at baseline and the dose had not been amended for at least 14 days. Initial immunosuppressive therapy included basiliximab (n = 37) plus corticosteroids, a calcineurin inhibitor, and an antimetabolite. Before 2009 most children started on cyclosporine and switched to tacrolimus after 3-6 months; an early-corticosteroid-withdrawal protocol was used in Rotterdam from July 2009 onward.",
    dose_range            = "Twice-daily oral tacrolimus titrated by therapeutic drug monitoring to a posttransplant year >=1 target trough (C0) of 4-8 ug/L; the cohort median daily dose was 0.12-0.15 mg/kg/day depending on CYP3A5 genotype (Prytula 2016 Section 3.1).",
    regions               = "Netherlands (Erasmus MC-Sophia Children's Hospital, Rotterdam, n = 45; Emma Children's Hospital, Amsterdam, n = 9).",
    cyp3a5_distribution   = "Of 49 genotyped patients (DNA unavailable for 5): *1/*1 n = 1 (2%), *1/*3 n = 12 (24.5%), *3/*3 n = 36 (73.5%); Hardy-Weinberg equilibrium confirmed.",
    abcb1_distribution    = "ABCB1 3435 C>T genotype (Table 1): T/C n = 32 (65.3%), C/C n = 9 (18.4%), T/T n = 8 (16.3%). ABCB1 polymorphism was tested as a covariate on CL/F but was not retained in the final model.",
    sampling_window       = "Abbreviated 4-h profiles (Rotterdam): predose (C0) and 10, 30, 90, 120, 240 min postdose. Abbreviated 2-h profiles (Amsterdam): C0 and C120 only. AUC12 was estimated by Bayesian analysis from the abbreviated profiles using the method previously validated for paediatric tacrolimus.",
    assay                 = "Whole-blood tacrolimus measured by LC-MS/MS (83% of samples; LLOQ 0.2 ug/L) or MEIA (remainder; LLOQ 1.0 ug/L). The assay-method covariate was tested in the model but not retained.",
    notes                 = "Retrospective cohort study; data collected April 1993 - June 2011 across two centres. The number of profiles per child ranged from 1 to 5 (median 2); 20 children had only one profile. External validation was performed on an independent cohort of 27 children (Section 2.1, Table 1 right columns). The simulated trough scenarios (Section 3.2.4 / Fig. 4) cover doses 0.025-0.2 mg/kg twice daily across body weights 10-90 kg and CYP3A5 strata, with target trough 4-8 ug/L."
  )

  ini({
    # Final-model fixed-effect estimates from Prytula 2016 Table 2 (p. 1136).
    # Reference subject: 70 kg paediatric kidney recipient, CYP3A5 *3/*3
    # nonexpresser, gGT 13 U/L, haematocrit 0.34. All apparent clearances
    # (CL/F, Q/F) are in L/h; apparent volumes (V1/F, V2/F) in L; ka in 1/h;
    # tlag in h.
    lka   <- log(0.23) ; label("Absorption rate constant ka at WT = 70 kg (1/h)")                                                  # Table 2 final ka = 0.23 1/h (SE 15%)
    ltlag <- log(0.47) ; label("Absorption lag time (h)")                                                                          # Table 2 final tlag = 0.47 h (SE 6%)
    lcl   <- log(35)   ; label("Apparent oral clearance CL/F at WT = 70 kg, CYP3A5 *3/*3, gGT 13 U/L, HCT 0.34 (L/h)")              # Table 2 final CL/F (CYP3A5 *3/*3) = 35 L/h/70 kg (SE 7%)
    lvc   <- log(12)   ; label("Apparent central volume V1/F at WT = 70 kg (L)")                                                   # Table 2 final V1/F = 12 L/70 kg (SE 12%)
    lq    <- log(68)   ; label("Apparent inter-compartmental clearance Q/F at WT = 70 kg (L/h)")                                   # Table 2 final Q/F = 68 L/h/70 kg (SE 17%)
    lvp   <- fixed(log(1090)) ; label("Apparent peripheral volume V2/F at WT = 70 kg (L; fixed during covariate analysis)")        # Table 2 V2/F = 1090 L/70 kg (FIX)

    # Covariate effects on CL/F. The Table 2 equations are:
    #   CYP3A5 = *3/*3: CL/F = 35 * (WT/70)^0.75 * (gGT/13)^-0.21 * (HCT/0.34)^-0.59 * exp(eta_CL + kappa_CL)
    #   CYP3A5 = *1 carrier: CL/F = 35 * (1 + 0.45) * (WT/70)^0.75 * (gGT/13)^-0.21 * (HCT/0.34)^-0.59 * exp(eta_CL + kappa_CL)
    # The CYP3A5*1 carrier effect is encoded as an additive shift on the
    # multiplier so the equation reads identically to the source paper for
    # CYP3A5_EXPR in {0, 1}.
    e_cyp3a5_expr_cl <- 0.45  ; label("CYP3A5*1-carrier additive-shift coefficient on CL/F multiplier (1 + e_cyp3a5_expr_cl * CYP3A5_EXPR); expressers have 45% higher CL/F")  # Table 2 theta(CYP3A5 *1/*1, *1/*3) = 0.45 (SE 38%)
    e_ggt_cl         <- -0.21 ; label("gGT power exponent on CL/F, centred at 13 U/L")                                                                                          # Table 2 theta_gGT = -0.21 (SE 31%)
    e_hct_cl         <- -0.59 ; label("Haematocrit power exponent on CL/F, centred at 0.34 fraction")                                                                           # Table 2 theta_Ht = -0.59 (SE 49%)

    # Allometric exponents (Prytula 2016 Section 3.2.1: "Allometric scaling
    # of the parameters with fixed exponents (values 0.75 [CL/F, Q/F], 1
    # [V1/F, V2/F] and -0.25 [ka])"). All five exponents are theory-based
    # and held fixed during estimation.
    e_wt_cl <- fixed(0.75)  ; label("Allometric exponent of (WT/70) on CL/F (unitless; fixed at theory value)")                    # Section 3.2.1
    e_wt_q  <- fixed(0.75)  ; label("Allometric exponent of (WT/70) on Q/F (unitless; fixed at theory value)")                     # Section 3.2.1
    e_wt_vc <- fixed(1)     ; label("Allometric exponent of (WT/70) on V1/F (unitless; fixed at theory value)")                    # Section 3.2.1
    e_wt_vp <- fixed(1)     ; label("Allometric exponent of (WT/70) on V2/F (unitless; fixed at theory value)")                    # Section 3.2.1
    e_wt_ka <- fixed(-0.25) ; label("Allometric exponent of (WT/70) on ka (unitless; fixed at theory value)")                      # Section 3.2.1

    # Q/F inter-individual variability is forced to perfect correlation with
    # CL/F's eta and is constructed as `q_eta_scale * etalcl` (Prytula 2016
    # Section 3.2.1: "The correlation between the inter-patient variability
    # of CL/F and Q/F was high and fixed to 1"). The single estimated
    # scaling parameter is "IIV-CL-Q" in Table 2. Stored as a structural
    # theta rather than as a separate `etalq` slot, because for a
    # correlation of exactly 1 the two forms are mathematically equivalent
    # and the deterministic form avoids the singular-covariance numerical
    # issue that a 1-correlation block triggers.
    q_eta_scale <- 2.0 ; label("Scaling factor relating eta_Q to eta_CL (correlation fixed to 1; eta_Q = q_eta_scale * etalcl)")    # Table 2 IIV-CL-Q = 2.0 (SE 17%)

    # Correlated 3 x 3 inter-individual variability block on (ka, CL/F, V1/F).
    # Table 2 reports IIV as %CV and the table footnote reports correlations
    # r(ka, CL) = 0.11, r(ka, V1) = 0.19, r(CL, V1) = 0.38 (rho_{ka.CL},
    # rho_{ka.V1}, rho_{CL.V1}). Convert %CV to log-scale variance via
    # omega^2 = log(1 + CV^2):
    #   ka   CV 58%  -> log(1 + 0.58^2)  = 0.28998
    #   CL   CV 45%  -> log(1 + 0.45^2)  = 0.18439
    #   V1   CV 170% -> log(1 + 1.70^2)  = 1.35841
    # Off-diagonals = r * sqrt(omega^2_i) * sqrt(omega^2_j):
    #   cov(ka, CL) = 0.11 * sqrt(0.28998) * sqrt(0.18439) = 0.02543
    #   cov(ka, V1) = 0.19 * sqrt(0.28998) * sqrt(1.35841) = 0.11925
    #   cov(CL, V1) = 0.38 * sqrt(0.18439) * sqrt(1.35841) = 0.19019
    etalka + etalcl + etalvc ~ c(0.28998,
                                  0.02543, 0.18439,
                                  0.11925, 0.19019, 1.35841)                                                                       # Table 2 IIV (ka, CL/F, V1/F) %CV + footnote correlations

    # Proportional residual error on the linear concentration scale (Prytula
    # 2016 Section 3.2.1: "The residual error was described with a
    # proportional error model").
    propSd <- 0.21 ; label("Proportional residual error (fraction)")                                                               # Table 2 final proportional RUV = 21% (SE 6%)
  })

  model({
    # Allometric scaling references a 70 kg adult.
    wt70 <- WT / 70

    # Individual PK parameters with Prytula 2016 covariate equations. eta_Q
    # is fully determined by etalcl via the perfect-correlation constraint
    # (no separate etalq is sampled).
    cl <- exp(lcl + etalcl) *
          (1 + e_cyp3a5_expr_cl * CYP3A5_EXPR) *
          wt70 ^ e_wt_cl *
          (GGT / 13) ^ e_ggt_cl *
          (HCT / 0.34) ^ e_hct_cl
    q    <- exp(lq + q_eta_scale * etalcl) * wt70 ^ e_wt_q
    vc   <- exp(lvc + etalvc) * wt70 ^ e_wt_vc
    vp   <- exp(lvp) * wt70 ^ e_wt_vp
    ka   <- exp(lka + etalka) * wt70 ^ e_wt_ka
    tlag <- exp(ltlag)

    # Two-compartment oral disposition. Dose lands in `depot`; bioavailability
    # is implicit in the apparent CL/F and V/F parameterisation.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    # Tacrolimus is reported in whole blood in ug/L (= ng/mL). Convert internal
    # mg/L (dose in mg, vc in L) to ng/mL by multiplying by 1000.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
