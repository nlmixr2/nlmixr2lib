Westra_2025_osimertinib_brownbase <- function() {
  description <- "Generalizability re-analysis of the Westra 2025 cobicistat boosting study: the previously published Brown 2017 joint osimertinib / AZ5104 population PK structure, with every Brown parameter held fixed and a single new multiplicative cobicistat factor of 0.678 estimated on osimertinib apparent clearance from the OSIBOOST data, a 32.2 percent reduction. First-order oral absorption feeds an osimertinib (parent) compartment, from which a serial AZ5104 metabolite compartment is formed at a rate fixed to 25 percent of the parent elimination rate constant. Body weight enters as estimated power exponents on parent clearance, parent volume and metabolite clearance, and serum albumin as a power exponent on parent volume, all inherited fixed from Brown 2017 at reference values 62 kg and 39 g/L. Companion to Westra_2025_osimertinib.R, which is the model developed de novo on the same data; this file exists because the paper reports it as a separate analysis establishing that the cobicistat effect reproduces in an independently developed model structure."
  reference <- "Westra N, Kruithof PD, Croes S, van Geel RMJM, Hendriks LEL, Touw DJ, Kosterink JGW, Stevens J, Oude Munnink TH, Mian P. Osimertinib Cost Minimization in Non-Small Cell Lung Cancer (NSCLC) Treatment: Hypothesis Generation for a Population Pharmacokinetic Approach for Equivalent Dose Optimization of Osimertinib in Combination with Cobicistat. J Clin Pharmacol. 2025;65(12):1687-1698. doi:10.1002/jcph.70085. Base structure and fixed parameters from Brown K, Comisar C, Witjes H, Maringwa J, de Greef R, Vishwanathan K, Cantarini M, Cox E. Population pharmacokinetics and exposure-response of osimertinib in patients with non-small cell lung cancer. Br J Clin Pharmacol. 2017;83(6):1216-1226. doi:10.1111/bcp.13223"
  vignette <- "Westra_2025_osimertinib"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight (baseline; reported in kg).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effects inherited fixed from Brown 2017 with reference weight 62 kg: exponent 0.56 on parent CL/F, 0.65 on parent V/F, and 0.99 on AZ5104 CL/F. AZ5104 V/F carries no weight effect. Supplementary Part SII NONMEM $PK block; $THETA 8, 9 and 10 are all flagged FIX.",
      source_name        = "BW"
    ),
    ALB = list(
      description        = "Baseline serum albumin concentration.",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on parent V/F with exponent 1.33 and reference albumin 39 g/L, inherited fixed from Brown 2017 (supplementary Part SII $PK 'V1 = (THETA(5)*((BW/62)**THETA(9))*((ALB/39)**THETA(11)))'; $THETA 11 is flagged FIX). Albumin was separately tested and NOT retained in the de novo Westra 2025 model, so it appears here only through the inherited Brown 2017 structure.",
      source_name        = "ALB"
    ),
    CONMED_COBICISTAT = list(
      description        = "Concomitant cobicistat 150 mg once-daily coadministration indicator (1 = boosted with cobicistat, 0 = osimertinib monotherapy).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (osimertinib monotherapy; the reference state is the pre-boost osimertinib 80 mg QD phase of the OSIBOOST trial).",
      notes              = "Multiplicative power-form effect on osimertinib apparent clearance only: cl is multiplied by 0.678^CONMED_COBICISTAT, i.e. a 32.2 percent reduction in CL/F when cobicistat is coadministered (Westra 2025 Results 'Generalizability'; supplementary Part SII $THETA 12, the ONLY theta in that control stream not flagged FIX). The paper reports the addition was highly significant (P < .0001, delta OFV = -63.4) and that extending the same covariate to AZ5104 CL/F did not further improve the fit (P > .05). This 32.2 percent reduction closely reproduces the 29.6 percent found in the de novo model Westra_2025_osimertinib.R, which is the point of the analysis. Time-varying within a subject, as in the companion model.",
      source_name        = "COBI"
    )
  )

  # Issue #482: verified against the supplementary Part SII NONMEM $MODEL
  # block (compartments ABSORB, PARENT, METABOLITE) and the $PK scaling
  # statements S2 = V1/1000 and S3 = V2/1000, which place both observed
  # concentrations in ug/L.
  compartmentData <- list(
    depot          = list(analyte = "osimertinib", units = "mg", specimen = "administration site", verified = TRUE),
    central        = list(analyte = "osimertinib", units = "mg", specimen = "plasma", verified = TRUE),
    central_az5104 = list(analyte = "AZ5104", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 11,
    n_studies      = 1,
    age_median     = "69.0 years",
    weight_median  = "78.5 kg",
    height_median  = "166.0 cm",
    bmi_median     = "23.6 kg/m^2",
    bsa_median     = "1.9 m^2",
    sex_female_pct = 63.6,
    race_ethnicity = c(Caucasian = 100.0),
    disease_state  = "Advanced non-small cell lung cancer, adenocarcinoma histology in 100 percent of the cohort, ECOG/WHO performance status 0-1 in 100 percent. All patients were on established osimertinib treatment and were selected for a relatively low steady-state trough concentration (Cmin,ss at or below 195 ug/L) under osimertinib 80 mg once daily.",
    dose_range     = "Oral osimertinib 80 mg once daily (one patient received an increased dose of 160 mg once daily), first as monotherapy and then with oral cobicistat 150 mg once daily added for at least 21 days to reach steady state.",
    regions        = "Netherlands (Maastricht University Medical Centre and the Antoni van Leeuwenhoek hospital).",
    n_observations = "88 osimertinib and AZ5104 plasma concentrations in total from the 11 patients, the same dataset as the de novo model.",
    notes          = "The cobicistat factor of this model was estimated on the 11-patient OSIBOOST cohort described above (Westra 2025 Table 1). Every OTHER parameter is inherited fixed from Brown 2017, whose own estimation population was 780 subjects (748 advanced-NSCLC patients from AURA and AURA2 plus 32 healthy volunteers, weight median 62 kg, albumin median 39 g/L). Users comparing this model against Westra_2025_osimertinib.R should note that its much larger between-subject variability is inherited from that broad Brown 2017 population and is not a property of the OSIBOOST cohort; see the vignette Errata for the omega-scale reading."
  )

  ini({
    # Structural parameters, inherited fixed from Brown 2017 Table 2.
    # Reference values are typical values at 62 kg body weight and
    # 39 g/L serum albumin (Brown 2017 Table 1 medians). Every value
    # below carries an explicit FIX flag in the supplementary Part SII
    # NONMEM $THETA block, so all are wrapped in fixed().
    lka        <- fixed(log(0.24)); label("First-order oral absorption rate constant (1/h)")                       # Part SII $THETA 3 '(0.24) FIX'; Brown 2017 Table 2 ka = 0.24 1/h
    lcl        <- fixed(log(14.2)); label("Apparent osimertinib clearance, CL/F, at 62 kg without cobicistat (L/h)")  # Part SII $THETA 4 '(14.2) FIX'; Brown 2017 Table 2 CLparent/F = 14.2 L/h
    lvc        <- fixed(log(986));  label("Apparent osimertinib central volume, V/F, at 62 kg and 39 g/L albumin (L)")  # Part SII $THETA 5 '(986) FIX'; Brown 2017 Table 2 Vparent/F = 986 L
    lcl_az5104 <- fixed(log(31.5)); label("Apparent AZ5104 clearance, CL/F, at 62 kg (L/h)")                       # Part SII $THETA 6 '(31.5) FIX'; Brown 2017 Table 2 CLmetabolite/F = 31.5 L/h
    lvc_az5104 <- fixed(log(207));  label("Apparent AZ5104 central volume, V/F (L)")                               # Part SII $THETA 7 '(207) FIX'; Brown 2017 Table 2 Vmetabolite/F = 207 L

    # Continuous-covariate power exponents, also inherited fixed.
    e_wt_cl        <- fixed(0.56); label("Power exponent for body weight on osimertinib CL/F (unitless)")  # Part SII $THETA 8 '(0.56) FIX'; Brown 2017 Table 2
    e_wt_vc        <- fixed(0.65); label("Power exponent for body weight on osimertinib V/F (unitless)")   # Part SII $THETA 9 '(0.65) FIX'; Brown 2017 Table 2
    e_wt_cl_az5104 <- fixed(0.99); label("Power exponent for body weight on AZ5104 CL/F (unitless)")       # Part SII $THETA 10 '(0.99) FIX'; Brown 2017 Table 2
    e_alb_vc       <- fixed(1.33); label("Power exponent for serum albumin on osimertinib V/F (unitless)") # Part SII $THETA 11 '(1.33) FIX'; Brown 2017 Table 2

    # The single parameter estimated in this model. It is the only theta
    # in the Part SII control stream without a FIX flag.
    e_cobi_cl <- 0.678; label("Multiplicative factor of concomitant cobicistat on osimertinib CL/F (unitless)")  # Part SII $THETA 12 '(0.678) ;12 cov cobi on CLp'; Westra 2025 Results 'Generalizability' reports the corresponding 32.2 percent CL/F reduction

    # Between-subject variability. The Part SII $OMEGA lower-triangular
    # block supplies, in ETA order, 0.46 (parent CL/F), covariance 0.44,
    # 0.52 (AZ5104 CL/F), then zero covariances and 0.89 (ka); the two
    # subsequent diagonal $OMEGA records supply 0.52 (parent V/F) and
    # 0.62 (AZ5104 V/F), both flagged FIX. Mapping each value to the ETA
    # index its $PK block references reproduces Brown 2017 Table 2
    # exactly - omega CL_parent 0.46, omega CL_metabolite 0.52, omega ka
    # 0.89, omega V_parent 0.52, omega V_metabolite 0.62 - which is how
    # the two trailing $OMEGA comments are known to be mis-numbered by
    # one relative to the ETAs their own $PK block uses. The code, not
    # the comment, is authoritative here.
    #
    # IMPORTANT scale note. Westra 2025 entered Brown 2017's reported
    # omega values into $OMEGA, where NONMEM reads them as VARIANCES;
    # the covariance 0.44 confirms this reading, because
    # 0.44 / sqrt(0.46 * 0.52) = 0.90, exactly the correlation Brown 2017
    # reports. The simulated spread of Westra 2025 Table S1 corroborates
    # it: a 90% interval of 10.2-97.4 for AUC0-144h implies a log-scale
    # SD of log(97.4/10.2)/(2*1.645) = 0.686, which matches sqrt(0.46) =
    # 0.678 and not 0.46. This is a larger between-subject variability
    # than the existing Brown_2017_osimertinib.R model file encodes,
    # because that file reads the same published numbers as omegas
    # (standard deviations) and squares them. Both readings are defensible
    # from Brown 2017 alone; this file reproduces what Westra 2025
    # actually ran, which is the variance reading. See the vignette
    # Errata.
    etalcl + etalcl_az5104 ~ c(0.46, 0.44, 0.52)                   # Part SII $OMEGA BLOCK, ETA(1) and ETA(2) with their covariance
    etalka        ~ 0.89                                           # Part SII $OMEGA BLOCK third diagonal, ETA(3), zero covariance with the two clearance etas
    etalvc        ~ fixed(0.52)                                    # Part SII $OMEGA '0.52 FIX', ETA(4), referenced by V1 in $PK
    etalvc_az5104 ~ fixed(0.62)                                    # Part SII $OMEGA '0.62 FIX', ETA(5), referenced by V2 in $PK

    # Residual error. The Part SII $ERROR block defines a single weight
    # W = IPRED*THETA(1)+THETA(2) with $SIGMA 1 FIX, i.e. a combined
    # proportional-plus-additive model applied to BOTH observed analytes
    # through the same two thetas, each flagged FIX. The additive term is
    # on the ug/L scale set by S2 = V1/1000 and S3 = V2/1000. Note that
    # the same numeric value 0.105 is reported by Brown 2017 on a nmol/L
    # scale; Westra 2025 carried the number across without a unit
    # conversion, so it is reproduced here as printed and flagged in the
    # vignette Errata.
    propSd        <- fixed(0.244); label("Proportional residual error on osimertinib (fraction)")                                         # Part SII $THETA 1 '(0.244) FIX'
    addSd         <- fixed(0.105); label("Additive residual error on osimertinib (ug/L)")                                                 # Part SII $THETA 2 '(0.105) FIX'
    propSd_az5104 <- fixed(0.244); label("Proportional residual error on AZ5104 (fraction; the same shared theta as osimertinib)")        # Part SII $THETA 1; $ERROR uses one W for both compartments
    addSd_az5104  <- fixed(0.105); label("Additive residual error on AZ5104 (ug/L; the same shared theta as osimertinib)")                # Part SII $THETA 2; $ERROR uses one W for both compartments
  })

  model({
    # Fraction of the osimertinib elimination flux appearing as AZ5104,
    # fixed at 0.25. Part SII $PK 'K23 = K20 * 0.25'; the $DES block
    # removes only K20*A(2) from the parent, so the metabolite formation
    # flux is a fraction OF the parent elimination flux and does not add
    # to it.
    fmet <- 0.25

    # Reference covariate values, Brown 2017 Table 1 medians, as used in
    # the Part SII $PK block.
    ref_wt  <- 62
    ref_alb <- 39

    # Individual parameters.
    ka <- exp(lka + etalka)

    cl <- exp(lcl + etalcl) *
          (WT / ref_wt)^e_wt_cl *
          e_cobi_cl^CONMED_COBICISTAT

    vc <- exp(lvc + etalvc) *
          (WT / ref_wt)^e_wt_vc *
          (ALB / ref_alb)^e_alb_vc

    cl_az5104 <- exp(lcl_az5104 + etalcl_az5104) *
                 (WT / ref_wt)^e_wt_cl_az5104

    vc_az5104 <- exp(lvc_az5104 + etalvc_az5104)

    # Micro-constants, named as in the Part SII $PK block.
    k20 <- cl / vc
    k30 <- cl_az5104 / vc_az5104
    k23 <- k20 * fmet

    # ODE system, amounts in mg. Reproduces the Part SII $DES block
    # exactly. As in the companion model, no molecular-weight correction
    # is applied between parent and metabolite, so the AZ5104 state and
    # concentration are in osimertinib mass equivalents.
    d/dt(depot)          <- -ka * depot
    d/dt(central)        <-  ka * depot - k20 * central
    d/dt(central_az5104) <-  k23 * central - k30 * central_az5104

    # Observations in ug/L, from the Part SII scaling S2 = V1/1000 and
    # S3 = V2/1000.
    Cc        <- 1000 * central / vc
    Cc_az5104 <- 1000 * central_az5104 / vc_az5104

    Cc        ~ prop(propSd) + add(addSd)
    Cc_az5104 ~ prop(propSd_az5104) + add(addSd_az5104)
  })
}
