Diao_2014_rFIXFc <- function() {
  description <- "Three-compartment population PK model for recombinant factor IX Fc fusion protein (rFIXFc, eftrenonacog alfa) in patients with severe to moderate haemophilia B aged 12-77 years (Diao 2014). Disposition is described by linear three-compartment kinetics with intravenous input and first-order elimination from the central compartment; body weight is the only retained covariate, scaling CL and V1 with estimated power exponents (not the canonical 0.75 / 1) and a reference weight of 73 kg."
  reference <- "Diao L, Li S, Ludden T, Gobburu J, Nestorov I, Jiang H. Population pharmacokinetic modelling of recombinant factor IX Fc fusion protein (rFIXFc) in patients with haemophilia B. Clin Pharmacokinet. 2014;53(5):467-477. doi:10.1007/s40262-013-0129-7"
  vignette <- "Diao_2014_rFIXFc"
  units <- list(time = "hour", dosing = "IU", concentration = "IU/dL")

  covariateData <- list(
    WT = list(
      description        = "Body weight (baseline)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference weight 73 kg (Diao 2014 Table 3, typical 73 kg patient). Body-weight power exponents were estimated (not fixed): 0.436 on CL and 0.396 on V1 (Diao 2014 Table 3 and Discussion p. 475). The paper found BW to be the only statistically significant covariate; its impact was limited (IIV on CL and V1 dropped by only 3.4 and 2.5 percentage points after including BW). BW does not enter Q2, V2, Q3, or V3 in the final model.",
      source_name        = "BW"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 135L,
    n_studies      = 2L,
    age_range      = "12.1-76.8 years",
    age_median     = "31.3 years",
    weight_range   = "45.0-186.7 kg",
    weight_median  = "73.3 kg",
    sex_female_pct = 0,
    race_ethnicity = c(White = 60.7, Asian = 22.2, Black = 8.9, Other = 7.4, AmericanIndianAlaskaNative = 0.74),
    disease_state  = "Previously treated patients with severe to moderate haemophilia B (endogenous FIX <= 2 IU/dL); HIV-positive in 3.7%, HCV-positive in 38.5%; FIX genotype 55.5% missense, 17.8% nonsense, 13.3% frameshift, 3.0% splice mutation, 10.4% other; haematocrit median 0.44, albumin median 46 g/L, IgG1 median 7.19 mg/mL.",
    dose_range     = "Phase 1/2a: 12.5, 25, 50, or 100 IU/kg rFIXFc single IV dose (~10-min infusion); phase 3 B-LONG: 50 or 100 IU/kg PK dose plus weekly (Arm 1, 20-100 IU/kg), individualised-interval (Arm 2, starting at 100 IU/kg), on-demand (Arm 3, 20-100 IU/kg), or perisurgical (Arm 4, 40-100 IU/kg) regimens.",
    regions        = "Multinational phase 3 (B-LONG, NCT01027364) plus single phase 1/2a study (NCT00716716).",
    n_observations = "1,400 FIX activity records from 135 baseline pharmacokinetic profiles plus 21 repeat profiles at week 26 in the Arm 1 sequential PK subgroup (modelling dataset); 1,027 trough/peak records (validation dataset).",
    reference_subject = "73 kg patient. (WT / 73)^0.436 for CL and (WT / 73)^0.396 for V1 both equal 1.",
    notes          = "Baseline characteristics in Diao 2014 Table 1; modelling/validation dataset summary in Table 2. Haemophilia B is X-linked and all subjects were male. Dependent variable was baseline- and residual-corrected FIX activity (Diao 2014 Eqs. 1-2), measured by one-stage aPTT clotting assay (LLOQ 1 IU/dL)."
  )

  ini({
    # Structural parameters - typical values for the reference 73 kg patient
    # (Diao 2014 Table 3). The paper reports CL, Q2, Q3 in dL/h and V1, V2, V3
    # in dL. Compartment 1 (central) holds V1 with clearance CL; compartments
    # 2 and 3 are the peripheral compartments, mapped to nlmixr2lib's
    # peripheral1 / peripheral2.
    lcl  <- log(2.39); label("Clearance for the reference 73 kg patient (CL, dL/h)")             # Diao 2014 Table 3: CL = 2.39 dL/h
    lvc  <- log(71.4); label("Central volume of distribution for the reference 73 kg patient (V1, dL)") # Diao 2014 Table 3: V1 = 71.4 dL
    lq   <- log(1.67); label("Intercompartmental clearance to peripheral 1 (Q2, dL/h)")          # Diao 2014 Table 3: Q2 = 1.67 dL/h
    lvp  <- log(87.0); label("Peripheral volume of distribution 1 (V2, dL)")                     # Diao 2014 Table 3: V2 = 87.0 dL
    lq2  <- log(39.3); label("Intercompartmental clearance to peripheral 2 (Q3, dL/h)")          # Diao 2014 Table 3: Q3 = 39.3 dL/h
    lvp2 <- log(39.9); label("Peripheral volume of distribution 2 (V3, dL)")                     # Diao 2014 Table 3: V3 = 39.9 dL

    # Body-weight power exponents - estimated by the paper rather than fixed to
    # the canonical 0.75 / 1 (Diao 2014 Discussion p. 475: "the exponents of BW
    # on CL and V1 were estimated during the modelling instead of being fixed
    # at presumed values, e.g. 0.75 for CL and 1 for V1"). Reference weight 73 kg.
    e_wt_cl <- 0.436; label("Power exponent of body weight on CL (unitless)")                    # Diao 2014 Table 3: BW exponent on CL = 0.436
    e_wt_vc <- 0.396; label("Power exponent of body weight on V1 (unitless)")                    # Diao 2014 Table 3: BW exponent on V1 = 0.396

    # Inter-individual variability. Diao 2014 Table 3 footnote: "IIV calculated
    # as sqrt(variance) * 100", i.e., the reported percentage equals the SD of
    # the log-scale eta times 100 and omega^2 = (IIV/100)^2 directly (no
    # log(CV^2 + 1) transform). Correlation between IIVs is reported as a
    # percentage on the variance/covariance scale, so cov = corr * sqrt(var1 * var2).
    # CL/V1 correlated block (Diao 2014 Table 3):
    #   omega^2_CL  = 0.177^2 = 0.031329
    #   cov_CL_V1   = 0.756 * 0.177 * 0.217 = 0.029042
    #   omega^2_V1  = 0.217^2 = 0.047089
    etalcl + etalvc ~ c(0.031329,
                        0.029042, 0.047089)   # Diao 2014 Table 3: IIV CL = 17.7%, IIV V1 = 21.7%, corr(CL,V1) = 75.6%

    # Independent IIVs on Q2, V2, V3 (Diao 2014 Table 3).
    # IIV on Q3 was not retained (high standard error 87%); see Diao 2014 Sec. 3.1.
    etalq   ~ 0.128164                        # Diao 2014 Table 3: IIV Q2 = 35.8%; omega^2 = 0.358^2
    etalvp  ~ 0.213444                        # Diao 2014 Table 3: IIV V2 = 46.2%; omega^2 = 0.462^2
    etalvp2 ~ 0.142129                        # Diao 2014 Table 3: IIV V3 = 37.7%; omega^2 = 0.377^2

    # Residual error - combined additive + proportional (Diao 2014 Table 3).
    propSd <- 0.106; label("Proportional residual error (fraction)")                             # Diao 2014 Table 3: proportional residual error = 10.6%
    addSd  <- 0.24;  label("Additive residual error on FIX activity (IU/dL)")                    # Diao 2014 Table 3: additive residual error = 0.24 IU/dL
  })

  model({
    # Individual PK parameters. BW enters CL and V1 only (Diao 2014 final
    # model); Q2, V2, Q3, V3 are not BW-scaled. Allometric exponents on CL and
    # V1 were estimated (e_wt_cl = 0.436, e_wt_vc = 0.396) - markedly lower
    # than the canonical 0.75 / 1, per the paper's Discussion.
    cl  <- exp(lcl  + etalcl)  * (WT / 73)^e_wt_cl
    vc  <- exp(lvc  + etalvc)  * (WT / 73)^e_wt_vc
    q   <- exp(lq   + etalq)
    vp  <- exp(lvp  + etalvp)
    q2  <- exp(lq2)
    vp2 <- exp(lvp2 + etalvp2)

    # Micro-constants (three-compartment, IV input to central, first-order
    # elimination from central; Diao 2014 Fig. 2).
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # ODE system. Dose enters central (rFIXFc is given IV); peripheral1 maps
    # to the paper's compartment 2 (Q2 / V2) and peripheral2 maps to compartment
    # 3 (Q3 / V3).
    d/dt(central)     <- -(kel + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # FIX activity in IU/dL. The model is fit to baseline- and residual-
    # corrected FIX activity (Diao 2014 Eqs. 1-2): observed activity minus
    # endogenous baseline (lowest activity per subject, set to 0 if <1 IU/dL)
    # minus decayed residual from prior FIX product. Cc therefore represents
    # the rFIXFc-attributable FIX activity above baseline.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
