Bosch_2022_liraglutide_qsp <- function() {
  description <- paste(
    "QSP. Integrated four-hormone-plus-glucose (4GI) model of human",
    "glucose regulation coupling glucose, insulin, GLP-1, glucagon,",
    "and glucose-dependent insulinotropic peptide (GIP) dynamics",
    "after meal and drug perturbation, with liraglutide (GLP-1RA)",
    "PK+PD baked in as the training drug. Glucose disposition is a",
    "two-compartment model with insulin-dependent and -independent",
    "elimination fed by an oral-meal dose compartment plus a buffer",
    "and three serial gut transit compartments; insulin, GLP-1,",
    "glucagon, and GIP each have one-compartment (GIP two-compartment)",
    "turnover with steady-state baseline production. Feedbacks: glucose",
    "stimulates insulin production (amplified by GLP-1 and GIP);",
    "insulin drives insulin-dependent glucose clearance through an",
    "effect compartment; glucose inhibits glucagon production (a",
    "below-baseline glucose-on-glucagon effect is estimable in healthy",
    "volunteers only); glucagon and GIP each stimulate glucagon and",
    "glucose production; GLP-1 inhibits gastric emptying and glucagon",
    "production, and stimulates insulin secretion. Liraglutide PK is",
    "the Watson 2010 body-weight-scaled one-compartment model (KAdrug",
    "0.154 /h, CL/F 0.013 L/h/kg, V/F 0.16 L/kg, fu 0.005) and drives",
    "the three GLP-1 receptor pathways via unbound-liraglutide EC50s",
    "derived from the in vitro EC50 ratio between endogenous GLP-1",
    "(1.919 pM) and liraglutide (6 pM). Population-mean fits only (no",
    "IIV per Bosch 2022 Methods; five proportional residual errors,",
    "one per biomarker output). Healthy volunteer (HV) vs type 2",
    "diabetes mellitus (T2DM) subjects differ in glucose clearance,",
    "insulin-dependent glucose clearance, the below-baseline",
    "glucose-on-glucagon exponent, and the GIP-on-insulin-secretion",
    "exponent, switched by the DIS_DIAB indicator (1 = T2DM, 0 = HV).",
    "18 ODE states; 5 outputs.",
    sep = " "
  )
  reference <- paste(
    "Bosch R, Petrone M, Arends R, Vicini P, Sijbrands EJG, Hoefman S,",
    "Snelder N. A novel integrated QSP model of in vivo human glucose",
    "regulation to support the development of a glucagon/GLP-1 dual",
    "agonist. CPT Pharmacometrics Syst Pharmacol. 2022;11(3):302-317.",
    "doi:10.1002/psp4.12752.",
    sep = " "
  )
  vignette <- "Bosch_2022_liraglutide_qsp"
  paper_specific_compartments <- c(
    "glc_dose", "glc_buffer", "glc_gut1", "glc_gut2", "glc_gut3",
    "glc_central", "glc_peripheral",
    "ins_central", "ins_effect",
    "glp_central", "glg_central",
    "gip_central", "gip_peripheral",
    "lira_depot", "lira_central"
  )

  units <- list(
    time          = "h",
    dosing        = "mmol (glucose meal) or pmol (liraglutide SC)",
    concentration = "mmol/L (glucose), pmol/L (insulin, GLP-1, glucagon, GIP, liraglutide)"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    glc_dose       = list(analyte = "glucose meal", units = NA_character_, specimen = "administration site", verified = FALSE),
    glc_buffer     = list(analyte = "glucose", units = NA_character_, specimen = "blood cell", verified = FALSE),
    glc_gut1       = list(analyte = "glucose", units = NA_character_, specimen = "administration site", verified = FALSE),
    glc_gut2       = list(analyte = "glucose", units = NA_character_, specimen = "administration site", verified = FALSE),
    glc_gut3       = list(analyte = "glucose", units = NA_character_, specimen = "administration site", verified = FALSE),
    glc_central    = list(analyte = "glucose", units = NA_character_, specimen = "plasma", verified = FALSE),
    glc_peripheral = list(analyte = "glucose", units = NA_character_, specimen = "tissue", verified = FALSE),
    ins_central    = list(analyte = "insulin", units = NA_character_, specimen = "plasma", verified = FALSE),
    ins_effect     = list(analyte = "insulin", units = NA_character_, specimen = "not applicable", verified = FALSE),
    glp_central    = list(analyte = "GLP-1", units = NA_character_, specimen = "plasma", verified = FALSE),
    glg_central    = list(analyte = "glucagon", units = NA_character_, specimen = "plasma", verified = FALSE),
    gip_central    = list(analyte = "GIP", units = NA_character_, specimen = "plasma", verified = FALSE),
    gip_peripheral = list(analyte = "GIP", units = NA_character_, specimen = "tissue", verified = FALSE),
    lira_depot     = list(analyte = "liraglutide", units = NA_character_, specimen = "administration site", verified = FALSE),
    lira_central   = list(analyte = "liraglutide", units = NA_character_, specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    DIS_DIAB = list(
      description        = "Type 2 diabetes mellitus indicator (1 = T2DM patient, 0 = healthy volunteer)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy volunteer)",
      notes              = paste(
        "Bosch 2022 uses a per-subject `PAT` flag in NONMEM (paper's",
        "convention: PAT=1 for healthy volunteer, PAT=0 for T2DM; see the",
        "supplementary NONMEM code Appendix S1 `IF (PAT.EQ.1.OR.STUD.EQ.2)",
        "THEN ; Healthy volunteers` block). Mapped to the canonical",
        "DIS_DIAB by inverting sign: DIS_DIAB = 0 (HV) uses the Silber",
        "2007 healthy-volunteer glucose disposition parameters (CLglc =",
        "5.36 L/h, CLglci = 0.072 (L/h)/(pmol/L)) plus the estimable",
        "below-baseline glucose-on-glucagon exponent POW_2L = 0.327 and",
        "the estimable GIP-on-insulin-secretion exponent POW_3 = 0.286;",
        "DIS_DIAB = 1 (T2DM) uses the Jauslin 2011 T2DM disposition",
        "parameters (CLglc = 1.72 L/h, CLglci = 0.0256 (L/h)/(pmol/L))",
        "with both POW_2L and POW_3 fixed to 0 (the paper Discussion:",
        "the below-baseline glucose-on-glucagon effect and the",
        "GIP-on-insulin-secretion amplification 'could only be estimated",
        "for the HV population' in T2DM). Time-fixed at study entry."
      ),
      source_name        = "PAT (inverted -- PAT=1 in the paper is HV, mapped here to DIS_DIAB=0)"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters as a linear multiplier on the Watson 2010 liraglutide",
        "PK parameters (CL/F = 0.013 L/h per kg body weight; V/F = 0.16",
        "L per kg body weight). Bosch 2022 quotes the Watson 2010 model",
        "in the supplementary NONMEM code as `CLdrug = 0.013 * BW0`,",
        "`VCdrug = 0.16 * BW0` -- i.e., a per-kg allometric-exponent-1",
        "linear form with no reference weight (the per-kg value IS the",
        "typical value). Only used when Cdrugf > 0 (liraglutide arm); a",
        "no-drug simulation can supply any positive WT without affecting",
        "the endogenous 4GI dynamics."
      ),
      source_name        = "BW0"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = 12L,
    age_range      = NA_character_,
    weight_range   = NA_character_,
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state  = paste(
      "Pooled healthy volunteers, healthy obese volunteers, and adults",
      "with type 2 diabetes mellitus from 12 published clinical studies",
      "used as development data (10 studies: Silber 2007 IVGTT +/-",
      "insulin in HV and T2DM; Jauslin 2011 24-h meal profiles in T2DM;",
      "Landersdorfer 2011 24-h meal profiles in T2DM; Tan 2012 IV GLP-1",
      "+/- glucagon in healthy obese; Edholm 2010 meal + IV GLP-1 or",
      "GIP in HV; Vilsboll 2006 IV GIP in HV and T2DM; Vilsboll 2002",
      "IV glucose + GLP-1 or GIP in HV and obese T2DM; Larsen 2001",
      "meal + IV GLP-1 in T2DM; Schneck 2013 meal in T2DM; Camastra",
      "2013 meal in morbidly-obese T2DM and non-diabetic controls) and",
      "3 long-term liraglutide studies (LEAD-3 52 weeks, LEAD-6 40",
      "weeks, AWARD-6 26 weeks); validation is against the dulaglutide",
      "arms of AWARD-6 and the dulaglutide + semaglutide arms of",
      "SUSTAIN-7 (also 40 weeks)."
    ),
    dose_range     = paste(
      "Nonpharmacological perturbations: 75 g oral glucose or standard",
      "meals (breakfast, lunch, dinner, snack); IV glucose infusions",
      "and clamps; IV endogenous-hormone infusions (GLP-1, glucagon,",
      "GIP at physiological doses). Pharmacological drug: subcutaneous",
      "liraglutide 1.2 or 1.8 mg once daily for up to 52 weeks in the",
      "LEAD-3 / LEAD-6 / AWARD-6 studies. External validation used",
      "dulaglutide 0.75 or 1.5 mg once weekly (AWARD-6, SUSTAIN-7) and",
      "semaglutide 0.5 or 1.0 mg once weekly (SUSTAIN-7)."
    ),
    regions        = NA_character_,
    notes          = paste(
      "Sample sizes are not tabulated in the aggregate paper because",
      "Bosch 2022 fits mean per-study time-course profiles digitized",
      "from published figures via DigitizeIt (Bosch 2022 Methods 'Data",
      "extraction'), not individual-level data -- hence the",
      "typical-value / no-IIV structure. Baselines (glucose, insulin,",
      "GLP-1, glucagon, GIP) were re-estimated per study; the packaged",
      "defaults are the HV-generic Jauslin-cohort fixed values plus a",
      "5 pmol/L GIP baseline (mid-range HV fasting GIP; the paper does",
      "not tabulate a single default). Users should override the",
      "baselines to match their target population and study protocol",
      "(the vignette walks the recommended HV vs T2DM defaults)."
    )
  )

  ini({
    # =====================================================================
    # GLUCOSE DISPOSITION (Bosch 2022 Table 3). Six parameters are fixed
    # to Silber 2007 (HV) / Jauslin 2011 (T2DM) IGI-model literature values
    # (T2DM values encoded as the reference DIS_DIAB=1 baseline; the HV
    # DIS_DIAB=0 CLglc / CLglci override is via the covariate effect).
    # Absorption rate constants (KAglc, Keglc, Kelglc) are estimated.
    # =====================================================================
    lclglc_t2dm  <- fixed(log(1.72))    ; label("Insulin-independent glucose clearance in T2DM CLglc (L/h)")                                     # Table 3, Fixed to Jauslin 2011
    lclglci_t2dm <- fixed(log(0.0256))  ; label("Insulin-dependent glucose clearance in T2DM CLglci ((L/h)/(pmol/L))")                            # Table 3, Fixed to Jauslin 2011
    lqglc        <- fixed(log(26.5))    ; label("Glucose intercompartmental clearance Qglc (L/h)")                                                # Table 3, Fixed to Silber 2007 / Jauslin 2011
    lvcglc       <- fixed(log(9.33))    ; label("Glucose central volume of distribution VCglc (L)")                                               # Table 3, Fixed to Silber 2007 / Jauslin 2011
    lvpglc       <- fixed(log(8.56))    ; label("Glucose peripheral volume of distribution VPglc (L)")                                            # Table 3, Fixed to Silber 2007 / Jauslin 2011
    lkaglc       <- log(0.853)          ; label("Glucose absorption rate from buffer to central KAglc (1/h)")                                     # Table 3, estimated (RSE 14.1%)
    lkeglc       <- log(0.281)          ; label("Glucose transit rate from buffer to gut1 Keglc (1/h)")                                           # Table 3, estimated (RSE 24.5%)
    lkelglc      <- log(1.93)           ; label("Glucose gut transit rate constant Kelglc (1/h)")                                                 # Table 3, estimated (RSE 15.3%)

    # HV-specific glucose clearance shifts (additive on log scale so
    # DIS_DIAB=0 gives the HV values CLglc=5.36 and CLglci=0.072).
    e_dis_diab0_lclglc  <- fixed(log(1.72 / 5.36))     ; label("Additive log-scale shift on CLglc for DIS_DIAB=0 (HV): CLglc_HV/CLglc_T2DM = 5.36/1.72")   # Table 3, Fixed to Silber 2007 HV CLglc = 5.36 L/h
    e_dis_diab0_lclglci <- fixed(log(0.0256 / 0.072))  ; label("Additive log-scale shift on CLglci for DIS_DIAB=0 (HV): CLglci_HV/CLglci_T2DM = 0.072/0.0256") # Table 3, Fixed to Silber 2007 HV CLglci = 0.072 (L/h)/(pmol/L)

    # =====================================================================
    # INSULIN DISPOSITION (Bosch 2022 Table 3). CLins and VCins fixed to
    # Silber 2007; Ke0ins (effect-compartment rate constant driving
    # insulin-dependent glucose clearance) estimated.
    # =====================================================================
    lclins <- fixed(log(73.2))          ; label("Insulin clearance CLins (L/h)")                                                                  # Table 3, Fixed to Silber 2007 / Jauslin 2011
    lvcins <- fixed(log(6.09))          ; label("Insulin central volume of distribution VCins (L)")                                               # Table 3, Fixed to Silber 2007 / Jauslin 2011
    lke0ins <- log(0.85)                ; label("Rate constant of insulin effect compartment on glucose clearance Ke0ins (1/h)")                  # Table 3, estimated (RSE 27.0%)

    # =====================================================================
    # GLP-1 DISPOSITION (Bosch 2022 Table 3). One-compartment with
    # Michaelis-Menten clearance (Vilsboll 2003). All three parameters
    # estimated; FACglp converts active GLP-1 (modelled state) to total
    # GLP-1 (for studies that measured total instead of active).
    # =====================================================================
    lvcglp <- log(16)                   ; label("GLP-1 central volume of distribution VCglp (L)")                                                 # Table 3, estimated (RSE 29.3%)
    lvmglp <- log(2893)                 ; label("GLP-1 Michaelis-Menten maximum clearance rate VMGLP (pmol/L*h)")                                 # Table 3, estimated (RSE 7.41%); NONMEM $THETA EXP(7.96952) = 2893
    lkmglp <- log(136)                  ; label("GLP-1 Michaelis-Menten half-saturation KMGLP (pmol/L)")                                          # Table 3, estimated (RSE 10.3%); NONMEM $THETA EXP(4.90603) = 135
    # NOTE: Bosch 2022 Table 3 also reports a total-to-active GLP-1
    # conversion factor `factor total GLP` = 3.8 (RSE 16.9%). This factor
    # is used in the paper's NONMEM $ERROR block ONLY for studies that
    # measured total (not active) GLP-1 -- specifically Tan 2012 (STUD.EQ.2),
    # where `IPRED = FACglp * A(4)/VCglp`. The packaged model observes the
    # active-GLP-1 concentration Cglp directly; users who want to compare
    # to a total-GLP-1 assay should multiply Cglp by 3.8. Not included as
    # an ini() parameter because it is not part of the ODE system and
    # is not referenced in model().

    # =====================================================================
    # GLUCAGON DISPOSITION (Bosch 2022 Table 3). One-compartment.
    # =====================================================================
    lclglg <- log(453)                  ; label("Glucagon clearance CLglg (L/h)")                                                                 # Table 3, estimated (RSE 5.55%)
    lvcglg <- log(64.6)                 ; label("Glucagon central volume of distribution VCglg (L)")                                              # Table 3, estimated (RSE 11.7%)

    # =====================================================================
    # GIP DISPOSITION (Bosch 2022 Table 3). Two-compartment; all four
    # parameters fixed to a re-fit of Vilsboll 2006 IV GIP data (per Bosch
    # 2022 Methods, kept fixed during subsequent 4GI optimization).
    # =====================================================================
    lclgip <- fixed(log(86.8))          ; label("GIP clearance CLgip (L/h)")                                                                      # Table 3, Fixed to a re-fit of Vilsboll 2006
    lvcgip <- fixed(log(9.21))          ; label("GIP central volume of distribution VCgip (L)")                                                   # Table 3, Fixed to a re-fit of Vilsboll 2006
    lqgip  <- fixed(log(49.4))          ; label("GIP intercompartmental clearance Qgip (L/h)")                                                    # Table 3, Fixed to a re-fit of Vilsboll 2006
    lvpgip <- fixed(log(22.8))          ; label("GIP peripheral volume of distribution VPgip (L)")                                                # Table 3, Fixed to a re-fit of Vilsboll 2006

    # =====================================================================
    # FOOD-EFFECT (KIN modulation) PARAMETERS (Bosch 2022 Table 4).
    # =====================================================================
    fdglp   <- 0.0102                   ; label("Food effect on GLP-1 production via glucose buffer FDGLP (1/mmol)")                              # Table 4, estimated (RSE 13.3%)
    fdglp_2 <- fixed(3.88)              ; label("Food effect on GLP-1 production via glucose gut3 FDGLP_2 (1/mmol) [HV only]")                    # Table 4, RSE dash (fixed). NB: NONMEM supplement Appendix S1 $THETA(23) shows `(0) FIX` in the MAXEVAL=0 rerun; Table 4 reports 3.88 (fixed) as the published final value -- trusting the paper table over the supplement code per policy.
    fdgip   <- 0.0343                   ; label("Food effect on GIP production via glucose buffer FDGIP (1/mmol)")                                # Table 4, estimated (RSE 19.4%); NONMEM code labels as (1/pmol) -- text and Table 4 label as (1/mmol) referring to buffer glucose amount (mmol)
    fdglg   <- 0.00329                  ; label("Food effect on glucagon production via glucose buffer FDGLG (1/mmol)")                           # Table 4, estimated (RSE 23.3%)

    # =====================================================================
    # GLUCOSE, GLP-1, GLUCAGON, GIP EFFECT PARAMETERS (Bosch 2022 Table 4).
    # =====================================================================
    glcins_s <- 2.46                    ; label("Power of glucose on insulin production GLCINS_S (unitless; per-mM)")                             # Table 4, estimated (RSE 7.56%); Cglc appears in eq. as Cglc**GLCINS_S so this is a power exponent
    pow_2h   <- 0.925                   ; label("Power of glucose-on-glucagon effect for glucose >= baseline POW_2H (unitless)")                  # Table 4, estimated (RSE 16.4%)
    pow_2l   <- fixed(0.327)            ; label("Power of glucose-on-glucagon effect for glucose < baseline in HV POW_2L (unitless) [HV only]")   # Table 4 not listed directly; NONMEM $THETA(48) 0.327 FIX (Bosch 2022 Discussion: 'could only be estimated for the HV population'); T2DM below-baseline effect = 0 (Table 4 GLCGLG_POWL, POW_2, T2DM = 0)

    # GLP-1 on glucose-dependent insulin secretion (Emax + Hill)
    emax_1  <- 10.7                     ; label("GLP-1 stimulation of insulin secretion Emax EMAX_1 (unitless)")                                  # Table 4, estimated (RSE 36.3%)
    ec50_1  <- 26.6                     ; label("GLP-1 half-max concentration for insulin secretion EC50_1 (pmol/L)")                             # Table 4, estimated (RSE 11.5%)
    hill_1  <- 1.79                     ; label("GLP-1 Hill coefficient for insulin secretion HILL_1 (unitless)")                                 # Table 4, estimated (RSE 18.4%)

    # GLP-1 on glucose absorption (gastric-emptying inhibition; Emax + Hill)
    emax_2  <- fixed(1)                 ; label("GLP-1 inhibition of glucose absorption Emax EMAX_2 (unitless)")                                  # Table 4, RSE dash (fixed to 1 = complete inhibition asymptote)
    ec50_2  <- 144                      ; label("GLP-1 half-max concentration for glucose absorption inhibition EC50_2 (pmol/L)")                 # Table 4, estimated (RSE 36.7%)
    hill_2  <- fixed(1)                 ; label("GLP-1 Hill coefficient for glucose absorption inhibition HILL_2 (unitless)")                     # Table 4, RSE dash (fixed to 1)

    # GLP-1 on glucagon production (inhibition; Emax + Hill)
    emax_3  <- fixed(1)                 ; label("GLP-1 inhibition of glucagon production Emax EMAX_3 (unitless)")                                 # Table 4, RSE dash (fixed to 1 = complete inhibition asymptote)
    ec50_3  <- 99.5                     ; label("GLP-1 half-max concentration for glucagon inhibition EC50_3 (pmol/L)")                           # Table 4, estimated (RSE 9.39%)
    hill_3  <- fixed(1)                 ; label("GLP-1 Hill coefficient for glucagon inhibition HILL_3 (unitless)")                               # Table 4, RSE dash (fixed to 1)

    # Glucagon on glucose production (stimulation; Emax + Hill)
    emax_4  <- 6.73                     ; label("Glucagon stimulation of glucose production Emax EMAX_4 (unitless)")                              # Table 4, estimated (RSE 31.8%)
    ec50_4  <- fixed(98.5)              ; label("Glucagon half-max concentration for glucose production EC50_4 (pmol/L)")                         # Table 4, fixed to Wendt 2007 (Bosch 2022 Results: 'EC50 for this effect was fixed to the estimate from Wendt et al. of 82 pM'; the paper's Table 4 reports 98.5 pM as the Wendt average used). NONMEM $THETA(41) = exp(4.59) = 98.5.
    hill_4  <- fixed(1)                 ; label("Glucagon Hill coefficient for glucose production HILL_4 (unitless)")                             # Table 4, RSE dash (fixed to 1)

    # GIP effects (power on insulin secretion for HV only; power on glucagon)
    pow_3   <- fixed(0.286)             ; label("GIP power on glucose-dependent insulin secretion POW_3 (unitless) [HV only]")                    # Table 4, estimated in HV subset (RSE 22.3%); T2DM POW_3 = 0 (Bosch 2022 Discussion: 'GIP effect on glucose-dependent insulin secretion could only be identified for HVs')
    pow_4   <- 0.109                    ; label("GIP power on glucagon production POW_4 (unitless)")                                              # Table 4, estimated (RSE 25.4%)

    # =====================================================================
    # LIRAGLUTIDE PK + POTENCY (Bosch 2022 Table 2 + supplementary NONMEM
    # code; Watson 2010 body-weight-scaled 1-compartment model).
    # =====================================================================
    lkadrug     <- fixed(log(0.154))    ; label("Liraglutide absorption rate constant KAdrug (1/h)")                                              # NONMEM supplement Appendix S1: 'PK Model Watson et al. 2010' with KAdrug = 0.154 (paper's Watson-2010 T2DM final value)
    lcldrug_kg  <- fixed(log(0.013))    ; label("Liraglutide apparent clearance per kg body weight CLdrug/BW (L/h/kg)")                           # NONMEM supplement Appendix S1: CLdrug = 0.013 * BW0 (Watson 2010)
    lvcdrug_kg  <- fixed(log(0.16))     ; label("Liraglutide apparent volume of distribution per kg body weight VCdrug/BW (L/kg)")                # NONMEM supplement Appendix S1: VCdrug = 0.16 * BW0 (Watson 2010)
    fu_lira     <- fixed(0.005)         ; label("Liraglutide unbound fraction in plasma fu (unitless)")                                           # NONMEM supplement Appendix S1: fu = 0.005; Bosch 2022 Table 2 free fraction = 0.51% = 0.0051 (Flint et al., ref 16 [47 in reprint])
    ec50_lira_in_vitro <- fixed(6)      ; label("Liraglutide in vitro EC50 at the GLP-1 receptor (pmol/L)")                                       # Bosch 2022 Table 2, in-house assay in absence of albumin
    ec50_endo_glp_in_vitro <- fixed(1.919) ; label("Endogenous GLP-1(7-36)NH2 in vitro EC50 at the GLP-1 receptor (pmol/L)")                      # Bosch 2022 Table 2, in-house assay in absence of albumin (NONMEM comment 'ECGLP = 1.919')

    # =====================================================================
    # BASELINES (Bosch 2022 Table S1 in supplement + NONMEM code $THETA).
    # Defaults are the Jauslin-cohort HV-generic fixed values (glucose 4.65
    # mM, insulin 49.1 pM, GLP-1 4.61 pM, glucagon 8.85 pM) and a
    # physiologically-typical HV fasting GIP baseline of 5 pmol/L. Users
    # simulating T2DM cohorts should override with study-specific baselines
    # (e.g., Landersdorfer glucose 10.8 mM, insulin 131 pM).
    # =====================================================================
    lbslglc <- log(4.65)                ; label("Baseline glucose concentration BSLglc (mmol/L)")                                                 # NONMEM $THETA(52) = 4.65 FIX (Jauslin HV), Bosch 2022 supplementary Table S1
    lbslins <- log(49.1)                ; label("Baseline insulin concentration BSLins (pmol/L)")                                                 # NONMEM $THETA(56) = 49.1 FIX (Jauslin HV), Bosch 2022 supplementary Table S1
    lbslglp <- log(4.61)                ; label("Baseline active GLP-1 concentration BSLglp (pmol/L)")                                            # NONMEM $THETA(59) = 4.61 FIX (Jauslin HV)
    lbslglg <- log(8.85)                ; label("Baseline glucagon concentration BSLglg (pmol/L)")                                                # NONMEM $THETA(61) = 8.85 FIX (Jauslin HV)
    lbslgip <- log(5)                   ; label("Baseline GIP concentration BSLgip (pmol/L)")                                                     # Physiological HV fasting GIP; Bosch 2022 uses study-specific IBGIP so no fixed default is published -- 5 pmol/L is mid-range for Vilsboll 2006 healthy volunteers

    # =====================================================================
    # RESIDUAL ERROR (Bosch 2022 Table 3; proportional per output).
    # =====================================================================
    propSd_Cglc <- 0.0211               ; label("Proportional residual SD for glucose (fraction)")                                                # Table 3, sigma1 = 0.0211 (RSE 27.5%); NONMEM $SIGMA linear-space proportional
    propSd_Cins <- 0.305                ; label("Proportional residual SD for insulin (fraction)")                                                # Table 3, sigma2 = 0.305 (RSE 38.4%)
    propSd_Cglp <- 0.0602               ; label("Proportional residual SD for GLP-1 (fraction)")                                                  # Table 3, sigma3 = 0.0602 (RSE 36.0%)
    propSd_Cglg <- 0.0348               ; label("Proportional residual SD for glucagon (fraction)")                                               # Table 3, sigma4 = 0.0348 (RSE 27.4%)
    propSd_Cgip <- 0.109                ; label("Proportional residual SD for GIP (fraction)")                                                    # Table 3, sigma5 = 0.109 (RSE 32.6%)
  })

  model({
    # ---------------------------------------------------------------------
    # Population-mean parameters (typical values; no IIV in this paper).
    # Structural parameters are back-transformed from their log-scale
    # ini() entries. DIS_DIAB switches HV vs T2DM population parameters:
    # DIS_DIAB = 0 (HV) -> CLglc_HV, CLglci_HV, plus POW_2L and POW_3
    #                      contribute their fixed values 0.327 and 0.286.
    # DIS_DIAB = 1 (T2DM) -> CLglc_T2DM, CLglci_T2DM, POW_2L and POW_3
    #                        forced to 0 (paper Discussion: those effects
    #                        could only be identified for HV).
    # ---------------------------------------------------------------------
    clglc  <- exp(lclglc_t2dm  + e_dis_diab0_lclglc  * (1 - DIS_DIAB))
    clglci <- exp(lclglci_t2dm + e_dis_diab0_lclglci * (1 - DIS_DIAB))
    qglc   <- exp(lqglc)
    vcglc  <- exp(lvcglc)
    vpglc  <- exp(lvpglc)
    kaglc  <- exp(lkaglc)
    keglc  <- exp(lkeglc)
    kelglc <- exp(lkelglc)

    clins  <- exp(lclins)
    vcins  <- exp(lvcins)
    ke0ins <- exp(lke0ins)

    vcglp  <- exp(lvcglp)
    vmglp  <- exp(lvmglp)
    kmglp  <- exp(lkmglp)

    clglg  <- exp(lclglg)
    vcglg  <- exp(lvcglg)

    clgip  <- exp(lclgip)
    vcgip  <- exp(lvcgip)
    qgip   <- exp(lqgip)
    vpgip  <- exp(lvpgip)

    # HV-only sub-effects: gated on (1 - DIS_DIAB) so T2DM sees exactly 0.
    pow_2l_i <- pow_2l * (1 - DIS_DIAB)
    pow_3_i  <- pow_3  * (1 - DIS_DIAB)
    # HV-only delayed GLP-1 food effect from the third gut transit
    # compartment (NONMEM code: `IF (A(15).GT.0.AND.PAT.EQ.1)`); PAT.EQ.1
    # is HV in the paper's convention, so this activates when DIS_DIAB=0.
    fdglp_2_i <- fdglp_2 * (1 - DIS_DIAB)

    # Baselines (mmol/L or pmol/L)
    bslglc <- exp(lbslglc)
    bslins <- exp(lbslins)
    bslglp <- exp(lbslglp)
    bslglg <- exp(lbslglg)
    bslgip <- exp(lbslgip)

    # Liraglutide PK: per-kg apparent clearance and volume.
    kadrug     <- exp(lkadrug)
    cldrug     <- exp(lcldrug_kg) * WT
    vcdrug     <- exp(lvcdrug_kg) * WT

    # Liraglutide in vivo EC50 (unbound) at each of the three GLP-1
    # pathways is derived from the endogenous in vivo EC50 by the ratio
    # of the drug's in vitro EC50 to the endogenous GLP-1 in vitro EC50
    # (Bosch 2022 Eqs. 28-30, assuming the in vitro:in vivo EC50 ratio
    # is the same for endogenous and exogenous agonists).
    ecglp1_lira <- ec50_lira_in_vitro / ec50_endo_glp_in_vitro * ec50_1
    ecglp2_lira <- ec50_lira_in_vitro / ec50_endo_glp_in_vitro * ec50_2
    ecglp3_lira <- ec50_lira_in_vitro / ec50_endo_glp_in_vitro * ec50_3

    # ---------------------------------------------------------------------
    # Compartmental concentrations. Amounts are held in the state; concen-
    # trations are algebraic. Glucose amounts in mmol -> mmol/L; peptide
    # amounts in pmol -> pmol/L; liraglutide amounts in pmol -> pmol/L.
    # ---------------------------------------------------------------------
    cglc  <- glc_central / vcglc
    cins  <- ins_central / vcins
    cglp  <- glp_central / vcglp
    cglg  <- glg_central / vcglg
    cgip  <- gip_central / vcgip
    cdrug <- lira_central / vcdrug

    # Unbound liraglutide (used in GLP-1-mediated drug effects only)
    cdrugf <- cdrug * fu_lira

    # ---------------------------------------------------------------------
    # Effect equations (Bosch 2022 Eqs. 2-24, 28-30).
    # ---------------------------------------------------------------------

    # GLP-1 (endogenous + liraglutide) on insulin secretion (stimulation)
    # Bosch 2022 Eq. 10-12 + Eq. 28: shared Emax / Hill; the drug is
    # added competitively with endogenous inside the Hill saturation.
    glpins_s_endo_num_bsl <- (bslglp / ec50_1) ^ hill_1
    glpins_s0 <- emax_1 * glpins_s_endo_num_bsl / (1 + glpins_s_endo_num_bsl)
    glpins_num_endo    <- (cglp   / ec50_1) ^ hill_1
    glpins_num_drug    <- (cdrugf / ecglp1_lira) ^ hill_1
    glpins_num_total   <- glpins_num_endo + glpins_num_drug
    glpins_s  <- emax_1 * glpins_num_total / (1 + glpins_num_total)

    # GIP on insulin secretion (power on Cgip; steady-state anchor)
    gipins_s  <- (cgip   ^ pow_3_i)
    gipins_s0 <- (bslgip ^ pow_3_i)

    # Combined incretin stimulation of insulin production (per Bosch 2022
    # Eq. 10-12: STglc = GLPINS_S + GIPINS_S adds the two effects).
    stglc  <- glpins_s  + gipins_s
    stglc0 <- glpins_s0 + gipins_s0

    # GLP-1 on gastric emptying / glucose absorption (inhibition; Bosch
    # 2022 Eqs. 1-2 + Eq. 28). GLP inhibits Kaglc2 = Kaglc * (1 - GLPGLU_AI).
    glpglu_ai_num_endo  <- (cglp   / ec50_2) ^ hill_2
    glpglu_ai_num_drug  <- (cdrugf / ecglp2_lira) ^ hill_2
    glpglu_ai_num_total <- glpglu_ai_num_endo + glpglu_ai_num_drug
    glpglu_ai <- emax_2 * glpglu_ai_num_total / (1 + glpglu_ai_num_total)
    kaglc2    <- kaglc * (1 - glpglu_ai)

    # GLP-1 on glucagon production (inhibition; Bosch 2022 Eqs. 21-23 +
    # Eq. 28). Steady-state anchor at BSLglp.
    glpglg_i0_num <- (bslglp / ec50_3) ^ hill_3
    glpglg_i0 <- emax_3 * glpglg_i0_num / (1 + glpglg_i0_num)
    glpglg_i_num_endo  <- (cglp   / ec50_3) ^ hill_3
    glpglg_i_num_drug  <- (cdrugf / ecglp3_lira) ^ hill_3
    glpglg_i_num_total <- glpglg_i_num_endo + glpglg_i_num_drug
    glpglg_i  <- emax_3 * glpglg_i_num_total / (1 + glpglg_i_num_total)
    # Multiplicative modulator, normalised so glpEFFglg = 1 at steady state
    glpEFFglg <- (1 - glpglg_i) / (1 - glpglg_i0)

    # Glucagon on glucose production (stimulation; Bosch 2022 Eqs. 5-7).
    glgglc_s0_num <- (bslglg / ec50_4) ^ hill_4
    glgglc_s0 <- emax_4 * glgglc_s0_num / (1 + glgglc_s0_num)
    glgglc_s_num <- (cglg / ec50_4) ^ hill_4
    glgglc_s  <- emax_4 * glgglc_s_num / (1 + glgglc_s_num)
    glgEFFglc <- (1 + glgglc_s) / (1 + glgglc_s0)

    # GIP on glucagon production (stimulation; Bosch 2022 Eq. 24). Power
    # model normalised at steady state (gipEFFglg = 1 at cgip = bslgip).
    gipEFFglg <- (cgip / bslgip) ^ pow_4

    # Glucose on glucagon production (Bosch 2022 Eq. 20 in the paper /
    # `glcEFFglg` in the NONMEM code): the exponent is POW_2H when
    # glucose is at or above baseline, and POW_2L when below (HV only;
    # T2DM POW_2L = 0 -> glcEFFglg = 1 below baseline).
    pow_2_i <- pow_2h
    if (cglc < bslglc) pow_2_i <- pow_2l_i
    glcEFFglg <- (bslglc / cglc) ^ pow_2_i

    # Food effects on baseline production rates (Bosch 2022 Eqs. 15-16,
    # 19, 27). Only active when the buffer/gut compartment has mass.
    fdglp_s   <- fdglp   * glc_buffer
    fdglp_s2  <- fdglp_2_i * glc_gut3
    fdgip_s   <- fdgip   * glc_buffer
    fdglg_s   <- fdglg   * glc_buffer

    # ---------------------------------------------------------------------
    # Baseline production fluxes derived from the steady-state anchors
    # (Bosch 2022 Eqs. 4, 9, 14, 18, 26).
    # ---------------------------------------------------------------------
    kinglc <- bslglc * (clglc + clglci * bslins)
    kinins <- bslins * clins / (1 + stglc0 * (bslglc ^ glcins_s))
    kinglp <- vmglp * (bslglp * vcglp) / (kmglp + bslglp)
    kinglg <- bslglg * clglg
    kingip <- bslgip * clgip

    # ---------------------------------------------------------------------
    # Micro-constants for two-compartment glucose and GIP distribution.
    # ---------------------------------------------------------------------
    k27  <- qglc / vcglc
    k72  <- qglc / vpglc
    k612 <- qgip / vcgip
    k126 <- qgip / vpgip

    # ---------------------------------------------------------------------
    # ODE system (Bosch 2022 Eqs. 1-3, 8, 13, 17, 25 + drug PK).
    # Amounts (mmol for glucose, pmol for peptides and drug).
    #
    # Compartment ordering here defines the rxode2 slot indices. Dosing
    # events must reference these names, NOT the algebraic concentrations
    # (Cglc, Cins, Cglp, Cglg, Cgip, Cdrug).
    # ---------------------------------------------------------------------

    # Glucose dose -> buffer -> central (with parallel gut transit loss)
    d/dt(glc_dose)   <- -kaglc2 * glc_dose
    d/dt(glc_buffer) <- kaglc2 * glc_dose - kaglc * glc_buffer - keglc * glc_buffer
    d/dt(glc_gut1)   <- keglc * glc_buffer - kelglc * glc_gut1
    d/dt(glc_gut2)   <- kelglc * (glc_gut1 - glc_gut2)
    d/dt(glc_gut3)   <- kelglc * (glc_gut2 - glc_gut3)

    # Glucose central: absorption + EGP*glucagon-mod - insulin-dep + insulin-indep clearance +/- peripheral exchange
    d/dt(glc_central)    <- kaglc * glc_buffer + kinglc * glgEFFglc - k27 * glc_central + k72 * glc_peripheral - (clglc / vcglc) * glc_central - (clglci * ins_effect / vcglc) * glc_central
    d/dt(glc_peripheral) <- k27 * glc_central - k72 * glc_peripheral

    # Insulin: production stimulated by glucose (incretin-amplified) - central elimination
    d/dt(ins_central) <- kinins * (1 + stglc * (cglc ^ glcins_s)) - (clins / vcins) * ins_central
    # Insulin effect compartment: first-order equilibration with Cins
    d/dt(ins_effect)  <- ke0ins * (cins - ins_effect)

    # GLP-1: baseline production x (1 + buffer effect) x (1 + gut3 effect) - Michaelis-Menten clearance
    d/dt(glp_central) <- kinglp * (1 + fdglp_s) * (1 + fdglp_s2) - vmglp * glp_central / (kmglp + cglp)

    # Glucagon: baseline production x (1 + buffer effect) x glcEFFglg x glpEFFglg x gipEFFglg - clearance
    d/dt(glg_central) <- kinglg * (1 + fdglg_s) * glcEFFglg * glpEFFglg * gipEFFglg - (clglg / vcglg) * glg_central

    # GIP: baseline production x (1 + buffer effect) - central clearance - central-to-peripheral + return
    d/dt(gip_central)    <- kingip * (1 + fdgip_s) - (clgip / vcgip) * gip_central - k612 * gip_central + k126 * gip_peripheral
    d/dt(gip_peripheral) <- k612 * gip_central - k126 * gip_peripheral

    # Liraglutide PK: first-order absorption + one-compartment elimination
    d/dt(lira_depot)   <- -kadrug * lira_depot
    d/dt(lira_central) <-  kadrug * lira_depot - (cldrug / vcdrug) * lira_central

    # ---------------------------------------------------------------------
    # Initial conditions for endogenous species: steady-state amounts.
    # Setting a steady-state initial condition here means simulations
    # without a dose start at the biological baseline instead of zero.
    # Users overriding the baselines above will get initial states
    # automatically re-computed from the new baseline values.
    # ---------------------------------------------------------------------
    glc_central(0)    <- bslglc * vcglc
    glc_peripheral(0) <- bslglc * vpglc
    ins_central(0)    <- bslins * vcins
    ins_effect(0)     <- bslins
    glp_central(0)    <- bslglp * vcglp
    glg_central(0)    <- bslglg * vcglg
    gip_central(0)    <- bslgip * vcgip
    gip_peripheral(0) <- bslgip * vpgip

    # ---------------------------------------------------------------------
    # Observables and residual error. Concentrations are already computed
    # above; declare each as an observable with a proportional residual.
    # ---------------------------------------------------------------------
    Cglc <- cglc
    Cins <- cins
    Cglp <- cglp
    Cglg <- cglg
    Cgip <- cgip
    Cdrug_lira <- cdrug

    Cglc ~ prop(propSd_Cglc)
    Cins ~ prop(propSd_Cins)
    Cglp ~ prop(propSd_Cglp)
    Cglg ~ prop(propSd_Cglg)
    Cgip ~ prop(propSd_Cgip)
  })
}
