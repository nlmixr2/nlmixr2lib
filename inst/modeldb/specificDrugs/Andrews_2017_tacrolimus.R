Andrews_2017_tacrolimus <- function() {
  description <- "Two-compartment population PK model with first-order absorption and an absorption lag time for twice-daily oral immediate-release tacrolimus (Prograft and Modigraf) in paediatric renal transplant recipients during the first 6 weeks post-transplantation (Andrews 2017). Apparent oral clearance CL/F and apparent inter-compartmental clearance Q/F scale allometrically with body weight at a fixed exponent of 0.75 referenced to a 70 kg adult; apparent central volume V1/F and apparent peripheral volume V2/F scale at a fixed exponent of 1.0; ka has no body-weight scaling. CL/F additionally varies with CYP3A5 expresser status (1.04 multiplier for *3/*3 or unknown genotype, 1.98 multiplier for *1/*1 or *1/*3 carriers; pooled with unknown because Andrews 2017 explicitly groups *3/*3 with unknown in the final equation), donor source (0.74 multiplier for living-donor recipients vs deceased-donor reference; equivalent to deceased-donor recipients having ~35% higher CL/F), eGFR (power exponent 0.19 centred at the cohort median 69 mL/min/1.73 m^2 of adapted-Schwartz eGFR), and a piecewise hematocrit effect (power exponent -0.44 centred at 0.3 L/L applied only when HCT < 0.3 L/L). Inter-individual variability is diagonal on ka, CL/F, V1/F, and V2/F. Residual error is a combined additive + proportional model with separate immunoassay and LC-MS/MS magnitudes selected by the per-sample IMMUNOASSAY indicator. Inter-occasion variability (IOV) on CL/F (18% CV) and V2/F (35% CV) reported by Andrews 2017 Table 2 is NOT encoded structurally here (per the Brooks 2021 tacrolimus precedent) -- the source paper does not define an operational occasion column for the model-library use case; downstream users who want to simulate IOV can add an OCC indicator and a per-occasion eta in rxode2."
  reference <- "Andrews LM, Hesselink DA, van Gelder T, Koch BCP, Cornelissen EAM, Bruggemann RJM, van Schaik RHN, de Wildt SN, Cransberg K, de Winter BCM. A Population Pharmacokinetic Model to Predict the Individual Starting Dose of Tacrolimus Following Pediatric Renal Transplantation. Clin Pharmacokinet. 2018;57(4):475-489. doi:10.1007/s40262-017-0567-8 (published online 5 July 2017)."
  vignette <- "Andrews_2017_tacrolimus"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying within subject during the first 6 weeks post-transplant. Allometric power scaling with reference weight 70 kg and theory-based fixed exponents: 0.75 on CL/F and Q/F, 1.0 on V1/F and V2/F; ka has no body-weight scaling. Andrews 2017 Section 3.1 confirms the fixed allometric exponents significantly improved the model and that estimating them did not (values not significantly different from theoretical). Model-building cohort median 28.4 kg, range 11.6-83.7 kg (Table 1).",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Adapted-Schwartz estimated glomerular filtration rate, BSA-normalized",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Stored under the canonical CRCL column per inst/references/covariate-columns.md (CRCL accepts BSA-normalised renal-function values regardless of the underlying derivation formula -- precedent: Cirincione 2017 MDRD eGFR, Li 2019 calculated GFR, Xu 2019 measured-CrCl BSA-normalised). Source name eGFR. Andrews 2017 Section 2.2 calculated eGFR using the adapted-Schwartz formula: eGFR = K * height(cm) / serum_creatinine(umol/L) with K = 36.5 (paediatric K constant). Reference value 69 mL/min/1.73 m^2 is the model-building cohort median (Table 1). Time-varying within subject during the first 6 weeks post-transplant; tacrolimus undergoes negligible renal elimination so the eGFR effect on CL/F is hypothesised to reflect intra-renal CYP3A5 metabolism rather than glomerular filtration of tacrolimus itself (Andrews 2017 Section 4 / Discussion).",
      source_name        = "eGFR"
    ),
    HCT = list(
      description        = "Hematocrit, expressed as a fraction of total blood volume (0-1, L/L)",
      units              = "fraction",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying during the first 6 weeks post-transplant. Andrews 2017 reports HCT as a fraction (L/L), NOT as percent. The canonical-register HCT entry's units (%) are explicitly overridden here so the centring value 0.3 L/L and the piecewise activation threshold of 0.3 L/L reproduce the paper's equation directly. To use a dataset that records HCT in percent, multiply the column by 0.01 before passing it to this model. Enters CL/F via the piecewise power form: f_HCT = (HCT/0.3)^(-0.44) if HCT < 0.3 L/L, else f_HCT = 1.0; the effect operates only on the anaemic subset because higher hematocrit binds more tacrolimus to erythrocytes (~70-80% of whole-blood tacrolimus is RBC-bound) and thus reduces whole-blood apparent CL/F. Model-building cohort median 0.29 L/L, range 0.16-0.43 L/L (Table 1).",
      source_name        = "hematocrit"
    ),
    CYP3A5_EXPR = list(
      description        = "CYP3A5 expresser indicator: 1 if the patient carries at least one functional CYP3A5*1 allele (genotype *1/*1 or *1/*3), 0 if homozygous *3/*3 or if the CYP3A5 genotype was not determined.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A5 *3/*3 nonexpresser, or genotype unknown -- Andrews 2017 pools these two groups in the final-model equation)",
      notes              = "Time-fixed germline genotype determined from rs776746 (CYP3A5 6986A>G). In the Andrews 2017 model-building cohort (n = 46), the CYP3A5 genotype distribution was: *1/*1 = 2 (4.3%), *1/*3 = 5 (10.9%), *3/*3 = 18 (39.1%), *3/*7 = 2 (4.3%), and Unknown = 19 (41.3%); the *3/*7 genotype was pooled with the *3/*3 nonexpresser category. The published final-model equation (Section 3.2) pools the Unknown stratum with *3/*3 nonexpressers under a single CL/F multiplier; this implementation follows that grouping (CYP3A5_EXPR = 0 means *3/*3 or unknown). Table 2 reports the final-model CL/F multipliers as 1.04 for the *3/*3-or-unknown reference stratum and 1.98 for the *1-carrier expresser stratum; the published equation rounds these to 1.0 and 2.0 respectively. This model file uses the Table 2 estimates (1.04 and 1.98) as the canonical reproducible parameter values. CL/F enters as f_CYP3A5 = e_cyp3a5_nonexpr_cl + (e_cyp3a5_expr_cl - e_cyp3a5_nonexpr_cl) * CYP3A5_EXPR, so expressers have approximately 90% higher apparent oral clearance than *3/*3-or-unknown subjects of the same weight, eGFR, hematocrit, and donor source.",
      source_name        = "CYP3A5"
    ),
    DONOR_DECEASED = list(
      description        = "Donor-source indicator: 1 if the recipient received the transplanted kidney from a deceased (post-mortal) donor; 0 if from a living donor.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (living-donor graft)",
      notes              = "Time-fixed per subject (assigned at transplantation date). Andrews 2017 Table 1 records 36 living-donor (78.3%) and 10 deceased-donor (21.7%) recipients in the model-building cohort. The published equation parameterises the donor effect as a 0.74 multiplier on CL/F for living-donor recipients (deceased-donor recipients are the reference at multiplier 1.0; Andrews 2017 Section 3.4 reports deceased-donor recipients having ~35% higher CL/F than living-donor recipients). Stored under the canonical DONOR_DECEASED column per inst/references/covariate-columns.md (1 = deceased-donor cohort, the smaller and higher-CL/F subgroup; 0 = living-donor cohort, the more common reference subgroup) so the value orientation matches the broader 'index condition = 1' register convention. CL/F enters as f_donor = e_donor_living_cl + (1 - e_donor_living_cl) * DONOR_DECEASED, equivalent to f_donor = 0.74 when DONOR_DECEASED = 0 (living donor) and f_donor = 1.0 when DONOR_DECEASED = 1 (deceased donor).",
      source_name        = "donor"
    ),
    IMMUNOASSAY = list(
      description        = "Per-sample bioanalytical assay indicator: 1 if the tacrolimus concentration was measured by a non-radioactive immunoassay (microparticle enzyme immunoassay MEIA or chemiluminescence immunoassay; pre-LC-MS/MS era); 0 if measured by LC-MS/MS.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (LC-MS/MS reference method)",
      notes              = "Per-sample (per-row) time-varying indicator: Andrews 2017 Section 2.3 reports that 9% of the model-building dataset's tacrolimus concentrations were measured before the laboratory's introduction of LC-MS/MS using an immunoassay (LLOQ 1.5 ng/mL), and the remaining 91% by a validated LC-MS/MS method (LLOQ 1.0 ng/mL). The assay-method indicator switches both the additive and the proportional residual-error magnitudes between the two analytical methods. The immunoassay arm uses additive SD 1.01 ng/mL and proportional SD 0.13 (final-model Table 2); the LC-MS/MS arm uses additive SD 0.28 ng/mL and proportional SD 0.21. In a pure-LC-MS/MS prospective dataset, set IMMUNOASSAY = 0 for every row.",
      source_name        = "assay"
    )
  )

  population <- list(
    species               = "human",
    n_subjects            = 46L,
    n_studies             = 1L,
    age_range             = "2.4-17.9 years",
    age_median            = "9.1 years",
    weight_range          = "11.6-83.7 kg",
    weight_median         = "28.4 kg",
    sex_female_pct        = 43.5,
    race_ethnicity        = c(Caucasian = 73.9, Black = 13.0, Asian = 4.3, Other = 8.7),
    disease_state         = "Paediatric kidney transplant recipients during the first 6 weeks post-transplantation. All children received twice-daily oral immediate-release tacrolimus (Prograft capsules, Modigraf granules for suspension, or extemporaneously prepared suspension) as part of the TWIST immunosuppressive protocol (basiliximab + tacrolimus + mycophenolic acid + a 5-day course of glucocorticoids). Initial tacrolimus dose was 0.3 mg/kg/day divided into two doses every 12 h; subsequent doses were adjusted by therapeutic drug monitoring to the institutional target pre-dose concentration (10-15 ng/mL during the first 3 weeks, then 7-12 ng/mL).",
    dose_range            = "0.3 mg/kg/day starting dose (twice-daily, every 12 h); subsequent doses titrated by therapeutic drug monitoring to a target pre-dose concentration of 10-15 ng/mL (first 3 weeks) then 7-12 ng/mL. Median initial daily dose 0.29 mg (range 0.23-0.39 mg, Table 1).",
    regions               = "Netherlands (model-building cohort: Erasmus Medical Center, Sophia Children's Hospital, Rotterdam; external validation cohort: Radboud University Medical Center, Amalia Children's Hospital, Nijmegen).",
    cyp3a5_distribution   = "Model-building cohort (n = 46): *1/*1 = 2 (4.3%), *1/*3 = 5 (10.9%), *3/*3 = 18 (39.1%), *3/*7 = 2 (4.3%), Unknown = 19 (41.3%). The 19 Unknown-genotype subjects were retrospectively pooled with the *3/*3 nonexpresser group in the final-model equation.",
    cyp3a4_distribution   = "Model-building cohort: *1/*1 = 10 (21.7%), *1/*1G = 4 (8.7%), *1G/*1G = 2 (4.3%), Unknown = 30 (65.2%). CYP3A4 genotype was tested as a covariate on V2/F (univariate OFV drop 5.0) but was not retained after backward elimination (Table 3).",
    donor_distribution    = "Model-building cohort: 36 living-donor recipients (78.3%) and 10 deceased-donor recipients (21.7%). Deceased-donor recipients had significantly higher tacrolimus CL/F than living-donor recipients (final-model multiplier 0.74 for the living-donor cohort vs the deceased-donor reference).",
    sampling_window       = "Tacrolimus concentrations collected during the first 6 weeks post-transplantation, total 722 blood samples (median 16 per subject, range 9-23). Sampling included routine pre-dose concentrations (C0) starting approximately 3 days post-transplantation, plus one abbreviated 4-h concentration-vs-time profile at approximately 2 weeks post-transplantation (samples pre-dose and at 10, 30, 90, 120, 240 min post-dose).",
    assay                 = "Whole-blood tacrolimus measured by validated LC-MS/MS (91% of model-building samples, LLOQ 1.0 ng/mL) or by immunoassay (9% of samples, pre-LC-MS/MS-introduction at the laboratory; LLOQ 1.5 ng/mL). The assay-method indicator (IMMUNOASSAY) was retained in the final residual-error model with separate additive + proportional magnitudes for each assay.",
    iov_structure         = "Inter-occasion variability (IOV) was identified on CL/F (18% CV) and V2/F (35% CV) in addition to the diagonal IIV (Table 2 final-model). This model file does NOT encode IOV structurally -- the source paper does not define an operational occasion column for the model-library use case, and the nlmixr2lib convention (Brooks 2021 precedent) is to omit IOV when no occasion mapping is defined; see vignette Assumptions and deviations.",
    external_validation   = "23 children from the Radboud University Medical Center (transplanted March 2012 - July 2015, selected with the same inclusion criteria) were used for external validation of the final model. The external VPC confirmed adequate median and variability prediction; the cohort contained insufficient CYP3A5 expressers to validate that subgroup.",
    notes                 = "Retrospective analysis of paediatric kidney transplant recipients transplanted November 2009 - April 2016 (model-building cohort). Routinely-collected clinical and demographic data abstracted from medical records. The model is intended for the FIRST 6 WEEKS post-transplantation in IMMEDIATE-RELEASE oral tacrolimus dosing in children aged 2-18 years; it is NOT validated for the once-daily extended-release formulation, for adults, or for the stable maintenance phase (>= 6 weeks post-transplantation). Andrews 2017 Section 4 / Discussion explicitly limits applicability."
  )

  ini({
    # Final-model fixed-effect parameter estimates from Andrews 2017 Table 2.
    # Reference subject for the typical-value structural parameters: WT = 70 kg
    # paediatric kidney transplant recipient, CYP3A5 *3/*3-or-unknown (multiplier
    # 1.04, the Table 2 final-model value), deceased-donor graft (multiplier
    # 1.0, the equation reference), eGFR 69 mL/min/1.73 m^2 (cohort median),
    # HCT >= 0.3 L/L (i.e., the piecewise hematocrit effect is inactive at the
    # reference). All apparent clearances (CL/F, Q/F) are in L/h; apparent
    # volumes (V1/F, V2/F) in L; ka in 1/h; tlag in h.
    lka   <- log(0.56) ; label("Absorption rate constant ka (1/h)")                                                         # Andrews 2017 Table 2 final ka = 0.56 (table column header reads 'L/h' which is a typo; the parameter is a first-order absorption rate constant with units 1/h)
    ltlag <- log(0.37) ; label("Absorption lag time tlag (h)")                                                              # Andrews 2017 Table 2 final tlag = 0.37 h
    lcl   <- log(50.5) ; label("Apparent oral clearance CL/F at the reference subject (L/h)")                               # Andrews 2017 Table 2 final CL/F = 50.5 L/h (at the reference covariate set defined above)
    lvc   <- log(206)  ; label("Apparent central volume V1/F at WT = 70 kg (L)")                                            # Andrews 2017 Table 2 final V1/F = 206 L
    lq    <- log(114)  ; label("Apparent inter-compartmental clearance Q/F at WT = 70 kg (L/h)")                            # Andrews 2017 Table 2 final Q/F = 114 L/h
    lvp   <- log(1520) ; label("Apparent peripheral volume V2/F at WT = 70 kg (L)")                                         # Andrews 2017 Table 2 final V2/F = 1520 L

    # Allometric exponents (Andrews 2017 Section 3.1: "Allometric scaling with
    # fixed exponents [0.75 (CL/F and Q/F) and 1 (V1/F and V2/F)] significantly
    # improved the model. Estimation of the exponents did not improve the model
    # and resulted in values not significantly different than the fixed
    # values."). Theory-based and held fixed during estimation. ka is NOT
    # allometrically scaled in Andrews 2017 (no exponent reported on ka,
    # unlike Prytula 2016 which uses -0.25).
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent of (WT/70) on CL/F (unitless; fixed at theory value)")              # Andrews 2017 Section 3.1
    e_wt_q  <- fixed(0.75) ; label("Allometric exponent of (WT/70) on Q/F (unitless; fixed at theory value)")               # Andrews 2017 Section 3.1
    e_wt_vc <- fixed(1)    ; label("Allometric exponent of (WT/70) on V1/F (unitless; fixed at theory value)")              # Andrews 2017 Section 3.1
    e_wt_vp <- fixed(1)    ; label("Allometric exponent of (WT/70) on V2/F (unitless; fixed at theory value)")              # Andrews 2017 Section 3.1

    # Covariate effects on CL/F. Andrews 2017 final-model equation (Section 3.2,
    # reproduced from the PDF layout):
    #   CL/F = 50.5 * (WT/70)^0.75 *
    #          [(1.0 if CYP3A5 *3/*3 or unknown) or (2.0 if CYP3A5 *1/*3 or *1/*1)] *
    #          (0.74 if living donor; else 1.0 deceased reference) *
    #          (eGFR/69)^0.19 *
    #          (HCT/0.3)^(-0.44) if HCT < 0.3 L/L, else 1.0
    # Table 2 reports the CYP3A5 *3/*3 multiplier as 1.04 (not 1.0) and the
    # CYP3A5 *1-carrier multiplier as 1.98 (not 2.0); the published equation
    # rounds these. This model file uses the Table 2 final-model estimates as
    # the canonical reproducible values. The "Unknown" stratum is pooled with
    # *3/*3 per the published equation.
    e_cyp3a5_nonexpr_cl <- 1.04  ; label("CYP3A5 *3/*3 or unknown multiplier on CL/F (vs the pooled-population reference)")  # Andrews 2017 Table 2 final CYP3A5 *3/*3 = 1.04
    e_cyp3a5_expr_cl    <- 1.98  ; label("CYP3A5 *1/*1 or *1/*3 expresser multiplier on CL/F")                               # Andrews 2017 Table 2 final CYP3A5 *1/*1 or *1/*3 = 1.98
    e_donor_living_cl   <- 0.74  ; label("Living-donor multiplier on CL/F (deceased-donor cohort as reference)")             # Andrews 2017 Table 2 final Donor living = 0.74
    e_egfr_cl           <- 0.19  ; label("eGFR power exponent on CL/F, centred at 69 mL/min/1.73 m^2")                       # Andrews 2017 Table 2 final eGFR exponent = 0.19
    e_hct_cl            <- -0.44 ; label("Hematocrit power exponent on CL/F, centred at 0.3 L/L; effective only when HCT < 0.3 L/L") # Andrews 2017 Table 2 final Hematocrit < 0.3 (L/L) = -0.44

    # Diagonal inter-individual variability on ka, CL/F, V1/F, V2/F. Andrews
    # 2017 Table 2 reports IIV as %CV with no inter-eta correlations
    # documented (no correlation matrix or off-diagonal entries reported in
    # the source). Variances on the internal log-scale are computed as
    # omega^2 = log(1 + CV^2):
    #   ka   CV 188% -> log(1 + 1.88^2) = log(4.5344) = 1.51178
    #   CL/F CV  25% -> log(1 + 0.25^2) = log(1.0625) = 0.06062
    #   V1/F CV  69% -> log(1 + 0.69^2) = log(1.4761) = 0.38944
    #   V2/F CV  62% -> log(1 + 0.62^2) = log(1.3844) = 0.32509
    etalka ~ 1.51178   # Andrews 2017 Table 2 IIV ka 188% CV
    etalcl ~ 0.06062   # Andrews 2017 Table 2 IIV CL/F 25% CV
    etalvc ~ 0.38944   # Andrews 2017 Table 2 IIV V1/F 69% CV
    etalvp ~ 0.32509   # Andrews 2017 Table 2 IIV V2/F 62% CV

    # Combined additive + proportional residual error with separate magnitudes
    # for the immunoassay (9% of model-building samples; LLOQ 1.5 ng/mL) and
    # LC-MS/MS (91% of samples; LLOQ 1.0 ng/mL) assays. Andrews 2017 Table 2
    # final-model values:
    #   Additive (ng/mL):     immunoassay 1.01;  LC-MS/MS 0.28
    #   Proportional (frac):  immunoassay 0.13;  LC-MS/MS 0.21
    # Switched per-row inside model() by the IMMUNOASSAY indicator.
    addSd_immuno <- 1.01 ; label("Additive residual SD for immunoassay samples (ng/mL)")                                     # Andrews 2017 Table 2 final immunoassay additive = 1.01 ng/mL
    addSd_lcms   <- 0.28 ; label("Additive residual SD for LC-MS/MS samples (ng/mL)")                                        # Andrews 2017 Table 2 final LC-MS/MS additive = 0.28 ng/mL
    propSd_immuno <- 0.13 ; label("Proportional residual SD for immunoassay samples (fraction)")                             # Andrews 2017 Table 2 final immunoassay proportional = 0.13
    propSd_lcms   <- 0.21 ; label("Proportional residual SD for LC-MS/MS samples (fraction)")                                # Andrews 2017 Table 2 final LC-MS/MS proportional = 0.21
  })

  model({
    # Body-weight scaling reference: 70 kg adult (Andrews 2017 Section 3.1).
    wt70 <- WT / 70

    # CYP3A5 multiplier on CL/F. Pools *3/*3 nonexpressers with unknown
    # genotype into CYP3A5_EXPR = 0, per Andrews 2017 final-model equation.
    f_cyp3a5 <- e_cyp3a5_nonexpr_cl + (e_cyp3a5_expr_cl - e_cyp3a5_nonexpr_cl) * CYP3A5_EXPR

    # Donor-source multiplier on CL/F. Living-donor recipients have a 0.74
    # multiplier relative to the deceased-donor reference (which itself has
    # multiplier 1.0). Encoded so that DONOR_DECEASED = 0 (living) gives
    # f_donor = e_donor_living_cl = 0.74, and DONOR_DECEASED = 1 (deceased)
    # gives f_donor = 1.0.
    f_donor <- e_donor_living_cl + (1 - e_donor_living_cl) * DONOR_DECEASED

    # eGFR power effect on CL/F, centred at the cohort median 69 mL/min/1.73 m^2.
    f_egfr <- (CRCL / 69) ^ e_egfr_cl

    # Piecewise hematocrit effect on CL/F. Active only when HCT < 0.3 L/L.
    # The boolean comparison (HCT < 0.3) returns 1 when active and 0
    # otherwise; raising the (HCT/0.3)^(-0.44) factor to that 0/1 indicator
    # gives the piecewise behaviour exactly: factor = 1 when HCT >= 0.3,
    # factor = (HCT/0.3)^(-0.44) when HCT < 0.3.
    hct_below <- (HCT < 0.3)
    f_hct <- ((HCT / 0.3) ^ e_hct_cl) ^ hct_below

    # Individual PK parameters with the Andrews 2017 covariate equations.
    ka   <- exp(lka + etalka)
    cl   <- exp(lcl + etalcl) * wt70 ^ e_wt_cl * f_cyp3a5 * f_donor * f_egfr * f_hct
    vc   <- exp(lvc + etalvc) * wt70 ^ e_wt_vc
    q    <- exp(lq)           * wt70 ^ e_wt_q
    vp   <- exp(lvp + etalvp) * wt70 ^ e_wt_vp
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

    # Tacrolimus whole-blood concentrations reported in ng/mL. Dose in mg, vc
    # in L, so central/vc is in mg/L = ug/mL; multiply by 1000 to convert to
    # ng/mL.
    Cc <- central / vc * 1000

    # Per-sample assay-method-conditional residual error. Following the
    # Cirincione_2017_exenatide precedent, switch between the immunoassay
    # and LC-MS/MS error magnitudes inside model() using the IMMUNOASSAY
    # binary indicator (1 = immunoassay, 0 = LC-MS/MS reference).
    addSd  <- addSd_immuno  * IMMUNOASSAY + addSd_lcms  * (1 - IMMUNOASSAY)
    propSd <- propSd_immuno * IMMUNOASSAY + propSd_lcms * (1 - IMMUNOASSAY)
    Cc ~ add(addSd) + prop(propSd)
  })
}
