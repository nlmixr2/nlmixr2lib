Martial_2017_micafungin <- function() {
  description <- "Two-compartment population PK model for IV micafungin in adult intensive-care-unit patients with suspected or proven fungal infection (Martial 2017). Body-weight allometric scaling (fixed exponents 0.75 on CL and Q, 1 on V1 and V2; 70 kg reference), log-normal IIV on CL and V1 (encoded as diagonal; the source reports a qualitative non-zero correlation but no numerical covariance), and a proportional residual error. No covariates were retained in the final model (only weight via the a-priori allometric structure)."
  reference <- paste(
    "Martial LC, ter Heine R, Schouten JA, Hunfeld NG, van Leeuwen HJ,",
    "Verweij PE, de Lange DW, Pickkers P, Bruggemann RJ; FANTASTIC",
    "Consortium. Population Pharmacokinetic Model and Pharmacokinetic",
    "Target Attainment of Micafungin in Intensive Care Unit Patients.",
    "Clin Pharmacokinet. 2017;56(10):1197-1206.",
    "doi:10.1007/s40262-017-0509-5.",
    sep = " "
  )
  vignette <- "Martial_2017_micafungin"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used for allometric scaling on CL and Q (exponent 0.75 fixed) and V1 and V2 (exponent 1 fixed) with reference weight 70 kg, per Martial 2017 Section 2.4 (Methods). Time-varying weight is permitted by the structural model; the ICU cohort had a median baseline weight of 76.5 kg (range 50-134 kg, Table 1).",
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 20L,
    n_studies      = 1L,
    age_range      = "20-84 years (median 68)",
    age_median     = "68 years",
    weight_range   = "50-134 kg (median 76.5)",
    weight_median  = "76.5 kg",
    sex_female_pct = 60,
    race_ethnicity = "Not stratified in the source (single-country Dutch ICU cohort).",
    disease_state  = "Adults admitted to the intensive care unit (ICU) and receiving micafungin for suspected or proven fungal infection. All 20 patients had pronounced hypoalbuminaemia (serum albumin <= 34 g/L); 5 (25%) received continuous veno-venous haemofiltration (CVVH) and 1 (5%) received intermittent haemodialysis at baseline (Martial 2017 Table 1).",
    dose_range     = "100 mg micafungin once daily as a ~1 h IV infusion; therapy continued as clinically indicated, capped at 14 days plus 3 follow-up days for sampling. No dose adaptations were performed during the study (Martial 2017 Sections 2.2 and 3.1).",
    bmi_range      = "16.3-47.5 kg/m^2 (median 24.6)",
    regions        = "Netherlands (multicentre FANTASTIC Consortium: Radboud University Medical Center, Canisius Wilhelmina Hospital, Erasmus Medical Center, Rijnstate Hospital, University Medical Centre Utrecht).",
    indication     = "Empiric or targeted treatment of Candida (and occasional Aspergillus) infections in the ICU; 100% of patients had Candida spp. infection.",
    n_observations = 356L,
    notes          = "Baseline demographics from Martial 2017 Table 1; n = 20 evaluable patients with a first PK curve on day 3. Hypoalbuminaemia categories: 25-34 g/L (10%), 15-24 g/L (65%), <15 g/L (25%). Renal replacement: CVVH 25%, intermittent haemodialysis 5%. The cohort was 60% female. No covariates beyond a-priori allometric weight scaling were retained in the final model; stepwise screening of albumin, CVVH, and SOFA score did not produce significant OFV drops (Section 3.2)."
  )

  ini({
    # Structural PK parameters from Martial 2017 Table 2 (final model NONMEM
    # estimates). Reported as point estimates with RSE; the reported RSE is
    # labelled "root square error" in the table caption but is the standard
    # NONMEM relative standard error from the covariance step. All four
    # parameters are allometrically scaled to a 70 kg reference subject
    # (Methods Section 2.4: "an allometric exponent of 0.75 for all flow
    # parameters and an exponent of 1 for all volume parameters standardised
    # to a 70-kg patient as proposed previously").
    lcl  <- log(1.10)  ; label("Clearance for a 70 kg subject (CL, L/h)")                                # Martial 2017 Table 2: CL = 1.10 [RSE 8%]
    lvc  <- log(17.6)  ; label("Central volume of distribution V1 for a 70 kg subject (L)")              # Martial 2017 Table 2: V1 = 17.6 [RSE 14%]
    lq   <- log(0.363) ; label("Inter-compartmental clearance Q for a 70 kg subject (L/h)")              # Martial 2017 Table 2: Q  = 0.363 [RSE 20%]
    lvp  <- log(3.63)  ; label("Peripheral volume of distribution V2 for a 70 kg subject (L)")           # Martial 2017 Table 2: V2 = 3.63 [RSE 8%]

    # Allometric exponents -- fixed a priori per Methods Section 2.4 and the
    # Discussion (page 1203): "Based on our data (n = 20 patients), the true
    # exponent is not identifiable and using an empirical estimate would not
    # allow for extrapolation beyond our dataset."
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent on CL (unitless; fixed)")                        # Martial 2017 Section 2.4 (Methods)
    e_wt_vc <- fixed(1)    ; label("Allometric exponent on V1 (unitless; fixed)")                        # Martial 2017 Section 2.4 (Methods)
    e_wt_q  <- fixed(0.75) ; label("Allometric exponent on Q (unitless; fixed)")                         # Martial 2017 Section 2.4 (Methods)
    e_wt_vp <- fixed(1)    ; label("Allometric exponent on V2 (unitless; fixed)")                        # Martial 2017 Section 2.4 (Methods)

    # Inter-individual variability (log-normal). The paper reports IIV as
    # %CV; the internal variance is omega^2 = log(CV^2 + 1). IIV was retained
    # on CL and V1 only (Martial 2017 Section 3.2: "data did not support the
    # addition of IIV on Q or V 2"). The text also notes that "Allowing a
    # correlation between IIV of CL and V 1 further improved the model
    # (difference OFV = 13.98)", but the numerical covariance / correlation
    # is not reported anywhere in the manuscript. Per nlmixr2lib precedent
    # for unreported correlations (Narwal_2013_sifalimumab,
    # Mould_2007_alemtuzumab, Chakraborty_2012_canakinumab,
    # Wang_2020_ontamalimab), the IIV is encoded as diagonal here and the
    # missing correlation is documented in the vignette Errata.
    etalcl ~ 0.14906   # Martial 2017 Table 2: IIV CL  40.1% CV (shrinkage 1.2%); omega^2 = log(1 + 0.401^2)
    etalvc ~ 0.42891   # Martial 2017 Table 2: IIV V1  73.2% CV (shrinkage 41%);   omega^2 = log(1 + 0.732^2)

    # Inter-occasion variability (IOV) on V1 was reported at 37.0% CV
    # (Martial 2017 Section 3.2 and Table 2 -- the Table 2 row label reads
    # "IOV V2" but the text describes IOV on V1 in five separate passages
    # and notes V2 carries no IIV either, so the Table 2 label is treated
    # as a typographical slip). The IOV is NOT encoded in this model
    # because the per-occasion structure (number and definition of
    # occasions) is not disclosed in the paper or any on-disk supplement
    # and the NONMEM control stream is not on disk; following the
    # Archary_2019_lamivudine precedent the IOV term is documented as a
    # deviation in the vignette Errata rather than encoded with a guessed
    # OCC partitioning.

    # Residual error -- proportional error model fit best per Methods
    # Section 2.4 (additive and combined additive + proportional were
    # also evaluated but not retained).
    propSd <- 0.17 ; label("Proportional residual error (fraction)")                                     # Martial 2017 Table 2: proportional residual error = 17% [RSE 24%]
  })

  model({
    # Individual PK parameters with body-weight allometric scaling
    # (reference 70 kg; fixed exponents per Methods Section 2.4). The final
    # model retained no other covariates -- albumin, CVVH, and SOFA score
    # were screened in stepwise forward addition / backward elimination
    # but none produced an OFV drop greater than 3.84 (Section 3.2).
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q  <- exp(lq)           * (WT / 70)^e_wt_q
    vp <- exp(lvp)          * (WT / 70)^e_wt_vp

    # Micro-constants for the 2-compartment IV ODE.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 2-compartment IV PK with first-order elimination from the central
    # compartment. Micafungin is given as a ~1 h IV infusion that is
    # encoded on the dose record (rate or dur) rather than as a depot
    # absorption process.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Plasma micafungin concentration. Dose in mg, Vc in L -> Cc in mg/L
    # (equivalent to ug/mL = mcg/mL, the units reported in the paper).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
