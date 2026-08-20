Xie_2025_midazolam <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for midazolam (MDZ) and its",
    "active metabolite 1-hydroxymidazolam (1-OH-MDZ) during continuous",
    "intravenous infusion in 61 mechanically ventilated Chinese adult ICU",
    "patients (Xie 2025). Midazolam: one-compartment IV disposition.",
    "1-OH-MDZ: a separate one-compartment model for formation and",
    "elimination, driven by the fraction of midazolam converted to the",
    "metabolite (FMET fixed at 0.6 from Seng et al); the metabolite CL and V",
    "are therefore apparent values conditional on that fixed fraction.",
    "Midazolam clearance decreases with aspartate aminotransferase (power",
    "form, reference 37 IU/L) and is 40.5% lower in carriers of the",
    "homozygous-variant NR1I2 (PXR) rs2461817 genotype (fractional-shift",
    "form). 1-OH-MDZ clearance increases with total body weight (power form,",
    "reference 62 kg). IIV on both clearances; midazolam residual error is",
    "combined proportional plus additive, 1-OH-MDZ residual error is",
    "additive only."
  )
  reference <- paste(
    "Xie H, Zheng Y, Zhang H, Guo Y, Liu M, Weng Q, Wu X.",
    "Association of NR1I2 Polymorphism with Midazolam Clearance in",
    "Mechanically Ventilated ICU Patients: A Population Pharmacokinetic and",
    "Pharmacogenetic Study.",
    "Drug Des Devel Ther. 2025;19:1527-1541.",
    "doi:10.2147/DDDT.S495647.",
    sep = " "
  )
  vignette <- "Xie_2025_midazolam"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Both states are verified against Xie 2025 Figure 2
  # (the published structural diagram) and the Blood Sampling section
  # (plasma, LC-MS/MS).
  compartmentData <- list(
    central      = list(analyte = "midazolam", units = "mg", specimen = "plasma", verified = TRUE),
    central_1ohm = list(analyte = "1-hydroxymidazolam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    AST = list(
      description        = "Serum aspartate aminotransferase activity",
      units              = "IU/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline hepatic-injury marker. Enters midazolam clearance as a",
        "power term (AST / 37)^-0.397, normalised to the population median",
        "of 37 IU/L (Xie 2025 Table 1: AST median 37.0, range 10.0-362.0).",
        "The paper does not print the covariate equation; the power form and",
        "the 37 IU/L reference are recovered exactly from the paper's own",
        "reported simulation clearances (Results 'Simulations': AST 22 ->",
        "27.8 L/h, AST 60 -> 18.7 L/h against the typical 22.6 L/h at the",
        "median) -- 22.6 * (22/37)^-0.397 = 27.78 and 22.6 * (60/37)^-0.397 =",
        "18.65. AST was collected once per patient from the electronic",
        "medical record and is treated here as time-fixed; the paper does not",
        "state that it was updated during the infusion. Units are reported by",
        "the paper as IU/L, used interchangeably with the register's",
        "canonical U/L."
      ),
      source_name        = "AST"
    ),
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters 1-OH-MDZ apparent clearance as a power term (WT / 62)^1.58,",
        "normalised to the population median of 62 kg (Xie 2025 Table 1:",
        "weight median 62.0 kg, range 38.0-87.6). The exponent is estimated",
        "(1.58, RSE 30%), not fixed at an allometric 0.75. The paper does not",
        "print the covariate equation; the power form and the 62 kg reference",
        "are recovered exactly from Results 'Simulations' (BW 54 kg -> 53.9",
        "L/h, BW 65 kg -> 72.3 L/h against the typical 67.1 L/h at the",
        "median) -- 67.1 * (54/62)^1.58 = 53.94 and 67.1 * (65/62)^1.58 =",
        "72.30. Body weight was NOT retained on midazolam clearance or on",
        "either volume; the Discussion argues total body weight tracks the",
        "metabolic (UGT-mediated) capacity that clears 1-OH-MDZ."
      ),
      source_name        = "BW"
    ),
    SNP_NR1I2_RS2461817_HOM = list(
      description        = paste(
        "Binary indicator for the homozygous-variant genotype of the NR1I2",
        "(PXR) rs2461817 single-nucleotide polymorphism"
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (wild-type homozygote or heterozygote)",
      notes              = paste(
        "Recessive-model encoding: 1 = homozygous mutant, 0 = the union of",
        "wild-type homozygotes and mutant heterozygotes. This pooling is the",
        "paper's own, stated verbatim in the Figure 6 caption: 'GENE = 1",
        "indicates that the NR1I2 rs2461817 genotype reflects a homozygous",
        "mutation. GENE = 0 indicates that the NR1I2 rs2461817 genotype is",
        "wild-type homozygous or mutant heterozygous.' Time-fixed per",
        "subject (germline genotype), determined by genotyping 23 loci across",
        "CYP3A4, CYP3A5, ABCB1 and NR1I2 (Methods 'Genotyping'). Enters",
        "midazolam clearance as a fractional shift (1 - 0.405 * GENE), not as",
        "an exponential factor: the paper's Results 'Simulations' report CL",
        "falling from 22.6 to 13.4 L/h for homozygous mutants, and",
        "22.6 * (1 - 0.405) = 13.45 whereas 22.6 * exp(-0.405) = 15.07."
      ),
      source_name        = "GENE (NR1I2 rs2461817)"
    )
  )

  # Covariates that Xie 2025 screened but did NOT retain in the final model.
  # Documentation only -- none is referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened in the stepwise covariate search (Methods 'Covariate Analysis'); not retained. Table 1 median 67 years (range 29-90)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened; not retained. Table 1: 42 male (68.9%) / 19 female (31.1%)."
    ),
    APACHE_II = list(
      description = "Acute Physiology and Chronic Health Evaluation II score",
      units       = "points",
      type        = "continuous",
      notes       = "Screened; not retained. Table 1 median 16 (range 2-24). The Discussion explicitly notes 'we did not observe a significant correlation between MDZ CL and APACHE II scores', in contrast to Swart 2004 where APACHE II >= 26 increased Vd."
    ),
    ALB = list(
      description = "Serum albumin concentration",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened; not retained. Table 1 median 31.7 g/L (range 21.7-65.1). The Discussion notes most patients had albumin in the normal range, unlike Franken 2017 / Vree where hypoalbuminaemia drove clearance."
    ),
    ALT = list(
      description = "Serum alanine aminotransferase activity",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Screened; not retained (AST was retained instead). Table 1 median 22.0 IU/L (range 2.0-391.0)."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened; not retained. Table 1 median 15.4 umol/L (range 3.5-169.7)."
    ),
    CRP = list(
      description = "C-reactive protein",
      units       = "mg/L",
      type        = "continuous",
      notes       = "Screened; not retained. Table 1 median 82.14 mg/L (range 1.95-295.74)."
    ),
    CRCL = list(
      description = "Creatinine clearance estimated by the Cockcroft-Gault equation",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened; not retained. Table 1 median 70.77 mL/min (range 8.97-203.4)."
    ),
    CONMED_PROPOFOL = list(
      description = "Concomitant propofol administration indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened (usage rate > 10% per Methods 'Covariate Analysis'); not retained. Table 1: 22 patients (36.1%). The Discussion notes the absence of a propofol effect on midazolam Vd despite a mechanistic expectation."
    ),
    CONMED_METHYLPREDNISOLONE = list(
      description = "Concomitant methylprednisolone administration indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened (usage rate > 10%); not retained. Table 1: 7 patients (11.5%)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 61L,
    n_studies      = 1L,
    age_range      = "29-90 years",
    age_median     = "67 years",
    weight_range   = "38.0-87.6 kg",
    weight_median  = "62.0 kg",
    sex_female_pct = 31.1,
    race_ethnicity = "Chinese (single-centre cohort, Fuzhou, Fujian Province)",
    disease_state  = paste(
      "Mechanically ventilated adult ICU patients requiring at least 24 h of",
      "mechanical ventilation and receiving continuous intravenous midazolam",
      "for sedation. Mixed admission diagnoses, mostly postoperative patients",
      "in relatively good general health (Discussion). Baseline severity was",
      "moderate: APACHE II median 16 (range 2-24). Exclusions were hepatic",
      "coma or cirrhosis, neurological inability to assess sedation,",
      "haemodynamic instability requiring frequent dose changes, and",
      "pregnancy / lactation / midazolam allergy."
    ),
    dose_range     = paste(
      "Continuous intravenous infusion of midazolam 1 mg/mL (50 mg diluted to",
      "40 mL of 0.9% saline or 5% glucose), pumped at 2-4 mL/h at initiation",
      "and titrated to the target Richmond Agitation-Sedation Scale score.",
      "Observed initial infusion rates ranged 2-6 mg/h. Observation window",
      "0-24 h."
    ),
    regions        = "China (Fujian Medical University Union Hospital, Fuzhou)",
    notes          = paste(
      "Prospective observational study, April 2020 - August 2022 (IRB",
      "2021YF003-01), reported per STROBE. 69 patients screened, 8 excluded",
      "(3 incomplete data, 5 not meeting inclusion criteria). 237 paired",
      "midazolam / 1-OH-MDZ plasma samples analysed; 3 midazolam and 3",
      "1-OH-MDZ values below the limit of quantitation (1.3% each) were",
      "excluded. Arterial sampling at t = 0 (pre-dose) and windows 0-0.5,",
      "1-3, 4-6 and 10-12 h. LC-MS/MS with midazolam-d4 internal standard;",
      "linear ranges 0.5-1000 ng/mL (midazolam) and 0.25-500 ng/mL",
      "(1-OH-MDZ); LLOQ 0.5 and 0.25 ng/mL. Observed plasma concentrations:",
      "midazolam median 149.26 ng/mL (range 1.53-1622.36), 1-OH-MDZ median",
      "18.66 ng/mL (range 0.28-166.19). Missing covariates imputed at the",
      "population median. Fitted in NONMEM 7.5.0 with PsN 5.0.0 and R 3.6.1;",
      "evaluated by goodness-of-fit plots, a 1000-replicate non-parametric",
      "bootstrap and a 1000-replicate VPC. All 22 evaluable SNP genotypes",
      "except NR1I2 rs1464603 satisfied Hardy-Weinberg equilibrium",
      "(Supplementary Table 1, not on disk)."
    )
  )

  ini({
    # ==================================================================
    # MIDAZOLAM DISPOSITION -- Xie 2025 Table 2, "Final Model" column.
    # One-compartment IV (Figure 2). The clearance value is the typical
    # value at the reference covariate state AST = 37 IU/L and
    # SNP_NR1I2_RS2461817_HOM = 0.
    # ==================================================================
    lcl <- log(22.6)
    label("Midazolam clearance at AST = 37 IU/L, non-homozygous rs2461817 (L/h)")   # Table 2 Final Model: CL MDZ = 22.6 L/h (RSE 16%; bootstrap median 22.15, 95% CI 12.55-30.99)

    lvc <- log(16.00)
    label("Midazolam volume of distribution (L)")                                   # Table 2 Final Model: V MDZ = 16.00 L (RSE 42%; bootstrap median 15.58, 95% CI 5.44-44.61)

    # ==================================================================
    # 1-OH-MIDAZOLAM DISPOSITION -- Xie 2025 Table 2, "Final Model".
    # A separate one-compartment model for formation and elimination
    # (Methods "Structural Model"; Figure 2). Both values are APPARENT:
    # Discussion "Many restrictions apply to this study" states "the
    # reported Vd and CL of 1-OH-MDZ are 'apparent' values because we
    # could not determine the proportion of MDZ metabolized to
    # 1-OH-MDZ", so they are conditional on FMET being fixed at 0.6.
    # The clearance is the typical value at the reference WT = 62 kg.
    # ==================================================================
    lcl_1ohm <- log(67.10)
    label("1-OH-midazolam apparent clearance at WT = 62 kg (L/h)")                  # Table 2 Final Model: CL 1-OH-MDZ = 67.10 L/h (RSE 10%; bootstrap median 65.79, 95% CI 52.43-79.36)

    lvc_1ohm <- log(86.20)
    label("1-OH-midazolam apparent volume of distribution (L)")                     # Table 2 Final Model: V 1-OH-MDZ = 86.20 L (RSE 35%; bootstrap median 86.32, 95% CI 1.39-149.375)

    # ==================================================================
    # FRACTION OF MIDAZOLAM CONVERTED TO 1-OH-MIDAZOLAM.
    # Fixed, not estimated: Table 2 reports F MET = 0.6 with RSE "/" in
    # both the base and final model columns, and Methods "Structural
    # Model" states "the conversion rate of MDZ to 1-OH-MDZ was fixed at
    # 0.6 based on prior literature" (Seng et al, reference 23).
    # This is a MOLAR fraction -- see the mw_mdz / mw_1ohm constants in
    # model() for the conversion into this file's mass units.
    # ==================================================================
    fm_1ohm <- fixed(0.6)
    label("Molar fraction of midazolam converted to 1-OH-midazolam (unitless)")     # Table 2 Final Model: F MET = 0.6 (fixed, RSE reported as "/"); Methods 'Structural Model and Choice of Statistical Model'

    # ==================================================================
    # COVARIATE EFFECTS -- Xie 2025 Table 2, "Final Model" column.
    #
    # Xie 2025 does not print the covariate equations anywhere in the
    # article, and no supplement carrying them is on disk. All three
    # functional forms and both reference values were recovered by
    # back-solving the paper's own reported simulation clearances
    # (Results "Simulations") and are exact to the reported precision:
    #
    #   cl      = 22.6 * (AST / 37)^-0.397 * (1 - 0.405 * GENE)
    #     AST 22, GENE 0 -> 27.78 L/h   (paper: 27.8)
    #     AST 37, GENE 0 -> 22.60 L/h   (paper: 22.6)
    #     AST 60, GENE 0 -> 18.65 L/h   (paper: 18.7)
    #     AST 37, GENE 1 -> 13.45 L/h   (paper: 13.4)
    #
    #   cl_1ohm = 67.1 * (WT / 62)^1.58
    #     WT 54 -> 53.94 L/h            (paper: 53.9)
    #     WT 62 -> 67.10 L/h            (paper: 67.1)
    #     WT 65 -> 72.30 L/h            (paper: 72.3)
    #
    # The exponential alternative for the genotype effect is excluded:
    # 22.6 * exp(-0.405) = 15.07 L/h, not the reported 13.4 L/h.
    # ==================================================================
    e_ast_cl <- -0.397
    label("AST power exponent on midazolam CL (unitless)")                          # Table 2 Final Model: 'ASTon CL MDZ' = -0.397 (RSE 47%; bootstrap median -0.385, 95% CI -0.783 to -0.0151)

    e_nr1i2_hom_cl <- -0.405
    label("Fractional shift in midazolam CL for rs2461817 homozygous mutants (unitless)")  # Table 2 Final Model: 'rs2461817 on CL MDZ' = -0.405 (RSE 29%; bootstrap median -0.389, 95% CI -0.591 to -0.0659)

    e_wt_cl_1ohm <- 1.58
    label("Body-weight power exponent on 1-OH-midazolam CL (unitless)")             # Table 2 Final Model: 'BWon CL 1-OH-MDZ' = 1.58 (RSE 30%; bootstrap median 1.54, 95% CI 0.52-2.58)

    # ==================================================================
    # INTER-INDIVIDUAL VARIABILITY.
    # Methods "Structural Model": "Individual variability (IIV) was
    # described using an exponential model". Table 2 reports IIV as CV%
    # and its footnote gives the conversion explicitly:
    # "CV% = sqrt(exp(OMEGA) - 1) x 100", i.e. OMEGA = log(1 + CV^2).
    # Only the two clearances carry IIV in the final model; neither
    # volume does.
    # ==================================================================
    etalcl      ~ 0.302203   # Table 2 Final Model IIV CL MDZ = 59.40 % CV (RSE 20%): omega^2 = log(1 + 0.594^2) = 0.302203
    etalcl_1ohm ~ 0.287379   # Table 2 Final Model IIV CL 1-OH-MDZ = 57.70 % CV (RSE 11%): omega^2 = log(1 + 0.577^2) = 0.287379

    # ==================================================================
    # RESIDUAL VARIABILITY -- Xie 2025 Table 2, "Final Model" column.
    # Results "Development of the Population Pharmacokinetic Model":
    # "for MDZ and 1-OH-MDZ, intra-individual variability was
    # characterized using a combination model and an additive model,
    # respectively."
    #
    # Both additive terms are tabulated in ng/mL, and that reading is
    # confirmed dimensionally: the model was fitted in molar units
    # (Figure 4 plots nmol/mL) where the entire observed midazolam range
    # tops out near 5 nmol/mL and 1-OH-MDZ near 0.49 nmol/mL, so
    # additive SDs of 3.43 and 0.6 could not possibly be nmol/mL. In
    # ng/mL they are 6.9x and 2.4x the respective LLOQs (0.5 and
    # 0.25 ng/mL), which is the expected magnitude.
    # ==================================================================
    propSd <- 0.400
    label("Midazolam proportional residual SD on Cc (fraction)")                    # Table 2 Final Model: Proportional (CV%) = 40.00 (RSE 38%; bootstrap median 41.14%, 95% CI 13.47-52.69%)

    addSd <- 3.43
    label("Midazolam additive residual SD on Cc (ng/mL)")                           # Table 2 Final Model: Additive (ng/mL) = 3.43 (RSE 56%; bootstrap median 3.22, 95% CI 1.15-5.02)

    addSd_1ohm <- 0.6
    label("1-OH-midazolam additive residual SD on Cc_1ohm (ng/mL)")                 # Table 2 Final Model: 1-OH-MDZ Additive (ng/mL) = 0.6 (RSE 22%; bootstrap median 0.58, 95% CI 0.41-0.76)
  })

  model({
    # ------------------------------------------------------------------
    # Reference covariate values (Xie 2025 Table 1 population medians;
    # the same values the Simulation section stratifies around).
    # ------------------------------------------------------------------
    ast_ref <- 37   # Table 1: AST median 37.0 IU/L; Simulation uses the 25th/50th/75th percentiles 22 / 37 / 60
    wt_ref  <- 62   # Table 1: weight median 62.0 kg; Simulation uses 54 / 62 / 65 kg

    # ------------------------------------------------------------------
    # Molecular weights, Methods "Structural Model": "The dose, infusion
    # rate, and plasma concentrations of MDZ (Mw = 325.77 g/mol) and
    # 1-OH-MDZ (Mw = 341.77 g/mol) were converted to molar equivalents".
    # FMET is therefore a MOLAR fraction. This file carries amounts in mg
    # rather than umol, so the metabolite formation flux is scaled by the
    # molecular-weight ratio: 1 mole of midazolam cleared through the
    # metabolic route yields FMET moles of 1-OH-MDZ, i.e.
    # FMET * (341.77 / 325.77) = 0.6294 mg of 1-OH-MDZ per mg of
    # midazolam. Retaining this 4.9% correction is what reproduces the
    # published Figure 7 plateaus (37.5 ng/mL at 62 kg); dropping it
    # would give 35.8 ng/mL.
    # ------------------------------------------------------------------
    mw_mdz  <- 325.77
    mw_1ohm <- 341.77

    # ------------------------------------------------------------------
    # Individual parameters.
    # Midazolam CL carries a power AST term and a fractional genotype
    # shift; neither volume carries IIV or a covariate.
    # ------------------------------------------------------------------
    ast_factor   <- (AST / ast_ref)^e_ast_cl
    nr1i2_factor <- 1 + e_nr1i2_hom_cl * SNP_NR1I2_RS2461817_HOM

    cl <- exp(lcl + etalcl) * ast_factor * nr1i2_factor
    vc <- exp(lvc)

    cl_1ohm <- exp(lcl_1ohm + etalcl_1ohm) * (WT / wt_ref)^e_wt_cl_1ohm
    vc_1ohm <- exp(lvc_1ohm)

    # Mass-unit formation fraction (see the molecular-weight note above).
    fm_mass <- fm_1ohm * mw_1ohm / mw_mdz

    # ------------------------------------------------------------------
    # ODE system, exactly as drawn in Xie 2025 Figure 2. Midazolam leaves
    # the central compartment at the total rate CL_MDZ / V_MDZ, split into
    # a non-metabolic arm (1 - FMET) and a metabolic arm (FMET) that feeds
    # the 1-OH-MDZ compartment; 1-OH-MDZ is then eliminated at
    # CL_1-OH-MDZ / V_1-OH-MDZ. Doses are continuous IV infusions into
    # `central` supplied by the event table.
    # ------------------------------------------------------------------
    d/dt(central)      <- -cl * central / vc
    d/dt(central_1ohm) <-  fm_mass * cl * central / vc -
                           cl_1ohm * central_1ohm / vc_1ohm

    # ------------------------------------------------------------------
    # Observations. States hold mg; volumes are L, so amount / volume is
    # mg/L = ug/mL and the factor 1000 converts to the ng/mL reporting
    # units used throughout Xie 2025 (Table 1, Figure 3, Figures 5-7).
    # ------------------------------------------------------------------
    Cc      <- central      / vc      * 1000
    Cc_1ohm <- central_1ohm / vc_1ohm * 1000

    Cc      ~ add(addSd) + prop(propSd)
    Cc_1ohm ~ add(addSd_1ohm)
  })
}
