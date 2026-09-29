Mohanan_2017_fludarabine <- function() {
  description <- paste(
    "Two-compartment IV population PK model for the circulating fludarabine nucleoside F-ara-A in",
    "patients with aplastic anemia (AA) or Fanconi anemia (FA) receiving fludarabine phosphate",
    "30 mg/m^2/day as a 1-h infusion for 6 days (day -7 to day -2) of conditioning prior to",
    "allogeneic hematopoietic stem cell transplantation (n = 53, age 3-57 years; Mohanan 2017).",
    "The model is parameterised on a body-surface-area-normalised basis exactly as published:",
    "clearance in L/h/m^2 and central volume in L/m^2, both multiplied by BSA inside model() so",
    "that doses may be supplied as absolute amounts in mg. Mohanan 2017 showed the BSA-normalised",
    "structure to be strongly preferred over a non-normalised one (-2LL 312.34 vs 376.03; Table 3).",
    "Clearance carries two binary covariates in the multiplicative exp(beta * covariate) form the",
    "paper specifies: Fanconi anemia (versus aplastic anemia) and carriage of the NT5E / CD73",
    "5'-UTR polymorphism rs2295890 (GC or CC versus GG wild-type). Together they reproduce all four",
    "clearance cells printed in Table 2: 7.12 (WT, AA), 5.03 (variant, AA), 2.90 (WT, FA) and",
    "2.05 (variant, FA) L/h/m^2. Central volume and the peripheral-to-central rate constant k21 each",
    "carry an uncentred exponential age effect; the central-to-peripheral rate constant k12 carries",
    "none. Inter-individual variability is log-normal on all four disposition parameters, and the",
    "paper additionally reports inter-day variability, encoded here as inter-occasion variability",
    "over six daily dosing occasions selected by the OCC column. Residual error is proportional",
    "(19% CV). NOTE: the sign of the age coefficient on volume is transcribed as POSITIVE here",
    "while Table 2 prints a negative sign; see the model's covariateData[[AGE]]$notes and the",
    "validation vignette for the arithmetic that falsifies the printed sign.",
    "SOLVE THIS MODEL WITH useLinCmt = FALSE. rxSolve.rxUi() defaults to useLinCmt = TRUE, whose",
    "ODE-to-linCmt auto-conversion silently discards the peripheral compartment of this",
    "micro-constant-parameterised two-compartment model: the solve returns no peripheral1 state",
    "and decays mono-exponentially at kel, giving a terminal half-life of 2.58 h instead of the",
    "correct 10.5 h at the cohort-median age. There is no error or warning, and AUC is unaffected",
    "because Dose/CL is preserved, so only the shape of the curve reveals it.",
    sep = " "
  )
  reference <- paste(
    "Mohanan E, Panetta JC, Lakshmi KM, Edison ES, Korula A, Fouzia NA, Abraham A, Viswabandya A,",
    "Mathews V, George B, Srivastava A, Balasubramanian P. Population pharmacokinetics of",
    "fludarabine in patients with aplastic anemia and Fanconi anemia undergoing allogeneic",
    "hematopoietic stem cell transplantation. Bone Marrow Transplant. 2017;52(7):977-983.",
    "doi:10.1038/bmt.2017.79.",
    "Correction: Bone Marrow Transplant. 2018;53(11):1490. doi:10.1038/s41409-018-0276-4",
    "(license change from CC BY-NC-ND 4.0 to CC BY 4.0 only; no scientific content was revised).",
    sep = " "
  )
  vignette <- "Mohanan_2017_fludarabine"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(
      analyte = "fludarabine",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "fludarabine",
      units = "mg",
      specimen = "plasma",
      verified = FALSE
    )
  )

  covariateData <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Mohanan 2017 Table 1: median 1.49 m^2, range 0.56-1.9 m^2 across the 53 patients. BSA is",
        "the size descriptor of the entire published model rather than an ordinary covariate:",
        "clearance is reported in L/h/m^2 and central volume in L/m^2, and the paper states that",
        "the BSA-normalised base model was significantly better than the non-normalised one",
        "(-2 log-likelihood 312.34 vs 376.03, a drop of 54.69 units; Results 'F-araA PK' and",
        "Table 3). Doses were likewise prescribed per square metre (30 mg/m^2/day). This file",
        "multiplies the published per-m^2 clearance and volume by BSA inside model() so that the",
        "event table can carry ABSOLUTE doses in mg. That rescaling is exactly equivalent to the",
        "paper's per-m^2 formulation and leaves every concentration unchanged: with a dose of",
        "D mg/m^2 the concentration is (D * BSA) / (V_per_m2 * BSA) = D / V_per_m2, and the",
        "elimination rate constant cl/vc = CL_per_m2 / V_per_m2 is BSA-free. The rate constants",
        "k12 and k21 are not BSA-scaled because they are first-order rate constants, not flows.",
        "The paper does not state which BSA formula was used.",
        sep = " "
      ),
      source_name = "BSA"
    ),
    AGE = list(
      description = "Age at the start of conditioning",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Mohanan 2017 Table 1: median 17 years, range 3-57 years. Enters the central volume and the",
        "peripheral-to-central rate constant k21 as an UNCENTRED exponential term,",
        "theta = theta_intercept * exp(beta * AGE), so the tabulated intercepts 21.25 L/m^2 and",
        "0.14 1/h are the values at AGE = 0 and not at a reference age. The covariate model form is",
        "given in Methods 'F-araA PK and Population PK modeling' as theta = theta_Base *",
        "exp(beta * covariate).",
        "SIGN OF THE VOLUME COEFFICIENT. Mohanan 2017 Table 2 prints the volume row as",
        "'21.25 x exp(-0.013 x age)', i.e. a volume that FALLS with age, but the Results text",
        "('Influence of genetic variants on F-araA PK') states that 'the parameters V and K21",
        "INCREASED significantly with respect to age'. This file follows the text and uses a",
        "POSITIVE coefficient +0.013, because the paper's own base model falsifies the printed",
        "sign arithmetically. Under an uncentred exponential model the base-model (covariate-free)",
        "typical value satisfies log(base) = log(intercept) + beta * mean(AGE). The k21 row, whose",
        "coefficient is printed unambiguously POSITIVE, gives log(0.19 / 0.14) / 0.016 = 19.1",
        "years, a wholly plausible mean age for a cohort of median 17 and range 3-57. Applying the",
        "identical calculation to the volume row gives log(27.56 / 21.25) / 0.013 = +20.0 years",
        "with a positive coefficient and -20.0 years -- a physically impossible mean age -- with",
        "the printed negative one. Two independent rows of the same table therefore agree on a",
        "cohort mean age of about 19-20 years only when both coefficients are positive. The",
        "printed minus sign is treated as a typesetting error; it is present in the published PDF",
        "and is not an artefact of text extraction. The validation vignette reproduces both the",
        "algebraic back-solve and a cohort simulation under each sign.",
        sep = " "
      ),
      source_name = "age"
    ),
    DIS_FANCONI = list(
      description = "Fanconi anemia versus aplastic anemia diagnosis indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = aplastic anemia (AA / SAA / VSAA)",
      notes = paste(
        "Mohanan 2017 Table 2 footnote and the Results caption define the covariate as",
        "'Diagnosis = (0 - AA/SAA/VSAA (aplastic anemia) or 1 - FA (Fanconi Anemia))'. Cohort",
        "composition: 40 aplastic anemia and 13 Fanconi anemia (Table 1). This is the single",
        "largest effect in the model: clearance in aplastic anemia is 7.12 / 2.90 = 2.46-fold that",
        "in Fanconi anemia, matching the 2.46x figure quoted in the Abstract (P < 1e-6). Both",
        "groups received the same 30 mg/m^2/day fludarabine dose, so the contrast is not",
        "confounded by dose. The paper offers no mechanism and reports (Discussion) that comparing",
        "demographic variables between aplastic-anemia patients above and below 3 L/h/m^2 found no",
        "significant difference, so the indicator should be read as a cohort label rather than an",
        "established physiological covariate. Reference category is aplastic anemia because that is",
        "the group the paper's printed typical value 7.12 L/h/m^2 belongs to; keeping the published",
        "orientation lets Table 2's values be quoted unchanged, following the DIS_OUD / DIS_CHB",
        "precedent in inst/references/covariate-columns.md.",
        sep = " "
      ),
      source_name = "Diagnosis"
    ),
    SNP_NT5E_RS2295890 = list(
      description = "NT5E (CD73) 5'-UTR rs2295890 G>C variant-allele carrier indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = homozygous wild-type GG",
      notes = paste(
        "Mohanan 2017 Results caption: 'rs2295890 = (0 - WT or 1 - heterozygous variant or",
        "Mutant)', i.e. 1 = GC or CC carrier and 0 = GG wild-type, the standard carrier",
        "(dominant-model) encoding of the SNP_<GENE>_<RSID> family. Table 1 genotype counts:",
        "GG 33, GC/CC 15, not available 5; the paper reports the variant-allele frequency as 0.18,",
        "matching the 1000 Genomes expectation. Carriers had 29% lower clearance than wild-type",
        "(7.12 vs 5.03 L/h/m^2, P = 0.038; Abstract and Table 2). NT5E encodes ecto-5'-nucleotidase",
        "/ CD73, which dephosphorylates nucleotides in the purine salvage pathway; the paper states",
        "this is the first report of an effect of NT5E genotype on fludarabine PK and that the",
        "functional relevance was still under evaluation. The SNP was found to be in complete",
        "linkage disequilibrium with rs9450278, rs9450279, rs4599602 and rs4458647 in the same",
        "region, so the column may equivalently be derived from any of those five markers in this",
        "cohort. Three further SNPs were screened (CNT3 rs7853758, NT5C2 rs4917996, hENT1",
        "rs747199) but none was retained in the final model; they are recorded in",
        "covariatesDataExcluded.",
        sep = " "
      ),
      source_name = "rs2295890"
    ),
    OCC = list(
      description = "Daily dosing-occasion index for inter-occasion (inter-day) variability",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Mohanan 2017 estimates an inter-day variability (IDV) alongside inter-individual",
        "variability for all four disposition parameters (Table 2, 'IIV, IDV' block), and Figure 1",
        "contrasts clearance on the first dose against the fifth dose, so the occasion is the",
        "daily dose. Fludarabine is given once daily for six days, day -7 to day -2 (Patients and",
        "Methods), so six occasions are encoded: OCC = 1 is the day -7 dose through OCC = 6, the",
        "day -2 dose. Decomposed inside model() into indicators occ1..occ6 that select per-occasion",
        "etas, following the Ding 2026 vancomycin and Chen 2023 nemonoxacin precedents; rxode2's",
        "'~ var | OCC' syntax parses but cannot be simulated from an rxUi. Only the first eta of",
        "each parameter's set carries the estimated variance and the remaining five are fixed to",
        "the same value, which is the nlmixr2lib encoding of a single shared IOV magnitude. PK was",
        "actually sampled on days -7, -4 and -2 in the first 7 patients and on days -7 and -3 in",
        "the remainder (Results 'F-araA PK'), so not every occasion is informed by data in the",
        "original fit; the six-occasion grid is the dosing schedule, not the sampling schedule. A",
        "user with a different number of doses can add or remove etaiov_*_<k> slots at the same",
        "fixed variances.",
        sep = " "
      ),
      source_name = "day of conditioning"
    )
  )

  covariatesDataExcluded <- list(
    SNP_SLC28A3_RS7853758 = list(
      description = "CNT3 (SLC28A3) exon 6 rs7853758 G>A variant-allele carrier indicator",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened by Mohanan 2017 (Table 1: GG 32, GA/AA 14, not available 7) as a candidate",
        "covariate on fludarabine PK but not retained in the final model; no coefficient is",
        "reported, so nothing can be encoded. CNT3 is a concentrative nucleoside transporter",
        "implicated in fludarabine uptake.",
        sep = " "
      )
    ),
    SNP_NT5C2_RS4917996 = list(
      description = "NT5C2 intronic rs4917996 A>C variant-allele carrier indicator",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened by Mohanan 2017 (Table 1: AA 17, AC/CC 30, not available 6) but not retained in",
        "the final PK model; no coefficient is reported. It was carried into the multivariate",
        "analysis of acute GvHD (spelled rs491799 there, an evident typographical contraction of",
        "rs4917996) but was not significant.",
        sep = " "
      )
    ),
    SNP_SLC29A1_RS747199 = list(
      description = "hENT1 (SLC29A1) exon 1 rs747199 G>C variant-allele carrier indicator",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened by Mohanan 2017 (Table 1: GG 37, GC/CC 12, not available 4) but not retained in",
        "the final model; no coefficient is reported. hENT1 is the equilibrative nucleoside",
        "transporter that mediates cellular uptake of fludarabine.",
        sep = " "
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 53L,
    n_studies = 1L,
    age_range = "3-57 years",
    age_median = "17 years",
    weight_range = "12-89 kg",
    weight_median = "50 kg",
    bsa_range = "0.56-1.9 m^2",
    bsa_median = "1.49 m^2",
    sex_female_pct = 34.0,
    race_ethnicity = "Not reported in source; single-centre Indian cohort",
    disease_state = paste(
      "Aplastic anemia (n = 40) or Fanconi anemia (n = 13) undergoing allogeneic hematopoietic",
      "stem cell transplantation. Conditioning regimens: fludarabine + cyclophosphamide (n = 29),",
      "fludarabine + cyclophosphamide + total body irradiation (n = 20), or fludarabine +",
      "cyclophosphamide + anti-thymocyte globulin (n = 4). Donor: matched sibling (n = 45) or",
      "alternate donor (n = 8). Graft: peripheral blood (n = 48) or bone marrow (n = 2).",
      sep = " "
    ),
    dose_range = paste(
      "Fludarabine phosphate 30 mg/m^2/day as a 1-h intravenous infusion once daily for 6 days,",
      "day -7 to day -2 (Patients and Methods). Identical dosing in both diagnosis groups. Note",
      "that Table 5, which compares this study against previous reports, instead lists the",
      "schedule as day -6 to day -2; the Patients and Methods statement is taken as authoritative.",
      sep = " "
    ),
    regions = "Single centre: Department of Hematology, Christian Medical College, Vellore, India",
    sampling_window = paste(
      "Blood drawn before the start of the infusion (0 h) and at 1, 2, 3, 5, 7 and 24 h after the",
      "start of the fludarabine infusion. PK days were -7, -4, -3 and -2 in the first 7 patients;",
      "subsequent patients were sampled on day -7 and day -3 only, to reduce blood-draw volume",
      "during conditioning.",
      sep = " "
    ),
    assay = paste(
      "F-ara-A quantified in plasma by LC-MS/MS (Shimadzu Prominence UFLC with API2000 triple",
      "quadrupole; MRM 286.0/154.0 for F-ara-A against 5-fluorocytidine internal standard",
      "262.1/130.0). LLOQ 1 ng/mL, linear over 7-7000 ng/mL (1-1000 uM), mean R^2 0.99 +/- 0.001.",
      sep = " "
    ),
    concomitant_medications = paste(
      "Cyclophosphamide 50 or 60 mg/kg/day for 2 days (aplastic anemia) or 10 mg/kg/day for 2 days",
      "(Fanconi anemia) on day -3 to -2; GvHD prophylaxis with ciclosporin 2.5 mg/kg BD plus",
      "methotrexate (n = 32) or post-transplant cyclophosphamide (n = 19); total body irradiation",
      "or anti-thymocyte globulin in a subset. None was tested as a PK covariate.",
      sep = " "
    ),
    notes = paste(
      "Prospective single-centre study of patients transplanted between January 2012 and December",
      "2014. Demographics are Mohanan 2017 Table 1. Non-linear mixed-effects estimation in Monolix",
      "4.3.3 by SAEM. Covariates were retained on a univariate -2 log-likelihood drop of at least",
      "3.84 units (P < 0.05, 1 df); the covariate-model build is tabulated in Table 3. The final",
      "model (diagnosis and rs2295890 on clearance, age on volume and on k21) explained 46% of the",
      "inter-individual variability in clearance. The paper also reports a limited-sampling model",
      "(best 4-point schedule 1, 5, 7 and 24 h) and an exposure-outcome analysis associating",
      "first-dose AUC above 29.4 uM*h with acute GvHD; neither is a pharmacokinetic structural",
      "model and neither is encoded here.",
      sep = " "
    )
  )

  ini({
    # --------------------------------------------------------------------
    # Structural disposition parameters. All values are the 'Final model'
    # column of Mohanan 2017 Table 2. Clearance and volume are reported
    # per square metre of body surface area and are multiplied by BSA in
    # model(); k12 and k21 are first-order rate constants and are not
    # BSA-scaled.
    # --------------------------------------------------------------------
    lcl <- log(7.12)
    label("BSA-normalized clearance in rs2295890 wild-type aplastic-anemia patients (L/h/m^2)")
    # Mohanan 2017 Table 2, final model, row 'WT, AA': CL = 7.12 L/h/m2 (RSE 10%)

    lvc <- log(21.25)
    label("BSA-normalized central volume intercept at age 0 (L/m^2)")
    # Mohanan 2017 Table 2, final model, row 'Age on V': 21.25 x exp(beta x age); RSE 11.8%

    lk12 <- log(0.36)
    label("Transfer rate constant central -> peripheral1 (1/h)")
    # Mohanan 2017 Table 2, final model, row 'k12': 0.36 1/h (RSE 7.6%)

    lk21 <- log(0.14)
    label("Transfer rate constant peripheral1 -> central, intercept at age 0 (1/h)")
    # Mohanan 2017 Table 2, final model, row 'Age on k21': 0.14 x exp(0.016 x age); RSE 10.1%

    # --------------------------------------------------------------------
    # Covariate effects, all in the paper's stated multiplicative form
    # theta = theta_Base * exp(beta * covariate) (Methods, 'F-araA PK and
    # Population PK modeling').
    #
    # Mohanan 2017 prints the four clearance CELLS rather than the two
    # clearance betas, so the two coefficients below are back-calculated
    # from those printed cells as beta = log(cell / 7.12). The
    # back-calculation is self-validating: the two betas are fitted from
    # the 'HET/MUT, AA' and 'WT, FA' cells, and their SUM then predicts the
    # fourth cell 'HET/MUT, FA' as 7.12 * exp(-0.3474877 - 0.8981970) =
    # 2.0487, which rounds to the printed 2.05. A multiplicative structure
    # is therefore confirmed, not assumed.
    # --------------------------------------------------------------------
    e_snp_nt5e_rs2295890_cl <- -0.3474877
    label("Effect of NT5E rs2295890 variant carriage on clearance (log scale)")
    # Mohanan 2017 Table 2: log(5.03 / 7.12); cells 'HET/MUT, AA' and 'WT, AA' (P = 3.8e-02)

    e_fanconi_cl <- -0.8981970
    label("Effect of Fanconi anemia (vs aplastic anemia) on clearance (log scale)")
    # Mohanan 2017 Table 2: log(2.90 / 7.12); cells 'WT, FA' and 'WT, AA' (P = 2.7e-07)

    e_age_vc <- 0.013
    label("Effect of age on central volume (1/year, log scale)")
    # Mohanan 2017 Table 2 row 'Age on V' prints -0.013 (RSE 40.5%, P = 1.3e-02); sign INVERTED to
    # positive per Results, which states that V and k21 INCREASED with age, and per the base-model
    # back-solve documented in covariateData[[AGE]]$notes

    e_age_k21 <- 0.016
    label("Effect of age on the peripheral -> central rate constant k21 (1/year, log scale)")
    # Mohanan 2017 Table 2, row 'Age on k21': 0.14 x exp(0.016 x age) (RSE 26.5%, P = 1.6e-04)

    # --------------------------------------------------------------------
    # Inter-individual variability. Mohanan 2017 Methods: 'The
    # inter-individual and inter-day variability of the parameters was
    # assumed to be log-normally distributed.' Table 2's 'IIV, IDV' block
    # is headed (CV%) and the Results quote the base-model entries as
    # '69 and 39% CV% for clearance and volume', confirming the tabulated
    # numbers are coefficients of variation rather than log-scale SDs.
    # Converted with omega^2 = log(CV^2 + 1).
    #   CL  51% CV -> 0.2311911
    #   V   36% CV -> 0.1218636
    #   k12 41% CV -> 0.1553785
    #   k21 14% CV -> 0.0194104
    # --------------------------------------------------------------------
    etalcl ~ 0.2311911 # Mohanan 2017 Table 2 final model, IIV on CL = 0.51 CV (RSE 11.0%)
    etalvc ~ 0.1218636 # Mohanan 2017 Table 2 final model, IIV on V = 0.36 CV (RSE 13.8%)
    etalk12 ~ 0.1553785 # Mohanan 2017 Table 2 final model, IIV on k12 = 0.41 CV (RSE 15.4%)
    etalk21 ~ 0.0194104 # Mohanan 2017 Table 2 final model, IIV on k21 = 0.14 CV (RSE 65.6%)

    # --------------------------------------------------------------------
    # Inter-day (inter-occasion) variability, same log-normal convention
    # and the same CV -> variance conversion. One shared magnitude per
    # parameter across six daily dosing occasions: the first slot carries
    # the estimated variance and slots 2-6 are fixed to it.
    #   CL  22% CV -> 0.0472652
    #   V   22% CV -> 0.0472652
    #   k12 19% CV -> 0.0354637
    #   k21 20% CV -> 0.0392207
    # --------------------------------------------------------------------
    etaiov_cl_1 ~ 0.0472652 # Mohanan 2017 Table 2 final model, IDV on CL = 0.22 CV (RSE 11.8%)
    etaiov_cl_2 ~ fixed(0.0472652)
    etaiov_cl_3 ~ fixed(0.0472652)
    etaiov_cl_4 ~ fixed(0.0472652)
    etaiov_cl_5 ~ fixed(0.0472652)
    etaiov_cl_6 ~ fixed(0.0472652)

    etaiov_vc_1 ~ 0.0472652 # Mohanan 2017 Table 2 final model, IDV on V = 0.22 CV (RSE 15.8%)
    etaiov_vc_2 ~ fixed(0.0472652)
    etaiov_vc_3 ~ fixed(0.0472652)
    etaiov_vc_4 ~ fixed(0.0472652)
    etaiov_vc_5 ~ fixed(0.0472652)
    etaiov_vc_6 ~ fixed(0.0472652)

    etaiov_k12_1 ~ 0.0354637 # Mohanan 2017 Table 2 final model, IDV on k12 = 0.19 CV (RSE 35.7%)
    etaiov_k12_2 ~ fixed(0.0354637)
    etaiov_k12_3 ~ fixed(0.0354637)
    etaiov_k12_4 ~ fixed(0.0354637)
    etaiov_k12_5 ~ fixed(0.0354637)
    etaiov_k12_6 ~ fixed(0.0354637)

    etaiov_k21_1 ~ 0.0392207 # Mohanan 2017 Table 2 final model, IDV on k21 = 0.20 CV (RSE 34.7%)
    etaiov_k21_2 ~ fixed(0.0392207)
    etaiov_k21_3 ~ fixed(0.0392207)
    etaiov_k21_4 ~ fixed(0.0392207)
    etaiov_k21_5 ~ fixed(0.0392207)
    etaiov_k21_6 ~ fixed(0.0392207)

    # --------------------------------------------------------------------
    # Residual error. Mohanan 2017 Methods: 'A proportional residual error
    # model was used with assumed normal distribution of the residuals.'
    # --------------------------------------------------------------------
    propSd <- 0.19
    label("Proportional residual error (fraction)")
    # Mohanan 2017 Table 2, final model, row 'sigma prop (CV%)': 0.19 (RSE 3.6%)
  })

  model({
    # 1. Inter-occasion (inter-day) variability terms. OCC indexes the
    #    daily fludarabine dose, 1 = day -7 through 6 = day -2. An
    #    observation whose OCC lies outside 1..6 receives no IOV
    #    contribution, which is the intended behaviour for records taken
    #    outside the conditioning window.
    occ1 <- (OCC == 1)
    occ2 <- (OCC == 2)
    occ3 <- (OCC == 3)
    occ4 <- (OCC == 4)
    occ5 <- (OCC == 5)
    occ6 <- (OCC == 6)

    iov_cl <- occ1 * etaiov_cl_1 + occ2 * etaiov_cl_2 + occ3 * etaiov_cl_3 +
      occ4 * etaiov_cl_4 + occ5 * etaiov_cl_5 + occ6 * etaiov_cl_6
    iov_vc <- occ1 * etaiov_vc_1 + occ2 * etaiov_vc_2 + occ3 * etaiov_vc_3 +
      occ4 * etaiov_vc_4 + occ5 * etaiov_vc_5 + occ6 * etaiov_vc_6
    iov_k12 <- occ1 * etaiov_k12_1 + occ2 * etaiov_k12_2 + occ3 * etaiov_k12_3 +
      occ4 * etaiov_k12_4 + occ5 * etaiov_k12_5 + occ6 * etaiov_k12_6
    iov_k21 <- occ1 * etaiov_k21_1 + occ2 * etaiov_k21_2 + occ3 * etaiov_k21_3 +
      occ4 * etaiov_k21_4 + occ5 * etaiov_k21_5 + occ6 * etaiov_k21_6

    # 2. Individual disposition parameters.
    #
    #    Clearance and central volume are published per square metre and
    #    are converted to absolute units by multiplying by BSA, so that
    #    doses can be given in mg. Because both are scaled by the same
    #    BSA, the elimination rate constant cl/vc is BSA-independent, as
    #    it must be.
    #
    #    Age enters V and k21 UNCENTRED: exp(beta * AGE), so exp(lvc) and
    #    exp(lk21) are the values at age 0 rather than at a reference age.
    cl <- exp(lcl + etalcl + iov_cl) * BSA *
      exp(e_snp_nt5e_rs2295890_cl * SNP_NT5E_RS2295890 + e_fanconi_cl * DIS_FANCONI)
    vc <- exp(lvc + etalvc + iov_vc) * BSA * exp(e_age_vc * AGE)
    k12 <- exp(lk12 + etalk12 + iov_k12)
    k21 <- exp(lk21 + etalk21 + iov_k21) * exp(e_age_k21 * AGE)

    # 3. Micro-constants.
    kel <- cl / vc

    # 4. Two-compartment IV disposition. Fludarabine phosphate is given as
    #    a 1-h infusion via the rate / dur column of the event table; the
    #    phosphate prodrug is dephosphorylated essentially immediately, so
    #    the infused dose enters the central F-ara-A compartment directly
    #    and no depot state is required.
    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 5. Observation. Amount in mg divided by volume in L gives mg/L,
    #    equivalently ug/mL. The paper reports plasma concentrations in
    #    ng/mL and exposures in uM*h; F-ara-A has a molar mass of
    #    285.23 g/mol, so 1 mg/L = 1000 ng/mL = 3.506 uM.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
