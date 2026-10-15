Dadkhah_2022_busulfan <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for intravenous busulfan",
    "and its metabolite sulfolane in adults with myelofibrosis receiving",
    "reduced-intensity busulfan/fludarabine conditioning before allogeneic",
    "haematopoietic stem cell transplantation (Dadkhah 2022). One compartment",
    "for each analyte with first-order elimination; a fixed metabolic",
    "fraction (MF = 0.0704) of total busulfan clearance forms sulfolane.",
    "Total body weight enters busulfan volume as a power function",
    "(exponent 0.854, 75 kg reference) and the GSTA1 -52G>A (rs3957356)",
    "SNP multiplies sulfolane clearance by exp(1.43). Log-normal IIV on all",
    "four clearances and volumes, inter-occasion variability on busulfan",
    "clearance (one occasion per sampled administration, up to three), and",
    "a proportional residual error per analyte.",
    sep = " "
  )
  reference <- paste(
    "Dadkhah A, Wicha SG, Kroger N, Muller A, Pfaffendorf C, Riedner M,",
    "Badbaran A, Fehse B, Langebrake C. Population Pharmacokinetics of",
    "Busulfan and Its Metabolite Sulfolane in Patients with Myelofibrosis",
    "Undergoing Hematopoietic Stem Cell Transplantation. Pharmaceutics.",
    "2022;14(6):1145. doi:10.3390/pharmaceutics14061145.",
    "Structure from Figure 1; covariate equations (1) and (2) from",
    "Section 3.3; parameter estimates from Table 2.",
    sep = " "
  )
  vignette <- "Dadkhah_2022_busulfan"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power function on busulfan volume normalised to 75 kg, the cohort",
        "median (Section 3.3 equation (1): VBu = VBu-typ * (TBW/75)^0.854;",
        "Table 1 weight median 75 kg, IQR 64.05-88.25). Weight does not",
        "enter busulfan clearance or either sulfolane parameter.",
        sep = " "
      ),
      source_name = "TBW"
    ),
    SNP_GSTA1_RS3957356 = list(
      description = "GSTA1 -52G>A (rs3957356) variant-carrier indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no GSTA1 -52G>A SNP detected)",
      notes = paste(
        "1 = the GSTA1 -52G>A SNP (the GSTA1*B promoter haplotype marker)",
        "was detected, 0 = not detected (Section 2.3; Section 2.5",
        "'Categorical covariates (GSTA1 SNP, ...) were coded as 0 or 1').",
        "The paper does not separate heterozygous from homozygous carriers.",
        "Exponential effect on sulfolane clearance only (Section 3.3",
        "equation (2): CLSu = CLSu-typ * exp(1.43) 'for GSTA1', CLSu =",
        "CLSu-typ 'for non-GSTA1'), so the Table 2 CLSu of 1.61 L/h is the",
        "non-carrier value and carriers clear sulfolane 4.18-fold faster.",
        "28 of 37 patients (75.7%) carried the SNP (Table 1).",
        sep = " "
      ),
      source_name = "GSTA1"
    ),
    OCC = list(
      description = "Integer-valued occasion index for inter-occasion variability on busulfan clearance",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Section 3.2: 'each administration, followed by blood sampling, was",
        "defined as a new occasion'. Q6H patients were sampled after the",
        "first, fifth (trough only) and ninth infusions (Section 2.2), giving",
        "up to three occasions; Q24H patients were sampled after the first",
        "infusion only. Values 1, 2 and 3 each select their own IOV eta on",
        "log busulfan clearance (one shared variance); any other value",
        "carries no IOV. Number occasions consecutively per subject by",
        "sampled administration, and keep each occasion's value from its",
        "dose until the next sampled dose.",
        sep = " "
      ),
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    JAK2_MUT = list(
      description = "JAK2 driver mutation",
      units = "(binary)",
      type = "binary",
      notes = "Significant on sulfolane clearance in forward inclusion (dOFV -12.6) but discarded because the parameter estimates became physiologically implausible (Section 3.3)."
    ),
    GSTM1_NULL = list(
      description = "GSTM1 gene deletion",
      units = "(binary)",
      type = "binary",
      notes = "Significant in forward inclusion (Section 3.3) but not retained in the final model."
    ),
    AST_ALT_RATIO = list(
      description = "De Ritis ratio (AST/ALT)",
      units = "(ratio)",
      type = "continuous",
      notes = "Reduced the OFV by 7.61 on busulfan clearance but caused high RSEs and was not retained (Section 3.3)."
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units = "U/L",
      type = "continuous",
      notes = "Power effect on busulfan clearance significant in forward inclusion (dOFV -5.56), removed in backward elimination (Section 3.3)."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Power effect on busulfan clearance significant in forward inclusion (dOFV -5.23), removed in backward elimination (Section 3.3). The paper reports bilirubin in mg/dL (Table 1 median 0.6 mg/dL = 10.3 umol/L)."
    )
  )

  compartmentData <- list(
    central = list(analyte = "busulfan", units = "mg", specimen = "plasma", verified = TRUE),
    central_sulfolane = list(analyte = "sulfolane", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 37L,
    n_studies = 1L,
    age_range = "median 60 years (IQR 53.5-65.5); adults >= 18 years",
    weight_range = "median 75 kg (IQR 64.05-88.25)",
    height = "median 174 cm (IQR 168-181)",
    bsa = "median 1.84 m^2 (IQR 1.75-2.07)",
    sex_female_pct = 51.4,
    race_ethnicity = "Not reported.",
    disease_state = paste(
      "Myelofibrosis scheduled for allogeneic haematopoietic stem cell",
      "transplantation after reduced-intensity busulfan/fludarabine",
      "conditioning (with anti-thymocyte globulin and levetiracetam",
      "prophylaxis): primary myelofibrosis 18 (48.7%), post-essential",
      "thrombocythaemia 9 (24.3%), post-polycythaemia vera 10 (27%).",
      "JAK2 mutation 26 (70.3%), CALR 7 (18.9%), MPL 1 (3%).",
      sep = " "
    ),
    dose_range = paste(
      "Q6H regimen (30 patients, 81%): 10 doses of 0.8 mg/kg intravenous",
      "busulfan every 6 h, each a 2 h infusion. Q24H regimen (7 patients,",
      "19%): 3.2 mg/kg as a 3 h infusion every 24 h for three doses, with",
      "doses adjusted where needed to a cumulative AUC of 50 mg*h/L.",
      sep = " "
    ),
    genotype = "GSTA1 -52G>A (rs3957356) SNP 28 (75.7%); GSTM1 deletion 19 (51.35%); both 10 (27%).",
    hepatic_function = "AST median 21 U/L, ALT median 21 U/L, alkaline phosphatase median 85 U/L, bilirubin median 0.6 mg/dL, albumin median 37.8 g/L (Table 1).",
    regions = "Germany (University Medical Center Hamburg-Eppendorf; prospective November 2018 - June 2020 plus retrospective October 2016 - October 2017).",
    n_observations = "523 plasma concentrations: 282 busulfan and 241 sulfolane. Sulfolane was measured only in the 30 prospectively enrolled patients; 70 sulfolane samples were below the 0.04 mg/L LLOQ and were kept in the fit.",
    notes = "Demographics from Table 1 and Section 3.1. NONMEM 7.4.3, FOCE-I; parameter uncertainty from SIR (M/m = 5000/1000)."
  )

  ini({
    # Busulfan (parent) -- Table 2 'Final Model' estimates.
    lcl <- log(16.3); label("Busulfan total clearance CLBu (L/h)") # Table 2 row 'CL Bu [L/h]' = 16.3 (RSE 3.6%, SIR 95% CI 15.18-17.35)
    lvc <- log(61.5); label("Busulfan volume of distribution VBu for a 75 kg patient (L)") # Table 2 row 'V Bu [L]' = 61.5 (RSE 2%, SIR 95% CI 59.37-63.78)
    e_wt_vc <- 0.854; label("Power exponent of total body weight on VBu (unitless)") # Table 2 row 'COV_V Bu _TBW [kg]' = 0.854 (RSE 11.6%); Section 3.3 equation (1)

    # Formation of sulfolane: Figure 1 arrow 'CLBu x MF' out of the busulfan
    # compartment; the caption states the formation is part of the total
    # busulfan clearance, so fm is a share of CLBu, not an extra arm.
    fm <- 0.0704; label("Metabolic fraction MF of busulfan clearance forming sulfolane, mass basis (fraction)") # Table 2 row 'MF' = 0.0704 (RSE 28.6%, SIR 95% CI 0.0463-0.1029)

    # Sulfolane (metabolite). The '_sulfolane' suffix is registered as a
    # metabolite-suffix entry in inst/references/compartment-names.md.
    lcl_sulfolane <- log(1.61); label("Sulfolane clearance CLSu in GSTA1 SNP non-carriers (L/h)") # Table 2 row 'CL Su [L/h]' = 1.61 (RSE 37%, SIR 95% CI 0.84-2.24)
    lvc_sulfolane <- log(48.8); label("Sulfolane volume of distribution VSu (L)") # Table 2 row 'V Su [L]' = 48.8 (RSE 35.2%, SIR 95% CI 30.75-78.46)
    e_gsta1_cl_sulfolane <- 1.43; label("Exponential effect of the GSTA1 -52G>A SNP on log CLSu (unitless)") # Table 2 row 'COV_CL Su _GSTA1' = 1.43 (RSE 43.6%); Section 3.3 equation (2)

    # IIV. Table 2 reports %CV with the footnote '%CV = sqrt(exp(OMEGA)-1) *
    # 100', so omega^2 = log(1 + CV^2).
    etalcl ~ 0.045171 # Table 2 row 'IIV CL Bu [CV%]' = 21.5; log(1 + 0.215^2) = 0.045171
    etalvc ~ 0.009950 # Table 2 row 'IIVV Bu [CV%]' = 10; log(1 + 0.10^2) = 0.009950
    etalcl_sulfolane ~ 0.820843 # Table 2 row 'IIV CL Su [CV%]' = 112.8; log(1 + 1.128^2) = 0.820843
    etalvc_sulfolane ~ 0.471361 # Table 2 row 'IIVV Su [CV%]' = 77.6; log(1 + 0.776^2) = 0.471361

    # IOV on log busulfan clearance, one shared variance across the (up to)
    # three sampled occasions.
    etaiov_cl_1 ~ 0.005759 # Table 2 row 'IOV CL Bu [CV%]' = 7.6; log(1 + 0.076^2) = 0.005759
    etaiov_cl_2 ~ fixed(0.005759) # same variance as occasion 1 (single IOV estimate in Table 2)
    etaiov_cl_3 ~ fixed(0.005759) # same variance as occasion 1 (single IOV estimate in Table 2)

    # Residual error: proportional per analyte (Section 3.2). Table 2 reports
    # each as a CV%, taken here as the proportional SD. The 11.8% correlation
    # between the two residual errors (NONMEM L2 sigma block, Section 3.2)
    # cannot be expressed in nlmixr2 and is omitted; see the vignette.
    propSd <- 0.071; label("Busulfan proportional residual error (fraction)") # Table 2 row 'Prop. sigma Bu [CV%]' = 7.1 (RSE 12.8%)
    propSd_sulfolane <- 0.362; label("Sulfolane proportional residual error (fraction)") # Table 2 row 'Prop. sigma Su [CV%]' = 36.2 (RSE 7.2%)
  })

  model({
    # Occasion indicators for the IOV etas on log busulfan clearance
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3

    # Busulfan: Section 3.3 equation (1), weight on volume only
    cl <- exp(lcl + etalcl + iov_cl)
    vc <- exp(lvc + etalvc) * (WT / 75)^e_wt_vc

    # Sulfolane: Section 3.3 equation (2), exponential GSTA1 SNP effect
    cl_sulfolane <- exp(lcl_sulfolane + etalcl_sulfolane + e_gsta1_cl_sulfolane * SNP_GSTA1_RS3957356)
    vc_sulfolane <- exp(lvc_sulfolane + etalvc_sulfolane)

    kel <- cl / vc
    kel_sulfolane <- cl_sulfolane / vc_sulfolane

    # Figure 1: total busulfan elimination CLBu, of which MF * CLBu forms
    # sulfolane. Amounts in mg on both sides (concentrations in mg/L; the
    # paper describes no molar conversion, so MF is a mass fraction).
    d/dt(central) <- -kel * central
    d/dt(central_sulfolane) <- fm * kel * central - kel_sulfolane * central_sulfolane

    Cc <- central / vc
    Cc_sulfolane <- central_sulfolane / vc_sulfolane

    Cc ~ prop(propSd)
    Cc_sulfolane ~ prop(propSd_sulfolane)
  })
}
