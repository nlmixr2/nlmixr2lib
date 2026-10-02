Li_2019_cyclosporine <- function() {
  description <- "One-compartment intravenous population PK model for cyclosporine A in Chinese children receiving allogeneic haematopoietic stem cell transplantation (Li 2019), with allometric body-weight scaling on CL and V, linear post-operative-day effects on CL and V, a power eGFR effect and a fractional triazole-antifungal effect on CL, CYP3A4*1G (rs2242480) genotype multipliers on CL, and combined proportional plus additive residual error."
  reference <- "Li TF, Hu L, Ma XL, Huang L, Liu XM, Luo XX, Feng WY, Wu CF. Population pharmacokinetics of cyclosporine in Chinese children receiving hematopoietic stem cell transplantation. Acta Pharmacol Sin. 2019;40(12):1603-1610. doi:10.1038/s41401-019-0277-x"
  vignette <- "Li_2019_cyclosporine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    # Li 2019 Abstract: "Whole blood samples were collected" and the CMIA
    # (Architect i2000SR) cyclosporine assay is a whole-blood assay, so the
    # modelled matrix is whole blood.
    central = list(analyte = "cyclosporine", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling (WT / 70)^0.75 on CL and (WT / 70)^1 on V, both",
        "exponents fixed (Li 2019 Eq. 3 and text: 'PWR is the allometric",
        "coefficient fixed at a value of 0.75 for clearance and a value of 1",
        "for distribution volume'). Cohort mean 31.93 kg, median 28.8 kg",
        "(Table 1). The Results text gives the range 6.5-69.0 kg while",
        "Table 1 gives 9.4-78.5 kg; the source does not resolve the",
        "discrepancy."
      ),
      source_name = "WT"
    ),
    POD = list(
      description = "Post-operative day (days since allogeneic HSCT)",
      units = "days",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying. Enters CL and V as the linear deviation",
        "(1 - (POD - 9) * theta) (Li 2019 Eq. 4 and Eq. 5), so both",
        "parameters decline linearly with POD; the centring value 9 days is",
        "printed in the equations. Table 1: mean 29.7, median 25, range",
        "0-67 days. Cyclosporine was started on Day -7 or -10; how the",
        "pre-transplant troughs were coded is not stated, but Table 1's",
        "lower bound of 0 suggests they carried POD = 0. The V factor reaches",
        "zero at POD = 9 + 1 / 0.0197 = 59.8 days and is negative beyond, so",
        "the published equation is only usable for POD below about 50 days",
        "(the model does not truncate it)."
      ),
      source_name = "POD"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate, BSA-normalised",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect (eGFR / 172.46)^0.545 on CL (Li 2019 Eq. 4); 172.46 is",
        "the Table 1 cohort median (mean 179.4, range 24.93-377.54). The",
        "estimating equation for eGFR is not stated in the paper (a",
        "creatinine-based paediatric equation such as bedside Schwartz is",
        "the usual choice; Table 1 reports serum creatinine in umol/L)."
      ),
      source_name = "eGFR"
    ),
    CONMED_AZOLE = list(
      description = "Concomitant triazole antifungal (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no triazole antifungal)",
      notes = paste(
        "Li 2019 'TAF' indicator ('patients were coadministered of TAF, the",
        "TAF = 1, otherwise TAF = 0'); fractional effect (1 - 0.36 * TAF) on",
        "CL (Eq. 4). The pooled triazoles were voriconazole (61 children),",
        "itraconazole (19), posaconazole (6) and fluconazole (1); only 4",
        "children never received one (Table 1). The paper does not say",
        "whether TAF was coded per record (time-varying) or per subject."
      ),
      source_name = "TAF"
    ),
    SNP_CYP3A4_RS2242480_VAR_COUNT = list(
      description = "CYP3A4*1G (rs2242480) variant (T) allele count: 0 = CC, 1 = CT, 2 = TT",
      units = "(count, 0/1/2 alleles per subject)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Li 2019 codes the genotype as 'Gene' = 1 for CC, 2 for TT and 3 for",
        "CT (Table 3 footnote) and multiplies CL by 0.984 for Gene = 1 and by",
        "1.22 for Gene = 2 or 3 (Results text below Eq. 5; Table 3 rows",
        "'Gene-1/2/3 ON CL'). Canonical mapping: CC -> 0, CT -> 1, TT -> 2.",
        "Table 2: CC 50 (61.0%), CT 26 (31.7%), TT 6 (7.3%); the T allele is",
        "the minor allele (frequency 23.2%). 82 of 86 children were",
        "genotyped; how the 4 ungenotyped children were coded is not",
        "reported, so no missing-genotype indicator is carried."
      ),
      source_name = "Gene"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Table 1 mean 8.38 (1.1-16.8). Screened; the age trend in CL disappeared after allometric weight scaling (Fig. 1)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "30 of 86 female (Table 1). Screened, not retained."
    ),
    HGB = list(
      description = "Haemoglobin",
      units = "g/L",
      type = "continuous",
      notes = "Table 1 mean 85.49 g/L. Screened, not retained."
    ),
    HCT = list(
      description = "Haematocrit",
      units = "%",
      type = "continuous",
      notes = "Table 1 mean 24.42%. Screened, not retained."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Table 1 mean 34.98 (printed with the unit 'U/L'). Screened, not retained."
    ),
    SNP_CYP3A5_RS776746 = list(
      description = "CYP3A5*3 (rs776746) genotype",
      units = "(genotype)",
      type = "categorical",
      notes = "Table 2: CC 42 / CT 34 / TT 5. Screened, not significant (Discussion)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 86L,
    n_studies = 1L,
    n_observations = 1010L,
    age_range = "1.1-16.8 years",
    age_mean = "8.38 years (SD 3.78; median 8.35)",
    weight_range = "9.4-78.5 kg (Table 1; the Results text gives 6.5-69.0 kg)",
    weight_median = "28.8 kg (mean 31.93, SD 16.75)",
    sex_female_pct = 34.9,
    race_ethnicity = "Chinese (single centre, Beijing)",
    disease_state = "Children with malignant haematological disorders (ALL 33.7%, AA 26.7%, AML 25.6%, NHL 5.8%, MDS 4.7%, other 3.5%) receiving allogeneic haematopoietic stem cell transplantation",
    dose_range = "Intravenous cyclosporine as a 2 h infusion every 12 h from Day -7 or -10, initial 2-3 mg/kg, then adjusted to a trough target of 150-250 ng/mL; per-dose amount mean 28.6 mg (median 25, range 5-125 mg)",
    regions = "China (Peking University People's Hospital, Beijing)",
    cyp3a4_1g_distribution = "rs2242480 CC 50 (61.0%), CT 26 (31.7%), TT 6 (7.3%) of 82 genotyped",
    pod_range = "0-67 days (mean 29.7, median 25)",
    egfr_median = "172.46 mL/min/1.73 m^2 (range 24.93-377.54)",
    notes = "Retrospective therapeutic-drug-monitoring data: 1010 whole-blood troughs (mean 12 per child) drawn before the morning intravenous infusion and measured by chemiluminescent microparticle immunoassay (LLOQ 30 ng/mL). Fit by NONMEM VII FOCE-I."
  )

  ini({
    # Structural parameters at the Eq. 4 / Eq. 5 reference subject: 70 kg,
    # POD 9 days, eGFR 172.46 mL/min/1.73 m^2, no triazole antifungal. The
    # CYP3A4*1G multiplier (0.984 CC / 1.22 T carrier) is applied on top of
    # the printed CL for every genotype, so no genotype reproduces CL = 42.3
    # exactly. Li 2019 Table 3 'Final model Estimate (RSE%)'.
    #
    # The data are trough-only, so V is identified from accumulation across
    # days; the resulting typical half-life is long (about 50 h at 70 kg),
    # as in the sibling Feng 2023 cyclosporine model.
    lcl <- log(42.3); label("Typical clearance CL at 70 kg, POD 9 d, eGFR 172.46, no triazole (L/h)") # Li 2019 Table 3: CL = 42.3 L/h (RSE 10.6%)
    lvc <- log(3100); label("Typical volume of distribution V at 70 kg, POD 9 d (L)") # Li 2019 Table 3: V = 3100 L (RSE 13.1%)

    # Allometric exponents fixed (Li 2019 Eq. 3 and text).
    e_wt_cl <- fixed(0.75); label("Allometric exponent of (WT / 70) on CL (unitless)") # Li 2019 Methods 'Covariate analysis': PWR fixed at 0.75 for clearance
    e_wt_vc <- fixed(1); label("Allometric exponent of (WT / 70) on V (unitless)") # Li 2019 Methods 'Covariate analysis': PWR fixed at 1 for distribution volume

    # Covariate effects, Li 2019 Eq. 4 and Eq. 5:
    #   CL = CLpop * (WT/70)^0.75 * (1 - theta_TAF * TAF) * (1 - (POD - 9) * theta_POD_CL)
    #        * (eGFR / 172.46)^theta_eGFR * exp(eta), then * 0.984 (CC) or * 1.22 (CT, TT)
    #   V  = Vpop * (WT/70) * (1 - (POD - 9) * theta_POD_V) * exp(eta)
    e_azole_cl <- 0.36; label("Fractional decrease in CL with a concomitant triazole antifungal (unitless)") # Li 2019 Table 3: theta TAF-CL = 0.36 (RSE 16.3%)
    e_pod_cl <- 0.00703; label("Fractional decrease in CL per post-operative day above 9 d (1/day)") # Li 2019 Table 3: theta POD-CL = 0.00703 (RSE 21.5%)
    e_pod_vc <- 0.0197; label("Fractional decrease in V per post-operative day above 9 d (1/day)") # Li 2019 Table 3: theta POD-V = 0.0197 (RSE 29.8%)
    e_crcl_cl <- 0.545; label("Power exponent of (eGFR / 172.46) on CL (unitless)") # Li 2019 Table 3: theta eGFR-CL = 0.545 (RSE 29.9%)
    e_cyp3a4_wild_cl <- 0.984; label("CL multiplier for CYP3A4*1G (rs2242480) CC (unitless)") # Li 2019 Table 3: Gene-1 ON CL = 0.984 (RSE 8.1%); Gene 1 = CC
    e_cyp3a4_varhom_cl <- 1.22; label("CL multiplier for CYP3A4*1G (rs2242480) TT (unitless)") # Li 2019 Table 3: Gene-2 ON CL = 1.22 (RSE 10.1%); Gene 2 = TT
    e_cyp3a4_het_cl <- 1.22; label("CL multiplier for CYP3A4*1G (rs2242480) CT (unitless)") # Li 2019 Table 3: Gene-3 ON CL = 1.22 (RSE 8.9%); Gene 3 = CT

    # Inter-individual variability, exponential (Li 2019 Eq. 1). Table 3
    # labels the rows 'omega^2', so they are variances:
    #   CL 0.0744 -> CV 27.8%;  V 0.454 -> CV 73.6%
    etalcl ~ 0.0744 # Li 2019 Table 3: omega^2 CL = 0.0744 (RSE 19.5%, shrinkage 8.4%)
    etalvc ~ 0.454 # Li 2019 Table 3: omega^2 V = 0.454 (RSE 24%, shrinkage 27.2%)

    # Residual error, Li 2019 Eq. 2: Cobs = Cpred * (1 + eps1) + eps2. Table 3
    # heads the rows 'sigma^2 1 (%)' = 0.154 and 'sigma^2 2 (ng/mL)' = 30.3.
    # The units printed (ng/mL, not (ng/mL)^2) and the Results text ('the
    # additive error was 30.3 ng/mL') identify these as standard deviations;
    # the text's 'proportional error was 22.9%' is the RSE column of the same
    # row. See the vignette for the scale decision and the Figure 2b check.
    propSd <- 0.154; label("Proportional residual error (fraction)") # Li 2019 Table 3: sigma 1 = 0.154 (RSE 22.9%); read as an SD
    addSd <- 30.3; label("Additive residual error (ng/mL)") # Li 2019 Table 3: sigma 2 = 30.3 ng/mL (RSE 21.7%); Results text 'additive error was 30.3 ng/mL'
  })
  model({
    # 1. CYP3A4*1G genotype multiplier (Li 2019 Results text below Eq. 5):
    # CC -> 0.984, CT or TT -> 1.22.
    cyp3a4_cl <- e_cyp3a4_wild_cl * (SNP_CYP3A4_RS2242480_VAR_COUNT == 0) +
      e_cyp3a4_het_cl * (SNP_CYP3A4_RS2242480_VAR_COUNT == 1) +
      e_cyp3a4_varhom_cl * (SNP_CYP3A4_RS2242480_VAR_COUNT == 2)

    # 2. Individual parameters (Li 2019 Eq. 4 and Eq. 5)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (1 - e_azole_cl * CONMED_AZOLE) *
      (1 - (POD - 9) * e_pod_cl) * (CRCL / 172.46)^e_crcl_cl * cyp3a4_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc * (1 - (POD - 9) * e_pod_vc)

    # 3. Micro-constant and ODE. All modelled troughs followed intravenous
    # 2 h infusions (Li 2019 'CsA administration'), so doses enter central.
    kel <- cl / vc
    d/dt(central) <- -kel * central

    # 4. Observation: mg / L * 1000 = ng/mL
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
