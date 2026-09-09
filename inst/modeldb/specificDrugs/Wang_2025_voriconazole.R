Wang_2025_voriconazole <- function() {
  description <- "One-compartment population pharmacokinetic model with linear elimination for intravenous voriconazole in critically ill adults with COVID-19-associated pulmonary aspergillosis (Wang 2025); clearance carries five covariates - platelet count, C-reactive protein, gamma-glutamyltransferase, aspartate aminotransferase and continuous renal replacement therapy - with continuous renal replacement therapy raising clearance 1.617-fold, and no interindividual variability on the volume of distribution"
  reference <- "Wang H, Shen Y, Luo X, Jin L, Zhu H, Wang J. Population pharmacokinetics and dose optimization of voriconazole in patients with COVID-19-associated pulmonary aspergillosis. Front Pharmacol. 2025;16:1554370. doi:10.3389/fphar.2025.1554370"
  vignette <- "Wang_2025_voriconazole"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Wang 2025 Section 2.1 (inclusion
  # criterion 3 restricts the cohort to patients "receiving an intravenous
  # infusion of at least 72 h of voriconazole", so there is no absorption
  # compartment) and Section 2.3 (plasma HPLC assay, UV detection at 262 nm,
  # calibration range 0.1-30 mg/L).
  compartmentData <- list(
    central = list(analyte = "voriconazole", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    PLT = list(
      description        = "Platelet count",
      units              = "10^9 cells/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters CL as (PLT/121)^0.248. The 121 x 10^9/L divisor appears ONLY inside the final-model equation printed in Wang 2025 Section 3.2 and does not equal the Table 1 cohort median of 126 x 10^9/L (IQR 73.63-196); the printed equation is used, matching the precedent already recorded under this canonical for Stitt 2026 (equation divisor 196 against a table median of 197) and for Wang 2024. The most likely explanation for the four-covariate-wide mismatch in this paper is that Table 1 summarises the 72 patients while the equation divisors are medians over the 150 concentration records. Wang 2025 Discussion reads the positive exponent as a liver-function marker rather than a platelet-mediated mechanism, citing Tang 2021 and noting that portal hypertension and reduced thrombopoietin lower the platelet count as liver function deteriorates.",
      source_name        = "PLT"
    ),
    CRP = list(
      description        = "C-reactive protein, standard (not high-sensitivity) clinical-chemistry assay",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters CL as (CRP/55.23)^-0.183. The 55.23 mg/L divisor appears ONLY inside the final-model equation in Wang 2025 Section 3.2 and does not equal the Table 1 cohort median of 60.4 mg/L (IQR 19.5-108.5). Beware a second, unrelated discrepancy: the Discussion states 'The median (IQR) of CRP in this study was 84 (56.08, 120.04) mg/L', but 84 (56.08, 120.04) is verbatim the Table 1 SCR (serum creatinine, umol/L) row, so that sentence mis-cites the creatinine row and is NOT a third candidate reference value. The negative exponent means clearance falls as inflammation rises; the paper attributes this to cytokine-mediated downregulation of CYP enzyme expression and argues at length that in this hypoalbuminaemic, severely infected cohort CRP affects clearance more than CYP2C19 genotype does.",
      source_name        = "CRP"
    ),
    GGT = list(
      description        = "Serum gamma-glutamyltransferase activity",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters CL as (GGT/68.72)^0.292. The 68.72 U/L divisor appears ONLY inside the final-model equation in Wang 2025 Section 3.2 and does not equal the Table 1 cohort median of 65.53 IU/L (IQR 41.79-136.32); the source reports the analyte in IU/L, which is used interchangeably with the canonical U/L. Wang 2025 groups GGT with AST as 'widely recognized biomarkers of liver function' and cites Li 2017 and Chantharit 2020 as prior voriconazole models retaining hepatic markers on clearance. Note that the Section 3.2 equation prints the exponent to four figures as 0.2928 while Table 2 gives 0.292; the Table 2 value is used here and the difference is immaterial (0.03% on clearance at a tenfold GGT deviation).",
      source_name        = "GGT"
    ),
    AST = list(
      description        = "Serum aspartate aminotransferase activity",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters CL as (AST/29.91)^-0.227. The 29.91 U/L divisor appears ONLY inside the final-model equation in Wang 2025 Section 3.2 and does not equal the Table 1 cohort median of 26.78 IU/L (IQR 19.2-43.68); the source reports the analyte in IU/L, which is used interchangeably with the canonical U/L. The negative exponent is the expected direction for a hepatically metabolised triazole - rising transaminase marks hepatocellular injury and slower clearance - and is the opposite sign to the GGT term retained in the same equation, which the paper does not comment on.",
      source_name        = "AST"
    ),
    RRT_CRRT_STATUS = list(
      description        = "Continuous renal replacement therapy received during voriconazole treatment (1) or not (0)",
      units              = "binary",
      type               = "binary",
      reference_category = "0 (no continuous renal replacement therapy)",
      notes              = "Subject-level rather than record-level, which is why the STATUS member of the RRT family is used rather than RRT_CRRT_ACTIVE: Wang 2025 Table 1 tabulates the covariate as 'CRRT during voriconazole therapy, n (%) of patients', 25 of 72 (34.7%), and the Monte Carlo simulations of Tables 3 and 4 stratify whole simulated patients into CRRT and non-CRRT arms. All 25 were treated in continuous veno-venous hemofiltration (CVVH) mode with blood flow 150-180 mL/min and dialysate flow 2 L/h (Wang 2025 Section 3.1), a continuous modality. Enters CL as the piecewise fractional-change form theta1 * (1 + theta2) that Wang 2025 Section 2.4.2 declares for categorical covariates, with theta2 = 0.617, i.e. a 1.617-fold multiplier while CRRT is in use. This is the paper's headline finding and the paper itself calls it 'unexpected'; the Discussion attributes it to convective solute removal across the hemofilter combined with hypoalbuminaemia in this cohort, and corroborates the direction against Wang 2024.",
      source_name        = "CRRT"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Tabulated in Wang 2025 Table 1 (median 66.5 kg, IQR 60-70) but explicitly NOT screened: Discussion limitation 2 states that most participants 'had been bedbound for an extended period, resulting in challenges in obtaining accurate weight data, with over 30% missing values. Consequently, weight was not incorporated into the modeling process.' This model therefore carries NO allometric term and both CL and V are absolute population values, not per-70-kg values - even though the paper's own Monte Carlo dose recommendations in Tables 3 and 4 are expressed in mg/kg."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as a continuous covariate (Wang 2025 Section 2.2; cohort median 77 years, IQR 69-84) but not retained in the final model. This is an unusually elderly critical-care cohort."
    ),
    SEXF = list(
      description = "Female sex",
      units       = "binary",
      type        = "binary",
      notes       = "Screened as a categorical covariate (Wang 2025 Section 2.2, recorded as 0 = female, 1 = male; 18 of 72 participants, 25.0%, were women) but not retained in the final model."
    ),
    CYP2C19_PHENOTYPE = list(
      description = "CYP2C19 metabolizer phenotype inferred from the *2 (rs4244285), *3 (rs4986893) and *17 (rs12248560) alleles",
      units       = "categorical",
      type        = "categorical",
      notes       = "Genotyped by Sanger sequencing in every participant and screened as a three-level covariate (UM/EM, IM, PM; observed counts 28 EM plus 1 RM, 39 IM, 4 PM per Wang 2025 Section 3.1) but NOT retained in the final model. The paper treats this as a substantive negative result rather than an omission: its Discussion argues that in this hypoalbuminaemic, severely infected cohort C-reactive protein masks the genotype effect, citing Li 2024 and Hao 2023 for the same masking in CRP-elevated populations. rs4244285 was the one SNP that departed from Hardy-Weinberg equilibrium (p < 0.05)."
    ),
    SNP_CYP3A4_RS4646437 = list(
      description = "CYP3A4 rs4646437 (c.671-202C>T) genotype",
      units       = "categorical",
      type        = "categorical",
      notes       = "Genotyped and screened as a two-level covariate (group 1 G/G, n = 59; group 2 G/A plus A/A, n = 13, per Wang 2025 Table 1) but not retained in the final model."
    ),
    CLCR = list(
      description = "Creatinine clearance by the Cockcroft-Gault equation",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened as a continuous covariate (Wang 2025 Section 2.2; Table 1 median 54.97 mL/min, IQR 36.04-96.6) but not retained. The estimated glomerular filtration rate (Table 1 median 78.2 mL/min/1.73 m^2) and serum creatinine (median 84 umol/L) were screened in the same step and likewise not retained; they are recorded here rather than as their own entries because this model uses none of the three. Note the contrast with the sibling Wang 2024 voriconazole model, which retained creatinine clearance on clearance."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as a continuous covariate (Wang 2025 Section 2.2; Table 1 median 21.4 IU/L) but not retained, unlike its transaminase partner AST. Wang 2025 Section 2.4.2 states that covariates correlated at r > 0.5 were not entered concurrently, which plausibly explains why only one of the ALT/AST pair survives."
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as a continuous covariate (Wang 2025 Section 2.2; Table 1 median 77.95 IU/L) but not retained, unlike its cholestatic partner GGT."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened as a continuous covariate (Wang 2025 Section 2.2; Table 1 median 9.68 umol/L) but not retained. Direct bilirubin (median 3.75 umol/L) was screened in the same step and likewise not retained."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened as a continuous covariate (Wang 2025 Section 2.2; Table 1 median 32.6 g/L) but not retained. The Discussion nevertheless makes hypoalbuminaemia load-bearing for its interpretation of two retained covariates: 45 of 72 patients (62.5%) were hypoalbuminaemic, which the paper invokes both to explain why CRP outweighs CYP2C19 genotype and to explain why CRRT increases the clearance of a drug that is only 58% protein-bound."
    ),
    CONMED_PPI = list(
      description = "Concomitant proton pump inhibitor",
      units       = "binary",
      type        = "binary",
      notes       = "Screened as a categorical covariate (Wang 2025 Table 1; 58 of 72 patients, 80.6%, split across omeprazole, esomeprazole, lansoprazole and pantoprazole) but not retained. The Discussion says the combination categories were too thinly populated to survive: 'only glucocorticoids and proton pump inhibitors were included, with few combinations in each category, which could not be integrated into the model after broad categorization.' Recorded despite the negative result because omeprazole-class agents have a well-known CYP2C19 interaction with voriconazole."
    ),
    CONMED_STEROID = list(
      description = "Concomitant systemic glucocorticoid",
      units       = "binary",
      type        = "binary",
      notes       = "Screened as a categorical covariate (Wang 2025 Table 1; 61 of 72 patients, 84.7%, across methylprednisolone, prednisolone, dexamethasone and hydrocortisone) but not retained, for the same thin-category reason given for CONMED_PPI. Glucocorticoids are CYP3A4 inducers with a documented voriconazole interaction, so the negative screen is informative rather than incidental."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 72L,
    n_studies      = 1L,
    n_observations = 150L,
    age_median     = "77 years (IQR 69-84)",
    weight_median  = "66.5 kg (IQR 60-70)",
    sex_female_pct = 25.0,
    race_ethnicity = c(Asian = 100),
    disease_state  = "Critically ill adults with COVID-19-associated pulmonary aspergillosis (CAPA) treated with intravenous voriconazole. 25 of 72 patients (34.7%) received continuous renal replacement therapy, uniformly in continuous veno-venous hemofiltration (CVVH) mode with blood flow 150-180 mL/min and dialysate flow 2 L/h. 45 of 72 (62.5%) were hypoalbuminaemic. The predominant pathogens were Aspergillus spp. (A. fumigatus, A. flavus, A. niger) and Candida spp.; 9 of 63 evaluable voriconazole courses were prophylactic. In-hospital mortality was 45.8%.",
    dose_range     = "Voriconazole given by intravenous infusion at approximately 4 mg/kg twice daily, with or without a 6 mg/kg loading dose; 50 of 72 patients received a loading dose. Administered doses ranged from 100 to 450 mg given once or twice daily, with the dosing interval, infusion rate and treatment course set by the treating team. Median treatment duration 9 days.",
    regions        = "Single center: Nanjing Drum Tower Hospital, Nanjing, Jiangsu, China.",
    renal_function = "Creatinine clearance (Cockcroft-Gault) median 54.97 mL/min, IQR 36.04-96.6; estimated glomerular filtration rate median 78.2 mL/min/1.73 m^2, IQR 48.39-125.6; serum creatinine median 84 umol/L, IQR 56.08-120.04. A third of the cohort was on continuous veno-venous hemofiltration.",
    notes          = "Retrospective single-center study of prospectively collected data, December 2022 to February 2023. Trough concentrations sampled at steady state after the fourth dose; most patients contributed two concentrations and 33 contributed one, for 150 concentration records in total. Observed troughs ranged from 0.15 to 11.0 mg/L, with 15.3% below 2 mg/L and 23.3% above 5 mg/L. Plasma HPLC with UV detection at 262 nm, calibration range 0.1-30 mg/L. NONMEM 7.3.0 with Pirana 2.9.0; final estimates and a 1000-sample nonparametric bootstrap (991 successful) per Table 2. CYP2C19 and CYP3A4 genotyping was performed on every participant but no genotype term survived covariate selection."
  )

  ini({
    # Structural parameters (Wang 2025 Table 2, "Final model" column).
    # These are ABSOLUTE population values: the paper carries no allometric
    # or other body-size term because over 30% of body weights were missing
    # and weight was never entered into the covariate screen (Discussion
    # limitation 2). The paper labels them CL/F and V/F, but the cohort
    # received intravenous voriconazole exclusively (Section 2.1 inclusion
    # criterion 3), so bioavailability is not identifiable here and the
    # "/F" is nominal - see the vignette Errata.
    lcl <- log(3.17); label("Clearance (L/h)")                      # Wang 2025 Table 2 (CL 3.17, RSE 6.1%, bootstrap median 3.176, 95% CI 2.76-3.57); base model 3.67
    lvc <- log(135);  label("Volume of distribution (L)")           # Wang 2025 Table 2 (V 135, RSE 22.1%, bootstrap median 130.919, 95% CI 54.37-200.50); base model 140

    # Covariate effects on clearance. All five come from the single
    # piecewise equation printed in Wang 2025 Section 3.2:
    #
    #   CL = 3.17 * (PLT/121)^0.248 * (CRP/55.23)^-0.183
    #             * (GGT/68.72)^0.2928 * (AST/29.91)^-0.227          , CRRT = 0
    #   CL = the same product * 1.617                                , CRRT = 1
    #
    # Two readings of that printed equation had to be settled.
    #
    # (1) The published equation carries a leading factor typeset as
    # e^0.0768, i.e. 3.17 * e^0.0768 * (the covariate product). Taken
    # literally that is a constant 1.0798-fold multiplier on the typical
    # value. It is NOT: 0.0768 is exactly the omega^2 for clearance in the
    # same Table 2, so the exponent slot holds the eta and the typeset
    # superscript has picked up the variance estimate. Three independent
    # statements in the paper falsify the literal reading and all agree
    # with exp(eta):
    #   - the Abstract: "The model estimated voriconazole's apparent
    #     clearance (CL/F) at 3.17 L/h ... for a standard patient";
    #   - the Discussion: "The standard values observed for the CL and V
    #     parameters in this investigation were 3.17 L/h and 135 L";
    #   - the Discussion again, decisively, on the CRRT arm: "The average
    #     CL of voriconazole under the influence of CRRT is approximately
    #     5.13 L/h". At median covariates every ratio is 1, so the two
    #     readings give 3.17 * 1.617 = 5.126 (matches 5.13) against
    #     3.17 * 1.0798 * 1.617 = 5.535 (does not).
    # The exp(eta) reading is therefore used and the e^0.0768 factor is not
    # carried as a constant.
    #
    # (2) Signs. Unlike the sibling Wang 2024 paper, Table 2 here prints
    # the coefficients WITH their signs (-0.183 for CRP, -0.227 for AST),
    # and the Section 3.2 equation agrees, so the two sources corroborate
    # rather than one supplying what the other omits. Both negative signs
    # also match the Discussion, which reads clearance as falling with
    # rising inflammation and rising transaminase.
    e_plt_cl  <-  0.248; label("Exponent of (PLT / 121 x 10^9/L) on clearance (unitless)")  # Wang 2025 Table 2, theta3 PLT on CL 0.248 (RSE 32.7%, bootstrap median 0.235, 95% CI 0.043-0.404); same value in the Section 3.2 equation
    e_crp_cl  <- -0.183; label("Exponent of (CRP / 55.23 mg/L) on clearance (unitless)")    # Wang 2025 Table 2, theta4 CRP on CL -0.183 (RSE 18.2%, bootstrap median -0.181, 95% CI -0.258 to -0.113); same value in the Section 3.2 equation
    e_ggt_cl  <-  0.292; label("Exponent of (GGT / 68.72 U/L) on clearance (unitless)")     # Wang 2025 Table 2, theta5 GGT on CL 0.292 (RSE 20.4%, bootstrap median 0.288, 95% CI 0.177-0.442); the Section 3.2 equation prints 0.2928
    e_ast_cl  <- -0.227; label("Exponent of (AST / 29.91 U/L) on clearance (unitless)")     # Wang 2025 Table 2, theta7 AST on CL -0.227 (RSE 34.2%, bootstrap median -0.224, 95% CI -0.419 to -0.078); same value in the Section 3.2 equation
    # CRRT enters as the fractional-change form theta1 * (1 + theta2) that
    # Wang 2025 Section 2.4.2 declares for categorical covariates, so the
    # stored coefficient is theta6 = 0.617 itself and the model block
    # applies (1 + 0.617 * CRRT). Section 3.2 writes out the resulting
    # 1.617 multiplier on the CRRT = 1 branch, which is the arithmetic
    # check on this encoding: 1 + 0.617 = 1.617.
    e_crrt_cl <-  0.617; label("Fractional increase in clearance while continuous renal replacement therapy is in use (unitless)") # Wang 2025 Table 2, theta6 CRRT on CL 0.617 (RSE 31.9%, bootstrap median 0.57, 95% CI 0.232-1.105); Section 3.2 equation multiplier 1.617

    # IIV. Wang 2025 Section 2.4.1 specifies exponential interindividual
    # variability. Table 2 heads the block "Interindividual variability"
    # and its footnote states that "omega^2CL, omega^2V are variance
    # estimates of the interindividual variability of CL, and V", so
    # 0.0768 is a variance (CV 28.5%) and is used directly, not squared.
    etalcl ~ 0.0768 # Wang 2025 Table 2, omega^2CL 0.0768 (RSE 34.4%, bootstrap median 0.0678, 95% CI 0.018-0.126); base model 0.212
    # No eta on V. Wang 2025 Table 2 reports omega^2V as "0 FIX" in both
    # the base and the final model, and Section 3.2 explains why: "In our
    # basic model, the shrinkage of V was high (58%), which means the V
    # parameter was not significant estimate inter individual variability.
    # Hence, we fixed the inter individual variability of V as zero." A
    # zero-variance eta is mechanically identical to no eta and would make
    # OMEGA singular, so it is omitted rather than written as ~ fixed(0).

    # Residual error. Wang 2025 Section 3.2 states that "the residual error
    # was best characterized by an addictive [sic] error model", so the
    # model is additive only with no proportional component. The Table 2
    # value 1.87 is a VARIANCE, not a standard deviation, so the standard
    # deviation nlmixr2 wants is sqrt(1.87) = 1.368 mg/L. Two independent
    # supports: the Table 2 footnote says "Additive error is the variance
    # estimate of the variance of the summed residual variance", and the
    # same table reports its interindividual variability rows as variances
    # under the omega^2 symbol, so the whole variability block is on the
    # variance scale. A magnitude check agrees: the 150 observed troughs
    # have median 3.6 mg/L and IQR 2.5-5 mg/L, an interquartile width
    # implying a total observation standard deviation near 1.85 mg/L, which
    # a 1.368 mg/L residual leaves room for once the 28.5% clearance IIV,
    # the covariate spread and the 100-450 mg dose range are added, but
    # which a 1.87 mg/L residual would already exceed on its own.
    addSd <- sqrt(1.87); label("Additive residual error (mg/L)")    # Wang 2025 Table 2, Additive error 1.87 (RSE 16.4%, bootstrap median 1.771, 95% CI 1.212-2.444), reported as a variance; base model 2.3
  })

  model({
    # Clearance: the Wang 2025 Section 3.2 equation. The four continuous
    # covariates keep the paper's power form on a median-normalised ratio;
    # CRRT keeps the paper's declared categorical fractional-change form.
    cl <- exp(lcl + etalcl) *
      (PLT / 121)^e_plt_cl *
      (CRP / 55.23)^e_crp_cl *
      (GGT / 68.72)^e_ggt_cl *
      (AST / 29.91)^e_ast_cl *
      (1 + e_crrt_cl * RRT_CRRT_STATUS)

    vc <- exp(lvc)

    kel <- cl / vc

    # One-compartment disposition with first-order elimination (Wang 2025
    # Section 2.4.1 and Section 3.2). Dosing is intravenous infusion
    # directly into central; there is no depot and no bioavailability term
    # because the cohort received no extravascular voriconazole.
    d/dt(central) <- -kel * central

    # Dose in mg divided by volume in L gives mg/L, the unit of the
    # reported plasma concentrations (equivalently ug/mL).
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
