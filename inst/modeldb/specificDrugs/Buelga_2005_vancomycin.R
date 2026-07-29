Buelga_2005_vancomycin <- function() {
  description <- "One-compartment IV intermittent-infusion population PK model for vancomycin in adult patients with hematological malignancies (Buelga 2005). CL is a purely multiplicative function of Cockcroft-Gault creatinine clearance (CL [L/h] = 1.08 x CLCR [L/h]) and V is a purely multiplicative function of total body weight (V [L] = 0.98 x TBW [kg]). Exponential inter-individual variability on CL and V with an estimated CL-V correlation; additive residual error in mg/L. The AML-1 and AML-2 subpopulation-specific models from the same paper are not packaged here; only the general final model (Table 4) is implemented."
  reference <- "Buelga DS, Fernandez de Gatta MM, Herrera EV, Dominguez-Gil A, Garcia MJ. Population pharmacokinetic analysis of vancomycin in patients with hematological malignancies. Antimicrob Agents Chemother. 2005;49(12):4934-4941. doi:10.1128/AAC.49.12.4934-4941.2005"
  vignette <- "Buelga_2005_vancomycin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Buelga 2005 Table 1: mean 64.7 kg, SD 11.3 (index set n = 215). Source paper symbol TBW (total body weight); stored under the canonical WT column. Buelga 2005 final model (Table 4): V (L) = theta2 * TBW with theta2 = 0.98 L/kg -- a purely multiplicative effect on V with no centering and no allometric exponent.",
      source_name        = "TBW"
    ),
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Buelga 2005 Table 1: mean 89.4 mL/min, SD 39.2 (index set n = 215). Estimated by the Cockcroft-Gault equation (Cockcroft and Gault 1976; Buelga 2005 Methods 'Data acquisition' and final-model paragraph in Results). The Cockcroft-Gault, Jelliffe, Tsubaki, and Levey formulas were compared; Cockcroft-Gault was selected for the final model. Stored under the canonical CRCL column with units mL/min (raw Cockcroft-Gault), matching the Goti_2018_vancomycin.R, Moore_2016_vancomycin.R, and Delattre_2010_amikacin.R precedents in inst/references/covariate-columns.md. Buelga 2005 expresses the CL covariate equation with CLCR in L/h (Buelga 2005 abstract and final-model paragraph: 'CL (liters/h) = 1.08 * CLCR(Cockcroft and Gault) (liters/h)'); the packaged model stores the column in mL/min (canonical) and converts to L/h inside model() by multiplying by 60/1000.",
      source_name        = "CLCR"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 215L,
    n_studies        = 1L,
    age_range        = ">=15 years (mean 51.5, SD 15.9)",
    age_median       = "51.5 years (mean, SD 15.9)",
    weight_range     = "mean 64.7 kg (SD 11.3)",
    weight_median    = "64.7 kg (mean, SD 11.3)",
    sex_female_pct   = 44.7,
    race_ethnicity   = "Not reported (single-center cohort at University Hospital of Salamanca, Spain; population presumed predominantly Spanish/European)",
    disease_state    = "Adult (>=15 years) inpatients with underlying hematological malignancy admitted to the Hematology Unit for suspected or documented gram-positive bacterial infection. Hematological diagnoses: AML 27.6%, NHL 30.8%, ALL 8.0%, CML 8.0%, CLL 4.2%, Hodgkin's 7.7%, MDS 4.5%, multiple myeloma 4.2%, other 4.9%. 15.7% had autologous bone marrow transplant; 43.7% had neutropenia (ANC < 500/mm^3); 45.1% had ECOG 0-2 and 12.9% had ECOG 3-4; 38.8% received concomitant amikacin and 21.0% concomitant amphotericin. ICU patients were excluded.",
    dose_range       = "Intravenous vancomycin by intermittent infusion over 0.5-1 h. Daily dose 200-3,900 mg/day (mean 1,535 +/- 280); dosing interval 6-48 h. Doses individualized by physician via a hematology-specific nomogram considering age, weight, and renal function.",
    regions          = "Spain (University Hospital of Salamanca)",
    renal_function   = "Cockcroft-Gault CRCL mean 89.4 mL/min (SD 39.2); serum creatinine mean 0.9 mg/dL (SD 0.4).",
    n_concentrations = 1004L,
    notes            = "Demographics from Buelga 2005 Table 1 (index set; n=215 patients, 1004 vancomycin serum concentrations). Retrospective collection of therapeutic-drug-monitoring (TDM) data from 1989-1999. 274 patients (348 treatment courses) were initially available; 11 excluded for incomplete data and 6 for ICU stay during therapy. The 215 index-set patients include 119 male / 96 female. Additional 59-patient validation cohort (40 male / 19 female; 124 concentrations) used for external validation; the validation set had a higher proportion of AML (40.3% vs 27.6%) and females (32% vs 44.7%) than the index set. Vancomycin assayed by fluorescence polarization immunoassay (TDx, Abbott; LLOQ 0.6 mg/L; CV <5% over 7-75 mg/L). Modeling done in NONMEM V (double precision, level 1.1) with first-order conditional estimation. Sampling was 'one or more post-distribution samples after the first doses or at steady state'; peak samples drawn at least 2 h after end of infusion, justifying the one-compartment model. Buelga 2005 also developed AML-specific subpopulation models (AML-1 with TBW, SCR, age, sex covariates on CL; AML-2 with the same CL = theta*CLCR / V = theta*TBW structure but with a 10%-higher CL coefficient of 1.17). The AML-specific variants are documented in the vignette's Assumptions and deviations section but not packaged as separate models; only the general final model from Buelga 2005 Table 4 is implemented here."
  )

  ini({
    # Structural parameters (Buelga 2005 Table 4 final-model column "Estimate (%estimation error)").
    # The general final model defines:
    #   CL (L/h) = theta1 * CLCR (L/h);  theta1 = 1.08 (%est. err. 2.12%)
    #   V (L)   = theta2 * TBW (kg);     theta2 = 0.98 L/kg (%est. err. 7.43%)
    # No centering; pure multiplicative covariate effects with no additional intercept.
    lcl <- log(1.08); label("CL coefficient on CLCR (L/h per L/h CLCR; dimensionless ratio)") # Buelga 2005 Table 4: theta1 = 1.08 (CL / CLCR)
    lvc <- log(0.98); label("V coefficient on TBW (L per kg TBW)")                            # Buelga 2005 Table 4: theta2 = 0.98 (V / TBW)

    # Inter-individual variability (Buelga 2005 Table 4 "omega" rows reported as CV in %).
    # Methods: 'theta_j = theta' * exp(eta_Theta_j)' (exponential IIV model); the FOCE
    # first-order method approximates the exponential error model as a proportional
    # error model. For log-normal eta with CV%, omega^2 = log(CV^2 + 1):
    #   etalcl: 28.16% CV -> omega^2 = log(0.2816^2 + 1) = log(1.07930) = 0.07631
    #   etalvc: 37.15% CV -> omega^2 = log(0.3715^2 + 1) = log(1.13800) = 0.12930
    # Buelga 2005 Table 4 also reports the CL-V random-effect correlation as
    # 'omega CL/V = 23.12%', which is interpreted here as the correlation coefficient
    # between etalcl and etalvc (rho = 0.2312); the off-diagonal covariance is
    #   cov = rho * sqrt(omega^2_lcl * omega^2_lvc)
    #       = 0.2312 * sqrt(0.07631 * 0.12930) = 0.2312 * 0.09931 = 0.02296
    # See vignette Assumptions and deviations for the discussion of this convention.
    # Block matrix order is var(etalcl), cov(etalcl, etalvc), var(etalvc).
    # Buelga 2005 Table 4: omega_CL = 28.16% CV; omega_V = 37.15% CV; omega_CL/V = 23.12% (correlation).
    etalcl + etalvc ~ c(0.07631,
                        0.02296, 0.12930)

    # Additive residual error (Buelga 2005 Methods: 'C_ij = C'_ij + epsilon' with
    # zero-mean variance sigma^2; Table 4: sigma = 3.52 mg/L, %est. err. 15.12%).
    addSd <- 3.52; label("Additive residual error (mg/L)") # Buelga 2005 Table 4: sigma = 3.52 mg/L
  })
  model({
    # Individual PK parameters.
    # CRCL is stored in mL/min (canonical CRCL units); Buelga 2005 expresses the CL
    # covariate in L/h. Convert with the factor 60/1000 (= 0.06) so the typical
    # subject with CRCL = 90 mL/min (5.4 L/h) and TBW = 65 kg has:
    #   CL = 1.08 * 5.4 = 5.832 L/h; V = 0.98 * 65 = 63.7 L; kel = 0.0916 1/h; t1/2 ~ 7.6 h
    crcl_Lh <- CRCL * 60 / 1000
    cl <- exp(lcl + etalcl) * crcl_Lh
    vc <- exp(lvc + etalvc) * WT

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, vc in L -> central/vc has units mg/L.
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
