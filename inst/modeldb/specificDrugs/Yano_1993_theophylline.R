Yano_1993_theophylline <- function() {
  description <- "One-compartment IV-infusion population PK model for theophylline (Yano 1993 Paper II) in 55 adult inpatients with stable chronic airway obstruction; clearance and volume of distribution are log-linear functions of arterial PaCO2 and a binary hepatic-dysfunction indicator."
  reference <- "Yano I, Tanigawara Y, Yasuhara M, Okumura K, Kawakatsu K, Nishimura K, Hori R. Population Pharmacokinetics of Theophylline. II: Intravenous Infusion to Patients with Stable Chronic Airway Obstruction. Biol Pharm Bull. 1993;16(5):501-505. doi:10.1248/bpb.16.501"
  vignette <- "Yano_1993_theophylline"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight (baseline)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Yano 1993 first assessed body-size scaling via SIZE = WT^theta1 * HT^theta2 (Eq. 3, page 503) and rejected H0: theta1 = 1 and H0: theta2 = 0 at p < 0.005 in favour of theta1 = 1, theta2 = 0. Thereafter CL (L/h/kg) and Vd (L/kg) are reported per kg of body weight (Table III, page 504), so simulation conversion to total clearance and total volume uses WT linearly with no allometric exponent (i.e., (WT/ref)^1 with no reference weight -- the per-kg parameterization is directly multiplicative).",
      source_name        = "WT"
    ),
    PACO2 = list(
      description        = "Arterial blood carbon dioxide tension (PaCO2)",
      units              = "mmHg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Measured at study admission by an automatic blood gas analyzer (ABL-3, Radiometer Copenhagen) on an arterial puncture sample. Cohort mean 42.5 +/- 6.8 mmHg, range 31-71 (Yano 1993 Table I, page 502). Enters both CL and Vd as a log-linear effect (Yano 1993 Eqs. 7-8, page 504). Treated as time-fixed at study admission.",
      source_name        = "PaCO2"
    ),
    HEPIMP = list(
      description        = "Binary hepatic-dysfunction indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Yano 1993 defines HF = 1 if a patient was classified as having hepatic dysfunction by their physician OR if serum glutamic pyruvic transaminase (ALT) or serum glutamic oxaloacetic transaminase (AST) exceeded 30 IU/L (Yano 1993 Patient Description, page 501). 5 of 55 patients had HF = 1. This is a paper-specific clinician-judgment + transaminase-cutoff scheme distinct from the NCI ODWG classification used by oncology-trial popPK papers (see inst/references/covariate-columns.md HEPIMP entry for the registered schemes). No patient had both hepatic dysfunction and a smoking habit, so HEPIMP is independent of the screened SMOKE covariate in this cohort.",
      source_name        = "HF"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age (years)",
      units       = "years",
      type        = "continuous",
      notes       = "Screened in Yano 1993 Table II forward selection (page 503): significant on CL at the first step (p < 0.025) and second step (p < 0.01) but lost significance by the third step (p = n.s.) once PaCO2 had entered. Not retained in the final regression model. Vignette narrative cites the discussion of age and smoking confounding (page 504)."
    ),
    SMOKE = list(
      description = "Current-smoker indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Yano 1993 SM column: 1 if patient smoked 20-60 cigarettes/day for 25-58 years (average 43 years) up to or within 6 months of the study (page 502). 11 of 55 patients were smokers. Screened on CL (1st step likelihood ratio not significant; -2 l.l.d. = 2.661) and not retained. The paper notes smoking is multivariate-confounded with age and the smokers in this cohort were geriatric; smoking did not reach significance even unadjusted (Discussion page 504)."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Yano 1993 Alb column. Screened on both CL and Vd at the first selection step (Table II); not significant on either. Cohort mean 4.2 +/- 0.5 g/dL, range 3.1-5.1 (Table I)."
    ),
    HCT = list(
      description = "Hematocrit",
      units       = "%",
      type        = "continuous",
      notes       = "Yano 1993 Hct column. Screened on both CL and Vd at the first selection step (Table II); not significant on either. Cohort mean 41.6 +/- 4.9 %, range 32-52 (Table I)."
    ),
    PaO2 = list(
      description = "Arterial blood oxygen tension",
      units       = "mmHg",
      type        = "continuous",
      notes       = "Yano 1993 PaO2 column. Screened on CL at the first selection step (Table II; -2 l.l.d. = 7.390, p < 0.01, but eliminated at the third step when PaCO2 entered). Not retained in the final regression model. Cohort mean 71.4 +/- 11.4 mmHg, range 52-108 (Table I). NOT a registered canonical covariate; documented here only because Yano 1993 screened it. Future popPK papers retaining PaO2 as a final covariate should propose a canonical (e.g., PAO2) at extraction time."
    ),
    BLOOD_PH = list(
      description = "Arterial blood pH",
      units       = "pH units",
      type        = "continuous",
      notes       = "Yano 1993 pH column. Screened on both CL and Vd in Table II (page 503): not significant on CL; significant on Vd at the first three selection steps (peaked at -2 l.l.d. = 8.539, p < 0.001) but lost significance at the fourth step when PaCO2 had entered Vd. Not retained in the final regression model. Cohort mean 7.402 +/- 0.029, range 7.300-7.483 (Table I). NOT a registered canonical covariate; documented here only because Yano 1993 screened it. Future popPK papers retaining blood pH as a final covariate should propose a canonical at extraction time."
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 55L,
    n_studies       = 1L,
    age_range       = "22-80 years",
    age_mean        = "59.0 years (SD 11.6)",
    weight_range    = "28-77 kg",
    weight_mean     = "49.3 kg (SD 11.6)",
    height_range    = "136-182 cm",
    height_mean     = "158.3 cm (SD 10.0)",
    sex_female_pct  = 36.4,
    race_ethnicity  = c(Japanese = 100),
    disease_state   = "Adult inpatients with stable chronic airway obstruction (chronic asthma 12, chronic pulmonary emphysema 20, diffuse panbronchiolitis 25, chronic bronchitis 2, other 1; 2 patients had two disease types) of mild to moderate severity. All were inpatients at the Chest Disease Research Institute, Kyoto University, initially treated with IV aminophylline for an asthmatic episode and continued on oral theophylline after responding (Yano 1993 Patient Description, page 501).",
    dose_range      = "Single IV short infusion of aminophylline 250 mg (= 200 mg theophylline) in 100 mL saline or 250 mL 5% glucose. Infusion duration 13-120 min (mean 0.517 +/- 0.250 h, range 0.22-2.0 h). Infusion rate 9.21 +/- 3.36 mg/h/kg, range 2.17-20.41 (Yano 1993 Table I, page 502).",
    regions         = "Japan (Kyoto University Hospital, Chest Disease Research Institute).",
    n_observations  = 276L,
    blood_gas       = "PaO2 71.4 +/- 11.4 mmHg, PaCO2 42.5 +/- 6.8 mmHg, blood pH 7.402 +/- 0.029, Hct 41.6 +/- 4.9 %, Alb 4.2 +/- 0.5 g/dL (Yano 1993 Table I, page 502).",
    hepatic_function = "5 of 55 patients with hepatic dysfunction by clinician judgment OR AST/ALT > 30 IU/L (Yano 1993 page 501).",
    smoking         = "11 of 55 current smokers; smokers were geriatric (mean smoking history 43 years).",
    notes           = "Serum samples drawn before infusion, immediately after infusion end, and up to 10 h after infusion end; 3-7 measurements/patient, 276 total. Theophylline assayed by fluorescence polarization immunoassay (TDx, Abbott) or HPLC (within-day and between-day CV < 3 % at 5, 10, 20 ug/mL; Kawakatsu 1989 reference 16). NONMEM III with likelihood-ratio test for stepwise forward covariate selection; final structural model was one-compartment IV-infusion because the data showed virtually mono-exponential decline (Figure 2, page 502) despite theophylline's previously-described two-compartment behavior (Mitenko & Ogilvie 1972 reference 18)."
  )

  ini({
    # Structural parameters -- log-scale clearance and volume at PaCO2 = 0 mmHg and HEPIMP = 0 (no hepatic dysfunction).
    # These are the published P1, P2 estimates from Yano 1993 Table III (page 504) substituted directly into Eqs. 7-8
    # (page 504): CL_j = exp(P1 + P3*HF + P5*PaCO2) [L/h/kg] and Vd_j = exp(P2 + P4*PaCO2) [L/kg].
    # The published "typical" values at PaCO2 = 42.5 mmHg, HEPIMP = 0 are
    #   CL = exp(-3.78 + 0.0233 * 42.5) = exp(-2.79) = 0.0614 L/h/kg = 61.4 mL/h/kg
    #   Vd = exp(-1.12 + 0.00934 * 42.5) = exp(-0.723) = 0.485 L/kg
    # which match the paper's reported typical values (Yano 1993 page 504).
    lcl <- -3.78; label("Log clearance per kg at PaCO2 = 0 mmHg, HEPIMP = 0 (log L/h/kg)")   # Yano 1993 Table III P1 = -3.78 (95% CI -4.13, -3.43)
    lvc <- -1.12; label("Log central volume of distribution per kg at PaCO2 = 0 mmHg (log L/kg)")  # Yano 1993 Table III P2 = -1.12 (95% CI -1.34, -0.896)

    # Covariate effects -- additive on the log scale because the published Eqs. 7-8 enter HEPIMP and PaCO2
    # inside the exponential. Direction: hepatic dysfunction lowers CL (P3 < 0); higher PaCO2 raises CL and Vd
    # (P5 > 0, P4 > 0).
    e_hepimp_cl <- -0.525;  label("Log-scale effect of hepatic dysfunction on CL")     # Yano 1993 Table III P3 = -0.525 (95% CI -0.897, -0.153)
    e_paco2_cl  <- 0.0233;  label("Log-scale effect of PaCO2 on CL (per mmHg)")        # Yano 1993 Table III P5 = 0.0233 (95% CI 0.0158, 0.0308)
    e_paco2_vc  <- 0.00934; label("Log-scale effect of PaCO2 on Vd (per mmHg)")        # Yano 1993 Table III P4 = 0.00934 (95% CI 0.00430, 0.0144)

    # IIV -- proportional error model on the individual parameters (Eq. 2, page 502): CL_j = CL_typ*(1 + eta_CL_j),
    # V_d,j = V_d_typ*(1 + eta_Vd_j). Yano 1993 reports inter-individual variabilities as CV %: 38.5 % on CL and
    # 12.5 % on Vd (Table III, page 504; 95% CIs 27.5-46.9 % and 8.7-15.4 % respectively). Translated to log-scale
    # variances via omega^2 = log(1 + CV^2): CL -> log(1 + 0.385^2) ~ 0.1383; Vd -> log(1 + 0.125^2) ~ 0.01552.
    # No correlation is reported between eta_CL and eta_Vd; encoded as diagonal.
    etalcl ~ 0.1383
    etalvc ~ 0.01552

    # Residual error -- proportional error model (Eq. 2, page 502): C_ij = C_ij_pred*(1 + eps_ij). Yano 1993
    # reports residual variability as 10.6 % CV (Table III, page 504; 95% CI 8.5-12.3 %).
    propSd <- 0.106; label("Proportional residual error (fraction)")                   # Yano 1993 Table III: residual variability 10.6 % CV (95% CI 8.5-12.3)
  })

  model({
    # Individual log-CL and log-Vd per kg of body weight (Yano 1993 Eqs. 7-8, page 504).
    cl_per_kg <- exp(lcl + e_hepimp_cl * HEPIMP + e_paco2_cl * PACO2 + etalcl)  # L/h/kg
    vc_per_kg <- exp(lvc + e_paco2_vc * PACO2 + etalvc)                          # L/kg

    # Convert to total clearance and total central volume by linear scaling on body weight
    # (Yano 1993 Eq. 4 with theta1 = 1, page 503).
    cl <- cl_per_kg * WT  # L/h
    vc <- vc_per_kg * WT  # L

    kel <- cl / vc

    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
