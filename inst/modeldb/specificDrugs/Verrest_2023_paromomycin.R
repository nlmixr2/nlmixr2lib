Verrest_2023_paromomycin <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption for",
    "intramuscular paromomycin in 265 Eastern African children and adults",
    "with visceral leishmaniasis (Verrest 2023), enrolled in the phase III",
    "randomized controlled trial NCT03129646 across six sites in Kenya,",
    "Sudan, Ethiopia and Uganda and treated with paromomycin sulphate",
    "20 mg/kg/day intramuscularly (equivalent to 15 mg/kg/day paromomycin",
    "base) for 14 days, combined with allometrically dosed oral miltefosine",
    "for 14 or 28 days. CL/F, Q/F, Vc/F and Vp/F are allometrically scaled",
    "on baseline body weight (exponents fixed at 0.75 and 1.0; reference",
    "27.5 kg, the cohort median). Relative bioavailability is fixed at 1.17",
    "for comparability with the earlier Kenyan and Sudanese monotherapy",
    "estimates. The clearance decrease observed over the treatment course",
    "is modelled mechanistically as a linear function of the time-varying",
    "absolute neutrophil count, which recovers from neutropenia as the",
    "disease resolves: clearance falls by 13% per 1 x 10^3 cells/uL rise in",
    "neutrophils above the population median of 0.98 x 10^3 cells/uL.",
    "Between-subject variability is log-normal on CL/F only (56.1% CV);",
    "residual error is proportional (53.4% CV). The companion miltefosine",
    "model from the same trial is Verrest_2023_miltefosine."
  )
  reference <- paste(
    "Verrest L, Roseboom IC, Wasunna M, Mbui J, Njenga S, Musa AM, Olobo J,",
    "Mohammed R, Ritmeijer K, Chu WY, Huitema ADR, Solomos A, Alves F,",
    "Dorlo TPC. Population pharmacokinetics of a combination of miltefosine",
    "and paromomycin in Eastern African children and adults with visceral",
    "leishmaniasis. J Antimicrob Chemother. 2023;78(11):2702-2714.",
    "doi:10.1093/jac/dkad286. ClinicalTrials.gov NCT03129646.",
    sep = " "
  )
  vignette <- "Verrest_2023_miltefosine_paromomycin"
  units    <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Amounts are mg of paromomycin BASE, matching the
  # trial dataset's AMT column ("PM: mg PM base (DOSE*375 mg/ml)",
  # supplementary NONMEM control stream $INPUT dataset description).
  compartmentData <- list(
    depot       = list(analyte = "paromomycin", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "paromomycin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "paromomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT_BASE = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric power scaling on CL/F and Q/F (exponent fixed 0.75)",
        "and on Vc/F and Vp/F (exponent fixed 1.00), normalized to the",
        "cohort median of 27.5 kg (Verrest 2023 Table 3 legend: 'WTmed,",
        "median population body weight (27.5 kg)'). Time-FIXED at the",
        "subject's baseline value: the supplementary NONMEM control",
        "stream writes the allometry against the dataset's BWT column",
        "('BWT: baseline body weight (kg)'), not the time-varying WT",
        "column that the same dataset also carries. The Table 3 equation",
        "prints the subscript as WT_i,t, which would imply a",
        "time-varying weight; the control stream is the executed model",
        "and is followed here. Cohort baseline weight was 32.9 kg mean",
        "(range 11.0-71.0) in the sparsely sampled group and 33.6 kg",
        "(20.5-54.0) in the intensively sampled group (Table 1)."
      ),
      source_name        = "BWT"
    ),
    NEUT = list(
      description        = "Absolute neutrophil count",
      units              = "cells/uL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "TIME-VARYING, and the mechanism by which this model captures the",
        "decrease in paromomycin clearance over the treatment course.",
        "Enters CL/F as the centred linear term",
        "(1 + e_neut_cl * (NEUT/1000 - 0.98)), where the median",
        "population neutrophil count is 0.98 x 10^3 cells/uL (Verrest",
        "2023 Table 3 legend). NOTE the unit change: the paper reports",
        "neutrophils in 10^3 cells/uL and its slope of -0.13 is 'per",
        "10^3 cells/uL', while the canonical NEUT column is in cells/uL",
        "(equivalently cells/mm^3); model() therefore divides by 1000",
        "before centring so that a standard clinical NEUT column can be",
        "supplied unchanged. Patients were neutropenic at treatment",
        "start (IQR 0.74-1.38 x 10^3 cells/uL) and recovered during",
        "treatment (Results, 'Patients and data'; Figure S3, where",
        "levels are interpolated between Day 0 and Day 28)."
      ),
      source_name        = "NEUTR"
    )
  )

  covariatesDataExcluded <- list(
    EGFR_ABS = list(
      description        = "Absolute (non-BSA-normalized) estimated glomerular filtration rate",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened on paromomycin clearance because paromomycin is",
        "cleared mainly renally, both as a fraction of eGFRabs and in",
        "combination with a non-renal clearance route, and NOT retained:",
        "'eGFRabs could not explain variability in clearance' (Verrest",
        "2023 Results). Derived per Eq. 3 as the BSA-unadjusted eGFR,",
        "with CKD-EPI (without the ethnicity adjustment) for patients",
        ">14 years and the Schwartz formula for children <=14 years. The",
        "Discussion attributes the null result to serum creatinine",
        "overestimating renal function in a malnourished, low-muscle-mass",
        "population. Cohort value 68.9 mL/min mean (range 14.9-288.1)."
      ),
      source_name        = "EGFRA"
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened for an inverse correlation with paromomycin volume of",
        "distribution, by analogy with other aminoglycosides in patients",
        "with haematological malignancies, and not retained in the final",
        "model (Verrest 2023 Methods 'Covariate analysis'; Table 3 lists",
        "no albumin term). Cohort value 28.3 g/L mean (range 1.5-54.6);",
        "levels were low at treatment start (IQR 23.4-32.4 g/L) and rose",
        "during treatment."
      ),
      source_name        = "ALB"
    ),
    AGE = list(
      description        = "Age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened and not retained: 'Age or country of origin could not",
        "explain remaining variability in any of the pharmacokinetic",
        "parameters on top of the identified covariates' (Verrest 2023",
        "Results). The supplementary control stream still multiplies CL",
        "by vestigial CLAGE and CLTIME terms that the stream never",
        "defines; both are inert (= 1) in the final model printed in",
        "Table 3. Cohort age 13.6 years mean (range 4-45)."
      ),
      source_name        = "AGE"
    ),
    CNTRY = list(
      description        = "Country of enrolment (Ethiopia, Kenya, Sudan, Uganda)",
      units              = "(categorical)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "A remaining population difference between countries was",
        "evaluated on paromomycin bioavailability, absorption rate,",
        "volume of distribution and clearance, and none was retained:",
        "'no significant pharmacokinetic differences between countries",
        "were identified, indicating that there are no geographical",
        "differences that are not already explained by demographic",
        "differences between populations or other covariates' (Verrest",
        "2023 Discussion)."
      ),
      source_name        = "CNTRY"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 265L,
    n_studies      = 1L,
    age_range      = "4-45 years",
    age_mean       = "13.6 years (sparse sampling group); 13.8 years (intensive sampling group)",
    weight_range   = "11.0-71.0 kg",
    weight_mean    = "32.9 kg (sparse sampling group); 33.6 kg (intensive sampling group)",
    sex_female_pct = 19.2,
    race_ethnicity = "Eastern African (Kenyan, Sudanese, Ethiopian and Ugandan cohorts; no further breakdown reported)",
    disease_state  = paste(
      "Eastern African children and adults with symptomatic,",
      "parasitologically confirmed visceral leishmaniasis; 59% were",
      "paediatric (<=12 years). Patients with relapse, severe",
      "malnutrition, severe VL, HIV co-infection or concomitant severe",
      "infection were excluded. The cohort was neutropenic (IQR",
      "0.74-1.38 x 10^3 cells/uL) and hypoalbuminaemic (IQR 23.4-32.4",
      "g/L) at treatment start, with both recovering during treatment."
    ),
    renal_function = paste(
      "Absolute eGFR 68.9 mL/min mean (range 14.9-288.1); serum",
      "creatinine 0.8 mg/dL mean (range 0.1-1.5). Two of 268 trial",
      "patients developed acute kidney injury; one reached a serum",
      "creatinine of 10.8 mg/dL by Day 14 and an AUC0-24 of 2388",
      "ug*h/mL, causing bilateral deafness."
    ),
    dose_range     = paste(
      "Paromomycin sulphate 20 mg/kg/day as a once-daily intramuscular",
      "injection, equivalent to 15 mg/kg/day paromomycin base, for 14",
      "days. Amounts in this model are mg of paromomycin BASE, matching",
      "the trial dataset. All patients also received oral miltefosine",
      "(allometric dose) twice daily for 14 or 28 days."
    ),
    regions        = "Eastern Africa (Kenya: Kacheliba; Uganda: Amudat; Sudan: Doka, Um El Kher; Ethiopia: Gondar, Abdurafi)",
    co_medication  = paste(
      "Oral miltefosine given simultaneously in every patient, for 14",
      "days (PM+MF14D arm) or 28 days (PM+MF28D arm). The exposures",
      "achieved for both drugs matched previous monotherapy studies,",
      "which the authors read as evidence against a drug-drug",
      "interaction; no interaction term is present in this model. The",
      "companion miltefosine model is Verrest_2023_miltefosine."
    ),
    samples        = paste(
      "229 paromomycin plasma concentrations from the 26 patients in the",
      "intensive-sampling cohort (Kenya and Sudan), sampled at 0, 1, 2,",
      "4 and 24 h or 0, 1, 2, 8 and 24 h on Day 1 and Day 14. Three",
      "observations were excluded as unreliable, leaving 232 in the",
      "analysis dataset per the Results text and 229 per Table 2. No",
      "observations were below the 5 ng/mL limit of quantification.",
      "Patients outside the intensive cohort contributed no paromomycin",
      "samples."
    ),
    notes          = paste(
      "NONMEM 7.5, ADVAN13, FOCE-I with interaction; parameter precision",
      "by sampling importance resampling. The starting point was the",
      "earlier paromomycin model in Kenyan and Sudanese VL patients",
      "(Verrest 2023 reference 8). Demographics are reported in Table 1",
      "as mean (range) at baseline; the sex percentage recorded here is",
      "derived from the 193 of 239 sparsely sampled patients reported as",
      "male. This study predominantly involved male patients, which the",
      "authors note may confound any sex-related effects."
    )
  )

  ini({
    # ============================================================
    # Structural PK parameters -- Verrest 2023 Table 3, 'Population
    # parameters' block. Values are apparent (CL/F etc.) and, per Table
    # 3 footnote b, are reported 'relative to a bioavailability of
    # 1.17'. The control stream (supplementary 'NONMEM control stream
    # paromomycine.docx') carries TIME in hours, so every rate constant
    # here is per hour; this matches the units printed in Table 3.
    # ============================================================
    lcl <- log(2.62)
    label("Apparent clearance CL/F at 27.5 kg and median neutrophils (L/h)") # Table 3 'CL (L/h)': 2.62 (95% CI 2.08-3.25); control stream $THETA(1)
    lvc <- log(9.17)
    label("Apparent central volume Vc/F at 27.5 kg (L)")                     # Table 3 'Vc (L)': 9.17 (95% CI 8.42-9.92); control stream $THETA(2)
    lq  <- log(0.26)
    label("Apparent inter-compartmental clearance Q/F at 27.5 kg (L/h)")     # Table 3 'Q': 0.26 (95% CI 0.22-0.31); units L/h from control stream $THETA(5) '0.262 ; 5 Q (L/hr)'
    lvp <- log(6.55)
    label("Apparent peripheral volume Vp/F at 27.5 kg (L)")                  # Table 3 'Vp (L)': 6.55 (95% CI 4.58-9.58); control stream $THETA(6)
    lka <- log(2.05)
    label("First-order absorption rate constant ka from the IM site (1/h)")  # Table 3 'ka (h-1)': 2.05 (95% CI 1.55-2.92); control stream $THETA(3)

    # Relative bioavailability was not estimated: it was FIXED to 1.17
    # so that the apparent parameters remain comparable with the earlier
    # Sudanese / Kenyan monotherapy estimates (Methods, 'Population
    # pharmacokinetic analysis').
    lfdepot <- fixed(log(1.17))
    label("Relative bioavailability F1 (unitless)")                   # Table 3 'F1': 1.17 (fixed); control stream $THETA(4) '(1.17) FIX'

    # ============================================================
    # Covariate effects.
    # ============================================================
    # Allometric exponents were not estimated; they were fixed at the
    # standard theory-based values (Methods, 'Covariate analysis':
    # 'fixed allometric exponents of 0.75 for clearance, and 1 for
    # volume of distribution'). The control stream applies the 0.75
    # exponent to Q as well as CL, and the 1.00 exponent to Vp as well
    # as Vc, which the Table 3 equation abbreviates as CL_TV and Vd_TV.
    e_wt_base_cl_q  <- fixed(0.75)
    label("Allometric exponent on baseline weight shared by CL/F and Q/F (unitless)")  # Table 3 equation exponent 0.75; control stream 'CL = TVCL *(BWT/27.5)**0.75', 'Q = TVQ *(BWT/27.5)**0.75'
    e_wt_base_vc_vp <- fixed(1.00)
    label("Allometric exponent on baseline weight shared by Vc/F and Vp/F (unitless)") # Table 3 equation exponent 1.00; control stream 'V2 = TVV2 *(BWT/27.5)**1.00', 'V3 = TVV3 *(BWT/27.5)**1.00'

    # Linear, centred neutrophil effect on clearance -- the only
    # covariate retained beyond allometry. Sign is negative: clearance
    # FALLS as neutrophils recover.
    e_neut_cl <- -0.13
    label("Fractional change in CL/F per 10^3 cells/uL rise in neutrophils above 0.98") # Table 3 'COV CL,neutr (fractional change/10^3 cells/uL)': -0.13 (95% CI -0.16 to -0.10); footnote c 'CL decreases by 13% per 1 x 10^3 cells/uL increase of neutrophils'; control stream $THETA(7)

    # ============================================================
    # Inter-individual variability -- Verrest 2023 Table 3
    # 'Between-subject variability'. Only CL/F carries BSV; the control
    # stream fixes the V2 and ka etas to zero ($OMEGA '0 FIX'), so no
    # eta is encoded for those. The paper reports BSV as CV%, and its
    # CV% is sqrt(omega^2) (not the log-normal sqrt(exp(omega^2)-1)
    # form): 56.1% squares to 0.3147, which reproduces the control
    # stream's $OMEGA of 0.314. omega^2 is therefore CV^2 directly.
    # ============================================================
    etalcl ~ 0.561^2
    # Table 3 BSV 'CL (CV%)': 56.1 (95% CI 44.1-76.6), shrinkage 0.0% -> omega^2 = 0.561^2 = 0.3147; control stream $OMEGA '0.314 ; 1 CL'

    # ============================================================
    # Residual error -- Verrest 2023 Table 3 'Residual variability'.
    # The control stream writes W = sqrt(THETA(8)^2*IPRED^2 +
    # THETA(9)^2) with THETA(9) fixed to zero, i.e. a pure
    # proportional error. The stream's THETA(8) of 0.285 is a stale
    # initial estimate; the final value is Table 3's 53.4% CV.
    # ============================================================
    propSd <- 0.534
    label("Proportional residual error on Cc (fraction)")                    # Table 3 'Proportional error (CV%)': 53.4 (95% CI 50.9-56.5), shrinkage 4.9%
  })

  model({
    # ------------------------------------------------------------
    # Individual PK parameters. Allometry on BASELINE body weight
    # normalized to the cohort median of 27.5 kg (Table 3 legend).
    # ------------------------------------------------------------
    # The paper's neutrophil slope is expressed per 10^3 cells/uL and
    # is centred on a median of 0.98 x 10^3 cells/uL, while the
    # canonical NEUT column is in cells/uL. Rescale first so that a
    # standard clinical neutrophil column can be supplied unchanged.
    neutK <- NEUT / 1000

    # CL_TV = CL_pop * (WT/WTmed)^0.75 * (1 + (NEUTR - NEUTRmed) * COV_CL,neutr)
    # (Verrest 2023 Table 3 equation). With e_neut_cl = -0.13, a
    # neutrophil count of 1.0 gives 2.61 L/h and 2.5 gives 2.10 L/h at
    # the reference weight, reproducing the worked example in Results.
    cl <- exp(lcl + etalcl) * (WT_BASE / 27.5)^e_wt_base_cl_q *
      (1 + e_neut_cl * (neutK - 0.98))
    vc <- exp(lvc) * (WT_BASE / 27.5)^e_wt_base_vc_vp
    q  <- exp(lq)  * (WT_BASE / 27.5)^e_wt_base_cl_q
    vp <- exp(lvp) * (WT_BASE / 27.5)^e_wt_base_vc_vp
    ka <- exp(lka)

    # ------------------------------------------------------------
    # Two-compartment disposition with first-order absorption from the
    # intramuscular injection site (control stream $MODEL / $DES).
    # ------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka  * depot
    d/dt(central)     <-  ka  * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot) <- exp(lfdepot)

    # ------------------------------------------------------------
    # Observation. Amounts are mg of paromomycin base and vc is in L,
    # so central / vc is mg/L = ug/mL, the unit in which the paper
    # reports exposure (AUC0-24 in ug*h/mL). The control stream's
    # S2 = V2/1000 exists only to match its ng/mL DV column.
    # ------------------------------------------------------------
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
