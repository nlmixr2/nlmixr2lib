vanderWalt_2025_avelumab <- function() {
  description <- "Bidirectional joint PK-ADA model for avelumab (anti-PD-L1 IgG1) in adults with solid tumors: two-compartment PK with sigmoid time-varying clearance, linked to a two-state continuous-time Markov model of anti-drug-antibody (ADA) status in which the probability of ADA positivity raises clearance and avelumab concentration lowers the ADA-negative to ADA-positive transition rate (van der Walt 2025)"
  reference <- "van der Walt JS, Wilkins J, Khandelwal A, Venkatakrishnan K, Gao W, Milenkovic-Grisic AM. Interplay between pharmacokinetics and immunogenicity of therapeutic proteins: stepwise development of a bidirectional joint pharmacokinetics-anti-drug antibodies model. J Pharmacokinet Pharmacodyn. 2025;52(3):33. doi:10.1007/s10928-025-09971-w"
  vignette <- "vanderWalt_2025_avelumab"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    central = list(analyte = "avelumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "avelumab", units = "mg", specimen = "serum", verified = TRUE),
    s_adaneg = list(
      analyte = "probability of being in the ADA-negative state",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    ),
    s_adapos = list(
      analyte = "probability of being in the ADA-positive state",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling with fixed exponents 0.75 on CL and Q and 1 on V1 and V2 (Supplementary Table 7). The paper does not state the reference weight; 70 kg is used, the reference of the standard allometry it cites (Anderson and Holford 2008).",
      source_name = "WT"
    ),
    ADA_POS_BASE = list(
      description = "ADA status at baseline (before the first avelumab dose)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = ADA-negative at baseline",
      notes = "Multiplies the ADA-negative to ADA-positive transition rate constant by (1 + 4.38). Also sets the initial state of the Markov chain: s_adapos(0) = ADA_POS_BASE.",
      source_name = "baseline ADA status"
    ),
    TUMTP_GASTRIC = list(
      description = "Gastric cancer or gastroesophageal junction cancer (GC/GEJC)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = not GC/GEJC; all TUMTP_* = 0 is non-small cell lung cancer (NSCLC)",
      notes = "JAVELIN Gastric 100 (1L maintenance) and JAVELIN Gastric 300 (3L), plus the gastric expansion of JAVELIN Solid Tumor JPN.",
      source_name = "tumor type (GC/GEJC)"
    ),
    TUMTP_BLADDER = list(
      description = "Urothelial carcinoma (UC)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = not UC; all TUMTP_* = 0 is NSCLC",
      notes = "UC cohorts of JAVELIN Solid Tumor (2L UC); the paper does not split bladder from upper-tract primaries.",
      source_name = "tumor type (UC)"
    ),
    TUMTP_MCC = list(
      description = "Merkel cell carcinoma (MCC)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = not MCC; all TUMTP_* = 0 is NSCLC",
      notes = "JAVELIN Merkel 200 (2L MCC).",
      source_name = "tumor type (MCC)"
    ),
    TUMTP_OTHER = list(
      description = "Other solid tumors",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = one of NSCLC, GC/GEJC, UC, MCC",
      notes = "Pool of the remaining solid tumors (lung, gastrointestinal tract, skin, breast, head and neck), mainly from JAVELIN Solid Tumor JPN; there were no ADA-positive patients of this type in the model development data set.",
      source_name = "tumor type (other solid tumors)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1513L,
    n_studies = 7L,
    age_range = "adults (not tabulated)",
    weight_range = "not reported",
    disease_state = "Advanced or metastatic solid tumors: NSCLC (1L and 2L), gastric / gastroesophageal junction cancer (1L maintenance and 3L), urothelial carcinoma (2L), Merkel cell carcinoma (2L) and other solid tumors",
    dose_range = "3, 10 or 20 mg/kg IV every 2 weeks; 10 mg/kg weekly for 12 weeks then every 2 weeks (JAVELIN Lung 100)",
    regions = "Global (multinational phase I-III trials, including a Japanese phase I study)",
    ada_incidence = "5.7-15.3% of patients ADA-positive at least once, by study (Table 1)",
    notes = "Model development subset (1513 of 1892 patients, about 80%, stratified by ever-ADA-positive status); 379 patients were held out for external validation. Studies: JAVELIN Solid Tumor (UC cohorts), Solid Tumor JPN, Merkel 200, Lung 200, Lung 100, Gastric 100, Gastric 300 (Table 1)."
  )

  ini({
    # Avelumab PK (Table 5, column 'ADA on CL + C1 on ADA'). Only CL was
    # re-estimated in the joint model; the other PK parameters were fixed to
    # the base population PK model.
    lcl <- log(0.0275); label("Clearance at time zero for a 70 kg patient with no ADA (L/h)") # Table 5 'ADA on CL + C1 on ADA': CL = 0.0275 L/h
    lvc <- fixed(log(3.52)); label("Central volume of distribution V1 for a 70 kg patient (L)") # Table 5: V1 = 3.52 L, fixed
    lvp <- fixed(log(0.582)); label("Peripheral volume of distribution V2 for a 70 kg patient (L)") # Table 5: peripheral volume = 0.582 L, fixed
    lq <- fixed(log(0.0128)); label("Intercompartmental clearance Q for a 70 kg patient (L/h)") # Table 5: intercompartmental clearance = 0.0128 L/h, fixed

    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL (unitless)") # Supplementary Table 7: weight on clearance = 0.75, fixed
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on V1 (unitless)") # Supplementary Table 7: weight on central volume = 1, fixed
    e_wt_vp <- fixed(1); label("Allometric exponent of body weight on V2 (unitless)") # Supplementary Table 7: weight on peripheral volume = 1, fixed
    e_wt_q <- fixed(0.75); label("Allometric exponent of body weight on Q (unitless)") # Supplementary Table 7: weight on intercompartmental clearance = 0.75, fixed

    # Sigmoid time-varying clearance, CL(t) = CL * exp(Imax * t^g / (T50^g + t^g))
    # (Methods equation), with additive IIV on Imax.
    cl_time_max <- fixed(-0.015); label("Maximal change in log-clearance as time >> T50, Imax (unitless)") # Table 5: Imax = -0.015, fixed
    lcl_t50 <- fixed(log(51.9)); label("Time to half of the maximal clearance change, T50 (day)") # Table 5: T50 = 51.9 days, fixed
    lcl_time_hill <- fixed(log(2.59)); label("Shape of the time effect on clearance, gamma (unitless)") # Table 5: gamma = 2.59, fixed

    # Two-state continuous-time Markov model of ADA status (Table 5)
    llambda01 <- log(0.00406); label("ADA-negative to ADA-positive transition rate constant, NSCLC, ADA-negative at baseline (1/h)") # Table 5 'ADA on CL + C1 on ADA': lambda01 = 0.00406 /h
    llambda10 <- fixed(log(0.146)); label("ADA-positive to ADA-negative transition rate constant (1/h)") # Table 5: lambda10 = 0.146 /h, fixed
    e_ada_pos_base_lambda01 <- fixed(4.38); label("Fractional change in lambda01 when ADA-positive at baseline (unitless)") # Table 5: 'Baseline ADA+: change in K34' = 4.38, fixed
    e_tumtp_gastric_lambda01 <- fixed(1.17); label("Fractional change in lambda01 for GC/GEJC vs NSCLC (unitless)") # Table 5: change in lambda01 relative to NSCLC, GC/GEJC = 1.17, fixed
    e_tumtp_bladder_lambda01 <- fixed(-0.404); label("Fractional change in lambda01 for UC vs NSCLC (unitless)") # Table 5: change in lambda01 relative to NSCLC, UC = -0.404, fixed
    e_tumtp_mcc_lambda01 <- fixed(-0.373); label("Fractional change in lambda01 for MCC vs NSCLC (unitless)") # Table 5: change in lambda01 relative to NSCLC, MCC = -0.373, fixed
    e_tumtp_other_lambda01 <- fixed(-0.99); label("Fractional change in lambda01 for other solid tumors vs NSCLC (unitless)") # Table 5: change in lambda01 relative to NSCLC, other solid tumors = -0.99, fixed

    # Link 1: probability of ADA+ on clearance (linear, P1eff = 1 + sl * p1)
    e_p_adapos_cl <- 0.15; label("Slope of the fractional clearance increase on the probability of ADA positivity (unitless)") # Table 5 'ADA on CL + C1 on ADA': slope of ADA effect on CL = 0.15

    # Link 2: avelumab concentration on lambda01 (sigmoid Emax)
    emax_lambda01 <- -0.374; label("Maximal fractional change in lambda01 at high avelumab concentration (unitless)") # Table 5 'ADA on CL + C1 on ADA': effect of exposure on ADA, maximum change = -0.374
    lec50_lambda01 <- log(349); label("Avelumab concentration giving half the maximal lambda01 change (ug/mL)") # Table 5: 50% of maximum change = 349; Results text: 349 ug/mL
    lhill_lambda01 <- log(2.58); label("Shape of the concentration effect on lambda01 (unitless)") # Table 5: shape of effect curve = 2.58

    # IIV (Table 5, column 'ADA on CL + C1 on ADA'). Percentages are the
    # log-normal CV, sqrt(exp(omega^2) - 1); omega^2 = log(1 + CV^2).
    etalcl ~ log(1 + 0.399^2) # Table 5: IIV on CL = 39.9%
    etalvc ~ log(1 + 0.246^2) # Table 5: IIV on V1 = 24.6%
    etallambda01 ~ log(1 + 3.29^2) # Table 5: IIV on lambda01 = 329%
    etacl_time_max ~ fixed(0.218^2) # Table 5: 'IIV on TMAX' = 21.8% (held at the base-model value); additive on Imax, SD 0.218

    propSd <- fixed(0.211); label("Proportional residual error (fraction)") # Table 5: proportional residual error PK = 21.1%, fixed
    addSd <- fixed(2.3786); label("Additive residual error (ug/mL)") # Table 5: additive residual error PK = 2378.6 ng/mL, fixed
  })
  model({
    # Allometric scaling (reference 70 kg)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vp <- exp(lvp) * (WT / 70)^e_wt_vp
    q <- exp(lq) * (WT / 70)^e_wt_q

    # Time-varying clearance; T50 is in days and model time in hours
    tday <- t / 24
    cl_time_max_i <- cl_time_max + etacl_time_max
    cl_t50 <- exp(lcl_t50)
    cl_time_hill <- exp(lcl_time_hill)
    td_cl <- exp(cl_time_max_i * tday^cl_time_hill / (cl_t50^cl_time_hill + tday^cl_time_hill))

    # Probability of ADA positivity and its effect on clearance
    p_adapos <- s_adapos / (s_adaneg + s_adapos)
    adapos_cl <- 1 + e_p_adapos_cl * p_adapos

    kel <- cl * td_cl * adapos_cl / vc
    k12 <- q / vc
    k21 <- q / vp

    Cc <- central / vc

    # Concentration effect on the ADA-negative to ADA-positive transition
    cpos <- Cc
    if (cpos < 0) cpos <- 0
    ec50_lambda01 <- exp(lec50_lambda01)
    hill_lambda01 <- exp(lhill_lambda01)
    conc_lambda01 <- 1 + emax_lambda01 * cpos^hill_lambda01 / (ec50_lambda01^hill_lambda01 + cpos^hill_lambda01)

    lambda01 <- exp(llambda01 + etallambda01) *
      (1 + e_ada_pos_base_lambda01 * ADA_POS_BASE) *
      (1 + e_tumtp_gastric_lambda01 * TUMTP_GASTRIC) *
      (1 + e_tumtp_bladder_lambda01 * TUMTP_BLADDER) *
      (1 + e_tumtp_mcc_lambda01 * TUMTP_MCC) *
      (1 + e_tumtp_other_lambda01 * TUMTP_OTHER) *
      conc_lambda01
    lambda10 <- exp(llambda10)

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(s_adaneg) <- -lambda01 * s_adaneg + lambda10 * s_adapos
    d/dt(s_adapos) <- lambda01 * s_adaneg - lambda10 * s_adapos
    s_adaneg(0) <- 1 - ADA_POS_BASE
    s_adapos(0) <- ADA_POS_BASE

    Cc ~ prop(propSd) + add(addSd)

    # Binary ADA status: category 1 = ADA-negative, 2 = ADA-positive
    p0_ada <- 1 - p_adapos
    ada_status ~ c(p0_ada)
  })
}
