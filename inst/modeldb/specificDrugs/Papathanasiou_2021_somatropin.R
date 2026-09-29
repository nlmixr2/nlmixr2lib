Papathanasiou_2021_somatropin <- function() {
  description <- "One-compartment population PK model with first-order absorption and zero-order endogenous GH production, linked to an indirect-response IGF-I model with an additive Emax stimulation of IGF-I production, for once-daily subcutaneous somatropin (Norditropin) in children and adults with growth hormone deficiency (Papathanasiou 2021)"
  reference <- "Papathanasiou T, Agerso H, Damholt BB, Rasmussen MH, Kildemoes RJ. Population Pharmacokinetics and Pharmacodynamics of Once-Daily Growth Hormone Norditropin(R) in Children and Adults. Clin Pharmacokinet. 2021;60(9):1217-1226. doi:10.1007/s40262-021-01011-3"
  vignette <- "Papathanasiou_2021_somatropin"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")

  # Issue #482: what each ODE state holds. Verified against Papathanasiou
  # 2021 Figure 1 (absorption compartment 'abs' -> 'central' via Ka, zero-order
  # endogenous input K_Endo into central, CL/F out of central; indirect
  # response from central GH onto the 'IGF-I' compartment). The IGF-I state
  # is carried directly as a serum concentration (Kin in ng/mL/h, Kout in 1/h,
  # Table 3), so its units are ng/mL rather than an amount.
  compartmentData <- list(
    depot = list(analyte = "somatropin", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "somatropin", units = "ug", specimen = "serum", verified = TRUE),
    igf1 = list(analyte = "insulin-like growth factor I", units = "ng/mL", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power covariate on Ka, CL/F and the endogenous baseline GH",
        "concentration with reference weight 70 kg (Methods section 2.6,",
        "first displayed equation, P_i = P_typ * (BW / 70 kg)^theta). On Emax",
        "the reference weight differs by age group: 85 kg for adults and 25 kg",
        "for children (Methods section 2.6, second and third displayed",
        "equations). V/F carries no weight effect. Table 1 does not state",
        "whether weight was time-varying; the trials were at most 4 weeks long,",
        "so it is treated as a baseline value."
      ),
      source_name = "BW"
    ),
    CHILD = list(
      description = "Child (pediatric GHD) age-cohort indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (adult with GHD)",
      notes = paste(
        "The paper codes the age group as 'adult' (1 = adult, 0 = child;",
        "Methods section 2.6, P_i = (P_typ,adult * adult + P_typ,child *",
        "(1 - adult)) * exp(eta)), so CHILD = 1 - adult. The age group switches",
        "(i) the typical Emax and its body-weight exponent and reference weight",
        "and (ii) the proportional residual error of GH (Table 2, separate",
        "AGHD and GHD rows). It also selects the Trial 3 IGF-I residual error,",
        "because Trial 3 enrolled only adults. Children were 6-11 years old",
        "and prepubertal; adults were 23-68 years old (Table 1)."
      ),
      source_name = "adult"
    ),
    STUDY_NCT00936403 = list(
      description = "Trial 2 (NCT00936403; NNC126-0083 pediatric trial) cohort indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Trial 1 NCT01973244 or Trial 3 NCT01706783)",
      notes = paste(
        "Selects the IGF-I proportional residual error estimated for Trial 2",
        "(Table 3, 'Proportional error Trial 2' = 8.7%). The separate",
        "per-trial IGF-I residual errors account for the different IGF-I",
        "assays (Siemens IMMULITE in Trial 2; IDS-iSYS in Trials 1 and 3;",
        "Methods section 2.3). Only meaningful for children (CHILD = 1); it is",
        "ignored when CHILD = 0 because Trial 3 is the only adult trial."
      ),
      source_name = "Trial"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 23,
    n_studies = 3,
    n_observations = "614 GH and 334 IGF-I concentrations",
    age_range = "6-68 years (children 6-11 years; adults 23-68 years)",
    age_mean = "23.4 years (children 8.2 years; adults 58.0 years)",
    weight_range = "17-102.2 kg (children 17-40.5 kg; adults 59.1-102.2 kg)",
    weight_mean = "44.0 kg (Trial 1 26.1 kg; Trial 2 29.7 kg; Trial 3 80.8 kg)",
    sex_female_pct = 17.4,
    race_ethnicity = "Not reported",
    disease_state = "Growth hormone deficiency (prepubertal children and adults), after a 7-14 day wash-out of prior GH treatment",
    dose_range = paste(
      "Once-daily subcutaneous Norditropin: 0.03 mg/kg for 7 days (Trial 1),",
      "0.035 mg/kg for 7 days (Trial 2), and the pre-trial dose (mean",
      "0.0042 mg/kg, 0.2-0.5 mg) for 4 weeks in adults (Trial 3)."
    ),
    regions = "Not reported in the article",
    notes = paste(
      "Pooled Norditropin comparator arms of three phase I trials: Trial 1",
      "(NCT01973244, somapacitan, prepubertal children, n = 8), Trial 2",
      "(NCT00936403, NNC126-0083, prepubertal children, n = 8) and Trial 3",
      "(NCT01706783, somapacitan, adults, n = 7 after exclusion of one adult",
      "with a suspected dosing error). Table 1 gives 8 + 8 + 7 = 23 subjects;",
      "the Results text says '15 children and eight adults', which conflicts",
      "with Table 1's per-trial counts (16 children) -- Table 1 is used here.",
      "Nineteen males and four females."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # PK (Table 2, final PK model). Typical values refer to a 70 kg subject.
    # The PK parameters were held fixed at these estimates while the PK/PD
    # model was estimated (Methods section 2.2); they are not fixed() here
    # because they were estimated in the PK step.
    # ---------------------------------------------------------------------
    lka <- log(0.122); label("Absorption rate constant Ka at 70 kg (1/h)")                          # Table 2 'Ka (1/h)' = 0.122 (95% CI 0.07-0.17), RSE 21.6%
    lvc <- log(28.2); label("Apparent volume of distribution V/F (L)")                              # Table 2 'V/F (L)' = 28.2 (95% CI 17.1-39.3), RSE 20.1%
    lcl <- log(24.7); label("Apparent clearance CL/F at 70 kg (L/h)")                               # Table 2 'CL/F (L/h)' = 24.7 (95% CI 16.2-33.2), RSE 17.6%
    lc0 <- log(0.188); label("Baseline GH concentration from endogenous production at 70 kg (ng/mL)") # Table 2 'GH Base (ng/mL)' = 0.188 (95% CI 0.03-0.35), RSE 43.0%

    e_wt_cl <- 0.982; label("Power exponent of (WT/70) on CL/F (unitless)")                         # Table 2 'theta BW CL/F' = 0.982 (95% CI 0.62-1.34)
    e_wt_ka <- -0.687; label("Power exponent of (WT/70) on Ka (unitless)")                          # Table 2 'theta BW Ka' = -0.687 (95% CI -1.17 to -0.20)
    e_wt_c0 <- -0.991; label("Power exponent of (WT/70) on the baseline GH concentration (unitless)") # Table 2 'theta BW GH Base' = -0.991 (95% CI -1.96 to -0.02)

    # ---------------------------------------------------------------------
    # PD (Table 3, final PK/PD model): indirect response, additive Emax
    # stimulation of the IGF-I production rate.
    # ---------------------------------------------------------------------
    lkin <- log(0.949); label("Zero-order IGF-I production rate at zero GH (ng/mL/h)")               # Table 3 'Kin (ng/mL/h)' = 0.949 (95% CI 0.411-1.49), RSE 28.9%
    lkout <- log(0.0262); label("First-order IGF-I elimination rate constant (1/h)")                 # Table 3 'Kout (1/h)' = 0.0262 (95% CI 0.0218-0.0306), RSE 8.6%
    lec50 <- log(2.13); label("GH concentration giving half-maximal IGF-I production stimulation (ng/mL)") # Table 3 'EC50 (ng/mL)' = 2.13 (95% CI 1.39-2.88), RSE 17.9%

    # Emax for adults was not identifiable and was fixed, with its weight
    # exponent, to the somapacitan meta-analysis estimate (Methods sections
    # 2.4 and 2.6; Table 3 '(fixed)').
    lemax_adult <- fixed(log(15.1)); label("Maximum increase in IGF-I production rate, 85 kg adult (ng/mL/h)") # Table 3 'Emax adult (ng/mL/h)' = 15.1 (fixed)
    e_wt_emax_adult <- fixed(0.46); label("Power exponent of (WT/85) on adult Emax (unitless)")           # Table 3 'theta BW Emax adult' = 0.46 (fixed)
    lemax_ped <- log(6.48); label("Maximum increase in IGF-I production rate, 25 kg child (ng/mL/h)")   # Table 3 'Emax child (ng/mL/h)' = 6.48 (95% CI 4.75-8.70), RSE 14.0%
    e_wt_emax_ped <- 1.94; label("Power exponent of (WT/25) on pediatric Emax (unitless)")             # Table 3 'theta BW Emax child' = 1.94 (95% CI 1.05-2.84), RSE 23.5%

    # ---------------------------------------------------------------------
    # IIV. Tables 2 and 3 report IIV as CV%; converted to log-normal
    # variances with omega^2 = log(CV^2 + 1). The Ka-CL/F and Ka-V/F
    # correlations that were estimated (ESM Table S1, dOFV -12.6) are not
    # reported anywhere in the article or ESM, so the etas are encoded as
    # independent (see vignette). Emax carries ONE random effect shared by
    # the adult and child typical values (Table 3 prints the same 33.1% CV
    # and 13.3% shrinkage on both rows). EC50 has no IIV.
    # ---------------------------------------------------------------------
    etalka ~ 0.5904   # Table 2 Ka IIV CV 89.7%: log(0.897^2 + 1)
    etalvc ~ 0.4704   # Table 2 V/F IIV CV 77.5%: log(0.775^2 + 1)
    etalcl ~ 0.1045   # Table 2 CL/F IIV CV 33.2%: log(0.332^2 + 1)
    etalc0 ~ 1.4190   # Table 2 GH Base IIV CV 177%: log(1.77^2 + 1)
    etalkin ~ 0.8624  # Table 3 Kin IIV CV 117%: log(1.17^2 + 1)
    etalkout ~ 0.02983 # Table 3 Kout IIV CV 17.4%: log(0.174^2 + 1)
    etalemax ~ 0.1040 # Table 3 Emax adult and child IIV CV 33.1%: log(0.331^2 + 1)

    # ---------------------------------------------------------------------
    # Residual error: proportional throughout (Methods section 2.3), with
    # separate magnitudes for adult and child GH data and for each trial's
    # IGF-I data.
    # ---------------------------------------------------------------------
    propSd_Cc_adult <- 0.414; label("Proportional residual SD of GH, adults (fraction)")         # Table 2 'Proportional error AGHD (%)' = 41.4
    propSd_Cc_ped <- 0.561; label("Proportional residual SD of GH, children (fraction)")         # Table 2 'Proportional error GHD (%)' = 56.1
    propSd_IGF1_trial1 <- 0.146; label("Proportional residual SD of IGF-I, Trial 1 (fraction)")  # Table 3 'Proportional error Trial 1 (%)' = 14.6
    propSd_IGF1_trial2 <- 0.087; label("Proportional residual SD of IGF-I, Trial 2 (fraction)")  # Table 3 'Proportional error Trial 2 (%)' = 8.7
    propSd_IGF1_trial3 <- 0.114; label("Proportional residual SD of IGF-I, Trial 3 (fraction)")  # Table 3 'Proportional error Trial 3 (%)' = 11.4
  })

  model({
    # Individual PK parameters (Methods section 2.6: P_i = P_typ * (BW/70)^theta * exp(eta)).
    ka <- exp(lka + etalka) * (WT / 70)^e_wt_ka
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    c0 <- exp(lc0 + etalc0) * (WT / 70)^e_wt_c0
    kel <- cl / vc

    # Zero-order endogenous GH production into central (Figure 1, K_Endo),
    # parameterised by the baseline concentration it sustains: at steady
    # state kendo / (CL/F) = c0. Written in the model's apparent (/F) amount
    # scale; the predicted concentration does not depend on F.
    kendo <- c0 * cl

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot + kendo - kel * central
    central(0) <- c0 * vc

    # ug / L = ng/mL
    Cc <- central / vc

    # Emax switches between the adult and child typical values and weight
    # relationships (Methods section 2.6, second to fourth displayed
    # equations), with one shared random effect.
    emax <- (exp(lemax_adult) * (WT / 85)^e_wt_emax_adult * (1 - CHILD) +
      exp(lemax_ped) * (WT / 25)^e_wt_emax_ped * CHILD) * exp(etalemax)
    kin <- exp(lkin + etalkin)
    kout <- exp(lkout + etalkout)
    ec50 <- exp(lec50)

    # Indirect response with an ADDITIVE Emax effect of total (endogenous +
    # exogenous) GH on the IGF-I production rate (Figure 1; Results section
    # 3.3.1 'best described as additive'; Emax carries production-rate
    # units, ng/mL/h). The IGF-I state starts at its endogenous steady state,
    # (Kin + Emax * c0 / (EC50 + c0)) / Kout, i.e. the Table 3 derived
    # 'Kin,endo' divided by Kout.
    d/dt(igf1) <- kin + emax * Cc / (ec50 + Cc) - kout * igf1
    igf1(0) <- (kin + emax * c0 / (ec50 + c0)) / kout

    IGF1 <- igf1

    propSd_Cc <- propSd_Cc_adult * (1 - CHILD) + propSd_Cc_ped * CHILD
    propSd_IGF1 <- propSd_IGF1_trial3 * (1 - CHILD) +
      CHILD * (propSd_IGF1_trial1 * (1 - STUDY_NCT00936403) +
        propSd_IGF1_trial2 * STUDY_NCT00936403)

    Cc ~ prop(propSd_Cc)
    IGF1 ~ prop(propSd_IGF1)
  })
}
