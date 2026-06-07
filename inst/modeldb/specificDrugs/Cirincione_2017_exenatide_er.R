Cirincione_2017_exenatide_er <- function() {
  description <- "Population PK model for extended-release (ER) microsphere SC exenatide in patients with type 2 diabetes (Cirincione 2017 AAPS J): two-compartment disposition with three parallel SC-absorption processes (initial first-order release plus two Savic 2007 analytical transit-compartment chains for the second- and third-phase microsphere release) and parallel linear plus saturable Michaelis-Menten elimination. Disposition parameters (CL, Q, Vc, Vp, Vmax, Km) and the eGFR-on-CL and WT-on-Vc covariate effects are fixed from the IR companion model (Cirincione 2017 BJCP)."
  reference <- paste(
    "Cirincione B, Edwards J, Mager DE.",
    "Population Pharmacokinetics of an Extended-Release Formulation of Exenatide Following",
    "Single- and Multiple-Dose Administration.",
    "AAPS J. 2017;19(2):487-496. doi:10.1208/s12248-016-9975-1.",
    "Disposition parameters (CL_int, Q, Vc_int, Vp, Vmax, Km) and the eGFR-on-CL and",
    "WT-on-Vc effects fixed from the IR companion model:",
    "Cirincione B, Mager DE. Population pharmacokinetics of exenatide.",
    "Br J Clin Pharmacol. 2017;83(3):517-526. doi:10.1111/bcp.13135;",
    "see modellib('Cirincione_2017_exenatide').",
    sep = " "
  )
  vignette <- "Cirincione_2017_exenatide_er"
  depends <- c("Cirincione_2017_exenatide")
  paper_specific_etas <- c("etalktr2", "etalf1", "etalfrelSD")
  paper_specific_residual_sds <- c("expSdSD", "expSdMD")

  units <- list(time = "day", dosing = "ug", concentration = "pg/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on central volume; reference weight 84.8 kg inherited from the IR Cirincione 2017 BJCP model. Table I baseline-only.",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Modification of Diet in Renal Disease (MDRD) estimated glomerular filtration rate (creatinine-based, BSA-normalized)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on linear clearance; reference 80 mL/min/1.73 m^2 inherited from the IR Cirincione 2017 BJCP model. Source column 'eGFR' (MDRD eGFR) maps to the canonical general-scope CRCL covariate.",
      source_name        = "eGFR"
    ),
    STUDY_MD = list(
      description        = "Multi-dose study cohort indicator (1 = phase II multi-dose study with weekly ER exenatide, 0 = phase II single-dose study)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (single-dose phase II study)",
      notes              = "Switches between study-specific relative bioavailability (f_rel(SD) vs f_rel(MD)) and between study-specific log-scale residual SDs (RV SD vs RV MD); see Table II. For the phase III external validation cohort (weekly 2 mg ER for 24 weeks), use STUDY_MD = 1.",
      source_name        = "derived"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 64L,
    n_studies      = 2L,
    age_range      = "30-72 years",
    age_mean       = "53.5 years (SD 9.78)",
    weight_range   = "59-155 kg",
    weight_mean    = "95.7 kg (SD 21)",
    sex_female_pct = 42,
    race_ethnicity = c(White = 37.5, Hispanic = 37.5, Black = 10.9, Asian = 7.8, Other = 6.3),
    renal_function = c(normal = 23, mild = 38, moderate = 3),
    disease_state  = "Type 2 diabetes mellitus.",
    dose_range     = "Single 2.5, 5, 7, or 10 mg SC (single-dose phase II study, n=41) or 0.8 mg or 2 mg SC once weekly for 15 weeks (multi-dose phase II study, n=23).",
    administration = "Subcutaneous, ER (PLGA microsphere) formulation",
    regions        = "Not reported in detail",
    notes          = "Demographics from Table I (combined single-dose + multi-dose cohorts). Mean eGFR 89.6 mL/min/1.73 m^2 (SD 23.2; range 56-169). External validation in a phase III multi-dose cohort (n=119, 2 mg weekly x 24 weeks) is described in the paper Results."
  )

  ini({
    # Estimated absorption parameters (Combined single- and multiple-dose model, Table II right block)
    lka     <- log(3.85);    label("First-order absorption rate constant (1/day)")                              # Table II combined: ka = 3.85 /day (RSE 11.2%)
    lktr1   <- log(0.105);   label("Transit-chain 1 first-order rate constant (1/day)")                         # Table II combined: ktr1 = 0.105 /day (RSE 9.50%)
    lntr1   <- log(0.570);   label("Transit-chain 1 number of transit compartments (dimensionless)")             # Table II combined: N1 = 0.570 (RSE 17.2%)
    lktr2   <- log(0.591);   label("Transit-chain 2 first-order rate constant (1/day)")                         # Table II combined: ktr2 = 0.591 /day (RSE 17.9%)
    lntr2   <- log(26.1);    label("Transit-chain 2 number of transit compartments (dimensionless)")             # Table II combined: N2 = 26.1 (RSE 19.6%)
    lf1     <- log(0.0118);  label("Fraction of bioavailable dose in initial first-order release (fraction)")    # Table II combined: f1 = 1.18% (RSE 8.56%)
    lf2     <- log(0.480);   label("Fraction of bioavailable dose absorbed via transit-chain 1 (fraction)")      # Table II combined: f2 = 48.0% (RSE 8.77%); f3 = 1 - f1 - f2 = 0.508 is derived per Methods (f3 = 1 - (f1 + f2))
    lfrelSD <- log(0.0886);  label("Relative bioavailability vs IR formulation, single-dose study (fraction)")   # Table II combined: f_rel(single-dose study) = 8.86% (RSE 7.56%)
    lfrelMD <- log(0.155);   label("Relative bioavailability vs IR formulation, multi-dose study (fraction)")    # Table II combined: f_rel(MD study) = 15.5% (RSE 5.93%)

    # Disposition parameters (fixed in the ER analysis from the IR Cirincione 2017 BJCP companion model;
    # see modellib('Cirincione_2017_exenatide')). ER values in day units = IR values in hour units * 24.
    lcl       <- fixed(log(110));    label("Linear clearance at reference CRCL (L/day)")                                          # Table II / IR Table 2: CL_int = 110 L/day (= 4.58 L/hr * 24)
    lq        <- fixed(log(89.3));   label("Intercompartmental clearance (L/day)")                                                # Table II / IR Table 2: CLd = 89.3 L/day (= 3.72 L/hr * 24)
    lvc       <- fixed(log(7.03));   label("Central volume at reference body weight (L)")                                          # Table II / IR Table 2: Vc_int = 7.03 L
    lvp       <- fixed(log(7.04));   label("Peripheral volume (L)")                                                                # Table II / IR Table 2: Vp = 7.04 L
    lvmax     <- fixed(log(37.2));   label("Maximum Michaelis-Menten elimination rate (ug/day)")                                   # Table II: Vmax = 0.0372 mg/day = 37.2 ug/day (matches IR 1.55 ug/hr * 24)
    lkm       <- fixed(log(0.567));  label("Michaelis-Menten constant for saturable elimination (ng/mL)")                          # Table II: Km = 567 pg/mL expressed as 0.567 ng/mL for unit-consistency with Cc_ngmL = central/vc (matches IR Km)
    e_crcl_cl <- fixed(0.838);       label("Power effect of CRCL (MDRD eGFR) on linear clearance (unitless)")                      # Table II / IR Table 2: CL_eGFR = 0.838
    e_wt_vc   <- fixed(2.67);        label("Power effect of body weight on central volume (unitless)")                             # Table II / IR Table 2: Vc_wtkg = 2.67

    # Inter-individual variability (combined model). Paper reports CV% on log-normal IIV;
    # nlmixr2 variance = log(1 + CV^2).
    etalka     ~ 0.2551  # ka IIV 53.9% (RSE 37.6%); log(1 + 0.539^2) = 0.2551 -- Table II combined right
    etalktr2   ~ 0.0292  # ktr2 IIV 17.2% (RSE 20.9%); log(1 + 0.172^2) = 0.0292 -- Table II combined right
    etalf1     ~ 0.3677  # f1 IIV 66.7% (RSE 23.4%); log(1 + 0.667^2) = 0.3677 -- Table II combined right
    etalfrelSD ~ 0.1755  # f_rel(single-dose study) IIV 43.8% (RSE 20.1%); log(1 + 0.438^2) = 0.1755 -- Table II combined right

    # Residual error -- log-normal on the observation, study-specific magnitude (Table II combined right).
    expSdSD <- 0.684; label("Log-scale residual SD, single-dose study (Cirincione 2017 RV SD Study)")  # Table II combined right: RV SD Study (Log SD) = 0.684 (RSE 4.58%)
    expSdMD <- 0.376; label("Log-scale residual SD, multi-dose study (Cirincione 2017 RV MD Study)")   # Table II combined right: RV MD Study (Log SD) = 0.376 (RSE 4.57%)
  })

  model({
    # Per-subject absorption parameters
    ka     <- exp(lka   + etalka)
    ktr1   <- exp(lktr1)
    ntr1   <- exp(lntr1)
    ktr2   <- exp(lktr2 + etalktr2)
    ntr2   <- exp(lntr2)
    f1     <- exp(lf1   + etalf1)
    f2     <- exp(lf2)
    # Third absorption fraction is derived structurally (Methods p. 488: f3 = 1 - (f1 + f2)),
    # so the total ER bioavailability across the three parallel processes equals f_rel.
    f3     <- 1 - f1 - f2

    # Study-specific relative bioavailability (single-dose study vs multi-dose study)
    frelSD <- exp(lfrelSD + etalfrelSD)
    frelMD <- exp(lfrelMD)
    frel   <- frelMD * STUDY_MD + frelSD * (1 - STUDY_MD)

    # Disposition parameters (fixed-from-IR with the published covariate effects on CL and Vc)
    cl   <- exp(lcl) * (CRCL / 80)^e_crcl_cl
    vc   <- exp(lvc) * (WT / 84.8)^e_wt_vc
    q    <- exp(lq)
    vp   <- exp(lvp)
    vmax <- exp(lvmax)
    km   <- exp(lkm)

    # Micro-constants. kel combines linear clearance and the saturable Michaelis-Menten arm;
    # the saturable arm is written with central (ug) and km*vc (ug/L * L = ug) so the
    # dimensions cancel and km can be expressed in concentration units (ng/mL) consistent
    # with Cc_ngmL = central / vc.
    k12 <- q / vc
    k21 <- q / vp
    kel <- cl / vc + vmax / (km * vc + central)

    # Three parallel SC absorption processes from a single microsphere dose (Cirincione 2017
    # Fig. 1 / Methods Eq. for the transit-compartment input rate). The Savic 2007 analytical
    # gamma-density input rate is evaluated explicitly from tad(depot) and podo(depot) so each
    # transit-chain contribution is independent of the depot mass:
    #   process 1 -- first-order from depot at rate ka, fraction f1 * frel of the administered
    #                dose loads into depot via f(depot);
    #   process 2 -- analytical transit chain 1 (ntr1, ktr1), fraction f2 * frel;
    #   process 3 -- analytical transit chain 2 (ntr2, ktr2), fraction f3 * frel.
    # NOTE on multi-dose: the analytical transit forms use the most-recent dose only
    # (rxode2 tad/podo are scalars). Process 1 (depot ODE) correctly superposes contributions
    # from earlier doses; processes 2 and 3 do not. Validation against the multi-dose VPC
    # is therefore approximate -- see the vignette Errata section.
    td <- tad(depot)
    pd <- podo(depot)

    in_t1 <- 0
    if (td > 0) {
      in_t1 <- (f2 * frel) * pd *
        exp((ntr1 + 1) * log(ktr1) + ntr1 * log(td) - ktr1 * td - lgamma(ntr1 + 1))
    }

    in_t2 <- 0
    if (td > 0) {
      in_t2 <- (f3 * frel) * pd *
        exp((ntr2 + 1) * log(ktr2) + ntr2 * log(td) - ktr2 * td - lgamma(ntr2 + 1))
    }

    # ODE system: 2-compartment disposition + linear + MM elimination + 3 parallel absorption
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot + in_t1 + in_t2 - kel * central -
                          k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Bioavailability: only f1 * frel fraction of the administered dose enters depot for
    # process 1; processes 2 and 3 deliver their (f2 * frel) and (f3 * frel) fractions
    # directly via in_t1 and in_t2.
    f(depot) <- f1 * frel

    # Observation: central (ug) / vc (L) = ug/L = ng/mL; multiply by 1000 to report in pg/mL
    # matching Table II's Km (567 pg/mL) and the published concentration figures.
    Cc <- (central / vc) * 1000

    # Residual error: log-normal observation, study-specific SD via STUDY_MD switch
    expSd <- expSdMD * STUDY_MD + expSdSD * (1 - STUDY_MD)
    Cc ~ lnorm(expSd)
  })
}
