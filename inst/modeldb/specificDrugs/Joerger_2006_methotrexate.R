Joerger_2006_methotrexate <- function() {
  description <- "Population PK model for methotrexate (MTX) and its principal circulating metabolite 7-hydroxy-methotrexate (7-OH-MTX) in adult cancer patients receiving high-dose intravenous MTX therapy (Joerger 2006). Joint parent + metabolite model: linear 3-compartment MTX (central + two peripheral compartments) with first-order elimination from the central compartment, feeding a linear 2-compartment 7-OH-MTX disposition through a fixed metabolic fraction of 10 percent. Additive-linear covariate effects of baseline creatinine clearance (Cockcroft-Gault, raw mL/min, truncated at 140), concurrent benzimidazole-class proton-pump-inhibitor comedication, and prior NSAID administration on both MTX and 7-OH-MTX total clearance."
  reference   <- "Joerger M, Huitema ADR, van den Bongard HJGD, Baas P, Schornagel JH, Schellens JHM, Beijnen JH. Determinants of the elimination of methotrexate and 7-hydroxy-methotrexate following high-dose infusional therapy to cancer patients. Br J Clin Pharmacol. 2006;62(1):71-80. doi:10.1111/j.1365-2125.2005.02513.x"
  vignette    <- "Joerger_2006_methotrexate"
  units       <- list(time = "hour", dosing = "umol", concentration = "umol/L")

  covariateData <- list(
    CRCL = list(
      description        = "Baseline creatinine clearance estimated from serum creatinine via the Cockcroft-Gault formula and truncated at 140 mL/min. Time-fixed per cycle; the baseline value before each high-dose MTX cycle is used.",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Raw Cockcroft-Gault (not BSA-normalized), to match Joerger 2006 Methods 'Population pharmacokinetic analysis' ('CL_CREA according to the Cockcroft-Gault formula, and with values > 140 mL/min truncated at this value'). Median 87.5 mL/min, range 40-140 mL/min in the source cohort (Table 1). Centered at 87 mL/min in the covariate equations (Joerger 2006 Eqs 1 and 2). Cockroft-Gault precedent for raw mL/min CRCL is established in Delattre_2010_amikacin.R.",
      source_name        = "CL_CREA"
    ),
    CONMED_PPI = list(
      description        = "Concurrent benzimidazole-class proton-pump-inhibitor (omeprazole 20-40 mg daily or lansoprazole 30 mg daily) at the time of high-dose MTX administration; 1 = on PPI, 0 = no PPI.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant benzimidazole-class PPI use)",
      notes              = "13 of 76 patients (17 percent) were taking benzimidazole-class PPIs concurrently with HDMTX (Joerger 2006 Results page 75: 10 on omeprazole 20-40 mg daily, 3 on lansoprazole 30 mg daily). No difference was demonstrable between omeprazole and lansoprazole effects, nor between dose strata (Results page 77 paragraph after Equations 1 and 2), so the indicator is a class-level binary. Joerger 2006 names the column PPI in Equations 1 and 2.",
      source_name        = "PPI"
    ),
    CONMED_NSAID = list(
      description        = "Prior NSAID exposure (diclofenac 75-300 mg daily in 5 patients, ibuprofen 100 mg three times daily in 1 patient) within the days leading up to high-dose MTX administration; 1 = prior NSAID, 0 = no prior NSAID. NSAIDs were discontinued on the day of HDMTX in all patients.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no prior NSAID use)",
      notes              = "6 of 76 patients (8 percent) had prior NSAID exposure (Joerger 2006 Results page 75). Joerger 2006 names the column NSAID in Equations 1 and 2. Despite being stopped on the day of HDMTX, the prior NSAID exposure still significantly impaired MTX and 7-OH-MTX clearance.",
      source_name        = "NSAID"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight (kg).",
      units       = "kg",
      type        = "continuous",
      notes       = "Tested as a candidate covariate on CL and Vd of MTX and 7-OH-MTX (Joerger 2006 Methods 'Population pharmacokinetic analysis'); not retained in the final model. Cohort median 78 kg (BSA median 1.94 m^2, range 1.56-2.45)."
    ),
    AGE = list(
      description = "Patient age (years).",
      units       = "years",
      type        = "continuous",
      notes       = "Tested as a candidate covariate; not retained. Cohort median 51.1 years, range 17.1-77.0 (Joerger 2006 Table 1)."
    ),
    SEXF = list(
      description = "Patient sex (1 = female, 0 = male).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested as a candidate covariate; not retained. Cohort 62 males / 14 females (Joerger 2006 Table 1). Joerger 2006 reports the source column as GEN with GEN = 0 for females and GEN = 1 for males; SEXF inverts the encoding to match the nlmixr2lib canonical (SEXF = 1 for female)."
    ),
    BSA = list(
      description = "Body surface area (m^2).",
      units       = "m^2",
      type        = "continuous",
      notes       = "Tested as a candidate covariate on Vd of MTX and 7-OH-MTX (Joerger 2006 Methods); not retained. Cohort median 1.94 m^2, range 1.56-2.45."
    ),
    AST = list(
      description = "Aspartate aminotransferase (U/L), surrogate marker for hepatic function.",
      units       = "U/L",
      type        = "continuous",
      notes       = "Tested as a candidate covariate; not retained. Cohort median 17 U/L (Joerger 2006 Table 1)."
    ),
    ALT = list(
      description = "Alanine aminotransferase (U/L), surrogate marker for hepatic function.",
      units       = "U/L",
      type        = "continuous",
      notes       = "Tested as a candidate covariate; not retained. Cohort median 20 U/L (Joerger 2006 Table 1)."
    ),
    BILIRUBIN = list(
      description = "Total bilirubin (umol/L), surrogate marker for hepatic function.",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Tested as a candidate covariate; not retained. Cohort median 6.0 umol/L (Joerger 2006 Table 1)."
    ),
    ALKPHOS = list(
      description = "Alkaline phosphatase (U/L), surrogate marker for hepatic function.",
      units       = "U/L",
      type        = "continuous",
      notes       = "Tested as a candidate covariate; not retained. Cohort median 92 U/L (Joerger 2006 Table 1)."
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 76,
    n_studies       = 1,
    n_cycles        = 304,
    age_range       = "17.1-77.0 years",
    age_median      = "51.1 years",
    sex_distribution = "62 males / 14 females",
    bsa_range       = "1.56-2.45 m^2",
    bsa_median      = "1.94 m^2",
    crcl_range      = "40-140 mL/min (raw Cockcroft-Gault, truncated at 140)",
    crcl_median     = "87.5 mL/min",
    serum_creatinine_range = "32-433 umol/L",
    disease_state   = "Solid tumours treated at The Netherlands Cancer Institute. Diagnoses: malignant pleural mesothelioma (n = 29), gastro-oesophageal cancer (n = 20), non-Hodgkin lymphoma (n = 12), head and neck cancer (n = 10), choriocarcinoma (n = 2), one each of acute lymphocytic leukaemia, trophoblastic tumour, and osteosarcoma.",
    dose_range      = "Intravenous methotrexate 300 mg/m^2 to 12 g/m^2 over 1-24 h infusion. 61 of 76 patients (80 percent) received MTX 1000-5000 mg/m^2 over a 1-6 h infusion. Routine HDMTX-supportive care (aggressive hydration, urine alkalinization to pH = 7, oral leucovorin rescue 15 mg every 6 h starting 24 h after MTX, escalated leucovorin if plasma MTX > 0.1 umol/L at 48 h).",
    co_medication   = "Concurrent benzimidazole-class proton-pump inhibitors (omeprazole or lansoprazole, n = 13) and prior NSAIDs (diclofenac or ibuprofen, n = 6) were retained as covariates on CL. Concurrent corticosteroids (all n = 12 NHL patients), 5-fluorouracil, vinca alkaloids, L-asparaginase, and intravenous doxorubicin (40 mg every 2 weeks; mesothelioma study) were tested but not retained.",
    regions         = "The Netherlands (single-centre, The Netherlands Cancer Institute / Antoni van Leeuwenhoek Hospital / Slotervaart Hospital, Amsterdam).",
    n_observations  = "21 of 76 patients had intensive sampling (end of infusion plus 3, 6, 8 h; 34 sampled cycles); remaining 55 had routine 24 h and 48 h samples only (270 cycles). HPLC LLOQ 0.04 umol/L for both MTX and 7-OH-MTX; within- and between-day CV <= 7 percent.",
    notes           = "Cohort baseline demographics in Joerger 2006 Table 1 (page 74). Two outliers (ID-36 and ID-63) with markedly decreased CL_MTX were retained per the source; ID-63 (29-year-old male, mesothelioma, 3 g 3-h MTX infusion plus doxorubicin) developed CTC grade 3 anuric renal failure on day 2, requiring continuous veno-venous haemofiltration, and died 20 days post-HDMTX (treatment-related mortality 1.3 percent)."
  )

  ini({
    # Structural PK - methotrexate (parent), 3-compartment linear disposition.
    # All point estimates from Joerger 2006 Table 2 'Full data set Estimate'
    # column (p. 76).
    lcl  <- log(8.85);  label("Typical MTX clearance at median CRCL (87 mL/min), no PPI, no NSAID (L/h)") # Joerger 2006 Table 2: CL_MTX = 8.85 L/h; Eqn 1 baseline at median CRCL with PPI = 0, NSAID = 0
    lvc  <- log(23.0);  label("MTX central volume of distribution (L)")                                   # Joerger 2006 Table 2: Volume of CENTRAL_MTX = 23.0 L
    lvp  <- log(185);   label("MTX first peripheral volume of distribution (L)")                          # Joerger 2006 Table 2: Volume of PERIPHERAL-1_MTX = 185 L
    lvp2 <- log(5.34);  label("MTX second peripheral volume of distribution (L)")                         # Joerger 2006 Table 2: Volume of PERIPHERAL-2_MTX = 5.34 L
    lq   <- log(0.444); label("MTX intercompartmental clearance to first peripheral compartment (L/h)")   # Joerger 2006 Table 2: Q1 = 0.444 L/h
    lq2  <- log(0.716); label("MTX intercompartmental clearance to second peripheral compartment (L/h)")  # Joerger 2006 Table 2: Q2 = 0.716 L/h

    # Structural PK - 7-OH-MTX (metabolite), 2-compartment linear disposition.
    lcl_7ohmtx <- log(2);     label("Typical 7-OH-MTX clearance at median CRCL (87 mL/min), no PPI, no NSAID (L/h)") # Joerger 2006 Table 2: CL_7-OH-MTX = 2 L/h; Eqn 2 baseline at median CRCL with PPI = 0, NSAID = 0
    lvc_7ohmtx <- log(21.6);  label("7-OH-MTX central volume of distribution (L)")                                    # Joerger 2006 Table 2: Volume of CENTRAL_7-OH-MTX = 21.6 L
    lvp_7ohmtx <- log(27.7);  label("7-OH-MTX peripheral volume of distribution (L)")                                 # Joerger 2006 Table 2: Volume of PERIPHERAL_7-OH-MTX = 27.7 L
    lq_7ohmtx  <- log(0.429); label("7-OH-MTX intercompartmental clearance (L/h)")                                    # Joerger 2006 Table 2: Q3 = 0.429 L/h

    # Fraction of MTX clearance converted to 7-OH-MTX. Fixed at 10 percent
    # per Joerger 2006 Results page 75: 'we assumed that 10 percent of MTX
    # was metabolized to 7-OH-MTX, in accordance with literature data.
    # Fixing the metabolic fraction to higher (up to 50 percent) or lower
    # (down to 2 percent) values resulted in a decreased fit with an
    # increased OFV.'
    lfm <- fixed(log(0.10)); label("Log fraction of MTX clearance converted to 7-OH-MTX (unitless)")     # Joerger 2006 Results p. 75: assumed 10 percent metabolic fraction; FIXED in NONMEM

    # Additive-linear covariate effects on parent CL_MTX (linear L/h per
    # unit covariate). Encoded with the standard centered form
    # (CRCL - 87) so the baseline lcl is the typical CL at the cohort
    # median CRCL of 87 mL/min. See Errata in the vignette for the sign
    # of the CRCL term: the published Equation 1 prints
    # '+ 0.0423 x (87 - CL_CREA)' which contradicts the abstract's
    # positive correlation between CL_CREA and CL_MTX; the operationally
    # consistent form '0.0423 x (CRCL - 87)' is used here so that higher
    # creatinine clearance produces higher MTX clearance.
    e_crcl_cl  <-  0.0423; label("Linear additive effect of (CRCL - 87) on MTX clearance (L/h per mL/min)") # Joerger 2006 Eq 1: 0.0423 L/h per mL/min CRCL deviation from median 87 mL/min
    e_ppi_cl   <- -2.45;   label("Linear additive effect of CONMED_PPI = 1 on MTX clearance (L/h)")         # Joerger 2006 Eq 1: -2.45 L/h decrement with concurrent benzimidazole PPI
    e_nsaid_cl <- -1.46;   label("Linear additive effect of CONMED_NSAID = 1 on MTX clearance (L/h)")       # Joerger 2006 Eq 1: -1.46 L/h decrement with prior NSAID exposure

    # Additive-linear covariate effects on metabolite CL_7-OH-MTX
    # (linear L/h per unit covariate). Sign on CRCL corrected as for
    # the parent (see vignette Errata).
    e_crcl_cl_7ohmtx  <-  0.0123; label("Linear additive effect of (CRCL - 87) on 7-OH-MTX clearance (L/h per mL/min)") # Joerger 2006 Eq 2: 0.0123 L/h per mL/min CRCL deviation from median 87 mL/min
    e_ppi_cl_7ohmtx   <- -0.369;  label("Linear additive effect of CONMED_PPI = 1 on 7-OH-MTX clearance (L/h)")         # Joerger 2006 Eq 2: -0.369 L/h decrement with concurrent benzimidazole PPI
    e_nsaid_cl_7ohmtx <- -0.357;  label("Linear additive effect of CONMED_NSAID = 1 on 7-OH-MTX clearance (L/h)")       # Joerger 2006 Eq 2: -0.357 L/h decrement with prior NSAID exposure

    # Interindividual variability. Joerger 2006 reports IIV as CV percent
    # in Table 2; converted here to log-normal variance via
    # omega^2 = log(1 + CV^2). IIV was reported as NOT included on
    # V_CENTRAL_MTX, V_PERIPHERAL-1_MTX, and Q2 (final paragraph,
    # Results page 75: 'did not improve the fit').
    # CL_MTX:                  log(1 + 0.196^2) = 0.03770
    # V_CENTRAL_7-OH-MTX:      log(1 + 0.0822^2) = 0.006735
    # V_PERIPHERAL-2_MTX:      log(1 + 0.310^2) = 0.09175
    # Q1:                      log(1 + 0.320^2) = 0.09748
    # CL_7-OH-MTX:             log(1 + 0.310^2) = 0.09175
    # V_PERIPHERAL_7-OH-MTX:   log(1 + 0.414^2) = 0.15824
    # Q3:                      log(1 + 0.271^2) = 0.07086
    etalcl        ~ 0.03770   # Joerger 2006 Table 2: 19.6 percent CV on CL_MTX
    etalvp2       ~ 0.09175   # Joerger 2006 Table 2: 31.0 percent CV on V_PERIPHERAL-2_MTX
    etalq         ~ 0.09748   # Joerger 2006 Table 2: 32.0 percent CV on Q1
    etalcl_7ohmtx ~ 0.09175   # Joerger 2006 Table 2: 31.0 percent CV on CL_7-OH-MTX
    etalvc_7ohmtx ~ 0.006735  # Joerger 2006 Table 2: 8.22 percent CV on V_CENTRAL_7-OH-MTX
    etalvp_7ohmtx ~ 0.15824   # Joerger 2006 Table 2: 41.4 percent CV on V_PERIPHERAL_7-OH-MTX
    etalq_7ohmtx  ~ 0.07086   # Joerger 2006 Table 2: 27.1 percent CV on Q3

    # Residual error. Joerger 2006 Methods page 73 specifies an additive
    # error model on the log-scale: log(C_ij) = log(C_ij_pred) + eps_ij.
    # In the linear-space convention used by nlmixr2lib this corresponds
    # to a proportional residual model y = ypred * (1 + eps); for small
    # eps the standard deviation on the log-scale approximates the
    # proportional CV in the linear scale (see vignette Assumptions and
    # deviations for the equivalence note).
    propSd           <- 0.523; label("MTX plasma concentration proportional residual error (fraction)")      # Joerger 2006 Table 2: 52.3 percent residual variability on MTX plasma concentration
    propSd_7ohmtx <- 0.571; label("7-OH-MTX plasma concentration proportional residual error (fraction)") # Joerger 2006 Table 2: 57.1 percent residual variability on 7-OH-MTX plasma concentration
  })

  model({
    # 1. Linear-additive covariate-adjusted typical clearances on the
    #    parent and metabolite. CRCL centred at 87 mL/min so the baseline
    #    log clearance corresponds to the cohort median CRCL with no PPI
    #    and no NSAID.
    tvcl <- exp(lcl) +
            e_crcl_cl  * (CRCL - 87) +
            e_ppi_cl   * CONMED_PPI +
            e_nsaid_cl * CONMED_NSAID
    tvcl_7ohmtx <- exp(lcl_7ohmtx) +
                   e_crcl_cl_7ohmtx  * (CRCL - 87) +
                   e_ppi_cl_7ohmtx   * CONMED_PPI +
                   e_nsaid_cl_7ohmtx * CONMED_NSAID

    # 2. Individual PK parameters. IIV applied multiplicatively around
    #    the linear-covariate-adjusted typical value (Joerger 2006 Eq.
    #    on page 73). Parameters with no reported IIV (V_CENTRAL_MTX,
    #    V_PERIPHERAL-1_MTX, Q2) carry the typical value only.
    cl  <- tvcl * exp(etalcl)
    vc  <- exp(lvc)
    vp  <- exp(lvp)
    vp2 <- exp(lvp2 + etalvp2)
    q   <- exp(lq  + etalq)
    q2  <- exp(lq2)

    cl_7ohmtx <- tvcl_7ohmtx * exp(etalcl_7ohmtx)
    vc_7ohmtx <- exp(lvc_7ohmtx + etalvc_7ohmtx)
    vp_7ohmtx <- exp(lvp_7ohmtx + etalvp_7ohmtx)
    q_7ohmtx  <- exp(lq_7ohmtx  + etalq_7ohmtx)

    # 3. Fixed metabolic fraction.
    fm <- exp(lfm)

    # 4. Micro-constants.
    kel       <- cl / vc                        # total MTX elimination rate from central
    k12       <- q  / vc
    k21       <- q  / vp
    k13       <- q2 / vc
    k31       <- q2 / vp2
    kel_7ohmtx <- cl_7ohmtx / vc_7ohmtx
    k12_7ohmtx <- q_7ohmtx  / vc_7ohmtx
    k21_7ohmtx <- q_7ohmtx  / vp_7ohmtx

    # 5. ODE system. IV MTX enters central; rate is supplied via the
    #    event table (RATE column). Of the total MTX elimination flux
    #    kel * central, a fraction fm appears as 7-OH-MTX in its
    #    central compartment.
    d/dt(central)         <- -kel * central -
                              k12 * central + k21 * peripheral1 -
                              k13 * central + k31 * peripheral2
    d/dt(peripheral1)     <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2)     <-  k13 * central - k31 * peripheral2
    d/dt(central_7ohmtx)  <-  fm * kel * central -
                              kel_7ohmtx * central_7ohmtx -
                              k12_7ohmtx * central_7ohmtx +
                              k21_7ohmtx * peripheral1_7ohmtx
    d/dt(peripheral1_7ohmtx) <- k12_7ohmtx * central_7ohmtx -
                                k21_7ohmtx * peripheral1_7ohmtx

    # 6. Observations (umol/L). Mole-equivalent stoichiometry is assumed
    #    (parent and metabolite molecular weights are within 4 percent:
    #    MTX 454.4 g/mol, 7-OH-MTX 470.4 g/mol).
    Cc        <- central        / vc
    Cc_7ohmtx <- central_7ohmtx / vc_7ohmtx

    Cc        ~ prop(propSd)
    Cc_7ohmtx ~ prop(propSd_7ohmtx)
  })
}
