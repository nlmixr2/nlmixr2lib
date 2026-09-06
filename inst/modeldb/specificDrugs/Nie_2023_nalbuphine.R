Nie_2023_nalbuphine <- function() {
  description <- paste(
    "Two-compartment IV population PK model of nalbuphine in adult patients undergoing general",
    "anaesthesia surgery (Nie 2023; 47 Chinese adults aged 21-78 years and 48-86 kg receiving a",
    "single 15 mg nalbuphine intravenous injection over 2-3 min for induction of anaesthesia,",
    "split into a 27-patient / 353-sample model-building set and a 20-patient / 100-sample",
    "external-validation set). Disposition is parameterized as CL, central volume V1,",
    "intercompartmental clearance Q and peripheral volume V2, with dosing directly into the",
    "central compartment (intravenous injection; no absorption phase). The single covariate",
    "retained by forward-inclusion / backward-elimination is the hourly net fluid volume infused",
    "during surgery, which enters Q as a power function; the authors attribute the effect to",
    "intra-operative changes in hepatic blood flow, nalbuphine having a hepatic extraction ratio",
    "of 0.5-0.7. Between-subject variability is exponential on all four disposition parameters",
    "and residual error is combined additive plus proportional.",
    sep = " "
  )
  reference <- paste(
    "Nie X, Gao X, Gao J, Heng T, Zhang Y, Sun Y, Feng Z, Jia L, Wang M.",
    "Population pharmacokinetics of nalbuphine in patients undergoing general anesthesia surgery.",
    "Front Pharmacol. 2023;14:1130287.",
    "doi:10.3389/fphar.2023.1130287.",
    sep = " "
  )
  vignette <- "Nie_2023_nalbuphine"

  # Nie 2023 reports plasma nalbuphine concentrations in ng/mL (LLOQ 0.1 ng/mL,
  # calibration range 0.1-500 ng/mL; Methods 2.2) and volumes in L, so the
  # amount unit that makes `central / vc` come out in ng/mL directly is the
  # microgram (ug/L == ng/mL). Doses are therefore expressed in ug: the 15 mg
  # clinical induction dose is 15000 ug and the 12 mg simulated fixed dose is
  # 12000 ug. Encoding the dose in mg instead would force the paper's printed
  # additive residual error (2.88 ng/mL) to be rewritten as 0.00288 mg/L, i.e.
  # a unit conversion applied to a published value, which this library avoids.
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")

  compartmentData <- list(
    central = list(
      analyte = "nalbuphine", units = "ug", specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "nalbuphine", units = "ug", specimen = "tissue", verified = TRUE
    )
  )

  covariateData <- list(
    PFA_NET_RATE = list(
      description        = "Hourly net fluid volume infused during surgery",
      units              = "mL/h",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column HNF. Nie 2023 Table 1 footnote defines it as",
        "HNF = (FVI + BVI - UVO) / OT, where FVI is the fluid volume infused, BVI the blood",
        "volume infused and UVO the urine volume output during surgery, and OT the operation",
        "time -- i.e. a net intra-operative fluid balance expressed as a rate. Time-fixed per",
        "subject (one scalar per operation). Model-building cohort 617.96 +/- 247.61 mL/h",
        "(median 563.56, range 234.26-1202.25); external-validation cohort 634.19 +/- 178.94",
        "mL/h (median 610.42, range 314.29-1047.24). The final model normalizes by 617.96 mL/h",
        "-- the model-building cohort MEAN, not the median. Nie 2023 Eq. 2 states the generic",
        "power form as theta1 * (cov_i / cov_median)^theta2, but the printed final-model",
        "equation (Results section 3.2) is Q = 245 * (HNF/617.96)^-0.58 and 617.96 is the mean",
        "from Table 1. The printed equation is used here per the standing rule that a printed",
        "equation outranks conflicting prose. Enters Q only; the negative exponent means a",
        "higher net intra-operative fluid load gives a LOWER intercompartmental clearance.",
        "One subject was dropped from the external-validation set for a missing HNF value",
        "(Results 3.1).",
        sep = " "
      ),
      source_name        = "HNF"
    )
  )

  # Screened during covariate model building (Nie 2023 Table 3) but NOT retained
  # in the final model, so they are documented rather than declared: none of
  # them is referenced in model(). ALT is the notable case -- it entered the
  # full model on CL (dOFV -7.39, p < 0.01) but was removed in backward
  # elimination because its dOFV of 7.40 did not exceed the 7.88 retention
  # criterion (Table 3 model 22, "> 0.005"), which the Discussion attributes to
  # the small sample and to most hepatic patients having ALT within 3x the
  # upper limit of normal.
  covariatesDataExcluded <- list(
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = paste(
        "Tested on CL as an exponential function (Nie 2023 Table 3 models 2 and 20, functional",
        "expression 'a' = Eq. 1 linear per the table footnote key). Entered the full model",
        "(OFV 2158.168 -> 2150.776 alone; 2148.529 -> 2141.133 on top of HNF-on-Q) but was",
        "eliminated in the backward step (dOFV 7.40 < 7.88). No point estimate is printed for",
        "the ALT coefficient anywhere in the paper, so the effect cannot be reconstructed.",
        sep = " "
      )
    ),
    GGT = list(
      description = "Gamma-glutamyltransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Tested on CL (Table 3 model 3, dOFV -4.90) and on V2 (model 13, dOFV -7.02); not retained. No coefficient printed."
    ),
    HR = list(
      description = "Heart rate",
      units       = "beats/min",
      type        = "continuous",
      notes       = "Tested on V1 (Table 3 model 5, dOFV -7.06) and on Q (model 8, dOFV -5.88); not retained. No coefficient printed."
    ),
    WT = list(
      description = "Total body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Tested on V2 only (Table 3 model 12, dOFV -4.07) and not retained; no allometric",
        "scaling appears anywhere in the final model. This is the paper's central dosing",
        "finding -- because body weight does not enter the PK, a fixed 12 mg dose and a",
        "0.2 mg/kg weight-based dose differ only by the cohort mean dose (bias < 6% at every",
        "sampled time, Table 6), and the fixed regimen carries the lower exposure variability",
        "(Figure 8). Weight is still needed as a COHORT attribute to simulate the weight-based",
        "arm, but it is not a model covariate.",
        sep = " "
      )
    ),
    UA = list(
      description = "Uric acid",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Tested on V2 (Table 3 model 14, dOFV -5.37); not retained. No coefficient printed."
    ),
    DDIMER = list(
      description = "Plasma D-dimer",
      units       = "mg/L",
      type        = "continuous",
      notes       = "Tested on V2 (Table 3 model 15, dOFV -4.79); not retained. No coefficient printed."
    ),
    SMOKE = list(
      description = "Current-smoker indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested on V2 (Table 3 model 19, dOFV -3.94); not retained. No coefficient printed. Table 1: 11/27 (40.7%) in the model-building cohort."
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 47L,
    n_studies       = 1L,
    age_range       = "21-78 years overall; model-building cohort 53.85 +/- 16.63 years (median 58, range 21-76), external-validation cohort 52.5 +/- 13.56 (median 54.5, range 27-78)",
    age_median      = "58 years (model-building cohort, Table 1)",
    weight_range    = "48-86 kg overall; model-building cohort 63.32 +/- 9.3 kg (median 62, range 48-82), external-validation cohort 62.25 +/- 9.67 (median 59.5, range 50.8-86)",
    weight_median   = "62 kg (model-building cohort, Table 1)",
    height_range    = "Model-building cohort 163.96 +/- 6.87 cm (median 165, range 153-175)",
    sex_female_pct  = 55.6,
    race_ethnicity  = "Not reported; single-centre Chinese cohort",
    disease_state   = "Adults scheduled for elective surgery under general anaesthesia. Liver resection 14.8%, cholecystectomy 12.8%, pancreatic resection 36.2%, other surgery 36.2%; 29.8% laparoscopic and 70.2% open. Model-building cohort: tumour 77.8%, hepatobiliary disease 55.6%, hypertension 25.9%, Child-Turcotte-Pugh class B 18.5% (the remainder class A).",
    dose_range      = "Single nalbuphine 15 mg intravenous injection over 2-3 min at induction (0.24 +/- 0.04 mg/kg), followed by midazolam 0.05 mg/kg, sufentanil 0.2 ug/kg, etomidate 0.03 mg/kg and cisatracurium 0.2 mg/kg; anaesthesia maintained with sevoflurane and remifentanil",
    regions         = "China (Shijiazhuang; the Fourth Hospital of Hebei Medical University, single centre, 2021)",
    hepatic_function = "55.6% of the model-building cohort had hepatobiliary disease; 18.5% were Child-Turcotte-Pugh class B and the rest class A. Excluded: known or suspected cardiopulmonary, renal or metabolic disease.",
    renal_function  = "Creatinine clearance (Cockcroft-Gault) 111.19 +/- 27.38 mL/min in the model-building cohort (median 105.14, range 53.4-160.98)",
    co_medication   = "Midazolam, sufentanil, etomidate, cisatracurium, sevoflurane and remifentanil, all part of the anaesthetic protocol. Nie 2023 Discussion notes as a limitation that nalbuphine, midazolam, sufentanil and sevoflurane are all CYP3A4 substrates and that the study design could not resolve any interaction between them.",
    notes           = "Baseline demographics per Nie 2023 Table 1. 458 concentrations were drawn from 48 patients; one subject (5 samples) was dropped for a missing HNF value, leaving 27 patients / 353 samples for model building (intensive sampling: pre-dose and 3, 5, 10, 15, 30, 45 min and 1, 1.5, 2, 3, 4, 5, 6, 12 h) and 20 patients / 100 samples for external validation (sparse sampling: pre-dose, during intubation, and 1, 3, 10 min after intubation). Every post-dose sample was above the 0.1 ng/mL LLOQ and none was flagged as an outlier."
  )

  ini({
    # Structural parameters, Nie 2023 Table 4 "Final model / Estimates" column.
    # The paper prints the final-model equations directly (Results 3.2):
    #   CL (L/h) = 32.9
    #   V1 (L)   = 32.5
    #   Q  (L/h) = 245 * (HNF/617.96)^-0.58
    #   V2 (L)   = 83.5
    # (Table 4's legend gives the units; the printed "Q (L)" line is a typo for
    # L/h, since Table 4 and the Discussion both state Q in L/h.)
    # All fixed-effect RSEs are < 30% and every estimate falls inside the
    # 2.5th-97.5th percentile of the 1000-sample bootstrap (984 successful
    # minimizations), so no parameter is poorly identified.
    lcl <- log(32.9); label("Clearance CL (L/h)")                                  # Nie 2023 Table 4 final model CL = 32.9 L/h (RSE 5.47%; bootstrap median 32.8, 95% CI 29.46-36.54)
    lvc <- log(32.5); label("Central volume V1 (L)")                               # Nie 2023 Table 4 final model V1 = 32.5 L (RSE 10.25%; bootstrap median 31.9, 95% CI 25.90-38.10)
    lq  <- log(245);  label("Intercompartmental clearance Q at PFA_NET_RATE = 617.96 mL/h (L/h)")  # Nie 2023 Table 4 final model Q = 245 L/h (RSE 13.99%; bootstrap median 247, 95% CI 216.63-297.38)
    lvp <- log(83.5); label("Peripheral volume V2 (L)")                            # Nie 2023 Table 4 final model V2 = 83.5 L (RSE 7.94%; bootstrap median 84.6, 95% CI 70.70-97.34)

    # Power-function exponent for the hourly net fluid volume infused on Q,
    # normalized to 617.96 mL/h. Applied in model() as
    # (PFA_NET_RATE / 617.96)^e_pfa_net_rate_q, matching the printed equation.
    e_pfa_net_rate_q <- -0.58; label("Power exponent of hourly net fluid volume infused on Q (unitless)")  # Nie 2023 Table 4 final model "HNF on Q" = -0.58 (RSE 14.11%; bootstrap median -0.557, 95% CI -0.802 to -0.084)

    # Between-subject variability, exponential on every disposition parameter
    # (Methods 2.3.1: "BSV ... was assumed to be log-normally distributed and
    # was applied by exponential model").
    #
    # Table 4 labels the BSV column "%CV", but the values printed there are the
    # NONMEM omega standard deviations x 100, not the exact log-normal CV. The
    # tell is the residual-error block of the SAME column: the final-model
    # proportional error prints as 0.139 while the basic model's (Table 2)
    # prints as 13.7 under the identical "%CV" header -- the software is
    # printing the raw parameter and the header is a mislabel. Reading the BSV
    # entries the same way, omega_CL = 0.277, omega_V1 = 0.401, omega_Q = 0.179
    # and omega_V2 = 0.408, so the variances below are those values squared.
    # This reading is corroborated by Results 3.2, which subtracts the numbers
    # arithmetically -- "a decrease in BSV from 35.5% to 17.9%, indicating that
    # 17.6% of BSV in Q was explained by HNF" (35.5 - 17.9 = 17.6) -- which is
    # only meaningful on the omega scale.
    etalcl ~ 0.076729  # Nie 2023 Table 4 final model omega CL = 27.7 -> 0.277^2 = 0.076729 (RSE 12.26%, shrinkage 4.30%)
    etalvc ~ 0.160801  # Nie 2023 Table 4 final model omega V1 = 40.1 -> 0.401^2 = 0.160801 (RSE 22.92%, shrinkage 13.10%)
    etalq  ~ 0.032041  # Nie 2023 Table 4 final model omega Q  = 17.9 -> 0.179^2 = 0.032041 (RSE 42.65%, shrinkage 30.72%; down from 35.5 in the basic model once HNF entered Q)
    etalvp ~ 0.166464  # Nie 2023 Table 4 final model omega V2 = 40.8 -> 0.408^2 = 0.166464 (RSE 13.84%, shrinkage 3.18%)

    # Residual unexplained variability: combined proportional + additive
    # (Methods 2.3.1 tested proportional, additive and combined; Results 3.2
    # selected the combination).
    propSd <- 0.139; label("Proportional residual error (fraction)")  # Nie 2023 Table 4 final model proportional error = 0.139 (RSE 8.02%, shrinkage 12.79%; the basic model prints the same quantity as 13.7 in Table 2)
    addSd  <- 2.88;  label("Additive residual error (ng/mL)")         # Nie 2023 Table 4 final model additive error = 2.88 ng/mL (RSE 19.81%; bootstrap median 2.9035, 95% CI 1.88-4.72)
  })

  model({
    # Individual disposition parameters. Only Q carries a covariate; CL, V1 and
    # V2 are covariate-free in the final model (weight and age were both
    # screened and rejected -- see covariatesDataExcluded and the Discussion).
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)

    # Q = 245 * (HNF/617.96)^-0.58 (Nie 2023 Results 3.2 printed equation).
    # 617.96 mL/h is the model-building cohort MEAN hourly net fluid volume
    # infused (Table 1), so exp(lq) is Q for a patient at that mean fluid load.
    q <- exp(lq + etalq) * (PFA_NET_RATE / 617.96)^e_pfa_net_rate_q

    # Two-compartment IV disposition micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Nalbuphine is given as an intravenous injection at induction, so doses
    # enter `central` directly from the event table; there is no depot and no
    # bioavailability term.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # Plasma nalbuphine concentration in ng/mL (amounts in ug, volumes in L).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
