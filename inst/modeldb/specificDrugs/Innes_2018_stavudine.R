Innes_2018_stavudine <- function() {
  description <- paste(
    "Two-compartment population PK model for intracellular stavudine",
    "triphosphate in peripheral blood mononuclear cells (PBMC) after oral",
    "stavudine, in 23 HIV-positive prepubertal children (4-11 years) on",
    "0.5 mg/kg twice daily and 24 HIV-positive adults on 20 mg twice daily",
    "(Cape Town, South Africa). Stavudine dosing is linked directly to the",
    "intracellular triphosphate by first-order appearance (ka) into a",
    "biphasic disposition model; plasma stavudine was not measured, so all",
    "clearances and volumes are apparent and expressed in cell units",
    "(10^12 cells/h and 10^12 cells), with relative bioavailability fixed",
    "to 1. Clearances and volumes are allometrically scaled on fat-free",
    "mass with fixed exponents 0.75 and 1, referenced to the adult median",
    "FFM of 47.2 kg. Between-subject variability on clearance;",
    "between-occasion variability on bioavailability and ka. No effect of",
    "age (maturation), sex or renal function was found."
  )
  reference <- paste(
    "Innes S, van der Laan L, Anderson PL, Cotton M, Denti P. Can we improve",
    "stavudine's safety profile in children? Pharmacokinetics of",
    "intracellular stavudine triphosphate with reduced dosing. Antimicrob",
    "Agents Chemother. 2018;62(11):e00761-18. doi:10.1128/AAC.00761-18.",
    "Parameter values from Table 2; structure from Results 'Pharmacokinetic",
    "model' and Methods 'Population pharmacokinetic modeling'.",
    sep = " "
  )
  vignette <- "Innes_2018_stavudine"
  units <- list(
    time = "h",
    dosing = "mg (stavudine)",
    concentration = "fmol/10^6 cells (intracellular stavudine triphosphate)"
  )

  covariateData <- list(
    FFM = list(
      description = "Fat-free mass, the body-size descriptor for allometric scaling of all clearances and volumes",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Reference value 47.2 kg, the adult cohort median (Table 2 footnote a:",
        "'the typical values reported here refer to an adult with a median",
        "fat-free mass of 47.2 kg'). Computed by the Al-Sallami et al. 2015",
        "equation (Methods reference 28; Clin Pharmacokinet 54:1169-1178),",
        "which applies a paediatric age correction to the Janmahasatian adult",
        "form. Fat-free mass was preferred over total body weight and body",
        "surface area (Methods 'Population pharmacokinetic modeling').",
        "Exponents fixed at 0.75 on CL and Q, and 1 on Vc and Vp. Children",
        "median FFM 19.1 kg.",
        sep = " "
      ),
      source_name = "FFM"
    ),
    OCC = list(
      description = "Occasion indicator for the between-occasion variability on bioavailability and ka",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Two occasions. Methods 'Population pharmacokinetic modeling': 'The",
        "pharmacokinetic samples collected predose were treated as a separate",
        "occasion in the model'. OCC = 1 is the dose preceding the day-7",
        "predose sample (the previous evening's dose); OCC = 2 is the",
        "observed day-7 morning dose whose post-dose samples (1, 2, 6 h or",
        "3, 4, 8 h) were collected. OCC is decomposed inside model() into",
        "binary indicators oc1 / oc2 that multiplex the per-occasion etas.",
        "Records with OCC outside {1, 2} carry no between-occasion",
        "variability.",
        sep = " "
      ),
      source_name = "OCC"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "stavudine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(
      analyte = "stavudine triphosphate",
      units = "mg (stavudine equivalents)",
      specimen = "blood cell",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "stavudine triphosphate",
      units = "mg (stavudine equivalents)",
      specimen = "blood cell",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 47,
    n_studies = 1,
    n_observations = 188,
    age_range = "children 4-11 years (median 8); adults 26-51 years (median 36)",
    weight_range = "children 18-31 kg (median 23); adults 62-158 kg (median 83)",
    ffm_range = "children 14-25 kg (median 19); adults 39-77 kg (median 47)",
    sex_female_pct = 72.3,
    race_ethnicity = "not reported; all participants enrolled in Cape Town, South Africa",
    disease_state = "HIV-1 infection, virally suppressed (>12 months) on first-line antiretroviral therapy; children prepubertal (Tanner stage 1 or 2)",
    dose_range = "stavudine 0.5 mg/kg twice daily (children, 1 mg/mL oral solution) or 20 mg twice daily (adults, capsules) for 7 days",
    regions = "South Africa",
    co_medication = "children: efavirenz, lamivudine and abacavir; adults: efavirenz, lamivudine and tenofovir",
    notes = paste(
      "Baseline characteristics in Table 1 (24 adults, 1 male; 23 children,",
      "12 male; one further child excluded for non-adherence). Sparse",
      "sampling on day 7 at predose and either 1, 2 and 6 h or 3, 4 and 8 h",
      "post-dose. Intracellular stavudine triphosphate measured by LC-MS/MS",
      "in PBMC; LLOQ 3 to 3.47 fmol/10^6 cells, 2 samples below the LLOQ",
      "handled by the M3 method. Fitted in Monolix 2016R1 (SAEM).",
      sep = " "
    )
  )

  ini({
    lcl <- log(454);   label("Apparent clearance of stavudine triphosphate at FFM 47.2 kg (10^12 cells/h)")                  # Innes 2018 Table 2 'CL (10^12 cells/h) 454 (9)'
    lvc <- log(2569);  label("Apparent central volume at FFM 47.2 kg (10^12 cells)")                                       # Innes 2018 Table 2 'Vc (10^12 cells) 2,569 (18)'
    lq  <- log(469);   label("Apparent intercompartmental clearance at FFM 47.2 kg (10^12 cells/h)")                      # Innes 2018 Table 2 'Q (10^12 cells/h) 469 (30)'
    lvp <- log(7500);  label("Apparent peripheral volume at FFM 47.2 kg (10^12 cells)")                                    # Innes 2018 Table 2 'Vp (10^12 cells) 7,500 (58)'
    lka <- log(1.17);  label("First-order appearance rate constant of stavudine triphosphate (1/h)")                       # Innes 2018 Table 2 'ka (1/h) 1.17 (46)', estimated with a weakly informative prior (footnote c)

    lfdepot <- fixed(log(1)); label("Relative bioavailability of a typical patient (fraction)")                             # Innes 2018 Table 2 'F 1 (fixed)'

    e_ffm_cl <- fixed(0.75); label("Allometric exponent of fat-free mass on CL and Q (unitless)")                           # Innes 2018 Methods: 'fixed coefficients of 3/4 and 1 was applied to all clearance (CL, Q) and volume of distribution (Vc, Vp) parameters'
    e_ffm_vc  <- fixed(1);    label("Allometric exponent of fat-free mass on Vc and Vp (unitless)")                         # Innes 2018 Methods, same sentence

    etalcl ~ 0.092416  # Innes 2018 Table 2 'IIV CL (%CV) 30.4 (21)', approximate CV = omega, so omega^2 = 0.304^2

    etaiov_fdepot_1 ~ 0.131044         # Innes 2018 Table 2 'IOV F (%CV) 36.2 (20)', omega^2 = 0.362^2
    etaiov_fdepot_2 ~ fixed(0.131044)  # one IOV variance shared by both occasions
    etaiov_ka_1     ~ 0.093636         # Innes 2018 Table 2 'IOV ka (%CV) 30.6 (61)', omega^2 = 0.306^2
    etaiov_ka_2     ~ fixed(0.093636)  # one IOV variance shared by both occasions

    addSd  <- fixed(0.6); label("Additive residual error (fmol/10^6 cells)")                                               # Innes 2018 Table 2 'Additive error (fmol/10^6 cells) 0.6 (fixed)', 20% of the median LLOQ of 3 fmol/10^6 cells (footnote d)
    propSd <- 0.274;      label("Proportional residual error (fraction)")                                                  # Innes 2018 Table 2 'Proportional error (%) 27.4 (8)'
  })

  model({
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2

    ffm_cl <- (FFM / 47.2)^e_ffm_cl
    ffm_v <- (FFM / 47.2)^e_ffm_vc

    cl <- exp(lcl + etalcl) * ffm_cl
    vc <- exp(lvc) * ffm_v
    q <- exp(lq) * ffm_cl
    vp <- exp(lvp) * ffm_v
    ka <- exp(lka + iov_ka)
    fdepot <- exp(lfdepot + iov_fdepot)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot

    # Dose in mg stavudine; amounts in mg stavudine equivalents; volumes in
    # 10^12 cells. Molar 1:1 conversion with the stavudine molecular weight
    # 224.21 g/mol (C10H12N2O4): 1 mg / 10^12 cells = 1e-3 / 224.21 * 1e15 fmol
    # per 10^6 * 10^6 cells = 1e6 / 224.21 fmol/10^6 cells.
    Cc <- 1e6 * central / (224.21 * vc)

    # Monolix combined error; the paper names a 'combined additive and
    # proportional model' without the form, so the Monolix default
    # combined1, y = f + (a + b * f) * e, is used.
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
