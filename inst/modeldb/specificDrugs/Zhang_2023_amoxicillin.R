Zhang_2023_amoxicillin <- function() {
  description <- paste(
    "Two-compartment intravenous population PK model for amoxicillin, as implemented by Zhang 2023 to",
    "simulate steady-state free-drug exposure in a 10,000-subject virtual population spanning creatinine",
    "clearance 10-150 mL/min/1.73 m2, for a pharmacodynamic analysis of mutant selection by the",
    "aztreonam/amoxicillin/clavulanate combination against NDM- and serine-beta-lactamase co-producing",
    "Escherichia coli and Klebsiella pneumoniae. Clearance is LINEARLY proportional to BSA-normalized",
    "creatinine clearance referenced to 102 mL/min/1.73 m2 (no power exponent); the volumes and",
    "intercompartmental clearance carry no covariates. Three algebraic observables are exposed: total",
    "plasma concentration Cc, unbound plasma concentration Ccfree = Cc * 0.83 (amoxicillin is 17% protein",
    "bound), and epithelial-lining-fluid concentration Celf = Cc * 0.30, which the paper treats as",
    "entirely unbound because no mucin binding was assumed. Every parameter is a FIXED prior transcribed",
    "from the upstream critically-ill analysis; Zhang 2023 estimated nothing and reported no residual",
    "error, so the residual SD is carried structurally at zero (see vignette Errata). The clavulanate",
    "co-formulant is a separate model, modellib('Zhang_2023_clavulanicAcid').",
    sep = " "
  )
  reference <- paste(
    "Zhang J, Wu M, Diao S, Zhu S, Song C, Yue J, Martins FS, Zhu P, Lv Z, Zhu Y, Yu M, Sy SKB.",
    "Pharmacokinetic/pharmacodynamic evaluation of aztreonam/amoxicillin/clavulanate combination against",
    "New Delhi metallo-beta-lactamase and serine-beta-lactamase co-producing Escherichia coli and",
    "Klebsiella pneumoniae. Pharmaceutics. 2023;15(1):251. doi:10.3390/pharmaceutics15010251",
    "(Section 2.4; Supplementary Materials 'Detailed description of pharmacokinetic models').",
    "Structural model and parameter estimates originate from Carlier M, Noe M, De Waele JJ, Stove V,",
    "Verstraete AG, Lipman J, Roberts JA. Population pharmacokinetics and dosing simulations of",
    "amoxicillin/clavulanic acid in critically ill patients. J Antimicrob Chemother.",
    "2013;68(11):2600-2608. doi:10.1093/jac/dkt240",
    sep = " "
  )
  vignette <- "Zhang_2023_aztreonam_amoxicillin_clavulanate"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central     = list(analyte = "amoxicillin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "amoxicillin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance, normalized to 1.73 m2 body surface area",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The ONLY covariate in this model, entering clearance as a LINEAR proportionality:",
        "CL = 10.3*(CRCL/102), with no power exponent. Both renderings of the source agree - Section 2.4",
        "prints 'CL = 10.3*(CrCL 102)' and the Supplementary Materials equation prints",
        "'CL = 10.3(CrCL/102) exp(omega)'. Neither carries an exponent, whereas the aztreonam equations",
        "in the same two places do print their (CrCL/100)^0.43 and (WT/70)^1.99 exponents, which",
        "confirms that the absence here is real and not a lost superscript. Zhang 2023 sampled CRCL from",
        "a uniform distribution over 10-150 mL/min normalized to 1.73 m2, drawn separately per renal",
        "function category (51-150, 31-50, 10-30) for the three dosing arms of Table 1.",
        "Source variable name CrCL.",
        sep = " "
      ),
      source_name        = "CrCL"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 10000,
    n_studies      = 1,
    sex_female_pct = 50,
    renal_function = "creatinine clearance 10-150 mL/min/1.73 m2, uniformly distributed",
    disease_state  = paste(
      "simulated adults; the source population PK analysis was conducted in critically ill patients",
      sep = " "
    ),
    dose_range     = paste(
      "amoxicillin/clavulanate as 3 h intravenous infusions by renal function category (Table 1):",
      "1/0.2 g q6h for CRCL > 50-150 and > 30-50 mL/min; 1/0.2 g followed by 500/100 mg q12h for",
      "CRCL 10-30 mL/min. Amoxicillin doses are therefore 1000 mg or 500 mg.",
      sep = " "
    ),
    notes          = paste(
      "The 10,000 subjects are a VIRTUAL population constructed by Zhang 2023 (Section 2.4 and the",
      "Supplementary Materials script), not an observed cohort: sex 50/50 and creatinine clearance",
      "sampled uniformly within each renal function category. Body weight was simulated for the",
      "aztreonam model but is not a covariate here, so it does not affect this model. No demographic",
      "table is reported because no patients were enrolled for the PK component. The parameter estimates",
      "come from the upstream critically-ill analysis cited in the reference field.",
      sep = " "
    )
  )

  ini({
    # Every parameter is FIXED: Zhang 2023 re-used the upstream population PK
    # estimates as simulation priors and estimated nothing.

    lcl <- fixed(log(10.3)); label("Clearance at CRCL 102 mL/min/1.73 m2 (L/h)")
    # Section 2.4 and Suppl. Materials: CL = 10.3*(CrCL/102)

    lvc <- fixed(log(13.5)); label("Central volume (L)")
    # Section 2.4 and Suppl. Materials: VC = 13.5, CV 38.7%

    lvp <- fixed(log(14.1)); label("Peripheral volume (L)")
    # Section 2.4 and Suppl. Materials: VP = 14.1, no variability reported

    lq  <- fixed(log(15.7)); label("Intercompartmental clearance (L/h)")
    # Section 2.4 and Suppl. Materials: Q = 15.7, no variability reported

    # Protein binding and lung penetration, used to build the free-drug and
    # epithelial-lining-fluid observables the pharmacodynamic analysis runs on.
    fu   <- fixed(0.83); label("Fraction of amoxicillin unbound in plasma (unitless)")
    # Section 2.4: amoxicillin plasma protein binding is 17%, so fu = 1 - 0.17

    relf <- fixed(0.30); label("Epithelial lining fluid to TOTAL plasma concentration ratio (unitless)")
    # Sections 2.4 and 3.4 both give a 30% ELF penetration rate for
    # amoxicillin/clavulanate. As for aztreonam, the ratio is applied to TOTAL
    # (not unbound) plasma concentration; see the vignette Errata.

    # IIV. The paper reports these as CV%; the Supplementary Materials write the
    # same numbers as omega ~ N(0, value) alongside the aztreonam script that
    # passes them to rnorm(sd = ...), so they are log-scale SDs and the omega
    # entries below are their squares.
    # Section 2.4 and Suppl. Materials: CL CV 39.9%, so 0.399^2 = 0.159201
    etalcl ~ fixed(0.159201)
    # Section 2.4 and Suppl. Materials: VC CV 38.7%, so 0.387^2 = 0.149769
    etalvc ~ fixed(0.149769)
    # No variability is reported on VP or Q anywhere in the paper or supplement;
    # carried structurally at zero.
    etalvp ~ fixed(0)
    etalq  ~ fixed(0)

    # Residual error is not reported: Zhang 2023 simulated concentration-time
    # profiles from between-subject variability only and never fitted data.
    propSd <- fixed(0); label("Proportional residual SD (fraction; not reported in the source)")
  })

  model({
    # Individual parameters. The CRCL term is linear, not a power function.
    cl <- exp(lcl + etalcl) * (CRCL / 102)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q  <- exp(lq  + etalq)

    # Two-compartment intravenous system, written as in the Supplementary
    # Materials ODEs (A1 = central, A2 = peripheral1)
    d/dt(central)     <- -central * (cl / vc + q / vc) + peripheral1 * (q / vp)
    d/dt(peripheral1) <-  central * (q / vc)            - peripheral1 * (q / vp)

    # Observables. Cc is total plasma concentration (A1 / VC). Ccfree is the
    # unbound plasma concentration that drives the plasma fT>MIC / fT>MPC
    # analysis. Celf is the epithelial-lining-fluid concentration; it scales
    # TOTAL plasma concentration and is itself treated as entirely unbound,
    # because Section 2.4 assumed no mucin binding in the ELF.
    Cc     <- central / vc
    Ccfree <- Cc * fu
    Celf   <- Cc * relf

    Cc ~ prop(propSd)
  })
}
