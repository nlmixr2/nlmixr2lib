Standing_2012_oseltamivir <- function() {
  description <- "Population PK model for oral oseltamivir and its active metabolite oseltamivir carboxylate in preterm and term neonates and infants (Standing 2012). One-compartment parent + one-compartment metabolite with first-order absorption, an empirical transit compartment delaying first-pass metabolite appearance, well-stirred-model hepatic first-pass conversion (FM derived from CLI / liver-blood-flow FQ), and physiologically scaled clearances combining (WT/70)^0.75 allometry with a Rhodin 2009 renal-maturation Hill sigmoid on CLU/CLM and a fitted HCE1 Hill sigmoid (PM50 86.1 wk, Hill 3.17) on intrinsic clearance CLI. Volumes (VD, VDM) and liver blood flow (FQ) fixed from external references."
  reference   <- "Standing JF, Nika A, Tsagris V, Kapetanakis I, Maltezou HC, Kafetzis DA, Tsolia MN. Oseltamivir pharmacokinetics and clinical experience in neonates and infants during an outbreak of H1N1 influenza A virus infection in a neonatal intensive care unit. Antimicrob Agents Chemother. 2012;56(7):3833-3840. doi:10.1128/AAC.00290-12"
  vignette    <- "Standing_2012_oseltamivir"
  units       <- list(time = "hour", dosing = "nmol", concentration = "nmol/L")

  # Standing 2012 Fig 1 compartment 3 is an empirical transit compartment that
  # delays first-pass metabolite appearance (causing flip-flop kinetics). It is
  # not a Savic absorption-chain transit; declare as paper-specific.
  paper_specific_compartments <- "transit_oc"

  covariateData <- list(
    WT = list(
      description        = "Total body weight (kg). Used for allometric size scaling on every disposition parameter: (WT/70)^0.75 on CLU, CLM, CLI, and liver blood flow FQ; (WT/70)^1.0 on VD and VDM. Reference weight 70 kg.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying in principle but constant within an individual over the 4-sample sparse sampling window of the source study. Source paper Methods (Pharmacokinetic modeling paragraph): 'all volume parameters were scaled with linear weight, and clearance of unchanged oseltamivir and oseltamivir carboxylate was scaled with a renal maturation model published by Rhodin et al. ... For the intrinsic clearance of oseltamivir to oseltamivir carboxylate, weight scaling used a power of 0.75.'",
      source_name        = "WT"
    ),
    PAGE = list(
      description        = "Postmenstrual age in months (gestational age in weeks / 4.348 + postnatal months). Drives two maturation Hill sigmoids: a Rhodin (2009) renal maturation function applied to CLU and CLM, and a Standing-fitted HCE1 maturation function applied to intrinsic clearance CLI.",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Range in the Standing 2012 cohort (Table 1) converted to months: ~6.4-12.0 months PAGE (Standing reports postmenstrual age in weeks; here it is rendered in months to match the canonical PAGE convention). The model converts PAGE back to weeks internally for the maturation Hill sigmoids because Standing 2012 and Rhodin 2009 published their TM50 values in weeks. Renal maturation: F_renal(PMA_wk) = PMA_wk^3.4 / (47.7^3.4 + PMA_wk^3.4) from Rhodin et al. 2009 doi:10.1007/s00228-008-0577-4 Table 1 (referenced by Standing 2012 ref 23; the TM50 and Hill values are not reproduced in the Standing text). HCE1 maturation: F_HCE1(PMA_wk) = PMA_wk^3.17 / (86.1^3.17 + PMA_wk^3.17), fitted by Standing 2012 to HCE1 expression data from Yang et al. 2009 (Standing 2012 Fig 2 caption: PM50 = 86.1 weeks; Hill = 3.17).",
      source_name        = "PMA (weeks)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 9,
    n_studies      = 1,
    age_range      = "PMA 28-52 weeks (postnatal age 2-86 days; gestational age 24-40 weeks)",
    age_median     = "PMA 36 weeks",
    weight_range   = "1.22-3.35 kg",
    weight_median  = "2.37 kg",
    sex_female_pct = NA_real_,
    disease_state  = "Preterm and term neonates / infants in a level-3 NICU during an A(H1N1) influenza outbreak. One subject had laboratory-confirmed H1N1 and received treatment dosing; the remaining eight received prophylaxis. Comorbidities included prematurity, respiratory distress syndrome, patent ductus arteriosus, necrotizing enterocolitis, chronic lung disease, intraventricular hemorrhage, and congenital heart disease.",
    dose_range     = "Oseltamivir suspension 1 mg/kg every 12 h for PMA <= 37 weeks; 3 mg/kg every 12 h (treatment) or every 24 h (prophylaxis) for PMA > 37 weeks. Administered via nasogastric tube (capsule contents suspended in water per manufacturer instructions).",
    regions        = "Greece (P. & A. Kyriakou Children's Hospital NICU, Athens).",
    n_observations = "72 plasma observations total (36 oseltamivir + 36 oseltamivir carboxylate). 6 of 36 oseltamivir samples (5 troughs + one 8-h sample) were below the LOQ of 1 ng/mL (~3.20 nM) and were handled by the NONMEM M3 method in the original fit. No carboxylate samples were below the carboxylate LOQ of 10 ng/mL (~35.2 nM).",
    sampling       = "Predose and 1, 4, and 8 h postdose at presumed steady state (variable numbers of preceding doses across subjects).",
    notes          = "Demographics and dosing summary from Standing 2012 Table 1. Three subjects had necrotizing enterocolitis and three were receiving total parenteral nutrition without enteral feeding; none showed compromised carboxylate exposure relative to the cohort."
  )

  ini({
    # Estimated structural parameters (Standing 2012 Table 2, "Estimate" column).
    # Adult typical reference: 70 kg body weight, infinite PMA (fully mature renal
    # function and HCE1 expression, both maturation functions -> 1).
    lka  <- log(0.22);   label("Absorption rate constant for oseltamivir and oseltamivir carboxylate (1/h)")            # Standing 2012 Table 2: Ka = 0.22 /h
    lcl  <- log(30.1);   label("Renal clearance of parent and carboxylate at 70 kg, CLU = CLM (L/h)")                   # Standing 2012 Table 2: CLU = CLM = 30.1 L/h/70 kg (set equal per the Results paragraph)
    lkam <- log(0.034);  label("Carboxylate appearance rate constant from transit compartment (1/h)")                   # Standing 2012 Table 2: Kam = 0.034 /h
    lcli <- log(3284);   label("Intrinsic hepatic clearance of oseltamivir to carboxylate at 70 kg, CLI (L/h)")         # Standing 2012 Table 2: CLI = 3284 L/h/70 kg

    # Fixed structural parameters (Standing 2012 Results paragraph and Methods).
    # VD and VDM were fixed to literature values (He et al. 2008, Standing ref 11)
    # to stabilize the fit given sparse data and largely flat profiles. FQ is the
    # adult typical liver blood flow used in the well-stirred hepatic model.
    lvc    <- fixed(log(91));    label("Volume of distribution of parent at 70 kg, VD (L) [fixed from He 2008]")        # Standing 2012 Results: VD fixed to 91 L/70 kg (ref 11)
    lvc_oc <- fixed(log(25.6));  label("Volume of distribution of carboxylate at 70 kg, VDM (L) [fixed from He 2008]")  # Standing 2012 Results: VDM fixed to 25.6 L/70 kg (ref 11)
    lfq    <- fixed(log(75));    label("Liver blood flow at 70 kg, FQ (L/h) [fixed adult value]")                       # Standing 2012 Methods: FQ adult value 75 L/h/70 kg (ref 21)

    # Fixed allometric exponents (Standing 2012 Methods, Tod et al. 2008 scaling).
    allo_cl <- fixed(0.75); label("Allometric exponent on clearances (CLU, CLM, CLI, FQ)")                              # Standing 2012 Methods: 0.75 on clearance (Tod et al.)
    allo_vc <- fixed(1.0);  label("Allometric exponent on volumes (VD, VDM)")                                           # Standing 2012 Methods: linear weight scaling on volumes (exponent 1.0)

    # HCE1 maturation Hill-sigmoid parameters fitted by Standing 2012 to the
    # Yang et al. 2009 HCE1 expression data (Standing 2012 Fig 2 caption).
    pma50_hce1 <- fixed(86.1); label("PMA for 50% mature HCE1 expression (weeks)")                                      # Standing 2012 Fig 2 caption: PM50 = 86.1 weeks
    hill_hce1  <- fixed(3.17); label("Hill coefficient for HCE1 maturation")                                            # Standing 2012 Fig 2 caption: Hill = 3.17

    # Renal maturation Hill-sigmoid parameters from Rhodin et al. 2009 (cited
    # in Standing 2012 as ref 23 but the TM50 / Hill values are not reproduced
    # in the Standing text; they come from Rhodin 2009 Table 1 directly).
    pma50_renal <- fixed(47.7); label("PMA for 50% mature renal function (weeks) [Rhodin 2009]")                        # Rhodin et al. 2009 doi:10.1007/s00228-008-0577-4 Table 1; referenced by Standing 2012 ref 23 but values not in Standing text
    hill_renal  <- fixed(3.4);  label("Hill coefficient for renal maturation [Rhodin 2009]")                            # Rhodin et al. 2009 doi:10.1007/s00228-008-0577-4 Table 1

    # IIV. Standing 2012 Table 2 reports BSV as %CV (log-normal). Convert via
    # omega^2 = log(1 + CV^2): 57.5% -> 0.2858, 59.4% -> 0.3023, 71.0% -> 0.4082.
    # Standing's eta on CLU/CLM is shared across both clearances (single THETA,
    # single ETA) so the eta on lcl drives both parent and metabolite renal CL.
    etalka  ~ 0.2858    # Standing 2012 Table 2: BSV(Ka) = 57.5% CV; omega^2 = log(1 + 0.575^2)
    etalcl  ~ 0.3023    # Standing 2012 Table 2: BSV(CLU and CLM) = 59.4% CV (single eta shared); omega^2 = log(1 + 0.594^2)
    etalcli ~ 0.4082    # Standing 2012 Table 2: BSV(CLI) = 71.0% CV; omega^2 = log(1 + 0.710^2)

    # Residual error (proportional on each output). Standing 2012 Table 2.
    propSd    <- 0.543; label("Proportional residual error on parent oseltamivir (fraction)")                           # Standing 2012 Table 2: 54.3%
    propSd_oc <- 0.232; label("Proportional residual error on oseltamivir carboxylate (fraction)")                      # Standing 2012 Table 2: 23.2%
  })

  model({
    # Reference body weight for allometric scaling (Standing 2012 used 70 kg adult).
    ref_wt <- 70

    # Convert PAGE (canonical postmenstrual age in months) to weeks for the
    # maturation Hill sigmoids, because Standing 2012 (HCE1) and Rhodin 2009
    # (renal) published their TM50 values in weeks. Conversion 1 month =
    # 30.4375 / 7 = 4.348 weeks (matches the canonical PAGE schema in
    # inst/references/covariate-columns.md ## PAGE entry).
    pma_weeks <- PAGE * 4.348125

    # Maturation fractions as functions of postmenstrual age in weeks (Hill
    # sigmoids). Both saturate to 1 as PMA -> infinity, so the adult typical
    # value reduces to no maturation scaling.
    fmat_renal <- pma_weeks^hill_renal / (pma50_renal^hill_renal + pma_weeks^hill_renal)
    fmat_hce1  <- pma_weeks^hill_hce1  / (pma50_hce1^hill_hce1  + pma_weeks^hill_hce1)

    # Individual structural parameters (typical * IIV * size * maturation).
    ka    <- exp(lka  + etalka)
    cl    <- exp(lcl  + etalcl)  * (WT / ref_wt)^allo_cl * fmat_renal
    cli   <- exp(lcli + etalcli) * (WT / ref_wt)^allo_cl * fmat_hce1
    fq    <- exp(lfq)            * (WT / ref_wt)^allo_cl
    vc    <- exp(lvc)            * (WT / ref_wt)^allo_vc
    vc_oc <- exp(lvc_oc)         * (WT / ref_wt)^allo_vc
    kam   <- exp(lkam)

    # Well-stirred hepatic model (Standing 2012 Methods):
    #   CLTM = (CLI * FQ) / (CLI + FQ)   = systemic parent -> metabolite CL
    #   FM   = CLI / (CLI + FQ)          = first-pass extraction ratio
    cltm <- (cli * fq) / (cli + fq)
    fm   <- cli / (cli + fq)

    # ODE system (Standing 2012 Fig 1). Compartment 1 = depot (oral dose,
    # oseltamivir as administered), compartment 2 = central parent (oseltamivir),
    # compartment 3 = transit_oc (empirical transit delaying first-pass
    # metabolite appearance), compartment 4 = central_oc (oseltamivir carboxylate).
    # The systemic CLTM route (parent -> carboxylate via hepatic well-stirred
    # extraction) feeds central_oc directly; only the first-pass FM fraction
    # passes through the transit_oc compartment (consistent with the Discussion
    # description of Kam as a 'mean absorption time' affected by cholestasis
    # and gut physiology).
    d/dt(depot)      <- -ka * depot
    d/dt(central)    <-  (1 - fm) * ka * depot - (cl / vc) * central - (cltm / vc) * central
    d/dt(transit_oc) <-  fm * ka * depot - kam * transit_oc
    d/dt(central_oc) <-  kam * transit_oc + (cltm / vc) * central - (cl / vc_oc) * central_oc

    # Observations and residual error. Both concentrations are in nM
    # (Standing 2012 transformed doses and observations to molar units using
    # MW oseltamivir = 312.40 g/mol and MW carboxylate = 284.35 g/mol).
    Cc    <- central    / vc
    Cc_oc <- central_oc / vc_oc

    Cc    ~ prop(propSd)
    Cc_oc ~ prop(propSd_oc)
  })
}
