# Simultaneous parent + active-metabolite (desbutyl-lumefantrine, DLF)
# population pharmacokinetic model for oral artemether-lumefantrine in
# pregnant women on the Thailand-Myanmar border with uncomplicated
# Plasmodium falciparum malaria (Kloprogge 2015, Antimicrob Agents
# Chemother 59:6375-6384; doi:10.1128/AAC.00267-15).

Kloprogge_2015_lumefantrine <- function() {
  description <- paste(
    "Simultaneous parent + active-metabolite (desbutyl-lumefantrine,",
    "DLF) population PK model for oral lumefantrine in 116 pregnant",
    "women (second or third trimester) with uncomplicated Plasmodium",
    "falciparum malaria on the Thailand-Myanmar border treated with",
    "the standard fixed-dose artemether-lumefantrine regimen",
    "(Kloprogge 2015). First-order absorption with lag time into a",
    "two-compartment LF disposition with relative bioavailability F",
    "fixed at 1 and Box-Cox-transformed IIV on F (Box-Cox shape -0.394",
    "not encoded -- see Errata); DLF is formed mole-for-mole from LF",
    "central elimination (linear drug-metabolite chain, fraction",
    "metabolised assumed = 1) and disposes through its own two-",
    "compartment chain with apparent CL/F = 197 L/h and Vc/F = 6,490",
    "L. Retained covariates: estimated gestational age (power on LF",
    "ka, linear-deviation on LF Q/F, both centered on the cohort",
    "median 22.8 weeks) and admission parasitaemia (log10-exponential",
    "on DLF CL/F, centered on cohort median log10(3,260) = 3.513).",
    "Venous-only residual error encoded; capillary-residual variance",
    "components and capillary conversion factors (LF 0.878, DLF",
    "0.464) NOT encoded -- see Errata. Time-to-event PD layer",
    "(Gompertz hazard with E_max LF effect on recrudescent malaria,",
    "Table 4) NOT encoded -- see Errata. Parameter values from",
    "Kloprogge 2015 Table 2.",
    sep = " "
  )
  reference <- paste(
    "Kloprogge F, McGready R, Hanpithakpong W, Blessborn D, Day NPJ,",
    "White NJ, Nosten F, Tarning J (2015). Lumefantrine and Desbutyl-",
    "Lumefantrine Population Pharmacokinetic-Pharmacodynamic",
    "Relationships in Pregnant Women with Uncomplicated Plasmodium",
    "falciparum Malaria on the Thailand-Myanmar Border. Antimicrobial",
    "Agents and Chemotherapy 59(10):6375-6384.",
    "doi:10.1128/AAC.00267-15.",
    sep = " "
  )
  vignette <- "Kloprogge_2015_lumefantrine"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  paper_specific_compartments <- c("central_dlf", "peripheral1_dlf")
  paper_specific_etas         <- c("etalcl_dlf", "etalvc_dlf")
  paper_specific_residual_sds <- c("propSd_dlf")

  covariateData <- list(
    GA = list(
      description        = paste(
        "Estimated gestational age at study admission. Time-fixed per",
        "subject in this study (the 3-day artemether-lumefantrine",
        "course is short enough that gestational age does not",
        "materially advance during PK sampling)."
      ),
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Kloprogge 2015 Table 1: all-cohort median 22.8 weeks, range",
        "13.1-39.0 weeks across the 116 pregnant women. Per operator",
        "sidecar (2026-06-07), the existing GA canonical (originally",
        "defined as gestational age at birth in pediatric popPK use",
        "cases like Hu 2026 clesrovimab and Clegg 2024 nirsevimab) is",
        "reused for this maternal at-admission semantic, with the",
        "distinction recorded here in the notes rather than via a",
        "separate EGA canonical. Kloprogge 2015's NONMEM column is",
        "EGA; the values map directly to GA (weeks, no transformation",
        "of the column value). Applied to two PK parameters: power",
        "form on LF ka (ka_indiv = TVka * (GA / 22.8)^e_ga_ka with",
        "e_ga_ka = -0.715, so ka decreases with later gestational age",
        "-- consistent with reduced gut motility in late pregnancy)",
        "AND linear-deviation form on LF Q/F (q_indiv = TVQ_LF * (1 +",
        "e_ga_q * (GA - 22.8)) per week with e_ga_q = -0.0271 per",
        "week, corresponding to the paper's Table 2 -2.71% per week",
        "centered at the cohort median).",
        sep = " "
      ),
      source_name        = "EGA"
    ),
    PARA = list(
      description        = "Plasmodium falciparum parasitaemia at admission",
      units              = "parasites/uL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Admission-only (time-fixed) asexual parasite count per",
        "microlitre. Kloprogge 2015 Table 1: all-cohort median 3,260",
        "parasites/uL, range 56.8-154,000 across the 116 pregnant",
        "women. Applied as a log10-transformed exponential effect on",
        "DLF elimination clearance:",
        "cl_dlf_indiv = TVcl_dlf * exp(e_para_cl_dlf *",
        "(log10(max(PARA, 1)) - log10(3260))) with",
        "e_para_cl_dlf = +0.133 per log10 parasitaemia. Centered on",
        "log10(3260) = 3.513 (the cohort median). The log10 transform",
        "inside model() follows the Mahidol-Oxford / Joel Tarning",
        "lab convention also used in Kloprogge_2014_quinine and",
        "Kloprogge_2018_lumefantrine (see inst/references/covariate-",
        "columns.md PARA entry); the max(PARA, 1) gate keeps the",
        "exponential finite for parasitaemia values below detection",
        "(PARA <= 0 collapses to log10(1) = 0, giving",
        "F_cl_dlf = exp(0.133 * (0 - 3.513)) = exp(-0.467) = 0.627",
        "as a fallback typical-low-parasitaemia value). For downstream",
        "simulations of no-parasite-effect typical PK, set PARA =",
        "3260 to recover the unmodified TVcl_dlf. Paper Discussion:",
        "'the mechanism underlying this covariate relationship is not",
        "known' -- the +13.3% per log10 increase in apparent CL_DLF",
        "is documented empirically without a proposed mechanism.",
        sep = " "
      ),
      source_name        = "PARA"
    )
  )

  population <- list(
    species             = "human",
    n_subjects          = 116L,
    n_studies           = 1L,
    n_pregnant          = 116L,
    n_dense_venous      = 13L,
    n_sparse_capillary  = 103L,
    pregnancy_trimester = "second or third trimester (median gestational age 22.8 weeks at admission; range 13.1-39.0 weeks)",
    age_range           = "14-42 years (all-cohort range, Table 1)",
    age_median          = "24 years (Table 1)",
    weight_range        = "35.0-65.0 kg (all-cohort range, Table 1)",
    weight_median       = "49.0 kg (Table 1)",
    sex_female_pct      = 100,
    body_temp_range     = "35.0-39.3 degC (all-cohort range, Table 1)",
    body_temp_median    = "36.7 degC (Table 1)",
    parasitemia_range   = "56.8-154,000 parasites/uL (all-cohort range, Table 1)",
    parasitemia_median  = "3,260 parasites/uL (Table 1)",
    primiparity_pct     = 23.4,
    disease_state       = paste(
      "Uncomplicated Plasmodium falciparum malaria. Primary",
      "presentations: 17 with PCR-confirmed recrudescent infections,",
      "22 with novel infections during follow-up, 77 with no parasite",
      "reappearance during 42-day follow-up or until delivery. Median",
      "time to recrudescent malaria 23 days (range 14-63); median",
      "time to novel infection 35 days (range 15-140). Inclusion was",
      "initially restricted to women who had reappearance of P.",
      "falciparum after quinine monotherapy (the first-line",
      "treatment at the time); after data-safety-monitoring-committee",
      "review of the first 20 delivered women, expanded to all",
      "patients with P. falciparum in the second or third trimester."
    ),
    dose_range          = paste(
      "Coartem (Novartis): 20 mg artemether + 120 mg lumefantrine per",
      "tablet; 4 tablets per dose (80 mg artemether + 480 mg",
      "lumefantrine) twice daily for 3 days (oral; dose times 0, 8,",
      "24, 36, 48, 60 hours) with 200-250 mL chocolate milk (6-7 g",
      "fat) to optimise lumefantrine bioavailability."
    ),
    regions             = "Thailand-Myanmar border (Shoklo Malaria Research Unit antenatal clinic, Mae Sot)",
    trial_registration  = "ISRCTN 86353884 (controlled-trials.com)",
    notes               = paste(
      "Pooled analysis combining a dense venous sampling sub-study",
      "(n = 13, frequent samples at pre-last-dose and 0.5, 1, 2, 4,",
      "6, 8, 12, 24, 48, 72, 96, 120, 144, 168 h after the last",
      "dose) with a sparse capillary sampling sub-study (n = 103,",
      "fingerprick samples at random windows 0-72 h, 72-96 h, 96-",
      "144 h, 144-336 h after first dose plus day 7), both within a",
      "larger drug-efficacy trial. PK assay: HPLC-UV for venous",
      "LF/DLF and capillary LF; LC-MS/MS reanalysis after 5 years",
      "for capillary DLF (LLOQ 9.7 ng/mL LF and 1.0 ng/mL DLF on 100",
      "uL plasma). Model built in NONMEM v7.2 with ADVAN6, FOCEI for",
      "PK and Laplacian for time-to-event PD. Total assay records:",
      "688 LF / 517 DLF (207 LF + 175 DLF venous; 481 LF + 342 DLF",
      "capillary). Paper Discussion includes a Gompertz hazard time-",
      "to-event PD layer (Table 4) with sigmoidal E_max LF",
      "concentration effect on recrudescent malaria (baseline hazard",
      "0.0845 events/week, hazard half-life 400 h, IC50 169 ng/mL);",
      "NOT encoded in this PK-only model file -- see vignette",
      "Errata. Eta-shrinkage 28.1-45.5%; epsilon-shrinkage 4.76%",
      "(venous LF) / 21.0% (capillary LF) / 3.90% (venous DLF) /",
      "28.8% (capillary DLF) per paper Results."
    )
  )

  ini({
    # ----- Structural PK parameters (Kloprogge 2015 Table 2 'Population
    # estimate (%RSE)' column with 95% bootstrap CI in column 2). The
    # paper reports linear-scale values; log() is applied here for the
    # nlmixr2 internal log scale. -----
    lka     <- log(0.0577)                                                           # Table 2: ka                = 0.0577    (RSE  7.93%; 95% CI 0.0526-0.0655)
    label("LF absorption rate constant ka (1/h, typical at GA = 22.8 wk)")
    ltlag   <- log(1.31)                                                             # Table 2: Lag time          = 1.31      (RSE 43.6%;  95% CI 0.0131-2.05)
    label("LF absorption lag time tlag (h)")

    # Relative bioavailability F is fixed at 1 by the source paper as a
    # structural anchor; all variability around F is captured via the
    # log-normal IIV (CV 51.2%) below. The paper additionally reports
    # a Box-Cox shape parameter on F = -0.394 that mildly distorts the
    # F IIV distribution away from strict log-normality; this Box-Cox
    # departure is NOT reproduced in the encoded log-normal IIV here
    # (matches the same-author Kloprogge_2018 file's precedent; see
    # vignette Errata).
    lfdepot <- fixed(log(1))                                                         # Table 2: F                 = 100 (fixed)
    label("LF relative bioavailability F (unitless, fixed at 1)")

    lcl     <- log(5.35)                                                             # Table 2: CL_LF/F           = 5.35      (RSE 12.9%;  95% CI 4.11-6.77)
    label("LF apparent elimination clearance CL/F (L/h)")
    lvc     <- log(28.4)                                                             # Table 2: VC_LF/F           = 28.4      (RSE 26.8%;  95% CI 17.3-46.8)
    label("LF apparent central volume of distribution Vc/F (L)")
    lq      <- log(1.55)                                                             # Table 2: Q_LF/F            = 1.55      (RSE 13.9%;  95% CI 1.17-2.00)
    label("LF apparent intercompartmental clearance Q/F (L/h, typical at GA = 22.8 wk)")
    lvp     <- log(147)                                                              # Table 2: VP_LF/F           = 147       (RSE 13.9%;  95% CI 110-187)
    label("LF apparent peripheral volume of distribution Vp/F (L)")

    lcl_dlf <- log(197)                                                              # Table 2: CL_DLF/F          = 197       (RSE 11.5%;  95% CI 156-245)
    label("DLF apparent elimination clearance CL/F (L/h, typical at log10(PARA) = 3.513)")
    lvc_dlf <- log(6490)                                                             # Table 2: VC_DLF/F          = 6,490     (RSE 21.1%;  95% CI 3,460-8,970)
    label("DLF apparent central volume of distribution Vc/F (L)")
    lq_dlf  <- log(250)                                                              # Table 2: Q_DLF/F           = 250       (RSE 19.1%;  95% CI 183-369)
    label("DLF apparent intercompartmental clearance Q/F (L/h)")
    lvp_dlf <- log(13200)                                                            # Table 2: VP_DLF/F          = 13,200    (RSE 12.6%;  95% CI 10,500-16,900)
    label("DLF apparent peripheral volume of distribution Vp/F (L)")

    # ----- Covariate-effect parameters (Kloprogge 2015 Table 2 lower
    # rows). The paper documents the three retained-covariate forms in
    # Table 2 footnote (e): 'Exponential covariate relationship is
    # determined as exp{theta * [covariate - median(covariate)]}; power
    # covariate relationship is determined as [covariate/
    # median(covariate)]^theta; linear covariate relationship is
    # determined as 1 + {theta * [covariate - median(covariate)]}.'
    # See model() for the per-parameter application. -----
    e_ga_ka       <- -0.715                                                          # Table 2: EGA power^e on ka                = -0.715    (RSE 19.1%; 95% CI -0.972 to -0.474)
    label("Power effect of GA on LF ka, centered at cohort median 22.8 wk: ka_indiv = TVka * (GA / 22.8)^e_ga_ka")
    e_ga_q        <- -0.0271                                                         # Table 2: EGA linear^e on Q LF (%)         = -2.71     (RSE 22.8%; 95% CI -3.59 to -1.34) -- '(%)' suffix in source means percentage points, divided by 100 here to give the per-week fractional change (-2.71% per week == -0.0271 per week)
    label("Linear-deviation effect of GA on LF Q/F per week from cohort median 22.8 wk: q_indiv = TVq * (1 + e_ga_q * (GA - 22.8))")
    e_para_cl_dlf <- 0.133                                                           # Table 2: Parasitemia exponentially^e on CL_DLF = 0.133  (RSE 31.8%; 95% CI 0.0455-0.210) -- log10 transform applied inside model() per same-author Kloprogge 2014 / 2018 lab convention (inst/references/covariate-columns.md PARA entry)
    label("Exponential effect of log10(PARA) on DLF CL/F per log10 parasites/uL, centered at log10(3260) = 3.513")

    # ----- Inter-individual variability (Kloprogge 2015 Table 2
    # 'IIV (%CV)' column). Table 2 footnote (b): the reported %CV
    # for IIV is computed as 100 * sqrt(exp(omega^2) - 1); the
    # internal log-normal variance is recovered by
    # omega^2 = log((CV/100)^2 + 1).
    #
    #   F      CV 51.2% -> log(0.512^2 + 1) = 0.23289
    #   CL_LF  CV 11.2% -> log(0.112^2 + 1) = 0.01246
    #   VC_LF  CV 119%  -> log(1.19^2  + 1) = 0.88213
    #   Q_LF   CV 23.9% -> log(0.239^2 + 1) = 0.05553
    #   CL_DLF CV 23.2% -> log(0.232^2 + 1) = 0.05241
    #   VC_DLF CV 90.5% -> log(0.905^2 + 1) = 0.59835
    #
    # The other structural parameters (ka, lag, Q_DLF, VP_LF, VP_DLF)
    # have dashes in the IIV columns of Table 2 -- no IIV is encoded
    # for them. -----
    etalfdepot ~ 0.23289   # Table 2: BSV on F      = 51.2% CV (RSE 14.3%; 95% CI 42.7-58.3) -- Box-Cox shape -0.394 NOT applied; see model description + Errata
    etalcl     ~ 0.01246   # Table 2: BSV on CL_LF  = 11.2% CV (RSE 49.5%; 95% CI 3.20-15.7)
    etalvc     ~ 0.88213   # Table 2: BSV on VC_LF  = 119%  CV (RSE 29.4%; 95% CI 78.1-178)
    etalq      ~ 0.05553   # Table 2: BSV on Q_LF   = 23.9% CV (RSE 45.2%; 95% CI 9.14-33.7)
    etalcl_dlf ~ 0.05241   # Table 2: BSV on CL_DLF = 23.2% CV (RSE 59.1%; 95% CI 6.93-36.5)
    etalvc_dlf ~ 0.59835   # Table 2: BSV on VC_DLF = 90.5% CV (RSE 93.9%; 95% CI 44.8-181)

    # ----- Residual error (Kloprogge 2015 Table 2 last 4 rows).
    # Table 2 footnote (a) declares sigma as the VARIANCE of the
    # unexplained residual variability on log-transformed
    # observations (Methods: 'Molar units of LF and DLF plasma
    # concentrations were transformed into their natural logarithms
    # and modeled simultaneously'). NONMEM additive-on-log-scale
    # error maps to proportional error in nlmixr2's linear-
    # concentration space; the proportional SD equals sqrt(variance).
    # Capillary-residual variances (sigma_capillary_LF = 0.0464,
    # sigma_capillary_DLF = 0.0326) and the matrix conversion
    # factors (LF 0.878, DLF 0.464) are documented in vignette
    # Errata but NOT encoded in this primary venous-only model file
    # (same simplification convention as Kloprogge_2013). -----
    propSd     <- sqrt(0.252)                                                        # Table 2: sigma_venous_LF  = 0.252 (RSE 54.0%; 95% CI 0.0208-0.550) -- variance on log scale; SD = sqrt(0.252) = 0.502
    label("LF venous-plasma proportional residual SD (SD on log scale, applied as proportional error in linear space)")
    propSd_dlf <- sqrt(0.335)                                                        # Table 2: sigma_venous_DLF = 0.335 (RSE 52.7%; 95% CI 0.0178-0.639) -- variance on log scale; SD = sqrt(0.335) = 0.579
    label("DLF venous-plasma proportional residual SD (SD on log scale, applied as proportional error in linear space)")
  })

  model({
    # Molecular weights for parent (lumefantrine) and active metabolite
    # (desbutyl-lumefantrine). Used to convert the parent molar
    # elimination flux into a DLF mass formation flux at the parent ->
    # metabolite step. The DLF molecule is the parent minus one
    # n-butyl group (CYP3A4-mediated N-debutylation).
    # PubChem CIDs 6437380 (lumefantrine: C30H32Cl3NO, 528.94 g/mol)
    # and 11281538 (desbutyl-lumefantrine: C26H24Cl3NO, 472.83 g/mol).
    # MW values are dimensionally g/mol = mg/mmol.
    MW_LF  <- 528.94
    MW_DLF <- 472.83

    # ----- Individual PK parameters. -----
    # LF absorption: gestational-age power effect on ka.
    ka      <- exp(lka)                  * (GA / 22.8)^e_ga_ka
    tlag    <- exp(ltlag)
    # LF disposition: body weight as power covariate on CL/Q/Vc/Vp
    # was tested and NOT retained in the paper's backward step (no
    # allometric scaling). Q/F carries the gestational-age linear-
    # deviation effect; the other LF parameters carry no covariate.
    cl      <- exp(lcl     + etalcl)
    vc      <- exp(lvc     + etalvc)
    q       <- exp(lq      + etalq)      * (1 + e_ga_q * (GA - 22.8))
    vp      <- exp(lvp)
    # DLF disposition: CL_DLF carries the admission-parasitaemia
    # log10-exponential effect; the other DLF parameters carry no
    # covariate.
    cl_dlf  <- exp(lcl_dlf + etalcl_dlf) * exp(e_para_cl_dlf * (log10(max(PARA, 1)) - log10(3260)))
    vc_dlf  <- exp(lvc_dlf + etalvc_dlf)
    q_dlf   <- exp(lq_dlf)
    vp_dlf  <- exp(lvp_dlf)

    # ----- Two-compartment disposition micro-constants for parent and
    # metabolite. -----
    kel     <- cl     / vc
    k12     <- q      / vc
    k21     <- q      / vp
    kel_dlf <- cl_dlf / vc_dlf
    k12_dlf <- q_dlf  / vc_dlf
    k21_dlf <- q_dlf  / vp_dlf

    # ----- ODE system. -----
    # LF: first-order absorption from depot (with absorption lag tlag)
    # into a 2-compartment disposition.
    d/dt(depot)           <- -ka * depot
    d/dt(central)         <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1)     <-                                k12 * central - k21 * peripheral1

    # DLF: formed mole-for-mole from LF central elimination (1:1
    # stoichiometric N-debutylation, implicit fm = 1), with the
    # mass-units flux scaled by MW_DLF / MW_LF (= 0.894) so that the
    # accumulated central_dlf is in DLF mass and Cc_dlf (= central_dlf
    # / vc_dlf) is in true DLF mass concentration. DLF disposition is
    # 2-compartment linear with apparent CL_DLF / Vc_DLF / Q_DLF /
    # Vp_DLF estimated by the paper's molar-units fit.
    d/dt(central_dlf)     <-  kel * central * (MW_DLF / MW_LF) - kel_dlf * central_dlf - k12_dlf * central_dlf + k21_dlf * peripheral1_dlf
    d/dt(peripheral1_dlf) <-                                                              k12_dlf * central_dlf - k21_dlf * peripheral1_dlf

    # ----- Bioavailability and absorption lag time on the depot. -----
    f(depot)    <- exp(lfdepot + etalfdepot)
    alag(depot) <- tlag

    # ----- Observations. Dose units are mg (1 tablet = 120 mg LF),
    # volumes are L, so central / vc returns concentration in mg/L =
    # ug/mL. Paper Tables 2-3 report day-7 and Cmax in ng/mL; the
    # validation vignette multiplies Cc and Cc_dlf by 1000 for ng/mL
    # comparison.
    Cc      <- central     / vc
    Cc_dlf  <- central_dlf / vc_dlf

    Cc      ~ prop(propSd)
    Cc_dlf  ~ prop(propSd_dlf)
  })
}
