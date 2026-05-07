# Joint parent-metabolite population PK model of oral artesunate and its
# active metabolite dihydroartemisinin (DHA) in pregnant and non-pregnant
# women with uncomplicated Plasmodium falciparum malaria, extracted from
# the DDMORE Foundation Model Repository entry DDMODEL00000297 (the
# Birgersson 2019 Wellcome Open Research paper).

Birgersson_2019_artesunate <- function() {
  description <- "Joint parent-metabolite population PK model for oral artesunate and dihydroartemisinin in pregnant and non-pregnant women with uncomplicated Plasmodium falciparum malaria, with a 3-compartment transit-absorption chain into a 1-compartment artesunate disposition model and complete in-vivo conversion to a 1-compartment DHA disposition model. Allometric body-weight scaling on CL (exponent 0.75) and V (exponent 1) is applied to both parent and metabolite. Covariate effects: pregnancy on DHA elimination clearance, and admission alanine aminotransferase and log-asexual-parasite-count on relative bioavailability of artesunate."
  reference <- paste(
    "Birgersson S, Valea I, Tinto H, Traore-Coulibaly M, Toe LC,",
    "Hoglund RM, Van Geertruyden JP, Ward SA, D'Alessandro U,",
    "Abelo A, Tarning J (2019).",
    "Population pharmacokinetics of artesunate and dihydroartemisinin in",
    "pregnant and non-pregnant women with uncomplicated",
    "Plasmodium falciparum malaria in Burkina Faso: an open label trial.",
    "Wellcome Open Res 4:45.",
    "doi:10.12688/wellcomeopenres.14849.2.",
    "DDMORE Foundation Model Repository: DDMODEL00000297.",
    sep = " "
  )
  vignette <- "Birgersson_2019_artesunate"
  ddmore_id <- "DDMODEL00000297"
  replicate_of <- NULL
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Allometric scaling with reference weight 52 kg (the population median in the Burkina Faso cohort) is applied to artesunate CL (exponent 0.75) and Vc (exponent 1.0) and to dihydroartemisinin CL (exponent 0.75) and Vc (exponent 1.0). Reference weight 52 kg is hard-coded in the source NONMEM .mod ($PK block: TVCLP = THETA(1)*((WT/52)**0.75) etc.).",
      source_name        = "WT"
    ),
    PREG = list(
      description        = "Pregnancy status",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "1 = pregnant, 0 = non-pregnant. The Birgersson 2019 cohort consists of 24 pregnant women in their second or third trimester paired with 24 non-pregnant women, all with uncomplicated Plasmodium falciparum malaria. The source NONMEM .mod uses pregnant women as the within-paper reference (CLMPREG = 1 when PREG = 1; CLMPREG = 1 + THETA(7) when PREG = 0). Verbatim source values are preserved by applying the effect parameter via (1 - PREG) in the model code: structural TVCLM = 190 L/h corresponds to pregnant women, and non-pregnant women have CLM scaled by (1 + e_preg_cl_dha) = 0.786.",
      source_name        = "PREG"
    ),
    ALT = list(
      description        = "Alanine aminotransferase activity at admission (baseline, time-fixed per subject)",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centered at the population median 20.75 U/L (hard-coded in the source NONMEM .mod). Linear-deviation effect on artesunate relative bioavailability F1: F1ALT = 1 + e_alt_fdepot * (ALT - 20.75) with e_alt_fdepot = +0.0215 per U/L. Birgersson 2019 includes ALT as a marker of acute liver-status change in symptomatic malaria.",
      source_name        = "ALT"
    ),
    LNPC = list(
      description        = "Natural logarithm of the asexual Plasmodium falciparum parasite count at admission (parasites per microlitre of blood)",
      units              = "log(parasites/uL)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at admission. Centered at the population median 5.88 log(parasites/uL) (hard-coded in the source NONMEM .mod). Linear-deviation effect on artesunate relative bioavailability F1: F1LNPC = 1 + e_lnpc_fdepot * (LNPC - 5.88) with e_lnpc_fdepot = +0.138 per log-unit increase in parasitaemia. The companion source column PARA carries the raw parasite count per microlitre; LNPC = log(PARA) is the active model covariate.",
      source_name        = "LNPC"
    )
  )

  population <- list(
    n_subjects     = 48,
    n_studies      = 1,
    n_pregnant     = 24,
    n_nonpregnant  = 24,
    age_range      = "Adult women",
    weight_range   = "Adult women (population median 52 kg, the WT allometric reference)",
    sex_female_pct = 100,
    disease_state  = "Uncomplicated Plasmodium falciparum malaria; pregnant women in the second or third trimester paired with matched non-pregnant women",
    dose_range     = "Oral artesunate (fixed-dose combination with mefloquine) once daily for three days at the standard adult dose. The DDMORE-shipped Simulated_run1.csv encodes a 520264 nmol single oral dose, equivalent to ~200 mg artesunate using the 384.42 g/mol molar mass.",
    regions        = "Burkina Faso (Bobo-Dioulasso and Nanoro)",
    trial_registration = "ClinicalTrials.gov NCT00701961",
    notes          = "Demographics summarized from the publication abstract (PMID 32025570) and the DDMORE bundle metadata (DDMODEL00000297). The publication PDF was not on disk during this extraction; the DDMORE bundle ships only Executable_run1.mod (a MAXEVAL=0 posthoc run with the published final estimates fixed as $THETA / $OMEGA / $SIGMA inputs), Output_simulated_run1.lst (a one-subject re-evaluation), Simulated_run1.csv, and DDMODEL00000297.rdf. No Output_real_*.lst or Model_Accomodations.text is included; final estimates were taken from the .mod file, which under MAXEVAL=0 functions as the authoritative source for the published estimates."
  )

  ini({
    # Structural-parameter values come directly from the Executable_run1.mod
    # $THETA / $OMEGA / $SIGMA blocks. Because the DDMORE bundle ships a
    # MAXEVAL=0 (POSTHOC, METHOD=1 LAPLACIAN INTER) run with no estimation
    # step, the .mod values are the published Birgersson 2019 final
    # estimates carried verbatim. The Output_simulated_run1.lst FINAL
    # PARAMETER ESTIMATE block is identical (TH 1 = 3.57E+03, etc.),
    # confirming that no further fitting was performed.

    lcl   <- log(3570)   ; label("Apparent artesunate elimination clearance, CLP/F at WT = 52 kg (L/h); equivalent to the artesunate-to-DHA conversion clearance under the source paper's assumption of complete in-vivo conversion")  # Executable_run1.mod $THETA TH 1 (CLP) = 3570 L/h
    lvc   <- log(1700)   ; label("Apparent artesunate central volume of distribution, V2/F at WT = 52 kg (L)")  # Executable_run1.mod $THETA TH 2 (V2) = 1700 L
    lcl_dha <- log(190)  ; label("Apparent dihydroartemisinin elimination clearance, CLM/F at WT = 52 kg in pregnant women (L/h); non-pregnant women have CLM reduced via the e_preg_cl_dha covariate effect")  # Executable_run1.mod $THETA TH 3 (CLM) = 190 L/h
    lvc_dha <- log(267)  ; label("Apparent dihydroartemisinin central volume of distribution, V3/F at WT = 52 kg (L)")  # Executable_run1.mod $THETA TH 4 (V3) = 267 L
    lmtt  <- log(0.832)  ; label("Mean transit time of the 3-compartment transit-absorption chain, MTT (h)")  # Executable_run1.mod $THETA TH 5 (MTT) = 0.832 h
    lfdepot <- fixed(log(1))   ; label("Reference relative bioavailability of artesunate, F1 (unitless); fixed at 1 (the source paper estimates F1 only via covariate effects of ALT and parasitaemia, with no absolute reference)")  # Executable_run1.mod $THETA TH 6 (F1) = 1 FIX

    # Covariate effects
    e_preg_cl_dha <- -0.214 ; label("Pregnancy-status effect on DHA elimination clearance, applied as (1 + e_preg_cl_dha * (1 - PREG)) so that pregnant women (PREG = 1) match the structural TVCLM = 190 L/h and non-pregnant women (PREG = 0) have CLM scaled by 0.786 (~21% lower CLM than pregnant women)")  # Executable_run1.mod $THETA TH 7 (CLMPREG1) = -0.214
    e_alt_fdepot  <-  0.0215 ; label("Linear-deviation effect of admission ALT (centered at 20.75 U/L) on artesunate relative bioavailability, per U/L")  # Executable_run1.mod $THETA TH 8 (F1ALT1) = 0.0215
    e_lnpc_fdepot <-  0.138  ; label("Linear-deviation effect of admission log-asexual-parasite-count (centered at 5.88 log(parasites/uL)) on artesunate relative bioavailability, per log-unit")  # Executable_run1.mod $THETA TH 9 (F1LNPC1) = 0.138

    # IIV (log-normal). The source $OMEGA carries six diagonal entries; the
    # entries for V2 and V3 are 0 FIXED (no IIV on artesunate or DHA volume),
    # so only four etas are carried into the model file.
    etalcl     ~ 0.0672    # Executable_run1.mod $OMEGA(1,1) (1.CL) = 0.0672 (variance on log-scale)
    etalcl_dha ~ 0.00809   # Executable_run1.mod $OMEGA(3,3) (3.CLM_) = 0.00809
    etalmtt    ~ 0.32      # Executable_run1.mod $OMEGA(5,5) (5.MTT_) = 0.32
    etalfdepot ~ 0.0887    # Executable_run1.mod $OMEGA(6,6) (8.F1) = 0.0887

    # Residual error. The source $ERROR block evaluates the M3 BQL
    # likelihood on the log-transformed observation: when not BQL,
    # Y = log(IPRED) + ERR(i), i.e., additive on the log scale. By the
    # convention documented in references/naming-conventions.md, NONMEM
    # additive-on-log-scale residual error maps to proportional residual
    # error in nlmixr2's linear space. The propSd value is the SD on log
    # scale, sqrt(variance). The M3 likelihood for BQL data is not
    # reproduced here -- the nlmixr2lib model is intended for
    # forward-simulation / typical-value use, where the BQL handling does
    # not affect the trajectory.
    propSd     <- sqrt(0.892) ; label("Proportional residual SD for artesunate plasma concentration (SD on log scale)")  # Executable_run1.mod $SIGMA(1,1) (RUV_ARS) = 0.892 (variance); SD = sqrt(0.892) ~= 0.945
    propSd_dha <- sqrt(0.66)  ; label("Proportional residual SD for dihydroartemisinin plasma concentration (SD on log scale)")  # Executable_run1.mod $SIGMA(2,2) (RUV_DHA) = 0.66 (variance); SD = sqrt(0.66) ~= 0.812
  })

  model({
    # Individual PK parameters. Allometric scaling on CL with exponent 0.75
    # and on V with exponent 1, applied per the source NONMEM .mod $PK
    # block (TVCLP = THETA(1)*((WT/52)**0.75); TVV2 = THETA(2)*((WT/52)**1)
    # etc.). Reference weight 52 kg is the cohort median.
    cl     <- exp(lcl     + etalcl)     * (WT / 52)^0.75
    vc     <- exp(lvc)                  * (WT / 52)
    cl_dha <- exp(lcl_dha + etalcl_dha) * (WT / 52)^0.75 *
              (1 + e_preg_cl_dha * (1 - PREG))
    vc_dha <- exp(lvc_dha)              * (WT / 52)

    # Mean transit time and chain rate. With NN = 3 transit compartments,
    # KTR = (NN + 1) / MTT = 4 / MTT (per the source NONMEM .mod $PK
    # block: NN = 3; KTR = (NN+1)/MT). The same KTR is used for K14, K45,
    # K56, K62, so the chain depot -> transit1 -> transit2 -> transit3 ->
    # central propagates with a single rate constant.
    mtt <- exp(lmtt + etalmtt)
    ktr <- 4 / mtt

    # Conversion / elimination micro-constants. Under the source paper's
    # assumption of complete in-vivo conversion of artesunate to DHA, all
    # artesunate clearance is metabolic conversion, so the rate from
    # central -> central_dha is cl/vc. DHA is then eliminated linearly at
    # cl_dha/vc_dha.
    kel_art <- cl     / vc
    kel_dha <- cl_dha / vc_dha

    # ODE system. Source NONMEM compartment indexing: 1 depot, 2 artesunate
    # central, 3 DHA central, 4-6 transit. Translated to canonical names
    # depot, central (artesunate), central_dha (DHA), transit1, transit2,
    # transit3 with the same connectivity (K14, K45, K56, K62, K23, K30 in
    # the source .mod).
    d/dt(depot)        <- -ktr * depot
    d/dt(transit1)     <-  ktr * depot    - ktr * transit1
    d/dt(transit2)     <-  ktr * transit1 - ktr * transit2
    d/dt(transit3)     <-  ktr * transit2 - ktr * transit3
    d/dt(central)      <-  ktr * transit3 - kel_art * central
    d/dt(central_dha)  <-  kel_art * central - kel_dha * central_dha

    # Bioavailability on the depot compartment. The source NONMEM .mod
    # parameterizes F1 = TVF1 * F1COV * EXP(ETA(6)) with TVF1 = 1 FIX and
    # F1COV = (1 + theta_alt * (ALT - 20.75)) * (1 + theta_lnpc * (LNPC -
    # 5.88)). Translated as a multiplicative chain on f(depot).
    f(depot) <- exp(lfdepot + etalfdepot) *
                (1 + e_alt_fdepot  * (ALT  - 20.75)) *
                (1 + e_lnpc_fdepot * (LNPC -  5.88))

    # Plasma concentrations in nmol/L. Source data file comment:
    # "DOSE(NMOL)CP(NMOL/L)" -- doses in nmol, concentrations in nmol/L.
    # Because both species are tracked on a molar basis under complete 1:1
    # conversion, no MW conversion is needed for the central -> central_dha
    # transfer.
    Cc     <- central     / vc
    Cc_dha <- central_dha / vc_dha

    # Residual error. NONMEM "additive on log scale" maps to proportional
    # error in linear space (see references/naming-conventions.md).
    Cc     ~ prop(propSd)
    Cc_dha ~ prop(propSd_dha)
  })
}
