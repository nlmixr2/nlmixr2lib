Urien_2004_cisplatin <- function() {
  description <- paste(
    "Integrated two-compartment population PK model for ultrafilterable",
    "(unbound) plasma platinum coupled to a metabolite compartment",
    "representing irreversibly protein-bound plasma platinum (Urien 2004",
    "BJCP). Fitted simultaneously to 396 unbound and 477 total plasma",
    "platinum concentration-time observations from 43 adult cancer",
    "patients receiving 30-min cisplatin infusions. Unbound clearance",
    "depends on body surface area and Cockcroft-Gault creatinine",
    "clearance; the unbound central volume depends on BSA. The bound",
    "formation parameter fm/Vm depends on dose per m^2 and total serum",
    "protein; its BSA and creatinine-clearance exponents are fixed at",
    "the negatives of the unbound CL exponents so the formation flux",
    "(fm/Vm)*CL*Cc is net BSA- and CLCr-neutral. The apparent metabolite",
    "volume Vm is not separately identifiable; the composite parameters",
    "fm/Vm (1/L) and CLm0/Vm (1/h) absorb Vm and the bound state is",
    "carried in concentration units."
  )
  reference <- "Urien S, Lokiec F. Population pharmacokinetics of total and unbound plasma cisplatin in adult patients. Br J Clin Pharmacol. 2004;57(6):756-763. doi:10.1111/j.1365-2125.2004.02082.x"
  vignette  <- "Urien_2004_cisplatin"

  # The protein-bound platinum state is a paper-mechanistic compartment
  # (irreversible plasma-protein-binding pathway from compartment 1, per
  # Urien 2004 Figure 1). Vm is structurally not identifiable; the bound
  # state is carried in concentration units (ug/L) with the identifiable
  # composite parameters fm/Vm and CLm0/Vm absorbing the volume.
  paper_specific_compartments <- c("bound")

  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling centred on the cohort median 1.62 m^2 (Urien 2004 Table 1). Exponent +1.60 on the unbound central volume Vc (Table 4) and +0.83 on unbound clearance CL (Table 4). The same BSA-on-CL effect also enters the bound-formation parameter fm/Vm as a FIXED exponent -0.83 (Table 5, opposite sign) so the bound-formation flux (fm/Vm)*CL*Cc is net BSA-neutral. Urien 2004 does not state which BSA formula was used; record 'unspecified' downstream.",
      source_name        = "BSA"
    ),
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling centred on the cohort median 81 mL/min (Urien 2004 Table 1). Exponent +0.36 on unbound clearance CL (Table 4). The same CLCr-on-CL effect enters fm/Vm as a FIXED exponent -0.36 (Table 5, opposite sign) so the bound-formation flux is net CLCr-neutral. The paper uses the raw Cockcroft-Gault CrCl in mL/min, NOT BSA-normalized to mL/min/1.73 m^2; document this in downstream simulations so body size is not double-corrected.",
      source_name        = "CLCr"
    ),
    TPRO = list(
      description        = "Total serum protein concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling centred on the cohort median 70 g/L (Urien 2004 Table 1). Exponent +1.33 on the bound-formation parameter fm/Vm (Table 5). Mechanism: higher serum protein concentration shifts the binding equilibrium toward protein-bound platinum.",
      source_name        = "PROT"
    ),
    DOSE = list(
      description        = "Cisplatin dose per infusion (mg) supplied as a per-record data column",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-record administered cisplatin dose in mg. Used inside model() to compute dose per BSA = DOSE / BSA (mg/m^2). Power scaling on fm/Vm centred on dose per BSA = 25 mg/m^2 (cohort median, Urien 2004 Table 1), exponent -0.48 (Table 5). The negative dose-per-area exponent reflects saturable plasma-protein binding at higher doses. For per-record simulation, LOCF the administered dose from the most recent infusion across subsequent observation records so dose_per_bsa stays defined between infusions.",
      source_name        = "Dose"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 43L,
    n_studies      = 2L,
    n_observations = "873 plasma platinum measurements (396 unbound + 477 total) across 146 30-min infusions",
    age_range      = "21-76 years",
    age_median     = "58 years",
    weight_range   = "40-102 kg",
    weight_median  = "64 kg",
    height_range   = "150-181 cm",
    height_median  = "168 cm",
    bsa_range      = "1.38-2.10 m^2",
    bsa_median     = "1.62 m^2",
    crcl_range     = "44-155 mL/min (Cockcroft-Gault)",
    crcl_median    = "81 mL/min",
    tpro_range     = "47-80 g/L",
    tpro_median    = "70 g/L",
    scr_range      = "43-120 umol/L",
    scr_median     = "76 umol/L",
    sex_female_pct = 41.9,
    disease_state  = "Adult patients with metastatic cancer receiving second- or third-line chemotherapy in two phase I studies (cisplatin combined with either irofulven 0.4 mg/kg or 5-fluorouracil 1 g/m^2/day continuous 120-h infusion).",
    dose_range     = "15-80 mg total per 30-min IV infusion (median 34.4 mg; median 25 mg/m^2); 5 consecutive daily infusions or twice-monthly schedules.",
    regions        = "France (Centre Rene Huguenin, Saint-Cloud).",
    notes          = "Demographics from Urien 2004 Table 1 (male/female 25/18). Concomitant chemotherapy (irofulven in 18 patients, 5-FU in 25 patients) was tested as a covariate and not retained. Infusion duration typically 30 min but varied 0.25-1 h across 146 infusions; one to five consecutive daily infusions per evaluation. Modelling software: MP2 nonlinear mixed-effect program (Urien 2004 reference 3)."
  )

  ini({
    # Unbound platinum, structural parameters - Urien 2004 Table 4 final-model
    # estimates. Typical values hold at the cohort median covariates: BSA = 1.62
    # m^2, CLCr = 81 mL/min, PROT = 70 g/L, dose per BSA = 25 mg/m^2. The
    # source's integrated three-step fit estimated these in Step 2 (unbound-only
    # 2-compartment fit) and then held them fixed in Step 3 (joint unbound +
    # bound fit). They are recorded here as the published Step-2 point estimates,
    # not wrapped in fixed(), so a downstream re-fitting workflow has the same
    # provenance signal as the rest of the popPK library.
    lvc <- log(23.4); label("Unbound central volume V1 (L)")                                              # Urien 2004 Table 4: V1 = 23.4 L
    lcl <- log(35.5); label("Unbound systemic clearance CL (L/h)")                                        # Urien 2004 Table 4: CL = 35.5 L/h
    lq  <- log(8.64); label("Unbound inter-compartmental clearance Q (L/h)")                              # Urien 2004 Table 4: Q  = 8.64 L/h
    lvp <- log(12.0); label("Unbound peripheral volume V2 (L)")                                           # Urien 2004 Table 4: V2 = 12.0 L

    # Unbound-layer covariate exponents - Urien 2004 Table 4.
    e_bsa_vc  <- 1.60; label("Power exponent of BSA/1.62 on Vc (unitless)")                               # Urien 2004 Table 4: theta_BSA on V1 = +1.60
    e_bsa_cl  <- 0.83; label("Power exponent of BSA/1.62 on CL (unitless)")                               # Urien 2004 Table 4: theta_BSA on CL = +0.83
    e_crcl_cl <- 0.36; label("Power exponent of CRCL/81 on CL (unitless)")                                # Urien 2004 Table 4: theta_CLCr on CL = +0.36

    # Bound-layer composite parameters - Urien 2004 Table 5. The metabolite
    # volume Vm is structurally not identifiable, so fm/Vm and CLm0/Vm
    # absorb Vm and parameterise the bound formation and elimination
    # respectively. The bound state is therefore carried in concentration
    # units (ug/L) directly.
    lfmVm   <- log(0.017); label("Bound-formation composite parameter fm/Vm (1/L)")                       # Urien 2004 Table 5: fm/Vm = 0.017 1/L
    lclm0Vm <- log(0.014); label("Bound-elimination rate constant CLm0/Vm (1/h)")                         # Urien 2004 Table 5: CLm0/Vm = 0.014 1/h

    # Bound-layer covariate effects - Urien 2004 Table 5. The BSA and CLCr
    # exponents on fm/Vm are FIXED at the negatives of their CL counterparts
    # (Table 5: "(fixed)") so the bound-formation flux (fm/Vm)*CL*Cc is net
    # BSA- and CLCr-neutral; only the dose-per-area and total-protein
    # exponents are estimated.
    e_dose_fmVm <- -0.48;        label("Power exponent of (DOSE/BSA)/25 on fm/Vm (unitless)")             # Urien 2004 Table 5: theta_DOSEm-2 on fm/Vm = -0.48
    e_tpro_fmVm <-  1.33;        label("Power exponent of TPRO/70 on fm/Vm (unitless)")                   # Urien 2004 Table 5: theta_PROT  on fm/Vm = +1.33
    e_bsa_fmVm  <- fixed(-0.83); label("Fixed power exponent of BSA/1.62 on fm/Vm (unitless)")            # Urien 2004 Table 5: theta_BSA   on fm/Vm = -0.83 (fixed)
    e_crcl_fmVm <- fixed(-0.36); label("Fixed power exponent of CRCL/81 on fm/Vm (unitless)")             # Urien 2004 Table 5: theta_CLCr on fm/Vm = -0.36 (fixed)

    # Inter-individual variability - log-normal IIV (omega^2 = log(1 + CV^2)).
    # Urien 2004 reports ISV as percent CV in Tables 4 and 5. Q and V2 ISVs
    # were fixed to zero in the source (Results paragraph beneath Table 4:
    # "The data did not allow reliable estimations of intersubject variabilities
    # for Q and V2 and fixing these parameters to zero did not increase the
    # OFV"); no eta is declared for lq or lvp.
    etalvc      ~ 0.05334  # log(1 + 0.234^2); ISV V1     = 23.4% CV (Urien 2004 Table 4)
    etalcl      ~ 0.02199  # log(1 + 0.149^2); ISV CL     = 14.9% CV (Urien 2004 Table 4)
    etalfmVm   ~ 0.04498  # log(1 + 0.214^2); ISV fm/Vm   = 21.4% CV (Urien 2004 Table 5)
    etalclm0Vm ~ 0.09617  # log(1 + 0.318^2); ISV CLm/Vm = 31.8% CV (Urien 2004 Table 5)

    # Residual error - additive on linear scale, separate values for unbound
    # and total plasma platinum observations. Urien 2004 reports residuals in
    # ug/L (= ng/mL); this file uses ug/mL = mg/L throughout, so 116 ug/L =
    # 0.116 ug/mL and 103 ug/L = 0.103 ug/mL. The integrated-model unbound
    # residual carries over from Step 2 (Table 4: 116 ug/L) and the total-
    # platinum residual is estimated in Step 3 (Table 5: 103 ug/L).
    addSd          <- 0.116; label("Unbound platinum additive residual SD (ug/mL)")                       # Urien 2004 Table 4: Res. Error = 116 ug/L = 0.116 ug/mL
    addSd_Cc_total <- 0.103; label("Total platinum additive residual SD (ug/mL)")                         # Urien 2004 Table 5: Res. Error = 103 ug/L = 0.103 ug/mL
  })

  model({
    # Derived covariate term: dose per body surface area (mg/m^2).
    dose_per_bsa <- DOSE / BSA

    # Individual unbound-layer PK parameters with the Urien 2004 Table 4
    # covariate equations (power form, centred on cohort medians).
    vc <- exp(lvc + etalvc) * (BSA / 1.62)^e_bsa_vc
    cl <- exp(lcl + etalcl) * (BSA / 1.62)^e_bsa_cl * (CRCL / 81)^e_crcl_cl
    q  <- exp(lq)
    vp <- exp(lvp)

    # Individual bound-layer composite parameters - Urien 2004 Table 5.
    fmVm <- exp(lfmVm + etalfmVm) *
             (dose_per_bsa / 25)^e_dose_fmVm *
             (TPRO / 70)^e_tpro_fmVm *
             (BSA / 1.62)^e_bsa_fmVm *
             (CRCL / 81)^e_crcl_fmVm
    clm0Vm <- exp(lclm0Vm + etalclm0Vm)

    # Unbound concentrations in the central and peripheral compartments.
    Cc <- central / vc
    C2 <- peripheral1 / vp

    # ODE system - Urien 2004 Figure 1 scheme. Compartments 1 and 2 carry
    # unbound platinum (amounts in mg); the bound state carries protein-bound
    # platinum concentration directly (ug/L), with Vm folded into fm/Vm and
    # CLm0/Vm so the apparent volume drops out.
    d/dt(central)     <- -cl * Cc - q * (Cc - C2)
    d/dt(peripheral1) <-  q * (Cc - C2)
    d/dt(bound)       <-  fmVm * cl * Cc - clm0Vm * bound

    # Observations. Total plasma platinum = unbound + protein-bound, both
    # measured in the same plasma sample. Dose in mg, vc in L, so Cc =
    # central/vc is in mg/L = ug/mL. The bound state was derived in matching
    # units (the (fm/Vm)*cl*Cc formation term is in ug/(mL*h)).
    Cc_total <- Cc + bound

    Cc       ~ add(addSd)
    Cc_total ~ add(addSd_Cc_total)
  })
}
