Allegaert_2015_paracetamol <- function() {
  description <- "Eight-compartment population PK model for IV propacetamol/paracetamol (APAP) and its glucuronide and sulphate metabolites in young women (Allegaert 2015), distributed in the DDMORE Foundation Model Repository as DDMODEL00000267. The structural model carries a three-compartment plasma disposition for parent APAP (central + two peripherals), two plasma metabolite compartments (APAP-glucuronide and APAP-sulphate, each with a metabolite-specific volume V_meta = 0.18 * V_central), and three cumulative-urine compartments (urine APAP, urine APAP-glucuronide, urine APAP-sulphate). Pregnancy state, time post partum, term-vs-preterm birth, oral-contraceptive use, and time-varying urine flow rate enter as covariates on the parent and metabolite-formation clearances; an OCC-conditional residual-error model gives non-pregnant volunteers on birth control a combined proportional + additive plasma error while every other occasion uses a proportional-only plasma error."
  reference <- paste(
    "Allegaert K, van der Marel CD, Debeer A, Pluim MAL, Van Lingen RA, Vanhole C,",
    "Tibboel D, Devlieger H. (2015).",
    "Pharmacokinetics of single-dose intravenous propacetamol in young women.",
    "BMC Anesthesiol 15:151.",
    "doi:10.1186/s12871-015-0144-3.",
    "DDMORE Foundation Model Repository: DDMODEL00000267.",
    sep = " "
  )
  vignette <- "Allegaert_2015_paracetamol"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")
  ddmore_id    <- "DDMODEL00000267"
  replicate_of <- NULL

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear allometric exponent of 1.0 on the inter-compartmental clearance Q1 (THETA(6) * (BW/70)^1 in the .mod). Reference weight 70 kg. Source column is `BW` (body weight, kg) -- same orientation, no transformation; rename to `WT` before passing the dataset to rxSolve.",
      source_name        = "BW"
    ),
    OCC = list(
      description        = "Five-level integer state column encoding pregnancy / postpartum / non-pregnant-volunteer status with a within-volunteer split for oral contraceptive use.",
      units              = "(categorical)",
      type               = "categorical",
      reference_category = "3 (one year post delivery; baseline value of CL_gluc and V_central)",
      notes              = "Allegaert 2015 NONMEM column with five mutually-exclusive levels: 1 = pregnant, 2 = 2 weeks postpartum, 3 = 1 year postpartum, 4 = non-pregnant volunteer not on birth control, 5 = non-pregnant volunteer on birth control. Decomposed inside `model()` into binary indicators `oc1` .. `oc5` that drive the per-occasion typical-value covariate effects on V_central (OCC = 1), CL_glucuronide (OCC = 1 and OCC = 2), Q2 (OCC = 2), and the residual-error switch on the plasma APAP observation (OCC = 5).",
      source_name        = "OCC"
    ),
    TERM_BIRTH = list(
      description        = "Term-vs-preterm birth indicator (1 = term birth, >= 37 weeks gestation; 0 = preterm).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (preterm)",
      notes              = "Selects between two typical-value sulphate-formation clearances via `CL_sulf = TERM_BIRTH * cl_sulf_term + (1 - TERM_BIRTH) * cl_sulf_preterm` (.mod $PK line `CLS=(TERM*THETA(3)+(1-TERM)*THETA(10))*EXP(ETA(3))`). Both arms are estimated parameters; neither category is the multiplicative reference. The .mod's $THETA comments label TH3 as 'CL formation APAP-S for TERM=0' and TH10 as 'CL formation APAP-S for TERM=1', but the $PK code applies TH3 when TERM=1 and TH10 when TERM=0 -- the code-level usage is taken as authoritative because `(TERM*TH3 + (1-TERM)*TH10)` evaluates that way.",
      source_name        = "TERM"
    ),
    BIRTHCONTROL = list(
      description        = "Oral hormonal contraceptive use indicator (1 = on oral contraceptive, 0 = not).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no oral contraceptive)",
      notes              = "Multiplicative scalar on the OCC-adjusted glucuronide-formation clearance (.mod $PK lines 39-43): `CL_gluc <- THETA(12) * CL_gluc_OCC` when BC == 1. Independent of OCC: a postpartum woman on contraceptives has both an OCC-2 typical-value adjustment AND the BIRTHCONTROL = 1 multiplier applied. The estrogen component of combined oral contraceptives is plausibly the driver via UGT2B7 induction.",
      source_name        = "BC"
    ),
    URINE_FLOW = list(
      description        = "Instantaneous urine flow rate over the urine-collection interval.",
      units              = "mL/h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear effect on renal clearance of unchanged paracetamol with centering at 100 mL/h: `CL_renal_apap = THETA(4) + THETA(16) * (URINE_FLOW - 100)` (.mod $PK lines 45-47), gated by `URINE_FLOW > 0` (urine-flow == 0 sets the slope contribution to zero, treated as 'no urine produced during the interval').",
      source_name        = "UF"
    )
  )

  population <- list(
    n_subjects     = 69L,
    n_studies      = 1L,
    age_range      = NA_character_,
    age_median     = NA_character_,
    weight_range   = NA_character_,
    weight_median  = NA_character_,
    sex_female_pct = 100,
    disease_state  = "Young women across five reproductive states (DDMODEL00000267 .mod $INPUT comment for OCC): pregnant (OCC = 1), 2 weeks postpartum (OCC = 2), 1 year postpartum (OCC = 3), non-pregnant volunteer not on birth control (OCC = 4), non-pregnant volunteer on birth control (OCC = 5). The intent of the analysis is to quantify how pregnancy, time post partum, and oral-contraceptive use modify paracetamol disposition in this population.",
    dose_range     = "Single-dose intravenous propacetamol (a paracetamol pro-drug; 1 g propacetamol is approximately equivalent to 0.5 g paracetamol on a molar basis), administered as a short infusion (the bundle's Simulated_APAP_YoungWomen.csv ships 1 g and 2 g propacetamol-equivalent doses delivered over ~20-30 min). Multiple-dose regimens are present in the simulated dataset for some subjects.",
    regions        = NA_character_,
    notes          = "n_subjects (69) and n_obs (1118) come from the Output_real_OriginalModelCode.lst run summary lines 'TOT. NO. OF INDIVIDUALS:' and 'TOT. NO. OF OBS RECS:' respectively. The Allegaert 2015 BMC Anesthesiol publication is not on disk in this worktree, so finer-grained demographics (age range, weight range, region, baseline-renal/hepatic profile) are recorded as NA. The .mod's covariate columns (OCC, BC -> BIRTHCONTROL, TERM -> TERM_BIRTH, UF -> URINE_FLOW, BW -> WT) reveal the variability axes the model resolves but do not by themselves constrain the underlying demographic distributions. See the validation vignette's Errata section for the full caveat list."
  )

  ini({
    # Final parameter estimates from
    # /home/bill/github/mab_human_consensus/literature/from_people/ddmore/ddmore_scraping/267/Output_real_OriginalModelCode.lst
    # FINAL PARAMETER ESTIMATE block (lst lines 519-567), captured after
    # `MINIMIZATION SUCCESSFUL` (lst line 452, OBJV 5286.743). NONMEM THETAs are
    # log-back-transformed values, so each `lX <- log(value)` wraps the lst
    # value to keep the internal scale log. NONMEM OMEGA / SIGMA entries are
    # variances on the relevant internal scale (log for log-normal IIV; linear
    # for proportional / additive residual error). nlmixr2 `prop()` / `add()`
    # SDs are sqrt of the NONMEM variances; the conversion is documented per
    # parameter.

    # Structural plasma-disposition parameters for parent APAP --------------
    lvc        <- log(18.5)   ; label("Central plasma volume of paracetamol V1 at OCC >= 2 (L)")               # TH8 FINAL = 1.85E+01
    lvp        <- log(19.7)   ; label("Plasma peripheral 1 volume V2 (L)")                                     # TH5 FINAL = 1.97E+01
    lvp2       <- log(23.9)   ; label("Plasma peripheral 2 volume V8 (L)")                                     # TH13 FINAL = 2.39E+01
    lq         <- log(1.29)   ; label("Inter-compartmental clearance Q1 (central <-> peripheral 1) at WT = 70 kg (L/h)") # TH6 FINAL = 1.29E+00
    lq2        <- log(61.1)   ; label("Inter-compartmental clearance Q2 (central <-> peripheral 2) at OCC != 2 (L/h)")   # TH14 FINAL = 6.11E+01

    # Metabolite-formation and renal-elimination clearances -----------------
    lcl_gluc      <- log(7.33)   ; label("CL of glucuronide formation at OCC >= 3 (L/h)")    # TH9 FINAL = 7.33E+00
    lcl_sulf      <- log(3.86)   ; label("CL of sulphate formation for preterm-birth subjects (TERM_BIRTH = 0) (L/h)") # TH10 FINAL = 3.86E+00
    lcl_renal     <- log(0.925)  ; label("Renal CL of unchanged paracetamol at urine flow 100 mL/h (L/h)")    # TH4 FINAL = 9.25E-01
    lf_meta_renal <- log(4.62)   ; label("k35 / k17: ratio of plasma-metabolite-to-urine elimination rate to plasma-APAP-to-urine elimination rate (unitless)") # TH7 FINAL = 4.62E+00 -- used as K35 = THETA(7) * K17, applied identically to glucuronide (k35) and sulphate (k46)

    # Covariate effects -----------------------------------------------------
    e_oc1_vc          <- 1.86   ; label("Multiplicative scalar on V_central for pregnancy (OCC = 1)") # TH1 FINAL = 1.86E+00
    e_oc1_cl_gluc     <- 2.03   ; label("Multiplicative scalar on glucuronide-formation CL for pregnancy (OCC = 1)") # TH2 FINAL = 2.03E+00
    e_oc2_cl_gluc     <- 0.547  ; label("Multiplicative scalar on glucuronide-formation CL for 2-weeks postpartum (OCC = 2)") # TH11 FINAL = 5.47E-01
    e_birthcontrol_cl_gluc  <- 1.46   ; label("Multiplicative scalar on glucuronide-formation CL for oral-contraceptive use (BIRTHCONTROL = 1, applied on top of the OCC effect)") # TH12 FINAL = 1.46E+00
    e_term_birth_cl_sulf <- 5.61 ; label("Sulphate-formation CL for term-birth subjects (TERM_BIRTH = 1) (L/h, on the linear scale; selector parameter rather than a multiplicative ratio)") # TH3 FINAL = 5.61E+00
    e_oc2_q2          <- 0.128  ; label("Multiplicative scalar on Q2 for 2-weeks postpartum (OCC = 2)") # TH15 FINAL = 1.28E-01
    e_urineflow_cl_renal     <- 0.00535 ; label("Slope of urine-flow effect on renal CL of paracetamol (L/h per mL/h, centered at 100 mL/h)") # TH16 FINAL = 5.35E-03

    # Inter-individual variability (IIV). NONMEM `$OMEGA` diagonal entries from
    # the FINAL PARAMETER ESTIMATE OMEGA block (lst lines 530-543). All etas
    # are Normal(0, var) on the log scale; the variance values are taken
    # verbatim. The sulphate-formation IIV slot was held FIX 0 in NONMEM (lst
    # OMEGA line 538-539, `0.00E+00`), so no eta is declared for cl_sulf.
    etalvc      ~ 0.0867  # OMEGA(1,1) FINAL = 8.67E-02 -- IIV V_central (log-normal variance)
    etalcl_gluc ~ 0.121   # OMEGA(2,2) FINAL = 1.21E-01 -- IIV CL_glucuronide
    etalcl_renal ~ 0.122  # OMEGA(4,4) FINAL = 1.22E-01 -- IIV CL_renal of unchanged paracetamol

    # Residual error. NONMEM `$SIGMA` diagonal entries from the FINAL
    # PARAMETER ESTIMATE SIGMA block (lst lines 549-567). The .mod $ERROR
    # uses linear-scale proportional and additive components: `Y1 = F*(1 +
    # ERR(1))` for plasma APAP at OCC != 5, `Y2 = F*(1 + ERR(5)) + ERR(6)`
    # for plasma APAP at OCC = 5, `Y5 = F*(1 + ERR(2))`, `Y6 = F*(1 +
    # ERR(3))`, `Y7 = F*(1 + ERR(4))` for the three urine outputs. NONMEM
    # SIGMA values are variances; nlmixr2 `prop()` / `add()` SDs are sqrt of
    # those variances. The OCC-conditional plasma error is implemented inside
    # `model()` by selecting between the OCC != 5 and OCC = 5 SD pairs based
    # on the `oc5` indicator; there is no canonical `<output>_<segment>_propSd`
    # name, so the OCC = 5 sub-arm parameters are named `Cc_oc5_propSd` /
    # `Cc_oc5_addSd` and a free-text justification appears in the vignette.
    CcpropSd        <- sqrt(0.0695)  ; label("Plasma APAP proportional residual SD at OCC != 5 (fraction)") # SIGMA(1,1) FINAL = 6.95E-02 -> SD 2.64E-01
    Cc_oc5_propSd   <- sqrt(0.0169)  ; label("Plasma APAP proportional residual SD at OCC = 5 (fraction)")  # SIGMA(5,5) FINAL = 1.69E-02 -> SD 1.30E-01
    Cc_oc5_addSd    <- sqrt(0.0160)  ; label("Plasma APAP additive residual SD at OCC = 5 (mg/L)")          # SIGMA(6,6) FINAL = 1.60E-02 -> SD 1.26E-01
    propSd_urineGluc <- sqrt(0.292)  ; label("Cumulative-urine APAP-glucuronide proportional residual SD (fraction)") # SIGMA(2,2) FINAL = 2.92E-01 -> SD 5.40E-01
    propSd_urineSulf <- sqrt(0.147)  ; label("Cumulative-urine APAP-sulphate proportional residual SD (fraction)")    # SIGMA(3,3) FINAL = 1.47E-01 -> SD 3.83E-01
    propSd_urineApap <- sqrt(0.152)  ; label("Cumulative-urine APAP proportional residual SD (fraction)")             # SIGMA(4,4) FINAL = 1.52E-01 -> SD 3.90E-01
  })
  model({
    # Decompose OCC into mutually-exclusive binary indicators. Used both for
    # typical-value covariate effects on V_central / CL_gluc / Q2 and for the
    # OCC = 5 conditional plasma residual-error switch (see end of this block).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc5 <- (OCC == 5)

    # Typical-value covariate multipliers. The `m_*` factors evaluate to 1.0
    # for the reference category (OCC >= 3 for V_central / CL_gluc, OCC != 2
    # for Q2, BIRTHCONTROL = 0 for CL_gluc) and to the corresponding `e_*`
    # parameter when the relevant indicator is 1. `m_occ_cl_gluc` and
    # `m_bc_cl_gluc` multiply each other as in the .mod $PK lines 39-43, where
    # BC = 1 scales the OCC-selected baseline by THETA(12) on top of the
    # OCC-specific adjustment (OCC = 1 -> TH2*TH9, OCC = 2 -> TH11*TH9,
    # OCC >= 3 -> TH9; then BC = 1 -> *TH12 applied uniformly).
    m_occ_vc       <- 1 + (e_oc1_vc - 1) * oc1
    m_occ_cl_gluc  <- 1 + (e_oc1_cl_gluc - 1) * oc1 + (e_oc2_cl_gluc - 1) * oc2
    m_bc_cl_gluc   <- 1 + (e_birthcontrol_cl_gluc - 1) * BIRTHCONTROL
    m_occ_q2       <- 1 + (e_oc2_q2 - 1) * oc2

    # Renal-CL covariate effect: linear in (URINE_FLOW - 100), with a
    # sentinel-zero rule -- when URINE_FLOW == 0 the slope contribution is
    # dropped (.mod $PK lines 45-47: `RCUF=(THETA(16)*(UF-100));
    # IF(UF.EQ.0) RCUF=0`; the source UF column is renamed to URINE_FLOW
    # canonical at data assembly).
    rcuf <- e_urineflow_cl_renal * (URINE_FLOW - 100) * (URINE_FLOW > 0)

    # Sulphate-formation CL: a binary selector between two typical values
    # rather than a multiplicative ratio (.mod $PK line 44:
    # `CLS=(TERM*THETA(3)+(1-TERM)*THETA(10))*EXP(ETA(3))`). Term-birth
    # subjects use `e_term_birth_cl_sulf` (TH3, log-untransformed because it
    # is not a typical structural CL with its own log slot), preterm-birth
    # subjects use `exp(lcl_sulf)` (TH10). The sulphate IIV slot was FIX 0
    # in NONMEM, so no eta multiplies this expression.
    cl_sulf <- TERM_BIRTH * e_term_birth_cl_sulf + (1 - TERM_BIRTH) * exp(lcl_sulf)

    # Individual PK parameters. WT enters Q1 with a fixed-1 allometric
    # exponent (.mod $PK line 33: `Q1=(THETA(6)*((BW/70)**1))`); no other
    # parameters carry an allometric term in this model.
    vc       <- exp(lvc + etalvc) * m_occ_vc
    vp       <- exp(lvp)
    vp2      <- exp(lvp2)
    q        <- exp(lq) * (WT / 70)
    q2       <- exp(lq2) * m_occ_q2
    cl_gluc  <- exp(lcl_gluc + etalcl_gluc) * m_occ_cl_gluc * m_bc_cl_gluc
    cl_renal <- (exp(lcl_renal) + rcuf) * exp(etalcl_renal)
    f_meta_renal <- exp(lf_meta_renal)

    # Metabolite plasma volumes set as 18% of central plasma volume (.mod
    # $PK lines 31-32: `V3=V1*0.18; V4=V3`).
    vc_gluc <- 0.18 * vc
    vc_sulf <- vc_gluc

    # Micro-constants reproducing the .mod $PK K12..K81 / K13..K46 lines.
    k12 <- q / vc
    k21 <- q / vp
    k18 <- q2 / vc
    k81 <- q2 / vp2
    k13 <- cl_gluc  / vc       # parent -> plasma APAP-glucuronide
    k14 <- cl_sulf  / vc       # parent -> plasma APAP-sulphate
    k17 <- cl_renal / vc       # parent -> urine APAP (unchanged-drug renal CL)
    k35 <- f_meta_renal * k17  # plasma APAP-glucuronide -> urine APAP-glucuronide
    k46 <- k35                 # plasma APAP-sulphate    -> urine APAP-sulphate (forced equal to k35 in the .mod)

    # ODEs (.mod $DES lines 73-80). Compartment names map verbatim to the
    # .mod's NCOMPARTMENTS=8 ordering: A(1)=central plasma APAP,
    # A(2)=peripheral1 plasma APAP, A(3)=plasma APAP-glucuronide,
    # A(4)=plasma APAP-sulphate, A(5)=urine APAP-glucuronide cumulative,
    # A(6)=urine APAP-sulphate cumulative, A(7)=urine APAP cumulative,
    # A(8)=peripheral2 plasma APAP. The `central_<metab>` plasma-metabolite
    # compartments use the canonical `<canonical>_<metab>` metabolite-suffix
    # convention with `gluc` and `sulf` newly registered in
    # R/conventions.R::registeredMetabolites for paracetamol Phase-II
    # conjugates. The three `urine_*` cumulative-elimination compartments are
    # paper-specific names that trigger `checkModelConventions()` warnings;
    # see the validation vignette's Assumptions and deviations section for
    # the justification (cumulative urine-elimination tracking is
    # mechanism-specific to this model and not yet a registered canonical).
    d/dt(central)      <- -k12 * central + k21 * peripheral1 -
                          k13 * central - k14 * central - k17 * central -
                          k18 * central + k81 * peripheral2
    d/dt(peripheral1)  <-  k12 * central - k21 * peripheral1
    d/dt(central_gluc) <-  k13 * central - k35 * central_gluc
    d/dt(central_sulf) <-  k14 * central - k46 * central_sulf
    d/dt(urine_gluc)   <-  k35 * central_gluc
    d/dt(urine_sulf)   <-  k46 * central_sulf
    d/dt(urine_apap)   <-  k17 * central
    d/dt(peripheral2)  <-  k18 * central - k81 * peripheral2

    # Plasma APAP concentration in the central compartment (mg/L). Dose units
    # are mg, V_central units are L, so central / vc is in mg/L.
    Cc <- central / vc
    # Cumulative urine outputs are reported as masses (mg) directly from the
    # cumulative compartments (.mod scale factors S5 = S6 = S7 = 1).
    urineGluc <- urine_gluc
    urineSulf <- urine_sulf
    urineApap <- urine_apap

    # Conditional plasma residual error: at OCC = 5 (non-pregnant volunteer
    # on birth control) the .mod $ERROR uses `Y2 = F*(1 + ERR(5)) + ERR(6)`,
    # i.e. combined proportional + additive; at every other OCC the model
    # uses `Y1 = F*(1 + ERR(1))`, proportional-only. The selector below
    # multiplexes the two SD pairs based on `oc5`. The additive SD evaluates
    # to zero for OCC != 5, so add(addSd_eff) collapses to no additive
    # contribution at those occasions.
    Cc_propSd_eff <- CcpropSd * (1 - oc5) + Cc_oc5_propSd * oc5
    Cc_addSd_eff  <- Cc_oc5_addSd * oc5
    Cc        ~ add(Cc_addSd_eff) + prop(Cc_propSd_eff)
    urineGluc ~ prop(propSd_urineGluc)
    urineSulf ~ prop(propSd_urineSulf)
    urineApap ~ prop(propSd_urineApap)
  })
}
