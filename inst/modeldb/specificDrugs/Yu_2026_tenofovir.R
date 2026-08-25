Yu_2026_tenofovir <- function() {
  description <- paste(
    "Eight-compartment semi-mechanistic population PK model describing every",
    "clinically relevant tenofovir moiety after either prodrug with a single",
    "unified parameter set: plasma tenofovir alafenamide (TAF), plasma tenofovir",
    "(TFV) arising from both TAF and tenofovir disoproxil fumarate (TDF), and the",
    "intracellular active anabolite tenofovir diphosphate (TFV-dp) in peripheral",
    "blood mononuclear cells (PBMCs) and in red cells sampled as dried blood spots",
    "(DBS). Plasma TAF is one-compartment; plasma TFV is two-compartment with a",
    "central and peripheral volume shared between the two prodrug routes. A",
    "fraction fm of the TAF eliminated from plasma converts to plasma TFV, of",
    "which a fraction frac converts immediately and the remainder passes through a",
    "transit compartment at rate ktr, producing the longer TFV half-life seen after",
    "TAF. Both intracellular pools are biophase compartments carried directly as",
    "concentrations: PBMC TFV-dp is fed by plasma TAF (the dominant route, via",
    "cathepsin A) and by plasma TFV, while red-cell TFV-dp is fed by plasma TFV",
    "alone because cathepsin A is absent from red cells. TDF is not carried as a",
    "state: its plasma half-life is about 24 seconds, so a 300 mg TDF dose is given",
    "directly into the tenofovir depot as its 136 mg TFV-equivalent. Pregnancy",
    "enters as second-trimester, third-trimester and postpartum fractional shifts",
    "on plasma TFV apparent clearance and on relative TAF bioavailability, carried",
    "over from the two companion pregnancy sub-models",
    "Yu_2026_tenofovir_pregnancy_tfv and Yu_2026_tenofovir_pregnancy_taf. Doses and",
    "amount states are in umol; all concentrations, including the two",
    "concentration-valued intracellular states, are in umol/L.",
    sep = " "
  )
  reference <- paste(
    "Yu Y, Brooks KM, Doncel GF, Best BM, Marzinke MA, Mirochnick M, Anderson P,",
    "Myer L, Celum C, Heffron R, Coleman J, Joseph Davey D, Hendrix CW, Momper JD,",
    "Bies R, Scott RK. Development of a Semi-Mechanistic Population Pharmacokinetic",
    "Model for Predicting Tenofovir Disoproxil Fumarate and Tenofovir Alafenamide",
    "Exposure in Plasma and Cellular Matrices During Pregnancy and Postpartum.",
    "Clin Pharmacokinet. 2026;65(1):133-148. doi:10.1007/s40262-025-01589-y.",
    "Structural parameters from Table 2; ordinary differential equations from",
    "Results section 3.2.1; pregnancy effects from Results section 3.2.2.",
    sep = " "
  )
  vignette <- "Yu_2026_tenofovir"
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  covariateData <- list(
    EGA = list(
      description        = "Estimated gestational age of the mother",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Used only to derive the paper's second- and third-trimester indicator",
        "variables TRI2 and TRI3; the model has no continuous EGA effect. Per the",
        "EGA register entry, a source paper that reports trimesters rather than",
        "weeks records the week values used per trimester here rather than",
        "introducing a trimester-indicator canonical. This model uses the standard",
        "obstetric boundaries TRI2 = 14 <= EGA < 28 weeks and TRI3 = EGA >= 28",
        "weeks. Yu 2026 does not print its own trimester cut-offs; the contributing",
        "IMPAACT P1026s protocol sampled the second trimester at 20-26 weeks and",
        "the third trimester at 30-38 weeks, both of which fall unambiguously",
        "inside these boundaries. Set EGA = 0 for non-pregnant and for postpartum",
        "subjects (the postpartum state is carried by TPP, not by EGA).",
        "No first-trimester data were available, so the model makes no prediction",
        "for 0 < EGA < 14 weeks; such a record is scored as non-pregnant.",
        sep = " "
      ),
      source_name        = "TRI2 / TRI3 (trimester indicator variables)"
    ),
    TPP = list(
      description        = "Time postpartum (time elapsed since delivery)",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Used only to derive the paper's postpartum indicator PP; the model has no",
        "continuous TPP effect and applies a single step change for any TPP > 0.",
        "Yu 2026 sampled the postpartum state at 6-12 weeks after delivery, so the",
        "estimated postpartum clearance and bioavailability describe that window",
        "and should not be read as immediately-post-delivery values. Set TPP = 0",
        "during pregnancy and for non-pregnant subjects.",
        sep = " "
      ),
      source_name        = "PP (postpartum indicator variable)"
    )
  )

  # Tenofovir alafenamide is the dosed parent for the TAF arm and is carried as
  # an ODE state, so it keeps the bare canonical depot / central / Cc names and
  # tenofovir takes the registered `_tfv` metabolite suffix (compartment-names.md
  # `tfv` entry, parent-wins rule; precedent Thoueille_2023_tenofovir_alafenamide.R).
  # `pbmc_tfvdp` and `rbc_tfvdp` are BIOPHASE states carried in concentration
  # units (umol/L), not amounts -- the paper writes their ODEs directly in
  # dC/dt form, so no volume divides them at the observation step.
  compartmentData <- list(
    depot           = list(analyte = "tenofovir alafenamide", units = "umol",   specimen = "administration site", verified = TRUE),
    central         = list(analyte = "tenofovir alafenamide", units = "umol",   specimen = "plasma",              verified = TRUE),
    transit1        = list(analyte = "tenofovir alafenamide-derived tenofovir precursor", units = "umol", specimen = "not applicable", verified = TRUE),
    depot_tfv       = list(analyte = "tenofovir",             units = "umol",   specimen = "administration site", verified = TRUE),
    central_tfv     = list(analyte = "tenofovir",             units = "umol",   specimen = "plasma",              verified = TRUE),
    peripheral1_tfv = list(analyte = "tenofovir",             units = "umol",   specimen = "plasma",              verified = TRUE),
    pbmc_tfvdp      = list(analyte = "tenofovir diphosphate", units = "umol/L", specimen = "blood cell", verified = TRUE),
    rbc_tfvdp       = list(analyte = "tenofovir diphosphate", units = "umol/L", specimen = "blood cell",          verified = TRUE)
  )

  population <- list(
    species       = "human",
    n_subjects    = 224L,
    n_studies     = 4L,
    sex_female_pct = 100,
    disease_state = paste(
      "Women receiving TDF 300 mg or TAF 10-25 mg once daily. The structural",
      "(non-pregnant) model was fitted to healthy HIV-negative non-pregnant women",
      "in CONRAD 137 (TDF 300 mg, TAF 10 mg and TAF 25 mg arms; 24 single-dose and",
      "73 multiple-dose participants) and in the directly-observed-therapy DOT-DBS",
      "study (TDF 300 mg, 28 participants). The pregnancy effects were identified",
      "in pregnant and postpartum women with HIV in IMPAACT P1026s (46 TDF-arm and",
      "25 TAF-arm participants) and IMPAACT 2026 (28 TAF-arm participants), all",
      "without a pharmacokinetic booster.",
      sep = " "
    ),
    dose_range = paste(
      "TDF 300 mg once daily (modelled as 136 mg = 474 umol tenofovir) or TAF",
      "25 mg once daily (52.5 umol); the CONRAD 137 multiple-dose phase also",
      "contributed a TAF 10 mg arm",
      sep = " "
    ),
    regions      = "USA and international IMPAACT network sites; CONRAD 137 and DOT-DBS conducted in the USA",
    co_medication = paste(
      "Participants co-administered TAF with a pharmacokinetic booster (cobicistat",
      "or ritonavir) were excluded, because P-glycoprotein inhibition raises TAF",
      "and TFV concentrations (Methods 2.1).",
      sep = " "
    ),
    notes = paste(
      "n_subjects is the sum of the per-arm participant counts in Table 1",
      "(24 + 73 + 28 + 46 + 25 + 28 = 224); the CONRAD 137 single-dose and",
      "multiple-dose phases may overlap, so the number of distinct women is",
      "somewhat lower. Age, weight, sex and race distributions are not reported in",
      "Yu 2026; the cohort is entirely female by design. Populations were",
      "deliberately not matched on age or race because covariate data were limited",
      "in the pooled datasets and neither factor had been identified as a",
      "significant covariate on TDF or TAF disposition (Methods 2.1). Plasma TAF",
      "was heavily censored (49.9-69.7% below the limit of quantification across",
      "the contributing arms), handled by the Beal M3 method; PBMC TFV-dp used the",
      "M1 method because each sample carried its own limit of quantification.",
      "TFV-dp in PBMCs was converted from fmol/million cells to umol/L using a",
      "single-PBMC volume of 282 fL, and TFV-dp in DBS from fmol/punch to umol/L",
      "using the red-cell volume and the average cell count of a 3-mm punch",
      "(Methods 2.4).",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Plasma tenofovir alafenamide (the dosed parent for the TAF arm --
    # bare canonical names). Yu 2026 Table 2.
    # ------------------------------------------------------------------
    lka <- log(2.44); label("First-order absorption rate constant Ka1 of oral TAF (1/h)")              # Table 2 (Ka1 = 2.44 /h, RSE 13%)
    lcl <- log(140);  label("Apparent clearance CL1/F of plasma TAF (L/h)")                            # Table 2 (CL1/F = 140 L/h, RSE 12%)
    lvc <- log(47.7); label("Apparent volume of distribution V2/F of plasma TAF (L)")                  # Table 2 (V2/F = 47.7 L, RSE 14%)

    # ------------------------------------------------------------------
    # Conversion of plasma TAF to plasma TFV (Results 3.2.1). Of the TAF
    # cleared from plasma, a fraction FTFV becomes plasma TFV; within
    # that, a fraction Ffast appears in plasma immediately and the rest
    # passes through the transit compartment at rate Kslow.
    # Canonical mapping (reuse-before-register; no new canonicals):
    #   FTFV  -> fm    (fraction metabolised; bare, because TAF is the
    #                   parent whose clearance is being split)
    #   Ffast -> frac  (registered generic bounded fraction)
    #   Kslow -> ktr   (registered transit-chain rate constant, the rate
    #                   out of the canonical `transit1` compartment)
    # ------------------------------------------------------------------
    lfm <- log(0.826); label("Fraction FTFV of eliminated plasma TAF that converts to plasma TFV (unitless)")            # Table 2 (FTFV = 82.6%, RSE 6%)
    lfrac   <- log(0.278); label("Fraction Ffast of the TAF-to-TFV conversion that reaches plasma immediately (unitless)")   # Table 2 (Ffast = 27.8%, RSE 13%)
    lktr    <- log(0.0212); label("Transit-compartment rate constant Kslow for delayed TAF-to-TFV conversion (1/h)")         # Table 2 (Kslow = 0.0212 /h, RSE 11%)

    # ------------------------------------------------------------------
    # Plasma tenofovir (`_tfv` metabolite suffix). Ka2 is the absorption
    # rate constant of the TFV-equivalent dose that stands in for oral
    # TDF; the central and peripheral volumes are shared by both prodrug
    # routes (Methods 2.4.2).
    # ------------------------------------------------------------------
    lka_tfv <- log(0.911); label("First-order absorption rate constant Ka2 of TFV after oral TDF (1/h)")   # Table 2 (Ka2 = 0.911 /h, RSE 26%)
    lcl_tfv <- log(52.7);  label("Apparent clearance CL2/F of plasma TFV, non-pregnant (L/h)")             # Table 2 (CL2/F = 52.7 L/h, RSE 4%)
    lvc_tfv <- log(425);   label("Apparent central volume of distribution V5/F of plasma TFV (L)")         # Table 2 (V5/F = 425 L, RSE 8%)
    lvp_tfv <- log(690);   label("Apparent peripheral volume of distribution V6/F of plasma TFV (L)")      # Table 2 (V6/F = 690 L, RSE 5%)
    lq_tfv  <- log(77.9);  label("Apparent intercompartmental clearance Q/F of plasma TFV (L/h)")          # Table 2 (Q/F = 77.9 L/h, RSE 11%)

    # ------------------------------------------------------------------
    # Intracellular TFV-dp biophase pools. The `_pbmc` influx / efflux
    # family was ratified for this extraction (operator sidecar
    # oare_PMC12783228 request-001 / response-001, q1 option A) as the
    # PBMC parallel of the registered `lkinf_rbc` / `lkeff_rbc` family;
    # the bare `lkinf_pbmc` is the parent-driven (plasma TAF) influx and
    # `lkinf_pbmc_tfv` the plasma-TFV-driven influx. Response q2 option B
    # keeps the red-cell pair on the registered bare names.
    # ------------------------------------------------------------------
    lkinf_pbmc     <- log(1.59);   label("Influx rate constant K(TAF-TFVdp-PBMC) from plasma TAF into PBMC TFV-dp (1/h)")   # Table 2 (K_TAF-TFV-dp-PBMC = 1.59 /h, RSE 11%)
    lkinf_pbmc_tfv <- log(0.0116); label("Influx rate constant K(TFV-TFVdp-PBMC) from plasma TFV into PBMC TFV-dp (1/h)")   # Table 2 (K_TFV-TFV-dp-PBMC = 0.0116 /h, RSE 9%)
    lkeff_pbmc     <- log(0.0127); label("Lumped loss rate constant Kel-PBMC out of the PBMC TFV-dp pool (1/h)")            # Table 2 (Kel-PBMC = 0.0127 /h, RSE 6%)
    lkinf_rbc      <- log(0.0068); label("Influx rate constant K(TFV-TFVdp-DBS) from plasma TFV into red-cell TFV-dp (1/h)") # Table 2 (K_TFV-TFV-dp-DBS = 0.0068 /h, RSE 5%)
    lkeff_rbc      <- log(0.0016); label("Lumped loss rate constant Kel-DBS out of the red-cell TFV-dp pool (1/h)")          # Table 2 (Kel-DBS = 0.0016 /h, RSE 2%)

    # ------------------------------------------------------------------
    # Pregnancy effects, carried into the semi-mechanistic model for the
    # clinical trial simulation (Methods 2.4.5). Results 3.2.2 states the
    # effects as fractional changes relative to the non-pregnant state,
    # which is exactly the additive covariate form of Methods 2.4.4
    # (TVP = th1 + TRI2*th2 + TRI3*th3 + PP*th4) written as a ratio.
    #
    # NOTE: the absolute clearances printed in ESM Table S1 imply smaller
    # trimester effects (+18.5% / +7.8% / -13.6%) than the percentages
    # quoted in the text. The text percentages are the ones the paper's
    # own simulations used -- see the vignette Errata, where the reported
    # steady-state DBS TFV-dp trough reductions arbitrate between them.
    # ------------------------------------------------------------------
    e_tri2_cl_tfv <- 0.249;  label("Fractional change in plasma TFV CL/F in the second trimester")            # Results 3.2.2 ("apparent clearance increased by 24.9% ... during the second ... trimester")
    e_tri3_cl_tfv <- 0.131;  label("Fractional change in plasma TFV CL/F in the third trimester")             # Results 3.2.2 ("... and 13.1% during the ... third trimester")
    e_pp_cl_tfv   <- -0.093; label("Fractional change in plasma TFV CL/F 6-12 weeks postpartum")              # Results 3.2.2 ("plasma TFV apparent clearance decreased by 9.3%")
    e_tri2_fdepot <- -0.173; label("Fractional change in relative TAF bioavailability in the second trimester") # Results 3.2.2 ("a decrease of 17.3% ... during the second ... trimester"); ESM Table S2 F = 0.827
    e_tri3_fdepot <- -0.051; label("Fractional change in relative TAF bioavailability in the third trimester")  # Results 3.2.2 ("... and 5.1% during the ... third trimester"); ESM Table S2 F = 0.949
    e_pp_fdepot   <- 0.18;   label("Fractional change in relative TAF bioavailability 6-12 weeks postpartum")   # Results 3.2.2 ("TAF bioavailability increased by 18% during the 6-12 week postpartum period"); ESM Table S2 F = 1.18

    # ------------------------------------------------------------------
    # Between-subject variability. Table 2 reports BSV as %CV for the
    # exponential model P = TVP * exp(eta) (Methods 2.4.1), so
    # omega^2 = log(CV^2 + 1). Q/F, FTFV, Kslow, V6/F and both red-cell
    # rate constants carry no BSV in Table 2.
    # ------------------------------------------------------------------
    etalcl            ~ 0.0878360  # CV 30.3%:  log(0.303^2 + 1)   (Table 2, BSV on CL1/F, RSE 37%, shrinkage 26%)
    etalvc            ~ 0.1106143  # CV 34.2%:  log(0.342^2 + 1)   (Table 2, BSV on V2/F, RSE 96%, shrinkage 50%)
    etalka            ~ 0.1100026  # CV 34.1%:  log(0.341^2 + 1)   (Table 2, BSV on Ka1, RSE 82%, shrinkage 44%)
    etalfrac          ~ 0.1881587  # CV 45.5%:  log(0.455^2 + 1)   (Table 2, BSV on Ffast, RSE 17%, shrinkage 30%)
    etalka_tfv        ~ 0.7610974  # CV 106.8%: log(1.068^2 + 1)   (Table 2, BSV on Ka2, RSE 37%, shrinkage 50%)
    etalcl_tfv        ~ 0.0786226  # CV 28.6%:  log(0.286^2 + 1)   (Table 2, BSV on CL2/F, RSE 17%, shrinkage 3%)
    etalvc_tfv        ~ 0.0759985  # CV 28.1%:  log(0.281^2 + 1)   (Table 2, BSV on V5/F, RSE 34%, shrinkage 33%)
    etalkinf_pbmc     ~ 0.0940330  # CV 31.4%:  log(0.314^2 + 1)   (Table 2, BSV on K_TAF-TFV-dp-PBMC, RSE 24%, shrinkage 30%)
    etalkinf_pbmc_tfv ~ 0.2167745  # CV 49.2%:  log(0.492^2 + 1)   (Table 2, BSV on K_TFV-TFV-dp-PBMC, RSE 31%, shrinkage 43%)
    etalkeff_pbmc     ~ 0.1009994  # CV 32.6%:  log(0.326^2 + 1)   (Table 2, BSV on Kel-PBMC, RSE 16%, shrinkage 18%)

    # ------------------------------------------------------------------
    # Residual variability, one model per analyte (Table 2). Plasma TFV
    # and red-cell TFV-dp carry a proportional term only.
    # ------------------------------------------------------------------
    propSd              <- 0.712;  label("Proportional residual error, plasma TAF (fraction)")           # Table 2 (PROP (TAF) = 71.2%, RSE 11%)
    addSd               <- 0.0007; label("Additive residual error, plasma TAF (umol/L)")                 # Table 2 (ADD (TAF) = 0.0007 umol/L, RSE 3%)
    propSd_tfv          <- 0.34;   label("Proportional residual error, plasma TFV (fraction)")           # Table 2 (PROP (TFV) = 34%, RSE 1%)
    propSd_Cpbmc_tfvdp  <- 0.306;  label("Proportional residual error, PBMC TFV-dp (fraction)")          # Table 2 (PROP (TFV-dp-PBMC) = 30.6%, RSE 2%)
    addSd_Cpbmc_tfvdp   <- 0.0272; label("Additive residual error, PBMC TFV-dp (umol/L)")                # Table 2 (ADD (TFV-dp-PBMC) = 0.0272 umol/L, RSE 2%)
    propSd_Crbc_tfvdp   <- 0.181;  label("Proportional residual error, red-cell (DBS) TFV-dp (fraction)") # Table 2 (PROP (TFV-dp-DBS) = 18.1%, RSE 2%)
  })

  model({
    # 1. Pregnancy state. The paper's TRI2 / TRI3 / PP indicators are
    #    derived here from the canonical EGA and TPP covariate columns.
    #    Order matters: postpartum is tested first so that a record with
    #    TPP > 0 is scored postpartum regardless of any residual EGA.
    TRI2 <- ifelse(TPP > 0, 0, ifelse(EGA >= 14 & EGA < 28, 1, 0))
    TRI3 <- ifelse(TPP > 0, 0, ifelse(EGA >= 28, 1, 0))
    PP   <- ifelse(TPP > 0, 1, 0)

    # 2. Individual parameters.
    ka     <- exp(lka + etalka)
    cl     <- exp(lcl + etalcl)
    vc     <- exp(lvc + etalvc)
    fm     <- exp(lfm)
    frac   <- exp(lfrac + etalfrac)
    ktr    <- exp(lktr)

    ka_tfv <- exp(lka_tfv + etalka_tfv)
    # Pregnancy raises plasma TFV apparent clearance in the second and
    # third trimesters and lowers it postpartum; the disposition
    # parameters are unaffected.
    cl_tfv <- exp(lcl_tfv + etalcl_tfv) *
      (1 + e_tri2_cl_tfv * TRI2 + e_tri3_cl_tfv * TRI3 + e_pp_cl_tfv * PP)
    vc_tfv <- exp(lvc_tfv + etalvc_tfv)
    vp_tfv <- exp(lvp_tfv)
    q_tfv  <- exp(lq_tfv)

    kinf_pbmc     <- exp(lkinf_pbmc + etalkinf_pbmc)
    kinf_pbmc_tfv <- exp(lkinf_pbmc_tfv + etalkinf_pbmc_tfv)
    keff_pbmc     <- exp(lkeff_pbmc + etalkeff_pbmc)
    kinf_rbc      <- exp(lkinf_rbc)
    keff_rbc      <- exp(lkeff_rbc)

    # 3. Micro-constants. Plasma TAF is cleared only into the TAF-derived
    #    pathways, so kel is the whole CL1/V2 flux out of `central`; the
    #    fraction fm of it becomes TFV and the remaining
    #    (1 - fm) leaves the system.
    kel     <- cl / vc
    kel_tfv <- cl_tfv / vc_tfv
    k12_tfv <- q_tfv / vc_tfv
    k21_tfv <- q_tfv / vp_tfv

    # 4. ODE system -- Yu 2026 Results 3.2.1, state by state.
    #    depot / central / transit1     = X1 / X2 / X3 (TAF side)
    #    depot_tfv / central_tfv /
    #      peripheral1_tfv              = X4 / X5 / X6 (TFV side)
    #    pbmc_tfvdp / rbc_tfvdp         = CPBMC / CDBS (concentration states)
    d/dt(depot)    <- -ka * depot
    d/dt(central)  <-  ka * depot - kel * central
    d/dt(transit1) <-  fm * (1 - frac) * kel * central - ktr * transit1

    d/dt(depot_tfv)   <- -ka_tfv * depot_tfv
    d/dt(central_tfv) <-  ka_tfv * depot_tfv +
                          fm * frac * kel * central +
                          ktr * transit1 -
                          kel_tfv * central_tfv -
                          k12_tfv * central_tfv + k21_tfv * peripheral1_tfv
    d/dt(peripheral1_tfv) <- k12_tfv * central_tfv - k21_tfv * peripheral1_tfv

    # Biophase pools. Both are driven by plasma CONCENTRATIONS and both
    # states are themselves concentrations, so the rate constants are in
    # 1/h and no volume term appears. Cathepsin A is absent from red
    # cells, so plasma TFV is the only source of red-cell TFV-dp.
    d/dt(pbmc_tfvdp) <- kinf_pbmc_tfv * (central_tfv / vc_tfv) +
                        kinf_pbmc * (central / vc) -
                        keff_pbmc * pbmc_tfvdp
    d/dt(rbc_tfvdp)  <- kinf_rbc * (central_tfv / vc_tfv) -
                        keff_rbc * rbc_tfvdp

    # 5. Relative bioavailability of oral TAF. Anchored at 1 in the
    #    non-pregnant state (Table 2 estimates CL1/F and V2/F, so F is
    #    already absorbed into them); pregnancy shifts it.
    f(depot) <- 1 + e_tri2_fdepot * TRI2 + e_tri3_fdepot * TRI3 + e_pp_fdepot * PP

    # 6. Observations. Plasma amounts are divided by their volumes; the
    #    two intracellular states are already concentrations.
    Cc           <- central / vc
    Cc_tfv       <- central_tfv / vc_tfv
    Cpbmc_tfvdp  <- pbmc_tfvdp
    Crbc_tfvdp   <- rbc_tfvdp

    Cc          ~ add(addSd) + prop(propSd)
    Cc_tfv      ~ prop(propSd_tfv)
    Cpbmc_tfvdp ~ add(addSd_Cpbmc_tfvdp) + prop(propSd_Cpbmc_tfvdp)
    Crbc_tfvdp  ~ prop(propSd_Crbc_tfvdp)
  })
}
