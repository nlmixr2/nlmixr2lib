# Joint parent-metabolite popPK model of tamoxifen and endoxifen at steady
# state in adult breast-cancer patients, extracted from the DDMORE Foundation
# Model Repository entry DDMODEL00000212 (the Ter Heine 2014 BJCP paper).

TerHeine_2014_tamoxifen <- function() {
  description <- "Joint parent-metabolite population PK model for tamoxifen and endoxifen at steady state in adult breast-cancer patients, with CYP2D6 and CYP3A4/5 individual-activity covariates on the endoxifen-formation clearance"
  reference <- paste(
    "Ter Heine R, Binkhorst L, de Graan AJM, et al. (2014).",
    "Population pharmacokinetic modelling to assess the impact of CYP2D6",
    "and CYP3A metabolic phenotypes on the pharmacokinetics of tamoxifen",
    "and endoxifen.",
    "Br J Clin Pharmacol 78(3):572-586.",
    "doi:10.1111/bcp.12388.",
    "DDMORE Foundation Model Repository: DDMODEL00000212.",
    sep = " "
  )
  vignette <- "TerHeine_2014_tamoxifen"
  ddmore_id <- "DDMODEL00000212"
  replicate_of <- NULL
  units <- list(time = "h", dosing = "mg", concentration = "nmol/L")

  covariateData <- list(
    CYP2D6 = list(
      description        = "CYP2D6 individual metabolic-activity score (dextromethorphan-probe model-based individual CL value)",
      units              = "ng/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-invariant per-subject covariate. Centered on the population median 1560 ng/L (.mdl GROUP_VARIABLES). Applied as `(CYP2D6 / 1560)^e_cyp2d6_cl_endox` on the endoxifen-formation clearance. Source: a separate dextromethorphan-probe popPK study run on the same patients (DDMORE RDF: 'CYP2D6 and CYP3A4/5 phenotypes (dextromethorphan model-based individual CL values)').",
      source_name        = "CYP2D6"
    ),
    CYP3A4 = list(
      description        = "CYP3A4 + CYP3A5 individual metabolic-activity score (dextromethorphan-probe model-based individual CL value)",
      units              = "ng/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-invariant per-subject covariate. Centered on the population median 44.7 ng/L (.mdl GROUP_VARIABLES). Applied as `(CYP3A4 / 44.7)^e_cyp3a4_cl_endox` on the endoxifen-formation clearance. The probe (dextromethorphan N-demethylation) cannot separate CYP3A4 from CYP3A5, so the column carries the combined CYP3A4 + CYP3A5 activity. The DDMORE source-data column name is `CYP3A4`; the source publication describes the same column as `CYP3A4/5`.",
      source_name        = "CYP3A4"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 39L,
    n_studies      = 1L,
    age_range      = "22 - 71 years",
    age_median     = "53 years",
    weight_range   = "48.5 - 114 kg",
    weight_median  = "72.7 kg",
    height_range   = "1.56 - 1.79 m",
    height_median  = "1.69 m",
    sex_female_pct = 100,
    disease_state  = paste(
      "Hormone-receptor-positive breast cancer on chronic oral tamoxifen at",
      "steady state. Patients on moderate or strong CYP3A inhibitors /",
      "inducers or ABCB1 / ABCG2 modulators were excluded."
    ),
    dose_range     = paste(
      "Tamoxifen 20 mg orally once daily in 70% of subjects and 40 mg",
      "orally once daily in 30% (Table 1 of the source publication;",
      "11 of 39 evaluable subjects on the 40 mg regimen). 349 tamoxifen",
      "and 331 endoxifen concentrations available."
    ),
    regions        = "The Netherlands (Erasmus MC Cancer Institute, Rotterdam; single-center)",
    notes          = paste(
      "Demographics from Table 1 of the source publication, which is now",
      "on disk (the earlier TODO markers from the DDMORE-only extraction",
      "have been resolved against the full PDF). 40 patients enrolled,",
      "1 dropped due to obstruction of a venous cannula; n = 39 evaluable.",
      "Cohort sex is 100% female because the indication is hormone-",
      "receptor-positive breast cancer on tamoxifen. Dutch Trial Registry",
      "NTR1751."
    )
  )

  ini({
    # Structural-parameter values come from the original Ter Heine 2014
    # publication's Table 2, embedded in the DDMORE bundle as a one-page PDF
    # `Output_real_tamoxifen.pdf`. The .mdl/.mod $THETA / $OMEGA / $SIGMA
    # blocks list the DDMORE convertor's INITIAL values, not the final
    # estimates, so the PDF is the authoritative source.
    lka <- log(1.90)
    label("First-order absorption rate constant from gut to hepatic compartment (1/h)")
    # PDF Table 2: k12 = 1.90 1/h, RSE 20.2%
    ltlag <- log(0.455)
    label("Absorption lag time (h)")
    # PDF Table 2: tlag = 0.455 h, RSE 10.4%
    lcl <- log(9.34)
    label("Apparent elimination clearance of tamoxifen, CL20/F (L/h)")
    # PDF Table 2: CL20/F = 9.34 L/h, RSE 6.20%
    lvc <- log(753)
    label("Apparent volume of distribution of tamoxifen, V2/F (L)")
    # PDF Table 2: V2/F = 753 L, RSE 9.00%
    lq <- log(61.8)
    label("Apparent hepatic flow Q connecting tamoxifen central and the algebraic hepatic compartment, Q/F (L/h)")
    # PDF Table 2: Q/F = 61.8 L/h, RSE 65.4%
    lcl_endox <- log(0.324)
    label("Apparent formation clearance of endoxifen from tamoxifen, CL23/F (L/h)")
    # PDF Table 2: CL23/F = 0.324 L/h, RSE 17.0%

    e_cyp2d6_cl_endox <- 0.262
    label("CYP2D6 power-law effect on endoxifen formation clearance (unitless)")
    # PDF Table 2: theta_CYP2D6 = 0.262, RSE 14.0%
    e_cyp3a4_cl_endox <- 0.157
    label("CYP3A4/5 power-law effect on endoxifen formation clearance (unitless)")
    # PDF Table 2: theta_CYP3A4/5 = 0.157, RSE 72.0%

    # IIV is log-normal. Paper reports CV%, which converts to internal (log)
    # variance via omega^2 = log(1 + (CV/100)^2). CL20 and V2 share a
    # correlated block (rho = 61.3% per PDF Table 2); covariance =
    # rho * sqrt(var_cl * var_vc) = 0.613 * sqrt(0.1336 * 0.0683) = 0.0586.
    etalcl + etalvc ~ c(0.1336, 0.0586, 0.0683)
    # PDF Table 2: omega CL20/F %CV = 37.8%, omega V2/F %CV = 26.6%, rho = 61.3%
    etalcl_endox ~ 0.0630
    # PDF Table 2: omega CL23/F %CV = 25.5%; var = log(1 + 0.255^2) = 0.0630

    # Residual error - proportional only. The original publication included a
    # residual-error correlation (rho = 62.2% per PDF Table 2) between the
    # parent and metabolite residuals; the DDMORE-implemented version drops
    # that correlation and parameterises both as standard deviations
    # (Model_Accommodations.txt notes the deviation). We follow the
    # DDMORE-implemented model (no residual correlation) so the model file
    # matches the .mdl / .mod structure shipped in the DDMORE bundle.
    propSd <- 0.138
    label("Proportional residual SD for tamoxifen plasma concentration (fraction)")
    # PDF Table 2: sigma_TAM %CV = 13.8%, RSE 11.3%
    propSd_endox <- 0.189
    label("Proportional residual SD for endoxifen plasma concentration (fraction)")
    # PDF Table 2: sigma_ENDX %CV = 18.9%, RSE 10.1%
  })

  model({
    # Centered phenotype log-ratios. The reference values 1560 ng/L and 44.7
    # ng/L are the population medians of the dextromethorphan-probe-derived
    # individual CL values, hard-coded in the .mdl GROUP_VARIABLES block.
    log_pheno_cyp2d6 <- log(CYP2D6 / 1560)
    log_pheno_cyp3a4 <- log(CYP3A4 / 44.7)

    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    cl_endx_form <- exp(lcl_endox + etalcl_endox +
                          e_cyp2d6_cl_endox * log_pheno_cyp2d6 +
                          e_cyp3a4_cl_endox * log_pheno_cyp3a4)

    # Fixed literature constants for endoxifen disposition. Source: Ahmad et
    # al. (2010) Clin Pharmacol Ther 88(6):814-817, doi:10.1038/clpt.2010.222
    # -- cited in the .mdl GROUP_VARIABLES block as `Ahmad et al, CPT Vol 88,
    # 2010`.
    cl_endx_elim <- 5.1
    vc_endox <- 400
    kel_endox <- cl_endx_elim / vc_endox
    # Molar masses for the compartment-mass-to-plasma-nM unit conversion. The
    # model tracks endoxifen amount in tamoxifen-mass-equivalents and converts
    # to nM via Mendx; because Mtam ~= Mendx within ~0.5%, the implicit
    # mass-conserving 1 mg parent -> 1 mg-equivalent metabolite assumption
    # introduces negligible error (matches the original .mdl scaling).
    mtam <- 371.51456
    mendx <- 373.48738

    # Algebraic hepatic compartment. The liver is at quasi-steady state: it
    # receives gut input (ratein) plus parent back-flow from systemic
    # (q * c2), exits via Q (return to systemic) and via the formation
    # clearance to endoxifen (cl_endx_form). At steady state:
    #   c_hep = (ratein + q * c2) / (q + cl_endx_form)
    # Rearranging into the parent-central and metabolite-central ODEs gives
    # the structure below.
    ratein <- ka * depot
    c2 <- central / vc
    c_hep <- (ratein + q * c2) / (q + cl_endx_form)

    d/dt(depot) <- -ratein
    d/dt(central) <- q * c_hep - q * c2 - cl * c2
    d/dt(central_endox) <- cl_endx_form * c_hep - kel_endox * central_endox

    # Standard NONMEM ALAG mechanism on the depot. The DDMORE-converted .mdl
    # implements the lag inside $DES via an `IF (T >= ALAG1)` switch which
    # gives slightly different ODE behaviour during 0 < t < lag (the gut
    # mass stays at the initial dose level rather than appearing at t = lag).
    # Both forms yield the same systemic exposure once t >= lag.
    lag(depot) <- exp(ltlag)

    # Plasma concentrations in nmol/L. With dose AMT in mg and Vc in L, the
    # compartment concentration in mg/L converts to nM as
    #   nM = (mg/L) * 1e6 / Mw,
    # so the scale factor is V * Mw / 1e6.
    Cc <- central / (vc * mtam / 1e6)
    Cc_endox <- central_endox / (vc_endox * mendx / 1e6)

    Cc ~ prop(propSd)
    Cc_endox ~ prop(propSd_endox)
  })
}
