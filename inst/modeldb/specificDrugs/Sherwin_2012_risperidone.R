# One-compartment mixture-model population PK of oral risperidone and its
# active metabolite (+/-)-9-hydroxyrisperidone in children and adolescents
# with neuropsychiatric disorders, with three CYP2D6 metabolizer
# subpopulations (PM / IM / EM) and allometric weight scaling on apparent
# clearance and volume (Sherwin 2012, Ther Drug Monit 34(5):535-544;
# doi:10.1097/FTD.0b013e318261c240).

Sherwin_2012_risperidone <- function() {
  description <- paste(
    "One-compartment parent-plus-metabolite population PK model for oral",
    "risperidone and its active metabolite (+/-)-9-hydroxyrisperidone in 45",
    "children and adolescents (aged 3-18.3 years, 16.8-110 kg) with",
    "neuropsychiatric disorders treated with maintenance oral risperidone",
    "(Sherwin 2012). First-order absorption (Ka fixed) into a single central",
    "compartment with first-order elimination; the fraction of risperidone",
    "metabolized to (+/-)-9-hydroxyrisperidone (KF) feeds a single",
    "metabolite compartment whose apparent volume of distribution is set",
    "equal to the parent apparent volume per Table 2 footnote (a).",
    "A mixture model with three CYP2D6 metabolizer subpopulations (poor",
    "PM, intermediate IM, extensive EM) yields subpopulation-specific",
    "apparent oral clearances and metabolite formation fractions; KF in IM",
    "subjects is fixed at 1 to stabilize the model per the paper's",
    "Mixture Model section. Allometric scaling (exponent 0.75 for CL/F and",
    "CLM/F, exponent 1 for Vd/F, reference 70 kg) is applied to all three",
    "subpopulations' clearance estimates and to the shared apparent volume.",
    "Inter-individual variability is reported separately for each",
    "subpopulation's CL/F (PM, IM, EM), for the metabolite CLM/F, and for",
    "the shared Vd/F; a combined additive-plus-proportional residual error",
    "is reported separately for risperidone and (+/-)-9-hydroxyrisperidone",
    "plasma concentrations.",
    sep = " "
  )
  reference <- paste(
    "Sherwin CMT, Saldana SN, Bies RR, Aman MG, Vinks AA (2012).",
    "Population pharmacokinetic modeling of risperidone and",
    "9-hydroxyrisperidone to estimate CYP2D6 subpopulations in children and",
    "adolescents. Ther Drug Monit 34(5):535-544.",
    "doi:10.1097/FTD.0b013e318261c240.",
    sep = " "
  )
  vignette <- "Sherwin_2012_risperidone"
  units    <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at baseline in Sherwin 2012; cohort range 16.8-110 kg",
        "(mean 43, SD 20.2; Table 1). Reference 70 kg with fixed allometric",
        "exponents 0.75 on apparent CL/F and CLM/F and 1.0 on apparent Vd/F",
        "(Methods, Equation 3 and surrounding text)."
      ),
      source_name        = "WT"
    ),
    CYP2D6_PM = list(
      description        = "CYP2D6 poor-metabolizer phenotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (intermediate or extensive metabolizer; both CYP2D6_PM and CYP2D6_EM = 0 indicates IM)",
      notes              = paste(
        "1 = subject is a CYP2D6 poor metabolizer, 0 otherwise. Paired with",
        "CYP2D6_EM to encode the three-level PM / IM / EM phenotype with two",
        "binary indicators on the SLCO1B1_HAP15_HET / SLCO1B1_HAP15_HOM",
        "pattern; intermediate metabolizer (IM) is the implicit reference",
        "(CYP2D6_PM = 0 and CYP2D6_EM = 0). In the Sherwin 2012 cohort the",
        "mixture-model assignment estimated 37.2% PM, 15.9% IM, and 46.9% EM",
        "(Table 2 P1 and P2; IM proportion = 1 - P1 - P2). For subjects with",
        "a genotyped CYP2D6 status, assign PM / IM / EM by the published",
        "genotype-to-phenotype mapping; for ungenotyped subjects (38% of the",
        "Sherwin 2012 cohort) the source paper inferred phenotype from the",
        "mixture-model posterior rather than fixing it from external data."
      ),
      source_name        = "P1 (mixture-model PM subpopulation fraction); Table 1 'CYP2D6 phenotype' for genotyped subjects"
    ),
    CYP2D6_EM = list(
      description        = "CYP2D6 extensive-metabolizer phenotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (intermediate or poor metabolizer; both CYP2D6_PM and CYP2D6_EM = 0 indicates IM)",
      notes              = paste(
        "1 = subject is a CYP2D6 extensive metabolizer, 0 otherwise. Paired",
        "with CYP2D6_PM; IM is the implicit reference (both indicators = 0).",
        "In the Sherwin 2012 cohort the mixture-model assignment estimated",
        "37.2% PM, 15.9% IM, and 46.9% EM (Table 2 P1 and P2; IM proportion",
        "= 1 - P1 - P2)."
      ),
      source_name        = "P2 (mixture-model EM subpopulation fraction); Table 1 'CYP2D6 phenotype' for genotyped subjects"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 45L,
    n_observations  = 497L,
    n_studies       = 3L,
    age_range       = "3-18.3 years (mean 9.6, SD 3.7; Table 1)",
    age_median      = "9.6 years (mean +/- SD reported, not median)",
    weight_range    = "16.8-110 kg (mean 43, SD 20.2; Table 1)",
    weight_median   = "43 kg (mean +/- SD reported, not median)",
    sex_female_pct  = 11.1,
    race_ethnicity  = c(White = 93.3, `White + Black (mixed)` = 6.7),
    ethnicity       = "100% non-Hispanic (Table 1)",
    disease_state   = paste(
      "Neuropsychiatric disorders treated with stable maintenance oral",
      "risperidone; autistic disorder was the predominant diagnosis",
      "(Results)."
    ),
    dose_range      = paste(
      "Oral risperidone tablet (n = 34) or liquid (n = 6) at the patient's",
      "clinician-prescribed dose; total daily dose 0.25-6.00 mg (mean 2.0,",
      "SD 1.5); most subjects (n = 39) dosed twice daily (Results)."
    ),
    regions         = "USA (Cincinnati OH, Columbus OH, Cleveland OH)",
    cyp2d6_distribution = paste(
      "Mixture-model assignment in the final model: 37.2% PM, 15.9% IM,",
      "46.9% EM (Table 2 P1 and P2 with IM = 1 - P1 - P2). For the 28",
      "subjects with confirmed CYP2D6 genotype, the observed phenotype",
      "distribution was 15 EM, 6 IM, 7 PM (Table 1)."
    ),
    notes           = paste(
      "Steady-state oral risperidone after at least 4 weeks at the same",
      "dose; samples taken pre-dose, 1, 2, 4, and 7 hours post-dose",
      "(Materials and Methods). One outlier subject was removed from the",
      "final dataset (Results). Concentrations of risperidone and the (+)",
      "and (-) enantiomers of 9-hydroxyrisperidone were measured separately",
      "by enantioselective LC-MS/MS (lower limit of quantification 0.2",
      "ng/mL); the model fits the combined racemic (+/-)-9-hydroxyrisperidone",
      "as a single 'metabolite' output following the paper's text. Source",
      "studies (n = 3) include the two open-label CCHMC / OSU / Rainbow",
      "Babies prior investigations enrolled June 2001-May 2003 and the",
      "December 2008-June 2010 enrichment protocol that prospectively",
      "recruited confirmed CYP2D6 PMs."
    )
  )

  ini({
    # Structural absorption -- Ka was fixed at 2.5 1/h in the final
    # mixture model (Table 2 superscript b: fixed). The paper's earlier
    # base-model section mentions 2.6 1/h; the Table 2 final mixture-
    # model value (2.5) is authoritative.
    lka <- fixed(log(2.5)); label("Absorption rate constant Ka (1/h); fixed per Table 2")  # Sherwin 2012 Table 2: Ka = 2.5 (Fixed)

    # Shared apparent volume of distribution. The paper constrains
    # Vd/F (risperidone) = VdM/F (9-OH risperidone) per Table 2
    # footnote (a). Allometric exponent on volume is 1 with reference
    # 70 kg (Methods Eq. 3 and surrounding text).
    lvc <- log(77.3); label("Apparent central volume of distribution Vd/F at 70 kg reference (L); shared with metabolite VdM/F (Table 2 footnote a)")  # Sherwin 2012 Table 2: Vd/F = VdM/F = 77.3 L (SE 12.6%; bootstrap 95% CI 55.3-101)

    # Subpopulation-specific apparent oral clearances. Each metabolizer
    # subpopulation has its own typical CL/F estimated by the mixture
    # model; allometric scaling (exponent 0.75, reference 70 kg) is
    # applied to all three. The PM / IM / EM suffix here marks the
    # CYP2D6 metabolizer subpopulation (not a metabolite), and is
    # selected at model() time by the paired CYP2D6_PM / CYP2D6_EM
    # binary indicators (both 0 = IM reference).
    lcl_pm <- log(9.38); label("Apparent oral clearance in CYP2D6 PMs, CL/F at 70 kg reference (L/h)")  # Sherwin 2012 Table 2: CL/F in PM = 9.38 L/h (SE 9.2%; bootstrap 95% CI 7.3-10.3)
    lcl_im <- log(29.2); label("Apparent oral clearance in CYP2D6 IMs, CL/F at 70 kg reference (L/h)")  # Sherwin 2012 Table 2: CL/F in IM = 29.2 L/h (SE 10.1%; bootstrap 95% CI 0.46-68.9 -- wide CI reflects small IM stratum)
    lcl_em <- log(37.4); label("Apparent oral clearance in CYP2D6 EMs, CL/F at 70 kg reference (L/h)")  # Sherwin 2012 Table 2: CL/F in EM = 37.4 L/h (SE 5.02%; bootstrap 95% CI 20.4-42.7)

    # Apparent clearance of (+/-)-9-hydroxyrisperidone (metabolite),
    # not subpopulation-specific. Allometric scaling 0.75 with
    # reference 70 kg per the paper's text ("9-OH risperidone
    # clearance was estimated to be 9.0 L/h/70 kg", Results).
    lcl_9oh <- log(9.0); label("Apparent clearance of (+/-)-9-hydroxyrisperidone CLM/F at 70 kg reference (L/h/70kg)")  # Sherwin 2012 Table 2: CLM/F = 9.0 L/h/70kg (SE 41.3%; bootstrap 95% CI 9-9.2)

    # Subpopulation-specific fraction of risperidone metabolized to
    # (+/-)-9-hydroxyrisperidone (KF). Linear-scale parameter (KF in
    # (0, 1]); estimated for PM and EM, fixed at 1 for IM per the
    # paper's Mixture Model section ("the fraction of risperidone to
    # (+/-) 9-OH risperidone for IMs was fixed to help stabilize the
    # model. When this fraction was unfixed, there was poor model fit
    # and the estimated relative clearances for CYP2D6 IMs were
    # unrealistically high.").
    kf_pm <- 0.16; label("Fraction of risperidone metabolized to (+/-)-9-hydroxyrisperidone in CYP2D6 PMs (unitless, 0-1)")  # Sherwin 2012 Table 2: KF-PM = 0.16 (SE 17.7%; bootstrap 95% CI 0.01-0.5)
    kf_em <- 0.13; label("Fraction of risperidone metabolized to (+/-)-9-hydroxyrisperidone in CYP2D6 EMs (unitless, 0-1)")  # Sherwin 2012 Table 2: KF-EM = 0.13 (SE 36.1%; bootstrap 95% CI 0.03-0.5)
    kf_im <- fixed(1); label("Fraction of risperidone metabolized to (+/-)-9-hydroxyrisperidone in CYP2D6 IMs (unitless, fixed at 1)")  # Sherwin 2012 Table 2: KF-IM = 1 (Fixed)

    # Inter-individual variability (NONMEM OMEGA, variance scale).
    # IIV is reported separately for each metabolizer subpopulation's
    # CL/F (Table 2), for the metabolite CLM/F, and for the shared
    # Vd/F. No IIV is reported on Ka (fixed parameter) or on KF.
    # The eta variances are taken from the "Parameter estimates"
    # column directly; the "CV%" column in Table 2 is sqrt(omega^2)
    # in percent, consistent with these values (e.g. sqrt(0.07) =
    # 26.5%, matching the table). Note: the BSV (omega)-CL EM
    # bootstrap distribution recovers a noticeably larger mean (0.65,
    # 95% CI 0.5-0.7) than the point estimate (0.06), reflecting a
    # poorly identifiable EM-stratum variance; the point estimate is
    # used here as the primary value per the paper's reporting.
    etalvc     ~ 0.20  # Sherwin 2012 Table 2: BSV(omega)-Vd = 0.20 (SE 7.9%; CV% 44.5%)
    etalcl_pm  ~ 0.07  # Sherwin 2012 Table 2: BSV(omega)-CL PM = 0.07 (SE 9.6%; CV% 26.5%)
    etalcl_im  ~ 0.03  # Sherwin 2012 Table 2: BSV(omega)-CL IM = 0.03 (SE 9.1%; CV% 18.7%)
    etalcl_em  ~ 0.06  # Sherwin 2012 Table 2: BSV(omega)-CL EM = 0.06 (SE 21.1%; CV% 25.5%; bootstrap mean 0.65, 95% CI 0.5-0.7 -- see vignette Assumptions and deviations)
    etalcl_9oh ~ 0.50  # Sherwin 2012 Table 2: BSV(omega)-CLM = 0.50 (SE 25.6%; CV% 70.7%)

    # Residual error -- combined additive + constant coefficient-of-
    # variation (proportional) error model for risperidone and (+/-)-9-
    # hydroxyrisperidone separately (Methods Equation 2 and surrounding
    # text). Table 2 reports the NONMEM sigma^2 variance directly in
    # the "Parameter estimates" column for both the proportional ("CV")
    # and additive ("SD") components; the "CV%" column shows
    # sqrt(sigma^2) in the original units (percent CV for proportional,
    # ng/mL for additive). nlmixr2 propSd / addSd are on the SD scale
    # so each is sqrt of the reported variance.
    propSd     <- sqrt(0.08); label("Proportional residual error for risperidone (fraction; SD = sqrt(NONMEM sigma^2))")  # Sherwin 2012 Table 2: RUV(sigma)-CV RSP = 0.08 (variance; CV% sqrt scale 29.4%; propSd = sqrt(0.08) = 0.283)
    addSd      <- sqrt(0.5);  label("Additive residual error for risperidone (ng/mL; SD = sqrt(NONMEM sigma^2))")  # Sherwin 2012 Table 2: RUV(sigma)-SD RSP = 0.5 (variance, ng/mL units; sqrt scale 0.71 ng/mL; addSd = sqrt(0.5) = 0.707)
    propSd_9oh <- sqrt(0.47); label("Proportional residual error for (+/-)-9-hydroxyrisperidone (fraction; SD = sqrt(NONMEM sigma^2))")  # Sherwin 2012 Table 2: RUV(sigma)-CV 9-OH = 0.47 (variance; CV% sqrt scale 68.6%; propSd_9oh = sqrt(0.47) = 0.686)
    addSd_9oh  <- sqrt(0.44); label("Additive residual error for (+/-)-9-hydroxyrisperidone (ng/mL; SD = sqrt(NONMEM sigma^2))")  # Sherwin 2012 Table 2: RUV(sigma)-SD 9-OH = 0.44 (variance, ng/mL units; sqrt scale 0.67 ng/mL; addSd_9oh = sqrt(0.44) = 0.663)
  })

  model({
    # CYP2D6 phenotype indicator-gating. PM and EM are explicit binary
    # covariates; IM is the implicit reference (both = 0). For any
    # subject only one of the three terms below evaluates to 1, so
    # only the corresponding subpopulation's structural CL/F and KF
    # (and the corresponding eta) contribute. nlmixr2 samples all
    # five etas for every subject; indicator gating zeroes the
    # contribution of the two non-active phenotypes, so the simulated
    # trajectory uses only the eta of the subject's assigned
    # phenotype -- consistent with the source NONMEM mixture model's
    # per-subject ETA assignment.
    ind_im <- (1 - CYP2D6_PM - CYP2D6_EM)

    # Allometric scaling. Exponent 0.75 on clearances and 1.0 on
    # volumes, reference body weight 70 kg, fixed exponents per
    # Methods Equation 3 ("by fixing the exponents in the allometric
    # model to 0.75 for CL/F and to 1 for Vd/F").
    wt_cl <- (WT / 70)^0.75
    wt_v  <- (WT / 70)

    # Subpopulation-specific apparent oral clearances of risperidone
    # with subpopulation-specific BSV. log-normal IIV applied on the
    # linear-scale typical value via cl_i = TVCL * exp(eta_i).
    cl_pm <- exp(lcl_pm + etalcl_pm) * wt_cl
    cl_im <- exp(lcl_im + etalcl_im) * wt_cl
    cl_em <- exp(lcl_em + etalcl_em) * wt_cl

    # Indicator-gated active CL and KF. Only one indicator is 1 per
    # subject, so cl equals the corresponding phenotype's cl_*; kf
    # equals the corresponding phenotype's kf_*.
    cl <- cl_pm * CYP2D6_PM + cl_em * CYP2D6_EM + cl_im * ind_im
    kf <- kf_pm * CYP2D6_PM + kf_em * CYP2D6_EM + kf_im * ind_im

    # Metabolite clearance (not phenotype-specific) and shared apparent
    # volume. VdM/F = Vd/F per Table 2 footnote (a).
    cl_9oh <- exp(lcl_9oh + etalcl_9oh) * wt_cl
    vc     <- exp(lvc     + etalvc)     * wt_v
    vc_9oh <- vc

    ka <- exp(lka)

    # Disposition. d/dt(depot): first-order absorption out of the gut.
    # d/dt(central): risperidone in the central compartment, eliminated
    # at rate (cl / vc). d/dt(central_9oh): metabolite formed at rate
    # kf * (cl / vc) * central (mass-fraction conversion as in the
    # source NONMEM $DES; the paper does not apply a molar correction
    # and the MWs of risperidone 410.5 and 9-hydroxyrisperidone 426.5
    # differ by only 4%) and eliminated at rate (cl_9oh / vc_9oh).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - (cl / vc) * central
    d/dt(central_9oh) <-  kf * (cl / vc) * central - (cl_9oh / vc_9oh) * central_9oh

    # Plasma concentrations in ng/mL. Dose in mg, volumes in L:
    # central / vc has units mg/L = ug/mL; multiply by 1000 to obtain
    # ng/mL, matching Table 1 observed concentration units (e.g.
    # risperidone 6.5 +/- 6.4 ng/mL; 9-hydroxyrisperidone enantiomers
    # 8.4 +/- 7.6 and 3.84 +/- 2.97 ng/mL).
    Cc     <- 1000 * central     / vc
    Cc_9oh <- 1000 * central_9oh / vc_9oh

    # Combined additive + proportional residual error, separate for
    # parent and metabolite per Methods Equation 2 and Table 2.
    Cc     ~ prop(propSd)     + add(addSd)
    Cc_9oh ~ prop(propSd_9oh) + add(addSd_9oh)
  })
}
