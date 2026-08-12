Mou_2024_aripiprazole <- function() {
  description <- "Compartmental reduction of a PK-Sim whole-body PBPK model for aripiprazole and its active metabolite dehydro-aripiprazole in adults, with CYP2D6 phenotype, mild hepatic impairment (Child-Pugh A) and risperidone coadministration effects (Mou 2024)"
  reference <- "Mou F, Huang Z, Cheng Y, Zhao X, Sun X, Li H, Yu S. Physiologically based pharmacokinetic modeling to predict the effect of risperidone on aripiprazole pharmacokinetics in subjects with different CYP2D6 genotypes and individuals with hepatic impairment. Ther Adv Drug Saf. 2024;15:20420986241303432. doi:10.1177/20420986241303432"
  vignette <- "Mou_2024_aripiprazole"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "aripiprazole", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "aripiprazole", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "aripiprazole", units = "mg", specimen = "plasma", verified = TRUE),
    central_dehydro = list(analyte = "dehydro-aripiprazole", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_dehydro = list(analyte = "dehydro-aripiprazole", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reference weight 70 kg. Mou 2024 Table 1 reports aripiprazole and dehydro-aripiprazole clearances",
        "in mL/h/kg, so every clearance is scaled linearly with WT (exponent 1, not allometric 0.75).",
        "Volumes are scaled linearly with WT as well; that is a reduction assumption, not a paper value",
        "(Mou 2024 reports no volume of distribution). Supplemental Table S1 gives the virtual population",
        "weight range 48.0-171.2 kg (median 78.0 kg); 70 kg is used as the reference because it is the only",
        "value at which the paper's own Table 1 clearance and its own reported AUC are mutually consistent",
        "with a bioavailability <= 1 (see the vignette Errata)."
      ),
      source_name        = "WT"
    ),
    CYP2D6_EM = list(
      description        = "CYP2D6 extensive (normal) metabolizer indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (intermediate metabolizer, when CYP2D6_PM = 0 and CYP2D6_UM = 0)",
      notes              = paste(
        "Mou 2024 calls this stratum the normal metabolizer (NM). Paired with CYP2D6_PM and CYP2D6_UM;",
        "the intermediate metabolizer (IM) is the implicit reference (all three indicators 0), following",
        "the CYP2D6_PM / CYP2D6_EM pattern established by Sherwin 2012."
      ),
      source_name        = "NM"
    ),
    CYP2D6_PM = list(
      description        = "CYP2D6 poor metabolizer indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (intermediate metabolizer, when CYP2D6_EM = 0 and CYP2D6_UM = 0)",
      notes              = paste(
        "Mou 2024 Table 2 sets both CYP2D6 Kcat values to 0 for the PM stratum (zero CYP2D6 activity) and",
        "additionally raises the two CYP3A4 Kcat values by ~1.34-fold (footnote a, parameter identification)."
      ),
      source_name        = "PM"
    ),
    CYP2D6_UM = list(
      description        = "CYP2D6 ultrarapid metabolizer indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (intermediate metabolizer, when CYP2D6_EM = 0 and CYP2D6_PM = 0)",
      notes              = paste(
        "Mou 2024 Table 2 assigns the ultrarapid metabolizer stratum exactly the same Kcat values as the",
        "normal metabolizer, so UM and NM are kinetically identical in this model. The paper's own Table S3",
        "reports identical AUC for NM and UM (748.25 umol*min/L), confirming the encoding."
      ),
      source_name        = "UM"
    ),
    HEPIMP_MILD = list(
      description        = "Mild hepatic impairment (Child-Pugh A) indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal hepatic function)",
      notes              = paste(
        "Child-Pugh A. Mou 2024 Table 3 lists the underlying PK-Sim system parameter changes (fractional",
        "liver mass 0.69, CYP2D6 abundance 0.76, CYP3A4 abundance 0.79, renal blood flow 0.88) and an",
        "aripiprazole unbound fraction of 1.23% versus 1.00%. Those factors are NOT composed mechanistically",
        "here: a well-stirred composition of them predicts roughly +50% AUC, whereas the paper's own",
        "Supplemental Table S3 reports +7.8%. The empirical Table S3 ratio is encoded instead; see the",
        "vignette Errata for the discrepancy."
      ),
      source_name        = "CP-A"
    ),
    CONMED_RISPERIDONE = list(
      description        = "Concomitant risperidone coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (aripiprazole monotherapy)",
      notes              = paste(
        "Mou 2024 models risperidone as a CYP2D6 inhibitor (inhibition constant Ki; the value is not",
        "reported anywhere in the paper or supplement). The predicted interaction is negligible: the",
        "aripiprazole AUC ratio was 1.00, 1.01 and 1.01 for risperidone 2, 4 and 6 mg. Encoded as a single",
        "flat multiplicative effect on total aripiprazole clearance taken from Supplemental Table S3",
        "(DDI 755.26 vs DGI NM 748.25 umol*min/L), i.e. dose-independent over 2-6 mg risperidone."
      ),
      source_name        = "risperidone"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 500,
    n_studies      = 4,
    age_range      = "18-65 years (virtual population); 19.0-77.0 years across the validation trials",
    weight_range   = "48.0-171.2 kg (virtual population; median 78.0 kg)",
    height_range   = "151.8-191.7 cm (virtual population; median 170.8 cm)",
    disease_state  = "healthy adults; one validation cohort with mild hepatic impairment (Child-Pugh A)",
    dose_range     = "10-30 mg oral aripiprazole, single and multiple dose",
    regions        = "virtual population built from NHANES; Mexican-American / White male demographics assumed where the source trials did not report them",
    notes          = paste(
      "Simulated virtual population of 500 individuals built in PK-Sim (Supplemental Table S1).",
      "Model validation used four clinical studies (Supplemental Table S2): Kneller 2021 CYP2D6-phenotyped",
      "single 10 mg dose (n = 28 NM, 21 IM, 1 PM, 5 UM), Mallikaarjun 2004 multiple 10 mg (n = 8) and 30 mg",
      "(n = 7), Waade 2009 multiple 20 mg with risperidone (n = 10), and Mallikaarjun 2008 single 15 mg in",
      "mild hepatic impairment (n = 8). The 500-subject count is the simulated population, not an observed",
      "sample size."
    )
  )

  ini({
    # =====================================================================
    # ABSORPTION AND DISTRIBUTION
    #
    # Mou 2024 reports NO volume of distribution, no half-life and no
    # absorption rate constant anywhere in the paper or the supplement:
    # distribution is a PK-Sim partition-coefficient database output. The
    # five parameters in this block were therefore recovered by digitising
    # the predicted mean concentration-time curves of Figure 2(a) (oral
    # single 10 mg aripiprazole, CYP2D6 normal metabolizer, aripiprazole
    # and dehydro-aripiprazole) and fitting this reduction to them with
    # every clearance held fixed at its published value. The fit reproduces
    # the digitised parent curve to 1.5% and the metabolite curve to 3.6%
    # (log RMSE), and reproduces the paper's own independently reported
    # AUC to 0.2%. See the vignette "Source trace" and "Errata" sections.
    # PROVENANCE: digitised from Mou 2024 Figure 2(a); not a paper value.
    # =====================================================================
    lka <- fixed(log(2.3638)); label("Aripiprazole first-order absorption rate constant (log 1/h)")  # digitised from Figure 2(a); Mou 2024 reports only Weibull dissolution (b = 0.91, T50 = 5 min)
    lfdepot <- fixed(log(0.974)); label("Aripiprazole oral bioavailability (log fraction)")  # derived: reconciles Table 1 CL_H (47.96 mL/h/kg at 70 kg) with the reported single-dose AUC of 773.85 umol*min/L after 20 mg
    lvc <- fixed(log(190.8)); label("Aripiprazole central volume of distribution (log L)")  # digitised from Figure 2(a)
    lvp <- fixed(log(32.0)); label("Aripiprazole peripheral volume of distribution (log L)")  # digitised from Figure 2(a)
    lq <- fixed(log(4.629)); label("Aripiprazole intercompartmental clearance (log L/h)")  # digitised from Figure 2(a)

    lvc_dehydro <- fixed(log(43.2)); label("Dehydro-aripiprazole central volume of distribution (log L)")  # digitised from Figure 2(a)
    lvp_dehydro <- fixed(log(193.3)); label("Dehydro-aripiprazole peripheral volume of distribution (log L)")  # digitised from Figure 2(a)
    lq_dehydro <- fixed(log(9.952)); label("Dehydro-aripiprazole intercompartmental clearance (log L/h)")  # digitised from Figure 2(a)

    # =====================================================================
    # CLEARANCE (all published)
    # Table 1 reports specific clearances in mL/h/kg; values below are at
    # the 70 kg reference weight. Supplemental section 1.3 independently
    # states the aripiprazole plasma clearance as 0.8 mL/min/kg, which is
    # 48.0 mL/h/kg and matches Table 1's 47.96 mL/h/kg.
    # =====================================================================
    lcl_nonren <- fixed(log(3.3572)); label("Aripiprazole hepatic (metabolic) clearance at 70 kg, CYP2D6 NM (log L/h)")  # Table 1 CL_H 47.96 mL/h/kg x 70 kg; cross-checked against Supplemental 1.3 "plasma clearance 0.8 mL/min/kg"
    lcl_renal <- fixed(log(0.0028)); label("Aripiprazole renal clearance at 70 kg (log L/h)")  # Table 1 CL_R 0.04 mL/h/kg x 70 kg; Supplemental 1.3 "renal clearance < 1%"
    lcl_nonren_dehydro <- fixed(log(2.80)); label("Dehydro-aripiprazole hepatic clearance at 70 kg (log L/h)")  # Table 1 dehydro-aripiprazole CL_H 40.00 mL/h/kg x 70 kg
    lcl_renal_dehydro <- fixed(log(0.0042)); label("Dehydro-aripiprazole renal clearance at 70 kg (log L/h)")  # Table 1 dehydro-aripiprazole CL_R 0.06 mL/h/kg x 70 kg

    # =====================================================================
    # METABOLIC PATHWAY SPLIT (Table 2 + Table 3)
    # Each pathway clearance is Kcat * [E] / Km (Supplemental equation 2),
    # with Table 2 Kcat / Km and Table 3 normal-liver enzyme abundances
    # ([CYP2D6] = 0.40 umol/L, [CYP3A4] = 4.32 umol/L). The four fractions
    # below are those four pathway clearances divided by their total, for
    # the CYP2D6 normal metabolizer. They sum to 1 and independently
    # reproduce the fractions stated in Supplemental section 1.3:
    # CYP2D6 0.0902 + 0.3439 = 0.434 (paper: 0.43), CYP3A4 0.3133 + 0.2527
    # = 0.566 (paper: 0.56), and dehydrogenation 0.0902 + 0.3133 = 0.403
    # (paper: "total dehydro-aripiprazole accounts for 40% of aripiprazole").
    # =====================================================================
    fmet2d6Dehydro <- fixed(0.09017); label("Fraction of CYP2D6-NM metabolic clearance via CYP2D6 dehydrogenation (unitless)")  # Table 2 Kcat 14.99 /min; 14.99*0.40/26.20 divided by the four-pathway total
    fmet2d6Other <- fixed(0.34387); label("Fraction of CYP2D6-NM metabolic clearance via CYP2D6 hydroxylation (unitless)")  # Table 2 Kcat 57.17 /min
    fmet3a4Dehydro <- fixed(0.31327); label("Fraction of CYP2D6-NM metabolic clearance via CYP3A4 dehydrogenation (unitless)")  # Table 2 Kcat 54.85 /min; 54.85*4.32/298.00 divided by the four-pathway total
    fmet3a4Other <- fixed(0.25268); label("Fraction of CYP2D6-NM metabolic clearance via CYP3A4 other routes (unitless)")  # Table 2 Kcat 44.24 /min

    # CYP2D6 phenotype activity relative to the normal metabolizer.
    # Table 2 gives IM Kcat 10.49 and 40.02 versus NM 14.99 and 57.17;
    # both ratios are 0.700. PM is 0 by construction (Table 2, reference 29:
    # "the enzyme activity of CYP2D6 poor metabolizer being 0"). UM equals NM.
    rel2d6Im <- fixed(0.6999); label("CYP2D6 activity in intermediate metabolizers relative to normal metabolizers (unitless)")  # Table 2: 10.49/14.99 = 0.6998 and 40.02/57.17 = 0.7000

    # CYP3A4 up-adjustment in CYP2D6 poor metabolizers (Table 2 footnote a,
    # "Calculated/parameter identification"). Applied to the CYP3A4
    # pathways only, and only in the PM stratum.
    rel3a4PmDehydro <- fixed(1.3453); label("CYP3A4 dehydrogenation activity in CYP2D6 poor metabolizers relative to other phenotypes (unitless)")  # Table 2: 73.79/54.85
    rel3a4PmOther <- fixed(1.3422); label("CYP3A4 other-route activity in CYP2D6 poor metabolizers relative to other phenotypes (unitless)")  # Table 2: 59.38/44.24

    # Molar mass conversion, aripiprazole -> dehydro-aripiprazole.
    mwRatioDehydro <- fixed(0.99554); label("Dehydro-aripiprazole / aripiprazole molar mass ratio (unitless)")  # Table 1: 446.4 / 448.4 g/mol

    # =====================================================================
    # COVARIATE EFFECTS ON CLEARANCE
    # Both are exposure ratios read off Supplemental Table S3 and inverted
    # to clearance multipliers, because the paper reports the outcome
    # (AUC change) rather than a clearance coefficient.
    # =====================================================================
    e_hepimp_cl <- fixed(0.9279); label("Child-Pugh A multiplicative effect on clearance (unitless)")  # Supplemental Table S3: DDZI AUC 806.39 vs DGI NM 748.25 umol*min/L -> AUC x 1.0777 -> CL x 1/1.0777
    e_risp_cl <- fixed(0.9907); label("Risperidone coadministration multiplicative effect on aripiprazole clearance (unitless)")  # Supplemental Table S3: DDI AUC 755.26 vs DGI NM 748.25 umol*min/L -> AUC x 1.00937 -> CL x 1/1.00937

    # =====================================================================
    # RESIDUAL ERROR
    # Mou 2024 report no residual-error model. The PBPK model is a forward
    # predictor assessed by MPE / MAPE against observed concentrations
    # (Methods "Evaluation of the PBPK model"), not by a fitted sigma.
    # Non-fitted placeholders per the in-repo PBPK convention
    # (Luo_2024_perindopril_pbpk.R, An_2012_mitoxantrone_human_pbpk.R).
    # =====================================================================
    propSd <- fixed(0.10); label("Proportional residual error placeholder, aripiprazole (fraction)")  # not reported by Mou 2024
    propSd_dehydro <- fixed(0.10); label("Proportional residual error placeholder, dehydro-aripiprazole (fraction)")  # not reported by Mou 2024
  })

  model({
    # ------------------------------------------------------------------
    # 1. CYP2D6 phenotype decode. CYP2D6_EM / CYP2D6_PM / CYP2D6_UM are
    #    mutually exclusive; all three 0 selects the intermediate
    #    metabolizer (IM) reference stratum. Written as a linear form in
    #    the indicators so that it reduces to the Table 2 activities:
    #    IM -> rel2d6Im, EM and UM -> 1, PM -> 0.
    # ------------------------------------------------------------------
    rel2d6 <- rel2d6Im +
      (1 - rel2d6Im) * (CYP2D6_EM + CYP2D6_UM) -
      rel2d6Im * CYP2D6_PM

    # CYP3A4 activity is unchanged except in the PM stratum (Table 2 footnote a).
    rel3a4Dehydro <- 1 + (rel3a4PmDehydro - 1) * CYP2D6_PM
    rel3a4Other <- 1 + (rel3a4PmOther - 1) * CYP2D6_PM

    # ------------------------------------------------------------------
    # 2. Clearances. Linear body-weight scaling because Table 1 reports
    #    every clearance per kg. The two covariate multipliers act on the
    #    whole aripiprazole clearance (Supplemental Table S3 reports them
    #    as total-exposure changes, not pathway-specific ones); the
    #    hepatic-impairment factor also acts on the metabolite, whose
    #    clearance is likewise hepatic.
    # ------------------------------------------------------------------
    covcl <- (1 + (e_hepimp_cl - 1) * HEPIMP_MILD) *
      (1 + (e_risp_cl - 1) * CONMED_RISPERIDONE)
    covclm <- 1 + (e_hepimp_cl - 1) * HEPIMP_MILD

    clmet <- exp(lcl_nonren) * (WT / 70) * covcl

    # Four Table 2 pathway clearances.
    cl2d6dehydro <- clmet * fmet2d6Dehydro * rel2d6
    cl2d6other <- clmet * fmet2d6Other * rel2d6
    cl3a4dehydro <- clmet * fmet3a4Dehydro * rel3a4Dehydro
    cl3a4other <- clmet * fmet3a4Other * rel3a4Other

    # Clearance forming dehydro-aripiprazole (both dehydrogenation routes).
    clform <- cl2d6dehydro + cl3a4dehydro

    cl <- cl2d6dehydro + cl2d6other + cl3a4dehydro + cl3a4other +
      exp(lcl_renal) * (WT / 70) * covcl

    cl_dehydro <- exp(lcl_nonren_dehydro) * (WT / 70) * covclm +
      exp(lcl_renal_dehydro) * (WT / 70)

    ka <- exp(lka)
    vc <- exp(lvc) * (WT / 70)
    vp <- exp(lvp) * (WT / 70)
    q <- exp(lq) * (WT / 70)
    vc_dehydro <- exp(lvc_dehydro) * (WT / 70)
    vp_dehydro <- exp(lvp_dehydro) * (WT / 70)
    q_dehydro <- exp(lq_dehydro) * (WT / 70)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    kelm <- cl_dehydro / vc_dehydro
    km12 <- q_dehydro / vc_dehydro
    km21 <- q_dehydro / vp_dehydro

    # ------------------------------------------------------------------
    # 3. Disposition. Dehydro-aripiprazole is formed from systemic
    #    aripiprazole at the two dehydrogenation clearances, converted
    #    from aripiprazole mass to dehydro-aripiprazole mass by the molar
    #    mass ratio.
    # ------------------------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    d/dt(central_dehydro) <- clform * (central / vc) * mwRatioDehydro -
      kelm * central_dehydro - km12 * central_dehydro + km21 * peripheral1_dehydro
    d/dt(peripheral1_dehydro) <- km12 * central_dehydro - km21 * peripheral1_dehydro

    f(depot) <- exp(lfdepot)

    # Dose in mg and volume in L give mg/L = ug/mL; x 1000 -> ng/mL.
    Cc <- 1000 * central / vc
    Cc_dehydro <- 1000 * central_dehydro / vc_dehydro

    Cc ~ prop(propSd)
    Cc_dehydro ~ prop(propSd_dehydro)
  })
}
