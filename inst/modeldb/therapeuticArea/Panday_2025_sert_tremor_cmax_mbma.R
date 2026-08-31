Panday_2025_sert_tremor_cmax_mbma <- function() {
  description <- "MBMA. Translational Emax model-based meta-analysis relating the incidence proportion of tremor (a characteristic manifestation of serotonin syndrome) to predicted brain serotonin-reuptake-transporter (SERT) target coverage, pooled across 20 SERT inhibitors from five drug classes (SSRI, SNRI, SMS, TCA, opioid). Consumes drug-level literature inputs supplied as data columns (steady-state MAXIMUM total plasma concentration, plasma fraction unbound, molecular weight, unbound brain-to-unbound-plasma partition coefficient and free-corrected SERT IC50) and returns a per-arm tremor incidence in percent. Suitable for simulating study-arm-level summary outcomes only; there is no PK time course, no dosing event and no individual-level prediction. This file carries the Cmax-based fit; see Panday_2025_sert_tremor_cavg_mbma for the Cavg-based fit, which the paper designates as preferred."
  reference <- "Panday SK, Lang BJ, Kapitanov GI, Subramanian K, Klopp-Schulze L, Venkatakrishnan K, Zutshi A, Alnaif AE. A Translational Model-Based Meta-Analysis to Predict Tremor Incidence Associated with Serotonin Reuptake Transporter Inhibition. Clin Pharmacol Ther. 2025. doi:10.1002/cpt.3696."
  vignette <- "Panday_2025_sert_tremor_mbma"

  # Algebraic MBMA: no rxode2 dose events are consumed (every input arrives as a
  # covariate data column) and the single output `tremor` is an incidence
  # proportion in percent rather than a drug concentration. The `units` entries
  # below carry placeholder strings chosen so that checkModelConventions() sees a
  # dimensionally consistent dose-vs-concentration pair; the same device is used
  # by Yoshioka_2018_FXa_inhibitors_mbma, the other algebraic MBMA in the library.
  units <- list(
    time          = "h",
    dosing        = "percent",
    concentration = "percent/percent"
  )

  covariateData <- list(
    CMAX = list(
      description        = "Study-arm-level steady-state MAXIMUM total (bound + unbound) plasma concentration of the SERT inhibitor at the dose for which tremor incidence was reported.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Total, not unbound: Eq. S1 applies FU internally. Collated per drug and dose from the primary literature and FDA labels (Panday 2025 Table S2, parameter rows 'plasma steady state total Cmax' / 'steady state plasma total Cmax' / 'plasma steady state cmax'). For drugs with a SERT-active metabolite (amitriptyline/nortriptyline, fluoxetine/norfluoxetine, venlafaxine/O-desmethylvenlafaxine) the paper computed target coverage separately for parent and metabolite and SUMMED the two coverages (Supplemental Methods, 'Estimating SERT target coverage for drugs with active metabolites'); that summation is NOT performed inside this model -- evaluate the model once per moiety and add the resulting coverages, or use the parameter set fitted excluding drugs with active metabolites (the default here), which is the paper's preferred relationship.",
      source_name        = "C_pl,ss (Panday 2025 Supplemental Methods Eq. S1); 'Cmax,ub,ss' after the unbound/brain corrections"
    ),
    FU = list(
      description        = "Fraction of the SERT inhibitor unbound in human plasma.",
      units              = "(unitless fraction)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Unitless fraction (e.g. 0.02 for sertraline, 0.7 for desvenlafaxine), NOT percent. Derived from analysis of human plasma samples (Panday 2025 Supplemental Methods, 'Model development for estimating the in vivo brain exposure'); values collated in Table S2 parameter rows 'plasma fu'. Enters Eq. S1 multiplicatively.",
      source_name        = "f_u,pl (Panday 2025 Supplemental Methods Eq. S1)"
    ),
    MW = list(
      description        = "Molecular weight of the SERT inhibitor, used to convert the reported mass-per-volume plasma concentration to a molar concentration so that it is commensurate with the nanomolar in vitro potency.",
      units              = "g/mol",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters Eq. S1 as the factor 1000/MW, which converts ng/mL to nM. Values in Panday 2025 Table S2 parameter rows 'MW' (e.g. 298.4 for vortioxetine, 306.2 for sertraline, 263.37 for desvenlafaxine).",
      source_name        = "MW (Panday 2025 Supplemental Methods Eq. S1)"
    ),
    KPUU_BRAIN = list(
      description        = "Unbound-brain to unbound-plasma partition coefficient (Kp,uu) of the SERT inhibitor.",
      units              = "(unitless ratio)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Typically derived from rat or mouse studies (Panday 2025 Supplemental Methods). The paper preferred the 'Kp,uu * f_u,pl' route (Eq. S2b, encoded here) over the 'Kp * f_u,br' route (Eq. S2a) because it considered plasma fraction-unbound measurements more reliable than brain fraction-unbound measurements. When only a total brain-to-plasma ratio Kp and a brain fraction unbound f_u,br are available, the two routes are algebraically identical for KPUU_BRAIN = Kp * f_u,br / f_u,pl, so supply that product. Brain fraction-unbound values measured in brain homogenate were first pH-partition corrected per Friden et al. (Eqs. S3-S6) before that conversion; the corrected values are tabulated in Panday 2025 Table S2 parameter rows 'brain fu corrected'.",
      source_name        = "K_p,uu (Panday 2025 Supplemental Methods Eq. S2b)"
    ),
    IC50_SERT = list(
      description        = "In vitro SERT inhibitory potency of the drug, corrected to the free (unbound) drug concentration in the assay.",
      units              = "nM",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Where an assay was reported to contain serum, the measured IC50 was corrected to a free-drug basis via Eq. S7 before use; where no serum was reported the measured potency was taken to be the free potency. Where several literature values existed for one drug the GEOMETRIC MEAN was used (Panday 2025 Supplemental Methods). Values in Panday 2025 Table S2 parameter rows 'SERT IC50' / 'SERT IC50 derived'. Substituting a free-corrected SERT K_D here instead of an IC50 is a supported variant of the paper's analysis, but it requires the K_D-based parameter set from Table S4 rather than the IC50-based default shipped in ini() -- see the vignette.",
      source_name        = "IC50_SERT (Panday 2025 Supplemental Methods Eqs. S7 and S8)"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 29677L,
    n_studies        = 33L,
    age_range        = "Adults; per-study demographics are not re-tabulated by the meta-analysis (each arm's source study is cited by PubMed ID or FDA label URL in Panday 2025 Table S2).",
    weight_range     = "Not a model covariate; one source study normalised tapentadol IR exposure to a 77 kg body weight (Panday 2025 Table S2).",
    sex_female_pct   = NA_real_,
    disease_state    = "Mixed: major depressive disorder, generalised anxiety disorder, fibromyalgia, obsessive-compulsive disorder, chronic pain and opioid dependence -- whichever indication the source study of each arm enrolled. The modelled endpoint is the treatment-emergent incidence proportion of tremor, used as a representative manifestation of serotonin syndrome.",
    dose_range       = "33 drug/dose arms across 20 SERT inhibitors (Panday 2025 Table S1): amitriptyline+nortriptyline 75 and 111 mg QD; citalopram 40 mg QD; desvenlafaxine 50, 100, 200 and 400 mg QD; duloxetine 60 and 90 mg QD; escitalopram 10 mg QD; fluoxetine+norfluoxetine 20, 25 and 40 mg QD; fluvoxamine 100 mg BID and 100 mg QD; imipramine 110 mg BID; methadone 100 mg QD; milnacipran 50 mg BID and 100 and 200 mg QD; nortriptyline 100 mg QD; paroxetine (and CR) 12.5, 22, 22.9, 25 and 30 mg QD; sertraline 50, 100, 144 and 200 mg QD; tapentadol (ER, IR) 75 mg Q6H and 175 mg BID; tramadol ER 200 mg QD; venlafaxine+O-desmethylvenlafaxine 50 and 75 mg BID and 75, 225 and 375 mg QD; venlafaxine XR 85, 150 and 225 mg QD; vilazodone 40 mg QD; vortioxetine 5 and 10 mg QD.",
    regions          = "Not reported; the arms are drawn from the published literature and from US FDA product labels.",
    drug_classes     = "Five SERT-inhibitor classes represented: selective serotonin reuptake inhibitors (SSRI), serotonin and norepinephrine reuptake inhibitors (SNRI), serotonin modulator and stimulator (SMS), tricyclic antidepressants (TCA) and opioids (Panday 2025 Figure 1b and Table S1).",
    notes            = "Summary-level MBMA: the modelled observations are per-arm tremor incidence proportions weighted by study-group size, NOT individual-patient data. n_subjects (29,677) and n_studies (33) are the sum and the count of the treatment-arm group sizes tabulated in Panday 2025 Table S2 ('tremor percent - treatment' rows); the paper itself does not print a pooled total. Placebo-arm group sizes total 20,381 records but double-count shared placebo groups across dose levels of the same trial, so they are not added here. Clinical tremor incidences and plasma fractions unbound are human; the brain distribution parameters (Kp, Kp,uu, brain fraction unbound) are predominantly rat or mouse, and the SERT potencies are in vitro -- this cross-species integration is the 'translational' element of the analysis."
  )

  ini({
    # Panday 2025 Eq. S9:  tremor = Emax * x / (EC50 + x) + B0,  x = SERT target coverage.
    # Emax, EC50 and B0 are assumed lognormal and are parameterised as
    # exp(theta_parameter) with theta_parameter normal (Panday 2025 Supplemental
    # Methods, final paragraph), so the theta values below ARE the log-scale
    # estimates printed in Table S2 and require no back-transformation.
    #
    # Default parameter set = Cmax,ub,ss / SERT IC50, data selection "Excluding
    # drugs with active metabolites" -- the Cmax counterpart of the paper's
    # preferred Cavg relationship. The paper drew "similar conclusions" from the
    # Cmax analyses and found the CI widths comparable, but concluded that no
    # determination of whether tremor tracks Cavg or Cmax could be made
    # (Discussion). The two other Cmax/IC50 data selections reported in Table S2
    # can be substituted with ini():
    #     All drugs                 lemax = 2.00  (SE 0.229), lec50 = 0.584  (SE 0.942)
    #     Subset to SSRIs and SNRIs lemax = 1.92  (SE 0.311), lec50 = 0.218  (SE 1.18)
    # and the SERT K_D-based Cmax fits from Table S4 are
    #     All drugs                 lemax = 1.71  (SE 0.691), lec50 = 1.42   (SE 2.34)
    #     Excluding active metab.   lemax = 1.66  (SE 0.739), lec50 = 1.42   (SE 2.71)
    #     Subset to SSRIs and SNRIs lemax = 1.74  (SE 1.1),   lec50 = 1.42   (SE 3.32)
    # le0 is 0.293 in every one of those fits.
    lemax <- 1.84
    label("Log maximum attributable tremor incidence, log(Emax) (Emax in percentage points)")            # Panday 2025 Table S2, Cmax,ub,ss / SERT IC50, excluding drugs with active metabolites, weighted estimate 1.84 (SE 0.257); back-transformed median exp(1.84) = 6.31% in Table S3

    lec50 <- 0.300
    label("Log SERT target coverage producing half-maximal tremor incidence, log(EC50) (EC50 unitless)") # Panday 2025 Table S2, Cmax,ub,ss / SERT IC50, excluding drugs with active metabolites, weighted estimate 0.300 (SE 1.27); back-transformed median exp(0.300) = 1.35 in Table S3

    # B0 was NOT estimated: it was set to the group-size-weighted arithmetic mean
    # of the placebo tremor incidences across every placebo arm in the dataset.
    # Main text Results: "In each analysis, a weighted arithmetic mean of placebo
    # tremor incidences (1.34%) was used as a baseline"; Supplemental Methods:
    # "A weighted average of placebo tremor incidences (1.34%) for all placebo
    # groups in the dataset was used as the B0." Consistently, Table S2 reports
    # theta_B0 = 0.293 identically for all twelve fits, and exp(0.293) = 1.340%.
    # The SE printed alongside it (1.20 for this fit) is a bootstrap dispersion,
    # not evidence of estimation -- hence fixed() here.
    le0 <- fixed(0.293)
    label("Log baseline (placebo) tremor incidence, log(B0) (B0 in percent)")                            # Panday 2025 Table S2 (theta_B0 = 0.293 in all fits) and Supplemental Methods Eq. S9; exp(0.293) = 1.34%, the group-size-weighted mean placebo tremor incidence
  })

  model({
    # Lognormal back-transforms (Panday 2025 Supplemental Methods: "parameter = exp(theta_parameter)")
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    b0   <- exp(le0)

    # Eq. S1 -- unbound steady-state plasma concentration, converted ng/mL -> nM
    cplssub <- CMAX * FU * 1000 / MW

    # Eq. S2b -- unbound steady-state brain concentration (the route the paper preferred)
    cbrssub <- cplssub * KPUU_BRAIN

    # Eq. S8 / main-text Eq. 1 -- SERT target coverage (unitless)
    sertcov <- cbrssub / IC50_SERT

    # Eq. S9 / main-text Eq. 2 -- Emax relationship for tremor incidence (percent)
    tremor <- emax * sertcov / (ec50 + sertcov) + b0
  })
}
