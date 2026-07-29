Zandvliet_2016_corifollitropin_alfa <- function() {
  description <- paste(
    "One-compartment subcutaneous population pharmacokinetic model for",
    "corifollitropin alfa (a long-acting recombinant gonadotrophin) in",
    "women undergoing controlled ovarian stimulation (Zandvliet 2016).",
    "Pooled analysis of 2557 evaluable women from five phase II and III",
    "trials (single SC doses of 60-180 ug). Corifollitropin alfa is",
    "absorbed first-order from a subcutaneous depot into a one-",
    "compartment central pool with first-order elimination; body weight",
    "is the major covariate, with allometric (WT/60)^exponent power",
    "effects on apparent clearance and apparent volume. Apparent",
    "bioavailability is modulated by body-mass index, race (Asian and",
    "Black indicators vs Caucasian reference), and remains anchored at",
    "F = 1. The model jointly describes total FSH immunoreactivity by",
    "adding an endogenous follicle stimulating hormone (FSH)",
    "compartment whose pre-dose steady-state baseline is FSHbaseline and",
    "whose synthesis is set to zero from corifollitropin administration",
    "onwards (per the paper's structural model). Total FSH",
    "immunoreactivity (IU/L) equals SCALE * corifollitropin concentration",
    "(ng/mL) plus the endogenous FSH compartment value, where SCALE is",
    "fixed at 6.11 IU/L per ng/mL from an upstream analysis. Trial-",
    "specific multiplicative effects on the FSH immunoreactivity",
    "prediction (1.26 for trial 06029 and 1.12 for trial 38825) are",
    "exposed via binary study indicators that default to zero for",
    "general simulation use."
  )
  reference <- paste(
    "Zandvliet AS, Prohn M, de Greef R, van Aarle F, McCrary Sisk C,",
    "Stegmann BJ. Impact of patient characteristics on the",
    "pharmacokinetics of corifollitropin alfa during controlled ovarian",
    "stimulation. Br J Clin Pharmacol. 2016;82(1):74-82.",
    "doi:10.1111/bcp.12939."
  )
  vignette <- "Zandvliet_2016_corifollitropin_alfa"
  paper_specific_compartments <- c("endo_fsh")
  units <- list(time = "hour", dosing = "ug", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Subject baseline body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Allometric power scaling on apparent CL",
        "((WT/60)^1.20) and apparent V ((WT/60)^1.23) with reference",
        "weight 60 kg (Zandvliet 2016 Table 3 theta_12 and theta_13).",
        "Also a power exponent (WT/60)^-0.832 on the endogenous-FSH",
        "elimination rate KeFSH (Zandvliet 2016 Table 3 theta_22).",
        "Cohort medians by trial range 54.0 to 68.8 kg",
        "(Table 2)."
      ),
      source_name        = "WT"
    ),
    BMI = list(
      description        = "Subject baseline body mass index.",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Power scaling (BMI/23.3)^-0.245 on",
        "apparent bioavailability F",
        "(Zandvliet 2016 Table 3 theta_16; reference 23.3 kg/m^2 from",
        "the cohort and the implicit normalisation of the power form).",
        "Higher BMI -> lower F. The paper text frames the same effect",
        "as ~14% higher dose-normalised exposure at BMI 18 vs 32",
        "(at matched body weight)."
      ),
      source_name        = "BMI"
    ),
    AGE = list(
      description        = "Subject age in years.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Power scaling (AGE/34)^-0.815 on the",
        "endogenous-FSH elimination rate KeFSH and",
        "(AGE/34)^0.423 on FSHbaseline",
        "(Zandvliet 2016 Table 3 theta_20 and theta_19;",
        "reference age 34 years from the implicit normalisation of the",
        "power form). AGE was not retained as a covariate on",
        "corifollitropin alfa CL, V, or F."
      ),
      source_name        = "AGE"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator (1 = Asian, 0 = otherwise).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian; Caucasian is the typical-value reference category in the source model)",
      notes              = paste(
        "Time-fixed per subject. Two effects:",
        "(i) multiplicative factor 0.843 on apparent bioavailability F",
        "for Asian vs non-Asian",
        "(Zandvliet 2016 Table 3 theta_17; 15.7% lower exposure in",
        "Asian subjects per Results / Race section);",
        "(ii) multiplicative factor 0.697 on the endogenous-FSH",
        "elimination rate KeFSH for Asian vs non-Asian",
        "(Zandvliet 2016 Table 3 theta_21).",
        "Cohort fractions ranged from 1.3% (phase II) to 45.1%",
        "(trial 107012 in Korea and Taiwan)."
      ),
      source_name        = "ASIAN"
    ),
    RACE_BLACK = list(
      description        = "Black race indicator (1 = Black, 0 = otherwise).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Black; Caucasian is the typical-value reference category in the source model)",
      notes              = paste(
        "Time-fixed per subject. Multiplicative factor 1.13 on apparent",
        "bioavailability F for Black vs non-Black",
        "(Zandvliet 2016 Table 3 theta_18; 13% higher exposure in Black",
        "subjects per Results / Race section).",
        "Cohort fractions ranged from 0.4% to 10%."
      ),
      source_name        = "BLACK"
    ),
    STUDY_06029 = list(
      description        = paste(
        "Trial 06029 (NCT01144416) indicator.",
        "1 = subject is from trial 06029 of the Zandvliet 2016",
        "integrated population PK analysis (phase III multiple-cycles",
        "non-inferiority trial; single 150 ug SC dose per cycle);",
        "0 = otherwise."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any non-06029 trial)",
      notes              = paste(
        "Subject-level / time-fixed. Multiplicative factor 1.26 on the",
        "total FSH immunoreactivity prediction in trial 06029",
        "(Zandvliet 2016 Table 3 theta_14). The trial effect was",
        "introduced to account for unexplained between-trial",
        "differences in the FSH immunoreactivity assay; it leaves",
        "corifollitropin alfa concentrations unaffected.",
        "Default 0 for general simulation use; set to 1 only when",
        "replicating Zandvliet 2016 trial 06029."
      ),
      source_name        = "STUD06029"
    ),
    STUDY_38825 = list(
      description        = paste(
        "Trial 38825 (NCT00696878) indicator.",
        "1 = subject is from trial 38825 of the Zandvliet 2016",
        "integrated population PK analysis (phase III open-label",
        "uncontrolled repeated-cycle trial; up to 3 cycles of 150 ug",
        "SC corifollitropin alfa); 0 = otherwise."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any non-38825 trial)",
      notes              = paste(
        "Subject-level / time-fixed. Multiplicative factor 1.12 on the",
        "total FSH immunoreactivity prediction in trial 38825",
        "(Zandvliet 2016 Table 3 theta_15). The trial effect was",
        "introduced to account for unexplained between-trial",
        "differences in the FSH immunoreactivity assay; it leaves",
        "corifollitropin alfa concentrations unaffected.",
        "Default 0 for general simulation use; set to 1 only when",
        "replicating Zandvliet 2016 trial 38825."
      ),
      source_name        = "STUD38825"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 2557L,
    n_studies       = 5L,
    age_range       = "phase II / III subjects undergoing controlled ovarian stimulation; cohort means 30.9 to 38.0 years (Table 2)",
    weight_range    = "cohort means 54.0 to 68.8 kg; population range broadly 50-90 kg per Results / Body weight section (Table 2)",
    sex_female_pct  = 100,
    race_ethnicity  = c(Caucasian = 84.5, Asian = 8.0, Black = 4.8, Other = 2.7),
    disease_state   = "Subfertility undergoing controlled ovarian stimulation in a GnRH antagonist protocol; switched to daily recombinant FSH (rFSH) from day 8 onwards.",
    dose_range      = "Single SC injection of 60, 100, 120, 150, or 180 ug corifollitropin alfa (Table 1). Dose 100 ug for body weight <= 60 kg and 150 ug for body weight > 60 kg in women aged <= 36 years; 150 ug irrespective of body weight in women > 36 years.",
    regions         = "Multinational (Europe, North America, Asia including Korea and Taiwan for trial 107012).",
    n_observations  = "Serum corifollitropin alfa concentrations (DVID = 3) and total FSH immunoreactivity levels (DVID = 4) measured at trial-specific visit schedules (Table 1).",
    notes           = paste(
      "Demographics from Zandvliet 2016 Table 2 (per-trial dose groups).",
      "Race-fraction aggregates above are coarse weighted averages",
      "across trials and dose groups; per-trial distributions vary",
      "(e.g., 45.1% Asian in trial 107012 vs 1-2% elsewhere). The",
      "scaling factor SCALE = 6.11 IU/L per ng/mL converting",
      "corifollitropin alfa concentration into FSH immunoreactivity",
      "was fixed from an upstream analysis (Table 3 theta_7 footnote).",
      "The covariance between CL/F and V/F random effects was",
      "estimated at 0.056 (95% CI 0.051, 0.061)."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # CORIFOLLITROPIN ALFA STRUCTURAL PARAMETERS
    # ------------------------------------------------------------------
    # Reference values 60 kg (body weight) and 23.3 kg/m^2 (BMI) are
    # the implicit normalisations in the source paper's power-form
    # covariate equations (Zandvliet 2016 Table 3 covariate-
    # relationship column).

    lka <- log(0.0436)
    label("Corifollitropin alfa absorption rate ka (1/h)")              # Zandvliet 2016 Table 3 theta_1

    lvc <- log(19.1)
    label("Corifollitropin alfa apparent central volume V/F (L)")       # Zandvliet 2016 Table 3 theta_2

    lcl <- log(0.19)
    label("Corifollitropin alfa apparent clearance CL/F (L/h)")         # Zandvliet 2016 Table 3 theta_3

    lfdepot <- fixed(log(1))
    label("Corifollitropin alfa bioavailability anchor F (FIXED at 1; the dose / V / CL scaling carries the apparent F)")  # Zandvliet 2016 Table 3 theta_6 = 1 FIX

    # ------------------------------------------------------------------
    # CORIFOLLITROPIN ALFA COVARIATE EFFECTS
    # ------------------------------------------------------------------

    e_wt_cl <- 1.20
    label("Power exponent of (WT/60) on apparent CL (unitless)")        # Zandvliet 2016 Table 3 theta_12

    e_wt_vc <- 1.23
    label("Power exponent of (WT/60) on apparent V (unitless)")         # Zandvliet 2016 Table 3 theta_13

    e_bmi_fdepot <- -0.245
    label("Power exponent of (BMI/23.3) on apparent F (unitless)")      # Zandvliet 2016 Table 3 theta_16 (negative: higher BMI -> lower exposure per Results / Body weight and BMI)

    e_race_asian_fdepot <- 0.843
    label("Multiplicative factor on apparent F for Asian vs non-Asian (unitless)")  # Zandvliet 2016 Table 3 theta_17

    e_race_black_fdepot <- 1.13
    label("Multiplicative factor on apparent F for Black vs non-Black (unitless)")  # Zandvliet 2016 Table 3 theta_18

    # ------------------------------------------------------------------
    # ENDOGENOUS FSH SUBMODEL STRUCTURAL PARAMETERS
    # ------------------------------------------------------------------
    # Pre-dose endogenous FSH is at steady state (zero-order
    # synthesis balanced against first-order elimination). After
    # corifollitropin administration, synthesis is set to zero
    # (Zandvliet 2016 Structural model paragraph). The model
    # therefore needs only the steady-state baseline FSHbaseline
    # (as the initial condition) and the first-order elimination
    # rate KeFSH.

    lkout <- log(0.0101)
    label("Endogenous FSH first-order elimination rate KeFSH (1/h)")    # Zandvliet 2016 Table 3 theta_4

    lrbase <- log(6.85)
    label("Endogenous FSH steady-state baseline FSHbaseline (IU/L)")    # Zandvliet 2016 Table 3 theta_5

    # ------------------------------------------------------------------
    # ENDOGENOUS FSH SUBMODEL COVARIATE EFFECTS
    # ------------------------------------------------------------------

    e_age_kout <- -0.815
    label("Power exponent of (AGE/34) on KeFSH (unitless)")             # Zandvliet 2016 Table 3 theta_20 (95% CI -1.206, -0.425; negative sign confirmed by CI bracket)

    e_race_asian_kout <- 0.697
    label("Multiplicative factor on KeFSH for Asian vs non-Asian (unitless)")  # Zandvliet 2016 Table 3 theta_21

    e_wt_kout <- -0.832
    label("Power exponent of (WT/60) on KeFSH (unitless)")              # Zandvliet 2016 Table 3 theta_22 (95% CI -1.37, -0.23; negative sign confirmed by CI bracket)

    e_age_rbase <- 0.423
    label("Power exponent of (AGE/34) on FSHbaseline (unitless)")       # Zandvliet 2016 Table 3 theta_19

    # ------------------------------------------------------------------
    # TRIAL EFFECTS ON FSH IMMUNOREACTIVITY PREDICTION
    # ------------------------------------------------------------------
    # The published model retained two trial-specific multiplicative
    # effects on the FSH immunoreactivity observation to account for
    # unexplained between-trial differences (Zandvliet 2016 Results /
    # Final model). For general simulation, both indicators default
    # to 0 and the trial effects are off.

    e_study_06029_fshir <- 1.26
    label("Multiplicative factor on total FSH immunoreactivity for trial 06029 (unitless)")  # Zandvliet 2016 Table 3 theta_14

    e_study_38825_fshir <- 1.12
    label("Multiplicative factor on total FSH immunoreactivity for trial 38825 (unitless)")  # Zandvliet 2016 Table 3 theta_15

    # ------------------------------------------------------------------
    # INTER-INDIVIDUAL VARIABILITY
    # ------------------------------------------------------------------
    # IIV is reported as CV% in Table 3 with the conversion
    # omega^2 = log(CV^2 + 1) (Table 3 footnote section sign).
    # Source CV% values:
    #   Ka          29.2%   -> omega^2 = log(0.292^2 + 1) = 0.08191
    #   V/F         29.3%   -> omega^2 = log(0.293^2 + 1) = 0.08243
    #   CL/F        22.1%   -> omega^2 = log(0.221^2 + 1) = 0.04764
    #   cov(CL,V)   0.056             (Table 3 footnote *)
    #   KeFSH       29%     -> omega^2 = log(0.29^2 + 1) = 0.08074
    #                          (see vignette Assumptions and deviations:
    #                          the published Table 3 IIV cell prints
    #                          "29%" with a 95% CI of 51.4-69.7% that
    #                          does not bracket the point estimate;
    #                          taken at face value here)
    #   FSHbaseline 26.7%   -> omega^2 = log(0.267^2 + 1) = 0.06877

    etalcl + etalvc ~ c(0.04764, 0.056, 0.08243)
    etalka          ~ 0.08191
    etalkout        ~ 0.08074
    etalrbase       ~ 0.06877

    # ------------------------------------------------------------------
    # RESIDUAL ERROR
    # ------------------------------------------------------------------
    # Corifollitropin alfa observation (Cc, ng/mL): proportional
    # 13.10% CV with additive component fixed at 0 (Zandvliet 2016
    # Table 3 theta_8 and theta_9).
    # FSH immunoreactivity observation (FSHir, IU/L): proportional
    # 4.86% CV with additive component 0.717 IU/L (Table 3 theta_10
    # and theta_11).

    propSd <- 0.131
    label("Proportional residual SD on corifollitropin alfa concentration (fraction)")  # Zandvliet 2016 Table 3 theta_8

    propSd_FSHir <- 0.0486
    label("Proportional residual SD on total FSH immunoreactivity (fraction)")          # Zandvliet 2016 Table 3 theta_10

    addSd_FSHir <- 0.717
    label("Additive residual SD on total FSH immunoreactivity (IU/L)")                  # Zandvliet 2016 Table 3 theta_11
  })

  model({
    # Reference values for the power-form covariate equations (Zandvliet
    # 2016 Table 3 covariate-relationship column; reference 60 kg, 23.3
    # kg/m^2, and 34 years are the implicit normalisations).
    ref_wt  <- 60     # kg
    ref_bmi <- 23.3   # kg/m^2
    ref_age <- 34     # years

    # Assay scaling factor converting corifollitropin alfa concentration
    # (ng/mL) into total FSH immunoreactivity (IU/L). Fixed at the
    # upstream-derived value 6.11 IU/L per ng/mL (Zandvliet 2016 Table 3
    # theta_7 footnote ** "scaling factor fixed to 6.11 as derived in a
    # previous analysis").
    scale_fsh <- 6.11

    # ------------------------------------------------------------------
    # Corifollitropin alfa individual parameters
    # ------------------------------------------------------------------
    cl <- exp(lcl + etalcl) * (WT / ref_wt)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / ref_wt)^e_wt_vc
    ka <- exp(lka + etalka)

    # Apparent bioavailability: anchored at exp(lfdepot) = 1; multiplied
    # by power-of-BMI factor and indicator-power factors for race. The
    # power-on-indicator forms theta_17^ASIAN and theta_18^BLACK reduce
    # to the multiplicative factor when the indicator equals 1 and to 1
    # otherwise (Zandvliet 2016 Table 3 theta_6 row).
    fdepot <- exp(lfdepot) *
              (BMI / ref_bmi)^e_bmi_fdepot *
              e_race_asian_fdepot^RACE_ASIAN *
              e_race_black_fdepot^RACE_BLACK

    # ------------------------------------------------------------------
    # Endogenous FSH submodel individual parameters
    # ------------------------------------------------------------------
    kout <- exp(lkout + etalkout) *
            (AGE / ref_age)^e_age_kout *
            e_race_asian_kout^RACE_ASIAN *
            (WT  / ref_wt)^e_wt_kout

    rbase <- exp(lrbase + etalrbase) *
             (AGE / ref_age)^e_age_rbase

    # ------------------------------------------------------------------
    # Micro-constants
    # ------------------------------------------------------------------
    kel <- cl / vc

    # ------------------------------------------------------------------
    # ODE system
    # Corifollitropin alfa: first-order absorption from depot into a
    # one-compartment central pool with first-order elimination
    # (Zandvliet 2016 Structural model paragraph and Figure 3).
    # Endogenous FSH: pre-dose steady state at rbase; zero-order
    # synthesis is set to zero from corifollitropin administration
    # onwards, leaving only first-order elimination from the
    # baseline value.
    # ------------------------------------------------------------------
    d/dt(depot)   <- -ka   * depot
    d/dt(central) <-  ka   * depot - kel * central
    d/dt(endo_fsh)    <- -kout * endo_fsh

    # Initial condition for the endogenous FSH compartment: pre-dose
    # steady-state baseline rbase (Zandvliet 2016 Structural model).
    endo_fsh(0) <- rbase

    # Bioavailability applied to the depot compartment.
    f(depot) <- fdepot

    # ------------------------------------------------------------------
    # Observations
    # ------------------------------------------------------------------
    # Corifollitropin alfa concentration (ng/mL).
    # Dose units = ug; vc in L gives ug/L = ng/mL directly when central
    # is in dose units of ug.
    Cc <- central / vc

    # Total FSH immunoreactivity (IU/L). Trial-specific multiplicative
    # effects modify only the FSH immunoreactivity prediction; the
    # power-on-indicator form leaves FSHir unchanged when the indicator
    # is 0 (Zandvliet 2016 Final model paragraph).
    FSHir <- (scale_fsh * Cc + endo_fsh) *
             e_study_06029_fshir^STUDY_06029 *
             e_study_38825_fshir^STUDY_38825

    Cc    ~ prop(propSd)
    FSHir ~ add(addSd_FSHir) + prop(propSd_FSHir)
  })
}
