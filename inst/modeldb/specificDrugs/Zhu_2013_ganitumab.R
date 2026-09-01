Zhu_2013_ganitumab <- function() {
  description <- "Two-compartment population PK model for intravenous ganitumab (AMG 479, anti-IGF1R IgG1 monoclonal antibody) in adults with metastatic pancreatic cancer or other advanced solid cancers, with pancreatic cancer type, body weight, serum albumin and serum creatinine on disposition (Zhu 2013)"
  reference <- "Zhu M, Gosselin NH, Kuchimanchi M, Johnson J, McCaffery I, Mouksassi MS, Loh E, Lu JF. Differential pharmacokinetics of ganitumab in patients with metastatic pancreatic cancer versus other advanced solid cancers. Clin Pharmacol Drug Dev. 2013;2(4):367-378. doi:10.1002/cpdd.48"
  vignette <- "Zhu_2013_ganitumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Ganitumab was quantified in SERUM by a validated double
  # anti-idiotypic antibody sandwich immunoassay (Zhu 2013 Methods,
  # "Bioanalytical Assays"), not in plasma.
  compartmentData <- list(
    central     = list(analyte = "ganitumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "ganitumab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    TUMTP_PANC = list(
      description        = "Pancreatic cancer type indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (other advanced solid cancers)",
      notes              = paste(
        "The single most significant covariate in the analysis (dMOF = -55.971, P < .001),",
        "entered on BOTH CL and Vc before any other covariate was formally tested",
        "(Zhu 2013 Results, 'Base PK Model Generation'). Zhu 2013 Table 1 reports the two",
        "typical values side by side rather than a multiplier: CL 0.679 L/day non-pancreatic",
        "vs 1.154 L/day pancreatic (1.7-fold), Vc 3.85 L vs 5.13 L (1.3-fold). This model file",
        "encodes the contrast as the exponential coefficients e_panc_cl = log(1.154/0.679) and",
        "e_panc_vc = log(5.13/3.85) applied to the non-pancreatic reference, which reproduces",
        "both published typical values exactly. Once pancreatic cancer type was in the model the",
        "effect of gemcitabine coadministration became non-significant (P > .2); see",
        "covariatesDataExcluded$CONMED_GEMCITABINE.",
        "Cohort split: 37 of 99 (37.4%) pancreatic, 62 of 99 (62.6%) non-pancreatic."
      ),
      source_name        = "cancer type (pancreatic vs non-pancreatic)"
    ),
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power-form effect on CL (exponent 0.984) and on Vc (exponent 0.559); Zhu 2013 Table 1.",
        "Reference 74.7 kg. Zhu 2013 equation (3) defines the reference as 'the median of the",
        "covariate among patients' but never tabulates the cohort medians, and the parenthetical",
        "'(e.g., 70 kg for body weight)' in that sentence is an illustration of what a reference",
        "value IS, not a statement of the value used. The reference weight was therefore recovered",
        "by inverting the paper's own published steady-state exposures (Zhu 2013 Results): at the",
        "median covariates every covariate term equals 1, so AUCss = Dose / CLtypical, and",
        "WTref = AUCss * CLtypical / (mg/kg dose). All four published AUCss values give the same",
        "answer to 3 significant figures -- non-pancreatic 1.32 and 2.20 mg*day/mL at 12 and",
        "20 mg/kg with CL 0.679 L/day, and pancreatic 0.78 and 1.29 mg*day/mL at the same doses",
        "with CL 1.154 L/day, implying 74.69, 74.69, 75.01 and 74.43 kg respectively. 74.7 kg",
        "reproduces all four published values exactly; 70 kg reproduces none of them",
        "(it predicts 1.24 / 2.06 / 0.73 / 1.21). See the vignette Errata."
      ),
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power-form effect on CL only (exponent -0.859); Zhu 2013 Table 1, whose footnote gives",
        "the unit as g/L, so no g/dL conversion is applied. Reference 40 g/L is ASSUMED -- Zhu 2013",
        "defines the reference as the cohort median (equation 3) but never reports the median or",
        "range of albumin, and unlike body weight there is no published quantity that identifies it",
        "(the AUCss inversion pins WT only). 40 g/L is a rounded standard adult serum albumin.",
        "The published exponent, not this centering constant, is the transferable part of the",
        "effect: a user working with their own cohort should re-centre on their own median.",
        "See the vignette Errata."
      ),
      source_name        = "ALB"
    ),
    CREAT = list(
      description        = "Baseline serum creatinine",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power-form effect on CL only (exponent -0.394); Zhu 2013 Table 1, whose footnote gives",
        "the unit as mg/dL. Reference 1 mg/dL is ASSUMED on the same basis as ALB above -- the",
        "cohort median creatinine is not reported anywhere in Zhu 2013 or its supplement.",
        "1 mg/dL is a rounded standard adult serum creatinine. Zhu 2013 Discussion notes that",
        "cachexia in pancreatic cancer lowers creatinine and thereby raises predicted ganitumab CL,",
        "so the true cohort median may sit below 1 mg/dL. Note also that creatinine clearance",
        "derived from WT and CR was deliberately excluded from the covariate screen as collinear",
        "(Zhu 2013 Methods). See the vignette Errata."
      ),
      source_name        = "CR"
    )
  )

  # Covariates Zhu 2013 screened but did NOT retain in the final model.
  # Documentation only -- none of these is referenced in model(). The
  # gemcitabine entry in particular is a headline negative finding of the
  # paper rather than an incidental screen failure.
  #
  # Zhu 2013 screened in TWO stages, and each entry below records which stage
  # rejected it, because the distinction is load-bearing provenance:
  #   Stage 1 (Methods, "Covariate Analysis and Final Model Development") --
  #     graphical exploration plus Pearson correlation against the posthoc CL
  #     and Vc. The full stage-1 list is: age, WT, BMI, sex, race, CR, ALB,
  #     total bilirubin, AST, ALT, ALP, disease stage at enrollment, ECOG
  #     performance status (0/1/2), platelet count, neutrophils, fasting blood
  #     glucose, BUN, LDH, circulating IGF-1, IGFBP3, gemcitabine
  #     coadministration, study, and MDRD-estimated GFR. CrCL derived from WT
  #     and CR was deliberately NOT screened, as collinear with its inputs.
  #   Stage 2 (Results, "Covariate Analysis and Final Model Generation") --
  #     formal stepwise forward addition (P < .01) then backward elimination
  #     (P < .005). Only the stage-1 survivors reached it: on CL, WT, sex, CR,
  #     ALB, AST, ALT, ALP, platelets, neutrophils, IGF-1, IGFBP3, ECOG,
  #     gemcitabine and GFR; on Vc, WT, sex, ALP and gemcitabine.
  # WT, ALB and CR are the only covariates that survived stage 2 (plus cancer
  # type, which was forced in before formal testing); everything else is here.
  #
  # Six screened concepts are recorded in this comment rather than as list
  # members, because no list entry can be written honestly for them:
  #   - platelet count, circulating IGF-1, circulating IGFBP3 (all stage 2, on
  #     CL) and disease stage at enrollment (stage 1) have no entry in the
  #     canonical covariate register. Registering a canonical for a covariate
  #     that no model actually references would put an unused name in the
  #     register, and none of the four carries a published point estimate.
  #   - race (stage 1) is registered only as per-group indicators
  #     (RACE_WHITE / RACE_BLACK / RACE_ASIAN / ...). Zhu 2013 publishes no
  #     baseline-demographics table and never states which race levels were
  #     present, so no RACE_<GROUP> indicator can be named without inventing
  #     the cohort composition.
  #   - study (stage 1) was screened as a categorical, but Zhu 2013 pooled
  #     Studies 1-3 and reports no per-study coefficient, so there is no
  #     STUDY_<id> contrast to record.
  covariatesDataExcluded <- list(
    CONMED_GEMCITABINE = list(
      description        = "Concomitant gemcitabine coadministration",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ganitumab without gemcitabine)",
      notes              = paste(
        "Screened on CL and on Vc and NOT retained -- a headline negative finding of Zhu 2013.",
        "Tested alone, gemcitabine coadministration was significant (dMOF = -34.848, P < .001),",
        "but it was confounded with pancreatic cancer type, which was the stronger effect",
        "(dMOF = -55.971); once cancer type was in the model gemcitabine became non-significant",
        "(P > .2). Corroborated non-compartmentally: mean (SD) CL in non-pancreatic patients was",
        "10.1 (3.8) mL/kg/day without gemcitabine (n = 33, Study 2) vs 11.1 (3.1) mL/kg/day with",
        "gemcitabine (n = 11, Study 3), P > .05; whereas among gemcitabine-treated patients CL was",
        "21.1 (6.4) mL/kg/day in pancreatic vs 11.1 (3.1) in non-pancreatic cancer, P < .001.",
        "46 of 99 patients (46.5%) received concomitant gemcitabine.",
        "No point estimate is published, so no coefficient can be encoded."
      ),
      source_name        = "gemcitabine coadministration"
    ),
    AGE = list(
      description = "Age", units = "years", type = "continuous",
      notes = "Stage 1 only: screened graphically against CL, did not advance to formal stepwise testing. No point estimate published."
    ),
    SEXF = list(
      description = "Sex", units = "(binary)", type = "binary",
      reference_category = "0 (male)",
      notes = "Reached stage 2 (formal stepwise testing) on BOTH CL and Vc and was eliminated from each. No point estimate published."
    ),
    BMI = list(
      description = "Body mass index", units = "kg/m^2", type = "continuous",
      notes = "Stage 1 only: screened graphically against CL, did not advance to formal stepwise testing. Body weight, which did advance and was retained, is the correlated size measure. No point estimate published."
    ),
    TBILI = list(
      description = "Total bilirubin", units = "umol/L", type = "continuous",
      notes = "Stage 1 only: screened graphically against CL, did not advance to formal stepwise testing. Unlike AST, ALT and ALP, bilirubin did not survive stage 1. No point estimate published."
    ),
    AST = list(
      description = "Aspartate aminotransferase", units = "U/L", type = "continuous",
      notes = "Reached stage 2 (formal stepwise testing) on CL and was eliminated. No point estimate published."
    ),
    ALT = list(
      description = "Alanine aminotransferase", units = "U/L", type = "continuous",
      notes = "Reached stage 2 (formal stepwise testing) on CL and was eliminated. No point estimate published."
    ),
    ALP = list(
      description = "Alkaline phosphatase", units = "U/L", type = "continuous",
      notes = "Reached stage 2 (formal stepwise testing) on BOTH CL and Vc -- the only laboratory covariate tested on Vc -- and was eliminated from each. No point estimate published."
    ),
    LDH = list(
      description = "Lactate dehydrogenase", units = "U/L", type = "continuous",
      notes = "Stage 1 only: screened graphically against CL, did not advance to formal stepwise testing. No point estimate published."
    ),
    BUN = list(
      description = "Blood urea nitrogen", units = "mmol/L", type = "continuous",
      notes = "Stage 1 only: screened graphically against CL, did not advance to formal stepwise testing. Note that serum creatinine, the other nitrogenous renal marker screened, WAS retained. No point estimate published."
    ),
    NEUT = list(
      description = "Absolute neutrophil count", units = "10^9/L", type = "continuous",
      notes = "Reached stage 2 (formal stepwise testing) on CL and was eliminated. No point estimate published."
    ),
    ECOG_GE1 = list(
      description = "ECOG performance status at or above 1", units = "(binary)", type = "binary",
      reference_category = "0 (ECOG 0)",
      notes = "Zhu 2013 screened ECOG performance status as a three-level categorical (0, 1, 2); it reached stage 2 (formal stepwise testing) on CL and was eliminated. The binary at-or-above-1 encoding here is the register's canonical form; the paper published no point estimate for any level."
    ),
    CRCL = list(
      description = "MDRD-estimated glomerular filtration rate", units = "mL/min/1.73 m^2", type = "continuous",
      notes = paste(
        "Reached stage 2 (formal stepwise testing) on CL and was eliminated. Zhu 2013 Methods",
        "computed it with the MDRD formula, which the CRCL register entry admits as a source alias.",
        "Its elimination is informative rather than incidental: serum creatinine was retained as a",
        "covariate on CL while the creatinine-based GFR estimate built from it was not, which fits",
        "the paper's cachexia interpretation -- creatinine acts here as a muscle-mass marker rather",
        "than a renal-function marker, and a ~150 kDa IgG1 is not renally cleared. Zhu 2013 also",
        "deliberately excluded WT-and-CR-derived CrCL from screening entirely, as collinear.",
        "No point estimate published."
      )
    ),
    FPG = list(
      description = "Fasting blood glucose", units = "mmol/L", type = "continuous",
      notes = paste(
        "Stage 1 only: screened graphically against CL, did not advance to formal stepwise testing.",
        "Glucose is mechanistically interesting for an IGF1R antagonist (IGF1R blockade causes",
        "hyperglycaemia), but Zhu 2013 screened it as a PK covariate only and reports no estimate.",
        "Units are not stated by Zhu 2013; mmol/L is the SI convention and is recorded here for",
        "completeness only, since no coefficient is encoded."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 99L,
    n_studies      = 3L,
    n_observations = 1227L,
    disease_state  = "Adults with metastatic pancreatic cancer or other advanced solid cancers. 37 of 99 (37.4%) had pancreatic cancer and 62 of 99 (62.6%) had non-pancreatic advanced solid tumours; 46 of 99 (46.5%) received concomitant gemcitabine.",
    dose_range     = "Ganitumab 1-20 mg/kg intravenously once every 2 weeks (Q2W) in all studies (Zhu 2013 Supplemental Table 1 footnote a).",
    notes          = paste(
      "Model-building dataset pooled Studies 1-3: Study 1 (NCT00630552, n = 35, 12 mg/kg,",
      "metastatic pancreatic cancer, ganitumab + gemcitabine), Study 2 (NCT00562380, n = 53,",
      "1/3/10/12/20 mg/kg, advanced solid tumours, ganitumab alone) and Study 3 (NCT00974896,",
      "n = 11, 6 or 12 mg/kg, advanced solid tumours, ganitumab + gemcitabine); see Zhu 2013",
      "Supplemental Table 1. 2% of samples were excluded as below the limit of quantification",
      "(LLOQ 20.2 ng/mL) or in error. Fitted in NONMEM VI level 1.1 with FOCE + INTERACTION.",
      "Study 4 (NCT00626106, n = 104, 12 mg/kg, advanced breast cancer with hormone therapy,",
      "802 concentrations) was held out and used only for external qualification by VPC and is",
      "NOT part of the 99-subject estimation dataset.",
      "IMPORTANT: Zhu 2013 reports no baseline-demographics table -- it states only that",
      "'individual demographic data and disease characteristics for each study were normally",
      "distributed and were similar among studies'. Age, sex, race, and the medians and ranges",
      "of weight, albumin and creatinine are therefore unavailable, which is why the covariate",
      "reference values had to be recovered or assumed; see covariateData notes and the",
      "vignette Errata."
    )
  )

  ini({
    # Structural parameters. Reference subject = NON-PANCREATIC cancer
    # (TUMTP_PANC = 0) at the cohort median covariates (WT 74.7 kg,
    # ALB 40 g/L, CREAT 1 mg/dL), for which every covariate term equals 1.
    lcl <- log(0.679); label("Clearance in non-pancreatic cancer at reference covariates (L/day)")             # Zhu 2013 Table 1: typical CL for non-pancreatic 0.679 L/day (RSE 3.47%)
    lvc <- log(3.85);  label("Central volume of distribution in non-pancreatic cancer at reference covariates (L)") # Zhu 2013 Table 1: typical Vc for non-pancreatic 3.85 L (RSE 2.91%)
    lq  <- log(0.626); label("Intercompartmental clearance (L/day)")                                           # Zhu 2013 Table 1: Q 0.626 L/day (RSE 9.73%)
    lvp <- log(3.41);  label("Peripheral volume of distribution (L)")                                          # Zhu 2013 Table 1: Vp 3.41 L (RSE 8.39%)

    # Continuous covariate effects, power form (COV / COVref)^exponent, per
    # Zhu 2013 equation (3). Neither Q nor Vp carries a covariate effect.
    #
    # Sign note: the ALB and CR exponents are NEGATIVE. Most markdown/text
    # extractions of this Wiley PDF render them as bare magnitudes, because the
    # journal encodes the minus sign as the C0 control byte 0x03 rather than as
    # an ASCII hyphen. Decoding the byte stream directly recovers the signs
    # unambiguously -- `pdftotext -layout ... | cat -A` on Table 1 gives:
    #     WT on CL     0.984 (12.8)      <- no control byte: POSITIVE
    #     ALB on CL   ^C0.859 (22.2)     <- ^C = 0x03 = minus: NEGATIVE
    #     CR on CL    ^C0.394 (24.7)     <- ^C = 0x03 = minus: NEGATIVE
    # The within-table contrast is what makes this conclusive: WT carries no
    # control byte while ALB and CR both do, so the missing signs are not a
    # blanket extraction artifact applied to every row. The same encoding
    # appears in the Results prose (dMOF "^C1255.55", "^C55.971", "^C34.848").
    # Corroborated independently by the paper's own prose in the same Results
    # paragraph: "higher baseline albumin and creatinine values were associated
    # with decreased ganitumab CL", and again in the Discussion: "baseline
    # albumin was negatively correlated with CL for all cancer types" and
    # "lower muscle mass ... may lead to lower creatinine levels and higher CL".
    e_wt_cl    <-  0.984; label("Power exponent of body weight on CL (unitless)")        # Zhu 2013 Table 1: WT on CL 0.984 (RSE 12.8%)
    e_alb_cl   <- -0.859; label("Power exponent of serum albumin on CL (unitless)")      # Zhu 2013 Table 1: ALB on CL -0.859 (RSE 22.2%)
    e_creat_cl <- -0.394; label("Power exponent of serum creatinine on CL (unitless)")   # Zhu 2013 Table 1: CR on CL -0.394 (RSE 24.7%)
    e_wt_vc    <-  0.559; label("Power exponent of body weight on Vc (unitless)")        # Zhu 2013 Table 1: WT on Vc 0.559 (RSE 21.6%)

    # Pancreatic cancer type, exponential form exp(coefficient * TUMTP_PANC).
    # Zhu 2013 Table 1 tabulates the two typical values rather than a ratio;
    # the coefficients below are the logs of the published ratios, so that
    # exp(lcl + e_panc_cl) = 1.154 L/day and exp(lvc + e_panc_vc) = 5.13 L
    # reproduce the published pancreatic typical values exactly.
    e_panc_cl <- log(1.154 / 0.679); label("Exponential coefficient of pancreatic cancer type on CL (unitless)") # Zhu 2013 Table 1: typical CL 1.154 (RSE 5.84%) vs 0.679 L/day => 1.7-fold
    e_panc_vc <- log(5.13 / 3.85);   label("Exponential coefficient of pancreatic cancer type on Vc (unitless)") # Zhu 2013 Table 1: typical Vc 5.13 (RSE 4.21%) vs 3.85 L => 1.3-fold

    # IIV. Zhu 2013 equation (1) is theta_i = theta_typical * exp(eta_i) with
    # "estimated variance v^2 (hereafter ... v is referred to as OMEGA)", and
    # Table 1 reports each IIV as a percentage. Under the paper's own notation
    # the tabulated percentage is OMEGA = v, i.e. the SD on the log scale, so
    # variance = (percent / 100)^2:
    #   CL 26.2% -> 0.0686,  Vc 21.7% -> 0.0471,
    #   Q  43.7% -> 0.1910,  Vp 60.2% -> 0.3624.
    # That reading is corroborated by the residual error, where Table 1's
    # "Proportional error (%) 16.4" matches the Results text "RV decreased from
    # 31.5% to 16.5%" -- an SD-scale percentage, not a variance.
    # Diagonal, per Zhu 2013 Methods: "IIV was implemented as a diagonal
    # variance matrix for random effects." IOV was evaluated and NOT included.
    etalcl ~ 0.068644  # Zhu 2013 Table 1: IIV CL 26.2% (RSE 18.8%); 0.262^2
    etalvc ~ 0.047089  # Zhu 2013 Table 1: IIV Vc 21.7% (RSE 17.7%); 0.217^2
    etalq  ~ 0.190969  # Zhu 2013 Table 1: IIV Q 43.7% (RSE 38.7%); 0.437^2
    etalvp ~ 0.362404  # Zhu 2013 Table 1: IIV Vp 60.2% (RSE 25.7%); 0.602^2

    # Residual error. Zhu 2013 equation (2) allows a combined additive plus
    # proportional error, but the final model is proportional only: the base
    # model "had a linear, time-invariant, 2-compartment structure with a
    # proportional error model" and Table 1 reports no additive term.
    propSd <- 0.164; label("Proportional residual error (fraction)")  # Zhu 2013 Table 1: proportional error 16.4% (RSE 9.26%)
  })
  model({
    # Individual disposition parameters. Continuous covariates enter as
    # (COV / COVref)^exponent (Zhu 2013 equation 3); pancreatic cancer type
    # enters exponentially. Reference covariates: WT 74.7 kg, ALB 40 g/L,
    # CREAT 1 mg/dL, TUMTP_PANC = 0.
    cl <- exp(lcl + etalcl) *
      (WT    / 74.7)^e_wt_cl *
      (ALB   / 40)^e_alb_cl *
      (CREAT / 1)^e_creat_cl *
      exp(e_panc_cl * TUMTP_PANC)

    vc <- exp(lvc + etalvc) *
      (WT / 74.7)^e_wt_vc *
      exp(e_panc_vc * TUMTP_PANC)

    q  <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    # Two-compartment micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # Dose in mg, volumes in L => central / vc has units mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
