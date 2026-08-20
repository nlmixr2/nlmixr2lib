Majid_2024_lenvatinib <- function() {
  description <- "Three-compartment population PK model for lenvatinib pooled across healthy subjects and patients with differentiated thyroid cancer (DTC), renal cell carcinoma (RCC), hepatocellular carcinoma (HCC), and other solid tumors (Majid 2024). Zero-order input into the absorption depot over a duration D1 followed by first-order absorption (ka) into the central compartment, linear elimination, and covariate effects of body weight (allometric exponent 0.75 on CL/F and Q/F, 1.0 on the volumes), concomitant CYP3A4 inhibitors (-10.4 percent on CL/F), serum albumin < 30 g/L (-10.0 percent), alkaline phosphatase > ULN (-9.0 percent), healthy-subject cohort (+19 percent), and tumor type (DTC -4.9 percent, HCC -13.8 percent, RCC -14.9 percent versus the other-solid-tumor reference), plus a capsule-versus-tablet relative bioavailability of 0.882. This is the updated successor to the Gupta 2016 lenvatinib popPK model; see modellib('Gupta_2016_lenvatinib')."
  reference <- paste(
    "Majid O, Hayato S, Sreerama Reddy SH, Lalovic B, Hihara T, Hoshi T,",
    "Funahashi Y, Aluri J, Takenaka O, Yasuda S, Hussein Z.",
    "Population pharmacokinetic-pharmacodynamic modeling of serum biomarkers as",
    "predictors of tumor dynamics following lenvatinib treatment in patients with",
    "radioiodine-refractory differentiated thyroid cancer (RR-DTC).",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(6):954-969.",
    "doi:10.1002/psp4.13130.",
    "Structural model updated from Gupta A, et al. Br J Clin Pharmacol.",
    "2016;81(6):1124-1133; see modellib('Gupta_2016_lenvatinib').",
    sep = " "
  )
  vignette <- "Majid_2024_lenvatinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric exponent fixed at 0.75 on CL/F, Q1/F and Q2/F and at 1.0 on V1/F, V2/F and V3/F, both normalised to the 73.2 kg PK-population median (Majid 2024 Table 1 footnote equations; Table S2 PK-population median weight 73.2 kg). Supplement Text S1 confirms THETA(12) = 0.75 FIX and THETA(13) = 1 FIX.",
      source_name        = "WGT"
    ),
    CONMED_CYP3A4_INH = list(
      description        = "Concomitant CYP3A4 inhibitor coadministration indicator (1 = concomitant CYP3A4 inhibitor; 0 = none).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant CYP3A4 inhibitor).",
      notes              = "Multiplicative power-form effect on CL/F: 0.896^CONMED_CYP3A4_INH (-10.4 percent when 1), Majid 2024 Table 1. 37 of 1921 PK-population subjects were positive (Table S2). Concomitant CYP3A4 inducers were carried in the Gupta 2016 predecessor model but are NOT in the Majid 2024 final model (only 6 of 1921 subjects were exposed; Table 1 lists no inducer term and Text S1 TVCL2 has no inducer factor).",
      source_name        = "INHIB"
    ),
    ALB = list(
      description        = "Serum albumin concentration.",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Binarized inline as alb_low <- (ALB < 30) per Majid 2024 Table 1 note ('ALB, 0 (>= ALB 30 g/L) or 1 (< ALB 30 g/L)') and Text S1 ('IF(ALB.LT.30) ALBC=1'). Multiplicative power-form effect on CL/F: 0.900^alb_low (-10.0 percent when 1). The Table 1 row is printed as 'Albumin (> ULN) ~ CL/F' but the table note and the Discussion ('CL/F decreased by 10.0 percent with albumin levels below 30 g/L') both make the < 30 g/L direction unambiguous.",
      source_name        = "ALB"
    ),
    ALP = list(
      description        = "Serum alkaline phosphatase activity.",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Binarized inline as alp_high <- (ALP > alp_uln) per Majid 2024 Table 1 note ('ALP, Alkaline phosphatase (IU/L) 0 (<= upper limit of normal) or 1 (> upper limit of normal value)'); Text S1 codes it from the ALP/ULN ratio column ALPR ('IF(ALPR.GT.1) ALPC=1'). Multiplicative power-form effect on CL/F: 0.910^alp_high (-9.0 percent when 1). The 120 U/L threshold used inline is a representative adult ULN and matches the choice made in Gupta_2016_lenvatinib.R; downstream users whose data carry an ALP/ULN ratio should scale the threshold to their own laboratory ULN.",
      source_name        = "ALP"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-subject cohort indicator (1 = healthy volunteer from a phase 1 clinical pharmacology study; 0 = cancer patient).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (cancer patient).",
      notes              = "Multiplicative power-form effect on CL/F: 1.19^DIS_HEALTHY (+19 percent when 1), Majid 2024 Table 1. Text S1 sets HV = 1 for studies numbered below 100. 199 of 1921 PK-population subjects were healthy volunteers (Table S2).",
      source_name        = "HV"
    ),
    TUMTP_DTC = list(
      description        = "Differentiated thyroid cancer tumor-type indicator (1 = DTC; 0 = other tumor type or healthy subject).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the reference stratum is cancer patients with a solid tumor other than DTC, RCC or HCC; Figure S2 legend).",
      notes              = "Multiplicative power-form effect on CL/F: 0.951^TUMTP_DTC (-4.9 percent when 1), Majid 2024 Table 1. Text S1 sets DTC = 1 for studies 303, 211 and 201. DTC, RCC, HCC and DIS_HEALTHY are mutually exclusive in the source dataset; setting all four to 0 recovers the other-solid-tumor reference.",
      source_name        = "DTC"
    ),
    TUMTP_HCC = list(
      description        = "Hepatocellular carcinoma tumor-type indicator (1 = HCC; 0 = other tumor type or healthy subject).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (cancer patients with a solid tumor other than DTC, RCC or HCC).",
      notes              = "Multiplicative power-form effect on CL/F: 0.862^TUMTP_HCC (-13.8 percent when 1), Majid 2024 Table 1. Text S1 sets HCC = 1 for studies 202 and 304.",
      source_name        = "HCC"
    ),
    TUMTP_RCC = list(
      description        = "Renal cell carcinoma tumor-type indicator (1 = RCC; 0 = other tumor type or healthy subject).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (cancer patients with a solid tumor other than DTC, RCC or HCC).",
      notes              = "Multiplicative power-form effect on CL/F: 0.851^TUMTP_RCC (-14.9 percent when 1), Majid 2024 Table 1. Text S1 sets RCC = 1 for studies 205, 218, 221 and 112. Concomitant everolimus (given with lenvatinib in the RCC studies) was tested and was NOT statistically significant on lenvatinib PK (Majid 2024 Results and Discussion), so no everolimus term appears in the final model.",
      source_name        = "RCC"
    ),
    FORM_CAPSULE = list(
      description        = "Capsule versus tablet formulation indicator (1 = capsule; 0 = tablet).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (tablet; relative bioavailability fixed at 1 for the tablet reference arm).",
      notes              = "Relative bioavailability F1 = 0.882^FORM_CAPSULE (Majid 2024 Table 1 footnote equation 'F1 = 1 * 0.882^FORM'; Text S1 sets FM = 1 when FORM < 4 and F1 = THETA(9)^FM). Unlike Gupta 2016, Majid 2024 estimates no inter-individual variability on F1. 1561 of 1921 PK-population subjects received capsules (Table S2).",
      source_name        = "FORM"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte  = "lenvatinib",
      units    = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte  = "lenvatinib",
      units    = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte  = "lenvatinib",
      units    = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral2 = list(
      analyte  = "lenvatinib",
      units    = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1921L,
    n_studies      = 19L,
    age_range      = "18-92 years",
    age_median     = "61.0 years",
    weight_range   = "32.6-190 kg",
    weight_median  = "73.2 kg",
    sex_female_pct = 32.2,
    race_ethnicity = c(
      White                          = 61.6,
      Japanese                       = 10.8,
      Asian_other_than_Chinese_Japanese = 9.0,
      Chinese                        = 8.8,
      Black_African_American         = 4.2,
      Other                          = 3.9,
      Missing                        = 1.7
    ),
    disease_state  = "Pooled cohort of 1921 subjects: HCC n = 534, DTC n = 542, RCC n = 491, other solid tumor n = 155, healthy subjects n = 199 (Table S2).",
    dose_range     = "3.2-32 mg oral lenvatinib as tablet or capsule; once-daily and twice-daily regimens across phase 1-3 studies. Subjects with a starting dose below 3.2 mg were excluded from the PK dataset (Table S1 footnote; Text S1 IGNORE=(DGRP.LT.3.2)).",
    regions        = "Multiregional; 19 pooled studies (Table S1).",
    n_observations = "17550 lenvatinib plasma concentrations from 1921 subjects across 19 studies.",
    lab_summary    = "Albumin 40.2 (4.9) g/L, median 40, range 19-58; ALP 121.5 (103.7) IU/L, median 88, range 19-1135; ALT 28.2 (27.4) IU/L; AST 33.3 (35.5) IU/L; bilirubin 9.8 (6.7) umol/L; creatinine clearance 93.2 (35.2) mL/min, median 89, range 17-305 (Table S2).",
    co_medication  = "Concomitant everolimus in 442 of 1921 subjects (RCC studies); CYP3A4 inhibitors in 37; CYP3A4 inducers in 6 (Table S2).",
    ecog_distribution = "ECOG 0: 1029; 1: 601; 2: 24; 3: 1; missing: 265 (Table S2).",
    notes          = "Demographics reproduced from Majid 2024 Table S2 (PK population, N = 1921). Sex 1302 male / 619 female gives 32.2 percent female. Race percentages computed from the Table S2 counts (White 1183, Japanese 207, Asian other than Chinese and Japanese 173, Chinese 170, Black/Afro-American 81, Others 75, Missing 32) over N = 1921. The PK/PD biomarker and tumor sub-populations are different, smaller cohorts; see Majid_2024_lenvatinib_biomarkers.R and Majid_2024_lenvatinib_tumor.R."
  )

  ini({
    # ---- Structural parameters: typical values at the 73.2 kg reference
    # weight, in a cancer patient with a solid tumor other than DTC / RCC /
    # HCC, no concomitant CYP3A4 inhibitor, albumin >= 30 g/L, ALP <= ULN,
    # and the tablet reference formulation (F1 = 1).
    lcl        <- log(6.28);   label("Apparent clearance CL/F at 73.2 kg (L/h)")                                 # Majid 2024 Table 1: CL/F = 6.28 L/h (%RSE 2.88; bootstrap 6.33 [5.93-6.63])
    lvc        <- log(46.0);   label("Apparent central volume V1/F at 73.2 kg (L)")                              # Majid 2024 Table 1: V1/F = 46.0 L (%RSE 2.61)
    lvp        <- log(28.3);   label("Apparent first peripheral volume V2/F at 73.2 kg (L)")                     # Majid 2024 Table 1: V2/F = 28.3 L (%RSE 5.90)
    lvp2       <- log(30.9);   label("Apparent second peripheral volume V3/F at 73.2 kg (L)")                    # Majid 2024 Table 1: V3/F = 30.9 L (%RSE 3.09)
    lq         <- log(3.57);   label("Apparent inter-compartmental clearance Q1/F (V1/F <-> V2/F) at 73.2 kg (L/h)")  # Majid 2024 Table 1: Q1 = 3.57 L/h (%RSE 2.89)
    lq2        <- log(0.688);  label("Apparent inter-compartmental clearance Q2/F (V1/F <-> V3/F) at 73.2 kg (L/h)")  # Majid 2024 Table 1: Q2 = 0.688 L/h (%RSE 3.15)
    lka        <- log(0.803);  label("First-order absorption rate constant ka (1/h)")                            # Majid 2024 Table 1: Ka = 0.803 1/h (%RSE 4.36)
    lduration  <- log(1.27);   label("Duration of the zero-order input into the depot, D1 (h)")                  # Majid 2024 Table 1: D1 = 1.27 h (%RSE 3.61)
    lfcap      <- log(0.882);  label("Relative bioavailability of capsule versus tablet (unitless)")             # Majid 2024 Table 1: F1 = 0.882 (%RSE 1.10)

    # ---- Allometric exponents. Both were held FIXED in the source control
    # stream (Text S1 $THETA slots 12 and 13 carry the FIX flag).
    e_wt_cl    <- fixed(0.75); label("Allometric exponent of body weight on CL/F, Q1/F and Q2/F (unitless)")     # Majid 2024 Text S1 $THETA 12: '0.75 FIX ;[WGT_CL]'; Table 1 footnote CL/F = 6.28 * (BW/73.2)^0.75
    allo_v     <- fixed(1.0);  label("Allometric exponent of body weight on V1/F, V2/F and V3/F (unitless)")     # Majid 2024 Text S1 $THETA 13: '1 FIX ;[WGT_Volume]'; Table 1 footnote V1/F = 46.0 * (BW/73.2)

    # ---- Covariate effects on CL/F. The source reports each as a
    # multiplicative ratio entering the power form theta^covariate; they are
    # stored on the log scale here so exp(e_* * cov) reproduces theta^cov.
    e_cyp3a4_inh_cl  <- log(0.896); label("Log-effect of concomitant CYP3A4 inhibitor on CL/F (unitless)")   # Majid 2024 Table 1: CYP3A4 inhibitors ~ CL/F ratio 0.896 (%RSE 5.30); Discussion: CL/F decreased by 10.4 percent
    e_alb_cl         <- log(0.900); label("Log-effect of serum albumin < 30 g/L on CL/F (unitless)")         # Majid 2024 Table 1: Albumin ~ CL/F ratio 0.900 (%RSE 1.98); Discussion: CL/F decreased by 10.0 percent with albumin below 30 g/L
    e_alp_cl         <- log(0.910); label("Log-effect of alkaline phosphatase > ULN on CL/F (unitless)")     # Majid 2024 Table 1: ALP (> ULN) ~ CL/F ratio 0.910 (%RSE 0.956); Discussion: CL/F decreased by 9 percent
    e_dis_healthy_cl <- log(1.19);  label("Log-effect of healthy-subject cohort on CL/F (unitless)")         # Majid 2024 Table 1: Healthy population ~ CL/F ratio 1.19 (%RSE 4.81)
    e_tumtp_dtc_cl   <- log(0.951); label("Log-effect of DTC tumor type on CL/F (unitless)")                 # Majid 2024 Table 1: DTC population ~ CL/F ratio 0.951 (%RSE 3.65)
    e_tumtp_hcc_cl   <- log(0.862); label("Log-effect of HCC tumor type on CL/F (unitless)")                 # Majid 2024 Table 1: HCC population ~ CL/F ratio 0.862 (%RSE 3.63)
    e_tumtp_rcc_cl   <- log(0.851); label("Log-effect of RCC tumor type on CL/F (unitless)")                 # Majid 2024 Table 1: RCC population ~ CL/F ratio 0.851 (%RSE 3.48)

    # ---- Inter-individual variability. Majid 2024 Table 1 reports IIV under
    # a 'CV%' header whose abbreviation footnote defines it explicitly as
    # 'CV%, square root of variance x 100' -- i.e. the printed number is
    # omega x 100 on the log scale, NOT a linear-space coefficient of
    # variation. The variances entered here are therefore (CV%/100)^2 with no
    # log-normal back-conversion. This reading is corroborated by every
    # $OMEGA initial estimate in Text S1, which is a variance on the same
    # scale and lands within a few percent of the final value: CL 0.15 vs
    # 0.127, V1/F 0.25 vs 0.183, V2/F 0.5 vs 0.465, ka 0.3 vs 0.346,
    # D1 1.0 vs 0.935, V3/F 0.15 vs 0.131.
    etalcl       ~ 0.126736;  label("IIV variance on log CL/F")                    # Majid 2024 Table 1: CL/F 35.6 CV% -> (0.356)^2; Text S1 $OMEGA init 0.15
    etalvc       ~ 0.183184;  label("IIV variance on log V1/F")                    # Majid 2024 Table 1: V1/F 42.8 CV% -> (0.428)^2; Text S1 $OMEGA init 0.25
    etalvp       ~ 0.465124;  label("IIV variance on log V2/F")                    # Majid 2024 Table 1: V2/F 68.2 CV% -> (0.682)^2; Text S1 $OMEGA init 0.5
    etalvp2      ~ 0.131044;  label("IIV variance on log V3/F")                    # Majid 2024 Table 1: V3/F 36.2 CV% -> (0.362)^2; Text S1 $OMEGA init 0.15
    etalka       ~ 0.345744;  label("IIV variance on log ka")                      # Majid 2024 Table 1: Ka 58.8 CV% -> (0.588)^2; Text S1 $OMEGA init 0.3
    etalduration ~ 0.935089;  label("IIV variance on log D1")                      # Majid 2024 Table 1: D1 96.7 CV% -> (0.967)^2; Text S1 $OMEGA init 1

    # ---- Residual error. Text S1 $ERROR implements three mutually exclusive
    # branches selected by study type and time after dose:
    #   Y = F * (1 + EPS(1))              healthy-subject studies, 17.6 CV%
    #   Y = F * (1 + EPS(2))              patient studies,         37.5 CV%
    #   Y = F * (1 + EPS(3)) + EPS(4)     time after dose <= 2 h,  42.9 CV% + 20.5 ng/mL
    # nlmixr2 has no covariate-conditional residual-error construct, so the
    # patient-study proportional arm is carried here: it is the branch that
    # applies to the RR-DTC cohort this paper's PK/PD models are built on.
    # The other two arms are quoted in the vignette 'Assumptions and
    # deviations' section.
    propSd <- 0.375;  label("Proportional residual error, patient studies (fraction)")  # Majid 2024 Table 1: Proportional (CV%; Patients studies) = 37.5 (%RSE 1.31); Text S1 $SIGMA init 0.2 for EPS(2)
  })

  model({
    # ---- Reference covariate values.
    ref_wt  <- 73.2  # kg; Majid 2024 Table 1 footnote normalises every allometric term to BW/73.2 (Table S2 PK-population median weight).
    alp_uln <- 120   # U/L; representative adult ALP upper limit of normal, see covariateData[[ALP]]$notes.

    # ---- Binary lab indicators, coded exactly as in Text S1 $PK.
    alb_low  <- (ALB < 30)
    alp_high <- (ALP > alp_uln)

    # ---- Individual parameters.
    ka       <- exp(lka + etalka)
    duration <- exp(lduration + etalduration)

    cl <- exp(lcl + etalcl) *
      (WT / ref_wt)^e_wt_cl *
      exp(e_cyp3a4_inh_cl  * CONMED_CYP3A4_INH) *
      exp(e_alb_cl         * alb_low) *
      exp(e_alp_cl         * alp_high) *
      exp(e_dis_healthy_cl * DIS_HEALTHY) *
      exp(e_tumtp_dtc_cl   * TUMTP_DTC) *
      exp(e_tumtp_hcc_cl   * TUMTP_HCC) *
      exp(e_tumtp_rcc_cl   * TUMTP_RCC)

    vc  <- exp(lvc  + etalvc)  * (WT / ref_wt)^allo_v
    vp  <- exp(lvp  + etalvp)  * (WT / ref_wt)^allo_v
    vp2 <- exp(lvp2 + etalvp2) * (WT / ref_wt)^allo_v

    # No IIV on the inter-compartmental clearances (Text S1 assigns no ETA
    # to Q3 or Q4, and Table 1 reports no IIV row for Q1 or Q2).
    q  <- exp(lq)  * (WT / ref_wt)^e_wt_cl
    q2 <- exp(lq2) * (WT / ref_wt)^e_wt_cl

    # ---- ODE system. Text S1 uses ADVAN12 TRANS4, whose compartment 1 is
    # the absorption depot, and defines a duration parameter D1 for that
    # compartment. The dose therefore enters the depot as a zero-order
    # process of duration D1 and is then absorbed first-order at ka into the
    # central compartment -- a sequential zero-order-then-first-order
    # absorption, which is what the paper's prose calls 'simultaneous first-
    # and zero-order absorption'. There is no parallel dose split: the
    # control stream has a single dosing compartment and a single F1.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - (cl + q + q2) / vc * central +
      q  / vp  * peripheral1 +
      q2 / vp2 * peripheral2
    d/dt(peripheral1) <-  q  / vc * central - q  / vp  * peripheral1
    d/dt(peripheral2) <-  q2 / vc * central - q2 / vp2 * peripheral2

    # Zero-order input duration into the depot. Dose records must carry
    # rate = -2 so rxode2 uses the modelled duration.
    dur(depot) <- duration

    # Relative bioavailability: 1 for the tablet reference arm, 0.882 for
    # capsules (Table 1 footnote F1 = 1 * 0.882^FORM).
    f(depot) <- exp(lfcap * FORM_CAPSULE)

    # Text S1 sets S2 = V2/1000, i.e. amounts in mg over volumes in L give
    # mg/L = ug/mL, and the factor 1000 converts to the reported ng/mL.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
