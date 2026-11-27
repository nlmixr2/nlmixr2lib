Goncalves_2017_primaquine <- function() {
  description <- "Semi-physiological population PK model for a single low oral dose of primaquine and its carboxyprimaquine metabolite in African children aged 2-14 years. The dose passes through a depot and two absorption transit compartments (ktr = 3/MAT) into a well-stirred liver compartment (volume from body surface area) that exchanges with a one-compartment primaquine plasma pool at the hepatic plasma flow. Hepatic intrinsic clearance is split into a monoamine-oxidase pathway that forms carboxyprimaquine (one-compartment disposition) and a CYP2D6 pathway proportional to the CYP2D6 activity score. Relative bioavailability matures with age as AGE/(AGE + F50). Flows and clearances scale allometrically to 70 kg with exponent 0.75 and volumes with exponent 1."
  reference <- paste(
    "Goncalves BP, Pett H, Tiono AB, Murry D, Sirima SB, Niemi M, Bousema T,",
    "Drakeley C, ter Heine R. Age, weight, and CYP2D6 genotype are major",
    "determinants of primaquine pharmacokinetics in African children.",
    "Antimicrob Agents Chemother. 2017;61(5):e02590-16.",
    "doi:10.1128/AAC.02590-16.",
    "Liver volume equation (cited by the paper but not printed in it):",
    "Johnson TN, Tucker GT, Tanner MS, Rostami-Hodjegan A. Changes in liver",
    "volume from birth to adulthood: a meta-analysis. Liver Transpl.",
    "2005;11(12):1481-1493. doi:10.1002/lt.20519.",
    sep = " "
  )
  vignette <- "Goncalves_2017_primaquine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline. Allometric scaling of every flow and clearance (exponent 0.75: QH, CLint,MAO, CLint,CYP2D6, CLCPQ) and every volume (exponent 1: VPQ, VCPQ) to a 70 kg reference (Methods 'Pharmacokinetic modeling'; Table 2 footnote). Also enters body surface area for the liver volume. Cohort medians 19.8 kg (0.25 mg/kg arm) and 24.9 kg (0.40 mg/kg arm), Table 1.",
      source_name = "body weight"
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline. Enters only the liver volume, through body surface area (Methods: 'liver volume (VL) was calculated from total body weight and height (44)'). Cohort medians 118 cm and 129 cm, Table 1.",
      source_name = "height"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline. Drives the maturation of relative oral bioavailability F = AGE / (AGE + F50) (Results 'Pharmacokinetics'); body weight alone does not carry it. Cohort 2-14 years (medians 8 and 10 years, Table 1).",
      source_name = "age"
    ),
    CYP2D6 = list(
      description = "CYP2D6 activity score (sum of the two allele activity values)",
      units = "(unitless)",
      type = "continuous",
      reference_category = NULL,
      notes = "Genotype-derived activity score (Methods 'CYP2D6 genotyping', reference 41). Multiplies the population CYP2D6 intrinsic clearance linearly: CLint,CYP2D6,i = AS * CLint,CYP2D6,pop, so a poor metabolizer (AS = 0) has no CYP2D6 clearance. Observed: AS 0 (n = 1), 0.5 or 1.0 (n = 10), 1.5 or 2.0 (n = 22), 3.0 (n = 3); the four children with missing or inconclusive genotype were set to 1.5 in the final model. No value transformation: the activity score is the column value.",
      source_name = "AS"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "primaquine", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "primaquine", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "primaquine", units = "mg", specimen = "administration site", verified = TRUE),
    liver = list(analyte = "primaquine", units = "mg", specimen = "tissue", verified = TRUE),
    central = list(analyte = "primaquine", units = "mg", specimen = "plasma", verified = TRUE),
    central_cpq = list(analyte = "carboxyprimaquine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 38,
    n_studies = 1,
    age_range = "2-14 years",
    age_median = "8 years (0.25 mg/kg arm) and 10 years (0.40 mg/kg arm)",
    weight_range = "IQR 16.2-26.3 kg (0.25 mg/kg arm) and 15.9-31.6 kg (0.40 mg/kg arm)",
    weight_median = "19.8 kg (0.25 mg/kg arm) and 24.9 kg (0.40 mg/kg arm)",
    height_median = "118 cm (0.25 mg/kg arm) and 129 cm (0.40 mg/kg arm)",
    sex_female_pct = 57.5,
    disease_state = "Afebrile children with asymptomatic Plasmodium falciparum parasitaemia (37 of 40 patent asexual parasites), G6PD-normal, receiving artemether-lumefantrine (six doses over 3 days) with a single low primaquine dose given with the fifth artemether-lumefantrine dose.",
    dose_range = "Single oral primaquine dose of 0.25 mg/kg (n = 20) or 0.40 mg/kg (n = 20), given as a 1 mg/mL solution of a crushed 15 mg tablet.",
    regions = "Burkina Faso (Balonghin)",
    cyp2d6_activity_score = "36 of 40 genotyped: AS 0 (poor metabolizer) n = 1; AS 0.5 or 1.0 (intermediate) n = 10; AS 1.5 or 2.0 (extensive) n = 22; AS 3.0 (ultrarapid) n = 3. Four missing or inconclusive genotypes were set to 1.5.",
    notes = "Baseline demographics from Table 1 (40 children, 55% and 60% female in the two arms). Two children were excluded from the PK analysis for inadequate sample volume, leaving 38 (Table 3: n = 18 and 20). Seven samples per child (one pre-dose, four in 0-12 h, two in 24-72 h; 274 samples). LC-MS, LLOQ 4 ng/mL; BLQ samples handled with Beal M6. Dose and concentrations were converted to molar units for the fit (Methods). NCT01935882."
  )

  ini({
    # Final estimates, Table 2. Flows and clearances are for a 70 kg
    # individual (allometric exponent 0.75) and volumes for 70 kg
    # (exponent 1), Table 2 footnote.
    lmat <- log(0.706); label("Mean absorption time MAT (h)") # Table 2 'Mean absorption time' = 0.706 h (RSE 12%)
    lvc <- log(127); label("Primaquine volume of distribution VPQ at 70 kg (L)") # Table 2 'VPQ (70 kg)' = 127 L (RSE 20%)
    lvc_cpq <- log(21.7); label("Carboxyprimaquine volume of distribution VCPQ at 70 kg (L)") # Table 2 'VCPQ (70 kg)' = 21.7 L (RSE 39%)
    lclint_mao <- log(7.35); label("Monoamine-oxidase intrinsic hepatic clearance CLint,MAO at 70 kg (L/h)") # Table 2 'CLint,MAO (70 kg)' = 7.35 L/h (RSE 41%)
    lclint_cyp2d6 <- log(6.70); label("Population CYP2D6 intrinsic hepatic clearance per unit activity score CLint,CYP2D6,pop at 70 kg (L/h)") # Table 2 'CLint,CYP2D6,pop (70 kg)' = 6.70 L/h (RSE 69%)
    lcl_cpq <- log(1.50); label("Carboxyprimaquine clearance CLCPQ at 70 kg (L/h)") # Table 2 'CLCPQ (70 kg)' = 1.50 L/h (RSE 35%)
    age50_fdepot <- 4.27; label("Age at which relative bioavailability is 50% of the mature value F50 (years)") # Table 2 'F50 (yrs)' = 4.27 (RSE 44%); Results F = age/(age + F50)

    # Liver plasma flow: adult liver blood flow 90 L/h times plasma
    # fraction 0.55 (haematocrit 45%), scaled allometrically like the
    # other flows (Methods 'Pharmacokinetic modeling').
    lqh <- fixed(log(49.5)); label("Liver plasma flow QH at 70 kg (L/h)") # Methods 'We assumed a liver plasma flow (QH) of 49.5 liters/h'

    # Allometric exponents fixed (Methods; Table 2 footnote).
    e_wt_cl <- fixed(0.75); label("Allometric exponent on all flows and clearances (QH, CLint,MAO, CLint,CYP2D6, CLCPQ)") # Methods 'fixed at 0.75 ... for flow'
    e_wt_vc <- fixed(1); label("Allometric exponent on all volumes (VPQ, VCPQ)") # Methods 'fixed at ... 1 for ... volume parameters'

    # Inter-individual variability, exponential model (Methods). Table 2
    # prints each as '%'; read as the log-normal CV and converted with
    # omega^2 = log(1 + CV^2). See the vignette Assumptions.
    etalvc ~ 0.523077 # Table 2 'VPQ interindividual variability (%)' 82.9 -> log(1 + 0.829^2)
    etalclint_mao ~ 0.355290 # Table 2 'CLint,MAO interindividual variability (%)' 65.3 -> log(1 + 0.653^2)

    # Residual error: combined proportional + additive for each analyte.
    # The additive SD is fixed at half the 4 ng/mL LLOQ (Methods, M6).
    propSd <- 0.327; label("Proportional residual error for primaquine (fraction)") # Table 2 'Residual error: PQ, Proportional (%)' = 32.7
    addSd <- fixed(2); label("Additive residual error for primaquine (ng/mL)") # Table 2 'Residual error: PQ, Additive (ng/ml)' = 2 (fixed)
    propSd_cpq <- 0.452; label("Proportional residual error for carboxyprimaquine (fraction)") # Table 2 'Residual error: CPQ, Proportional (%)' = 45.2
    addSd_cpq <- fixed(2); label("Additive residual error for carboxyprimaquine (ng/mL)") # Table 2 'Residual error: CPQ, Additive (ng/ml)' = 2 (fixed)
  })

  model({
    # 1. Body surface area and liver volume. The paper computes VL from
    #    weight and height with Johnson 2005 (its reference 44), whose
    #    equation is VL (L) = 0.722 * BSA^1.176. The BSA formula is not
    #    printed; the DuBois (< 15 kg) / Haycock (>= 15 kg) pair that
    #    Brussee_2018b_midazolam_children_pbpk uses with the same
    #    liver-volume equation is used here. VL only sets the liver transit time (VL/QH, about
    #    2 min), which leaves the plasma profiles unchanged; the vignette
    #    shows this.
    bsa_lo <- 0.007184 * HT^0.725 * WT^0.425
    bsa_hi <- 0.024265 * HT^0.3964 * WT^0.5378
    bsa <- (WT < 15) * bsa_lo + (WT >= 15) * bsa_hi
    vliver <- 0.722 * bsa^1.176

    # 2. Individual parameters with allometric scaling to 70 kg.
    mat <- exp(lmat)
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vc_cpq <- exp(lvc_cpq) * (WT / 70)^e_wt_vc
    clint_mao <- exp(lclint_mao + etalclint_mao) * (WT / 70)^e_wt_cl
    # CLint,CYP2D6,i = AS * CLint,CYP2D6,pop (Methods).
    clint_cyp2d6 <- CYP2D6 * exp(lclint_cyp2d6) * (WT / 70)^e_wt_cl
    cl_cpq <- exp(lcl_cpq) * (WT / 70)^e_wt_cl
    qh <- exp(lqh) * (WT / 70)^e_wt_cl

    # 3. Transit absorption, ktr = (n + 1)/MAT with n = 2 transit
    #    compartments (Methods; Figure 2: k12 = k23 = k34 = 3/MAT).
    ktr <- 3 / mat

    # 4. Well-stirred liver (Methods): EH = CLint / (QH + CLint), with
    #    CLint the sum of the two pathway intrinsic clearances, and
    #    CLH = EH * QH split between the pathways by their share of
    #    CLint. Figure 2: k40 = CLH,CYP2D6/VL, k45 = CLH,MAO/VL,
    #    k46 = QH (1 - EH)/VL, k64 = QH/VPQ.
    clint <- clint_mao + clint_cyp2d6
    eh <- clint / (qh + clint)
    clh_mao <- qh * clint_mao / (qh + clint)
    clh_cyp2d6 <- qh * clint_cyp2d6 / (qh + clint)

    # 5. Relative bioavailability matures with age (Results).
    fdepot <- AGE / (AGE + age50_fdepot)

    # 6. Molar conversion of the MAO metabolic flux. The paper fitted in
    #    molar units; amounts here are in mg, so the primaquine mass
    #    converted to carboxyprimaquine is multiplied by the molecular
    #    weight ratio (carboxyprimaquine C15H18N2O3, 274.32 g/mol;
    #    primaquine base C15H21N3O, 259.35 g/mol; the paper's monitored
    #    [M+H]+ ions are 275.25 and 260.30 m/z).
    mw_ratio <- 274.32 / 259.35

    # 7. ODEs (Figure 2).
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(liver) <- ktr * transit2 + qh * central / vc -
      qh * (1 - eh) * liver / vliver -
      (clh_mao + clh_cyp2d6) * liver / vliver
    d/dt(central) <- qh * (1 - eh) * liver / vliver - qh * central / vc
    d/dt(central_cpq) <- mw_ratio * clh_mao * liver / vliver -
      cl_cpq * central_cpq / vc_cpq

    f(depot) <- fdepot

    # 8. Observations: mg/L * 1000 = ng/mL.
    Cc <- 1000 * central / vc
    Cc_cpq <- 1000 * central_cpq / vc_cpq

    Cc ~ add(addSd) + prop(propSd)
    Cc_cpq ~ add(addSd_cpq) + prop(propSd_cpq)
  })
}
