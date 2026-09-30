Martial_2021_tacrolimus <- function() {
  description <- "Two-compartment population pharmacokinetic model for once-daily oral meltdose tacrolimus (Envarsus) in stable adult liver transplant recipients converted from prolonged-release tacrolimus (Martial 2021). Delayed absorption uses the Savic transit-compartment input (1.58 transit compartments, mean transit time 3.39 h) feeding a first-order absorption compartment. Oral bioavailability is fixed at 0.23 with log-normal between-subject and between-occasion variability, and the peripheral volume is fixed at 500 L. Body weight enters by fixed allometry (exponent 0.75 on CL and Q, 1 on both volumes, reference 70 kg). Log-normal IIV on CL, Vc, Q and ka. Proportional residual error differs between venous whole-blood and dried-blood-spot samples."
  reference <- "Martial LC, Biewenga M, Ruijter BN, Keizer R, Swen JJ, van Hoek B, Moes DJAR. Population pharmacokinetics and genetics of oral meltdose tacrolimus (Envarsus) in stable adult liver transplant recipients. Br J Clin Pharmacol. 2021;87(11):4262-4272. doi:10.1111/bcp.14842"
  vignette <- "Martial_2021_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")
  # Oral doses only; the dose record targets depot (bolus suppressed by
  # f(depot) <- 0, drug delivered by the Savic transit input).
  dosing <- "depot"

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling on CL and Q (exponent fixed 0.75) and on Vc and Vp (exponent fixed 1), standardised to 70 kg (Martial 2021 Methods 2.6; supplementary file S2 NONMEM $PK). Body weight was included a priori and was the only covariate retained in the final model. Cohort median 81.5 kg (range 54-133 kg; Table 1).",
      source_name = "WT"
    ),
    OCC = list(
      description = "Integer occasion index for the between-occasion variability on oral bioavailability",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "Martial 2021 Methods 2.4 defines occasions 1, 2 and 3 as the first (full, 2 weeks after conversion), second (abbreviated, 3 months after conversion) and optional third (clinical-care) AUC measurements, and Table 3 prints one shrinkage per occasion (28, 40, 98%). Decomposed into indicators oc1..oc3 inside model(); any other value (e.g. 0) carries no IOV. The InsightRX re-implementation in supplementary file S2 instead bins occasions by 24-h clock windows, which is a forecasting device and not the estimation design.",
      source_name = "occasion (AUC1 / AUC2 / AUC3)"
    ),
    SAMPLE_DBS = list(
      description = "Per-observation sampling-matrix indicator (1 = dried-blood-spot sample; 0 = venous whole-blood sample)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (venous whole blood)",
      notes = "Selects the matrix-specific proportional residual error only; Martial 2021 estimated no systematic DBS-to-venous conversion factor. Tacrolimus is quantified in whole blood, so the reference level for this model is venous whole blood (there are no plasma records). In the study the 0-6 h samples of the full AUC were venous whole blood and the 8, 12 and 24 h samples and the whole abbreviated AUC were dried blood spots (Methods 2.2). The supplementary file S2 control stream codes the same switch as DBS = 1 (whole blood, ERR(1)) and DBS = 2 (dried blood spot, ERR(2)); SAMPLE_DBS = DBS - 1.",
      source_name = "DBS"
    )
  )

  covariatesDataExcluded <- list(
    HCT = list(
      description = "Hematocrit",
      units = "L/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL and V1 (linear, hockey-stick, exponential and power forms). Significant in univariate analysis (lower hematocrit, higher CL) but not significant in the stepwise covariate model; not retained (Results 3.2). Cohort median 0.41 L/L (range 0.30-0.51; Table 1).",
      source_name = "Haematocrit"
    ),
    CYP3A5_EXPR = list(
      description = "Recipient CYP3A5 expresser status (1 = at least one CYP3A5*1 allele at rs776746)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CYP3A5*3/*3 nonexpresser)",
      notes = "Screened on CL and F. Recipient expressers had on average 43% higher CL in univariate analysis, but the effect was not significant in the stepwise covariate model and was not retained (Results 3.2). 12 of 54 recipients were expressers (Table 2).",
      source_name = "Recipient CYP3A5*3"
    ),
    CYP3A5_EXPR_DONOR = list(
      description = "Donor (graft) CYP3A5 expresser status",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CYP3A5*3/*3 nonexpresser graft)",
      notes = "Screened on CL and F, alone and as the recipient + donor combination C1-C4; not retained (Results 3.2). 12 of 49 donors were expressers (Table 2).",
      source_name = "Donor CYP3A5*3"
    ),
    SNP_CYP3A4_RS35599367 = list(
      description = "Recipient CYP3A4*22 (rs35599367) carrier indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (C/C noncarrier)",
      notes = "Screened on CL and F (recipient, donor and combination); not retained (Results 3.2, Figure 3). 4 of 54 recipients and 3 of 49 donors were C/T carriers (Table 2). IL-6 (rs1800796), IL-10 (rs1800871) and IL-18 (rs5744247) genotypes of recipient and donor were also screened on CL by stepwise covariate modelling and not retained.",
      source_name = "CYP3A4*22"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 53L,
    n_enrolled = 55L,
    n_observations = 748L,
    age_range = "21-70 years",
    age_median = "57 years",
    weight_range = "54-133 kg",
    weight_median = "81.5 kg",
    height_median = "175 cm (range 151-189 cm)",
    sex_female_pct = 34.5,
    race_ethnicity = "Caucasian 87% (48 of 55); others not stratified.",
    disease_state = "Stable adult liver transplant recipients at least 6 months after transplantation (median 67 months, range 6-240), on a stable prolonged-release tacrolimus (Advagraf) regimen, converted to meltdose tacrolimus (Envarsus) at a 1:0.7 dose ratio. Indications: primary sclerosing cholangitis 21.8%, hepatocellular carcinoma 21.8%, alcoholic liver disease 12.7%, hepatitis C 7.3%, polycystic liver disease 7.3%, other 29%; 9.1% retransplantations.",
    dose_range = "Once-daily oral Envarsus, median 2 mg (range 0.75-6 mg) (Table 1), adjusted to individual whole-blood trough targets.",
    regions = "The Netherlands (Leiden University Medical Center).",
    co_medication = "Mycophenolate mofetil 49%, prednisone 9.1%, everolimus 5.5%, sirolimus 3.6%, azathioprine 3.6% (Table 1).",
    sampling_design = "Full AUC (0, 1, 2, 3, 4, 6, 8, 12, 24 h) 2 weeks after conversion, abbreviated AUC (0, 4, 8, 12 h) 3 months after conversion, a trough sample, and an optional third clinical-care AUC. The 0-6 h samples of the full AUC were venous whole blood; the 8-24 h samples and the abbreviated AUC were dried blood spots. Tacrolimus was measured by LC-MS/MS.",
    notes = "Baseline demographics from Martial 2021 Table 1 (all 55 enrolled). Two patients were excluded from the PK analysis for data inconsistencies, leaving 53 patients with 748 concentrations. Table 1 median Envarsus AUC was 144 ug*h/L (range 25-323)."
  )

  ini({
    # Martial 2021 Table 3 'Final model' column. Bioavailability is fixed at
    # 0.23 and applied to the dose inside the transit input (supplementary
    # file S2 NONMEM: BIO = THETA(6)*EXP(ETA(5)+IOV) multiplying PODO), so the
    # clearances and volumes are systemic values even though Table 3 labels
    # them CL/F, V1/F, Q/F and V2/F. The Discussion confirms this: it
    # back-calculates the apparent clearance as 3.27/0.23 = 14.2 L/h, and the
    # Methods define the model AUC as DOSE*F/CL.
    lcl <- log(3.27); label("Clearance CL at 70 kg (L/h)") # Table 3 'CL/F (L/h)' = 3.27 (RSE 8%); S2 $THETA(1) = 3.27
    lvc <- log(94.9); label("Central volume of distribution Vc at 70 kg (L)") # Table 3 'V1/F (L)' = 94.9 (RSE 29%); S2 $THETA(2) = 94.9
    lq <- log(9.62); label("Intercompartmental clearance Q at 70 kg (L/h)") # Table 3 'Q/F (L/h)' = 9.62 (RSE 14%); S2 $THETA(4) = 9.62
    lvp <- fixed(log(500)); label("Peripheral volume of distribution Vp at 70 kg (L)") # Table 3 'V2/F (L) (fixed)' = 500; Results 3.2 fixed 'based on literature'; S2 $THETA(5) 500 FIX
    lka <- log(2.97); label("First-order absorption rate constant ka out of the absorption compartment (1/h)") # Table 3 'Ka (h-1)' = 2.97 (RSE 112%); S2 $THETA(3) prints 2.96
    lmtt <- log(3.39); label("Mean transit time of the Savic transit input (h)") # Table 3 'MTT (h)' = 3.39 (RSE 12%); S2 $THETA(7) = 3.39
    lntr <- log(1.58); label("Number of transit compartments of the Savic transit input (unitless)") # Table 3 'Ntrans' = 1.58 (RSE 18%); S2 $THETA(8) = 1.58
    lfdepot <- fixed(log(0.23)); label("Oral bioavailability F (fraction)") # Table 3 'F (fixed)' = 0.23; S2 $THETA(6) 0.23 FIX

    # Allometric exponents fixed a priori, reference weight 70 kg
    # (Methods 2.6; S2 $PK: CL and Q *(WT/70)**0.75, V2 and V3 *(WT/70)).
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent on CL and Q (unitless)") # Methods 2.6 'power exponent of 0.75'; S2 $PK
    e_wt_vc_vp <- fixed(1); label("Allometric exponent on Vc and Vp (unitless)") # Methods 2.6 'power exponent of 1.0'; S2 $PK

    # IIV. Table 3 prints CV% = 100*sqrt(omega^2); the omega^2 values are the
    # S2 $OMEGA estimates (0.116 -> 34%, 1.99 -> 141%, 0.0571 -> 24%,
    # 3.02 -> 174%, 0.131 -> 36%). The CL/V1 $OMEGA BLOCK(2) has a zero
    # covariance, so the etas are independent.
    etalcl ~ 0.116 # Table 3 'CL/F (CV%)' = 34; S2 $OMEGA 'IIV CL' 0.116
    etalvc ~ 1.99 # Table 3 'V1/F (CV%)' = 141; S2 $OMEGA (second element of BLOCK(2)) 1.99
    etalq ~ 0.0571 # Table 3 'Q/F (CV%)' = 24; S2 $OMEGA 'IIV Q' 0.0571
    etalka ~ 3.02 # Table 3 'Ka (CV%)' = 174; S2 $OMEGA 'IIV KA' 3.02
    etalfdepot ~ 0.131 # Table 3 'F (CV%)' = 36; S2 $OMEGA 'IIV F' 0.131

    # Between-occasion variability on F, one eta per AUC occasion with a
    # shared variance (S2 $OMEGA BLOCK(1) 0.0388 then BLOCK(1) SAME).
    etaiov_fdepot_1 ~ 0.0388 # Table 3 'IOV F (block) (CV%)' = 19.7; S2 $OMEGA BLOCK(1) 0.0388
    etaiov_fdepot_2 ~ fixed(0.0388) # BLOCK(1) SAME as occasion 1
    etaiov_fdepot_3 ~ fixed(0.0388) # BLOCK(1) SAME as occasion 1

    # Proportional residual error per sampling matrix. Table 3 prints
    # 100*sqrt(sigma^2); S2 $SIGMA gives 0.0111 (whole blood) and 0.062 (DBS).
    propSd <- 0.10536; label("Proportional residual error, venous whole blood (fraction)") # Table 3 'Whole blood (%)' = 10.5; S2 $SIGMA 0.0111, sqrt = 0.10536
    propSdDbs <- 0.24900; label("Proportional residual error, dried blood spot (fraction)") # Table 3 'DBS (V%)' = 24.9; S2 $SIGMA 0.062, sqrt = 0.24900
  })

  model({
    # Occasion indicators for the IOV on F (OCC = 1, 2, 3)
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 + oc3 * etaiov_fdepot_3

    # Individual PK parameters (S2 $PK)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q <- exp(lq + etalq) * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp) * (WT / 70)^e_wt_vc_vp
    ka <- exp(lka + etalka)
    mtt <- exp(lmtt)
    ntr <- exp(lntr)
    fdepot <- exp(lfdepot + etalfdepot + iov_fdepot)

    # Savic transit rate constant, KTR = (NN+1)/MTT (S2 $PK)
    ktr <- (ntr + 1) / mtt
    # log(ntr!) by the Stirling approximation exactly as coded in the
    # estimation model (S2 $PK LNFAC; 2.5066 = sqrt(2*pi)). It differs from
    # lgamma(ntr + 1) by 0.07% at ntr = 1.58.
    lnfac <- log(2.5066) + (ntr + 0.5) * log(ntr) - ntr + log(1 + 1 / (12 * ntr))

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # The dose record targets depot with its bolus suppressed (S2: F1 = 0);
    # the bioavailable amount F*PODO arrives through the Savic gamma-density
    # input driven by the time since the last dose (S2 $DES DADT(1)).
    # podo(depot) is not bioavailability-adjusted, so f(depot) <- 0 does not
    # remove the dose from this term.
    d/dt(depot) <- exp(log(fdepot * podo(depot)) + log(ktr) + ntr * log(ktr * tad(depot)) -
      ktr * tad(depot) - lnfac) - ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- 0

    # Whole-blood tacrolimus in ug/L: mg / L * 1000 (S2 $PK S2 = V2/1000)
    Cc <- central / vc * 1000
    propSdCc <- propSd * (1 - SAMPLE_DBS) + propSdDbs * SAMPLE_DBS
    Cc ~ prop(propSdCc)
  })
}
