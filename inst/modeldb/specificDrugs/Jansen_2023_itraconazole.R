Jansen_2023_itraconazole <- function() {
  description <- "Semi-mechanistic population PK model for intravenous itraconazole nanocrystal formulation (NCF) and its active metabolite hydroxy-itraconazole in allogeneic haematopoietic cell transplant recipients (Jansen 2023). A nanocrystal-bound itraconazole compartment receives the infusion and dissolves by a fixed first-order rate constant into a two-compartment dissolved-itraconazole disposition model; all eliminated itraconazole (fraction metabolised fixed to 1) enters a one-compartment hydroxy-itraconazole model. Observed itraconazole is the sum of nanocrystal-bound and dissolved concentrations. Allometric weight scaling with fixed exponents 0.75 on clearances and 1 on volumes. All amounts and concentrations are molar, as the source data were converted to molar equivalents before fitting."
  reference   <- "Jansen AME, Ter Heine R, Donnelly JP, Blijlevens N, Bruggemann RJM. Repurposing antifungals: population pharmacokinetics of itraconazole and hydroxy-itraconazole following administration of a nanocrystal formulation. J Antimicrob Chemother. 2023;78(5):1172-1178. doi:10.1093/jac/dkad072"
  vignette    <- "Jansen_2023_itraconazole"
  units       <- list(time = "h", dosing = "mmol", concentration = "mmol/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against the supplementary final model control
  # stream ($MODEL block: CEN / PERI / MET / NANO) and the Methods paragraph
  # describing the nanocrystal dissolution compartment.
  compartmentData <- list(
    central_np   = list(analyte = "itraconazole (nanocrystal-bound)", units = "mmol", specimen = "plasma", verified = TRUE),
    central      = list(analyte = "itraconazole (dissolved)",         units = "mmol", specimen = "plasma", verified = TRUE),
    peripheral1  = list(analyte = "itraconazole (dissolved)",         units = "mmol", specimen = "plasma", verified = TRUE),
    central_ohi  = list(analyte = "hydroxy-itraconazole",             units = "mmol", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight recorded at baseline.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "A priori allometric scaling to a 70 kg reference: (WT/70)^0.75 on itraconazole CL, intercompartmental Q and hydroxy-itraconazole CL; (WT/70)^1 on every volume, including the nanocrystal-bound volume, which the control stream sets equal to the itraconazole central volume. Exponents were fixed, not estimated. Body weight is the only covariate in the model; the study was explicitly not powered to screen covariates.",
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 10,
    n_studies      = 1,
    age_range      = "22.0-59.0 years",
    age_median     = "47.5 years",
    weight_range   = "60.0-92.5 kg",
    weight_median  = "78.3 kg",
    height_range   = "164-190 cm",
    bmi_range      = "20.6-29.4 kg/m2",
    sex_female_pct = 50,
    disease_state  = "Adults (18-65 years) receiving a matched allogeneic bone marrow transplant after conditioning with idarubicin, cyclophosphamide and total body irradiation; underlying disease acute lymphatic leukaemia (30%), acute myeloid leukaemia (20%), chronic myelomonocytic leukaemia (20%), non-Hodgkin lymphoma (20%), myelofibrosis (10%). No signs or symptoms of fungal infection at inclusion.",
    dose_range     = "Itraconazole nanocrystal formulation 200 mg IV as a 2 h infusion twice daily for 2 days, then 200 mg once daily until day 14.",
    regions        = "The Netherlands (Radboud University Medical Center, Nijmegen).",
    n_observations = "471 itraconazole and 471 paired hydroxy-itraconazole plasma concentrations, all above the LLOQ (0.002 mg/L itraconazole, 0.005 mg/L hydroxy-itraconazole). Observed ranges 0.06-6.96 mg/L (itraconazole) and 0.09-1.86 mg/L (hydroxy-itraconazole).",
    notes          = "Baseline characteristics from Jansen 2023 Table 1. Prospective open-label Phase II study; identification of covariates was not part of the study design. Full PK curves on days 7 and 14, pre- and post-infusion samples until day 6, pre-infusion on days 10 and 12, and washout samples on days 16, 17, 18, 19 and 28."
  )

  ini({
    # Structural parameters - dissolved itraconazole, reference 70 kg.
    # Final estimates from Jansen 2023 Table 2, identical to the $THETA block
    # of the supplementary final model control stream.
    lcl <- log(4.29);  label("Itraconazole clearance, CL_P, at 70 kg (L/h)")                                     # Table 2: CL_P 4.29 L/h (95% CI 3.85-4.80); $THETA TVCLP
    lvc <- log(14.1);  label("Itraconazole central volume of distribution, V_P1, at 70 kg (L)")                  # Table 2: V_P1/V_N 14.1 L (95% CI 11.8-17.2); $THETA TVVP1/TVVN
    lq  <- log(53.0);  label("Itraconazole intercompartmental clearance, Q_P, at 70 kg (L/h)")                   # Table 2: Q_P 53.0 L/h (95% CI 45.3-63.4); $THETA TVQP
    lvp <- log(1660);  label("Itraconazole peripheral volume of distribution, V_P2, at 70 kg (L)")               # Table 2: V_P2 1660 L (95% CI 1558-1775); $THETA TVVP2

    # Structural parameters - hydroxy-itraconazole. Both are apparent to the
    # fraction metabolised, which the paper fixed to 1 (Methods, Population
    # pharmacokinetic analysis), so no separate fm parameter is carried.
    lcl_ohi <- log(2.86);  label("Hydroxy-itraconazole clearance, CL_M/fm, at 70 kg (L/h)")                      # Table 2: CL_M 2.86 L/h (95% CI 2.43-3.33); $THETA TVCLM
    lvc_ohi <- log(43.1);  label("Hydroxy-itraconazole volume of distribution, V_M/fm, at 70 kg (L)")            # Table 2: V_M 43.1 L (95% CI 36.9-50.4); $THETA TVVM

    # Dissolution of itraconazole from the nanocrystals. Not estimated: the
    # paper fixes it from a 9 min dissolution half-life in human plasma
    # (log(2)/0.15 h = 4.62 /h) reported for the same formulation.
    lkdiss <- fixed(log(4.62));  label("Nanocrystal dissolution rate constant, kN (1/h)")                        # Methods: "The rate of dissolution from the nanocrystals (kN) was fixed at 4.62 h-1 based on a dissolution half-life of 9 min in human plasma"; $PK KN = 4.62

    # Allometric exponents, fixed a priori (Methods, Population
    # pharmacokinetic analysis) and applied in $PK as (WT/70)**0.75 and
    # (WT/70)**1.
    e_wt_cl <- fixed(0.75);  label("Allometric exponent on clearances (unitless)")                               # Methods: "allometrically scaled to a total body weight of 70 kg with a fixed exponent of 3/4 for (intercompartmental) clearance"
    e_wt_vc <- fixed(1.0);   label("Allometric exponent on volumes (unitless)")                                  # Methods: "and 1 for volumes of distribution"

    # IIV. Table 2 reports %CV back-transformed from the log-normal variance
    # via sqrt(exp(omega^2) - 1); the variances themselves are printed in the
    # supplementary $OMEGA block and are used directly here.
    #   CL_P: 0.0126 -> sqrt(exp(0.0126) - 1) = 11.3% (Table 2 IIV CL_P)
    #   CL_M: 0.0506 -> sqrt(exp(0.0506) - 1) = 22.8% (Table 2 IIV CL_M)
    # The control stream also declares OMEGA entries for V_P1/V_N, Q_P, V_P2
    # and V_M, but all four are "0 FIX"; they are omitted here because a
    # zero-variance diagonal makes the OMEGA matrix singular for simulation.
    etalcl     ~ 0.0126                                                                                          # $OMEGA 0.0126 ; IIV CLP
    etalcl_ohi ~ 0.0506                                                                                          # $OMEGA 0.0506 ; IIV CLM

    # Residual error. The paper fitted an additive error on the log scale for
    # both analytes ($ERROR: IPRED = LOG(...); Y = IPRED + EPS), which is the
    # nlmixr2 lnorm() structure; expSd is the SD on the log scale. Table 2
    # reports these as 48.8% and 19.1% after the same sqrt(exp(s^2) - 1)
    # back-transformation applied to the IIV.
    expSd     <- sqrt(0.214);   label("Itraconazole residual SD on the log scale (log units)")                   # $SIGMA BLOCK(2) diagonal 0.214; Table 2 Error_P 48.8% = sqrt(exp(0.214) - 1)
    expSd_ohi <- sqrt(0.0358);  label("Hydroxy-itraconazole residual SD on the log scale (log units)")           # $SIGMA BLOCK(2) diagonal 0.0358; Table 2 Error_M 19.1% = sqrt(exp(0.0358) - 1)
    # The published model additionally correlates the two residuals
    # ($SIGMA off-diagonal 0.0185, i.e. correlation 0.211; Table 2 reports
    # 48.5% after the same back-transformation). nlmixr2 has no syntax for a
    # correlated residual-error block across endpoints, so the correlation is
    # not reproduced here - see the vignette Assumptions and deviations.
  })

  model({
    # Reference body weight for the a priori allometric scaling.
    ref_wt <- 70

    # Individual parameters - dissolved itraconazole ($PK CLP, VP1, QP, VP2).
    cl <- exp(lcl + etalcl) * (WT / ref_wt)^e_wt_cl
    vc <- exp(lvc)          * (WT / ref_wt)^e_wt_vc
    q  <- exp(lq)           * (WT / ref_wt)^e_wt_cl
    vp <- exp(lvp)          * (WT / ref_wt)^e_wt_vc

    # Individual parameters - hydroxy-itraconazole ($PK CLM, VM).
    cl_ohi <- exp(lcl_ohi + etalcl_ohi) * (WT / ref_wt)^e_wt_cl
    vc_ohi <- exp(lvc_ohi)              * (WT / ref_wt)^e_wt_vc

    # The nanocrystal-bound volume is set equal to the itraconazole central
    # volume ($PK VN = THETA(2)*((WT/70)**1)*EXP(ETA(2))); the paper states
    # this is an assumption made in the absence of nanocrystal-bound
    # concentration measurements.
    kdiss <- exp(lkdiss)

    # Micro-constants, named as in the control stream's $PK rate-constant
    # block. There is no K10: with fm fixed to 1 every unit of itraconazole
    # leaving the central compartment by clearance appears as
    # hydroxy-itraconazole.
    kpm     <- cl / vc          # K13, itraconazole -> hydroxy-itraconazole
    k12     <- q / vc           # K12
    k21     <- q / vp           # K21
    kel_ohi <- cl_ohi / vc_ohi  # K30

    # ODE system (ADVAN5: NANO -> CEN <-> PERI, CEN -> MET -> out).
    # The infusion is administered into the nanocrystal-bound compartment.
    d/dt(central_np)  <- -kdiss * central_np
    d/dt(central)     <-  kdiss * central_np - kpm * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(central_ohi) <-  kpm * central - kel_ohi * central_ohi

    # Observations. $ERROR fits the itraconazole record against the SUM of
    # the nanocrystal-bound and nanocrystal-unbound concentrations
    # (IPRED = LOG(C1 + C4)), both divided by the same volume. Cc_np is the
    # nanocrystal-bound term on its own; it is not measured separately and so
    # carries no residual error, but the paper plots it in Figure 1.
    Cc     <- (central + central_np) / vc
    Cc_np  <- central_np / vc
    Cc_ohi <- central_ohi / vc_ohi

    Cc     ~ lnorm(expSd)
    Cc_ohi ~ lnorm(expSd_ohi)
  })
}
