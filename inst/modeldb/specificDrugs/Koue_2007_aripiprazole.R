Koue_2007_aripiprazole <- function() {
  description <- "Two-compartment population PK model for oral aripiprazole in healthy Japanese male volunteers (Koue 2007), with first-order absorption, an absorption lag time, body-weight linear scaling on Vc/F, Vp/F, Q/F, and CL/F, additive linear-deviation CYP2D6 intermediate- and poor-metabolizer effects on CL/F (Group 1 = extensive metabolizer reference), an additive linear-deviation itraconazole-coadministration (CYP3A4 inhibitor) effect on CL/F, independent inter-individual variability on every structural parameter, and a log-normal (exponential) residual error."
  reference <- "Koue T, Kubo M, Funaki T, Fukuda T, Azuma J, Takaai M, Kayano Y, Hashimoto Y. Nonlinear mixed effects model analysis of the pharmacokinetics of aripiprazole in healthy Japanese males. Biol Pharm Bull. 2007;30(11):2154-2158. doi:10.1248/bpb.30.2154"
  vignette <- "Koue_2007_aripiprazole"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed in Koue 2007 (single-, multi-, and itraconazole-coadministration trials in healthy Japanese male volunteers; weight range 48.2-75.2 kg, mean 60.7 kg per the Pharmacokinetic Data section). Enters every clearance and volume parameter as a multiplicative scaling per kg with no reference-weight normalisation, i.e. Vc/F = theta3 * WT, Vp/F = theta8 * WT, Q/F = theta9 * WT, CL/F = (theta4 - theta5*G2 - theta6*G3 - theta7*ITZ) * WT. Equivalent to a fixed allometric exponent of 1.0 on every dispositional parameter.",
      source_name        = "WT"
    ),
    CYP2D6_PM = list(
      description        = "CYP2D6 poor-metabolizer phenotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (extensive or intermediate metabolizer)",
      notes              = "1 = subject is a CYP2D6 poor metabolizer (genotype encoding no functional CYP2D6 activity; Koue 2007 Group 3 = CYP2D6 *5/*10, *10/*10, or *41/*41), 0 otherwise. Paired with CYP2D6_EM to encode the three-level EM/IM/PM phenotype using two binary indicators; in this paper Group 1 (EM) is the explicit reference rather than IM, so within model() the IM indicator is derived as (1 - CYP2D6_EM - CYP2D6_PM). Maps to the paper's binary G3 (Eq. 4): G3 = 1 for Group 3, 0 otherwise.",
      source_name        = "G3 (paper Eq. 4)"
    ),
    CYP2D6_EM = list(
      description        = "CYP2D6 extensive-metabolizer phenotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (intermediate or poor metabolizer)",
      notes              = "1 = subject is a CYP2D6 extensive metabolizer (Koue 2007 Group 1 = CYP2D6 *1/*1, *1/*2, or *2/*2), 0 otherwise. Paired with CYP2D6_PM to encode the three-level EM/IM/PM phenotype; Koue 2007 takes Group 1 (EM) as the model reference and Group 2 (IM) and Group 3 (PM) as additive linear-deviation shifts on CL/F. Within model() the implicit IM stratum is derived as (1 - CYP2D6_EM - CYP2D6_PM), so this column is required even though it does not appear in the paper's equation as a named indicator.",
      source_name        = "Group 1 indicator (paper Table 1 / Eq. 4); derived as 1 - G2 - G3"
    ),
    CONMED_AZOLE = list(
      description        = "Concomitant azole antifungal coadministration (itraconazole) indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant itraconazole)",
      notes              = "1 = aripiprazole observation falls within the itraconazole-coadministration period of the protocol-3 trial (itraconazole 100 mg once daily for 21 consecutive days, aripiprazole 3 mg on day 7 after itraconazole start; Koue 2007 Pharmacokinetic Data section), 0 otherwise. Only itraconazole was studied in Koue 2007; the canonical CONMED_AZOLE column pools any systemic azole antifungal, but here only itraconazole was administered. Maps to the paper's binary ITZ (Eq. 4): ITZ = 1 for the coadministration period, 0 otherwise.",
      source_name        = "ITZ (paper Eq. 4)"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 68L,
    n_studies        = 3L,
    n_observations   = NA_integer_,
    age_range        = "20-32 years (mean 23.1)",
    weight_range     = "48.2-75.2 kg (mean 60.7)",
    sex_female_pct   = 0,
    race_ethnicity   = "100% Japanese (per Pharmacokinetic Data section)",
    disease_state    = "Healthy adult Japanese male volunteers across three clinical trials: single-dose (protocol 1, n=26 -- the abstract reports 24 contributing subjects), multiple-dose (protocol 2, n=15), and itraconazole-coadministration (protocol 3, n=27).",
    dose_range       = "Oral aripiprazole (ABILIFY tablet): 6 mg single dose (protocol 1), 3 mg once daily for 14 d (protocol 2), and 3 mg single dose with and without itraconazole 100 mg once daily for 21 d (protocol 3).",
    regions          = "Japan",
    notes            = "CYP2D6 genotype distribution (Koue 2007 Table 1): Group 1 (EM = *1/*1, *1/*2, *2/*2) n=14; Group 2 (IM = *1/*5, *1/*10, *2/*5, *2/*10, *2/*41) n=39; Group 3 (PM = *5/*10, *10/*10, *41/*41) n=15. Genotyping was performed by PCR-RFLP and long-PCR for CYP2D6 *2, *4, *5, *10, *14, *18, *36, *41; CYP3A5 was not genotyped. All food effects on absorption were judged not significant (Results and Discussion). The parent dataset comes from two previously published clinical-pharmacology studies cited as references 13 and 14 of Koue 2007."
  )

  ini({
    # Structural parameters -- Koue 2007 Table 2 "Population Pharmacokinetic
    # Parameters of Aripiprazole" (point estimates) and Equations 1-7 (model
    # equations). Parameters in /F form are the apparent (oral-clearance and
    # oral-volume) values; absolute bioavailability F = 0.87 (Results and
    # Discussion) is NOT used inside this model so all CL and V quantities
    # are /F throughout.

    ltlag <- log(0.805) ; label("Absorption lag time ALAG (h)")                                                # Koue 2007 Table 2 theta1 = 0.805 h (Eq. 1)
    lka   <- log(2.65)  ; label("Absorption rate constant KA (1/h)")                                           # Koue 2007 Table 2 theta2 = 2.65 1/h (Eq. 2)
    lvc   <- log(3.84)  ; label("Apparent central volume Vc/F per kg body weight (L/kg)")                      # Koue 2007 Table 2 theta3 = 3.84 L/kg (Eq. 3)
    lcl   <- log(0.0645); label("Apparent oral clearance CL/F per kg body weight, Group 1 reference (L/h/kg)") # Koue 2007 Table 2 theta4 = 0.0645 L/h/kg (Eq. 4); reference = Group 1 = CYP2D6 extensive metabolizer
    lvp   <- log(1.54)  ; label("Apparent peripheral volume Vp/F per kg body weight (L/kg)")                   # Koue 2007 Table 2 theta8 = 1.54 L/kg (Eq. 5)
    lq    <- log(0.168) ; label("Apparent inter-compartmental clearance Q/F per kg body weight (L/h/kg)")      # Koue 2007 Table 2 theta9 = 0.168 L/h/kg (Eq. 6)

    # Covariate effects on apparent oral clearance -- Koue 2007 Eq. 4:
    #   CL/F_i = (theta4 - theta5*G2 - theta6*G3 - theta7*ITZ) * WT_i * exp(eta_CL/F)
    # The paper reports the e_* values as positive decrements that are
    # SUBTRACTED from theta4; we store them as negative additive shifts
    # (linear-deviation L/h/kg) so the model() expression follows the
    # standard "(typical_value + sum_of_effects)" pattern.
    e_2d6im_cl  <- -0.0135 ; label("CYP2D6 intermediate-metabolizer additive shift on CL/F (L/h/kg)")  # Koue 2007 Table 2 theta5 = 0.0135 L/h/kg (paper reports as a positive decrement: CL/F decreases by 0.0135 L/h/kg in Group 2)
    e_2d6pm_cl  <- -0.0293 ; label("CYP2D6 poor-metabolizer additive shift on CL/F (L/h/kg)")          # Koue 2007 Table 2 theta6 = 0.0293 L/h/kg (paper reports as a positive decrement: CL/F decreases by 0.0293 L/h/kg in Group 3)
    e_azole_cl  <- -0.0181 ; label("Itraconazole coadministration additive shift on CL/F (L/h/kg)")    # Koue 2007 Table 2 theta7 = 0.0181 L/h/kg (paper reports as a positive decrement: CL/F decreases by 0.0181 L/h/kg during itraconazole coadministration)

    # Inter-individual variability -- Koue 2007 Table 2 reports omega
    # (standard deviation on the internal log scale, NOT omega^2). The
    # paper text immediately above Eq. 1 reads "h is a random variable
    # distributed with a mean of zero and variance of omega^2_ALAG" --
    # i.e., omega^2 is the variance and omega (Table 2 column) is the SD.
    # nlmixr2's `~ value` syntax declares the variance, so we square each
    # SD here.
    etaltlag ~ 0.0121                       # Koue 2007 Table 2 omega_ALAG  = 0.110; variance = 0.110^2 = 0.0121
    etalka   ~ 0.8281                       # Koue 2007 Table 2 omega_KA    = 0.910; variance = 0.910^2 = 0.8281
    etalvc   ~ 0.1452                       # Koue 2007 Table 2 omega_V1/F  = 0.381; variance = 0.381^2 = 0.1452
    etalcl   ~ 0.1576                       # Koue 2007 Table 2 omega_CL/F  = 0.397; variance = 0.397^2 = 0.1576
    etalvp   ~ 0.1332                       # Koue 2007 Table 2 omega_V2/F  = 0.365; variance = 0.365^2 = 0.1332
    etalq    ~ 0.0686                       # Koue 2007 Table 2 omega_Q/F   = 0.262; variance = 0.262^2 = 0.0686

    # Residual unexplained variability -- Koue 2007 Eq. 7 writes the
    # residual error as C_ij = C*_ij * exp(eps_ij) with eps ~ N(0, sigma^2),
    # i.e. log-normal / exponential residual on the linear-scale prediction.
    # Table 2 reports sigma = 0.166 (the SD on the log scale).
    expSd <- 0.166 ; label("Log-normal residual error SD (log scale)")  # Koue 2007 Table 2 sigma = 0.166; Eq. 7 exponential residual error model
  })

  model({
    # CYP2D6 intermediate-metabolizer indicator derived from the canonical
    # paired CYP2D6_PM / CYP2D6_EM register: a subject is IM exactly when
    # they are neither PM nor EM. Koue 2007 Eq. 4 takes Group 1 (EM) as
    # the reference category so the IM indicator G2 corresponds to "not
    # EM and not PM".
    ind_cyp2d6_im <- 1 - CYP2D6_EM - CYP2D6_PM

    # Typical-value structural parameters. Body weight enters every
    # clearance and volume as a per-kg multiplier with no reference-weight
    # normalisation (Koue 2007 Equations 3, 4, 5, 6). CYP2D6 IM / PM and
    # itraconazole effects on CL/F are linear-deviation additive shifts
    # in L/h/kg added to theta4 before the WT scaling.
    tvtlag <- exp(ltlag)
    tvka   <- exp(lka)
    tvvc   <- exp(lvc)  * WT
    tvcl   <- (exp(lcl) + e_2d6im_cl * ind_cyp2d6_im + e_2d6pm_cl * CYP2D6_PM + e_azole_cl * CONMED_AZOLE) * WT
    tvvp   <- exp(lvp)  * WT
    tvq    <- exp(lq)   * WT

    # Individual PK parameters -- all etas are log-scale (parameter_i =
    # typical_value * exp(eta)).
    tlag <- tvtlag * exp(etaltlag)
    ka   <- tvka   * exp(etalka)
    vc   <- tvvc   * exp(etalvc)
    cl   <- tvcl   * exp(etalcl)
    vp   <- tvvp   * exp(etalvp)
    q    <- tvq    * exp(etalq)

    # Micro-constants for the two-compartment oral-absorption ODE system.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment oral PK with first-order absorption and an
    # absorption lag time. Dose lands in `depot`. Bioavailability is
    # implicit in every /F parameter -- the popPK model was developed
    # from oral-only data and cannot separate F from CL or volumes.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    alag(depot)       <- tlag

    # Aripiprazole plasma concentration. The dose is in mg, central
    # amount is in mg, vc is in L; multiplying by 1000 gives ng/mL.
    Cc <- central / vc * 1000
    Cc ~ lnorm(expSd)
  })
}
