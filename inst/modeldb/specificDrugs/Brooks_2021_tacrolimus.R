Brooks_2021_tacrolimus <- function() {
  description <- paste0(
    "Two-compartment population pharmacokinetic model for IV continuous-infusion ",
    "tacrolimus in pediatric and young adult patients undergoing allogeneic ",
    "hematopoietic cell transplantation (Brooks 2021). Allometric weight scaling ",
    "on all PK parameters with fixed theoretic exponents (0.75 on CL and Q, 1.0 ",
    "on V and V2; reference weight 70 kg); a structural ratio Fact fixed at 2.0 ",
    "links Q to CL and V2 to V; and a multiplicative azole-antifungal ",
    "(voriconazole or posaconazole) factor of 0.8 on CL captures the CYP3A4/5 ",
    "inhibitor co-treatment effect."
  )
  reference <- paste0(
    "Brooks JT, Keizer RJ, Long-Boyle JR, Kharbanda S, Dvorak CC, Friend BD. ",
    "Population Pharmacokinetic Model Development of Tacrolimus in Pediatric ",
    "and Young Adult Patients Undergoing Hematopoietic Cell Transplantation. ",
    "Front Pharmacol. 2021;12:750672. doi:10.3389/fphar.2021.750672."
  )
  vignette <- "Brooks_2021_tacrolimus"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight (actual body weight, ABW)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Allometric power scaling on all PK parameters with reference 70 kg ",
        "(Brooks 2021 Results / final-model equation block following Figure 4). ",
        "Theoretic exponents fixed at 0.75 on CL and Q and 1.0 on V and V2; ",
        "the paper estimated the exponents (0.73 on CL and 0.83 on V) and fixed ",
        "them at the canonical theoretic values because the estimates were close ",
        "to those values. Allometry on fat-free mass was also evaluated but did ",
        "not improve fit. Cohort median 23.9 kg (range 5.5-155.5 kg, Table 1)."
      ),
      source_name        = "WT"
    ),
    CONMED_AZOLE = list(
      description        = "Concomitant azole antifungal (voriconazole or posaconazole) indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant voriconazole or posaconazole)",
      notes              = paste0(
        "Time-varying. 1 = patient coadministered voriconazole or posaconazole ",
        "during the observation, 0 = otherwise (Brooks 2021 Methods / Results). ",
        "Voriconazole and posaconazole were prescribed for prophylaxis or ",
        "treatment of fungal infections; pooled into a single CYP3A4/5-inhibitor ",
        "covariate. 193 of 1648 concentration samples (11.7%) carried ",
        "CONMED_AZOLE = 1 in the model-building dataset. Effect on CL is a ",
        "multiplicative factor of 0.8 (Table 2 theta_INH = 0.8, RSE 6.97%; 20% ",
        "reduction in CL with 95% CI 10-33%)."
      ),
      source_name        = "CYP3A4/5 inhibitor (voriconazole or posaconazole)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 111L,
    n_studies      = 1L,
    age_range      = "0.5-25 years",
    age_median     = "7.3 years",
    weight_range   = "5.5-155.5 kg",
    weight_median  = "23.9 kg",
    sex_female_pct = 39.0,
    race_ethnicity = c(
      Caucasian_NonHispanic = 34.2,
      Caucasian_Hispanic    = 34.2,
      Asian                 = 16.2,
      African_American      = 5.4,
      Asian_Caucasian_Hispanic = 3.6,
      Other_Declined        = 2.7,
      American_Indian_Alaskan = 1.8,
      Multi_ancestry        = 1.8
    ),
    disease_state  = paste0(
      "Pediatric and young adult patients undergoing allogeneic hematopoietic ",
      "cell transplantation (HCT) for graft-versus-host disease prophylaxis. ",
      "58.6% malignant diagnoses (acute lymphoblastic leukemia 30.6%, acute ",
      "myeloid leukemia 13.5%, juvenile myelomonocytic leukemia 7.2%, others); ",
      "41.4% non-malignant (primary immunodeficiencies 17.1%, aplastic anemia / ",
      "bone-marrow-failure syndromes 14.4%, inborn errors of metabolism 6.3%, ",
      "hemoglobinopathies 3.6%). Conditioning regimens included busulfan / ",
      "fludarabine / clofarabine (34.3%), busulfan / fludarabine (12.6%), ",
      "cyclophosphamide / fludarabine (12.6%), cyclophosphamide / TBI (9.9%), ",
      "melphalan / fludarabine (9.0%), and other regimens (21.6%). Donor source: ",
      "bone marrow 49.5%, peripheral blood stem cells 44.1%, umbilical cord ",
      "blood 6.3%. Serotherapy: rabbit antithymocyte globulin 57.7%, alemtuzumab ",
      "37.8%, none 4.5%."
    ),
    dose_range     = paste0(
      "IV continuous infusion. Recommended starting dose 1.25 mcg/kg/h with a ",
      "goal therapeutic trough range of 7-10 ng/mL; observed dosing rate range ",
      "0.45-1.25 mcg/kg/h, median 1.25 mcg/kg/h (Brooks 2021 Results). ",
      "Continuous infusions were assumed to run over exactly 24 h with samples ",
      "drawn the 15 min before the next infusion start."
    ),
    regions        = "United States (UCSF Benioff Children's Hospital, San Francisco, CA)",
    co_medication  = paste0(
      "Voriconazole or posaconazole concomitant with tacrolimus in 193/1648 ",
      "(11.7%) of plasma samples. Ursodeoxycholic acid was given to nearly all ",
      "patients as standard practice and was therefore not modeled as a ",
      "covariate."
    ),
    n_concentrations = 1648L,
    sampling_design  = paste0(
      "1,648 steady-state trough plasma concentrations over a median of 14 days ",
      "after starting tacrolimus continuous IV infusion; troughs measured every ",
      "24-48 h after initiation and every 24 h after a dosing change. Median ",
      "initial trough 10.2 ng/mL (range 1.8-24.2 ng/mL). 929 (56.4%) samples ",
      "were outside the 7-10 ng/mL target window."
    ),
    outcomes       = paste0(
      "Acute GVHD 17.1%; chronic GVHD 14.4%; deceased at time of data collection ",
      "18.9% (Table 1). None of the transplant-specific covariates (donor source, ",
      "HLA mismatch, aGVHD, cGVHD, survival) reached significance in the ",
      "covariate analysis."
    ),
    notes          = paste0(
      "Single-center retrospective chart review (February 2016 to July 2020). ",
      "All concentrations are steady-state troughs from continuous IV infusion; ",
      "no oral / SC / IM dosing data are in the model-development dataset. ",
      "Software: NONMEM v7.4, PsN v4.8.1, PiranaJS; ggplot2 / vpc for diagnostics. ",
      "Baseline demographics from Brooks 2021 Table 1; concentration counts and ",
      "covariate evaluations from the Results section."
    )
  )

  ini({
    # ----- Structural PK (Brooks 2021 Table 2 final-model column) -----
    # Two-compartment IV continuous-infusion disposition. CL and V are the
    # primary estimated parameters; Q and V2 are linked to CL and V by the
    # fixed multiplication factor Fact = 2.0 (Table 2 footnote: Fact "fixed"),
    # so the typical Q = 2.0 * 4.2 = 8.4 L/h and typical V2 = 2.0 * 61.9 =
    # 123.8 L. Allometric exponents are fixed at the theoretic values 0.75 on
    # CL/Q and 1.0 on V/V2 (Brooks 2021 Results: estimated exponents 0.73 and
    # 0.83 were close to the theoretic 0.75 and 1.0 and were fixed there).
    # Reference weight is 70 kg (final-model equation block following Figure 4).
    lcl       <- log(4.2);   label("Clearance (CL, L/h) at reference weight 70 kg, no concomitant azole antifungal")  # Brooks 2021 Table 2 theta_CL = 4.2 L/h (RSE 2.95%)
    lvc       <- log(61.9);  label("Central volume of distribution (V, L) at reference weight 70 kg")                  # Brooks 2021 Table 2 theta_V = 61.9 L (RSE 5.98%)

    # CYP3A4/5-inhibitor (voriconazole or posaconazole) multiplicative factor
    # on CL. The paper's final-model equation has the form
    #   CL = theta_CL * (WT/70)^0.75 * theta_INH^CONMED_AZOLE * exp(eta_CL + kappa)
    # so theta_INH is applied as a power-of-binary-indicator: when
    # CONMED_AZOLE = 0 the factor is 1; when CONMED_AZOLE = 1 the factor is
    # theta_INH = 0.8 (20% reduction in CL, Results / Table 2 narrative
    # "estimate 20%, 95% CI 10-33%").
    e_azole_cl <- 0.8;       label("Multiplicative effect of concomitant azole antifungal on CL (CL_typical * e_azole_cl^CONMED_AZOLE)")  # Brooks 2021 Table 2 theta_INH = 0.8 (RSE 6.97%); equivalently a 20% reduction in CL

    # Structural ratio Fact: Q = Fact * CL and V2 = Fact * V. Brooks 2021
    # fixed Fact at 2.0 ("Fact was determined based on the average of several
    # published population PK studies, and it was assumed to be the same for
    # both Q and V2", citing Xue 2011, Kassir 2014, Moes 2016, Andrews 2018,
    # Andrews 2020). Encoded as a unitless fixed parameter so the
    # provenance of the linked Q/V2 typical values is explicit in ini().
    fact_q_vp <- fixed(2.0); label("Structural ratio linking Q to CL and V2 to V (Fact = Q/CL = V2/V; unitless, fixed)")  # Brooks 2021 Table 2 Fact = 2.0 (fixed)

    # Allometric weight exponents (fixed at theoretic values). Reference 70 kg.
    # Shared exponent on CL and Q (Q = Fact * CL inherits the CL scaling) and
    # on V and V2 (V2 = Fact * V inherits the V scaling).
    e_wt_cl_q  <- fixed(0.75); label("Allometric weight exponent on CL and Q (unitless, fixed theoretic)")  # Brooks 2021 Results: theoretic 0.75; estimated 0.73 then fixed
    e_wt_vc_vp <- fixed(1.00); label("Allometric weight exponent on V and V2 (unitless, fixed theoretic)")  # Brooks 2021 Results: theoretic 1.0; estimated 0.83 then fixed

    # ----- Inter-individual variability (Brooks 2021 Table 2) -----
    # IIV on CL only. Reported as 26.1% CV; omega^2 = log(CV^2 + 1)
    # = log(0.261^2 + 1) = log(1.06812) = 0.06591.
    # The paper's Results narrative also says BSV on V was evaluated and
    # included in the final model, but Table 2 reports no IIV magnitude for
    # V (column entry "-"); the in-file model therefore follows Table 2 and
    # carries an IIV term on CL alone. See vignette's Assumptions and
    # deviations section for this point.
    etalcl ~ 0.06591  # Brooks 2021 Table 2 IIV CL = 26.1% CV -> log(1 + 0.261^2) = 0.06591

    # Inter-occasion variability (IOV) on CL was reported in Brooks 2021
    # Table 2 (28.7% CV) but is NOT encoded structurally here: nlmixr2lib
    # model files target a single subject-level eta per parameter, and the
    # paper does not define what "occasion" means operationally (no IOV
    # column / per-dose-change occasion definition is given). Downstream
    # users who want to simulate IOV can add an OCC indicator and a per-
    # occasion eta in rxode2. See vignette's Assumptions and deviations
    # section for this point.

    # ----- Residual error (Brooks 2021 Table 2 / Results) -----
    # Proportional residual error 17.9% (additive component could not be
    # identified from the data and was dropped, Results).
    propSd <- 0.179; label("Proportional residual error (fraction)")  # Brooks 2021 Results: proportional residual = 17.9% (additive could not be identified and was dropped)
  })

  model({
    # ----- 1. Derived covariate terms -----
    # Allometric weight scaling on all PK parameters with reference 70 kg.
    # Shared exponents: 0.75 on CL and Q; 1.0 on V and V2.
    wt_cl_q  <- (WT / 70)^e_wt_cl_q
    wt_vc_vp <- (WT / 70)^e_wt_vc_vp

    # CYP3A4/5-inhibitor effect on CL (power-of-binary-indicator form).
    azole_cl <- e_azole_cl^CONMED_AZOLE

    # ----- 2. Individual PK parameters -----
    cl <- exp(lcl + etalcl) * wt_cl_q * azole_cl
    vc <- exp(lvc) * wt_vc_vp
    q  <- fact_q_vp * exp(lcl + etalcl) * wt_cl_q * azole_cl
    vp <- fact_q_vp * exp(lvc) * wt_vc_vp

    # ----- 3. Micro-constants -----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ----- 4. ODE system (2-cmt IV continuous-infusion disposition) -----
    # No depot: dose is given as IV continuous infusion directly into the
    # central compartment.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ----- 5. Observation and error -----
    # Dose in mg, central amount in mg, vc in L -> mg/L; multiply by 1000 to
    # report tacrolimus whole-blood concentration in ng/mL, the units used
    # throughout Brooks 2021.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
