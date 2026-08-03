Rodjun_2023_colistin <- function() {
  description <- "Population PK model for colistimethate sodium (CMS, the inactive prodrug) and formed colistin (the active polymyxin) in critically ill adults, as encoded by Rodjun 2023 for Monte Carlo probability-of-target-attainment simulation against carbapenem-, multidrug- and colistin-resistant Acinetobacter baumannii. CMS follows a two-compartment model with an intravenous 30-min infusion; its total clearance is the sum of a renal component proportional to creatinine clearance and a non-renal component. The entire non-renal CMS clearance is the formation pathway for colistin, which then follows a one-compartment model whose total clearance is likewise the sum of a CrCL-proportional renal component and a non-renal component. Unbound colistin concentration is returned as Ccu_col using the reported unbound fraction of 0.49, because the paper's PK/PD target is the unbound AUC ratio fAUC/MIC >= 7.4."
  reference <- paste(
    "Rodjun V, Montakantikul P, Houngsaitong J, Jitaree K, Nosoongnoen W.",
    "Pharmacokinetic/pharmacodynamic (PK/PD) simulation for dosage optimization of",
    "colistin and sitafloxacin, alone and in combination, against carbapenem-,",
    "multidrug-, and colistin-resistant Acinetobacter baumannii.",
    "Front Microbiol. 2023;14:1275909. doi:10.3389/fmicb.2023.1275909.",
    "Parameter values (Table 1) and the unbound fraction are reproduced by Rodjun 2023",
    "from Nation RL, Garonzik SM, Thamlikitkul V, Giamarellos-Bourboulis EJ, Forrest A,",
    "Paterson DL, et al. Dosing guidance for intravenous colistin in critically-ill",
    "patients. Clin Infect Dis. 2017;64(5):565-571. doi:10.1093/cid/ciw839.",
    "The differential equations (Rodjun 2023 Eqs. 1-3) are stated by Rodjun 2023 to be",
    "modified from Garonzik SM, Li J, Thamlikitkul V, Paterson DL, Shoham S, Jacob J,",
    "et al. Antimicrob Agents Chemother. 2011;55(7):3284-3294. doi:10.1128/AAC.01733-10.",
    "See also modellib('Rodjun_2023_sitafloxacin') for the companion agent.",
    sep = " "
  )
  vignette <- "Rodjun_2023_colistin_sitafloxacin"

  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance (raw, NOT BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters the model as a linear multiplier on BOTH the renal clearance of CMS (Rodjun 2023 Table 1: CLR = 0.0340 L/h per mL/min of CrCL) and the renal clearance of formed colistin (CLRC = 0.00834 L/h per mL/min of CrCL). No centering or normalization is applied - the Table 1 units 'L/h/CrCL' make the slope a clearance per mL/min, so CRCL is used in raw mL/min. Rodjun 2023 simulated the four discrete values 90, 50, 30 and 10 mL/min (Materials and methods, Simulated dosage regimens); the underlying Nation 2017 estimation cohort spanned 0-236 mL/min (Materials and methods, Pharmacokinetic model / Colistin). Stored under canonical CRCL per inst/references/covariate-columns.md, which accepts raw mL/min when the source paper does not apply BSA normalization (same convention as Karaiskos_2015_colistin.R).",
      source_name        = "CrCL"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_simulated    = 10000L,
    n_studies      = 1L,
    age_range      = "Not reported on disk. Rodjun 2023 reproduces only the parameter table of Nation 2017 and does not restate that study's baseline demographics.",
    weight_range   = "Not reported on disk. The colistin model carries no body-weight covariate, so weight does not enter the simulation.",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported on disk.",
    disease_state  = "Critically ill adults receiving intravenous colistimethate sodium (the Nation 2017 estimation cohort, described by Rodjun 2023 only as 'critically ill patients with a CrCL of 0-236 mL/min'). The simulation target population is patients with carbapenem-resistant (CRAB), multidrug-resistant (MDR-AB) or colistin-resistant (CoR-AB) Acinetobacter baumannii infection.",
    dose_range     = "Intravenous loading dose of 300 mg or 450 mg colistin base activity (CBA) infused over 30 min, followed by maintenance doses from 50 mg q48h to 450 mg q12h according to creatinine clearance (Materials and methods, Simulated dosage regimens / Colistin). All doses are expressed in mg of colistin base activity.",
    regions        = "Thailand (simulation study, Mahidol University, Bangkok); the underlying Nation 2017 PK cohort was multinational.",
    renal_function = "Simulated at creatinine clearance 90, 50, 30 and 10 mL/min. The Nation 2017 estimation cohort spanned CrCL 0-236 mL/min.",
    notes          = "Monte Carlo simulation of 10,000 virtual subjects per dosage regimen (Crystal Ball version 2017, Decisioneering Inc.), with log-normal between-patient variability on every PK parameter except the unbound fraction of colistin, which was drawn from a uniform distribution (Materials and methods, Monte Carlo simulation). Rodjun 2023 does not report the size or demographics of the Nation 2017 estimation cohort, so n_subjects is NA; n_simulated records the virtual cohort size instead. The PK/PD target is fAUC/MIC >= 7.4 (Cheah 2015 murine thigh-infection model), applied to the unbound 24-h colistin AUC. See the validation vignette for the extent to which the published PTA tables reproduce."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters - Rodjun 2023 Table 1 ("Population
    # pharmacokinetic parameters of colistin (Nation et al., 2017)").
    # The ODE system is Rodjun 2023 Eqs. 1-3.
    # ------------------------------------------------------------------

    # CMS (prodrug) disposition. Table 1 names V1 / V2 / CLD1; the
    # canonical nlmixr2lib names for a two-compartment central /
    # peripheral pair are lvc / lvp / lq.
    lvc <- log(12.9); label("CMS central volume of distribution V1 (L)")               # Rodjun 2023 Table 1: CMS V1 = 12.9 L (%SE not reported; %IIV 40.4)
    lvp <- log(16.1); label("CMS peripheral volume of distribution V2 (L)")            # Rodjun 2023 Table 1: CMS V2 = 16.1 L (%SE not reported; %IIV 70.9)
    lq  <- log(9.57); label("CMS inter-compartmental clearance CLD1 (L/h)")            # Rodjun 2023 Table 1: CMS CLD1 = 9.57 (%SE 10.5; %IIV 80.1)

    # CMS clearance is split into a CrCL-proportional renal arm and a
    # non-renal arm. Rodjun 2023 Table 1 footnote defines CLTCMS as the
    # "total intrinsic clearance for CMS"; it is not tabulated separately
    # because it is the sum CLR * CrCL + CLNR,CMS, which is what Eq. 1
    # uses. Only the non-renal arm forms colistin (Eq. 3).
    lcl_renal  <- log(0.0340); label("CMS renal clearance slope CLR (L/h per mL/min of CrCL)")  # Rodjun 2023 Table 1: CMS CLR = 0.0340 L/h/CrCL (%SE 6.85; %IIV 75.2)
    lcl_nonren <- log(2.52);   label("CMS non-renal clearance CLNR,CMS (L/h); the colistin formation pathway")  # Rodjun 2023 Table 1: CMS CLNR CMS = 2.52 L/h (%SE 3.71; %IIV 39.8)

    # Formed colistin disposition (one compartment, Eq. 3).
    lvc_col <- log(57.2); label("Colistin volume of distribution V3 (L)")              # Rodjun 2023 Table 1: Colistin V3 = 57.2 L (%SE 5.13; %IIV 43.5)

    # Colistin total clearance CLTC (Eq. 3). Table 1 reports the typical
    # total (CLTC = 3.59 L/h, %IIV 37.9) AND its two components
    # (CLRC = 0.00834 L/h per mL/min of CrCL, CLNRC = 3.11 L/h, neither
    # carrying its own %IIV). This extraction builds CLTC from the two
    # components so that colistin elimination retains the creatinine-
    # clearance dependence the simulation requires across CrCL 10-90
    # mL/min, and applies the single tabulated 37.9% IIV to the resulting
    # total. The tabulated typical CLTC of 3.59 L/h is recovered at
    # CrCL = (3.59 - 3.11) / 0.00834 = 57.6 mL/min, consistent with a
    # cohort-median creatinine clearance. See the validation vignette,
    # which shows this encoding holds the average steady-state colistin
    # concentration flat (CV 3.1%) across the four reference regimens
    # Rodjun 2023 attributes to Nation et al., whereas holding CLTC fixed
    # at 3.59 L/h does not (CV 7.3%).
    lcl_renal_col  <- log(0.00834); label("Colistin renal clearance slope CLRC (L/h per mL/min of CrCL)")  # Rodjun 2023 Table 1: Colistin CLR C = 0.00834 L/h/CrCL (%SE 27.7; no %IIV reported)
    lcl_nonren_col <- log(3.11);    label("Colistin non-renal clearance CLNRC (L/h)")                      # Rodjun 2023 Table 1: Colistin CLNR C = 3.11 L/h (%SE 4.38; no %IIV reported)

    # Unbound fraction of colistin in plasma, used to convert total
    # colistin concentration to the unbound concentration that drives the
    # paper's fAUC/MIC target. Rodjun 2023 drew this from a UNIFORM
    # distribution over 0.49 +/- 0.11 in the Monte Carlo simulation
    # (Materials and methods, Monte Carlo simulation); nlmixr2 random
    # effects are log-normal, so the typical value is fixed here and the
    # uniform draw is reproduced explicitly in the validation vignette.
    fu_col <- fixed(0.49); label("Colistin plasma unbound fraction (unitless)")  # Rodjun 2023 Materials and methods, Pharmacokinetic model / Colistin: unbound fraction 0.49 +/- 0.11 by ultracentrifugation (Nation 2017)

    # ------------------------------------------------------------------
    # Inter-individual variability. Rodjun 2023 Table 1 reports %IIV as a
    # coefficient of variation and its footnote states "standard
    # deviation, SD, were calculated from %IVV x mean" - i.e. the %IIV is
    # an ARITHMETIC CV, which Crystal Ball then used as the spread of a
    # log-normal draw. The equivalent log-scale variance is
    # omega^2 = log(CV^2 + 1).
    # ------------------------------------------------------------------
    etalvc        ~ 0.15118  # log(0.404^2 + 1); Rodjun 2023 Table 1 CMS V1 %IIV = 40.4
    etalvp        ~ 0.40725  # log(0.709^2 + 1); Rodjun 2023 Table 1 CMS V2 %IIV = 70.9
    etalq         ~ 0.49567  # log(0.801^2 + 1); Rodjun 2023 Table 1 CMS CLD1 %IIV = 80.1
    etalcl_renal  ~ 0.44802  # log(0.752^2 + 1); Rodjun 2023 Table 1 CMS CLR %IIV = 75.2
    etalcl_nonren ~ 0.14703  # log(0.398^2 + 1); Rodjun 2023 Table 1 CMS CLNR CMS %IIV = 39.8
    etalvc_col    ~ 0.17330  # log(0.435^2 + 1); Rodjun 2023 Table 1 Colistin V3 %IIV = 43.5

    # Colistin clearance IIV. Table 1 reports a single 37.9% IIV on the
    # TOTAL colistin clearance CLTC and none on its renal / non-renal
    # components, so model() applies this one eta to BOTH arms. Because
    # (CLNRC + CLRC * CrCL) * exp(eta) == CLNRC * exp(eta) + CLRC * CrCL *
    # exp(eta), the shared eta is mathematically identical to placing the
    # tabulated 37.9% IIV on the reconstructed total. It is named for the
    # non-renal arm to satisfy the eta / fixed-effect pairing convention;
    # the same shared-eta device is used in Karaiskos_2015_colistin.R.
    etalcl_nonren_col ~ 0.13421  # log(0.379^2 + 1); Rodjun 2023 Table 1 Colistin CLT C %IIV = 37.9 (shared across both colistin clearance arms)

    # ------------------------------------------------------------------
    # Residual error. Rodjun 2023 reports none: the Monte Carlo
    # simulation is deterministic once each virtual subject's PK
    # parameters and unbound fraction have been drawn (Materials and
    # methods, Monte Carlo simulation). Both residual SDs are therefore
    # fixed at 0 so the packaged model reproduces the paper's noise-free
    # simulation; a user fitting real data should re-estimate them.
    # ------------------------------------------------------------------
    propSd     <- fixed(0); label("CMS proportional residual error (fraction; 0 - not reported)")       # Rodjun 2023 Materials and methods (no residual error model reported)
    propSd_col <- fixed(0); label("Colistin proportional residual error (fraction; 0 - not reported)")  # Rodjun 2023 Materials and methods (no residual error model reported)
  })

  model({
    # 1. Individual CMS parameters. Doses are administered as a 30-min
    #    intravenous infusion of colistin base activity into `central`;
    #    the infusion rate R1 of Eq. 1 is supplied by the event table
    #    (rate / dur on the dose record) rather than by the model.
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q  <- exp(lq + etalq)

    cl_renal  <- exp(lcl_renal + etalcl_renal) * CRCL
    cl_nonren <- exp(lcl_nonren + etalcl_nonren)
    cl_cms    <- cl_renal + cl_nonren   # CLTCMS, the total intrinsic CMS clearance of Eq. 1

    # 2. Individual colistin parameters. The tabulated 37.9% IIV on total
    #    colistin clearance is applied to the reconstructed
    #    CrCL-dependent total (see the ini() comment on lcl_renal_col).
    vc_col <- exp(lvc_col + etalvc_col)
    cl_col <- (exp(lcl_nonren_col + etalcl_nonren_col) +
                 exp(lcl_renal_col + etalcl_nonren_col) * CRCL)

    # 3. ODE system - Rodjun 2023 Eqs. 1-3 verbatim.
    #    central       = CMS in the central compartment (Eq. 1, CMSc)
    #    peripheral1   = CMS in the peripheral compartment (Eq. 2, CMSp)
    #    central_col   = formed colistin in its single compartment (Eq. 3)
    #    Only the NON-RENAL CMS clearance forms colistin; the renal
    #    component leaves the system as unchanged CMS.
    d/dt(central)     <- -q * (central / vc - peripheral1 / vp) - cl_cms * (central / vc)
    d/dt(peripheral1) <-  q * (central / vc - peripheral1 / vp)
    d/dt(central_col) <-  cl_nonren * (central / vc) - cl_col * (central_col / vc_col)

    # 4. Observations. Cc is the CMS concentration in the central
    #    compartment, Cc_col the total formed-colistin concentration, and
    #    Ccu_col the unbound colistin concentration that drives the
    #    paper's fAUC/MIC >= 7.4 target.
    Cc     <- central / vc
    Cc_col <- central_col / vc_col
    Ccu_col <- Cc_col * fu_col

    Cc     ~ prop(propSd)
    Cc_col ~ prop(propSd_col)
  })
}
