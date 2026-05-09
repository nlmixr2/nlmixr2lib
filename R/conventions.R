#' nlmixr2lib convention standards
#'
#' Internal register of canonical parameter, compartment, covariate, and
#' residual-error names used by [checkModelConventions()]. The canonical
#' covariate list is parsed at runtime from
#' `inst/references/covariate-columns.md` (installed as
#' `system.file("references", "covariate-columns.md", package = "nlmixr2lib")`),
#' so that register remains the single authoritative source. Remaining fields
#' mirror the `extract-literature-model` skill's `naming-conventions.md` and
#' `vignettes/create-model-library.Rmd`.
#'
#' @keywords internal
#' @noRd
.nlmixr2libConventionsStatic <- list(
  pkParams = c(
    "lka", "lcl", "lvc", "lvp", "lvp2", "lq", "lq2", "lfdepot",
    "lvmax", "lcl_ss", "lcl_time", "lcl_renal", "lcl_nonren"
  ),
  pkBareParams = c(
    "ka", "cl", "vc", "vp", "vp2", "q", "q2", "kel",
    "k12", "k21", "k13", "k31", "fdepot",
    "vmax", "cl_ss", "cl_time", "cl_renal", "cl_nonren"
  ),
  compartments = c(
    "depot", "central", "peripheral1", "peripheral2", "effect",
    "target", "complex", "total_target",
    # Semi-physiological liver compartment used by paper-specific
    # extraction-ratio first-pass models (Xie_2019_agomelatine).
    "liver",
    # Cumulative-hazard state used by time-to-event / dropout sub-models
    # (Girard_2012_pimasertib). The state integrates the instantaneous
    # hazard so that survival = exp(-cumhaz); the source NONMEM idiom is
    # `$MODEL COMP=(CUMHAZ)` with `DADT(<cumhaz>) = HAZARD`.
    "cumhaz",
    # Renal-cortex accumulation compartment used by aminoglycoside
    # nephrotoxicity models (Llanos-Paez_2017_gentamicin). Tracks drug
    # amount sequestered in the renal cortex via saturable uptake from
    # the central compartment plus first-order tubular reabsorption back
    # out (Rougier 2003 / Croes 2011 mechanism).
    "renal_cortex"
  ),
  # Bare numbered chains (transit / effect / precursor / lat / dar /
  # depot) and metabolite-suffixed compartments are validated
  # separately via .matchesCompartment() so that the registered
  # metabolite list can be honored at runtime; this static regex
  # covers only the numbered-chain patterns. `depot[0-9]+` accommodates
  # parallel-absorption models with two or more depots.
  compartmentRegex = "^(transit|effect|precursor|lat|depot)[0-9]+$",
  darCompartmentRegex = "^dar[0-9]+_(central|peripheral[0-9]?)$",
  observationVar = "Cc",
  residualError = c("propSd", "addSd"),
  transformPrefixes = c("l", "logit", "probit"),
  # Covariate-effect names match e_<cov>(_<continuation>)+_<param>
  # (canonical) or have an additional trailing token for metabolite /
  # shared / CL-component (e_<cov>_<param>_<suffix>). The pattern
  # accepts up to 6 underscore-separated tokens after the leading "e_"
  # to accommodate compound covariates (RACE_BLACK, ADA_POSITIVE,
  # FORM_CHO_PHASE2). Semantic interpretation of the trailing tokens
  # is handled by .classifyCovEffect().
  covEffectPattern = "^e_[A-Za-z0-9]+(_[A-Za-z0-9]+){1,5}$",
  # Lowercase paper / payload names allowed as a third-token suffix on
  # parameters and compartments for a non-parent species. The set is
  # extended whenever a new ADC payload, target binding partner, or
  # secondary-analyte appears in a model. Keep mutually disjoint from
  # pkBareParams and clComponents to preserve covariate-effect
  # disambiguation.
  registeredMetabolites = c(
    "mmae", "dxd", "sn38", "dm4", "medm4", "mcmmaf",
    "complex", "ige", "tab", "nab",
    "dar0", "dar1", "dar2", "dar3", "dar4", "dar5", "dar6", "dar7", "dar8",
    # Small-molecule metabolites of agomelatine (Xie 2019): 3-hydroxy
    # and 7-desmethyl. Suffixes start with a digit; this is fine
    # because the convention check matches on `endsWith(name, "_<metab>")`
    # rather than treating the metabolite name itself as an R identifier.
    "3oh", "7dm",
    # N-desmethyl-bedaquiline metabolite (M2) of bedaquiline
    # (Svensson 2016 DDMODEL00000219).
    "m2",
    # Endoxifen (4-hydroxy-N-desmethyltamoxifen), the major active
    # metabolite of tamoxifen -- Ter Heine 2014.
    "endx",
    # Lidocaine sequential metabolites (DDMODEL00000281, NA_NA_lidocaine):
    # MEGX = monoethylglycinexylidide (LID -> MEGX via CYP1A2/3A4),
    # GX   = glycinexylidide (MEGX -> GX), and 2,6-XYL = 2,6-xylidide
    # (LID -> 2,6-XYL minor pathway). Each metabolite is a separate
    # central compartment with its own apparent volume in the source's
    # ADVAN5 parent + 3-metabolite structure.
    "megx", "gx", "xyl",
    # Morphine-3-glucuronide and morphine-6-glucuronide, the two major
    # glucuronide metabolites of morphine -- Knibbe 2009 DDMODEL00000248.
    "m3g", "m6g",
    # Phase-II conjugates: glucuronide (gluc) and sulphate (sulf).
    # Used for paracetamol-glucuronide / paracetamol-sulphate plasma
    # metabolite compartments in Allegaert 2015 (DDMODEL00000267).
    "gluc", "sulf",
    # Paracetamol (APAP) phase-II conjugate metabolites -- APAP-glucuronide
    # ("apapg") and APAP-sulphate ("apaps") used in the Cook 2016 newborn
    # model (DDMODEL00000271).
    "apapg", "apaps",
    # Colistin, the active polymyxin generated in vivo by hydrolysis of
    # the prodrug colistimethate sodium (CMS). Used as a metabolite
    # suffix in parent-prodrug CMS / metabolite-active-drug colistin
    # popPK models (Leuppi-Taegtmeyer 2019 DDMODEL00000295).
    "col",
    # Dihydroartemisinin, the active metabolite of artesunate
    # (Birgersson 2019 DDMODEL00000297).
    "dha",
    # Hydroxy-itraconazole (OH-ITZ), the major active metabolite of
    # itraconazole produced by CYP3A4 hydroxylation. Used as a metabolite
    # suffix in parent + metabolite simultaneous popPK models
    # (Hennig 2006 Clin Pharmacokinet 45(11):1099-1114; Hennig 2007 BJCP
    # 63(4):438-450).
    "ohi",
    # Doxorubicinol, the C-13 alcohol metabolite of doxorubicin
    # (Kunarajah 2017 paediatric oncology popPK/PD model).
    "doxol",
    # 25-O-desacetyl rifabutin, the primary active metabolite of
    # rifabutin formed by arylacetamide deacetylase (Hennig 2015
    # AAC doi:10.1128/AAC.01195-15).
    "desrbn",
    # AZ5104 (N-desmethyl osimertinib), an active EGFR-inhibitor
    # metabolite of osimertinib formed predominantly via CYP3A4/5
    # (Brown 2017 BJCP 83(6):1216-1226 doi:10.1111/bcp.13223).
    "az5104",
    # Capecitabine sequential metabolites (Urien 2005
    # doi:10.1007/s10928-005-0018-2): 5'-DFCR (5'-deoxy-5-
    # fluorocytidine, formed in the liver by carboxylesterase from
    # capecitabine), 5'-DFUR (5'-deoxy-5-fluorouridine, formed from
    # 5'-DFCR by cytidine deaminase in liver and tumour cells), and
    # 5-FU (5-fluorouracil, formed from 5'-DFUR by thymidine
    # phosphorylase preferentially in tumour tissue). Each metabolite
    # has its own central compartment with apparent volume fixed to
    # 1 L (only output rate constants K23, K34, K40 are identifiable
    # in the source NONMEM ADVAN6 fit).
    "dfcr", "dfur", "5fu",
    # AS(N-1)3' truncated antisense strand of GalNAc-conjugated
    # siRNAs (givosiran and other galnac-siRNA conjugates), formed
    # by removal of the 3'-terminal nucleotide from the antisense
    # strand. Treated as the active metabolite that is equipotent
    # with the parent in terms of RISC loading and target mRNA
    # silencing (Ayyar 2024 doi:10.1016/j.xphs.2023.10.026).
    "asn1"
  ),
  # Suffixes allowed for multi-component CL parameters. `_ss` denotes
  # the steady-state arm; `_time` denotes the time-varying decay arm.
  # `_renal` denotes the glomerular-filtration / tubular-secretion arm,
  # `_nonren` the non-renal (hepatic / metabolic / extra-renal) arm,
  # used by additive renal-plus-non-renal popPK models for renally
  # cleared small molecules (e.g. Jonckheere 2019 cefepime, where
  # CL_total = CL_renal + CL_nonren is the structural form).
  clComponents = c("ss", "time", "renal", "nonren"),
  # Paper-named mechanistic parameters that don't fit any canonical PK
  # naming pattern but recur across published models. Treated as
  # acceptable bare names (with the usual `l<name>` convention if
  # log-transformed). Add to this list rather than introducing a new
  # ad-hoc pattern.
  paperNamedParams = c(
    "kd", "kd0", "kdes", "kdecay", "krel", "kss", "kint",
    "frac", "alfm", "ksyn", "p", "vd", "kcat", "kpro", "krmr",
    # Transit-absorption naming used by published popPK models that
    # parameterize via mean-absorption-time / fraction-of-MAT
    # (Svensson 2016 bedaquiline DDMODEL00000219, Kovalenko 2020
    # dupilumab, etc.). `mat` = mean absorption time (hours / days);
    # `mtt` = mean transit time; `fr` = fraction of MAT in the transit
    # delay; `ktr` = first-order transit rate constant (= n_transit /
    # MTT for a chain of length n).
    "mat", "mtt", "fr", "ktr"
  ),
  requiredUnits = c("time", "dosing", "concentration"),
  requiredMetadata = c("description", "reference", "units"),
  deprecatedResidualError = c(
    "prop.err", "add.err", "propErr", "addErr",
    "err.prop", "err.add"
  ),
  deprecatedIivPrefixes = c("iiv_", "IIV_", "bsv_", "BSV_"),
  # Bare volume names that should be replaced with vc / vp / vp2.
  deprecatedVolumeNames = c("v", "v1", "v2", "v3", "lv", "lv1", "lv2", "lv3"),
  # Deprecated Michaelis-Menten Vmax names.
  deprecatedVmaxNames = c("vm", "lvm"),
  # Deprecated parent-suffix marker. A model that names a parent-side
  # parameter `<base>_adc` should drop the `_adc` suffix; the parent
  # uses the canonical name unsuffixed.
  deprecatedParentSuffix = "_adc"
)

.covariateRegisterCache <- new.env(parent = emptyenv())

.covariateColumnsPath <- function() {
  p <- system.file("references", "covariate-columns.md",
                   package = "nlmixr2lib")
  if (nzchar(p)) return(p)
  stop("Could not locate inst/references/covariate-columns.md in the ",
       "nlmixr2lib package. Check the installation.", call. = FALSE)
}

#' Parse the canonical covariate register from covariate-columns.md.
#'
#' Walks the Markdown register and extracts one entry per H3 heading
#' (`### NAME (**...**)`). For each entry, captures the `Units`, `Type`,
#' `Scope`, `Source aliases`, and `Example models` fields. Aliases whose
#' backticked content is not a bare R identifier (e.g. `DVID = "study1"`)
#' are skipped. Example-model tokens are accepted as backticked file names
#' ending in `.R`; the `.R` suffix is stripped so the value matches the
#' bare model function name used throughout the rest of the package.
#'
#' @param path Path to the markdown file.
#' @return A named list keyed by canonical name. Each entry is a list with
#'   `units`, `type`, `scope` (one of `"general"` / `"specific"` / `NA`),
#'   `aliases` (character vector of alias names), and `example_models`
#'   (character vector of model function names).
#' @keywords internal
#' @noRd
.parseCovariateColumns <- function(path) {
  lines <- readLines(path, warn = FALSE)
  entries <- list()
  current <- NULL
  state <- "idle"

  flush <- function() {
    if (is.null(current)) return(invisible())
    for (nm in current$names) {
      entries[[nm]] <<- list(
        units = current$units %||% "",
        type = current$type %||% "",
        scope = current$scope %||% NA_character_,
        aliases = current$aliases %||% character(),
        example_models = current$example_models %||% character()
      )
    }
  }

  aliasRegex <- "^\\s*-\\s*`([^`]+)`"
  identRegex <- "^[A-Za-z_][A-Za-z0-9_]*$"
  modelFileRegex <- "^[A-Za-z_][A-Za-z0-9_-]*\\.R$"

  extractBacktickedModels <- function(text) {
    toks <- regmatches(text, gregexpr("`([^`]+)`", text))[[1]]
    models <- character()
    for (tok in toks) {
      inner <- gsub("`", "", tok)
      if (grepl(modelFileRegex, inner)) {
        models <- c(models, sub("\\.R$", "", inner))
      }
    }
    models
  }

  for (line in lines) {
    if (startsWith(line, "## ") && !startsWith(line, "### ")) {
      flush()
      current <- NULL
      state <- "idle"
      next
    }
    if (startsWith(line, "### ")) {
      flush()
      heading <- sub("^###\\s+", "", line)
      heading <- sub("\\s*\\(\\*\\*.*\\*\\*\\)\\s*$", "", heading)
      nms <- trimws(strsplit(heading, ",")[[1]])
      nms <- nms[grepl(identRegex, nms)]
      current <- list(names = nms, aliases = character(),
                      example_models = character())
      state <- "header"
      next
    }
    if (is.null(current)) next

    m <- regmatches(line, regexec("^- \\*\\*Units:\\*\\*\\s*(.*)$", line))[[1]]
    if (length(m) == 2) {
      current$units <- trimws(m[[2]])
      state <- "header"
      next
    }
    m <- regmatches(line, regexec("^- \\*\\*Type:\\*\\*\\s*(.*)$", line))[[1]]
    if (length(m) == 2) {
      current$type <- trimws(m[[2]])
      state <- "header"
      next
    }
    m <- regmatches(line, regexec("^- \\*\\*Scope:\\*\\*\\s*(.*)$", line))[[1]]
    if (length(m) == 2) {
      scope_raw <- tolower(trimws(sub("\\.$", "", m[[2]])))
      if (scope_raw %in% c("general", "specific")) {
        current$scope <- scope_raw
      }
      state <- "header"
      next
    }
    if (grepl("^- \\*\\*Source aliases:\\*\\*", line)) {
      state <- "aliases"
      after <- sub("^- \\*\\*Source aliases:\\*\\*\\s*", "", line)
      # "none", "none known", "none;"-style declarations have no aliases.
      if (grepl("^none\\b", after, ignore.case = TRUE)) next
      # Capture inline aliases up to the first em-dash prose separator.
      after <- strsplit(after, "\\s+\u2014\\s+", perl = TRUE)[[1]][1]
      inline <- regmatches(after, gregexpr("`([^`]+)`", after))[[1]]
      for (tok in inline) {
        inner <- gsub("`", "", tok)
        if (grepl(identRegex, inner)) {
          current$aliases <- c(current$aliases, inner)
        }
      }
      next
    }
    if (grepl("^- \\*\\*Example models:\\*\\*", line)) {
      state <- "example_models"
      after <- sub("^- \\*\\*Example models:\\*\\*\\s*", "", line)
      current$example_models <- c(current$example_models,
                                  extractBacktickedModels(after))
      next
    }

    if (state == "aliases") {
      m <- regmatches(line, regexec(aliasRegex, line))[[1]]
      if (length(m) == 2) {
        inner <- m[[2]]
        if (grepl(identRegex, inner)) {
          current$aliases <- c(current$aliases, inner)
        }
        next
      }
      if (grepl("^- \\*\\*", line)) {
        state <- "header"
      }
    }

    if (state == "example_models") {
      # Continuation bullet lines in a multi-line Example-models list.
      if (grepl("^\\s+-\\s", line)) {
        current$example_models <- c(current$example_models,
                                    extractBacktickedModels(line))
        next
      }
      if (grepl("^- \\*\\*", line)) {
        state <- "header"
      }
    }
  }
  flush()
  entries
}

.loadCanonicalCovariates <- function(force = FALSE) {
  if (!force && !is.null(.covariateRegisterCache$canonical)) {
    return(.covariateRegisterCache$canonical)
  }
  entries <- .parseCovariateColumns(.covariateColumnsPath())
  .covariateRegisterCache$canonical <- entries
  entries
}

.nlmixr2libConventions <- function() {
  out <- .nlmixr2libConventionsStatic
  out$canonicalCovariates <- .loadCanonicalCovariates()
  out
}

.nlmixr2libCovariateAliasMap <- function() {
  conv <- .loadCanonicalCovariates()
  out <- character()
  for (canon in names(conv)) {
    for (a in conv[[canon]]$aliases) {
      out[a] <- canon
    }
  }
  out
}

# Return TRUE when `name` is a canonical log-transformed PK parameter or
# a metabolite-suffixed PK parameter (`l<base>_<metab>`).
.isPkParam <- function(name, conv) {
  if (name %in% conv$pkParams) return(TRUE)
  for (metab in conv$registeredMetabolites) {
    suf <- paste0("_", metab)
    if (endsWith(name, suf)) {
      base <- substr(name, 1, nchar(name) - nchar(suf))
      if (base %in% conv$pkParams) return(TRUE)
    }
  }
  FALSE
}

# Return TRUE when `name` is a canonical bare PK parameter or a
# metabolite-suffixed bare PK parameter (`<base>_<metab>`).
.isPkBareParam <- function(name, conv) {
  if (name %in% conv$pkBareParams) return(TRUE)
  for (metab in conv$registeredMetabolites) {
    suf <- paste0("_", metab)
    if (endsWith(name, suf)) {
      base <- substr(name, 1, nchar(name) - nchar(suf))
      if (base %in% conv$pkBareParams) return(TRUE)
    }
  }
  FALSE
}

# Compartment name validator. Recognizes:
#   - canonical names from conv$compartments
#   - numbered chains via conv$compartmentRegex (transit/effect/precursor/lat)
#   - DAR-numbered ADC isoforms via conv$darCompartmentRegex
#   - metabolite-suffixed compartments: <canonical>_<metab>
.matchesCompartment <- function(name, conv) {
  if (name %in% conv$compartments) return(TRUE)
  if (grepl(conv$compartmentRegex, name)) return(TRUE)
  if (grepl(conv$darCompartmentRegex, name)) return(TRUE)
  for (metab in conv$registeredMetabolites) {
    suf <- paste0("_", metab)
    if (endsWith(name, suf)) {
      base <- substr(name, 1, nchar(name) - nchar(suf))
      if (base %in% conv$compartments) return(TRUE)
    }
  }
  FALSE
}

# TRUE when `name` ends with `_<metab>` for any registered metabolite.
.endsWithMetabolite <- function(name, conv) {
  for (metab in conv$registeredMetabolites) {
    if (endsWith(name, paste0("_", metab))) return(TRUE)
  }
  FALSE
}

# TRUE when `name` ends with `_<component>` for any registered CL component.
.endsWithClComponent <- function(name, conv) {
  for (comp in conv$clComponents) {
    if (endsWith(name, paste0("_", comp))) return(TRUE)
  }
  FALSE
}

# TRUE when `name` ends with `_<param>` for any bare PK parameter.
# Used to detect shared-exponent covariate effects like e_wt_cl_q.
.endsWithBarePkParam <- function(name, conv) {
  for (p in conv$pkBareParams) {
    if (endsWith(name, paste0("_", p))) return(TRUE)
  }
  FALSE
}

# Classify a covariate-effect name by its trailing suffix. Returns one of:
#   "two_token"  - matches e_<cov>_<param> with no third-token suffix
#   "metabolite" - matches e_<cov>_<param>_<metab>
#   "shared"     - matches e_<cov>_<param>_<param2> (shared exponent)
#   "component"  - matches e_<cov>_<param>_<component> (multi-CL arm)
#   "unknown"    - has a third-token suffix that doesn't match any
#                  registered category
.classifyCovEffect <- function(name, conv) {
  if (!startsWith(name, "e_")) return(NA_character_)
  if (!grepl(conv$covEffectPattern, name)) return(NA_character_)
  if (.endsWithMetabolite(name, conv)) return("metabolite")
  if (.endsWithClComponent(name, conv)) return("component")
  if (.endsWithBarePkParam(name, conv)) return("shared")
  # Strip the leading `e_` and check whether the rest is a single
  # `<cov>_<param>` pair (no third-token suffix). If yes, two_token;
  # otherwise the name has an unrecognized trailing suffix.
  rest <- substr(name, 3, nchar(name))
  parts <- strsplit(rest, "_", fixed = TRUE)[[1]]
  if (length(parts) == 2) return("two_token")
  "unknown"
}
