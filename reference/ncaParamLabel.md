# Map PKNCA parameter codes to friendly human-readable labels

Vectorized recoder that maps the PKNCA parameter codes (\`PPTESTCD\`
values such as \`"cmax"\`, \`"aucinf.obs"\`, \`"half.life"\`) to the
friendly pharmacometric labels readers expect in tables and prose
(\`"Cmax"\`, \`"AUC0-Inf (obs)"\`, \`"t1/2"\`). Codes that are not in
the mapping table are returned unchanged with a warning, so unfamiliar
parameters never get silently dropped.

## Usage

``` r
ncaParamLabel(code, units = NULL)
```

## Arguments

- code:

  Character vector of PKNCA parameter codes. \`NA\` values are passed
  through unchanged.

- units:

  Optional named character vector keyed by PKNCA code. For each code
  with a matching entry, the friendly label is suffixed with the unit
  string in parentheses (e.g. \`"Cmax (ng/mL)"\`).

## Value

Character vector the same length as \`code\`.

## Author

Bill Denney

## Examples

``` r
ncaParamLabel(c("cmax", "tmax", "auclast", "aucinf.obs", "half.life"))
#> [1] "Cmax"         "Tmax"         "AUClast"      "AUC0-∞ (obs)" "t½"          
ncaParamLabel(
  c("cmax", "auclast"),
  units = c(cmax = "ng/mL", auclast = "ng*h/mL")
)
#> [1] "Cmax (ng/mL)"      "AUClast (ng*h/mL)"
```
