# Force the second absorption delay longer than the first

Reparameterizes a double-absorption model so the second path's lag is
estimated as a positive increment over the first path's lag (\`Tlag2 =
Tlag1 + diffTlag2\`), matching Monolix's "force delay2 longer than
delay1" models. Both paths must have lag times and the model must not
already be sequential.

## Usage

``` r
convertAbsForceLongerDelay(
  ui,
  central = "central",
  depot = "depot",
  depot2 = "depot2"
)
```

## Arguments

- ui:

  The model as a function (or something convertible to an rxUi object)

- central:

  central compartment name

- depot:

  depot compartment name

- depot2:

  name of the second depot compartment

## Value

a model with the second lag tied above the first

## See also

Other absorption:
[`addSecondAbsorption()`](https://nlmixr2.github.io/nlmixr2lib/reference/addSecondAbsorption.md),
[`addTransit()`](https://nlmixr2.github.io/nlmixr2lib/reference/addTransit.md),
[`addWeibullAbs()`](https://nlmixr2.github.io/nlmixr2lib/reference/addWeibullAbs.md),
[`addZeroOrderAbs()`](https://nlmixr2.github.io/nlmixr2lib/reference/addZeroOrderAbs.md),
[`convertAbsSequential()`](https://nlmixr2.github.io/nlmixr2lib/reference/convertAbsSequential.md),
[`removeTransit()`](https://nlmixr2.github.io/nlmixr2lib/reference/removeTransit.md),
[`removeZeroOrderAbs()`](https://nlmixr2.github.io/nlmixr2lib/reference/removeZeroOrderAbs.md)

## Author

Matthew L. Fidler
