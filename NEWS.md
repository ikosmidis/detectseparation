# detectseparation 0.4

## New functionality

+ `detect_separation()` can now distinguish between complete and
  quasi-complete separation when `separation_type = TRUE` in
  `detect_separation_control()`.

## Other improvements, updates and additions

+ Documentation improvements.

+ Updated CI workflow.

+ `detect_separation_control()` does not infer the `separator`
  function, which is now done internally in
  `.detect_infinite_estimates()`.

+ Fixed typos in some internal variable names (thanks to [Kasper
  Daniel Hansen](https://github.com/kasperdanielhansen)).

# detectseparation 0.3

## New functionality

* New `detect_infinite_estimates()` method for the `glm` function:
  `detect_infinite_estimates()` detects infinite components in the
  maximum likelihood estimates of generalized linear models with
  binomial responses, using the methods of Schwendinger et
  al. (2021) for `"log"` links and the methods in Konis (2007) for
  '"logit"', '"probit"', and '"cauchit"' links. See
  `?detect_infinite_estimates` for details.

* `detect_separation()` is now a wrapper of `detect_infinite_estimates`.

## Other improvements, updates and additions

* Major rewrite and enrichment of help pages and examples.

* New vignette "Detecting separation and infinite estimates in log-binomial regression".


# detectseparation 0.2

## Other improvements, updates and additions

* `detect_separation()` returns a warnings if the link function is not
  one of `logit`, `probit`, `cloglog`, `cauchit`.
  
* Documentation and citations updates.

* Added new tests.

# detectseparation 0.1

* First public release.
