# Matrix utilities for Chez Scheme

Implementation of common matrix usages / operations for Chez 9.5 as a
way to go through the rite of passage and learn some scheme.

## Disclaimer
1. Currently very barebones! And not really structured well! 
2. I've never used scheme before so this will be potentially
not-very-scheme-y. 
3. I don't know the quirks on scheme function implementations and
their computational complexity. In light of this, I've tried to do
implement slightly-better-than-naiive solutions rather than highly
optimized with the hope of going back later after benchmarking

## TODOs that should be in the "Issues" tab
### Structural
* [ ] exceptions for mismatched dimensions or values
* [ ] tests (need test library?)
* [ ] benchmarking common operations (need benchmarking library?)
* [ ] think harder about naming convention for more famous functions
      (e.g the frobenius norm: should it be `frobenius-norm`?
      `array-frobenius-norm`? kind of long-winded!)
* [ ] use library or module? are these different in R6RS? i'm new here
* [ ] Separate modules for matrix functionalities (base, exactly
      solving equations, approximately solving equations, operations
      on matrices, linear approximations, etc.)

### Implementation
easier:
* [ ] conjugate transpose (dagger)
* [ ] SVD decomposition
* [ ] solutions to linear equation
* [ ] solutions to normal equations
* [ ] regression with regularizers
* [ ] read matrix from CSV
* [ ] fixed value types? no clue how this goes in scheme


harder:
* [ ] Generalize to higher dimensions (a la numpy arrays)
* [ ] slices
* [ ] broadcasting to slices
