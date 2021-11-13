# Matrix utilities for Chez Scheme

Implementation of common matrix usages / operations for Chez 9.5 as a
way to go through the rite of passage and learn some scheme. For
better or worse, I've decided to represent matrices/tensors as
expansions of vector operations. E.g., if you were to index into a 2x2
matrix, `(matrix-ref A i j)`, this would expand to `(vector-ref
(vector-ref A i) j)`. Therefore matrices are simply nested vectors of
uniform size across each axis. Because of this, these may be useful in
constructing a data-frame-esque record type with the underlying
operations being `any`-typed 2x2 matrix operations.

## Disclaimer
1. Currently very barebones! The structure isn't great!
2. I've never used scheme before so this will be potentially
not-very-scheme-like.
3. I don't know the quirks on scheme function implementations and
their computational complexity. In light of this, I've tried to 
implement slightly-better-than-naiive solutions rather than highly
optimized with the hope of going back later after benchmarking

# Lamentations
I think the scheme vector function interface `vector-*` is clunky for
matrices. What I would *like* to do is implement matrices and tensors as
functions acting on, e.g., four numbers: their indices. This is consistent with
the mathematical definition of index-based structures as simply being a map `f:
Natural Number -> Value` so that `f_i := f(i)`. This also gives a much
more natural syntax for array references, `(A i j)`. Unfortunately, this would
result in extensive computation for even very simple operations (such as
obtaining the length along each axis of `A`). This also has the downside of
obscuring the structure of a matrix from the user: a lambda could be anything!
I'm not adept enough with scheme (alt: don't know enough SRFI 

# Running
If you have [akku](https://akkuscm.org/), life is easy:

```bash
akku install
.akku/env
./tests/test-chez-matrices.sps
```

This repo only depends on [ SRFI 64
](https://srfi.schemers.org/srfi-64/srfi-64.html) for the test suite. A life
without a testing library is a sad one after-all.

# Outstanding issues
Check out the issues tab to see what isn't here yet. Most of the usual matrix
operations are implemented (although probably not very efficiently).

Complete implementations:
- Common matrix operations (`make-matrix, matrix-ref, matrix-set!, matrix-rows,
  matrix-cols, matrix-copy, matrix?, matrix-ref-row ...`)
- Common matrix operations (`matrix-min, matrix-max, mul, T
  (transpose), tr (trace), euclidean-norm ...`)
- Functions/macros to iterate over a matrix (`do-matrix, matrix-fold, matrix-map,
  matrix-contract (only 1D contractions)`).
  - `do-matrix` is a notable macro to iterate over a matrix. It allows syntax of
    the form:
    ```scheme
    (do-matrix m (i j () l) ...)
    ```
    Where `m` is the matrix to iterate over, `i, j, _, l` are indices to loop
    over. Indices can be skipped with an empty list `()` and the iteration goes
    over the entire length of an axis. 
