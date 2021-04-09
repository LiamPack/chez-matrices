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
1. Currently very barebones! And not really structured well!
2. I've never used scheme before so this will be potentially
not-very-scheme-y.
3. I don't know the quirks on scheme function implementations and
their computational complexity. In light of this, I've tried to do
implement slightly-better-than-naiive solutions rather than highly
optimized with the hope of going back later after benchmarking

# Lamentations
I think the scheme vector function interface `vector-*` is clunky for
matrices. What I would *like* to do is implement matrices and tensors
as functions acting on four numbers, their indices. This is consistent
with the mathematical definition of index-based structures as simply
being a map `f: Natural Number -> Value` so that `f_i defined as
f(i)`. This also gives a much more natural syntax for array
references, `(A i j)`. Unfortunately, this would result in extensive
computation for even very simple operations (such as obtaining the
length along each axis of `A`). If you know of a way to obtain this
kind of interface while still keeping the benefits of fixed-size
vectors, please lord let me know I'm fiending for it.

# TODOs that should be in the "Issues" tab
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
macro transition:
* [-] matrix map (this would be do-matrix over all indices)
  * 2D case OK
* [-] matrix fold
  * 2D case OK
* [x] matrix?
* [x] syntax for expanding ranges in matrix-ref (i 10 20) or something

easier:
* [ ] matrix-cols to handle single-vector as a row vector (length of
vector = cols)
* [ ] conjugate transpose (dagger)
* [ ] determinant
* [x] inverse
* could do this quick or just use the gauss-jordan reduction for
left inverse if it exists
* [ ] SVD decomposition
* [ ] generate matrix with values given by a lambda
* [x] solutions to linear equation
* [ ] read matrix from CSV
* [ ] matrix/tensor contractions (contract (i j k) (M i j) (D j k)) or something


harder:
* [ ] Generalize to higher dimensions (a la numpy arrays)
* [ ] slices
* [ ] broadcasting to slices
* [ ] fixed value types? no clue how this goes in scheme. Currently
the data type is "anything goes". probably need byte vectors.
