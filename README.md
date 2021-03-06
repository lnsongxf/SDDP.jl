# SDDP

[![Build Status](https://travis-ci.org/odow/SDDP.jl.svg?token=BjRx6YCjMdN19LP812Rj&branch=master)](https://travis-ci.org/odow/SDDP.jl)

## Installation
This package is unregistered so you will need to `Pkg.clone` it as follows:
```julia
Pkg.clone("https://github.com/odow/SDDP.jl.git")
```

## Documentation

The documentation is still very incomplete, and the internals of the library need a tidy and a refactor, however the user-facing API from the examples should be stable enough to use.

However, you can find the documentation at https://odow.github.io/SDDP.jl/build/index.html

## SDDiP

[@lkapelevich](https://github.com/lkapelevich) wrote an extension for SDDP.jl to solve multi-stage stochastic programs with binary state variables. Check it out at https://github.com/lkapelevich/SDDiP.jl!

## Quick Start Guide
For now the best documentation is probably contained in the examples. There is
quite a few and they provide a fairly comprehensive overview of the library.

### A Note on Value Functions

You may notice we parameterise the SDDPModel by the DefaultValueFunction. Although
this is the only value function provided in this package, it enables extensibility
for some of our research codes that are not yet at the point for public release.
