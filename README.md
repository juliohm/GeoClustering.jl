# GeoClustering.jl

[![][build-img]][build-url] [![][codecov-img]][codecov-url]

Geostatistical clustering methods for the [GeoStats.jl](https://github.com/JuliaEarth/GeoStats.jl) framework.

### SLIC

Simple Linear Iterative Clustering (SLIC) produces clusters
that are spatially connected based on a geospatial distance `dₛ`.
The samples in these clusters are similar in terms of their features
according to a distance `dᵥ`. The tradeoff is controlled with a
hyperparameter parameter `m` in an additive model `dₜ = √(dᵥ² + m²(dₛ/s)²)`.
The original method developed for images is described in
[Achanta et al. 2011](https://ieeexplore.ieee.org/document/6205760).
It has been generalized in this package for any geospatial data set
(e.g. point sets).

### GHC

Geostatistical Hierarchical Clustering (GHC) produces clusters
based on (cross-)variograms between covariates and on a kernel
function between geospatial coordinates. The method is described
in [Fouedjio, F. 2016](https://www.sciencedirect.com/science/article/abs/pii/S2211675316300367).

### GSC

Geostatistical Spectral Clustering (GSC) produces clusters
based on the spectral decomposition of the graph Laplacian
constructed with weights that are a function of features and
locations. The method is described in
[Romary et al. 2015](https://www.sciencedirect.com/science/article/pii/S0098300415001314).

## Installation

Get the latest stable release with Julia's package manager:

```julia
] add GeoClustering
```

## Usage

This package is part of the [GeoStats.jl](https://github.com/JuliaEarth/GeoStats.jl) framework.

For a simple example of usage, please check the main documentation.

## Asking for help

If you have any questions, don't hesitate to ask in our community channels:

[![GITTER][gitter-img]][gitter-url]
[![ZULIP][zulip-img]][zulip-url]

[build-img]: https://img.shields.io/github/workflow/status/JuliaEarth/GeoClustering.jl/CI?style=flat-square
[build-url]: https://github.com/JuliaEarth/GeoClustering.jl/actions

[codecov-img]: https://img.shields.io/codecov/c/github/JuliaEarth/GeoClustering.jl?style=flat-square
[codecov-url]: https://codecov.io/gh/JuliaEarth/GeoClustering.jl

[gitter-img]: https://img.shields.io/badge/chat-on%20gitter-bc0067?style=flat-square
[gitter-url]: https://gitter.im/JuliaEarth/GeoStats.jl

[zulip-img]: https://img.shields.io/badge/chat-on%20zulip-9cf?style=flat-square
[zulip-url]: https://julialang.zulipchat.com/#narrow/stream/276201-geostats.2Ejl
