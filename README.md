# QG4L
Quasi Geostrophic 4-layer toy model

Code developed for the simulations presented in the paper "Semi-analytical approach to study the role of abyssal stratification in the propagation of potential vorticity in a four-layer ocean basin" submitted to OS (Ocean Science). The code describes potential vorticity advection using Quasi Geostrophic (QG) formalism [1, 2, 3], with a 4-layer scheme based on in-situ observation (data are available at http://doi.org/10.5281/zenodo.7871735) [4]. The model equations are derived with a z-coordinate scheme for the derivative discretization instead of the density coordinate as is done typically in QG layered models, to account for the observed density jumps [5, 6].

Contents:
- QG4L.m			main 
- QG4L_evolution.m	core function
- QGeqs.m			QG equations 
- inverse.m		Gauss-Schiedel for laplacian inversion
- arakawa.m		Arakawa scheme for Jacobian
- laplaciancr.m		Laplacian discretization
- diricheletbc.m		Boundary Condition
- periodic.m		Boundary Condition
- periodicx.m		Boundary Condition

Needed for running as it is:
- TEOS-10 functions http://www.teos-10.org/
- Crameri's Scientific colour maps https://www.fabiocrameri.ch/colourmaps/

References:
[1] Pedlosky, J. (2013). Geophysical fluid dynamics. Springer Science & Business Media. 
[2] Cushman-Roisin, B., & Beckers, J. M. (2011). Introduction to geophysical fluid dynamics: physical and numerical aspects. Academic press.
[3] Vallis, G. K. (2017). Atmospheric and oceanic fluid dynamics. Cambridge University Press.
[4] Giambenedetti, B., Lo Bue, N., Kokoszka, F., Artale, V., and Falcini, F. (2023). Multiapproach analysis of baroclinic internal tide perturbation in the Ionian Sea abyssal layer (Mediterranean Sea). Geophysical Research Letters, 50, e2023GL104311. https://doi.org/10.1029/2023GL104311
[5] Benzi, R., Pierini, S., Vulpiani, A., & Salusti, E. (1982). On nonlinear hydrodynamic stability of planetary vortices. Geophysical & Astrophysical Fluid Dynamics, 20(3-4), 293-306. https://doi.org/10.1080/03091928208213657
[6] Griffies et al., 2000).
