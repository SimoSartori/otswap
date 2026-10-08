# The method

## The problem

A catalogue of $N$ tracers is observed at positions $\mathbf{x}_i$, the
Eulerian positions. A set of $N$ points, the randoms $\mathbf{y}_j$, is drawn
from a random catalogue that samples the same volume uniformly; they stand for
the homogeneous initial conditions. Pairing each tracer with one random, so that the sum of
the squared distances

$$
C(\sigma) = \sum_{i=1}^{N} \left| \mathbf{x}_i - \mathbf{y}_{\sigma(i)} \right|^2
$$

is minimal over the permutations $\sigma$, is the discrete form of the optimal
transport problem with the Monge cost function (Monge 1781): the total squared
displacement under one-to-one assignment (Brenier et al. 2003). Each pair
defines a displacement, from the tracer's observed position to its partner.
The tracer's displacement is the mean of these vectors over independent
realizations (below).

This is a reconstruction because of the Lagrangian picture of structure
formation. Within the Zel'dovich approximation (Zel'dovich 1970) the
displacement from the initial positions is irrotational and derives from the
gradient of a scalar potential, which guarantees a unique mapping before
shell crossing (Brenier 1991; Brenier et al. 2003). The optimal transport
problem with quadratic cost is a convex optimization problem with a unique
solution under mild regularity assumptions (Brenier 1991), and the
reconstruction holds as long as the large-scale motions are predominantly
irrotational and single-stream. Villani (2008) reviews optimal transport;
Brenier (1991), Benamou & Brenier (2000), Frisch et al. (2002), Brenier et
al. (2003) and Nikakhtar et al. (2022, 2023, 2024) introduce it and some of
its cosmological applications.

## How otswap solves it

Solving the full variational problem is demanding. otswap builds instead on
the pairwise-swapping concept of the Particle Interchange Zel'dovich
Approximation, PIZA (Croft & Gaztañaga 1997), further refined by Elyiv et al.
(2015): pairs are swapped so as to reduce the total cost. The scheme of
Sartori et al. (2026), developed independently and implemented in
CosmoBolognaLib (Marulli et al. 2016), swaps quartets rather than pairs. It
works in two phases.

1. **First guess.** A local first guess matches each tracer to a nearby free
   random: each starting tracer first with its nearest free random, so that
   every tracer leaves this phase with a partner, then up to 31 pairs per
   starter, within 8 mean particle separations.
2. **Local search.** Each sweep visits every tracer, takes three of its 113
   nearest tracers, tries the 24 permutations of the four partners and keeps
   the best. The sweeps continue until the fraction of successful swaps in a
   sweep, swaps per tracer visited, no longer exceeds a chosen threshold,
   `convergence`.

This is 4-opt local search on an assignment problem: a heuristic, not an
exact solver. The result is a good assignment, not a certified optimum. The
cost is fixed to the squared Euclidean distance in three dimensions; otswap is
not a general optimal transport solver and takes no cost matrix.

## Realizations

The reconstruction is repeated on independent sets of $N$ points extracted
from the random catalogue. A tracer's displacement is the mean, over the
realizations valid for it, of the vector from the tracer to its partner; the
mask filter below can make a realization invalid for a tracer. As the number
of realizations grows, the mean Lagrangian position tends to the centroid of
the equal-volume Laguerre cell that optimal transport assigns to the tracer:
the Lagrangian volume from which the observed object has evolved.

## The lightcone

A survey gives sky coordinates. They are converted to comoving Cartesian
coordinates with a redshift to distance table. The mean particle separation
changes with redshift, so otswap measures it from the tracers, in redshift
bins, and scales the seeding and the search with it.

The survey's footprint is a HEALPix mask. Tracers and randoms on unobserved
pixels are left out, and a displacement whose great-circle arc, from the
tracer to its random, crosses unobserved pixels is rejected: the arc is
bisected adaptively, and the displacement is rejected when it crosses more
than a set number of distinct unobserved pixels.

## Cosmological use

A tracer's mean displacement points to the centroid of the Lagrangian volume
it evolved from, so the displacement field estimates the large-scale motion
that carried each tracer from the initial state to where it is observed.
otswap provides
two uses of it.

### The input of a void finder

The displacement field is the input of the Back-in-time Void Finder, which
identifies voids dynamically from it (Sartori et al. 2026).

### Redshift-space correction

A tracer's redshift includes its peculiar velocity along the line of sight
$\hat{\mathbf{l}}$, so its redshift-space position $\mathbf{s}$ differs from
its real-space one $\mathbf{r}$ (Kaiser 1987). In linear theory the velocity
is $\mathbf{v} = a f H \boldsymbol{\Psi}$, with $\boldsymbol{\Psi}$ the
displacement from the Lagrangian position and $f$ the linear growth rate, and

$$
\mathbf{r} = \mathbf{s} - f \, \Psi_\parallel \, \hat{\mathbf{l}},
\qquad \Psi_\parallel = \boldsymbol{\Psi} \cdot \hat{\mathbf{l}} .
$$

The reconstruction does not give $\boldsymbol{\Psi}$. In first-order
Lagrangian theory $\nabla \cdot \boldsymbol{\Psi} = -\delta_m$, so in Fourier
space, with $\mu = \hat{\mathbf{k}} \cdot \hat{\mathbf{l}}$,

$$
\boldsymbol{\Psi}(\mathbf{k}) = \frac{i \mathbf{k}}{k^2} \, \delta_m(\mathbf{k}),
\qquad
\Psi_\parallel(\mathbf{k}) = \frac{i \mu}{k} \, \delta_m(\mathbf{k}) .
$$

The tracers are observed in redshift space and are biased, with linear bias
$b$: their density contrast is $\delta_g^s(\mathbf{k}) = (b + f\mu^2)\,
\delta_m(\mathbf{k})$ (Kaiser 1987). The optimal transport displacement is
irrotational (Brenier 1991), and in the linear regime its divergence is set
by the density of the tracers it moves, $\nabla \cdot
\boldsymbol{\Psi}_\mathrm{OT} = -\delta_g^s$ (Sartori et al. 2026,
Appendix C). Along the line of sight, mode by mode,

$$
\Psi_{\mathrm{OT},\parallel}(\mathbf{k}) = (b + f\mu^2) \, \Psi_\parallel(\mathbf{k}),
\qquad
f \, \Psi_\parallel(\mathbf{k}) = \frac{f}{b + f\mu^2} \, \Psi_{\mathrm{OT},\parallel}(\mathbf{k}) .
$$

The correction is therefore mode dependent. Sartori et al. (2026) take
$f/(b + f)$ (their eq. C.22, eq. 21 in the text), which corresponds to
$\mu = 1$. otswap uses a single factor
averaged over the modes. The line-of-sight displacement of a mode is
proportional to $\mu$, so its power to $\mu^2$; the factor that relates the
two line-of-sight components on average is the ratio of their
cross-correlation to the power of $\Psi_\parallel$, which for an isotropic
field reduces to an average over directions:

$$
\frac{\langle \Psi_{\mathrm{OT},\parallel} \, \Psi_\parallel^* \rangle}
     {\langle |\Psi_\parallel|^2 \rangle}
= b + f \, \frac{\langle \mu^4 \rangle}{\langle \mu^2 \rangle}
= b + f \, \frac{1/5}{1/3}
= b + \frac{3}{5} f .
$$

This refines eq. C.22 and approximates the mode-dependent correction. The
shift of tracer $i$ along its line of sight is

$$
s_i = \frac{f(z_i)}{b(z_i) + 3 f(z_i)/5} \,
      \langle \boldsymbol{\Psi}_\mathrm{OT} \cdot \hat{\mathbf{r}} \rangle_i ,
$$

where otswap's displacement points from the observed to the reconstructed
position, backward in time, which accounts for the sign, and
$\langle\,\rangle_i$ is a gaussian average of width $\sigma$ over the tracers
around $i$, which keeps the correction on the scales where linear theory
holds. In a lightcone the line of sight is each tracer's own direction, and
the tracer moves to the redshift of comoving distance $d(z_i) + s_i$; in a box
it is a Cartesian axis. $f$ comes from the distance table and $b(z)$ from a
table of the tracers' bias.

## Parameters

The choice of the parameters is discussed in Appendix B of Sartori et al.
(2026). The guide lists the ones that matter and what each one changes.

## References

- Benamou, J.-D. & Brenier, Y. 2000, Numerische Mathematik, 84, 375
- Brenier, Y. 1991, Communications on Pure and Applied Mathematics, 44, 375
- Brenier, Y., Frisch, U., Hénon, M., et al. 2003, MNRAS, 346, 501
- Croft, R. A. C. & Gaztañaga, E. 1997, MNRAS, 285, 793
- Elyiv, A., Marulli, F., Pollina, G., et al. 2015, MNRAS, 448, 642
- Frisch, U., Matarrese, S., Mohayaee, R., & Sobolevski, A. 2002, Nature, 417,
  260
- Kaiser, N. 1987, MNRAS, 227, 1
- Marulli, F., Veropalumbo, A., & Moresco, M. 2016, Astronomy and Computing,
  14, 35
- Monge, G. 1781, Histoire de l'Académie Royale des Sciences de Paris, 666
- Nikakhtar, F., Sheth, R. K., Lévy, B., & Mohayaee, R. 2022, Phys. Rev. Lett.,
  129, 251101
- Nikakhtar, F., Padmanabhan, N., Lévy, B., Sheth, R. K., & Mohayaee, R. 2023,
  Phys. Rev. D, 108, 083534
- Nikakhtar, F., Sheth, R. K., Padmanabhan, N., Lévy, B., & Mohayaee, R. 2024,
  Phys. Rev. D, 109, 123512
- Sartori, S., et al. 2026, *The Back-in-time Void Finder: dynamical
  identification of cosmic voids through optimal transport reconstruction*,
  [arXiv:2601.15378](https://arxiv.org/abs/2601.15378)
- Villani, C. 2008, *Optimal Transport: Old and New*, Grundlehren der
  mathematischen Wissenschaften, 338 (Berlin: Springer)
- Zel'dovich, Y. B. 1970, A&A, 5, 84
