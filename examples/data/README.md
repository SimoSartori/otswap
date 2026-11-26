# Example data

These files are provided for running the examples in this directory's sibling
folders, and may be redistributed unchanged together with them. They are not
covered by the licence of the otswap code, and may not be used for any other
purpose.

| File | Content | Rows | Columns and units |
|---|---|---|---|
| `box_halos_real_space.dat` | halos at redshift 1 in a cube of side 400 Mpc/h, [269, 669) on each axis, real space | 50,547 | x y z, Mpc/h |
| `box_halos_redshift_space.dat` | the same halos, row by row, in redshift space with the line of sight along y, in the same cube | 50,547 | x y z, Mpc/h |
| `box_randoms.dat` | randoms uniform in the cube, 8 per halo | 404,376 | x y z, Mpc/h |
| `periodic_halos.dat` | halos at redshift 0 in a periodic simulation cube of side 1000 Mpc/h, [0, 1000) on each axis: the real-space position, then z in redshift space with the line of sight along z, wrapped into the box (flat ΛCDM, omega_b = 0.02244, omega_c = 0.1206, h = 0.67; linear bias 1.14) | 50,000 | x y z z_s, Mpc/h |
| `lightcone_tracers.dat` | tracers in a lightcone, redshift 0.885 to 1.10 | 29,579 | RA Dec (degrees), z |
| `lightcone_randoms.dat` | randoms with the same footprint and redshift distribution | 242,548 | RA Dec (degrees), z |
| `lightcone_tracers_cartesian.dat` | the same tracers, row by row, in comoving Cartesian coordinates (flat ΛCDM, w0 = −1, wa = 0, Omega_m = 0.286, h = 0.7; observer at the origin) | 29,579 | x y z, Mpc/h |
| `lightcone_randoms_cartesian.dat` | the same randoms, in the same coordinates | 242,548 | x y z, Mpc/h |
| `lightcone_mask.fits` | the footprint: a full-sky HEALPix map, NSIDE 1024, RING, 1 observed and 0 unobserved | 12,582,912 pixels | |
| `lightcone_bias.dat` | the tracers' linear bias, five redshift nodes from 0.600 to 1.114 | 5 | z b |

The ASCII files are space separated, with one header line starting with `#`.
