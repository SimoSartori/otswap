# Example data

These files are provided only for running the examples in this directory's
sibling folders. They are not covered by the licence of the otswap code, and
they may not be used, copied or redistributed for any other purpose.

| File | Content | Rows | Columns and units |
|---|---|---|---|
| `box_halos_real_space.dat` | halos at redshift 1 in a cube of side 400 Mpc/h, [269, 669) on each axis, real space | 50,547 | x y z, Mpc/h |
| `box_halos_redshift_space.dat` | the same halos, row by row, in redshift space with the line of sight along y, in the same cube | 50,547 | x y z, Mpc/h |
| `box_randoms.dat` | randoms uniform in the cube, 8 per halo | 404,376 | x y z, Mpc/h |
| `lightcone_tracers.dat` | tracers in a lightcone, redshift 0.885 to 1.10 | 29,579 | RA Dec (degrees), z |
| `lightcone_randoms.dat` | randoms with the same footprint and redshift distribution | 242,548 | RA Dec (degrees), z |
| `lightcone_mask.fits` | the footprint: a full-sky HEALPix map, NSIDE 1024, RING, 1 observed and 0 unobserved | 12,582,912 pixels | |
| `lightcone_bias.dat` | the tracers' linear bias, five redshift nodes from 0.600 to 1.114 | 5 | z b |

The ASCII files are space separated, with one header line starting with `#`.

SHA-256:

```
ab1e40e53b9a7d595bbff84073648d294b8bc7e618b3f104bd438905c168b688  box_halos_real_space.dat
66168cb28e590fc2c98a2cda23c0308af111b3bf9784d533253e11e6ca62032e  box_halos_redshift_space.dat
c92ea58393b7836a1ac7a48f60dcd9d9a414994ffa649c4f249815a60bb06512  box_randoms.dat
7c4e77be0f10af7087ba4c3f7866beac5e3074fd48311a6243a6753f7b7f1a2d  lightcone_tracers.dat
83de91bb94588a2cbc3cdad589773b8f33cf948938bbc3319a7b5ad639e880e3  lightcone_randoms.dat
d96195ceb66fc59ecd3820ce2f36c56c3347feefbe1f5baeac78429ef8f1892f  lightcone_mask.fits
a8567446f0b99aff6b482ea757e59169fe7901041510f9a048560544a59bd5cb  lightcone_bias.dat
```
