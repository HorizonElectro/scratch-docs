# scratch-docs

## Introduction

This is public documentation for ```scratch```, an efficient EFIE solver accelerated by MLFMA for arbitrary PEC geometries with a convenient programmable interface.

## Example

Define your geometry using pyvista. Using standard pyvista objects, define the incident E fields, and then use ```solve``` to solve for the unkown surface currents.

```
from scratch import GaussianBeam, Mirror, solve
from scipy.constants import milli

lam = 3 * milli

gb = GaussianBeam(w0=10 * milli, lam=lam, num_waist=1.5, num_lam=2)

m1 = Mirror(L=15 * milli, d=1 * milli, f=20 * milli)
m1 = m1.rotate_x(135)
m1 = m1.translate([0, 0, 50 * milli])

tri = solve(lam=lam, pec=m1 + gb)
```

the solved triangulation of surface currents can be used to solve for probe fields in free-space

```
from scratch import ds2grid, post

vol = ds2grid(tri, d=1 / 2 * milli, scale=1.5)

post(lam=lam, source=tri, probe=vol)
```

the resulting obects can be visualized

```
vol.plot(volume=True)
```

## Contact

For questions and licensing options, contact alex@horizonelectro.com