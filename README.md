# Introduction

```
Surface Current Radiation Analysis Toolkit; Computational Helmholtz
```

![](./assets/iso.PNG)

`scratch` solves electromagnetic scattering and radiation from arbitrary PEC surfaces, two ways:

- **EFIE / Method of Moments** — Rao–Wilton–Glisson (RWG) edge basis functions, solving the
  electric field integral equation directly for surface currents [1].
- **Cascading physical optics** — fast, approximate surface-to-surface propagation of induced
  currents via the Stratton-Chu equivalence-principle radiation integral [2, 3], for
  quasi-optical systems (mirror/launcher trains) where a full MoM solve isn't needed.

# References

1. S. M. Rao, D. R. Wilton, and A. W. Glisson, "Electromagnetic scattering by surfaces of
   arbitrary shape," *IEEE Trans. Antennas Propag.*, AP-30(3), 1982.
2. J. A. Kong, *Electromagnetic Wave Theory*. EMW Publishing, 2000.
3. C. A. Balanis, *Antenna Theory: Analysis and Design*, 4th ed. Wiley, 2016.

# Licensing

Contact alex@horizonelectro.com for licensing options or a trial.
