# EPPIonization.jl

Lookup table for ionization from energetic particle precipitation from Earth's radiation belts.

## Usage

ionizationprofile, neutralprofiles, chargeprofiles

`ionizationprofile(altitude, energy, energydis, pitchangle, pitchdis, massdensity)`

`neutralprofiles(lat, lon, z, dt::DateTime)`

```julia
using Dates, EPPIonization
using LMPTools  # for isday, zenithangle

z = 0:110  # altitudes in km
dt = DateTime(2020, 1, 1, 2, 30)
lat, lon = 60, 258
flux = 1e5  # electron flux in el/cm²/s

# Compute table of profiles `np` with columns of height `h` in kilometers,
# electron density `Ne` in m⁻³ from [FIRITools.jl](https://github.com/fgasdia/FIRITools.jl),
# neutral temperature `Tn` in Kelvin, and densities of neutral `O`, `N2`, and `O2` in m⁻³,
# both from NRLMSISE-00.
np = neutralprofiles(lat, lon, z, dt)

# Without specifying flux, `chargeprofiles` returns only the [GPI](https://github.com/fgasdia/GPILowerIonosphere.jl)-derived background profiles.
bp = chargeprofiles(lat, lon, z, dt)
bp = chargeprofiles(np, z, daytime)

# When specifying the flux, we get the [GPI](https://github.com/fgasdia/GPILowerIonosphere.jl)-derived background and EPP-perturbed profiles.
bp, epp = chargeprofiles(flux, lat, lon, z, dt)
bp, epp = chargeprofiles(flux, np, z, daytime)
```

## License

This code is currently not distributed with a license and _cannot be distributed under other licenses_. This is a nearly direct Julia implementation of Wei Xu's `ionprofile.m` and related files. Contact Wei at: wexu6668@colorado.edu.
