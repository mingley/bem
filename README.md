# bem

A Rust port of the [EnergyPlus](https://github.com/NatLabRockies/EnergyPlus) whole-building energy simulation engine. This is not an official EnergyPlus release. The EnergyPlus name is a U.S. Department of Energy trademark; the program in this repository is `bem`. Portions are derived from EnergyPlus `develop` at `ae0eaba` and keep that license. See `LICENSE.txt`.

The port is incomplete. What runs today is enough to simulate a zone with ideal loads. The rest of the engine (air loops, plant, airflow networks, the full window stack, CTF generators, EMS) is still ahead.

## Run

```bash
cargo run --release -- run examples/box.idf examples/weather.epw --out results.csv
cargo run --release -- run examples/box.epJSON examples/weather.epw --out results.csv
```

`examples/weather.epw` is three days of synthetic January weather, not a measured climate file. The CSV is hourly. The process also prints heating and cooling in kWh.

```bash
cargo test
```

## What is ported

Translated from the EnergyPlus sources, with the same equations:

- Psychrometrics (`Psychrometrics.cc` / `.hh`), including Hyland–Wexler saturation pressure and the wet-bulb iterator
- `General::Iterate`
- Daily solar coefficients and `SUN4` sun direction
- IDF field grammar, epJSON (current IDD field names), and EPW hourly records
- Relative-to-world vertex transform
- Third-order zone air temperature and humidity correctors, and the ideal-loads setpoint load
- Schedule:Constant and Schedule:Compact

Filled in so a zone can run, and labeled where they are not yet the EnergyPlus model:

- Opaque walls use fully implicit finite difference, not conduction transfer functions
- Inside convection is the TARP natural-convection model
- Outside convection is the DOE-2 model (TARP natural plus MoWiTT forced, with the EnergyPlus roughness multipliers)
- Windows use the simple-glazing block model. Transmitted beam is absorbed by floors; the rest of the transmitted solar is spread by area times inside absorptance
- Solar is isotropic beam plus diffuse and a 0.2 ground albedo, with no shading
- Opaque surfaces exchange long-wave radiation with the Script F method and approximate view factors
- Ideal loads meet sensible setpoint only. Dehumidification, outdoor-air minimums, and flow limits are not applied
- `ZoneHVAC:Baseboard:Convective:Electric` meets the remaining heating-setpoint load up to its capacity. Electricity is delivered heat divided by efficiency. Equipment runs in SequentialLoad heating-sequence order. Capacity per floor area is applied. Autosize is rejected
- `Output:Variable` for zone temperature, humidity ratio, ideal-loads heating and cooling rates, baseboard heating and electricity rates, and outdoor dry-bulb. With no output request, the CSV is the full hourly summary

## Layout

```text
src/baseboard.rs          convective electric baseboard
src/convection.rs        TARP interior and DOE-2 exterior films
src/psychrometrics.rs   moist-air properties
src/solar.rs            declination, equation of time, sun direction
src/idf.rs              IDF reader
src/epw.rs              EPW reader
src/model.rs            IDF objects this slice understands
src/simulate.rs         timestep, conduction, ideal loads
src/main.rs             bem run
```
