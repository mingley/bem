# bem

`bem` is a Rust program that runs a whole-building energy simulation from EnergyPlus input. It reads an IDF or epJSON model and an EPW weather file, steps through the run period, and writes heating and cooling results.

This is an independent port. It is not an official EnergyPlus release. EnergyPlus is a trademark of the U.S. Department of Energy; the program in this repository is named `bem`. Derived numerical code follows the EnergyPlus license in `LICENSE.txt`. The upstream snapshot used for translation is `develop` at `ae0eaba`.

The port is in progress. A single zone with ideal loads, or ideal loads plus an electric baseboard, can run today. Air loops, plant loops, and most HVAC systems are still ahead.

## Requirements

- [Rust](https://www.rust-lang.org/tools/install) with Cargo (edition 2021)
- No extra services, accounts, or data files

## Run the example

From the repository root:

```bash
cargo run --release -- run examples/box.idf examples/weather.epw --out results.csv
```

The same model as epJSON:

```bash
cargo run --release -- run examples/box.epJSON examples/weather.epw --out results.csv
```

`examples/box.idf` is one 5 m × 5 m × 2.7 m room in Chicago, with a south window, two people, lights, infiltration, a 21°C heating setpoint, a 24°C cooling setpoint, and ideal loads. `examples/weather.epw` is three synthetic January days, not a measured climate file.

That command prints a short summary on stderr. With the default seven warmup days, the example currently reports:

```text
heating 81.239 kWh
cooling 0.000 kWh
timesteps 432
```

## Commands

```text
bem run <model.idf|model.epJSON> <weather.epw> [--out results.csv] [--warmup-days 7]
bem help
```

| Flag | Default | Meaning |
| --- | --- | --- |
| `--out` | print CSV to stdout | Write the results file |
| `--warmup-days` | `7` | Repeat the first day this many times before recording results |

Warmup lets surface and zone history settle. Use `--warmup-days 0` for a short check. A normal annual run should keep the default.

Exit code `0` means the run finished. Exit code `1` means the input could not be read or a required object is missing. The error is printed to stderr.

## Results

If the model has no `Output:Variable` objects, `bem` writes one row per hour:

| Column | Unit |
| --- | --- |
| Date | `YYYY-MM-DD` |
| Hour | `1`–`24` |
| Outdoor Drybulb | °C |
| `<Zone>` Temperature | °C |
| `<Zone>` Humidity Ratio | kg/kg |
| `<Zone>` Heating | W |
| `<Zone>` Cooling | W |

Heating in that table is all sensible heat added to the zone, including baseboard heat. Cooling is the ideal-loads cooling rate.

If the model requests `Output:Variable`, the CSV uses the EnergyPlus column style:

```text
Date/Time,KEY:Variable Name [unit](Frequency)
```

The timestamp is the end of the interval. Hour 24 is written `24:00:00`. Rates are averaged over the interval.

| Variable | Key | Unit |
| --- | --- | --- |
| Site Outdoor Air Drybulb Temperature | `Environment` or `*` | C |
| Zone Mean Air Temperature | zone name or `*` | C |
| Zone Air Humidity Ratio | zone name or `*` | kgWater/kgDryAir |
| Zone Ideal Loads Supply Air Sensible Heating Rate | zone name or `*` | W |
| Zone Ideal Loads Supply Air Sensible Cooling Rate | zone name or `*` | W |
| Baseboard Total Heating Rate | baseboard name or `*` | W |
| Baseboard Electricity Rate | baseboard name or `*` | W |

`*` expands to every matching zone or baseboard. Recognized frequencies are `Timestep`, `Hourly` (the default), `Daily`, `Monthly`, and `RunPeriod`. `Detailed` is treated as `Timestep`. A variable this version does not produce becomes a warning, and the run continues.

## What a model can contain

`bem` reads IDF and epJSON. epJSON field names follow the current EnergyPlus IDD: lowercase, with spaces and punctuation turned into `_`.

| Group | Objects |
| --- | --- |
| Run setup | `Version`, `SimulationControl`, `Building`, `Site:Location`, `RunPeriod`, `Timestep`, `GlobalGeometryRules` |
| Schedules | `ScheduleTypeLimits`, `Schedule:Constant`, `Schedule:Compact` |
| Envelope | `Material`, `WindowMaterial:SimpleGlazingSystem`, `Construction`, `Zone`, `BuildingSurface:Detailed`, `FenestrationSurface:Detailed` |
| Gains | `People`, `Lights`, `ElectricEquipment`, `GasEquipment`, `OtherEquipment`, `ZoneInfiltration:DesignFlowRate` |
| Controls | `ZoneControl:Thermostat`, `ThermostatSetpoint:DualSetpoint` |
| HVAC | `ZoneHVAC:IdealLoadsAirSystem`, `ZoneHVAC:Baseboard:Convective:Electric`, `ZoneHVAC:EquipmentList`, `ZoneHVAC:EquipmentConnections` |
| Shading | `Shading:Overhang` |
| Reports | `Output:Variable` |

An equipment type outside that HVAC list stops the run with an error. Other unrecognized objects are skipped and named in a warning.

Geometry may be in world coordinates. `BuildingSurface:Detailed` uses the current field order, including the blank Space Name between the zone and the outside boundary. A fenestration’s vertices live under `extensions` in epJSON.

`Schedule:Compact` accepts `Through`, `For`, and `Until`, including `AllDays`, `Weekdays`, `AllOtherDays`, and a named weekday. Weekdays take priority over `AllDays`.

## How the timestep works

Each zone timestep does four things:

1. Place the sun and split solar onto windows and opaque surfaces.
2. Update opaque-wall conduction with conduction transfer functions, then the inside and outside surface temperatures.
3. Solve the zone air temperature and humidity with the EnergyPlus third-order corrector.
4. Meet the thermostat with zone equipment, in `SequentialLoad` order.

Ideal loads add or remove sensible heat so the zone can hold its setpoint, up to the heating and cooling capacity entered on the object. A blank limit means unlimited capacity. Flow limits, outdoor-air minimums, and dehumidification are still ahead.

`ZoneHVAC:Baseboard:Convective:Electric` runs on its heating sequence. When its availability schedule is positive and the remaining load to the heating setpoint is at least 1 W, it delivers heat up to its capacity. Electricity is that heat divided by efficiency. A lower heating sequence runs first. Ideal loads later in the sequence cover the remaining heating, and they still do the cooling.

Enter baseboard capacity as watts (`HeatingDesignCapacity`) or as W/m² times the zone floor area (`CapacityPerFloorArea`). Autosize is rejected, because zone sizing is not in this version yet.

## Modeling choices in this version

These pieces follow the EnergyPlus equations:

- Moist-air properties, including Hyland–Wexler saturation pressure and the wet-bulb solver
- Daily solar coefficients and the `SUN4` sun direction
- Inside convection with the TARP model
- Outside convection with the DOE-2 model
- Long-wave exchange among opaque surfaces and the inside face of simple glazing, with Hottel’s Script F method
- The simple-glazing window model. Beam solar uses the angular transmittance and reflectance curves. The glazing resistance is what remains after the standard films are removed, and the glass face temperatures are solved with the simulation's films. Transmitted beam is absorbed by the floor; the rest of the transmitted solar is spread by area and inside absorptance. The inside glass face joins the long-wave enclosure
- The third-order zone air temperature and humidity updates
- Opaque-wall conduction transfer functions from the state-space method, in one dimension and without an internal heat source
- Design-flow infiltration, `Idesign * schedule * (A + B*|Tzone - Toutdoor| + C*wind + D*wind²)`, with outdoor, indoor, or standard density. Wind speed is taken from the weather file

These pieces are simpler stand-ins, so annual results will move as the port grows:

- `Shading:Overhang` shades the window or door it references. Beam solar uses the sunlit area of that surface.
- Isotropic sky diffuse uses the same shadow on each patch of a 6 by 24 sky dome. The horizon-band ratio is stored for a later anisotropic sky. Ground-reflected diffuse is not shaded.
- Fins, detached shades, and window shades are still ahead. Solar is isotropic beam plus diffuse, with a ground reflectance of 0.2.

## Tests

```bash
cargo test
```

The tests check psychrometrics, sun position, convection coefficients, long-wave balance, glazing, overhang beam and sky-diffuse shadows, schedules, ideal loads, the electric baseboard, and that the example IDF and epJSON produce the same heating and cooling.

## Layout

| File | Role |
| --- | --- |
| `src/main.rs` | `bem run` |
| `src/model.rs` | Turn IDF and epJSON into a model |
| `src/idf.rs` | IDF field reader |
| `src/epjson.rs` | epJSON field reader |
| `src/epw.rs` | EPW weather reader |
| `src/ctf.rs` | Conduction transfer functions |
| `src/simulate.rs` | Timestep, walls, zone air, and equipment |
| `src/psychrometrics.rs` | Moist-air properties |
| `src/solar.rs` | Sun position |
| `src/convection.rs` | TARP and DOE-2 surface films |
| `src/radiosity.rs` | Script F long-wave exchange |
| `src/simple_glazing.rs` | Simple window solar split |
| `src/shading.rs` | Overhang geometry, beam shadows, and sky-diffuse ratios |
| `src/baseboard.rs` | Convective electric baseboard |
| `src/report.rs` | `Output:Variable` CSV |
| `examples/` | One-zone model and three-day weather |

## License

BSD-3-Clause, with the EnergyPlus terms in `LICENSE.txt`. Keep those terms with any redistribution. Do not use the EnergyPlus name for this program.
