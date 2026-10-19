//! Convective electric baseboard from `BaseboardElectric::SimElectricConvective`.
//!
//! The unit meets the remaining sensible load to the heating setpoint, up to
//! its nominal capacity, when the availability schedule is positive and that
//! load is at least `HVAC::SmallLoad` (1 W). Electricity is the delivered heat
//! divided by efficiency. Efficiency does not reduce the heat added to the zone.

/// `HVAC::SmallLoad` in DataHVACGlobals.hh, watts.
pub const SMALL_LOAD_W: f64 = 1.0;

/// Returns `(delivered_heat_w, electricity_w)`.
pub fn convective_electric(
    load_to_heating_setpoint_w: f64,
    capacity_w: f64,
    efficiency: f64,
    available: bool,
) -> (f64, f64) {
    if available && load_to_heating_setpoint_w >= SMALL_LOAD_W && efficiency > 0.0 {
        let power = load_to_heating_setpoint_w.min(capacity_w.max(0.0));
        (power, power / efficiency)
    } else {
        (0.0, 0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn capacity_limits_heat_and_efficiency_sets_electricity() {
        let (heat, elec) = convective_electric(800.0, 500.0, 0.8, true);
        assert!((heat - 500.0).abs() < 1e-12);
        assert!((elec - 625.0).abs() < 1e-12);

        let (heat, elec) = convective_electric(100.0, 500.0, 0.5, true);
        assert!((heat - 100.0).abs() < 1e-12);
        assert!((elec - 200.0).abs() < 1e-12);
    }

    #[test]
    fn small_load_and_unavailable_schedule_stay_off() {
        let (heat, elec) = convective_electric(0.999, 500.0, 1.0, true);
        assert_eq!(heat, 0.0);
        assert_eq!(elec, 0.0);

        let (heat, elec) = convective_electric(1.0, 500.0, 1.0, true);
        assert!((heat - 1.0).abs() < 1e-12);
        assert!((elec - 1.0).abs() < 1e-12);

        let (heat, elec) = convective_electric(400.0, 500.0, 1.0, false);
        assert_eq!(heat, 0.0);
        assert_eq!(elec, 0.0);
    }
}
