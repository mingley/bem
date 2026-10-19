//! EnergyPlus epJSON reader.
//!
//! An epJSON file is a JSON object whose keys are IDD object types. Each type
//! maps instance name to fields. Field names follow the schema generator in
//! `idd/schema/idd_parser.py`: lowercase, then non-alphanumeric runs become `_`.
//! Extensible groups use the names assigned in `modify_schema.py` (`vertices`,
//! `data`, `equipment`, or the default `extensions`).
//!
//! The reader emits the same positional [`IdfObject`]s the IDF path uses, in
//! current IDD order, so blank omitted fields stay aligned.

use crate::idf::IdfObject;
use serde_json::Value;

struct Spec {
    /// When true, the JSON object key is the Name field.
    named: bool,
    fields: &'static [&'static str],
    /// JSON array key for the extensible group. The first match that is an array is used.
    extensions: &'static [&'static str],
    extension_fields: &'static [&'static str],
}

pub fn parse_epjson(text: &str) -> Result<Vec<IdfObject>, String> {
    let root: Value = serde_json::from_str(text).map_err(|err| format!("epJSON: {err}"))?;
    let Some(root) = root.as_object() else {
        return Err("epJSON root must be an object".to_string());
    };
    let mut objects = Vec::new();
    for (object_type, body) in root {
        if object_type.starts_with('$') || object_type.eq_ignore_ascii_case("epJSON") {
            continue;
        }
        let Some(instances) = body.as_object() else {
            return Err(format!(
                "epJSON '{object_type}' must be an object of named instances"
            ));
        };
        let spec = spec_for(object_type);
        for (instance_name, instance) in instances {
            let Some(fields_obj) = instance.as_object() else {
                return Err(format!(
                    "epJSON {object_type} '{instance_name}' must be an object"
                ));
            };
            objects.push(match spec {
                Some(spec) => object_from_spec(object_type, instance_name, fields_obj, spec),
                None => IdfObject {
                    object_type: object_type.clone(),
                    fields: vec![instance_name.clone()],
                },
            });
        }
    }
    Ok(objects)
}

fn object_from_spec(
    object_type: &str,
    instance_name: &str,
    fields_obj: &serde_json::Map<String, Value>,
    spec: &Spec,
) -> IdfObject {
    let mut fields = Vec::new();
    if spec.named {
        fields.push(instance_name.to_string());
    }
    for key in spec.fields {
        fields.push(lookup(fields_obj, key));
    }
    if let Some(group) = extension_array(fields_obj, spec.extensions) {
        for item in group {
            let Some(item) = item.as_object() else {
                continue;
            };
            if spec.extension_fields.len() == 1 {
                fields.push(lookup(item, spec.extension_fields[0]));
            } else {
                for key in spec.extension_fields {
                    fields.push(lookup(item, key));
                }
            }
        }
    }
    IdfObject {
        object_type: object_type.to_string(),
        fields,
    }
}

fn extension_array<'a>(
    fields_obj: &'a serde_json::Map<String, Value>,
    keys: &[&str],
) -> Option<&'a Vec<Value>> {
    for key in keys {
        if let Some(value) = fields_obj
            .iter()
            .find(|(name, _)| name.eq_ignore_ascii_case(key))
            .and_then(|(_, value)| value.as_array())
        {
            return Some(value);
        }
    }
    None
}

fn lookup(object: &serde_json::Map<String, Value>, key: &str) -> String {
    object
        .iter()
        .find(|(name, _)| name.eq_ignore_ascii_case(key))
        .map(|(_, value)| scalar(value))
        .unwrap_or_default()
}

fn scalar(value: &Value) -> String {
    match value {
        Value::Null => String::new(),
        Value::Bool(flag) => flag.to_string(),
        Value::Number(number) => number.to_string(),
        Value::String(text) => text.clone(),
        Value::Array(_) | Value::Object(_) => String::new(),
    }
}

fn spec_for(object_type: &str) -> Option<&'static Spec> {
    const ZONE: Spec = Spec {
        named: true,
        fields: &[
            "direction_of_relative_north",
            "x_origin",
            "y_origin",
            "z_origin",
            "type",
            "multiplier",
            "ceiling_height",
            "volume",
            "floor_area",
        ],
        extensions: &[],
        extension_fields: &[],
    };
    const SURFACE: Spec = Spec {
        named: true,
        fields: &[
            "surface_type",
            "construction_name",
            "zone_name",
            "space_name",
            "outside_boundary_condition",
            "outside_boundary_condition_object",
            "sun_exposure",
            "wind_exposure",
            "view_factor_to_ground",
            "number_of_vertices",
        ],
        extensions: &["vertices"],
        extension_fields: &[
            "vertex_x_coordinate",
            "vertex_y_coordinate",
            "vertex_z_coordinate",
        ],
    };
    const WINDOW: Spec = Spec {
        named: true,
        fields: &[
            "surface_type",
            "construction_name",
            "building_surface_name",
            "outside_boundary_condition_object",
            "view_factor_to_ground",
            "frame_and_divider_name",
            "multiplier",
            "number_of_vertices",
        ],
        extensions: &["vertices", "extensions"],
        extension_fields: &[
            "vertex_x_coordinate",
            "vertex_y_coordinate",
            "vertex_z_coordinate",
        ],
    };
    const MATERIAL: Spec = Spec {
        named: true,
        fields: &[
            "roughness",
            "thickness",
            "conductivity",
            "density",
            "specific_heat",
            "thermal_absorptance",
            "solar_absorptance",
            "visible_absorptance",
        ],
        extensions: &[],
        extension_fields: &[],
    };
    const GLAZING: Spec = Spec {
        named: true,
        fields: &[
            "u_factor",
            "solar_heat_gain_coefficient",
            "visible_transmittance",
        ],
        extensions: &[],
        extension_fields: &[],
    };
    const CONSTRUCTION: Spec = Spec {
        named: true,
        fields: &[
            "outside_layer",
            "layer_2",
            "layer_3",
            "layer_4",
            "layer_5",
            "layer_6",
            "layer_7",
            "layer_8",
            "layer_9",
            "layer_10",
        ],
        extensions: &[],
        extension_fields: &[],
    };
    const SCHEDULE_CONSTANT: Spec = Spec {
        named: true,
        fields: &["schedule_type_limits_name", "hourly_value"],
        extensions: &[],
        extension_fields: &[],
    };
    const SCHEDULE_COMPACT: Spec = Spec {
        named: true,
        fields: &["schedule_type_limits_name"],
        extensions: &["data"],
        extension_fields: &["field"],
    };
    const PEOPLE: Spec = Spec {
        named: true,
        fields: &[
            "zone_or_zonelist_or_space_or_spacelist_name",
            "number_of_people_schedule_name",
            "number_of_people_calculation_method",
            "number_of_people",
            "people_per_floor_area",
            "floor_area_per_person",
            "fraction_radiant",
            "sensible_heat_fraction",
            "activity_level_schedule_name",
        ],
        extensions: &[],
        extension_fields: &[],
    };
    const LIGHTS: Spec = Spec {
        named: true,
        fields: &[
            "zone_or_zonelist_or_space_or_spacelist_name",
            "schedule_name",
            "design_level_calculation_method",
            "lighting_level",
            "watts_per_floor_area",
            "watts_per_person",
            "return_air_fraction",
            "fraction_radiant",
            "fraction_visible",
        ],
        extensions: &[],
        extension_fields: &[],
    };
    const EQUIPMENT: Spec = Spec {
        named: true,
        fields: &[
            "zone_or_zonelist_or_space_or_spacelist_name",
            "schedule_name",
            "design_level_calculation_method",
            "design_level",
            "watts_per_floor_area",
            "watts_per_person",
            "fraction_latent",
            "fraction_radiant",
            "fraction_lost",
        ],
        extensions: &[],
        extension_fields: &[],
    };
    const GAS: Spec = Spec {
        named: true,
        fields: &[
            "zone_or_zonelist_or_space_or_spacelist_name",
            "schedule_name",
            "design_level_calculation_method",
            "design_level",
            "power_per_floor_area",
            "power_per_person",
            "fraction_latent",
            "fraction_radiant",
            "fraction_lost",
        ],
        extensions: &[],
        extension_fields: &[],
    };
    const OTHER: Spec = Spec {
        named: true,
        fields: &[
            "fuel_type",
            "zone_or_zonelist_or_space_or_spacelist_name",
            "schedule_name",
            "design_level_calculation_method",
            "design_level",
            "power_per_floor_area",
            "power_per_person",
            "fraction_latent",
            "fraction_radiant",
            "fraction_lost",
        ],
        extensions: &[],
        extension_fields: &[],
    };
    const INFILTRATION: Spec = Spec {
        named: true,
        fields: &[
            "zone_or_zonelist_or_space_or_spacelist_name",
            "schedule_name",
            "design_flow_rate_calculation_method",
            "design_flow_rate",
            "flow_rate_per_floor_area",
            "flow_rate_per_exterior_surface_area",
            "air_changes_per_hour",
            "constant_term_coefficient",
            "temperature_term_coefficient",
            "velocity_term_coefficient",
            "velocity_squared_term_coefficient",
            "density_basis",
        ],
        extensions: &[],
        extension_fields: &[],
    };
    const THERMOSTAT: Spec = Spec {
        named: true,
        fields: &[
            "zone_or_zonelist_name",
            "control_type_schedule_name",
            "control_1_object_type",
            "control_1_name",
        ],
        extensions: &[],
        extension_fields: &[],
    };
    const DUAL: Spec = Spec {
        named: true,
        fields: &[
            "heating_setpoint_temperature_schedule_name",
            "cooling_setpoint_temperature_schedule_name",
        ],
        extensions: &[],
        extension_fields: &[],
    };
    const IDEAL: Spec = Spec {
        named: true,
        fields: &[
            "availability_schedule_name",
            "zone_supply_air_node_name",
            "zone_exhaust_air_node_name",
            "system_inlet_air_node_name",
            "maximum_heating_supply_air_temperature",
            "minimum_cooling_supply_air_temperature",
            "maximum_heating_supply_air_humidity_ratio",
            "minimum_cooling_supply_air_humidity_ratio",
            "heating_limit",
            "maximum_heating_air_flow_rate",
            "maximum_sensible_heating_capacity",
            "cooling_limit",
            "maximum_cooling_air_flow_rate",
            "maximum_total_cooling_capacity",
        ],
        extensions: &[],
        extension_fields: &[],
    };
    const BASEBOARD: Spec = Spec {
        named: true,
        fields: &[
            "availability_schedule_name",
            "heating_design_capacity_method",
            "heating_design_capacity",
            "heating_design_capacity_per_floor_area",
            "fraction_of_autosized_heating_design_capacity",
            "efficiency",
        ],
        extensions: &[],
        extension_fields: &[],
    };
    const EQUIP_LIST: Spec = Spec {
        named: true,
        fields: &["load_distribution_scheme"],
        extensions: &["equipment"],
        extension_fields: &[
            "zone_equipment_object_type",
            "zone_equipment_name",
            "zone_equipment_cooling_sequence",
            "zone_equipment_heating_or_no_load_sequence",
            "zone_equipment_sequential_cooling_fraction_schedule_name",
            "zone_equipment_sequential_heating_fraction_schedule_name",
        ],
    };
    const CONNECTIONS: Spec = Spec {
        named: false,
        fields: &[
            "zone_name",
            "zone_conditioning_equipment_list_name",
            "zone_air_inlet_node_or_nodelist_name",
            "zone_air_exhaust_node_or_nodelist_name",
            "zone_air_node_name",
            "zone_return_air_node_or_nodelist_name",
        ],
        extensions: &[],
        extension_fields: &[],
    };
    const BUILDING: Spec = Spec {
        named: true,
        fields: &["north_axis"],
        extensions: &[],
        extension_fields: &[],
    };
    const LOCATION: Spec = Spec {
        named: true,
        fields: &["latitude", "longitude", "time_zone", "elevation"],
        extensions: &[],
        extension_fields: &[],
    };
    const RUN_PERIOD: Spec = Spec {
        named: true,
        fields: &[
            "begin_month",
            "begin_day_of_month",
            "begin_year",
            "end_month",
            "end_day_of_month",
            "end_year",
        ],
        extensions: &[],
        extension_fields: &[],
    };
    const TIMESTEP: Spec = Spec {
        named: false,
        fields: &["number_of_timesteps_per_hour"],
        extensions: &[],
        extension_fields: &[],
    };
    const GEOMETRY: Spec = Spec {
        named: false,
        fields: &[
            "starting_vertex_position",
            "vertex_entry_direction",
            "coordinate_system",
        ],
        extensions: &[],
        extension_fields: &[],
    };
    const OVERHANG: Spec = Spec {
        named: true,
        fields: &[
            "window_or_door_name",
            "height_above_window_or_door",
            "tilt_angle_from_window_door",
            "left_extension_from_window_door_width",
            "right_extension_from_window_door_width",
            "depth",
        ],
        extensions: &[],
        extension_fields: &[],
    };
    const OUTPUT: Spec = Spec {
        named: false,
        fields: &[
            "key_value",
            "variable_name",
            "reporting_frequency",
            "schedule_name",
        ],
        extensions: &[],
        extension_fields: &[],
    };
    const VERSION: Spec = Spec {
        named: false,
        fields: &["version_identifier"],
        extensions: &[],
        extension_fields: &[],
    };
    const SIMCONTROL: Spec = Spec {
        named: false,
        fields: &[
            "do_zone_sizing_calculation",
            "do_system_sizing_calculation",
            "do_plant_sizing_calculation",
            "run_simulation_for_sizing_periods",
            "run_simulation_for_weather_file_run_periods",
        ],
        extensions: &[],
        extension_fields: &[],
    };
    Some(match object_type {
        "Zone" => &ZONE,
        "BuildingSurface:Detailed" => &SURFACE,
        "FenestrationSurface:Detailed" => &WINDOW,
        "Material" => &MATERIAL,
        "WindowMaterial:SimpleGlazingSystem" => &GLAZING,
        "Construction" => &CONSTRUCTION,
        "Schedule:Constant" => &SCHEDULE_CONSTANT,
        "Schedule:Compact" => &SCHEDULE_COMPACT,
        "People" => &PEOPLE,
        "Lights" => &LIGHTS,
        "ElectricEquipment" => &EQUIPMENT,
        "GasEquipment" => &GAS,
        "OtherEquipment" => &OTHER,
        "ZoneInfiltration:DesignFlowRate" => &INFILTRATION,
        "ZoneControl:Thermostat" => &THERMOSTAT,
        "ThermostatSetpoint:DualSetpoint" => &DUAL,
        "ZoneHVAC:IdealLoadsAirSystem" => &IDEAL,
        "ZoneHVAC:Baseboard:Convective:Electric" => &BASEBOARD,
        "ZoneHVAC:EquipmentList" => &EQUIP_LIST,
        "ZoneHVAC:EquipmentConnections" => &CONNECTIONS,
        "Building" => &BUILDING,
        "Site:Location" => &LOCATION,
        "RunPeriod" => &RUN_PERIOD,
        "Timestep" => &TIMESTEP,
        "GlobalGeometryRules" => &GEOMETRY,
        "Shading:Overhang" => &OVERHANG,
        "Output:Variable" => &OUTPUT,
        "Version" => &VERSION,
        "SimulationControl" => &SIMCONTROL,
        _ => return None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::load_idf;

    #[test]
    fn epjson_surface_matches_idf_geometry() {
        let epjson = r#"{
          "GlobalGeometryRules": {
            "GlobalGeometryRules 1": {
              "starting_vertex_position": "UpperLeftCorner",
              "vertex_entry_direction": "Counterclockwise",
              "coordinate_system": "World"
            }
          },
          "Zone": {
            "Box": { "volume": 30 }
          },
          "Material": {
            "Ins": {
              "roughness": "Rough",
              "thickness": 0.1,
              "conductivity": 0.04,
              "density": 30,
              "specific_heat": 1400
            }
          },
          "Construction": { "Wall": { "outside_layer": "Ins" } },
          "BuildingSurface:Detailed": {
            "SouthWall": {
              "surface_type": "Wall",
              "construction_name": "Wall",
              "zone_name": "Box",
              "outside_boundary_condition": "Outdoors",
              "vertices": [
                {"vertex_x_coordinate": 0, "vertex_y_coordinate": 0, "vertex_z_coordinate": 3},
                {"vertex_x_coordinate": 0, "vertex_y_coordinate": 0, "vertex_z_coordinate": 0},
                {"vertex_x_coordinate": 10, "vertex_y_coordinate": 0, "vertex_z_coordinate": 0},
                {"vertex_x_coordinate": 10, "vertex_y_coordinate": 0, "vertex_z_coordinate": 3}
              ]
            }
          }
        }"#;
        let model = load_idf(epjson).unwrap();
        assert!((model.surfaces[0].area - 30.0).abs() < 1e-6);
        assert!(model.surfaces[0].normal[1] < -0.99);
        assert!((model.surfaces[0].azimuth_deg - 180.0).abs() < 1e-6);
        assert!((model.zones[0].volume - 30.0).abs() < 1e-9);
    }

    #[test]
    fn compact_schedule_data_array_is_positional() {
        let objects = parse_epjson(
            r#"{
              "Schedule:Compact": {
                "Occ": {
                  "schedule_type_limits_name": "Fraction",
                  "data": [
                    {"field": "Through: 12/31"},
                    {"field": "For: AllDays"},
                    {"field": "Until: 24:00"},
                    {"field": 1}
                  ]
                }
              }
            }"#,
        )
        .unwrap();
        assert_eq!(objects[0].field(0), "Occ");
        assert_eq!(objects[0].field(3), "For: AllDays");
        assert_eq!(objects[0].field(5), "1");
    }

    #[test]
    fn example_epjson_matches_idf_results() {
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
        let idf =
            load_idf(&std::fs::read_to_string(root.join("examples/box.idf")).unwrap()).unwrap();
        let epjson =
            load_idf(&std::fs::read_to_string(root.join("examples/box.epJSON")).unwrap()).unwrap();
        assert_eq!(idf.surfaces.len(), epjson.surfaces.len());
        let mut idf_surfaces = idf.surfaces.clone();
        let mut epjson_surfaces = epjson.surfaces.clone();
        idf_surfaces.sort_by(|left, right| left.name.cmp(&right.name));
        epjson_surfaces.sort_by(|left, right| left.name.cmp(&right.name));
        for (left, right) in idf_surfaces.iter().zip(&epjson_surfaces) {
            assert_eq!(left.name, right.name);
            assert!(
                (left.area - right.area).abs() < 1e-6,
                "{} area {} vs {}",
                left.name,
                left.area,
                right.area
            );
            assert_eq!(left.boundary, right.boundary);
        }
        let weather = crate::epw::parse_epw(
            &std::fs::read_to_string(root.join("examples/weather.epw")).unwrap(),
        )
        .unwrap();
        let options = crate::simulate::SimOptions { warmup_days: 0 };
        let from_idf = crate::simulate::simulate(&idf, &weather, &options).unwrap();
        let from_epjson = crate::simulate::simulate(&epjson, &weather, &options).unwrap();
        let scale = from_idf.heating_j.abs().max(1.0);
        assert!(
            (from_idf.heating_j - from_epjson.heating_j).abs() / scale < 1e-8,
            "heating {} vs {}",
            from_idf.heating_j,
            from_epjson.heating_j
        );
        assert!(
            (from_idf.cooling_j - from_epjson.cooling_j).abs() < 1e-6,
            "cooling {} vs {}",
            from_idf.cooling_j,
            from_epjson.cooling_j
        );
    }
}
