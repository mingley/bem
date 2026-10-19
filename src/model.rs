//! Translation of an EnergyPlus IDF into the objects this port can simulate.
//!
//! Field order follows Energy+.idd: zones, schedules, materials, constructions,
//! surfaces, internal gains, infiltration, ideal loads, and convective electric baseboard.

use crate::convection::{azimuth_from_normal, Roughness};
use crate::idf::{parse_f64, parse_idf, IdfObject};
use crate::simple_glazing::simple_glazing;
use crate::DEG_TO_RAD;

#[derive(Clone, Debug)]
pub struct Model {
    pub location: Option<Location>,
    pub north_axis: f64,
    pub world_coordinates: bool,
    pub timesteps_per_hour: usize,
    pub run_period: Option<RunPeriod>,
    pub schedules: Vec<Schedule>,
    pub zones: Vec<Zone>,
    pub surfaces: Vec<Surface>,
    pub outputs: Vec<OutputVariable>,
    pub warnings: Vec<String>,
}

#[derive(Clone, Debug)]
pub struct Location {
    pub name: String,
    pub latitude: f64,
    pub longitude: f64,
    pub time_zone: f64,
    pub elevation: f64,
}

#[derive(Clone, Debug)]
pub struct RunPeriod {
    pub begin_month: u32,
    pub begin_day: u32,
    pub begin_year: Option<i32>,
    pub end_month: u32,
    pub end_day: u32,
    pub end_year: Option<i32>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReportFrequency {
    Timestep,
    Hourly,
    Daily,
    Monthly,
    RunPeriod,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OutputVariable {
    pub key: String,
    pub variable: String,
    pub frequency: ReportFrequency,
}

#[derive(Clone, Debug)]
pub enum Schedule {
    Constant(f64),
    Compact(CompactSchedule),
}

#[derive(Clone, Debug)]
pub struct CompactSchedule {
    pub through: Vec<Through>,
}

#[derive(Clone, Debug)]
pub struct Through {
    pub month: u32,
    pub day: u32,
    pub rules: Vec<DayRule>,
}

#[derive(Clone, Debug)]
pub struct DayRule {
    pub kind: DayKind,
    pub untils: Vec<(u32, f64)>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DayKind {
    All,
    Weekday,
    Weekend,
    Dow(u32),
    Unused,
}

#[derive(Clone, Debug)]
pub struct Zone {
    pub name: String,
    pub volume: f64,
    pub floor_area: f64,
    pub heating_schedule: Option<usize>,
    pub cooling_schedule: Option<usize>,
    pub heating_capacity_w: f64,
    pub cooling_capacity_w: f64,
    /// Zone equipment in `ZoneHVAC:EquipmentList` order. Empty means this zone
    /// still uses the single implicit ideal-loads unit and the capacity fields.
    pub equipment: Vec<ZoneEquipment>,
    pub gains: Vec<Gain>,
    pub infiltration_m3_s: Vec<Infiltration>,
}

/// One piece of zone equipment. Heating and cooling sequence numbers follow
/// `ZoneHVAC:EquipmentList`: a sequence of 0 keeps that mode off, and a lower
/// positive sequence runs first under SequentialLoad.
#[derive(Clone, Debug)]
pub enum ZoneEquipment {
    IdealLoads {
        heating_capacity_w: f64,
        cooling_capacity_w: f64,
        heating_sequence: i32,
        cooling_sequence: i32,
    },
    ElectricBaseboard {
        name: String,
        capacity_w: f64,
        efficiency: f64,
        availability: Option<usize>,
        heating_sequence: i32,
    },
}

#[derive(Clone, Debug)]
pub enum Gain {
    Power {
        schedule: usize,
        convective_w: f64,
        radiant_w: f64,
        latent_w: f64,
    },
    People {
        occupancy: usize,
        activity: Option<usize>,
        count: f64,
        fraction_radiant: f64,
        sensible_fraction: f64,
    },
}

#[derive(Clone, Debug)]
pub struct Infiltration {
    pub schedule: Option<usize>,
    pub flow_m3_s: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Boundary {
    Outdoors,
    Ground,
    Adiabatic,
    /// Far side is the air of another zone. Partition should be entered once.
    Zone(usize),
}

#[derive(Clone, Debug)]
pub struct Layer {
    pub k: f64,
    pub dx: f64,
    pub n: usize,
    pub rho: f64,
    pub cp: f64,
}

#[derive(Clone, Debug)]
pub enum SurfaceKind {
    Opaque {
        layers: Vec<Layer>,
    },
    Window {
        u_factor: f64,
        transmittance: f64,
        inward_absorptance: f64,
    },
}

#[derive(Clone, Debug)]
pub struct Surface {
    pub name: String,
    pub zone: usize,
    pub area: f64,
    /// Outward unit normal, EnergyPlus world axes: x east, y north, z up.
    pub normal: [f64; 3],
    /// Outward azimuth, degrees clockwise from north.
    pub azimuth_deg: f64,
    pub roughness: Roughness,
    pub kind: SurfaceKind,
    pub boundary: Boundary,
    pub solar_absorptance: f64,
    /// Inside-face solar absorptance, from the innermost material layer.
    pub solar_absorptance_int: f64,
    pub thermal_absorptance_ext: f64,
    pub thermal_absorptance_int: f64,
    pub active: bool,
}

impl Schedule {
    pub fn value(&self, year: i32, month: u32, day: u32, minutes: u32) -> f64 {
        match self {
            Schedule::Constant(value) => *value,
            Schedule::Compact(compact) => compact.value(year, month, day, minutes),
        }
    }
}

impl CompactSchedule {
    fn value(&self, year: i32, month: u32, day: u32, minutes: u32) -> f64 {
        let through = self
            .through
            .iter()
            .find(|block| block.month > month || (block.month == month && block.day >= day))
            .or_else(|| self.through.last());
        let Some(through) = through else {
            return 0.0;
        };
        let dow = crate::solar::weekday(year, month, day);
        let weekend = dow == 0 || dow == 6;
        let mut best: Option<(u8, &DayRule)> = None;
        for rule in &through.rules {
            let rank = match rule.kind {
                DayKind::Dow(day) if day == dow => 3,
                DayKind::Weekday if !weekend => 2,
                DayKind::Weekend if weekend => 2,
                DayKind::All => 1,
                _ => 0,
            };
            if rank == 0 {
                continue;
            }
            if best.map(|(have, _)| rank > have).unwrap_or(true) {
                best = Some((rank, rule));
            }
        }
        let Some((_, rule)) = best else {
            return 0.0;
        };
        for (until, value) in &rule.untils {
            if minutes <= *until {
                return *value;
            }
        }
        rule.untils.last().map(|(_, value)| *value).unwrap_or(0.0)
    }
}

pub fn load_idf(text: &str) -> Result<Model, String> {
    load_input(text)
}

/// Load an EnergyPlus IDF or epJSON file. epJSON is detected by a leading `{`.
pub fn load_input(text: &str) -> Result<Model, String> {
    let trimmed = text.trim_start();
    if trimmed.starts_with('{') {
        interpret(&crate::epjson::parse_epjson(text)?)
    } else {
        interpret(&parse_idf(text))
    }
}

fn interpret(objects: &[IdfObject]) -> Result<Model, String> {
    let mut warnings = Vec::new();
    let mut outputs = Vec::new();
    let mut north_axis = 0.0;
    let mut world_coordinates = false;
    let mut timesteps_per_hour = 6usize;
    let mut location = None;
    let mut run_period = None;
    let mut schedules: Vec<(String, Schedule)> = Vec::new();
    let mut materials: Vec<Material> = Vec::new();
    let mut glazing: Vec<SimpleGlazing> = Vec::new();
    let mut constructions: Vec<(String, Vec<String>)> = Vec::new();
    let mut zones_raw: Vec<ZoneDraft> = Vec::new();
    let mut surface_drafts: Vec<SurfaceDraft> = Vec::new();
    let mut window_drafts: Vec<WindowDraft> = Vec::new();
    let mut gains_raw: Vec<GainDraft> = Vec::new();
    let mut inf_raw: Vec<InfDraft> = Vec::new();
    let mut thermostats: Vec<ThermostatDraft> = Vec::new();
    let mut dual_setpoints: Vec<DualDraft> = Vec::new();
    let mut ideal: Vec<IdealDraft> = Vec::new();
    let mut baseboards: Vec<BaseboardDraft> = Vec::new();
    let mut equipment_lists: Vec<EquipList> = Vec::new();
    let mut connections: Vec<EquipConn> = Vec::new();

    for object in objects {
        if object.is_type("Building") {
            north_axis = parse_f64(object.field(1)).unwrap_or(0.0);
        } else if object.is_type("GlobalGeometryRules") {
            world_coordinates = object.field(2).eq_ignore_ascii_case("World");
        } else if object.is_type("Timestep") {
            let n = parse_f64(object.field(0)).unwrap_or(6.0) as usize;
            if n == 0 {
                return Err("Timestep must be at least 1".to_string());
            }
            timesteps_per_hour = n;
        } else if object.is_type("Site:Location") {
            location = Some(Location {
                name: object.field(0).to_string(),
                latitude: required_f64(object, 1, "Site:Location latitude")?,
                longitude: required_f64(object, 2, "Site:Location longitude")?,
                time_zone: required_f64(object, 3, "Site:Location time zone")?,
                elevation: parse_f64(object.field(4)).unwrap_or(0.0),
            });
        } else if object.is_type("RunPeriod") {
            run_period = Some(RunPeriod {
                begin_month: required_f64(object, 1, "RunPeriod begin month")? as u32,
                begin_day: required_f64(object, 2, "RunPeriod begin day")? as u32,
                begin_year: parse_f64(object.field(3)).map(|y| y as i32),
                end_month: required_f64(object, 4, "RunPeriod end month")? as u32,
                end_day: required_f64(object, 5, "RunPeriod end day")? as u32,
                end_year: parse_f64(object.field(6)).map(|y| y as i32),
            });
        } else if object.is_type("Schedule:Constant") {
            schedules.push((
                object.field(0).to_string(),
                Schedule::Constant(parse_f64(object.field(2)).unwrap_or(0.0)),
            ));
        } else if object.is_type("Schedule:Compact") {
            schedules.push((object.field(0).to_string(), parse_compact(object)?));
        } else if object.is_type("Material") {
            materials.push(parse_material(object)?);
        } else if object.is_type("WindowMaterial:SimpleGlazingSystem") {
            glazing.push(SimpleGlazing {
                name: object.field(0).to_string(),
                u_factor: required_f64(object, 1, "simple glazing U-factor")?,
                shgc: required_f64(object, 2, "simple glazing SHGC")?,
            });
        } else if object.is_type("Construction") {
            let layers = object
                .fields
                .iter()
                .skip(1)
                .cloned()
                .filter(|f| !f.is_empty())
                .collect();
            constructions.push((object.field(0).to_string(), layers));
        } else if object.is_type("Zone") {
            zones_raw.push(ZoneDraft {
                name: object.field(0).to_string(),
                rel_north: parse_f64(object.field(1)).unwrap_or(0.0),
                origin: [
                    parse_f64(object.field(2)).unwrap_or(0.0),
                    parse_f64(object.field(3)).unwrap_or(0.0),
                    parse_f64(object.field(4)).unwrap_or(0.0),
                ],
                multiplier: parse_f64(object.field(6)).unwrap_or(1.0).max(1.0),
                ceiling_height: parse_f64(object.field(7)),
                volume: parse_f64(object.field(8)),
                floor_area: parse_f64(object.field(9)),
            });
        } else if object.is_type("BuildingSurface:Detailed") {
            surface_drafts.push(parse_surface(object)?);
        } else if object.is_type("FenestrationSurface:Detailed") {
            window_drafts.push(parse_window(object)?);
        } else if object.is_type("People") {
            gains_raw.push(parse_people(object)?);
        } else if object.is_type("Lights") {
            gains_raw.push(parse_lights(object)?);
        } else if object.is_type("ElectricEquipment") || object.is_type("GasEquipment") {
            gains_raw.push(parse_equipment(object, 1)?);
        } else if object.is_type("OtherEquipment") {
            // Fuel Type sits in front of the zone name in the current IDD.
            gains_raw.push(parse_equipment(object, 2)?);
        } else if object.is_type("ZoneInfiltration:DesignFlowRate") {
            inf_raw.push(InfDraft {
                zone: object.field(1).to_string(),
                schedule: optional_name(object.field(2)),
                method: object.field(3).to_string(),
                flow: parse_f64(object.field(4)),
                per_area: parse_f64(object.field(5)),
                per_exterior: parse_f64(object.field(6)),
                ach: parse_f64(object.field(7)),
            });
        } else if object.is_type("ThermostatSetpoint:DualSetpoint") {
            dual_setpoints.push(DualDraft {
                name: object.field(0).to_string(),
                heating: object.field(1).to_string(),
                cooling: object.field(2).to_string(),
            });
        } else if object.is_type("ZoneControl:Thermostat") {
            thermostats.push(ThermostatDraft {
                zone: object.field(1).to_string(),
                control_type: object.field(3).to_string(),
                control_name: object.field(4).to_string(),
            });
        } else if object.is_type("ZoneHVAC:IdealLoadsAirSystem") {
            ideal.push(parse_ideal(object));
        } else if object.is_type("ZoneHVAC:Baseboard:Convective:Electric") {
            baseboards.push(parse_baseboard(object)?);
        } else if object.is_type("ZoneHVAC:EquipmentList") {
            equipment_lists.push(parse_equip_list(object)?);
        } else if object.is_type("ZoneHVAC:EquipmentConnections") {
            connections.push(EquipConn {
                zone: object.field(0).to_string(),
                list: object.field(1).to_string(),
            });
        } else if object.is_type("Output:Variable") {
            match parse_output_variable(object) {
                Ok(request) => outputs.push(request),
                Err(message) => warnings.push(message),
            }
        } else if object.is_type("Version")
            || object.is_type("SimulationControl")
            || object.is_type("ScheduleTypeLimits")
            || object.is_type("Output:Meter")
            || object.is_type("OutputControl:Table:Style")
            || object.is_type("Output:Table:SummaryReports")
            || object.is_type("Output:SQLite")
            || object.is_type("Output:Constructions")
            || object.is_type("Output:Surfaces:List")
            || object.is_type("Output:Diagnostics")
        {
            // Bookkeeping that does not change the heat balance.
        } else {
            warnings.push(format!(
                "ignored {} '{}'",
                object.object_type,
                object.field(0)
            ));
        }
    }

    if zones_raw.is_empty() {
        return Err("IDF has no Zone objects".to_string());
    }

    let schedule_index = |name: &str| -> Result<usize, String> {
        schedules
            .iter()
            .position(|(n, _)| n.eq_ignore_ascii_case(name))
            .ok_or_else(|| format!("schedule '{name}' was not found"))
    };

    let mut zones = Vec::new();
    for draft in &zones_raw {
        zones.push(Zone {
            name: draft.name.clone(),
            volume: 0.0,
            floor_area: 0.0,
            heating_schedule: None,
            cooling_schedule: None,
            heating_capacity_w: f64::INFINITY,
            cooling_capacity_w: f64::INFINITY,
            equipment: Vec::new(),
            gains: Vec::new(),
            infiltration_m3_s: Vec::new(),
        });
    }
    let zone_index = |name: &str| -> Result<usize, String> {
        zones_raw
            .iter()
            .position(|zone| zone.name.eq_ignore_ascii_case(name))
            .ok_or_else(|| format!("zone '{name}' was not found"))
    };

    let mut surfaces = Vec::new();
    for draft in &surface_drafts {
        let zone = zone_index(&draft.zone)?;
        let (kind, solar_abs, solar_int, thermal_ext, thermal_int, roughness) =
            construction_kind(&draft.construction, &constructions, &materials, &glazing)?;
        let world = to_world(
            &draft.vertices,
            &zones_raw[zone],
            north_axis,
            world_coordinates,
        );
        let (area, normal) = polygon_area_normal(&world);
        surfaces.push(Surface {
            name: draft.name.clone(),
            zone,
            area,
            normal,
            azimuth_deg: azimuth_from_normal(normal),
            roughness,
            kind,
            boundary: Boundary::Outdoors,
            solar_absorptance: solar_abs,
            solar_absorptance_int: solar_int,
            thermal_absorptance_ext: thermal_ext,
            thermal_absorptance_int: thermal_int,
            active: area > 1.0e-8,
        });
        // Boundary is filled after every surface exists, so pairs can be matched.
        let _ = draft;
    }
    // Second pass for boundaries, now that names exist.
    for (i, draft) in surface_drafts.iter().enumerate() {
        surfaces[i].boundary = boundary_of(
            &draft.boundary,
            &draft.boundary_object,
            &surface_drafts,
            &zones_raw,
        )?;
    }
    for draft in &window_drafts {
        let parent = surface_drafts
            .iter()
            .position(|surf| surf.name.eq_ignore_ascii_case(&draft.parent))
            .ok_or_else(|| {
                format!(
                    "window '{}' parent surface '{}' was not found",
                    draft.name, draft.parent
                )
            })?;
        let zone = zone_index(&surface_drafts[parent].zone)?;
        let (kind, solar_abs, solar_int, thermal_ext, thermal_int, roughness) =
            construction_kind(&draft.construction, &constructions, &materials, &glazing)?;
        let SurfaceKind::Window { .. } = kind else {
            return Err(format!(
                "fenestration '{}' construction '{}' is not a simple glazing system",
                draft.name, draft.construction
            ));
        };
        let world = to_world(
            &draft.vertices,
            &zones_raw[zone],
            north_axis,
            world_coordinates,
        );
        let (mut area, normal) = polygon_area_normal(&world);
        area *= draft.multiplier.max(1.0);
        let parent_normal = surfaces[parent].normal;
        let normal = if dot(normal, parent_normal) < 0.0 {
            scale(normal, -1.0)
        } else {
            normal
        };
        surfaces.push(Surface {
            name: draft.name.clone(),
            zone,
            area,
            normal,
            azimuth_deg: azimuth_from_normal(normal),
            roughness,
            kind,
            boundary: surfaces[parent].boundary,
            solar_absorptance: solar_abs,
            solar_absorptance_int: solar_int,
            thermal_absorptance_ext: thermal_ext,
            thermal_absorptance_int: thermal_int,
            active: area > 1.0e-8,
        });
    }

    // Vertex order follows the EnergyPlus right-hand rule (outward normal).

    // A paired interzone surface is solved from one side so the wall is not counted twice.
    let names: Vec<String> = surfaces.iter().map(|s| s.name.clone()).collect();
    for i in 0..surfaces.len() {
        if let Boundary::Zone(_) = surfaces[i].boundary {
            if let Some(partner) = surface_drafts.iter().enumerate().find_map(|(_idx, draft)| {
                if draft.name.eq_ignore_ascii_case(&names[i]) {
                    let other = draft.boundary_object.as_str();
                    names
                        .iter()
                        .position(|name| name.eq_ignore_ascii_case(other))
                } else {
                    None
                }
            }) {
                if partner < i {
                    surfaces[i].active = false;
                }
            }
        }
    }

    for (zi, draft) in zones_raw.iter().enumerate() {
        let floor_from_geom: f64 = surfaces
            .iter()
            .filter(|s| s.zone == zi && s.active && s.normal[2] < -0.7)
            .map(|s| s.area)
            .sum();
        let floor_area = draft.floor_area.unwrap_or(floor_from_geom);
        let volume = match draft.volume {
            Some(v) => v,
            None => {
                let height = draft.ceiling_height.unwrap_or(0.0);
                if floor_area > 0.0 && height > 0.0 {
                    floor_area * height
                } else {
                    return Err(format!(
                        "zone '{}' needs a Volume, or Floor Area and Ceiling Height",
                        draft.name
                    ));
                }
            }
        };
        let mult = draft.multiplier;
        zones[zi].volume = volume * mult;
        zones[zi].floor_area = floor_area * mult;
        for surface in surfaces.iter_mut().filter(|s| s.zone == zi) {
            surface.area *= mult;
        }
    }

    let exterior_area = |zone: usize| -> f64 {
        surfaces
            .iter()
            .filter(|s| s.zone == zone && s.active && matches!(s.boundary, Boundary::Outdoors))
            .map(|s| s.area)
            .sum()
    };

    for gain in gains_raw {
        let zi = zone_index(&gain.zone)?;
        let floor = zones[zi].floor_area / zones_raw[zi].multiplier;
        let mult = zones_raw[zi].multiplier;
        match gain.spec {
            GainSpec::People {
                occupancy,
                activity,
                count,
                per_area,
                area_per,
                method,
                fraction_radiant,
                sensible_fraction,
            } => {
                let n = people_count(&method, count, per_area, area_per, floor)?;
                zones[zi].gains.push(Gain::People {
                    occupancy: schedule_index(&occupancy)?,
                    activity: match activity {
                        Some(name) => Some(schedule_index(&name)?),
                        None => None,
                    },
                    count: n * mult,
                    fraction_radiant,
                    sensible_fraction,
                });
            }
            GainSpec::Power {
                schedule,
                convective_w,
                radiant_w,
                latent_w,
            } => {
                let schedule = schedule_index(&schedule)?;
                zones[zi].gains.push(Gain::Power {
                    schedule,
                    convective_w: convective_w * mult,
                    radiant_w: radiant_w * mult,
                    latent_w: latent_w * mult,
                });
            }
        }
    }

    for inf in inf_raw {
        let zi = zone_index(&inf.zone)?;
        let floor = zones[zi].floor_area / zones_raw[zi].multiplier.max(1.0);
        let volume = zones[zi].volume / zones_raw[zi].multiplier.max(1.0);
        let exterior = exterior_area(zi) / zones_raw[zi].multiplier.max(1.0);
        let flow = infiltration_flow(&inf, floor, volume, exterior)?;
        let schedule = match inf.schedule {
            Some(name) => Some(schedule_index(&name)?),
            None => None,
        };
        zones[zi].infiltration_m3_s.push(Infiltration {
            schedule,
            flow_m3_s: flow * zones_raw[zi].multiplier,
        });
    }

    for conn in &connections {
        let zi = zone_index(&conn.zone)?;
        let list = equipment_lists
            .iter()
            .find(|list| list.name.eq_ignore_ascii_case(&conn.list))
            .ok_or_else(|| format!("equipment list '{}' was not found", conn.list))?;
        zones[zi].heating_capacity_w = 0.0;
        zones[zi].cooling_capacity_w = 0.0;
        zones[zi].equipment.clear();
        for item in &list.items {
            if item
                .kind
                .eq_ignore_ascii_case("ZoneHVAC:IdealLoadsAirSystem")
            {
                let unit = ideal
                    .iter()
                    .find(|unit| unit.name.eq_ignore_ascii_case(&item.name))
                    .ok_or_else(|| format!("ideal loads system '{}' was not found", item.name))?;
                zones[zi].heating_capacity_w = unit.heating_capacity_w;
                zones[zi].cooling_capacity_w = unit.cooling_capacity_w;
                zones[zi].equipment.push(ZoneEquipment::IdealLoads {
                    heating_capacity_w: unit.heating_capacity_w,
                    cooling_capacity_w: unit.cooling_capacity_w,
                    heating_sequence: item.heating_sequence,
                    cooling_sequence: item.cooling_sequence,
                });
            } else if item
                .kind
                .eq_ignore_ascii_case("ZoneHVAC:Baseboard:Convective:Electric")
            {
                let unit = baseboards
                    .iter()
                    .find(|unit| unit.name.eq_ignore_ascii_case(&item.name))
                    .ok_or_else(|| format!("electric baseboard '{}' was not found", item.name))?;
                let capacity_w = match unit.capacity {
                    BaseboardCapacity::Watts(watts) => watts,
                    BaseboardCapacity::PerFloorArea(per_area) => per_area * zones[zi].floor_area,
                };
                let availability = match &unit.availability {
                    Some(name) => Some(schedule_index(name)?),
                    None => None,
                };
                zones[zi].equipment.push(ZoneEquipment::ElectricBaseboard {
                    name: unit.name.clone(),
                    capacity_w,
                    efficiency: unit.efficiency,
                    availability,
                    heating_sequence: item.heating_sequence,
                });
            } else {
                return Err(format!(
                    "{} '{}' is not ported yet. This slice simulates ZoneHVAC:IdealLoadsAirSystem and ZoneHVAC:Baseboard:Convective:Electric",
                    item.kind, item.name
                ));
            }
        }
    }

    for stat in &thermostats {
        let zi = zone_index(&stat.zone)?;
        if !stat
            .control_type
            .eq_ignore_ascii_case("ThermostatSetpoint:DualSetpoint")
        {
            return Err(format!(
                "thermostat control '{}' is not ported yet; use ThermostatSetpoint:DualSetpoint",
                stat.control_type
            ));
        }
        let dual = dual_setpoints
            .iter()
            .find(|dual| dual.name.eq_ignore_ascii_case(&stat.control_name))
            .ok_or_else(|| format!("dual setpoint '{}' was not found", stat.control_name))?;
        zones[zi].heating_schedule = Some(schedule_index(&dual.heating)?);
        zones[zi].cooling_schedule = Some(schedule_index(&dual.cooling)?);
    }

    let schedules = schedules
        .into_iter()
        .map(|(_, schedule)| schedule)
        .collect();
    Ok(Model {
        location,
        north_axis,
        world_coordinates,
        timesteps_per_hour,
        run_period,
        schedules,
        zones,
        surfaces,
        outputs,
        warnings,
    })
}

fn boundary_of(
    boundary: &str,
    object: &str,
    surfaces: &[SurfaceDraft],
    zones: &[ZoneDraft],
) -> Result<Boundary, String> {
    if boundary.eq_ignore_ascii_case("Outdoors") || boundary.is_empty() {
        Ok(Boundary::Outdoors)
    } else if boundary.eq_ignore_ascii_case("Ground")
        || boundary.eq_ignore_ascii_case("GroundFCfactorMethod")
        || boundary.eq_ignore_ascii_case("GroundSlabPreprocessorAverage")
    {
        Ok(Boundary::Ground)
    } else if boundary.eq_ignore_ascii_case("Adiabatic") {
        Ok(Boundary::Adiabatic)
    } else if boundary.eq_ignore_ascii_case("Surface") {
        let partner = surfaces
            .iter()
            .find(|surf| surf.name.eq_ignore_ascii_case(object))
            .ok_or_else(|| format!("outside boundary surface '{object}' was not found"))?;
        let zone = zones
            .iter()
            .position(|zone| zone.name.eq_ignore_ascii_case(&partner.zone))
            .ok_or_else(|| format!("partner surface '{}' has no zone", partner.name))?;
        Ok(Boundary::Zone(zone))
    } else if boundary.eq_ignore_ascii_case("Zone") {
        let zone = zones
            .iter()
            .position(|zone| zone.name.eq_ignore_ascii_case(object))
            .ok_or_else(|| format!("outside boundary zone '{object}' was not found"))?;
        Ok(Boundary::Zone(zone))
    } else {
        Err(format!(
            "outside boundary condition '{boundary}' is not ported yet"
        ))
    }
}

fn people_count(
    method: &str,
    count: f64,
    per_area: f64,
    area_per: f64,
    floor: f64,
) -> Result<f64, String> {
    if method.eq_ignore_ascii_case("People/Area") {
        Ok(per_area * floor)
    } else if method.eq_ignore_ascii_case("Area/Person") {
        if area_per <= 0.0 {
            return Err("Area/Person must be positive".to_string());
        }
        Ok(floor / area_per)
    } else {
        Ok(count)
    }
}

fn infiltration_flow(
    inf: &InfDraft,
    floor: f64,
    volume: f64,
    exterior: f64,
) -> Result<f64, String> {
    let method = inf.method.as_str();
    if method.eq_ignore_ascii_case("Flow/Zone") {
        Ok(inf.flow.unwrap_or(0.0))
    } else if method.eq_ignore_ascii_case("Flow/Area") {
        Ok(inf.per_area.unwrap_or(0.0) * floor)
    } else if method.eq_ignore_ascii_case("Flow/ExteriorArea")
        || method.eq_ignore_ascii_case("Flow/ExteriorWallArea")
    {
        Ok(inf.per_exterior.unwrap_or(0.0) * exterior)
    } else if method.eq_ignore_ascii_case("AirChanges/Hour")
        || method.is_empty() && inf.ach.is_some()
    {
        Ok(inf.ach.unwrap_or(0.0) * volume / 3600.0)
    } else if method.is_empty() {
        if let Some(flow) = inf.flow {
            Ok(flow)
        } else if let Some(per_area) = inf.per_area {
            Ok(per_area * floor)
        } else {
            Ok(0.0)
        }
    } else {
        Err(format!("infiltration method '{method}' is not ported yet"))
    }
}

fn construction_kind(
    name: &str,
    constructions: &[(String, Vec<String>)],
    materials: &[Material],
    glazing: &[SimpleGlazing],
) -> Result<(SurfaceKind, f64, f64, f64, f64, Roughness), String> {
    let layers = constructions
        .iter()
        .find(|(n, _)| n.eq_ignore_ascii_case(name))
        .ok_or_else(|| format!("construction '{name}' was not found"))?;
    if layers.1.is_empty() {
        return Err(format!("construction '{name}' has no layers"));
    }
    if let Some(glass) = glazing
        .iter()
        .find(|g| g.name.eq_ignore_ascii_case(&layers.1[0]))
    {
        let u_factor = glass.u_factor;
        let optics = simple_glazing(u_factor, glass.shgc);
        return Ok((
            SurfaceKind::Window {
                u_factor,
                transmittance: optics.transmittance,
                inward_absorptance: optics.inward_absorptance(),
            },
            0.0,
            0.84,
            0.84,
            0.0,
            Roughness::VerySmooth,
        ));
    }
    let mut built = Vec::new();
    let mut solar = 0.6;
    let mut solar_int = 0.6;
    let mut thermal_ext = 0.9;
    let mut thermal_int = 0.9;
    let mut roughness = Roughness::MediumRough;
    for (i, layer_name) in layers.1.iter().enumerate() {
        let mat = materials
            .iter()
            .find(|mat| mat.name.eq_ignore_ascii_case(layer_name))
            .ok_or_else(|| {
                format!("material '{layer_name}' in construction '{name}' was not found")
            })?;
        if mat.thickness <= 0.0 || mat.conductivity <= 0.0 {
            return Err(format!(
                "material '{}' needs thickness and conductivity",
                mat.name
            ));
        }
        let n = ((mat.thickness / 0.05).round() as usize).clamp(1, 20);
        let dx = mat.thickness / n as f64;
        if i == 0 {
            solar = mat.solar_absorptance;
            thermal_ext = mat.thermal_absorptance;
            roughness = mat.roughness;
        }
        thermal_int = mat.thermal_absorptance;
        solar_int = mat.solar_absorptance;
        built.push(Layer {
            k: mat.conductivity,
            dx,
            n,
            rho: mat.density.max(1.0),
            cp: mat.specific_heat.max(1.0),
        });
    }
    Ok((
        SurfaceKind::Opaque { layers: built },
        solar,
        solar_int,
        thermal_ext,
        thermal_int,
        roughness,
    ))
}

fn to_world(
    vertices: &[[f64; 3]],
    zone: &ZoneDraft,
    north_axis: f64,
    world: bool,
) -> Vec<[f64; 3]> {
    if world {
        return vertices.to_vec();
    }
    let cos_z = (-zone.rel_north * DEG_TO_RAD).cos();
    let sin_z = (-zone.rel_north * DEG_TO_RAD).sin();
    let cos_b = (-north_axis * DEG_TO_RAD).cos();
    let sin_b = (-north_axis * DEG_TO_RAD).sin();
    vertices
        .iter()
        .map(|v| {
            let xb = v[0] * cos_z - v[1] * sin_z + zone.origin[0];
            let yb = v[0] * sin_z + v[1] * cos_z + zone.origin[1];
            [
                xb * cos_b - yb * sin_b,
                xb * sin_b + yb * cos_b,
                v[2] + zone.origin[2],
            ]
        })
        .collect()
}

fn polygon_area_normal(vertices: &[[f64; 3]]) -> (f64, [f64; 3]) {
    let mut n = [0.0, 0.0, 0.0];
    if vertices.len() < 3 {
        return (0.0, [0.0, 0.0, 1.0]);
    }
    for i in 0..vertices.len() {
        let a = vertices[i];
        let b = vertices[(i + 1) % vertices.len()];
        n[0] += (a[1] - b[1]) * (a[2] + b[2]);
        n[1] += (a[2] - b[2]) * (a[0] + b[0]);
        n[2] += (a[0] - b[0]) * (a[1] + b[1]);
    }
    let mag = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    if mag < 1.0e-12 {
        return (0.0, [0.0, 0.0, 1.0]);
    }
    (mag * 0.5, [n[0] / mag, n[1] / mag, n[2] / mag])
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn parse_compact(object: &IdfObject) -> Result<Schedule, String> {
    let mut through = Vec::new();
    let mut i = 2;
    while i < object.fields.len() {
        let field = object.field(i);
        if !field.to_ascii_lowercase().starts_with("through:") {
            return Err(format!(
                "Schedule:Compact '{}': expected Through: at '{}'",
                object.field(0),
                field
            ));
        }
        let date = field
            .split_once(':')
            .map(|(_, rest)| rest.trim())
            .unwrap_or("");
        let (month, day) = parse_month_day(date)
            .ok_or_else(|| format!("Schedule:Compact '{}': bad date '{date}'", object.field(0)))?;
        i += 1;
        let mut rules = Vec::new();
        while i < object.fields.len()
            && !object.field(i).to_ascii_lowercase().starts_with("through:")
        {
            let for_field = object.field(i);
            if !for_field.to_ascii_lowercase().starts_with("for:") {
                return Err(format!(
                    "Schedule:Compact '{}': expected For: at '{}'",
                    object.field(0),
                    for_field
                ));
            }
            let kind = day_kind(
                for_field
                    .split_once(':')
                    .map(|(_, r)| r.trim())
                    .unwrap_or(""),
            );
            i += 1;
            let mut untils = Vec::new();
            while i < object.fields.len()
                && object.field(i).to_ascii_lowercase().starts_with("until:")
            {
                let time = object
                    .field(i)
                    .split_once(':')
                    .map(|(_, r)| r.trim())
                    .unwrap_or("");
                let minutes = parse_until(time).ok_or_else(|| {
                    format!("Schedule:Compact '{}': bad time '{time}'", object.field(0))
                })?;
                i += 1;
                let value = parse_f64(object.field(i)).ok_or_else(|| {
                    format!(
                        "Schedule:Compact '{}': value missing after {time}",
                        object.field(0)
                    )
                })?;
                untils.push((minutes, value));
                i += 1;
            }
            rules.push(DayRule { kind, untils });
        }
        through.push(Through { month, day, rules });
    }
    Ok(Schedule::Compact(CompactSchedule { through }))
}

fn parse_month_day(text: &str) -> Option<(u32, u32)> {
    let (m, d) = text.split_once('/')?;
    Some((m.trim().parse().ok()?, d.trim().parse().ok()?))
}

fn parse_until(text: &str) -> Option<u32> {
    let (h, m) = text.split_once(':')?;
    let hour: u32 = h.trim().parse().ok()?;
    let minute: u32 = m.trim().parse().ok()?;
    Some(hour * 60 + minute)
}

fn day_kind(text: &str) -> DayKind {
    match text.to_ascii_lowercase().as_str() {
        "alldays" | "allotherdays" => DayKind::All,
        "weekdays" => DayKind::Weekday,
        "weekends" => DayKind::Weekend,
        "sunday" => DayKind::Dow(0),
        "monday" => DayKind::Dow(1),
        "tuesday" => DayKind::Dow(2),
        "wednesday" => DayKind::Dow(3),
        "thursday" => DayKind::Dow(4),
        "friday" => DayKind::Dow(5),
        "saturday" => DayKind::Dow(6),
        _ => DayKind::Unused,
    }
}

fn parse_material(object: &IdfObject) -> Result<Material, String> {
    Ok(Material {
        name: object.field(0).to_string(),
        thickness: required_f64(object, 2, "material thickness")?,
        conductivity: required_f64(object, 3, "material conductivity")?,
        density: required_f64(object, 4, "material density")?,
        specific_heat: required_f64(object, 5, "material specific heat")?,
        thermal_absorptance: parse_f64(object.field(6)).unwrap_or(0.9),
        solar_absorptance: parse_f64(object.field(7)).unwrap_or(0.6),
        roughness: Roughness::parse(object.field(1)),
    })
}

fn parse_surface(object: &IdfObject) -> Result<SurfaceDraft, String> {
    // Current IDD inserts Space Name after Zone Name, so vertices start at field 11.
    let vertices = vertices_from(object, 11)?;
    Ok(SurfaceDraft {
        name: object.field(0).to_string(),
        construction: object.field(2).to_string(),
        zone: object.field(3).to_string(),
        boundary: object.field(5).to_string(),
        boundary_object: object.field(6).to_string(),
        vertices,
    })
}

fn parse_window(object: &IdfObject) -> Result<WindowDraft, String> {
    Ok(WindowDraft {
        name: object.field(0).to_string(),
        construction: object.field(2).to_string(),
        parent: object.field(3).to_string(),
        multiplier: parse_f64(object.field(7)).unwrap_or(1.0),
        vertices: vertices_from(object, 9)?,
    })
}

fn vertices_from(object: &IdfObject, start: usize) -> Result<Vec<[f64; 3]>, String> {
    let rest = &object.fields[start.min(object.fields.len())..];
    if rest.len() % 3 != 0 || rest.len() < 9 {
        return Err(format!(
            "surface '{}' vertices are incomplete ({} numbers after field {start})",
            object.field(0),
            rest.len()
        ));
    }
    let mut vertices = Vec::new();
    for chunk in rest.chunks(3) {
        let x = chunk[0]
            .parse::<f64>()
            .map_err(|_| format!("bad X on '{}'", object.field(0)))?;
        let y = chunk[1]
            .parse::<f64>()
            .map_err(|_| format!("bad Y on '{}'", object.field(0)))?;
        let z = chunk[2]
            .parse::<f64>()
            .map_err(|_| format!("bad Z on '{}'", object.field(0)))?;
        vertices.push([x, y, z]);
    }
    Ok(vertices)
}

fn parse_people(object: &IdfObject) -> Result<GainDraft, String> {
    Ok(GainDraft {
        zone: object.field(1).to_string(),
        spec: GainSpec::People {
            occupancy: object.field(2).to_string(),
            method: object.field(3).to_string(),
            count: parse_f64(object.field(4)).unwrap_or(0.0),
            per_area: parse_f64(object.field(5)).unwrap_or(0.0),
            area_per: parse_f64(object.field(6)).unwrap_or(0.0),
            fraction_radiant: parse_f64(object.field(7)).unwrap_or(0.3),
            sensible_fraction: parse_f64(object.field(8)).unwrap_or(0.6),
            activity: optional_name(object.field(9)),
        },
    })
}

fn parse_lights(object: &IdfObject) -> Result<GainDraft, String> {
    let level = design_level(object, 3, 4, 5, 6)?;
    let ret = parse_f64(object.field(7)).unwrap_or(0.0);
    let rad = parse_f64(object.field(8)).unwrap_or(0.0);
    let vis = parse_f64(object.field(9)).unwrap_or(0.0);
    let conv = (1.0 - ret - rad - vis).clamp(0.0, 1.0);
    // Return-air heat has nowhere else to go until plenums exist, so it stays in the zone air.
    Ok(GainDraft {
        zone: object.field(1).to_string(),
        spec: GainSpec::Power {
            schedule: object.field(2).to_string(),
            convective_w: level * (conv + ret),
            radiant_w: level * (rad + vis),
            latent_w: 0.0,
        },
    })
}

fn parse_equipment(object: &IdfObject, zone_index: usize) -> Result<GainDraft, String> {
    let method = zone_index + 2;
    let level_index = zone_index + 3;
    let level = design_level(
        object,
        method,
        level_index,
        level_index + 1,
        level_index + 2,
    )?;
    let lat = parse_f64(object.field(level_index + 3)).unwrap_or(0.0);
    let rad = parse_f64(object.field(level_index + 4)).unwrap_or(0.0);
    let lost = parse_f64(object.field(level_index + 5)).unwrap_or(0.0);
    let conv = (1.0 - lat - rad - lost).clamp(0.0, 1.0);
    Ok(GainDraft {
        zone: object.field(zone_index).to_string(),
        spec: GainSpec::Power {
            schedule: object.field(zone_index + 1).to_string(),
            convective_w: level * conv,
            radiant_w: level * rad,
            latent_w: level * lat,
        },
    })
}

fn design_level(
    object: &IdfObject,
    method_i: usize,
    level_i: usize,
    per_area_i: usize,
    per_person_i: usize,
) -> Result<f64, String> {
    let method = object.field(method_i);
    if method.eq_ignore_ascii_case("Watts/Area") || method.eq_ignore_ascii_case("Power/Area") {
        Err(format!(
            "{} '{}' uses {} which is applied later from floor area; enter a design level for this slice or use Watts",
            object.object_type,
            object.field(0),
            method
        ))
    } else if method.eq_ignore_ascii_case("Watts/Person") {
        let _ = per_person_i;
        Err(format!(
            "{} '{}' Watts/Person is not ported yet",
            object.object_type,
            object.field(0)
        ))
    } else {
        let _ = per_area_i;
        Ok(parse_f64(object.field(level_i)).unwrap_or(0.0))
    }
}

fn parse_ideal(object: &IdfObject) -> IdealDraft {
    let heat_limit = object.field(9);
    let cool_limit = object.field(12);
    let heat_cap = if heat_limit.eq_ignore_ascii_case("LimitCapacity")
        || heat_limit.eq_ignore_ascii_case("LimitFlowRateAndCapacity")
    {
        parse_f64(object.field(11)).unwrap_or(f64::INFINITY)
    } else {
        f64::INFINITY
    };
    let cool_cap = if cool_limit.eq_ignore_ascii_case("LimitCapacity")
        || cool_limit.eq_ignore_ascii_case("LimitFlowRateAndCapacity")
    {
        parse_f64(object.field(14)).unwrap_or(f64::INFINITY)
    } else {
        f64::INFINITY
    };
    IdealDraft {
        name: object.field(0).to_string(),
        heating_capacity_w: heat_cap,
        cooling_capacity_w: cool_cap,
    }
}

fn parse_baseboard(object: &IdfObject) -> Result<BaseboardDraft, String> {
    let name = object.field(0).to_string();
    let method = {
        let text = object.field(2).trim();
        if text.is_empty() {
            "HeatingDesignCapacity"
        } else {
            text
        }
    };
    let capacity = if method.eq_ignore_ascii_case("HeatingDesignCapacity") {
        let watts = parse_f64(object.field(3)).ok_or_else(|| {
            format!(
                "ZoneHVAC:Baseboard:Convective:Electric '{name}' needs a numeric Heating Design Capacity. Autosize is not ported yet"
            )
        })?;
        if watts < 0.0 {
            return Err(format!(
                "ZoneHVAC:Baseboard:Convective:Electric '{name}' heating capacity {watts} W is negative"
            ));
        }
        BaseboardCapacity::Watts(watts)
    } else if method.eq_ignore_ascii_case("CapacityPerFloorArea") {
        let per_area = parse_f64(object.field(4)).ok_or_else(|| {
            format!(
                "ZoneHVAC:Baseboard:Convective:Electric '{name}' needs Heating Design Capacity Per Floor Area"
            )
        })?;
        if per_area <= 0.0 {
            return Err(format!(
                "ZoneHVAC:Baseboard:Convective:Electric '{name}' capacity per floor area must be greater than 0"
            ));
        }
        BaseboardCapacity::PerFloorArea(per_area)
    } else if method.eq_ignore_ascii_case("FractionOfAutosizedHeatingCapacity") {
        return Err(format!(
            "ZoneHVAC:Baseboard:Convective:Electric '{name}' uses FractionOfAutosizedHeatingCapacity, which needs zone sizing. Enter HeatingDesignCapacity or CapacityPerFloorArea"
        ));
    } else {
        return Err(format!(
            "ZoneHVAC:Baseboard:Convective:Electric '{name}' heating capacity method '{method}' is not recognized"
        ));
    };
    let efficiency = match object.field(6).trim() {
        "" => 1.0,
        text => {
            let efficiency = parse_f64(text).ok_or_else(|| {
                format!(
                    "ZoneHVAC:Baseboard:Convective:Electric '{name}' efficiency '{text}' is not a number"
                )
            })?;
            if efficiency <= 0.0 || efficiency > 1.0 {
                return Err(format!(
                    "ZoneHVAC:Baseboard:Convective:Electric '{name}' efficiency {efficiency} must be greater than 0 and at most 1"
                ));
            }
            efficiency
        }
    };
    Ok(BaseboardDraft {
        name,
        availability: optional_name(object.field(1)),
        capacity,
        efficiency,
    })
}

fn parse_equip_list(object: &IdfObject) -> Result<EquipList, String> {
    // Load distribution scheme is field 2 when present. Equipment groups follow.
    let mut start = 1;
    let mut scheme = "SequentialLoad".to_string();
    if object.fields.len() > 2 {
        let maybe = object.field(1);
        if maybe.eq_ignore_ascii_case("SequentialLoad")
            || maybe.eq_ignore_ascii_case("UniformLoad")
            || maybe.eq_ignore_ascii_case("UniformPLR")
            || maybe.eq_ignore_ascii_case("SequentialUniformPLR")
        {
            scheme = maybe.to_string();
            start = 2;
        }
    }
    if !scheme.eq_ignore_ascii_case("SequentialLoad") {
        return Err(format!(
            "equipment list '{}' uses load distribution '{scheme}'. Only SequentialLoad is ported",
            object.field(0)
        ));
    }
    let mut items = Vec::new();
    let mut i = start;
    while i < object.fields.len() {
        if object.field(i).is_empty() {
            i += 1;
            continue;
        }
        if !object.field(i).contains(':') {
            return Err(format!(
                "equipment list '{}' has '{}' where an equipment type was expected",
                object.field(0),
                object.field(i)
            ));
        }
        if i + 1 >= object.fields.len() || object.field(i + 1).is_empty() {
            return Err(format!(
                "equipment list '{}' is missing a name for {}",
                object.field(0),
                object.field(i)
            ));
        }
        let kind = object.field(i).to_string();
        let name = object.field(i + 1).to_string();
        let end = (i + 2..object.fields.len())
            .find(|j| object.field(*j).contains(':'))
            .unwrap_or(object.fields.len());
        let extras: Vec<&str> = (i + 2..end).map(|j| object.field(j)).collect();
        let cooling_sequence = sequence_field(extras.first().copied().unwrap_or(""), &kind)?;
        let heating_sequence = sequence_field(extras.get(1).copied().unwrap_or(""), &kind)?;
        items.push(EquipItem {
            kind,
            name,
            cooling_sequence,
            heating_sequence,
        });
        i = end;
    }
    if items.is_empty() {
        return Err(format!("equipment list '{}' is empty", object.field(0)));
    }
    Ok(EquipList {
        name: object.field(0).to_string(),
        items,
    })
}

fn sequence_field(text: &str, equipment: &str) -> Result<i32, String> {
    let text = text.trim();
    if text.is_empty() {
        return Ok(1);
    }
    let value = parse_f64(text).ok_or_else(|| {
        format!("{equipment} sequence '{text}' is not a number. Cooling and heating sequences come before the fraction schedules")
    })?;
    Ok(value.round() as i32)
}

fn required_f64(object: &IdfObject, index: usize, what: &str) -> Result<f64, String> {
    parse_f64(object.field(index))
        .ok_or_else(|| format!("{what} is missing on '{}'", object.field(0)))
}

fn parse_output_variable(object: &IdfObject) -> Result<OutputVariable, String> {
    let variable = object.field(1).trim().to_string();
    if variable.is_empty() {
        return Err("Output:Variable is missing a variable name".to_string());
    }
    if !output_variable_known(&variable) {
        return Err(format!(
            "Output:Variable '{variable}' is not produced by this port yet"
        ));
    }
    let key = {
        let raw = object.field(0).trim();
        if raw.is_empty() {
            "*".to_string()
        } else {
            raw.to_string()
        }
    };
    Ok(OutputVariable {
        key,
        variable,
        frequency: parse_report_frequency(object.field(2)),
    })
}

fn parse_report_frequency(text: &str) -> ReportFrequency {
    match text.trim().to_ascii_lowercase().as_str() {
        "detailed" | "timestep" => ReportFrequency::Timestep,
        "daily" => ReportFrequency::Daily,
        "monthly" => ReportFrequency::Monthly,
        "runperiod" | "environment" | "annual" => ReportFrequency::RunPeriod,
        _ => ReportFrequency::Hourly,
    }
}

pub fn output_variable_known(name: &str) -> bool {
    matches!(
        name.trim().to_ascii_lowercase().as_str(),
        "site outdoor air drybulb temperature"
            | "zone mean air temperature"
            | "zone air humidity ratio"
            | "zone ideal loads supply air sensible heating rate"
            | "zone ideal loads supply air sensible cooling rate"
            | "baseboard total heating rate"
            | "baseboard electricity rate"
    )
}

fn optional_name(text: &str) -> Option<String> {
    if text.trim().is_empty() {
        None
    } else {
        Some(text.to_string())
    }
}

struct Material {
    name: String,
    thickness: f64,
    conductivity: f64,
    density: f64,
    specific_heat: f64,
    thermal_absorptance: f64,
    solar_absorptance: f64,
    roughness: Roughness,
}

struct SimpleGlazing {
    name: String,
    u_factor: f64,
    shgc: f64,
}

struct ZoneDraft {
    name: String,
    rel_north: f64,
    origin: [f64; 3],
    multiplier: f64,
    ceiling_height: Option<f64>,
    volume: Option<f64>,
    floor_area: Option<f64>,
}

struct SurfaceDraft {
    name: String,
    construction: String,
    zone: String,
    boundary: String,
    boundary_object: String,
    vertices: Vec<[f64; 3]>,
}

struct WindowDraft {
    name: String,
    construction: String,
    parent: String,
    multiplier: f64,
    vertices: Vec<[f64; 3]>,
}

struct GainDraft {
    zone: String,
    spec: GainSpec,
}

enum GainSpec {
    People {
        occupancy: String,
        activity: Option<String>,
        count: f64,
        per_area: f64,
        area_per: f64,
        method: String,
        fraction_radiant: f64,
        sensible_fraction: f64,
    },
    Power {
        schedule: String,
        convective_w: f64,
        radiant_w: f64,
        latent_w: f64,
    },
}

struct InfDraft {
    zone: String,
    schedule: Option<String>,
    method: String,
    flow: Option<f64>,
    per_area: Option<f64>,
    per_exterior: Option<f64>,
    ach: Option<f64>,
}

struct ThermostatDraft {
    zone: String,
    control_type: String,
    control_name: String,
}

struct DualDraft {
    name: String,
    heating: String,
    cooling: String,
}

struct IdealDraft {
    name: String,
    heating_capacity_w: f64,
    cooling_capacity_w: f64,
}

enum BaseboardCapacity {
    Watts(f64),
    PerFloorArea(f64),
}

struct BaseboardDraft {
    name: String,
    availability: Option<String>,
    capacity: BaseboardCapacity,
    efficiency: f64,
}

struct EquipList {
    name: String,
    items: Vec<EquipItem>,
}

struct EquipItem {
    kind: String,
    name: String,
    cooling_sequence: i32,
    heating_sequence: i32,
}

struct EquipConn {
    zone: String,
    list: String,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compact_schedule_weekday_and_all_day() {
        let idf = r#"
Schedule:Compact,
  Occ, Any Number,
  Through: 12/31,
  For: Weekdays,
  Until: 08:00, 0,
  Until: 18:00, 1,
  Until: 24:00, 0,
  For: AllOtherDays,
  Until: 24:00, 0;
"#;
        let model_schedules = {
            let objects = crate::idf::parse_idf(idf);
            parse_compact(&objects[0]).unwrap()
        };
        // Friday 1 Jan 2021, 09:00 -> 9*60 = 540, within 08:00-18:00.
        assert!((model_schedules.value(2021, 1, 1, 9 * 60) - 1.0).abs() < 1e-9);
        assert!((model_schedules.value(2021, 1, 1, 7 * 60) - 0.0).abs() < 1e-9);
        // Saturday
        assert!((model_schedules.value(2021, 1, 2, 12 * 60) - 0.0).abs() < 1e-9);
    }

    #[test]
    fn south_wall_normal_points_south() {
        let idf = r#"
GlobalGeometryRules, UpperLeftCorner, CounterClockwise, World;
Zone, Box, 0, 0, 0, 0, 1, 1, 3, 30, 10;
Material, Ins, Rough, 0.1, 0.04, 30, 1400, 0.9, 0.6, 0.6;
Construction, Wall, Ins;
BuildingSurface:Detailed,
  SouthWall, Wall, Wall, Box, , Outdoors, , SunExposed, WindExposed, 0.5, 4,
  0, 0, 3,
  0, 0, 0,
  10, 0, 0,
  10, 0, 3;
"#;
        let model = load_idf(idf).unwrap();
        let wall = &model.surfaces[0];
        assert!((wall.area - 30.0).abs() < 1e-6, "area {}", wall.area);
        assert!(wall.normal[1] < -0.99, "normal {:?}", wall.normal);
        assert!((model.zones[0].volume - 30.0).abs() < 1e-6);
    }

    #[test]
    fn output_variable_request_is_kept() {
        let idf = r#"
Zone, Box, 0, 0, 0, 0, 1, 1, 3, 30, 10;
Output:Variable, *, Zone Mean Air Temperature, Timestep;
Output:Variable, *, Not A Real Variable, Hourly;
"#;
        let model = load_idf(idf).unwrap();
        assert_eq!(model.outputs.len(), 1);
        assert_eq!(model.outputs[0].variable, "Zone Mean Air Temperature");
        assert!(model
            .warnings
            .iter()
            .any(|warning| warning.contains("Not A Real Variable")));
    }

    #[test]
    fn electric_baseboard_capacity_per_area_and_heating_sequence() {
        let idf = r#"
Schedule:Constant, On, Any Number, 1;
Zone, Box, 0, 0, 0, 0, 1, 1, 3, 30, 10;
ZoneHVAC:Baseboard:Convective:Electric,
  BB, On, CapacityPerFloorArea, , 20, , 0.8;
ZoneHVAC:IdealLoadsAirSystem, Ideal;
ZoneHVAC:EquipmentList,
  List, SequentialLoad,
  ZoneHVAC:Baseboard:Convective:Electric, BB, 1, 1, , ,
  ZoneHVAC:IdealLoadsAirSystem, Ideal, 1, 2;
ZoneHVAC:EquipmentConnections,
  Box, List;
"#;
        let model = load_idf(idf).unwrap();
        let zone = &model.zones[0];
        assert!(zone.equipment.len() == 2, "{:?}", zone.equipment.len());
        match &zone.equipment[0] {
            ZoneEquipment::ElectricBaseboard {
                name,
                capacity_w,
                efficiency,
                availability,
                heating_sequence,
            } => {
                assert_eq!(name, "BB");
                assert!((capacity_w - 200.0).abs() < 1e-9, "capacity {capacity_w}");
                assert!((efficiency - 0.8).abs() < 1e-12);
                assert_eq!(*availability, Some(0));
                assert_eq!(*heating_sequence, 1);
            }
            other => panic!("expected baseboard, got {other:?}"),
        }
        match &zone.equipment[1] {
            ZoneEquipment::IdealLoads {
                heating_sequence,
                cooling_sequence,
                heating_capacity_w,
                ..
            } => {
                assert_eq!(*heating_sequence, 2);
                assert_eq!(*cooling_sequence, 1);
                assert!(heating_capacity_w.is_infinite());
            }
            other => panic!("expected ideal loads, got {other:?}"),
        }
    }

    #[test]
    fn electric_baseboard_autosize_is_rejected() {
        let idf = r#"
Zone, Box, 0, 0, 0, 0, 1, 1, 3, 30, 10;
ZoneHVAC:Baseboard:Convective:Electric, BB, , HeatingDesignCapacity, autosize;
ZoneHVAC:EquipmentList, List, SequentialLoad, ZoneHVAC:Baseboard:Convective:Electric, BB, 1, 1;
ZoneHVAC:EquipmentConnections, Box, List;
"#;
        let err = load_idf(idf).unwrap_err();
        assert!(err.contains("Autosize"), "{err}");
    }
}
