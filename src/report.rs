//! Requested `Output:Variable` table, in the EnergyPlus CSV column style.
//!
//! Rates and state variables are averaged over the reporting interval. The
//! Date/Time column is the end of that interval, with hour 24 written as 24:00:00.

use crate::model::{OutputVariable, ReportFrequency};
use crate::simulate::{BaseboardStep, Simulation, TimestepResult, ZoneStep};

pub fn write_report(requests: &[OutputVariable], simulation: &Simulation) -> String {
    if simulation.steps.is_empty() || requests.is_empty() {
        return String::new();
    }
    let columns = expand_columns(requests, &simulation.steps[0]);
    let mut csv = String::from("Date/Time");
    for column in &columns {
        csv.push(',');
        csv.push_str(&column.header);
    }
    csv.push('\n');
    for bucket in group_steps(&simulation.steps, &columns) {
        csv.push_str(&format_timestamp(
            bucket.month,
            bucket.day,
            bucket.end_minute,
        ));
        for value in bucket.values {
            csv.push(',');
            csv.push_str(&format_value(value));
        }
        csv.push('\n');
    }
    csv
}

struct Column {
    header: String,
    frequency: ReportFrequency,
    sample: Sample,
}

enum Sample {
    Outdoor,
    Zone {
        index: usize,
        field: ZoneField,
    },
    Baseboard {
        zone: usize,
        index: usize,
        field: BaseboardField,
    },
}

#[derive(Clone, Copy)]
enum ZoneField {
    Temperature,
    Humidity,
    Heating,
    Cooling,
}

#[derive(Clone, Copy)]
enum BaseboardField {
    Heating,
    Electric,
}

struct Bucket {
    month: u32,
    day: u32,
    end_minute: u32,
    values: Vec<f64>,
}

fn expand_columns(requests: &[OutputVariable], first: &TimestepResult) -> Vec<Column> {
    let mut columns = Vec::new();
    for request in requests {
        let frequency = request.frequency;
        let variable = request.variable.trim();
        let lower = variable.to_ascii_lowercase();
        if lower == "site outdoor air drybulb temperature" {
            let key = if request.key == "*" || request.key.eq_ignore_ascii_case("environment") {
                "ENVIRONMENT"
            } else {
                request.key.as_str()
            };
            columns.push(Column {
                header: header(key, variable, "C", frequency),
                frequency,
                sample: Sample::Outdoor,
            });
            continue;
        }
        if let Some(field) = baseboard_field(&lower) {
            for (zone_index, zone) in first.zones.iter().enumerate() {
                for (index, baseboard) in zone.baseboards.iter().enumerate() {
                    if request.key != "*" && !baseboard.name.eq_ignore_ascii_case(&request.key) {
                        continue;
                    }
                    columns.push(Column {
                        header: header(&baseboard.name, variable, "W", frequency),
                        frequency,
                        sample: Sample::Baseboard {
                            zone: zone_index,
                            index,
                            field,
                        },
                    });
                }
            }
            continue;
        }
        let field = match lower.as_str() {
            "zone mean air temperature" => ZoneField::Temperature,
            "zone air humidity ratio" => ZoneField::Humidity,
            "zone ideal loads supply air sensible heating rate" => ZoneField::Heating,
            "zone ideal loads supply air sensible cooling rate" => ZoneField::Cooling,
            _ => continue,
        };
        let zones: Vec<usize> = if request.key == "*" {
            (0..first.zones.len()).collect()
        } else {
            first
                .zones
                .iter()
                .enumerate()
                .filter(|(_, zone)| zone.name.eq_ignore_ascii_case(&request.key))
                .map(|(index, _)| index)
                .collect()
        };
        for index in zones {
            columns.push(Column {
                header: header(&first.zones[index].name, variable, unit(field), frequency),
                frequency,
                sample: Sample::Zone { index, field },
            });
        }
    }
    columns
}

fn header(key: &str, variable: &str, unit: &str, frequency: ReportFrequency) -> String {
    let frequency = match frequency {
        ReportFrequency::Timestep => "Timestep",
        ReportFrequency::Hourly => "Hourly",
        ReportFrequency::Daily => "Daily",
        ReportFrequency::Monthly => "Monthly",
        ReportFrequency::RunPeriod => "RunPeriod",
    };
    format!("{key}:{variable} [{unit}]({frequency})")
}

fn baseboard_field(name: &str) -> Option<BaseboardField> {
    match name {
        "baseboard total heating rate" => Some(BaseboardField::Heating),
        "baseboard electricity rate" => Some(BaseboardField::Electric),
        _ => None,
    }
}

fn unit(field: ZoneField) -> &'static str {
    match field {
        ZoneField::Temperature => "C",
        ZoneField::Humidity => "kgWater/kgDryAir",
        ZoneField::Heating | ZoneField::Cooling => "W",
    }
}

fn group_steps(steps: &[TimestepResult], columns: &[Column]) -> Vec<Bucket> {
    if columns.is_empty() {
        return Vec::new();
    }
    // Columns can request different frequencies. Each output row follows the finest
    // frequency present, and coarser columns repeat their interval average.
    let finest = columns
        .iter()
        .map(|column| column.frequency)
        .min_by_key(|frequency| rank(*frequency))
        .unwrap();
    let mut rows = Vec::new();
    let mut index = 0;
    while index < steps.len() {
        let start = index;
        index = next_group_end(steps, index, finest);
        let group = &steps[start..index];
        let end = group.last().unwrap();
        let mut values = Vec::with_capacity(columns.len());
        for column in columns {
            let sample = match column.frequency {
                freq if freq == finest => group,
                _ => interval_slice(steps, end, column.frequency),
            };
            values.push(average(sample, &column.sample));
        }
        rows.push(Bucket {
            month: end.month,
            day: end.day,
            end_minute: end_minute_for(end, finest),
            values,
        });
    }
    rows
}

fn rank(frequency: ReportFrequency) -> u8 {
    match frequency {
        ReportFrequency::Timestep => 0,
        ReportFrequency::Hourly => 1,
        ReportFrequency::Daily => 2,
        ReportFrequency::Monthly => 3,
        ReportFrequency::RunPeriod => 4,
    }
}

fn next_group_end(steps: &[TimestepResult], start: usize, frequency: ReportFrequency) -> usize {
    let first = &steps[start];
    let mut end = start + 1;
    while end < steps.len() && same_interval(first, &steps[end], frequency) {
        end += 1;
    }
    end
}

fn same_interval(
    left: &TimestepResult,
    right: &TimestepResult,
    frequency: ReportFrequency,
) -> bool {
    match frequency {
        ReportFrequency::Timestep => false,
        ReportFrequency::Hourly => {
            left.year == right.year
                && left.month == right.month
                && left.day == right.day
                && left.hour == right.hour
        }
        ReportFrequency::Daily => {
            left.year == right.year && left.month == right.month && left.day == right.day
        }
        ReportFrequency::Monthly => left.year == right.year && left.month == right.month,
        ReportFrequency::RunPeriod => true,
    }
}

fn interval_slice<'a>(
    steps: &'a [TimestepResult],
    end: &TimestepResult,
    frequency: ReportFrequency,
) -> &'a [TimestepResult] {
    let start = steps
        .iter()
        .position(|step| same_interval(step, end, frequency))
        .unwrap_or(0);
    let stop = steps
        .iter()
        .rposition(|step| same_interval(step, end, frequency))
        .unwrap_or(start)
        + 1;
    &steps[start..stop]
}

fn average(steps: &[TimestepResult], sample: &Sample) -> f64 {
    if steps.is_empty() {
        return 0.0;
    }
    let sum: f64 = steps
        .iter()
        .map(|step| match sample {
            Sample::Outdoor => step.outdoor_db,
            Sample::Zone { index, field } => zone_value(&step.zones[*index], *field),
            Sample::Baseboard { zone, index, field } => step.zones[*zone]
                .baseboards
                .get(*index)
                .map(|baseboard| baseboard_value(baseboard, *field))
                .unwrap_or(0.0),
        })
        .sum();
    sum / steps.len() as f64
}

fn zone_value(zone: &ZoneStep, field: ZoneField) -> f64 {
    match field {
        ZoneField::Temperature => zone.temperature,
        ZoneField::Humidity => zone.humidity_ratio,
        ZoneField::Heating => zone.ideal_heating_w,
        ZoneField::Cooling => zone.cooling_w,
    }
}

fn baseboard_value(baseboard: &BaseboardStep, field: BaseboardField) -> f64 {
    match field {
        BaseboardField::Heating => baseboard.heating_w,
        BaseboardField::Electric => baseboard.electric_w,
    }
}

fn end_minute_for(step: &TimestepResult, frequency: ReportFrequency) -> u32 {
    match frequency {
        ReportFrequency::Timestep => step.end_minute,
        ReportFrequency::Hourly => step.hour * 60,
        ReportFrequency::Daily | ReportFrequency::Monthly | ReportFrequency::RunPeriod => 24 * 60,
    }
}

fn format_timestamp(month: u32, day: u32, end_minute: u32) -> String {
    let hour = end_minute / 60;
    let minute = end_minute % 60;
    format!("{month:02}/{day:02}  {hour:02}:{minute:02}:00")
}

fn format_value(value: f64) -> String {
    if value.abs() >= 100.0 || value == 0.0 {
        format!("{value:.4}")
    } else {
        format!("{value:.6}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{OutputVariable, ReportFrequency};
    use crate::simulate::{BaseboardStep, Simulation, TimestepResult, ZoneStep};

    fn step(hour: u32, temp: f64, heat: f64) -> TimestepResult {
        TimestepResult {
            year: 2021,
            month: 1,
            day: 1,
            hour,
            end_minute: hour * 60,
            outdoor_db: -2.0,
            zones: vec![ZoneStep {
                name: "Living".into(),
                temperature: temp,
                humidity_ratio: 0.004,
                heating_w: heat,
                cooling_w: 0.0,
                ideal_heating_w: heat,
                baseboards: vec![],
            }],
        }
    }

    #[test]
    fn hourly_request_averages_the_zone_timestep_values() {
        let simulation = Simulation {
            steps: vec![
                TimestepResult {
                    end_minute: 10,
                    ..step(1, 20.0, 100.0)
                },
                TimestepResult {
                    end_minute: 20,
                    ..step(1, 22.0, 300.0)
                },
            ],
            heating_j: 0.0,
            cooling_j: 0.0,
        };
        let csv = write_report(
            &[OutputVariable {
                key: "*".into(),
                variable: "Zone Mean Air Temperature".into(),
                frequency: ReportFrequency::Hourly,
            }],
            &simulation,
        );
        assert!(csv.starts_with("Date/Time,Living:Zone Mean Air Temperature [C](Hourly)\n"));
        assert!(csv.contains("01/01  01:00:00,21.000000\n"), "{csv}");
    }

    #[test]
    fn baseboard_electricity_is_keyed_by_the_equipment_name() {
        let mut row = step(1, 20.0, 0.0);
        row.end_minute = 10;
        row.zones[0].baseboards = vec![BaseboardStep {
            name: "BB".into(),
            heating_w: 500.0,
            electric_w: 1000.0,
        }];
        let csv = write_report(
            &[OutputVariable {
                key: "*".into(),
                variable: "Baseboard Electricity Rate".into(),
                frequency: ReportFrequency::Timestep,
            }],
            &Simulation {
                steps: vec![row],
                heating_j: 0.0,
                cooling_j: 0.0,
            },
        );
        assert!(
            csv.contains("BB:Baseboard Electricity Rate [W](Timestep)"),
            "{csv}"
        );
        assert!(csv.contains(",1000"), "{csv}");
    }
}
