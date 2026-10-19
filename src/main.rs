use std::env;
use std::fs;
use std::process::ExitCode;

use bem::epw::parse_epw;
use bem::model::load_idf;
use bem::simulate::{simulate, SimOptions};

fn main() -> ExitCode {
    if let Err(err) = run() {
        eprintln!("error: {err}");
        return ExitCode::from(1);
    }
    ExitCode::SUCCESS
}

fn run() -> Result<(), String> {
    let mut args = env::args().skip(1);
    let command = args.next().unwrap_or_else(|| "help".to_string());
    match command.as_str() {
        "run" => {
            let idf_path = args
                .next()
                .ok_or("usage: bem run <file.idf|file.epJSON> <file.epw> [--out results.csv]")?;
            let epw_path = args
                .next()
                .ok_or("usage: bem run <file.idf|file.epJSON> <file.epw> [--out results.csv]")?;
            let mut out_path = None;
            let mut warmup_days = 7u32;
            while let Some(flag) = args.next() {
                match flag.as_str() {
                    "--out" => {
                        out_path = Some(args.next().ok_or("--out needs a path")?);
                    }
                    "--warmup-days" => {
                        let text = args.next().ok_or("--warmup-days needs an integer")?;
                        warmup_days = text
                            .parse()
                            .map_err(|_| format!("bad warmup days '{text}'"))?;
                    }
                    other => return Err(format!("unknown flag '{other}'")),
                }
            }
            let idf = fs::read_to_string(&idf_path)
                .map_err(|err| format!("reading {idf_path}: {err}"))?;
            let epw = fs::read_to_string(&epw_path)
                .map_err(|err| format!("reading {epw_path}: {err}"))?;
            let model = load_idf(&idf)?;
            for warning in &model.warnings {
                eprintln!("warning: {warning}");
            }
            let weather = parse_epw(&epw)?;
            let simulation = simulate(&model, &weather, &SimOptions { warmup_days })?;
            let csv = if model.outputs.is_empty() {
                hourly_csv(&simulation)
            } else {
                bem::report::write_report(&model.outputs, &simulation)
            };
            match out_path {
                Some(path) => {
                    fs::write(&path, &csv).map_err(|err| format!("writing {path}: {err}"))?;
                    println!("wrote {path}");
                }
                None => print!("{csv}"),
            }
            let hours = simulation.heating_j / 3_600_000.0;
            let cool = simulation.cooling_j / 3_600_000.0;
            eprintln!("heating {hours:.3} kWh");
            eprintln!("cooling {cool:.3} kWh");
            eprintln!("timesteps {}", simulation.steps.len());
            Ok(())
        }
        "help" | "--help" | "-h" => {
            println!(
                "bem — Rust port of the EnergyPlus engine\n\n\
                 usage:\n  \
                 bem run <file.idf|file.epJSON> <file.epw> [--out results.csv] [--warmup-days 7]\n\n\
                 This is not EnergyPlus. Ideal loads, convective electric baseboard,\n\
                 zone heat balance, psychrometrics, and a finite-difference wall are\n\
                 implemented. Input may be IDF or epJSON. Other HVAC systems are rejected."
            );
            Ok(())
        }
        other => Err(format!("unknown command '{other}'. Try `bem help`.")),
    }
}

fn hourly_csv(simulation: &bem::simulate::Simulation) -> String {
    if simulation.steps.is_empty() {
        return String::new();
    }
    let zone_names: Vec<&str> = simulation.steps[0]
        .zones
        .iter()
        .map(|z| z.name.as_str())
        .collect();
    let mut header = "Date,Hour,Outdoor Drybulb (C)".to_string();
    for name in &zone_names {
        header.push_str(&format!(
            ",{name} Temperature (C),{name} Humidity Ratio (kg/kg),{name} Heating (W),{name} Cooling (W)"
        ));
    }
    header.push('\n');
    let mut body = header;
    let mut index = 0;
    while index < simulation.steps.len() {
        let step = &simulation.steps[index];
        let mut count = 0.0;
        let mut outdoor = 0.0;
        let mut temp = vec![0.0; zone_names.len()];
        let mut hum = vec![0.0; zone_names.len()];
        let mut heat = vec![0.0; zone_names.len()];
        let mut cool = vec![0.0; zone_names.len()];
        while index < simulation.steps.len()
            && simulation.steps[index].year == step.year
            && simulation.steps[index].month == step.month
            && simulation.steps[index].day == step.day
            && simulation.steps[index].hour == step.hour
        {
            let row = &simulation.steps[index];
            outdoor += row.outdoor_db;
            for (zi, zone) in row.zones.iter().enumerate() {
                temp[zi] += zone.temperature;
                hum[zi] += zone.humidity_ratio;
                heat[zi] += zone.heating_w;
                cool[zi] += zone.cooling_w;
            }
            count += 1.0;
            index += 1;
        }
        body.push_str(&format!(
            "{:04}-{:02}-{:02},{},{:.4}",
            step.year,
            step.month,
            step.day,
            step.hour,
            outdoor / count
        ));
        for zi in 0..zone_names.len() {
            body.push_str(&format!(
                ",{:.4},{:.6},{:.2},{:.2}",
                temp[zi] / count,
                hum[zi] / count,
                heat[zi] / count,
                cool[zi] / count
            ));
        }
        body.push('\n');
    }
    body
}
