//! EnergyPlus weather (EPW) reader.
//!
//! The first eight lines are the header. Data rows follow the EPW data dictionary.

#[derive(Clone, Debug)]
pub struct WeatherHour {
    pub year: i32,
    pub month: u32,
    pub day: u32,
    pub hour: u32,
    pub dry_bulb: f64,
    pub dew_point: f64,
    pub relative_humidity: f64,
    pub pressure: f64,
    pub horizontal_ir: f64,
    pub direct_normal: f64,
    pub diffuse_horizontal: f64,
    pub wind_direction: f64,
    pub wind_speed: f64,
}

#[derive(Clone, Debug)]
pub struct WeatherFile {
    pub city: String,
    pub latitude: f64,
    pub longitude: f64,
    pub time_zone: f64,
    pub elevation: f64,
    pub hours: Vec<WeatherHour>,
}

pub fn parse_epw(text: &str) -> Result<WeatherFile, String> {
    let lines: Vec<&str> = text
        .lines()
        .filter(|line| !line.trim().is_empty())
        .collect();
    if lines.len() < 9 {
        return Err("EPW file needs an 8-line header and at least one data row".to_string());
    }
    let header: Vec<&str> = lines[0].split(',').map(str::trim).collect();
    if !header
        .first()
        .is_some_and(|key| key.eq_ignore_ascii_case("LOCATION"))
    {
        return Err("EPW first line must be LOCATION".to_string());
    }
    let city = header.get(1).unwrap_or(&"Unknown").to_string();
    let latitude = header
        .get(6)
        .and_then(|v| v.parse().ok())
        .ok_or("EPW LOCATION latitude is missing")?;
    let longitude = header
        .get(7)
        .and_then(|v| v.parse().ok())
        .ok_or("EPW LOCATION longitude is missing")?;
    let time_zone = header.get(8).and_then(|v| v.parse().ok()).unwrap_or(0.0);
    let elevation = header.get(9).and_then(|v| v.parse().ok()).unwrap_or(0.0);

    let mut hours = Vec::new();
    for line in lines.iter().skip(8) {
        let col: Vec<&str> = line.split(',').map(str::trim).collect();
        if col.len() < 22 {
            return Err(format!(
                "EPW data row has {} fields, need at least 22",
                col.len()
            ));
        }
        let num = |index: usize, name: &str| -> Result<f64, String> {
            col[index]
                .parse::<f64>()
                .map_err(|_| format!("EPW {name} is not a number: '{}'", col[index]))
        };
        hours.push(WeatherHour {
            year: num(0, "year")? as i32,
            month: num(1, "month")? as u32,
            day: num(2, "day")? as u32,
            hour: num(3, "hour")? as u32,
            dry_bulb: num(6, "dry bulb")?,
            dew_point: num(7, "dew point")?,
            relative_humidity: num(8, "relative humidity")?,
            pressure: num(9, "pressure")?,
            horizontal_ir: num(12, "horizontal IR")?,
            direct_normal: num(14, "direct normal")?,
            diffuse_horizontal: num(15, "diffuse horizontal")?,
            wind_direction: num(20, "wind direction")?,
            wind_speed: num(21, "wind speed")?,
        });
    }
    if hours.is_empty() {
        return Err("EPW file has no data rows".to_string());
    }
    Ok(WeatherFile {
        city,
        latitude,
        longitude,
        time_zone,
        elevation,
        hours,
    })
}

impl WeatherFile {
    pub fn hour(&self, year: i32, month: u32, day: u32, hour: u32) -> Option<&WeatherHour> {
        self.hours.iter().find(|row| {
            row.year == year && row.month == month && row.day == day && row.hour == hour
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_location_and_one_hour() {
        let text = "\
LOCATION,Denver,CO,USA,TMY,999,39.83,-104.65,-7.0,1608.0
DESIGN CONDITIONS,0
TYPICAL/EXTREME PERIODS,0
GROUND TEMPERATURES,0
HOLIDAYS/DAYLIGHT SAVINGS,No,0,0,0
COMMENTS 1,x
COMMENTS 2,y
DATA PERIODS,1,1,Data,Sunday,1/1,1/1
2021,1,1,1,60,flags,-8.0,-10.0,80,81100,0,0,200,0,0,0,0,0,0,0,0,2.5,0,0,9999,99999,9,999999999,0,0.04,0,88,0.2,0,0
";
        let weather = parse_epw(text).unwrap();
        assert_eq!(weather.city, "Denver");
        assert!((weather.latitude - 39.83).abs() < 1e-9);
        let hour = weather.hour(2021, 1, 1, 1).unwrap();
        assert!((hour.dry_bulb + 8.0).abs() < 1e-9);
        assert!((hour.wind_speed - 2.5).abs() < 1e-9);
    }
}
