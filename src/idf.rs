//! EnergyPlus IDF reader.
//!
//! Strips `!` comments, then splits objects on `;` and fields on `,`.
//! This is the same field grammar as the EnergyPlus input processor.

#[derive(Clone, Debug)]
pub struct IdfObject {
    pub object_type: String,
    pub fields: Vec<String>,
}

impl IdfObject {
    pub fn is_type(&self, name: &str) -> bool {
        self.object_type.eq_ignore_ascii_case(name)
    }

    pub fn field(&self, index: usize) -> &str {
        self.fields.get(index).map(String::as_str).unwrap_or("")
    }
}

pub fn parse_idf(text: &str) -> Vec<IdfObject> {
    let mut stripped = String::with_capacity(text.len());
    for line in text.split_inclusive(['\n', '\r']) {
        let content = match line.find('!') {
            Some(i) => &line[..i],
            None => line,
        };
        stripped.push_str(content);
        if line.ends_with('\n') || line.ends_with('\r') {
            stripped.push('\n');
        }
    }
    let mut objects = Vec::new();
    for chunk in stripped.split(';') {
        let mut fields: Vec<String> = chunk
            .split(',')
            .map(|field| field.split_whitespace().collect::<Vec<_>>().join(" "))
            .collect();
        while fields.last().is_some_and(|f| f.is_empty()) {
            fields.pop();
        }
        if fields.is_empty() || fields.iter().all(|f| f.is_empty()) {
            continue;
        }
        let object_type = fields.remove(0);
        if object_type.is_empty() {
            continue;
        }
        objects.push(IdfObject {
            object_type,
            fields,
        });
    }
    objects
}

pub fn parse_f64(text: &str) -> Option<f64> {
    let trimmed = text.trim();
    if trimmed.is_empty() || trimmed.eq_ignore_ascii_case("autocalculate") {
        return None;
    }
    trimmed.parse::<f64>().ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_comments_and_blank_fields() {
        let idf = r#"
Version, 25.1;
Zone,
  Zone One, !- Name
  0,        !- north
  ,         !- x
  1.5;
"#;
        let objects = parse_idf(idf);
        assert_eq!(objects.len(), 2);
        assert!(objects[0].is_type("Version"));
        assert_eq!(objects[0].field(0), "25.1");
        assert!(objects[1].is_type("Zone"));
        assert_eq!(objects[1].field(0), "Zone One");
        assert_eq!(objects[1].field(1), "0");
        assert_eq!(objects[1].field(2), "");
        assert_eq!(objects[1].field(3), "1.5");
    }
}
