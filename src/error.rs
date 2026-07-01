pub type PrimerResult<T> = Result<T, String>;

#[derive(Debug)]
pub struct PrimerError;

impl PrimerError {
    pub fn unknown_chemistry(raw: &str) -> String {
        format!("unknown chemistry: {raw}")
    }

    pub fn invalid_grammar(message: impl AsRef<str>) -> String {
        format!("invalid grammar: {}", message.as_ref())
    }

    pub fn invalid_coordinates(message: impl AsRef<str>) -> String {
        format!("invalid sequence coordinates: {}", message.as_ref())
    }

    pub fn rhapsody(message: impl AsRef<str>) -> String {
        format!("rhapsody whitelist error: {}", message.as_ref())
    }
}
