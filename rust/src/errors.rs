/// Errors that can occur during radar calculations.
#[derive(Debug, Clone, PartialEq)]
pub enum RadarError {
    /// Invalid input parameter (e.g., negative value where positive expected).
    InvalidParameter(String),
    /// Division by zero or other mathematical error.
    MathError(String),
}

impl std::fmt::Display for RadarError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RadarError::InvalidParameter(msg) => write!(f, "Invalid parameter: {}", msg),
            RadarError::MathError(msg) => write!(f, "Math error: {}", msg),
        }
    }
}

impl std::error::Error for RadarError {}

/// Result type for radar calculations.
pub type Result<T> = std::result::Result<T, RadarError>;
