# Potential Features

This document captures feature ideas that are not implemented today but would be
valuable in future updates.

## Core Calculation Enhancements
- Add atmospheric loss, rain attenuation, and propagation models (ITU-R based).
- Support bistatic and multistatic radar range equations alongside monostatic.
- Provide configurable noise temperature and receiver noise figure modeling.
- Introduce probability of detection/false alarm calculations with Swerling models.
- Add pulse integration and coherent processing interval utilities for common radar modes.

## Data & Scenario Modeling
- Build a scenario API for multiple targets, clutter sources, and jamming emitters.
- Provide time-series simulation helpers for range/Doppler tracks.
- Add parameter sweeps and sensitivity analysis utilities for design trade studies.

## Units, Validation, and UX
- Offer unit-aware calculations (e.g., Pint/Dart extensions/Rust newtypes) to reduce errors.
- Provide structured input validation with clear error reporting in all languages.
- Add configuration loaders for JSON/YAML radar system definitions.

## Visualization & Reporting
- Export plots and summaries in standardized report formats (PDF/HTML).
- Provide interactive notebooks or dashboard templates for common radar workflows.

## Cross-Language Consistency
- Generate language bindings from a shared equation metadata catalog.
- Ensure equation coverage parity tests across Python, Dart, and Rust.
