# Changelog

All notable changes to this wrapper will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2025-04-10

### Added
- First stable release of the Picard AddOrReplaceReadGroups wrapper
- Basic read group manipulation functionality
- Support for all Picard read group parameters:
  - Required: ID, LB, PL, PU, SM
  - Optional: CN, DS, DT, FO, KS, PG, PI, PM
- Additional Picard parameters:
  - validation_stringency: Validation stringency control
  - create_index: BAM index creation
  - compression_level: BAM compression control
- Version checking functionality
- Input file validation
- Comprehensive error handling
- Detailed documentation in README.md
- Type hints for better code maintainability

### Changed
- Improved command line argument handling
- Better organized code structure
- Enhanced parameter handling with sensible defaults
- Removed unnecessary shlex.quote usage
- Added explicit type annotations

### Deprecated
- N/A

### Removed
- N/A

### Fixed
- N/A

### Security
- N/A

## [0.0.0] - 2025-04-10

### Added
- Initial development version
- Basic structure and functionality
- Preliminary implementation of AddOrReplaceReadGroups wrapper 