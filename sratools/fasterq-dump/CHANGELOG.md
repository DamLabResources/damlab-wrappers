# Changelog

All notable changes to the SRA Toolkit Fasterq-dump Wrapper will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [1.0.1] - 2025-04-21
 - Added support for downloading SRA accessions directly without requiring `prefetch`

## [1.0.0] - 2025-04-17
### Added
- Initial release
- Support for converting SRA files to FASTQ format
- Comprehensive parameter support including:
  - Output format specification
  - Buffer and cache size tuning
  - Memory limits
  - Thread count control
  - Progress monitoring
  - Split options (spot/files/3)
  - Technical read handling
  - Minimum read length filtering
- Performance optimization options
- Test suite with different output format scenarios 