# ADIAS - Automated Detection, Identification, and Analysis System

ADIAS is a comprehensive Python-based pipeline for automated processing of astronomical CCD images, specifically designed for astrometric measurements of celestial objects such as asteroids, comets, and other moving targets.

## Overview

This system provides end-to-end processing capabilities for astronomical observations, from raw CCD images to precise astrometric measurements. The pipeline handles image calibration, stellar detection, catalog matching, and statistical analysis of observation residuals.

## Key Features

### Image Processing Pipeline
- **Automated Image Calibration**: Supports both standard (bias/dark/flat) and algorithmic background subtraction methods
- **Stellar Detection**: Advanced star detection using connected component analysis with configurable thresholds
- **Catalog Matching**: Precise stellar pattern matching using GAIA star catalogs with proper motion corrections
- **Astrometric Solution**: Plate constant determination using polynomial models of varying complexity

### Statistical Analysis
- **Outlier Detection**: Iterative sigma-clipping algorithms for robust statistical analysis
- **Residual Analysis**: Comprehensive O-C (Observed minus Calculated) residual computation
- **Quality Assessment**: Automated evaluation of observation quality and reliability

### Multi-Object Support
- **Simultaneous Processing**: Handles multiple target objects in a single observation
- **Ephemeris Integration**: Supports both IMCCE and JPL ephemeris formats
- **Automated Reporting**: Generates detailed processing reports and visualization plots

## System Architecture

The system is organized into four main processing modules:

### 1. Preprocessing Module (`utils/preprocess.py`)
Handles image calibration and background estimation:
- Standard calibration using bias, dark, and flat field corrections
- Algorithmic calibration using median filtering for background subtraction
- Robust background estimation using iterative sigma-clipping

### 2. Detection Module (`utils/detect.py`)
Implements stellar detection algorithms:
- Connected component labeling for star identification
- Flux-weighted centroid calculation for precise positioning
- Signal-to-noise ratio computation for quality assessment
- Overflow detection for saturated stars

### 3. Matching Module (`utils/match.py`)
Performs astrometric calibration and catalog matching:
- GAIA catalog loading with magnitude filtering
- Proper motion corrections for accurate stellar positions
- Ephemeris interpolation for target object predictions
- Plate constant solution using least-squares fitting
- Coordinate transformations between pixel and celestial coordinates

### 4. Analysis Module (`utils/comoc.py`)
Conducts statistical analysis of observation results:
- Iterative outlier removal using configurable sigma thresholds
- Comprehensive residual statistics computation
- Quality-based data filtering and reporting
- Monthly and summary report generation

## Installation and Dependencies

### Required Python Packages
```bash
pip install numpy pandas scipy astropy astroquery matplotlib pathlib
```

### Key Dependencies
- **NumPy**: Numerical computations and array operations
- **Pandas**: Data manipulation and analysis
- **SciPy**: Scientific computing and statistical functions
- **Astropy**: Astronomical data handling and coordinate transformations
- **Astroquery**: Online astronomical catalog access
- **Matplotlib**: Plotting and visualization

## Configuration

The system uses a comprehensive configuration system defined in `config.py`. The configuration is organized into several categories:

### Preprocessing Configuration
Controls image calibration parameters:
- Calibration method selection (standard vs algorithmic)
- Median filter window sizes for background estimation
- Background subtraction modes

### Detection Configuration
Manages stellar detection parameters:
- Background threshold multipliers for star detection
- Signal-to-noise ratio thresholds
- Centroid calculation methods

### Matching Configuration
Defines astrometric solution parameters:
- Telescope focal length and pixel scale
- Field of view specifications
- Catalog matching tolerances
- Plate model complexity levels

### Analysis Configuration
Sets statistical analysis parameters:
- Outlier detection thresholds
- Residual limits for quality assessment
- Output formatting options

## Usage Examples

### Basic Image Processing
```python
from main import process_single_image
from config import config

# Process a single FITS image
results = process_single_image(
    "path/to/image.fits",
    calibration_data,
    catalogs,
    ephemerides,
    "output_directory"
)
```

### Batch Processing
```python
# Process all images in configured directories
python main.py --output results --debug
```

### Custom Configuration
```python
# Modify detection parameters
config.detect.background_threshold = 3
config.detect.signal_noise_threshold = 2.0

# Adjust matching parameters
config.match.match_limit = 3.0
config.match.model_type = 12  # Use 12-parameter plate model
```

## Data Formats

### Input Requirements
- **FITS Images**: Standard astronomical FITS format with proper header information
- **Star Catalogs**: GAIA format with RA, DEC, proper motions, and magnitudes
- **Ephemerides**: IMCCE or JPL format with time-series position data

### Output Products
- **Observation Files**: Standardized format with astrometric measurements
- **Statistical Reports**: Comprehensive analysis of observation quality
- **Summary Plots**: Visual representations of processing results
- **JSON Results**: Machine-readable processing metadata

## Processing Workflow

The complete processing workflow follows these sequential steps:

1. **Image Loading**: FITS files are loaded with header parsing for observation metadata
2. **Calibration**: Images undergo bias/dark/flat correction or algorithmic background subtraction
3. **Detection**: Stars are identified using threshold-based connected component analysis
4. **Catalog Preparation**: Reference catalogs are loaded and corrected for proper motion
5. **Pattern Matching**: Detected stars are matched against catalog positions
6. **Astrometric Solution**: Plate constants are determined through least-squares fitting
7. **Target Positioning**: Object positions are computed using ephemeris predictions
8. **Quality Assessment**: Results undergo statistical validation and outlier detection
9. **Report Generation**: Comprehensive reports and visualizations are produced

## Quality Control

The system implements multiple quality control mechanisms:

### Automated Validation
- Minimum number of matched stars required for reliable solutions
- Maximum residual thresholds for outlier detection
- Signal-to-noise ratio filtering for detection quality
- Iterative sigma-clipping for robust statistics

### Statistical Metrics
- Root-mean-square residuals in both coordinate directions
- Number of iterations required for convergence
- Percentage of outliers removed during processing
- Overall success rates for object detection and matching

## File Organization

```
ADIAS/
├── main.py                 # Main processing script
├── config.py              # Configuration management
├── utils/                 # Core processing modules
│   ├── preprocess.py      # Image calibration
│   ├── detect.py          # Stellar detection
│   ├── match.py           # Catalog matching
│   └── comoc.py           # Statistical analysis
├── output/                # Processing results
├── demo_output/           # Example outputs
└── __pycache__/           # Python cache files
```

## Performance Considerations

The system is designed for efficient processing of large image datasets:

### Processing Times
- Typical processing time: 10-30 seconds per image depending on complexity
- Scalable to hundreds of images per night
- Parallel processing capabilities for multiple objects

### Memory Requirements
- Moderate memory usage suitable for standard desktop computers
- Efficient array operations using NumPy optimizations
- Configurable processing parameters for different hardware capabilities

## Troubleshooting

### Common Issues and Solutions

**Insufficient Matched Stars**: Verify that star catalogs cover the observation field and check magnitude limits in configuration.

**High Residuals**: Review telescope parameters (focal length, pixel scale) and consider increasing the plate model complexity.

**Detection Problems**: Adjust background threshold parameters or examine image quality for focus and tracking issues.

**Coordinate System Errors**: Ensure proper epoch handling and verify ephemeris file formats and time systems.

## Contributing

This project follows standard Python development practices:
- Code documentation using docstrings
- Type hints for function parameters and returns
- Modular design for easy extension and maintenance
- Comprehensive error handling and logging

## Future Enhancements

Potential areas for system expansion include:
- GPU acceleration for image processing operations
- Machine learning integration for improved stellar classification
- Real-time processing capabilities for automated observatories
- Extended support for additional catalog formats and coordinate systems
- Web-based interface for remote processing and monitoring

## Support and Documentation

For detailed technical information about specific algorithms and implementation details, refer to the inline documentation within each module. The code includes comprehensive comments explaining the mathematical foundations and processing steps involved in each operation.