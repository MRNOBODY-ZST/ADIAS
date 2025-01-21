from dataclasses import dataclass, field
from typing import List


@dataclass
class PreprocessConfig:
    bias: int = 0
    dark: int = 0
    flat: int = 0
    super: int = 1
    med_length: int = 35
    med_width: int = 35
    background: int = 2
    enhance: int = 1


@dataclass
class DetectConfig:
    background_threshold: int = 5
    signal_noise_threshold: float = 1.0
    position_method: int = 1


@dataclass
class MatchConfig:
    telescope_label: str = "km100B"
    telescope_focal: float = 13300.0
    ccd_scale: float = 0.0135
    ccd_field_size: int = 15
    model_type: int = 6
    plate_angle: int = 181
    object_eph_source: str = "IMCCE"
    gaia_min_mag: int = 5
    gaia_max_mag: int = 17
    match_limit: float = 5.0
    delta_t: float = 0.0
    obj_total: int = 8
    eph_files: List[str] = field(default_factory=lambda: [
        "EPH_S1_202311.DAT",
        "EPH_S2_202311.DAT",
        "EPH_S3_202311.DAT",
        "EPH_S4_202311.DAT",
        "EPH_S5_202311.DAT",
        "EPH_S6_202311.DAT",
        "EPH_S7_202311.DAT",
        "EPH_S8_202311.DAT"
    ])
    gaia_files: List[str] = field(default_factory=lambda: ["GAIA3_S_202311.DAT"] * 7)


@dataclass
class ObservationCalculationConfig:
    specified_output: str = "data"
    eps: float = 0.5
    standard_limit: int = 3
    mean_limit: int = 0
    delete_flag: int = 0



@dataclass
class Config:
    preprocess: PreprocessConfig = field(default_factory=PreprocessConfig)
    detect: DetectConfig = field(default_factory=DetectConfig)
    match: MatchConfig = field(default_factory=MatchConfig)
    observation_calculation: ObservationCalculationConfig = field(default_factory=ObservationCalculationConfig)
    fits_directory = ["data/img/2023/11/05"]
    dark_file = "data/img/2023/11/05/dark.fits"
    flat_file = "data/img/2023/11/05/flat.fits"
    bias_file = "data/img/2023/11/05/bias.fits"
    


# Create a global config instance
config = Config()