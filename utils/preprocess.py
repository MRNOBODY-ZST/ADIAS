import os
from typing import Dict, List, Tuple, Optional
import numpy as np
from astropy.io import fits
import configparser
from pathlib import Path

def load_config(config_path:str):
    config = configparser.ConfigParser()
    config.read(config_path)
    sections = config.sections()
    return sections