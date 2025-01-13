import argparse
from utils import preprocess

if __name__ == '__main__':
    fits_path_file = "data/fitspath.ini"
    config_path = "data/adias.cfg"

    config = preprocess.load_config(config_path)
    print(config)
