import argparse
from utils import preprocess, detect
import logging
import numpy as np
from config import config
from astropy.io import fits
import os
from matplotlib import pyplot as plt
import time


def fits_open(file_name: str):
    try:
        fits_file = fits.open(file_name)
    except FileNotFoundError:
        fits_file = None
    return fits_file


if __name__ == "__main__":
    for fits_folder in config.fits_directory:
        image_dir = os.path.join(fits_folder, "fits")
        fits_file_dirs = os.listdir(image_dir)
        fits_file_list = []
        fits_header_list = []
        dark_file = fits_open(config.dark_file)
        flat_file = fits_open(config.flat_file)
        bias_file = fits_open(config.bias_file)
        for fits_file_dir in fits_file_dirs:
            print(f"Processing FITS file: {os.path.join(image_dir,fits_file_dir)}")
            fits_file = fits.open(os.path.join(image_dir, fits_file_dir))

            # fits_file_list.append(np.array(fits_file[0].data))
            # fits_header_list.append(fits_file.info())
            # fits_file.close()
            if config.preprocess.super == 0:
                calibrated_image = preprocess.standard_calibration(
                    fits_file[0].data,
                    dark_file,
                    flat_file,
                    bias_file,
                )
            else:
                calibrated_image = preprocess.algorithm_calibration(
                    fits_file[0].data,
                    config.preprocess.med_length,
                    config.preprocess.med_width,
                    config.preprocess.background,
                )
            
            # star_list = detect.stellar_locate(fits_file[0].data, 5)

            plt.figure()
            plt.subplot(1, 2, 1)
            plt.imshow(
                fits_file[0].data,
                cmap="gray",
                vmin=fits_file[0].data.min(),
                vmax=fits_file[0].data.max(),
                origin="lower",
            )
            plt.subplot(1, 2, 2)
            plt.imshow(
                calibrated_image,
                cmap="gray",
                vmin=calibrated_image.min(),
                vmax=calibrated_image.max(),
                origin="lower",
            )
            plt.show()
