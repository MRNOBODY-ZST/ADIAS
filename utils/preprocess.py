from typing import Optional, overload, Any
import numpy as np
from numpy import floating
from scipy import signal, stats


def enhance(fits_image: np.ndarray): ...


def calculate_background(
    data: np.ndarray, sigma_factor: float = 2.6
) -> tuple[floating[Any], floating[Any]]:
    flat_data = data.ravel()
    avg_value = np.mean(flat_data)
    sigma = np.std(flat_data, ddof=1)  # ddof=1 for sample standard deviation

    while True:
        mask = np.abs(flat_data - avg_value) <= sigma_factor * sigma
        new_avg = np.mean(flat_data[mask])
        new_sigma = np.std(flat_data[mask], ddof=1)
        if new_sigma <= 1e-6:
            return new_avg, new_sigma
        if abs(new_sigma - sigma) < 0.01 * new_sigma:
            return new_avg, new_sigma
        avg_value = new_avg
        sigma = new_sigma


def calculate_background_robust(
    data: np.ndarray, sigma_factor: float = 2.6
) -> tuple[floating[Any], float]:
    flat_data = data.ravel()
    median = np.median(flat_data)
    mad = stats.median_abs_deviation(flat_data)
    sigma = mad * 1.4826
    while True:
        mask = np.abs(flat_data - median) <= sigma_factor * sigma
        new_median = np.median(flat_data[mask])
        new_mad = stats.median_abs_deviation(flat_data[mask])
        new_sigma = new_mad * 1.4826
        if new_sigma <= 1e-6:
            return new_median, new_sigma
        if abs(new_sigma - sigma) < 0.01 * new_sigma:
            return new_median, new_sigma
        median = new_median
        sigma = new_sigma


def algorithm_calibration(
    fits_image: np.ndarray,
    med_length: int,
    med_width: int,
    background_mode: int = 1,
) -> tuple[np.ndarray, np.ndarray, floating[Any], floating[Any]]:
    """
    Process FITS image background using median filtering.

    Args:
        fits_image: Input FITS image data
        med_length: Length of median filter window
        med_width: Width of median filter window
        background_mode: Background subtraction mode (1: division, 2: subtraction)

    Returns:
        Tuple containing:
        - processed_image: Background-corrected image
        - background_map: Estimated background map
        - final_background: Final background value
        - final_sigma: Final background sigma
    """
    # 确保输入是numpy数组
    fits_image = np.asarray(fits_image)

    # 计算原始图像的背景和sigma
    background0, background_sigma = calculate_background(fits_image)
    print(
        f"Before processing - Background: {background0:.3f}, Sigma: {background_sigma:.3f}"
    )

    # 创建输出数组
    processed_image = fits_image.copy()
    background_map = fits_image.copy()

    # 应用中值滤波
    if med_length != med_width:
        # 如果长宽不同，进行两次滤波
        print(f"First median filter window (length, width): {med_width}, {med_length}")
        temp = signal.medfilt2d(
            fits_image.astype(float), kernel_size=(med_width, med_length)
        )

        print(f"Second median filter window (length, width): {med_length}, {med_width}")
        background_map = signal.medfilt2d(temp, kernel_size=(med_length, med_width))
    else:
        # 长宽相同，一次滤波
        print(f"Median filter window (length, width): {med_length}, {med_width}")
        background_map = signal.medfilt2d(
            fits_image.astype(float), kernel_size=(med_length, med_width)
        )

    # 根据background_mode处理图像
    mask = np.abs(fits_image) >= 2  # 只处理绝对值>=2的像素

    if background_mode == 1:
        # 使用除法模式
        processed_image = np.where(
            mask, fits_image / background_map * background0, background_map
        )
    elif background_mode == 2:
        # 使用减法模式
        processed_image = np.where(
            mask, fits_image - background_map + background0, background_map
        )
    else:
        raise ValueError("background_mode must be 1 or 2")

    # 计算处理后的背景和sigma
    final_background, final_sigma = calculate_background(processed_image)
    print(
        f"After processing - Background: {final_background:.3f}, Sigma: {final_sigma:.3f}"
    )

    return processed_image


def standard_calibration(
    fits_image: np.ndarray,
    flat_image: Optional[np.ndarray] = None,
    bias_image: Optional[np.ndarray] = None,
    dark_image: Optional[np.ndarray] = None,
) -> np.ndarray:
    if all(x is None for x in (flat_image, bias_image, dark_image)):
        return fits_image / np.mean(fits_image)
    if not flat_image:
        raise ValueError("For full calibration, flat field must be provided")
    return (
        fits_image - bias_image
        if bias_image
        else (
            np.zeros(fits_image.size) - dark_image
            if dark_image
            else np.zeros(fits_image.size)
        )
    ) / flat_image
