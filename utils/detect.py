import numpy as np
from scipy.ndimage import label
from typing import Tuple, Optional
from utils.preprocess import calculate_background

def stellar_locate(fits_image: np.ndarray, background_threshold: int = 5):
    """
    Detect and locate stars in astronomical images.

    Args:
        fits_image: Input image as numpy array
        background_threshold: Threshold factor for background detection (default: 5)

    Returns:
        tuple containing:
        - x_centers: Array of x coordinates for detected stars
        - y_centers: Array of y coordinates for detected stars
        - fluxes: Array of star fluxes
        - signal_noise_ratio_list: Array of signal-to-noise ratios
        - num_pixels: Array of pixel counts per star
    """
    # 计算背景和噪声
    background, background_sigma = calculate_background(fits_image)

    # 创建二值化图像用于标记
    threshold = background + background_threshold * background_sigma
    binary_image = fits_image > threshold

    # 使用8邻域连通性进行标记
    structure = np.ones((3, 3), dtype=np.int32)  # 8邻域
    labeled_image, num_features = label(binary_image, structure=structure)

    # 初始化结果数组
    x_centers = []
    y_centers = []
    fluxes = []
    signal_noise_ratio_list = []
    num_pixels = []

    # 处理每个检测到的目标
    for i in range(1, num_features + 1):
        # 获取当前目标的掩码
        mask = labeled_image == i

        # 跳过边缘目标
        if np.any(mask[0, :]) or np.any(mask[-1, :]) or np.any(mask[:, 0]) or np.any(mask[:, -1]):
            continue

        # 获取目标区域的坐标和像素值
        y_coords, x_coords = np.where(mask)
        pixel_values = fits_image[mask]

        # 计算总流量
        flux = np.sum(pixel_values - background)

        # 检查目标大小
        if len(pixel_values) < 3 or len(pixel_values) > 700:  # 对应原代码的minpix和maxpix
            continue

        # 计算质心
        x_center = np.sum(x_coords * pixel_values) / np.sum(pixel_values)
        y_center = np.sum(y_coords * pixel_values) / np.sum(pixel_values)

        # 计算信噪比
        noise = np.sqrt(flux + len(pixel_values) * background_sigma ** 2)
        snr = flux / noise if noise > 0 else 0

        # 存储结果
        x_centers.append(x_center)
        y_centers.append(y_center)
        print("Find star at: ", x_center, y_center)
        fluxes.append(flux)
        signal_noise_ratio_list.append(snr)
        num_pixels.append(len(pixel_values))

    # 转换为numpy数组
    results = map(np.array, [x_centers, y_centers, fluxes, signal_noise_ratio_list, num_pixels])
    arrays = list(results)

    # 按流量排序
    if len(arrays[0]) > 0:
        sort_idx = np.argsort(arrays[2])[::-1]  # 按流量降序排序
        arrays = [arr[sort_idx] for arr in arrays]

    return tuple(arrays)


def stellar_locate_enhanced(
        fits_image: np.ndarray,
        background_threshold: int = 5,
        pos_method: int = 2,
        flat_flag: bool = False,
        super_flag: bool = False
) -> Tuple[np.ndarray, ...]:
    """
    Enhanced stellar detection with equivalence table processing.

    Args:
        fits_image: Input image as numpy array
        background_threshold: Threshold factor for background detection
        pos_method: Power factor for weighted centroid calculation
        flat_flag: Flag for flat field correction
        super_flag: Flag for super star processing

    Returns:
        tuple containing:
        - x_centers: Array of x coordinates
        - y_centers: Array of y coordinates
        - fluxes: Array of star fluxes
        - signal_noise_ratio_list: Array of SNRs
        - star_ids: Array of star IDs
        - num_pixels: Array of pixel counts per star
        - overflow_flags: Array of overflow flags
    """
    # 计算背景和噪声
    background, background_sigma = calculate_background(fits_image)

    # 创建背景减除后的图像
    thresh_image = fits_image - (background + background_threshold * background_sigma)
    thresh_image[thresh_image < 0] = 0

    # 创建二值化图像
    binary_image = thresh_image > 0

    # 初始标记
    height, width = fits_image.shape
    labeled_image = np.zeros_like(binary_image, dtype=np.int32)
    next_label = 1

    # 第一遍标记 - 8邻域连通域标记
    def get_neighbors(y: int, x: int) -> list:
        """获取像素的8邻域标记"""
        neighbors = []
        for dy, dx in [(-1, 1), (-1, 0), (-1, -1), (0, -1)]:  # 右上、正上、左上、左
            ny, nx = y + dy, x + dx
            if 0 <= ny < height and 0 <= nx < width:
                if labeled_image[ny, nx] > 0:
                    neighbors.append(labeled_image[ny, nx])
        return neighbors

    # 第一遍扫描
    equivalence_dict = {}  # 等价表
    for i in range(1, height - 1):
        for j in range(1, width - 1):
            if binary_image[i, j]:
                neighbors = get_neighbors(i, j)
                if not neighbors:  # 没有邻居
                    labeled_image[i, j] = next_label
                    next_label += 1
                else:
                    min_label = min(neighbors)
                    labeled_image[i, j] = min_label
                    # 更新等价表
                    for n in neighbors:
                        if n != min_label:
                            if n not in equivalence_dict:
                                equivalence_dict[n] = set()
                            equivalence_dict[n].add(min_label)

    # 解析等价表
    final_labels = {}
    for label in range(1, next_label):
        if label in equivalence_dict:
            min_equiv = min(equivalence_dict[label])
            final_labels[label] = min_equiv
        else:
            final_labels[label] = label

    # 第二遍扫描 - 统一标签
    for i in range(1, height - 1):
        for j in range(1, width - 1):
            if labeled_image[i, j] > 0:
                labeled_image[i, j] = final_labels[labeled_image[i, j]]

    # 初始化结果数组
    x_centers = []
    y_centers = []
    fluxes = []
    signal_noise_ratio_list = []
    star_ids = []
    num_pixels = []
    overflow_flags = []

    # 处理每个检测到的目标
    unique_labels = np.unique(labeled_image[labeled_image > 0])

    for label in unique_labels:
        mask = labeled_image == label

        # 跳过边缘10个像素的目标
        if (np.any(mask[:10, :]) or np.any(mask[-10:, :]) or
                np.any(mask[:, :10]) or np.any(mask[:, -10:])):
            continue

        y_coords, x_coords = np.where(mask)
        pixel_values = thresh_image[mask]
        original_values = fits_image[mask]

        # 大小检查
        if len(pixel_values) < 3 or len(pixel_values) > 700:
            continue

        # 检查是否溢出
        if np.issubdtype(fits_image.dtype, np.integer):
            # 对于整数类型，检查是否达到最大值
            overflow = np.any(original_values >= np.iinfo(fits_image.dtype).max)
        else:
            # 对于浮点数类型，使用一个合理的阈值来判断溢出
            # 通常CCD图像的饱和值在40000-65000之间
            saturation_threshold = 60000.0
            overflow = np.any(original_values >= saturation_threshold)

        # 使用指定方法计算质心
        weights = pixel_values ** pos_method
        x_center = np.sum(x_coords * weights) / np.sum(weights)
        y_center = np.sum(y_coords * weights) / np.sum(weights)

        # 计算总流量和信噪比
        flux = np.sum(original_values)
        noise = np.sqrt(flux + len(pixel_values) * background_sigma ** 2)
        snr = flux / noise if noise > 0 else 0

        if not overflow:  # 只保存未溢出的星象
            x_centers.append(x_center)
            y_centers.append(y_center)
            fluxes.append(flux)
            signal_noise_ratio_list.append(snr)
            star_ids.append(label)
            num_pixels.append(len(pixel_values))
            overflow_flags.append(overflow)

    # 转换为numpy数组并按流量排序
    results = map(np.array, [
        x_centers, y_centers, fluxes, signal_noise_ratio_list,
        star_ids, num_pixels, overflow_flags
    ])
    arrays = list(results)

    if len(arrays[0]) > 0:
        sort_idx = np.argsort(arrays[2])[::-1]
        arrays = [arr[sort_idx] for arr in arrays]

    return tuple(arrays)

