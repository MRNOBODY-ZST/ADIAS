import argparse
from utils import preprocess, detect, match, comoc
import logging
import numpy as np
from config import config
from astropy.io import fits
from astropy.time import Time
import os
from matplotlib import pyplot as plt
import time
from pathlib import Path
from typing import Dict, List, Any
import json
from datetime import datetime

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def parse_fits_header(header: fits.Header) -> Dict[str, Any]:
    """
    解析FITS头文件获取观测信息
    
    Args:
        header: FITS文件头
        
    Returns:
        观测信息字典
    """
    obs_info = {}
    
    # 尝试不同的时间关键字
    time_keywords = ['DATE-OBS', 'DATE-STA', 'DATE', 'OBSDATE']
    for keyword in time_keywords:
        if keyword in header:
            obs_info['obs_time_str'] = header[keyword]
            break
    else:
        logger.warning("No observation time found in header")
        obs_info['obs_time_str'] = '2023-11-05T00:00:00'
    
    # 解析时间
    try:
        if 'T' in obs_info['obs_time_str']:
            obs_time = Time(obs_info['obs_time_str'])
        else:
            # 如果只有日期，假设是UTC午夜
            obs_time = Time(obs_info['obs_time_str'] + 'T00:00:00')
            
        obs_info['obs_time'] = obs_time
        obs_info['mjd'] = obs_time.mjd
        obs_info['year'] = obs_time.datetime.year
        obs_info['month'] = obs_time.datetime.month
        obs_info['day'] = obs_time.datetime.day
        obs_info['hour'] = obs_time.datetime.hour
        obs_info['minute'] = obs_time.datetime.minute
        obs_info['second'] = obs_time.datetime.second
        
        # 计算天的小数（用于历表内插）
        obs_info['day_fraction'] = (obs_info['day'] + 
                                   (obs_info['hour'] * 3600 + obs_info['minute'] * 60 + obs_info['second']) / 86400.0)
        
    except Exception as e:
        logger.error(f"Error parsing observation time: {e}")
        obs_info['day_fraction'] = 1.0
    
    # 获取曝光时间
    exp_keywords = ['EXPTIME', 'EXPOSURE', 'ITIME']
    for keyword in exp_keywords:
        if keyword in header:
            obs_info['exposure_time'] = header[keyword]
            break
    else:
        obs_info['exposure_time'] = 1.0
        
    # 获取其他信息
    obs_info['object'] = header.get('OBJECT', 'Unknown')
    obs_info['telescope'] = header.get('TELESCOP', config.match.telescope_label)
    obs_info['instrument'] = header.get('INSTRUME', 'Unknown')
    
    return obs_info


def fits_open(file_name: str):
    """安全打开FITS文件"""
    try:
        fits_file = fits.open(file_name)
        return fits_file
    except FileNotFoundError:
        logger.warning(f"File not found: {file_name}")
        return None
    except Exception as e:
        logger.error(f"Error opening FITS file {file_name}: {e}")
        return None


def load_calibration_data():
    """加载标定数据"""
    dark_file = fits_open(config.dark_file)
    flat_file = fits_open(config.flat_file) 
    bias_file = fits_open(config.bias_file)
    
    calibration_data = {
        'dark': dark_file[0].data if dark_file else None,
        'flat': flat_file[0].data if flat_file else None,
        'bias': bias_file[0].data if bias_file else None
    }
    
    # 关闭文件
    if dark_file: dark_file.close()
    if flat_file: flat_file.close() 
    if bias_file: bias_file.close()
    
    return calibration_data


def load_reference_data():
    """加载参考数据（星表和历表）"""
    logger.info("Loading reference catalogs and ephemerides...")
    
    catalogs = []
    ephemerides = []
    
    # 加载多个目标的星表和历表
    for i in range(config.match.obj_total):
        # 星表文件
        if i < len(config.match.gaia_files):
            gaia_file = config.match.gaia_files[i]
        else:
            gaia_file = config.match.gaia_files[0]  # 使用第一个作为默认
            
        catalog = match.CatalogManager.load_gaia_catalog(
            gaia_file, 
            config.match.gaia_min_mag,
            config.match.gaia_max_mag
        )
        catalogs.append(catalog)
        
        # 历表文件
        if i < len(config.match.eph_files):
            eph_file = config.match.eph_files[i]
            ephemeris = match.EphemerisManager.load_ephemeris(
                eph_file, 
                config.match.object_eph_source
            )
            ephemerides.append(ephemeris)
        else:
            # 如果没有对应的历表文件，创建空历表
            ephemerides.append({
                "time": np.array([]),
                "ra": np.array([]),
                "dec": np.array([])
            })
    
    logger.info(f"Loaded {len(catalogs)} catalogs and {len(ephemerides)} ephemerides")
    return catalogs, ephemerides


def save_results(results: List[Dict[str, Any]], output_dir: str, image_name: str):
    """保存匹配结果"""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # 保存为JSON格式
    json_file = output_path / f"{Path(image_name).stem}_results.json"
    
    # 转换numpy数组为列表以便JSON序列化
    def convert_numpy(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.integer):
            return int(obj)
        return obj
    
    serializable_results = []
    for result in results:
        serializable_result = {}
        for key, value in result.items():
            if isinstance(value, dict):
                serializable_result[key] = {k: convert_numpy(v) for k, v in value.items()}
            else:
                serializable_result[key] = convert_numpy(value)
        serializable_results.append(serializable_result)
    
    with open(json_file, 'w') as f:
        json.dump(serializable_results, f, indent=2)
    
    logger.info(f"Results saved to {json_file}")
    
    # 保存观测数据格式（类似原始Fortran程序的object.out）
    obs_file = output_path / f"{Path(image_name).stem}_observations.txt"
    
    with open(obs_file, 'w') as f:
        f.write("# Year Month Day_fraction RA_obs DEC_obs RA_eph DEC_eph Res_RA Res_DEC Sigma N_matched Angle Image\n")
        
        for i, result in enumerate(results):
            if result['success']:
                # 这里需要从结果中提取观测时间信息
                f.write(f"{2023} {11} {5.0:.6f} {result['target_ra_obs']:.8f} {result['target_dec_obs']:.8f} "
                       f"{result['target_ra_eph']:.8f} {result['target_dec_eph']:.8f} "
                       f"{result['residual_ra']:.4f} {result['residual_dec']:.4f} "
                       f"{result['sigma']:.3f} {result['n_matched']} {result['plate_angle']:.2f} "
                       f"{Path(image_name).name}\n")
    
    logger.info(f"Observations saved to {obs_file}")


def create_summary_plots(image_data: np.ndarray, calibrated_image: np.ndarray,
                        detected_stars: tuple, match_results: List[Dict[str, Any]],
                        output_dir: str, image_name: str):
    """创建汇总图表"""
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle(f"Processing Results: {Path(image_name).name}", fontsize=16)
    
    # 原始图像
    im1 = axes[0, 0].imshow(image_data, cmap='gray', origin='lower',
                           vmin=np.percentile(image_data, 1),
                           vmax=np.percentile(image_data, 99))
    axes[0, 0].set_title('Original Image')
    axes[0, 0].set_xlabel('X (pixels)')
    axes[0, 0].set_ylabel('Y (pixels)')
    
    # 标定后图像
    im2 = axes[0, 1].imshow(calibrated_image, cmap='gray', origin='lower',
                           vmin=np.percentile(calibrated_image, 1),
                           vmax=np.percentile(calibrated_image, 99))
    axes[0, 1].set_title('Calibrated Image')
    axes[0, 1].set_xlabel('X (pixels)')
    axes[0, 1].set_ylabel('Y (pixels)')
    
    # 检测到的星象 - 正确解包7个元素
    x_det, y_det, flux_det, snr_det, star_ids, npix_det, overflow_flags = detected_stars
    
    if len(x_det) > 0:
        scatter = axes[1, 0].scatter(x_det, y_det, c=np.log10(flux_det + 1), 
                                    s=30, cmap='viridis', alpha=0.7)
        axes[1, 0].set_title(f'Detected Stars ({len(x_det)} stars)')
        axes[1, 0].set_xlabel('X (pixels)')
        axes[1, 0].set_ylabel('Y (pixels)')
        plt.colorbar(scatter, ax=axes[1, 0], label='log10(Flux)')
    else:
        axes[1, 0].set_title('No Stars Detected')
        axes[1, 0].set_xlabel('X (pixels)')
        axes[1, 0].set_ylabel('Y (pixels)')
    
    # 匹配结果统计
    n_successful = sum(1 for r in match_results if r['success'])
    n_matched_total = sum(r['n_matched'] for r in match_results if r['success'])
    avg_sigma = np.mean([r['sigma'] for r in match_results if r['success']]) if n_successful > 0 else 0
    
    # 创建匹配统计图
    axes[1, 1].axis('off')
    stats_text = f"""Matching Statistics:
    
Objects processed: {len(match_results)}
Successful matches: {n_successful}
Total matched stars: {n_matched_total}
Average sigma: {avg_sigma:.3f}"
Total detected stars: {len(x_det)}
Overflow stars: {np.sum(overflow_flags) if len(overflow_flags) > 0 else 0}

Successful Objects:"""
    
    for i, result in enumerate(match_results):
        if result['success']:
            stats_text += f"\nObj {i+1}: {result['n_matched']} stars, σ={result['sigma']:.3f}\""
            stats_text += f"\n   RA residual: {result['residual_ra']:.2f}\""
            stats_text += f"\n   DEC residual: {result['residual_dec']:.2f}\""
    
    axes[1, 1].text(0.05, 0.95, stats_text, transform=axes[1, 1].transAxes,
                   fontsize=10, verticalalignment='top', fontfamily='monospace')
    
    plt.tight_layout()
    
    # 保存图像
    plot_file = Path(output_dir) / f"{Path(image_name).stem}_summary.png"
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Summary plot saved to {plot_file}")


def process_single_image(fits_file_path: str, calibration_data: Dict[str, Any],
                        catalogs: List[Dict[str, np.ndarray]], 
                        ephemerides: List[Dict[str, np.ndarray]],
                        output_dir: str) -> Dict[str, Any]:
    """处理单个FITS图像"""
    
    logger.info(f"Processing: {fits_file_path}")
    
    # 打开FITS文件
    fits_file = fits.open(fits_file_path)
    image_data = fits_file[0].data.astype(np.float64)
    header = fits_file[0].header
    
    # 解析头文件
    obs_info = parse_fits_header(header)
    
    # 预处理
    start_time = time.time()
    if config.preprocess.super == 0:
        calibrated_image = preprocess.standard_calibration(
            image_data,
            calibration_data['flat'],
            calibration_data['bias'], 
            calibration_data['dark']
        )
    else:
        calibrated_image = preprocess.algorithm_calibration(
            image_data,
            config.preprocess.med_length,
            config.preprocess.med_width,
            config.preprocess.background
        )
    
    preprocess_time = time.time() - start_time
    
    # 星象检测
    start_time = time.time()
    detected_stars = detect.stellar_locate_enhanced(
        calibrated_image,
        config.detect.background_threshold,
        config.detect.position_method
    )
    
    detect_time = time.time() - start_time
    
    # detected_stars包含7个元素：x, y, flux, snr, star_ids, npix, overflow_flags
    x_det, y_det, flux_det, snr_det, star_ids, npix_det, overflow_flags = detected_stars
    
    logger.info(f"Detected {len(x_det)} stars")
    
    # 星表匹配
    start_time = time.time()
    
    telescope_params = {
        "focal_length": config.match.telescope_focal,
        "pixel_scale": config.match.ccd_scale * 1000,  # 转换为毫角秒/像素
        "field_size": config.match.ccd_field_size
    }
    
    match_config = {
        "model_type": config.match.model_type,
        "plate_angle": config.match.plate_angle,
        "match_limit": config.match.match_limit,
        "field_size": config.match.ccd_field_size
    }
    
    match_results = match.match_multiple_objects(
        detected_stars, catalogs, ephemerides,
        obs_info['day_fraction'], telescope_params, match_config
    )
    
    match_time = time.time() - start_time
    
    # 保存结果
    save_results(match_results, output_dir, fits_file_path)
    
    # 创建汇总图表
    create_summary_plots(image_data, calibrated_image, detected_stars,
                        match_results, output_dir, fits_file_path)
    
    # 关闭FITS文件
    fits_file.close()
    
    # 计算统计信息
    n_successful = sum(1 for r in match_results if r['success'])
    
    processing_stats = {
        'image_name': Path(fits_file_path).name,
        'preprocess_time': preprocess_time,
        'detect_time': detect_time,
        'match_time': match_time,
        'total_time': preprocess_time + detect_time + match_time,
        'n_detected_stars': len(x_det),
        'n_successful_matches': n_successful,
        'n_total_objects': len(match_results),
        'obs_info': obs_info,
        'match_results': match_results
    }
    
    logger.info(f"Processing completed in {processing_stats['total_time']:.2f}s")
    logger.info(f"  - Preprocessing: {preprocess_time:.2f}s")
    logger.info(f"  - Detection: {detect_time:.2f}s") 
    logger.info(f"  - Matching: {match_time:.2f}s")
    logger.info(f"Successfully matched {n_successful}/{len(match_results)} objects")
    
    return processing_stats


def main():
    """主程序"""
    parser = argparse.ArgumentParser(description='CCD图像天体测量处理流程')
    parser.add_argument('--config', type=str, help='配置文件路径')
    parser.add_argument('--output', type=str, default='output', help='输出目录')
    parser.add_argument('--debug', action='store_true', help='启用调试模式')
    parser.add_argument('--skip-analysis', action='store_true', help='跳过统计分析步骤')
    
    args = parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # 加载标定数据
    logger.info("Loading calibration data...")
    calibration_data = load_calibration_data()
    
    # 加载参考数据
    catalogs, ephemerides = load_reference_data()
    
    # 创建输出目录
    output_dir = args.output
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # 处理所有目录中的图像
    all_stats = []
    all_object_files = []  # 收集所有目标文件用于后续统计分析
    
    for fits_folder in config.fits_directory:
        image_dir = os.path.join(fits_folder, "fits")
        
        if not os.path.exists(image_dir):
            logger.warning(f"Image directory does not exist: {image_dir}")
            continue
            
        fits_files = [f for f in os.listdir(image_dir) if f.lower().endswith(('.fits', '.fit'))]
        fits_files.sort()
        
        logger.info(f"Found {len(fits_files)} FITS files in {image_dir}")
        
        # 处理当前目录的所有图像
        day_stats = []
        day_object_files = []
        
        for fits_file in fits_files:
            fits_file_path = os.path.join(image_dir, fits_file)
            
            try:
                stats = process_single_image(
                    fits_file_path, calibration_data, 
                    catalogs, ephemerides, output_dir
                )
                day_stats.append(stats)
                all_stats.append(stats)
                
                # 收集观测数据文件路径
                for obj_id in range(1, config.match.obj_total + 1):
                    obj_file = os.path.join(os.path.dirname(fits_file_path), f"object_{obj_id}.out")
                    if os.path.exists(obj_file):
                        day_object_files.append(obj_file)
                        all_object_files.append(obj_file)
                
            except Exception as e:
                logger.error(f"Error processing {fits_file}: {e}")
                continue
        
        # 为当前目录生成日度汇总
        if day_stats:
            generate_daily_summary(day_stats, fits_folder, output_dir)
    
    # 生成总结报告
    generate_summary_report(all_stats, output_dir)
    
    # 执行统计分析和野值剔除
    if not args.skip_analysis and all_object_files:
        logger.info("Starting observation data analysis...")
        perform_observation_analysis(all_object_files, output_dir)
    
    logger.info("Processing complete!")


def perform_observation_analysis(object_files: List[str], output_dir: str):
    """
    执行观测数据统计分析和野值剔除
    
    Args:
        object_files: 观测数据文件列表
        output_dir: 输出目录
    """
    try:
        # 按目标ID分组文件
        object_groups = {}
        for file_path in object_files:
            # 从文件名提取目标ID
            filename = Path(file_path).name
            if 'object_' in filename:
                try:
                    obj_id = int(filename.split('object_')[1].split('.')[0])
                    if obj_id not in object_groups:
                        object_groups[obj_id] = []
                    object_groups[obj_id].append(file_path)
                except (ValueError, IndexError):
                    logger.warning(f"Could not extract object ID from {filename}")
                    continue
        
        if not object_groups:
            logger.warning("No valid object files found for analysis")
            return
        
        logger.info(f"Found observation data for {len(object_groups)} objects")
        
        # 为每个目标创建合并的观测文件
        merged_files = []
        object_ids = []
        
        for obj_id, files in object_groups.items():
            merged_file = Path(output_dir) / f"merged_object_{obj_id}.out"
            merge_observation_files(files, merged_file)
            merged_files.append(str(merged_file))
            object_ids.append(obj_id)
        
        # 配置统计分析参数
        analysis_config = {
            'eps': config.observation_calculation.eps,
            'standard_limit': config.observation_calculation.standard_limit,
            'mean_limit': config.observation_calculation.mean_limit,
            'delete_flag': config.observation_calculation.delete_flag,
            'specified_output': output_dir
        }
        
        # 执行统计分析
        results = comoc.process_observation_results(
            analysis_config, merged_files, object_ids
        )
        
        logger.info("Observation data analysis completed")
        logger.info(f"Analysis results saved to {output_dir}")
        
        # 清理临时合并文件
        for merged_file in merged_files:
            try:
                os.remove(merged_file)
            except OSError:
                pass
        
        return results
        
    except Exception as e:
        logger.error(f"Error during observation analysis: {e}")
        return None


def merge_observation_files(file_paths: List[str], output_file: Path):
    """
    合并多个观测数据文件
    
    Args:
        file_paths: 观测文件路径列表
        output_file: 输出文件路径
    """
    all_observations = []
    
    for file_path in file_paths:
        try:
            if os.path.exists(file_path):
                with open(file_path, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if line and not line.startswith('#'):
                            all_observations.append(line)
        except Exception as e:
            logger.warning(f"Error reading {file_path}: {e}")
            continue
    
    # 按时间排序观测数据
    def extract_time(line):
        try:
            parts = line.split()
            if len(parts) >= 3:
                year = int(parts[0])
                month = int(parts[1])
                day_fraction = float(parts[2])
                return year * 10000 + month * 100 + day_fraction
        except:
            pass
        return 0
    
    all_observations.sort(key=extract_time)
    
    # 写入合并文件
    with open(output_file, 'w') as f:
        f.write("# Merged observation data\n")
        f.write("# Year Month Day_fraction RA_obs DEC_obs RA_eph DEC_eph Res_RA Res_DEC Sigma N_matched Angle Image\n")
        for obs in all_observations:
            f.write(f"{obs}\n")
    
    logger.debug(f"Merged {len(all_observations)} observations from {len(file_paths)} files")


def generate_daily_summary(day_stats: List[Dict[str, Any]], 
                          fits_folder: str, output_dir: str):
    """
    生成日度汇总报告
    
    Args:
        day_stats: 当日处理统计列表
        fits_folder: FITS文件目录
        output_dir: 输出目录
    """
    date_str = Path(fits_folder).name
    summary_file = Path(output_dir) / f"daily_summary_{date_str}.txt"
    
    with open(summary_file, 'w') as f:
        f.write(f"Daily Processing Summary - {date_str}\n")
        f.write("=" * 50 + "\n\n")
        
        total_images = len(day_stats)
        total_time = sum(s['total_time'] for s in day_stats)
        total_detected = sum(s['n_detected_stars'] for s in day_stats)
        total_successful = sum(s['n_successful_matches'] for s in day_stats)
        
        f.write(f"Images processed: {total_images}\n")
        f.write(f"Total processing time: {total_time:.2f} seconds\n")
        f.write(f"Average time per image: {total_time/total_images:.2f} seconds\n") 
        f.write(f"Stars detected: {total_detected}\n")
        f.write(f"Successful matches: {total_successful}\n")
        f.write(f"Success rate: {total_successful/(total_images*config.match.obj_total)*100:.1f}%\n\n")
        
        f.write("Image Details:\n")
        f.write("-" * 20 + "\n")
        
        for stats in day_stats:
            f.write(f"\n{stats['image_name']}:\n")
            f.write(f"  Processing time: {stats['total_time']:.2f}s\n")
            f.write(f"  Stars detected: {stats['n_detected_stars']}\n")
            f.write(f"  Successful matches: {stats['n_successful_matches']}/{stats['n_total_objects']}\n")
    
    logger.info(f"Daily summary saved to {summary_file}")


def generate_summary_report(all_stats: List[Dict[str, Any]], output_dir: str):
    """生成总结报告"""
    
    if not all_stats:
        logger.warning("No processing statistics to summarize")
        return
    
    report_file = Path(output_dir) / "processing_summary.txt"
    
    with open(report_file, 'w') as f:
        f.write("CCD Image Astrometry Processing Summary\n")
        f.write("=" * 50 + "\n\n")
        
        total_images = len(all_stats)
        total_time = sum(s['total_time'] for s in all_stats)
        total_detected = sum(s['n_detected_stars'] for s in all_stats)
        total_successful = sum(s['n_successful_matches'] for s in all_stats)
        
        f.write(f"Total images processed: {total_images}\n")
        f.write(f"Total processing time: {total_time:.2f} seconds\n")
        f.write(f"Average time per image: {total_time/total_images:.2f} seconds\n")
        f.write(f"Total stars detected: {total_detected}\n")
        f.write(f"Average stars per image: {total_detected/total_images:.1f}\n")
        f.write(f"Total successful matches: {total_successful}\n")
        f.write(f"Success rate: {total_successful/(total_images*config.match.obj_total)*100:.1f}%\n\n")
        
        f.write("Individual Image Results:\n")
        f.write("-" * 30 + "\n")
        
        for stats in all_stats:
            f.write(f"\nImage: {stats['image_name']}\n")
            f.write(f"  Processing time: {stats['total_time']:.2f}s\n")
            f.write(f"  Stars detected: {stats['n_detected_stars']}\n")
            f.write(f"  Successful matches: {stats['n_successful_matches']}/{stats['n_total_objects']}\n")
            
            for result in stats['match_results']:
                if result['success']:
                    f.write(f"    Object {result['object_id']}: {result['n_matched']} stars, "
                           f"σ={result['sigma']:.3f}\", "
                           f"residuals=({result['residual_ra']:.2f}\", {result['residual_dec']:.2f}\")\n")
    
    logger.info(f"Summary report saved to {report_file}")


if __name__ == "__main__":
    main()