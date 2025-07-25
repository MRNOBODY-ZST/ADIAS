import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, Any
from pathlib import Path
import json
import logging
from datetime import datetime
import shutil
import os

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class ObservationStatistics:
    """观测数据统计分析类"""
    
    def __init__(self, config: Dict[str, Any]):
        """
        初始化统计分析器
        
        Args:
            config: 配置参数字典
        """
        self.config = config
        self.eps = config.get('eps', 0.5)  # 绝对残差限制
        self.std_limit = config.get('standard_limit', 2.6)  # 标准偏差倍数
        self.mean_limit = config.get('mean_limit', 0.03)  # 均值限制
        self.delete_flag = config.get('delete_flag', 0)  # 是否删除临时文件
        self.output_dir = config.get('specified_output', 'output')
        
        # 创建输出目录
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)
        
        # 统计信息
        self.statistics = {}
        
    def process_observation_data(self, object_files: List[str], 
                               object_ids: List[int]) -> Dict[str, Any]:
        """
        处理观测数据文件
        
        Args:
            object_files: object.out文件路径列表
            object_ids: 目标对象ID列表
            
        Returns:
            处理结果统计
        """
        results = {}
        
        for obj_id, obj_file in zip(object_ids, object_files):
            logger.info(f"Processing object {obj_id}: {obj_file}")
            
            # 读取观测数据
            obs_data = self._read_observation_file(obj_file)
            
            if obs_data is None or len(obs_data) == 0:
                logger.warning(f"No valid data found in {obj_file}")
                continue
                
            # 统计分析
            stats = self._analyze_residuals(obs_data, obj_id)
            
            # 剔除野值
            filtered_data, final_stats = self._remove_outliers(obs_data, stats)
            
            # 保存结果
            self._save_results(filtered_data, final_stats, obj_file, obj_id)
            
            results[obj_id] = {
                'original_stats': stats,
                'final_stats': final_stats,
                'original_count': len(obs_data),
                'final_count': len(filtered_data),
                'file_path': obj_file
            }
            
        # 生成汇总报告
        self._generate_summary_report(results)
        
        # 删除临时文件（如果设置）
        if self.delete_flag == 1:
            self._clean_temporary_files(object_files)
            
        return results
    
    def _read_observation_file(self, file_path: str) -> Optional[pd.DataFrame]:
        """
        读取观测数据文件
        
        Args:
            file_path: 文件路径
            
        Returns:
            观测数据DataFrame
        """
        try:
            if not os.path.exists(file_path):
                logger.warning(f"File not found: {file_path}")
                return None
                
            # 读取固定格式的观测数据
            data = []
            with open(file_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                        
                    try:
                        # 解析观测数据行
                        # 格式：Year Month Day_fraction RA_obs DEC_obs RA_eph DEC_eph Res_RA Res_DEC Sigma N_matched Angle Image
                        parts = line.split()
                        if len(parts) >= 12:
                            data.append({
                                'year': int(parts[0]),
                                'month': int(parts[1]), 
                                'day_fraction': float(parts[2]),
                                'ra_obs': float(parts[3]),
                                'dec_obs': float(parts[4]),
                                'ra_eph': float(parts[5]),
                                'dec_eph': float(parts[6]),
                                'res_ra': float(parts[7]),
                                'res_dec': float(parts[8]),
                                'sigma': float(parts[9]),
                                'n_matched': int(parts[10]),
                                'angle': float(parts[11]),
                                'image': parts[12] if len(parts) > 12 else '',
                                'original_line': line
                            })
                    except (ValueError, IndexError) as e:
                        logger.debug(f"Skipping invalid line: {line}")
                        continue
            
            if not data:
                return None
                
            df = pd.DataFrame(data)
            
            # 过滤掉明显异常的数据
            df = df[
                (np.abs(df['res_ra']) < 100.0) & 
                (np.abs(df['res_dec']) < 100.0) &
                (df['sigma'] > 0.001) &
                (df['sigma'] < 10.0)
            ]
            
            return df
            
        except Exception as e:
            logger.error(f"Error reading observation file {file_path}: {e}")
            return None
    
    def _analyze_residuals(self, data: pd.DataFrame, obj_id: int) -> Dict[str, Any]:
        """
        分析残差统计
        
        Args:
            data: 观测数据
            obj_id: 目标ID
            
        Returns:
            统计结果
        """
        n_obs = len(data)
        
        if n_obs < 1:
            return {
                'n_obs': 0,
                'mean_ra': 0.0,
                'mean_dec': 0.0,
                'std_ra': 0.0,
                'std_dec': 0.0
            }
        
        # 计算均值
        mean_ra = np.mean(data['res_ra'])
        mean_dec = np.mean(data['res_dec'])
        
        # 计算标准偏差
        if n_obs > 1:
            std_ra = np.std(data['res_ra'], ddof=1)
            std_dec = np.std(data['res_dec'], ddof=1)
        else:
            std_ra = 0.0
            std_dec = 0.0
        
        stats = {
            'n_obs': n_obs,
            'mean_ra': mean_ra,
            'mean_dec': mean_dec,
            'std_ra': std_ra,
            'std_dec': std_dec,
            'min_ra': np.min(data['res_ra']),
            'max_ra': np.max(data['res_ra']),
            'min_dec': np.min(data['res_dec']),
            'max_dec': np.max(data['res_dec'])
        }
        
        logger.info(f"Object {obj_id} original statistics:")
        logger.info(f"  Observations: {n_obs}")
        logger.info(f"  Mean RA residual: {mean_ra:.4f}\"")
        logger.info(f"  Mean DEC residual: {mean_dec:.4f}\"")
        logger.info(f"  Std RA residual: {std_ra:.4f}\"")
        logger.info(f"  Std DEC residual: {std_dec:.4f}\"")
        
        return stats
    
    def _remove_outliers(self, data: pd.DataFrame, 
                        original_stats: Dict[str, Any]) -> Tuple[pd.DataFrame, Dict[str, Any]]:
        """
        迭代剔除野值
        
        Args:
            data: 观测数据
            original_stats: 原始统计
            
        Returns:
            (filtered_data, final_stats): 过滤后数据和最终统计
        """
        current_data = data.copy()
        iteration = 0
        max_iterations = 10
        
        while iteration < max_iterations:
            iteration += 1
            n_before = len(current_data)
            
            if n_before < 2:
                break
                
            # 计算当前统计
            mean_ra = np.mean(current_data['res_ra'])
            mean_dec = np.mean(current_data['res_dec'])
            std_ra = np.std(current_data['res_ra'], ddof=1)
            std_dec = np.std(current_data['res_dec'], ddof=1)
            
            # 应用剔除条件
            if self.std_limit > 0:
                # 使用标准偏差模式
                mask = (
                    (np.abs(current_data['res_ra'] - mean_ra) < self.std_limit * std_ra) &
                    (np.abs(current_data['res_dec'] - mean_dec) < self.std_limit * std_dec) &
                    (np.abs(current_data['res_ra']) < self.eps) &
                    (np.abs(current_data['res_dec']) < self.eps)
                )
            else:
                # 使用均值限制模式
                mask = (
                    (np.abs(current_data['res_ra'] - mean_ra) < self.mean_limit) &
                    (np.abs(current_data['res_dec'] - mean_dec) < self.mean_limit)
                )
            
            filtered_data = current_data[mask]
            n_after = len(filtered_data)
            
            # 如果没有剔除任何数据，停止迭代
            if n_after == n_before:
                break
                
            current_data = filtered_data
            
            logger.debug(f"Iteration {iteration}: {n_before} -> {n_after} observations")
        
        # 计算最终统计
        final_stats = self._analyze_residuals(current_data, 0)
        final_stats['iterations'] = iteration
        final_stats['removed_count'] = len(data) - len(current_data)
        
        logger.info(f"Outlier removal completed after {iteration} iterations")
        logger.info(f"Removed {final_stats['removed_count']} outliers")
        logger.info(f"Final statistics:")
        logger.info(f"  Observations: {final_stats['n_obs']}")
        logger.info(f"  Mean RA residual: {final_stats['mean_ra']:.4f}\"")
        logger.info(f"  Mean DEC residual: {final_stats['mean_dec']:.4f}\"")
        logger.info(f"  Std RA residual: {final_stats['std_ra']:.4f}\"")
        logger.info(f"  Std DEC residual: {final_stats['std_dec']:.4f}\"")
        
        return current_data, final_stats
    
    def _save_results(self, filtered_data: pd.DataFrame, 
                     final_stats: Dict[str, Any],
                     original_file: str, obj_id: int):
        """
        保存处理结果
        
        Args:
            filtered_data: 过滤后的数据
            final_stats: 最终统计
            original_file: 原始文件路径
            obj_id: 目标ID
        """
        file_stem = Path(original_file).stem
        dir_path = Path(original_file).parent
        
        # 1. 保存到原始目录
        final_object_file = dir_path / f"final_object_{obj_id}.out"
        final_oc_file = dir_path / f"final_oc_{obj_id}.out"
        
        # 保存过滤后的观测数据
        with open(final_object_file, 'w') as f:
            f.write("# Filtered observation data after outlier removal\n")
            f.write("# Year Month Day_fraction RA_obs DEC_obs RA_eph DEC_eph Res_RA Res_DEC Sigma N_matched Angle Image\n")
            for _, row in filtered_data.iterrows():
                f.write(f"{row['original_line']}\n")
        
        # 保存统计结果
        with open(final_oc_file, 'w') as f:
            f.write("=================Residual Data Statistics=================\n")
            f.write(f"  Original data count: {len(filtered_data) + final_stats['removed_count']}\n")
            f.write(f"  After outlier removal: {final_stats['n_obs']}\n")
            f.write(f"  Removed outliers: {final_stats['removed_count']}\n")
            f.write(f"  Iterations: {final_stats['iterations']}\n")
            f.write(f"  Final mean RA residual: {final_stats['mean_ra']:.4f}\"\n")
            f.write(f"  Final mean DEC residual: {final_stats['mean_dec']:.4f}\"\n")
            f.write(f"  Final std RA residual: {final_stats['std_ra']:.4f}\"\n")
            f.write(f"  Final std DEC residual: {final_stats['std_dec']:.4f}\"\n")
            
            if self.std_limit > 0:
                f.write(f"*********Outlier removal criteria (std_limit): {self.std_limit}*********\n")
            else:
                f.write(f"*********Outlier removal criteria (mean_limit): {self.mean_limit}*********\n")
        
        # 2. 保存到指定输出目录
        date_str = self._extract_date_from_path(original_file)
        output_file = Path(self.output_dir) / f"{date_str}_{obj_id}.out"
        oc_output_file = Path(self.output_dir) / f"{date_str}_oc_{obj_id}.out"
        obsdata_file = Path(self.output_dir) / f"{date_str}_obsdata_{obj_id}.out"
        
        # 只保存质量好的数据
        if (final_stats['std_ra'] < 0.3 and final_stats['std_dec'] < 0.3 and 
            final_stats['n_obs'] > 1):
            
            # 完整数据
            shutil.copy2(final_object_file, output_file)
            
            # O-C统计
            with open(oc_output_file, 'w') as f:
                f.write(f"{date_str} & {len(filtered_data) + final_stats['removed_count']} & "
                       f"{final_stats['mean_ra']:.3f} & {final_stats['mean_dec']:.3f} & "
                       f"{final_stats['std_ra']:.3f} & {final_stats['std_dec']:.3f} & "
                       f"{final_stats['n_obs']} & {final_stats['mean_ra']:.3f} & "
                       f"{final_stats['mean_dec']:.3f} & {final_stats['std_ra']:.3f} & "
                       f"{final_stats['std_dec']:.3f}\n")
            
            # IMCCE格式数据
            with open(obsdata_file, 'w') as f:
                f.write("# IMCCE format observation data\n")
                for _, row in filtered_data.iterrows():
                    # 简化的IMCCE格式
                    f.write(f"{row['year']} {row['month']:02d} {row['day_fraction']:08.5f} "
                           f"{row['ra_obs']:.8f} {row['dec_obs']:.8f}\n")
        
        logger.info(f"Results saved for object {obj_id}")
    
    def _extract_date_from_path(self, file_path: str) -> str:
        """从文件路径提取日期字符串"""
        path = Path(file_path)
        
        # 尝试从路径中提取日期
        for part in path.parts:
            if len(part) == 8 and part.isdigit():  # YYYYMMDD格式
                return part
            elif len(part) == 10 and part[4] == '-' and part[7] == '-':  # YYYY-MM-DD格式
                return part.replace('-', '')
                
        # 如果找不到，使用文件名
        return path.stem.replace('object_', '').replace('_', '')
    
    def _generate_summary_report(self, results: Dict[int, Dict[str, Any]]):
        """
        生成汇总报告
        
        Args:
            results: 处理结果字典
        """
        summary_file = Path(self.output_dir) / "observation_summary.txt"
        
        with open(summary_file, 'w') as f:
            f.write("Observation Data Processing Summary\n")
            f.write("=" * 50 + "\n\n")
            
            total_original = sum(r['original_count'] for r in results.values())
            total_final = sum(r['final_count'] for r in results.values())
            total_removed = total_original - total_final
            
            f.write(f"Total objects processed: {len(results)}\n")
            f.write(f"Total original observations: {total_original}\n")
            f.write(f"Total final observations: {total_final}\n")
            f.write(f"Total removed outliers: {total_removed}\n")
            f.write(f"Overall retention rate: {total_final/total_original*100:.1f}%\n\n")
            
            f.write("Individual Object Results:\n")
            f.write("-" * 30 + "\n")
            
            for obj_id, result in results.items():
                orig_stats = result['original_stats']
                final_stats = result['final_stats']
                
                f.write(f"\nObject {obj_id}:\n")
                f.write(f"  File: {result['file_path']}\n")
                f.write(f"  Original: {result['original_count']} observations\n")
                f.write(f"    Mean residuals: ({orig_stats['mean_ra']:.4f}\", {orig_stats['mean_dec']:.4f}\")\n")
                f.write(f"    Std residuals: ({orig_stats['std_ra']:.4f}\", {orig_stats['std_dec']:.4f}\")\n")
                f.write(f"  Final: {result['final_count']} observations\n")
                f.write(f"    Mean residuals: ({final_stats['mean_ra']:.4f}\", {final_stats['mean_dec']:.4f}\")\n")
                f.write(f"    Std residuals: ({final_stats['std_ra']:.4f}\", {final_stats['std_dec']:.4f}\")\n")
                f.write(f"  Removed: {final_stats['removed_count']} outliers\n")
                f.write(f"  Quality: {'Good' if final_stats['std_ra'] < 0.3 and final_stats['std_dec'] < 0.3 else 'Poor'}\n")
        
        # 生成JSON格式的详细结果
        json_summary_file = Path(self.output_dir) / "observation_summary.json"
        
        # 转换numpy类型为Python原生类型
        def convert_numpy(obj):
            if isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, np.integer):
                return int(obj)
            return obj
        
        json_results = {}
        for obj_id, result in results.items():
            json_results[str(obj_id)] = {
                key: convert_numpy(value) if not isinstance(value, dict) 
                else {k: convert_numpy(v) for k, v in value.items()}
                for key, value in result.items()
            }
        
        with open(json_summary_file, 'w') as f:
            json.dump(json_results, f, indent=2)
        
        logger.info(f"Summary reports saved to {summary_file} and {json_summary_file}")
    
    def _clean_temporary_files(self, object_files: List[str]):
        """
        清理临时文件
        
        Args:
            object_files: 对象文件列表
        """
        if self.delete_flag != 1:
            return
            
        logger.info("Cleaning temporary files...")
        
        for obj_file in object_files:
            dir_path = Path(obj_file).parent
            
            try:
                # 删除处理过程文件
                for pattern in ['*_n.fit', '*.reg', '*.ref.reg']:
                    for file_path in dir_path.glob(pattern):
                        file_path.unlink()
                        logger.debug(f"Deleted: {file_path}")
                        
            except Exception as e:
                logger.warning(f"Error cleaning files in {dir_path}: {e}")
        
        logger.info("Temporary file cleanup completed")


def create_monthly_summary(results_dir: str, output_file: str = None) -> Dict[str, Any]:
    """
    创建月度汇总报告
    
    Args:
        results_dir: 结果目录
        output_file: 输出文件路径（可选）
        
    Returns:
        月度统计结果
    """
    results_path = Path(results_dir)
    
    if not results_path.exists():
        logger.error(f"Results directory does not exist: {results_dir}")
        return {}
    
    # 收集所有结果文件
    summary_files = list(results_path.glob("*_oc_*.out"))
    
    if not summary_files:
        logger.warning("No O-C summary files found")
        return {}
    
    monthly_data = {}
    
    for file_path in summary_files:
        try:
            with open(file_path, 'r') as f:
                content = f.read()
                
            # 提取统计信息
            # 这里需要根据实际的O-C文件格式来解析
            # 简化处理，假设每行包含日期和统计数据
            
        except Exception as e:
            logger.warning(f"Error reading {file_path}: {e}")
            continue
    
    # 生成月度报告
    if output_file is None:
        output_file = results_path / "monthly_summary.txt"
    
    with open(output_file, 'w') as f:
        f.write("Monthly Observation Summary\n")
        f.write("=" * 30 + "\n\n")
        
        # 这里添加月度统计的具体内容
        f.write("Monthly statistics will be implemented based on specific requirements.\n")
    
    logger.info(f"Monthly summary saved to {output_file}")
    return monthly_data


def process_observation_results(config: Dict[str, Any], 
                              object_files: List[str],
                              object_ids: List[int] = None) -> Dict[str, Any]:
    """
    处理观测结果的主函数
    
    Args:
        config: 配置参数
        object_files: 观测数据文件列表
        object_ids: 目标ID列表（可选）
        
    Returns:
        处理结果
    """
    if object_ids is None:
        object_ids = list(range(1, len(object_files) + 1))
    
    # 创建统计分析器
    analyzer = ObservationStatistics(config)
    
    # 处理数据
    results = analyzer.process_observation_data(object_files, object_ids)
    
    # 创建月度汇总
    create_monthly_summary(config.get('specified_output', 'output'))
    
    return results


# 为兼容性提供的函数接口
def analyze_residuals(residuals_ra: np.ndarray, residuals_dec: np.ndarray,
                     std_limit: float = 2.6, eps: float = 0.5) -> Dict[str, Any]:
    """
    分析残差数据（兼容性函数）
    
    Args:
        residuals_ra: RA残差数组
        residuals_dec: DEC残差数组  
        std_limit: 标准偏差限制倍数
        eps: 绝对值限制
        
    Returns:
        统计结果
    """
    # 创建临时DataFrame
    data = pd.DataFrame({
        'res_ra': residuals_ra,
        'res_dec': residuals_dec
    })
    
    # 创建临时配置
    config = {
        'standard_limit': std_limit,
        'eps': eps,
        'mean_limit': 0,
        'delete_flag': 0,
        'specified_output': 'temp'
    }
    
    # 创建分析器
    analyzer = ObservationStatistics(config)
    
    # 分析原始统计
    original_stats = analyzer._analyze_residuals(data, 0)
    
    # 剔除野值
    filtered_data, final_stats = analyzer._remove_outliers(data, original_stats)
    
    return {
        'original': original_stats,
        'final': final_stats,
        'filtered_indices': filtered_data.index.tolist()
    }