import numpy as np
from typing import Tuple, List, Optional, Dict, Any
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.gaia import Gaia
from astroquery.jplhorizons import Horizons
import pandas as pd
from scipy.optimize import leastsq
from scipy.interpolate import interp1d
from datetime import datetime, timedelta
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class CoordinateTransform:
    """坐标变换相关函数"""
    
    @staticmethod
    def radec_to_tangent(ra: np.ndarray, dec: np.ndarray, 
                        ra0: float, dec0: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        将赤经赤纬转换为切平面坐标（理想坐标）
        
        Args:
            ra: 赤经数组 (度)
            dec: 赤纬数组 (度) 
            ra0: 参考点赤经 (度)
            dec0: 参考点赤纬 (度)
            
        Returns:
            xi, eta: 切平面坐标 (弧度)
        """
        ra_rad = np.radians(ra)
        dec_rad = np.radians(dec)
        ra0_rad = np.radians(ra0)
        dec0_rad = np.radians(dec0)
        
        cos_dec = np.cos(dec_rad)
        sin_dec = np.sin(dec_rad)
        cos_dec0 = np.cos(dec0_rad)
        sin_dec0 = np.sin(dec0_rad)
        cos_dra = np.cos(ra_rad - ra0_rad)
        
        denominator = sin_dec * sin_dec0 + cos_dec * cos_dec0 * cos_dra
        
        xi = cos_dec * np.sin(ra_rad - ra0_rad) / denominator
        eta = (sin_dec * cos_dec0 - cos_dec * sin_dec0 * cos_dra) / denominator
        
        return xi, eta
    
    @staticmethod
    def tangent_to_radec(xi: np.ndarray, eta: np.ndarray,
                        ra0: float, dec0: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        将切平面坐标转换为赤经赤纬
        
        Args:
            xi, eta: 切平面坐标 (弧度)
            ra0: 参考点赤经 (度)
            dec0: 参考点赤纬 (度)
            
        Returns:
            ra, dec: 赤经赤纬 (度)
        """
        ra0_rad = np.radians(ra0)
        dec0_rad = np.radians(dec0)
        
        cos_dec0 = np.cos(dec0_rad)
        sin_dec0 = np.sin(dec0_rad)
        tan_dec0 = np.tan(dec0_rad)
        
        # 计算赤经
        ra_offset = np.arctan(xi / (cos_dec0 - eta * tan_dec0))
        ra = np.degrees(ra0_rad + ra_offset)
        
        # 计算赤纬
        dec_term = (eta + tan_dec0) * np.cos(ra_offset) / (1 - eta * tan_dec0)
        dec = np.degrees(np.arctan(dec_term))
        
        return ra, dec
    
    @staticmethod
    def xy_to_radec(x: np.ndarray, y: np.ndarray, params: np.ndarray,
                   ra0: float, dec0: float, focal_length: float, 
                   pixel_scale: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        将像素坐标转换为赤经赤纬
        
        Args:
            x, y: 像素坐标
            params: 底片常数
            ra0, dec0: 参考点赤经赤纬 (度)
            focal_length: 焦距 (mm)
            pixel_scale: 像素尺度 (弧秒/像素)
            
        Returns:
            ra, dec: 赤经赤纬 (度)
        """
        # 转换为角坐标 (弧度)
        x_ang = x * np.radians(pixel_scale / 3600) / focal_length
        y_ang = y * np.radians(pixel_scale / 3600) / focal_length
        
        # 应用底片常数
        if len(params) == 6:
            xi = params[0] * x_ang + params[1] * y_ang + params[2]
            eta = params[3] * x_ang + params[4] * y_ang + params[5]
        elif len(params) == 12:
            xi = (params[0] * x_ang + params[1] * y_ang + params[2] +
                  params[3] * x_ang**2 + params[4] * x_ang * y_ang + params[5] * y_ang**2)
            eta = (params[6] * x_ang + params[7] * y_ang + params[8] +
                   params[9] * x_ang**2 + params[10] * x_ang * y_ang + params[11] * y_ang**2)
        else:
            raise ValueError(f"Unsupported parameter count: {len(params)}")
        
        # 转换为赤经赤纬
        return CoordinateTransform.tangent_to_radec(xi, eta, ra0, dec0)


class PlateModel:
    """底片模型和常数求解"""
    
    @staticmethod
    def solve_plate_constants(x_obs: np.ndarray, y_obs: np.ndarray,
                            ra_cat: np.ndarray, dec_cat: np.ndarray,
                            ra0: float, dec0: float, model_type: int = 6,
                            max_iterations: int = 30, sigma_factor: float = 2.6) -> Tuple[np.ndarray, float, int, np.ndarray]:
        """
        求解底片常数
        
        Args:
            x_obs, y_obs: 观测的像素坐标
            ra_cat, dec_cat: 星表赤经赤纬 (度)
            ra0, dec0: 参考点赤经赤纬 (度)  
            model_type: 模型类型 (6, 12, 20, 30)
            max_iterations: 最大迭代次数
            sigma_factor: 剔野值因子
            
        Returns:
            params: 底片常数
            sigma: 拟合误差 (角秒)
            n_used: 使用的星数
            residuals: 残差
        """
        n_stars = len(x_obs)
        if n_stars < model_type:
            raise ValueError(f"Not enough stars ({n_stars}) for model type {model_type}")
        
        # 转换星表位置为理想坐标
        xi_cat, eta_cat = CoordinateTransform.radec_to_tangent(ra_cat, dec_cat, ra0, dec0)
        
        # 设置初始权重
        weights = np.ones(n_stars, dtype=bool)
        
        for iteration in range(max_iterations):
            # 选择参与计算的星
            n_used = np.sum(weights)
            if n_used < model_type:
                break
                
            x_use = x_obs[weights]
            y_use = y_obs[weights]
            xi_use = xi_cat[weights]
            eta_use = eta_cat[weights]
            
            # 构建设计矩阵
            A = PlateModel._build_design_matrix(x_use, y_use, model_type)
            b = np.concatenate([xi_use, eta_use])
            
            # 最小二乘求解
            params, residuals, rank, s = np.linalg.lstsq(A, b, rcond=None)
            
            # 计算残差
            xi_calc, eta_calc = PlateModel._apply_model(x_obs, y_obs, params, model_type)
            res_xi = xi_calc - xi_cat
            res_eta = eta_calc - eta_cat
            res_total = np.sqrt(res_xi**2 + res_eta**2)
            
            # 计算sigma
            if n_used > model_type:
                sigma = np.sqrt(np.sum(res_total[weights]**2) / (n_used - model_type))
            else:
                sigma = 0.0
                
            # 更新权重（剔除outliers）
            new_weights = res_total < sigma_factor * sigma
            if np.array_equal(weights, new_weights):
                break
            weights = new_weights
            
        # 转换sigma为角秒
        sigma_arcsec = np.degrees(sigma) * 3600
        
        return params, sigma_arcsec, n_used, res_total * 206265  # 转换为角秒
    
    @staticmethod
    def _build_design_matrix(x: np.ndarray, y: np.ndarray, model_type: int) -> np.ndarray:
        """构建设计矩阵"""
        n = len(x)
        
        if model_type == 6:
            A = np.zeros((2 * n, 6))
            # xi 方程
            A[0::2, 0] = x
            A[0::2, 1] = y  
            A[0::2, 2] = 1
            # eta 方程
            A[1::2, 3] = x
            A[1::2, 4] = y
            A[1::2, 5] = 1
            
        elif model_type == 12:
            A = np.zeros((2 * n, 12))
            # xi 方程
            A[0::2, 0] = x
            A[0::2, 1] = y
            A[0::2, 2] = 1
            A[0::2, 3] = x**2
            A[0::2, 4] = x * y
            A[0::2, 5] = y**2
            # eta 方程  
            A[1::2, 6] = x
            A[1::2, 7] = y
            A[1::2, 8] = 1
            A[1::2, 9] = x**2
            A[1::2, 10] = x * y
            A[1::2, 11] = y**2
            
        else:
            raise ValueError(f"Model type {model_type} not implemented")
            
        return A
    
    @staticmethod
    def _apply_model(x: np.ndarray, y: np.ndarray, params: np.ndarray, 
                    model_type: int) -> Tuple[np.ndarray, np.ndarray]:
        """应用底片模型"""
        if model_type == 6:
            xi = params[0] * x + params[1] * y + params[2]
            eta = params[3] * x + params[4] * y + params[5]
        elif model_type == 12:
            xi = (params[0] * x + params[1] * y + params[2] +
                  params[3] * x**2 + params[4] * x * y + params[5] * y**2)
            eta = (params[6] * x + params[7] * y + params[8] +
                   params[9] * x**2 + params[10] * x * y + params[11] * y**2)
        else:
            raise ValueError(f"Model type {model_type} not implemented")
            
        return xi, eta


class CatalogManager:
    """星表管理"""
    
    @staticmethod
    def load_gaia_catalog(catalog_file: str, min_mag: float = 5.0, 
                         max_mag: float = 17.0) -> Dict[str, np.ndarray]:
        """
        加载GAIA星表
        
        Args:
            catalog_file: 星表文件路径
            min_mag: 最小星等
            max_mag: 最大星等
            
        Returns:
            包含ra, dec, pmra, pmdec, mag的字典
        """
        try:
            # 尝试读取本地文件
            if Path(catalog_file).exists():
                return CatalogManager._load_local_gaia(catalog_file, min_mag, max_mag)
            else:
                logger.warning(f"Local catalog file {catalog_file} not found")
                return {"ra": np.array([]), "dec": np.array([]), 
                       "pmra": np.array([]), "pmdec": np.array([]), "mag": np.array([])}
        except Exception as e:
            logger.error(f"Error loading GAIA catalog: {e}")
            return {"ra": np.array([]), "dec": np.array([]), 
                   "pmra": np.array([]), "pmdec": np.array([]), "mag": np.array([])}
    
    @staticmethod
    def _load_local_gaia(catalog_file: str, min_mag: float, max_mag: float) -> Dict[str, np.ndarray]:
        """加载本地GAIA星表文件"""
        ra_list, dec_list, pmra_list, pmdec_list, mag_list = [], [], [], [], []
        
        with open(catalog_file, 'r') as f:
            # 跳过前60行注释
            for _ in range(60):
                f.readline()
                
            for line in f:
                if line.strip().startswith('#END') or not line.strip():
                    break
                    
                try:
                    # 解析格式: 21 17 38.5772643376 -15 30 01.349582279 319.41072917498 -15.50039051142 -2.110 -3.630 14.8833
                    parts = line.split()
                    if len(parts) >= 11:
                        ra = float(parts[2])  # 已经是度
                        dec = float(parts[3])  # 已经是度
                        pmra = float(parts[4])  # mas/yr
                        pmdec = float(parts[5])  # mas/yr
                        mag = float(parts[6])
                        
                        if min_mag <= mag <= max_mag:
                            ra_list.append(ra)
                            dec_list.append(dec)
                            pmra_list.append(pmra)
                            pmdec_list.append(pmdec)
                            mag_list.append(mag)
                            
                except (ValueError, IndexError):
                    continue
        
        return {
            "ra": np.array(ra_list),
            "dec": np.array(dec_list), 
            "pmra": np.array(pmra_list),
            "pmdec": np.array(pmdec_list),
            "mag": np.array(mag_list)
        }
    
    @staticmethod
    def apply_proper_motion(catalog: Dict[str, np.ndarray], 
                           epoch: float, reference_epoch: float = 2016.0) -> Dict[str, np.ndarray]:
        """
        应用自行修正
        
        Args:
            catalog: 星表数据
            epoch: 目标历元
            reference_epoch: 参考历元
            
        Returns:
            修正后的星表数据
        """
        dt = epoch - reference_epoch
        
        # 修正赤经赤纬
        ra_corrected = catalog["ra"] + (dt * catalog["pmra"] / 1000.0 / 3600.0 / np.cos(np.radians(catalog["dec"])))
        dec_corrected = catalog["dec"] + (dt * catalog["pmdec"] / 1000.0 / 3600.0)
        
        result = catalog.copy()
        result["ra"] = ra_corrected
        result["dec"] = dec_corrected
        
        return result


class EphemerisManager:
    """历表管理"""
    
    @staticmethod
    def load_ephemeris(eph_file: str, source: str = "IMCCE") -> Dict[str, np.ndarray]:
        """
        加载历表数据
        
        Args:
            eph_file: 历表文件路径
            source: 历表来源 ("IMCCE" 或 "JPL")
            
        Returns:
            包含time, ra, dec的字典
        """
        if source.upper() == "IMCCE":
            return EphemerisManager._load_imcce_ephemeris(eph_file)
        elif source.upper() == "JPL":
            return EphemerisManager._load_jpl_ephemeris(eph_file)
        else:
            raise ValueError(f"Unknown ephemeris source: {source}")
    
    @staticmethod
    def _load_imcce_ephemeris(eph_file: str) -> Dict[str, np.ndarray]:
        """加载IMCCE格式历表"""
        times, ra_list, dec_list = [], [], []
        
        try:
            with open(eph_file, 'r') as f:
                # 跳过前10行
                for _ in range(10):
                    f.readline()
                    
                for line in f:
                    if line.startswith('---'):
                        continue
                        
                    try:
                        parts = line.split()
                        if len(parts) >= 8:
                            year, month, day = int(parts[0]), int(parts[1]), int(parts[2])
                            hour, minute, second = int(parts[3]), int(parts[4]), float(parts[5])
                            ra, dec = float(parts[6]) * 15.0, float(parts[7])  # RA转为度
                            
                            # 转换为天的小数
                            time_decimal = day + (hour * 3600 + minute * 60 + second) / 86400.0
                            
                            times.append(time_decimal)
                            ra_list.append(ra)
                            dec_list.append(dec)
                            
                    except (ValueError, IndexError):
                        continue
                        
        except FileNotFoundError:
            logger.warning(f"Ephemeris file {eph_file} not found")
            
        return {
            "time": np.array(times),
            "ra": np.array(ra_list),
            "dec": np.array(dec_list)
        }
    
    @staticmethod
    def _load_jpl_ephemeris(eph_file: str) -> Dict[str, np.ndarray]:
        """加载JPL格式历表"""
        times, ra_list, dec_list = [], [], []
        
        month_map = {
            'Jan': 1, 'Feb': 2, 'Mar': 3, 'Apr': 4, 'May': 5, 'Jun': 6,
            'Jul': 7, 'Aug': 8, 'Sep': 9, 'Oct': 10, 'Nov': 11, 'Dec': 12
        }
        
        try:
            with open(eph_file, 'r') as f:
                # 跳过前40行
                for _ in range(40):
                    f.readline()
                    
                for line in f:
                    if line.startswith('$$EOE'):
                        break
                        
                    try:
                        # 解析JPL格式
                        # 格式示例: 2006-Oct-08 00:00     m  02 59 46.96 +15 41 36.2
                        if len(line) > 50:
                            date_part = line[1:12].strip()
                            time_part = line[13:18].strip()
                            ra_part = line[25:35].strip()
                            dec_part = line[37:47].strip()
                            
                            # 解析日期
                            year = int(date_part[:4])
                            month_str = date_part[5:8]
                            day = int(date_part[9:11])
                            month = month_map[month_str]
                            
                            # 解析时间
                            hour, minute = map(int, time_part.split(':'))
                            
                            # 解析赤经 (时分秒)
                            ra_parts = ra_part.split()
                            ra_hour = int(ra_parts[0])
                            ra_min = int(ra_parts[1])
                            ra_sec = float(ra_parts[2])
                            ra = (ra_hour + ra_min/60.0 + ra_sec/3600.0) * 15.0
                            
                            # 解析赤纬 (度分秒)
                            dec_parts = dec_part.split()
                            dec_sign = 1 if dec_parts[0][0] == '+' else -1
                            dec_deg = int(dec_parts[0][1:])
                            dec_min = int(dec_parts[1])
                            dec_sec = float(dec_parts[2])
                            dec = dec_sign * (dec_deg + dec_min/60.0 + dec_sec/3600.0)
                            
                            # 转换为天的小数
                            time_decimal = day + (hour * 3600 + minute * 60) / 86400.0
                            
                            times.append(time_decimal)
                            ra_list.append(ra)
                            dec_list.append(dec)
                            
                    except (ValueError, IndexError):
                        continue
                        
        except FileNotFoundError:
            logger.warning(f"Ephemeris file {eph_file} not found")
            
        return {
            "time": np.array(times),
            "ra": np.array(ra_list), 
            "dec": np.array(dec_list)
        }
    
    @staticmethod
    def interpolate_position(ephemeris: Dict[str, np.ndarray], 
                           target_time: float) -> Tuple[float, float]:
        """
        内插指定时间的天体位置
        
        Args:
            ephemeris: 历表数据
            target_time: 目标时间 (天的小数)
            
        Returns:
            ra, dec: 内插得到的赤经赤纬 (度)
        """
        if len(ephemeris["time"]) == 0:
            return 0.0, 0.0
            
        # 创建内插函数
        ra_interp = interp1d(ephemeris["time"], ephemeris["ra"], 
                            kind='linear', fill_value='extrapolate')
        dec_interp = interp1d(ephemeris["time"], ephemeris["dec"],
                             kind='linear', fill_value='extrapolate')
        
        return float(ra_interp(target_time)), float(dec_interp(target_time))


def stellar_match(detected_stars: Tuple[np.ndarray, ...], 
                 catalog: Dict[str, np.ndarray],
                 ephemeris: Dict[str, np.ndarray],
                 observation_time: float,
                 telescope_params: Dict[str, Any],
                 match_config: Dict[str, Any]) -> Dict[str, Any]:
    """
    主要的星表匹配函数
    
    Args:
        detected_stars: 检测到的星象 (x, y, flux, snr, npix)
        catalog: 星表数据
        ephemeris: 历表数据  
        observation_time: 观测时间 (天的小数)
        telescope_params: 望远镜参数
        match_config: 匹配配置
        
    Returns:
        匹配结果字典
    """
    x_det, y_det, flux_det, snr_det, npix_det = detected_stars
    
    if len(x_det) == 0:
        logger.warning("No detected stars provided")
        return {"success": False, "message": "No detected stars"}
    
    # 内插目标位置
    target_ra, target_dec = EphemerisManager.interpolate_position(ephemeris, observation_time)
    if target_ra == 0 and target_dec == 0:
        logger.warning("Could not interpolate target position")
        return {"success": False, "message": "Could not interpolate target position"}
    
    # 应用自行修正
    epoch = 2020.0 + observation_time / 365.25  # 简化的历元计算
    catalog_corrected = CatalogManager.apply_proper_motion(catalog, epoch)
    
    # 选择视场内的星
    field_size = match_config.get("field_size", 15) / 60.0  # 转换为度
    field_mask = ((np.abs(catalog_corrected["ra"] - target_ra) < field_size) &
                  (np.abs(catalog_corrected["dec"] - target_dec) < field_size))
    
    if not np.any(field_mask):
        logger.warning("No catalog stars in field")
        return {"success": False, "message": "No catalog stars in field"}
    
    field_catalog = {key: val[field_mask] for key, val in catalog_corrected.items()}
    logger.info(f"Found {len(field_catalog['ra'])} catalog stars in field")
    
    # 尝试不同的旋转角度找到最佳匹配
    plate_angle = match_config.get("plate_angle", 0)
    if plate_angle == 0:
        # 自动确定旋转角度
        best_result = _determine_rotation_angle(
            x_det, y_det, field_catalog, target_ra, target_dec,
            telescope_params, match_config
        )
    else:
        # 使用指定角度
        best_result = _match_with_angle(
            x_det, y_det, field_catalog, target_ra, target_dec,
            plate_angle, telescope_params, match_config
        )
    
    if best_result["n_matched"] < 3:
        return {"success": False, "message": "Insufficient matched stars"}
    
    # 计算目标星位置
    target_result = _find_target_position(
        x_det, y_det, flux_det, snr_det,
        target_ra, target_dec, best_result,
        telescope_params, match_config
    )
    
    result = {
        "success": True,
        "n_matched": best_result["n_matched"],
        "plate_angle": best_result["plate_angle"],
        "plate_constants": best_result["plate_constants"],
        "sigma": best_result["sigma"],
        "target_x": target_result["target_x"],
        "target_y": target_result["target_y"],
        "target_ra_obs": target_result["target_ra_obs"],
        "target_dec_obs": target_result["target_dec_obs"],
        "target_ra_eph": target_ra,
        "target_dec_eph": target_dec,
        "residual_ra": target_result["residual_ra"],
        "residual_dec": target_result["residual_dec"],
        "matched_stars": best_result["matched_stars"]
    }
    
    return result


def _determine_rotation_angle(x_det: np.ndarray, y_det: np.ndarray,
                            catalog: Dict[str, np.ndarray],
                            target_ra: float, target_dec: float,
                            telescope_params: Dict[str, Any],
                            match_config: Dict[str, Any]) -> Dict[str, Any]:
    """确定最佳旋转角度"""
    best_match = {"n_matched": 0, "sigma": float('inf')}
    
    # 粗匹配：3度步长
    for angle in range(0, 360, 3):
        result = _match_with_angle(
            x_det, y_det, catalog, target_ra, target_dec,
            angle, telescope_params, match_config, coarse=True
        )
        
        if (result["n_matched"] > best_match["n_matched"] and 
            result["sigma"] < 1.0 and result["sigma"] > 0.001):
            best_match = result
    
    if best_match["n_matched"] > 6:
        # 细匹配：在最佳角度附近1度步长
        base_angle = best_match["plate_angle"]
        for offset in range(-3, 4):
            angle = base_angle + offset
            result = _match_with_angle(
                x_det, y_det, catalog, target_ra, target_dec,
                angle, telescope_params, match_config, coarse=False
            )
            
            if (result["n_matched"] >= best_match["n_matched"] and 
                result["sigma"] < best_match["sigma"]):
                best_match = result
    
    return best_match


def _match_with_angle(x_det: np.ndarray, y_det: np.ndarray,
                     catalog: Dict[str, np.ndarray],
                     target_ra: float, target_dec: float,
                     plate_angle: float,
                     telescope_params: Dict[str, Any],
                     match_config: Dict[str, Any],
                     coarse: bool = False) -> Dict[str, Any]:
    """使用指定角度进行匹配"""
    
    # 初始化旋转矩阵参数
    angle_rad = np.radians(plate_angle)
    cos_a, sin_a = np.cos(angle_rad), np.sin(angle_rad)
    
    initial_params = np.array([cos_a, 0, 0, sin_a, 1, 0])
    
    focal_length = telescope_params["focal_length"]
    pixel_scale = telescope_params["pixel_scale"]
    match_limit = match_config["match_limit"] * (1.5 if coarse else 1.0)
    
    # 只使用前100颗最亮的检测星
    n_use = min(100, len(x_det))
    x_use = x_det[:n_use]
    y_use = y_det[:n_use]
    
    best_match = {"n_matched": 0, "sigma": float('inf')}
    
    # 尝试每颗检测星作为可能的目标星
    for i in range(min(80, n_use)):  # 假设目标在前80颗中
        x_center, y_center = x_use[i], y_use[i]
        
        # 计算相对坐标
        x_rel = (x_use - x_center) * pixel_scale / focal_length * np.pi / 180 / 3600
        y_rel = (y_use - y_center) * pixel_scale / focal_length * np.pi / 180 / 3600
        
        # 应用初始变换得到大致的天球坐标
        xi_approx = initial_params[0] * x_rel + initial_params[1] * y_rel
        eta_approx = initial_params[3] * x_rel + initial_params[4] * y_rel
        
        ra_approx, dec_approx = CoordinateTransform.tangent_to_radec(
            xi_approx, eta_approx, target_ra, target_dec
        )
        
        # 与星表匹配
        matched_indices = []
        matched_det_indices = []
        
        for j, (ra_det, dec_det) in enumerate(zip(ra_approx, dec_approx)):
            # 计算与星表中每颗星的角距离
            sep = np.sqrt((catalog["ra"] - ra_det)**2 * np.cos(np.radians(catalog["dec"]))**2 +
                         (catalog["dec"] - dec_det)**2) * 3600  # 转换为角秒
            
            min_idx = np.argmin(sep)
            if sep[min_idx] < match_limit:
                matched_indices.append(min_idx)
                matched_det_indices.append(j)
        
        n_matched = len(matched_indices)
        
        if n_matched >= 3:
            try:
                # 求解底片常数
                x_matched = x_rel[matched_det_indices]
                y_matched = y_rel[matched_det_indices]
                ra_matched = catalog["ra"][matched_indices]
                dec_matched = catalog["dec"][matched_indices]
                
                params, sigma, n_used, residuals = PlateModel.solve_plate_constants(
                    x_matched, y_matched, ra_matched, dec_matched,
                    target_ra, target_dec, match_config.get("model_type", 6)
                )
                
                if n_matched > best_match["n_matched"] and sigma < 1.0 and sigma > 0.001:
                    best_match = {
                        "n_matched": n_matched,
                        "sigma": sigma,
                        "plate_angle": plate_angle,
                        "plate_constants": params,
                        "target_candidate_idx": i,
                        "matched_stars": {
                            "x_det": x_use[matched_det_indices],
                            "y_det": y_use[matched_det_indices],
                            "ra_cat": catalog["ra"][matched_indices],
                            "dec_cat": catalog["dec"][matched_indices],
                            "mag_cat": catalog["mag"][matched_indices]
                        }
                    }
                    
            except Exception as e:
                logger.debug(f"Failed to solve plate constants: {e}")
                continue
    
    return best_match


def _find_target_position(x_det: np.ndarray, y_det: np.ndarray,
                         flux_det: np.ndarray, snr_det: np.ndarray,
                         target_ra_eph: float, target_dec_eph: float,
                         match_result: Dict[str, Any],
                         telescope_params: Dict[str, Any],
                         match_config: Dict[str, Any]) -> Dict[str, Any]:
    """寻找目标星位置"""
    
    params = match_result["plate_constants"]
    focal_length = telescope_params["focal_length"]
    pixel_scale = telescope_params["pixel_scale"]
    
    # 使用历表位置预测像素位置
    matched_stars = match_result["matched_stars"]
    x_matched = matched_stars["x_det"]
    y_matched = matched_stars["y_det"]
    
    # 计算参考星中心
    x_center = np.mean(x_matched)
    y_center = np.mean(y_matched)
    
    # 预测目标位置
    xi_target, eta_target = CoordinateTransform.radec_to_tangent(
        np.array([target_ra_eph]), np.array([target_dec_eph]),
        target_ra_eph, target_dec_eph
    )
    
    # 反解像素坐标（简化处理）
    x_pred = x_center  # 简化假设目标在参考星中心附近
    y_pred = y_center
    
    # 在预测位置附近寻找最近的检测星
    distances = np.sqrt((x_det - x_pred)**2 + (y_det - y_pred)**2)
    closest_idx = np.argmin(distances)
    
    if distances[closest_idx] < 3.0:  # 3像素内
        target_x = x_det[closest_idx]
        target_y = y_det[closest_idx]
        target_flux = flux_det[closest_idx]
        target_snr = snr_det[closest_idx]
    else:
        target_x = x_pred
        target_y = y_pred
        target_flux = 0.0
        target_snr = 0.0
    
    # 计算观测位置
    x_rel = (target_x - x_center) * pixel_scale / focal_length * np.pi / 180 / 3600
    y_rel = (target_y - y_center) * pixel_scale / focal_length * np.pi / 180 / 3600
    
    xi_obs, eta_obs = PlateModel._apply_model(
        np.array([x_rel]), np.array([y_rel]), params, len(params)
    )
    
    target_ra_obs, target_dec_obs = CoordinateTransform.tangent_to_radec(
        xi_obs, eta_obs, target_ra_eph, target_dec_eph
    )
    
    # 计算残差 (角秒)
    residual_ra = (target_ra_obs[0] - target_ra_eph) * 3600 * np.cos(np.radians(target_dec_eph))
    residual_dec = (target_dec_obs[0] - target_dec_eph) * 3600
    
    return {
        "target_x": target_x,
        "target_y": target_y,
        "target_ra_obs": target_ra_obs[0],
        "target_dec_obs": target_dec_obs[0],
        "residual_ra": residual_ra,
        "residual_dec": residual_dec,
        "target_flux": target_flux,
        "target_snr": target_snr
    }


def match_multiple_objects(detected_stars: Tuple[np.ndarray, ...],
                          catalogs: List[Dict[str, np.ndarray]],
                          ephemerides: List[Dict[str, np.ndarray]], 
                          observation_time: float,
                          telescope_params: Dict[str, Any],
                          match_config: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    多目标匹配
    
    Args:
        detected_stars: 检测到的星象
        catalogs: 多个星表数据
        ephemerides: 多个历表数据
        observation_time: 观测时间
        telescope_params: 望远镜参数
        match_config: 匹配配置
        
    Returns:
        多个目标的匹配结果列表
    """
    results = []
    
    for i, (catalog, ephemeris) in enumerate(zip(catalogs, ephemerides)):
        logger.info(f"Processing object {i+1}/{len(catalogs)}")
        
        result = stellar_match(
            detected_stars, catalog, ephemeris, observation_time,
            telescope_params, match_config
        )
        
        result["object_id"] = i + 1
        results.append(result)
    
    return results