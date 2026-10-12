"""Wavelength-resolved diffraction for the calculated thermal disk image.

At each wavelength this module evaluates the *scalar* observed Planck
radiance, averages detector subpixels in linear light, convolves it with the
corresponding circular-pupil Airy PSF, and only then integrates CIE XYZ.
It never assigns an effective wavelength to an RGB channel. The user supplies
T_obs = g*f_col*T_eff and dilution = Upsilon(mu)/f_col**4, so applying another
g**3 or g**4 here would count the relativistic transfer twice.

Sources: Airy (1835), circular-pupil diffraction; official CIE 1931 2-degree
observer data, CIE (2019), doi:10.25039/CIE.DS.xvudnb9b; and ICC's sRGB
interpretation, section 7. The atmosphere and relativistic transfer belong
to the existing emission/optics modules and their separately cited papers.

This is a hypothetical diffraction-limited camera. It adds no gravitational
halo, atmospheric seeing, detector noise, arbitrary bloom, or lens flare.
Finite-support kernels are normalized to unit sum. Inputs outside the supplied
field are zero, so render an overscan margin before cropping the final frame.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from PIL import Image

from emission import (C, H, KB, CIE_FILENAME, CIE_MD5, XYZ_TO_LINEAR_SRGB,
                      ThermalSpectrum, planck_lambda)
from optics import bessel_j1, _next_fast_length


AIRY_FIRST_ZERO = 3.8317059702075125
FIRST_NULL_FACTOR = AIRY_FIRST_ZERO / math.pi
WIEN_WAVELENGTH_T = 2.897771955185172e-3  # B_lambda maximum, m K


@lru_cache(maxsize=4)
def _thermal_basis(wavelength_tuple, rank=18):
    """A measured low-rank factorization of scalar Planck spectra.

    Normalization is only an internal conditioning operation: the absolute
    Planck amplitude is restored for every source subpixel before averaging.
    Interpolate smooth normalized coefficients, never temperature or RGB.
    """
    wavelengths = np.asarray(wavelength_tuple)*1e-9
    # This is the existing planck_lambda exponent cutoff, not a new cold-light
    # threshold. Below it every sampled visible wavelength is exactly zero in
    # that implementation. Stay infinitesimally above it when normalizing to
    # avoid dividing by zero at the boundary itself.
    zero_temperature = H*C/(wavelengths[-1]*KB*700.)
    lower_temperature = zero_temperature*(1.+1e-12)
    training_t = np.geomspace(lower_temperature, 1e5, 4096)
    def normalized(t):
        peak_wavelength = np.clip(WIEN_WAVELENGTH_T/t, wavelengths[0], wavelengths[-1])
        scale = planck_lambda(peak_wavelength, t)
        return planck_lambda(wavelengths[None, :], t[:, None])/scale[:, None]
    _, singular, vectors = np.linalg.svd(normalized(training_t), full_matrices=False)
    rank = min(rank, len(wavelengths))
    basis = vectors[:rank]
    log_grid = np.linspace(math.log(lower_temperature), math.log(1e5), 32769)
    coefficients = normalized(np.exp(log_grid)) @ basis.T
    # Half-grid samples independently measure both SVD and interpolation error.
    validation_log = (log_grid[:-1]+log_grid[1:])/2
    exact = normalized(np.exp(validation_log))
    reconstructed = ((coefficients[:-1]+coefficients[1:])/2) @ basis
    peak_error = float(np.max(np.abs(reconstructed-exact)))
    return basis, log_grid, coefficients, {
        'method': 'normalized scalar Planck SVD; absolute amplitude restored before BOX average; RGB-weighted spectral OTF sums',
        'rank': rank,
        'temperature_domain_k': [zero_temperature, 100000.],
        'exact_zero_temperature_k': zero_temperature,
        'exact_zero_rule': 'same exponent >= 700 cutoff already used by direct planck_lambda; no additional light truncation',
        'temperature_training_samples': len(training_t),
        'temperature_lookup_samples': len(log_grid),
        'independent_half_grid_validation_samples': len(validation_log),
        'max_monochromatic_error_relative_to_continuous_Planck_peak': peak_error,
        'discarded_singular_value': float(singular[rank]) if rank < len(singular) else 0.,
        'out_of_domain_policy': 'direct wavelength-resolved fallback; zero sentinel retained',
        'cold_scene_policy': 'direct fallback when all emitting temperatures are below 2000 K',
        'FFT_count_per_frame': {'forward': rank, 'inverse': 3},
        'optical_model_changed': False,
    }


def spectral_grid(step_nm=10.0):
    """Wavelengths, CIE CMFs and trapezoidal weights (metres), 360--830 nm.

    Endpoint weights are half their adjoining interval. For a step that does
    not divide 470 nm, the final shorter interval has its actual width.
    The official 1 nm data are linearly interpolated at noninteger wavelengths.
    """
    if not math.isfinite(step_nm) or not 0 < step_nm <= 470:
        raise ValueError('wavelength_step_nm must be in (0, 470]')
    raw = CIE_FILENAME.read_bytes()
    if hashlib.md5(raw).hexdigest() != CIE_MD5:
        raise ValueError('CIE table checksum disagrees with official dataset')
    table = np.loadtxt(CIE_FILENAME, delimiter=',')
    if table.shape != (471, 4) or not np.array_equal(table[:, 0], np.arange(360, 831)):
        raise ValueError('unexpected CIE table wavelength grid')
    wavelengths_nm = np.r_[np.arange(360., 830., step_nm), 830.]
    wavelengths_m = wavelengths_nm * 1e-9
    intervals = np.diff(wavelengths_m)
    weights = np.r_[intervals[0]/2, (intervals[:-1]+intervals[1:])/2, intervals[-1]/2]
    cmf = np.column_stack([np.interp(wavelengths_nm, table[:, 0], table[:, c])
                           for c in (1, 2, 3)])
    return wavelengths_nm, weights, cmf


@lru_cache(maxsize=4)
def _bessel_lookup(maximum_argument):
    """Independently evaluated J1 table; linear interpolation error is tested.

    J1's defining quadrature is shared with optics.py. A <=0.001 argument
    spacing resolves the oscillations well, without reevaluating the integral
    for every pixel in every wavelength's kernel.
    """
    count = max(65537, int(math.ceil(maximum_argument/.001))+1)
    x = np.linspace(0., maximum_argument, count)
    y = bessel_j1(x)
    y[0] = 0.
    return x, y


def spectral_airy_kernel(first_null_pixels, support_first_nulls=24., subsamples=3, *, _use_lookup=True):
    """Pixel-integrated unit-sum Airy PSF, evaluated with a checked J1 LUT.

    Gauss--Legendre quadrature integrates intensity over each detector pixel.
    The minimum nodes per axis are ``subsamples``; undersampled PSFs use more.
    The
    retained circular support has radius ceil(first_null*support) pixels;
    its omitted Airy wings are an explicit small finite-support approximation.
    """
    if first_null_pixels <= 0 or support_first_nulls < 2 or subsamples < 1:
        raise ValueError('invalid Airy kernel parameters')
    radius = max(2, int(math.ceil(first_null_pixels*support_first_nulls)))
    y, x = np.mgrid[-radius:radius+1, -radius:radius+1]
    kernel = np.zeros(x.shape, dtype=np.float64)
    node_count = max(subsamples, int(math.ceil(3/first_null_pixels)))
    offsets, quadrature_weights = np.polynomial.legendre.leggauss(node_count)
    offsets, quadrature_weights = offsets/2, quadrature_weights/2
    # Use one support-derived lookup for every wavelength; radius rounding
    # can add one pixel, hence this conservative argument allowance.
    limit = AIRY_FIRST_ZERO * (support_first_nulls + 2/first_null_pixels)
    # Round upward to a common coarse bound so nearby kernels reuse the table.
    table_limit = float(math.ceil(limit/32)*32)
    tx, ty = _bessel_lookup(table_limit)
    for dy, wy in zip(offsets, quadrature_weights):
        for dx, wx in zip(offsets, quadrature_weights):
            distance = np.hypot(x+dx, y+dy)
            arg = AIRY_FIRST_ZERO*distance/first_null_pixels
            # Values outside the circular cutoff do not contribute. Clip their
            # lookup argument only to avoid unnecessary tables for the corners.
            j1 = (np.interp(np.minimum(arg, tx[-1]), tx, ty) if _use_lookup
                  else bessel_j1(arg))
            ratio = np.ones_like(arg)
            np.divide(2*j1, arg, out=ratio, where=arg > 1e-10)
            intensity = ratio**2
            intensity[distance > radius+.5] = 0.
            kernel += intensity*wx*wy
    kernel /= kernel.sum()
    return kernel


class SpectralCamera:
    """Reusable wavelength-dependent clear-aperture camera.

    ``width,height`` describe detector pixels, including any overscan supplied
    by the caller. ``fov_rg`` is the corresponding *full* horizontal field.
    Preserve fov/width when expanding an image by ``padding_pixels`` on each
    edge. ``apply`` accepts larger same-aspect-ratio source arrays: monochromatic
    radiance, not temperature, is then BOX-averaged into detector pixels.

    The cached Fourier kernels are real because the discrete PSFs are centrally
    symmetric. Float32 storage saves memory; radiance and color sums use float64.
    This finite precision is checked against float64 direct convolution.
    """
    def __init__(self, width, height, fov_rg, aperture_m=10.,
                 observer_distance_rg=1e6, wavelength_step_nm=10.,
                 support_first_nulls=24., kernel_subsamples=3):
        if not isinstance(width, int) or not isinstance(height, int) or min(width, height) < 8:
            raise ValueError('detector dimensions must be integers >=8')
        numbers = (fov_rg, aperture_m, observer_distance_rg, support_first_nulls)
        if not all(math.isfinite(x) and x > 0 for x in numbers) or support_first_nulls < 2:
            raise ValueError('invalid field or circular-aperture parameters')
        if not isinstance(kernel_subsamples, int) or kernel_subsamples < 1:
            raise ValueError('kernel_subsamples must be a positive integer')
        self.width, self.height = width, height
        self.fov_rg = float(fov_rg)
        self.aperture_m = float(aperture_m)
        self.observer_distance_rg = float(observer_distance_rg)
        self.wavelengths_nm, weights, self.cmf = spectral_grid(wavelength_step_nm)
        self.wavelengths_m = self.wavelengths_nm*1e-9
        self.weights_m = weights
        self.xyz_weights = weights[:, None]*self.cmf
        self.rgb_weights = self.xyz_weights @ XYZ_TO_LINEAR_SRGB.T
        self.first_nulls_pixels = FIRST_NULL_FACTOR*self.wavelengths_m/aperture_m / (
            fov_rg/(observer_distance_rg*width))
        self.padding_pixels = max(2, int(math.ceil(self.first_nulls_pixels.max()*support_first_nulls)))
        self.fft_shape = (_next_fast_length(height+2*self.padding_pixels),
                          _next_fast_length(width+2*self.padding_pixels))
        self.kernels = []
        self._transfers = []
        started = time.perf_counter()
        for first_null in self.first_nulls_pixels:
            kernel = spectral_airy_kernel(float(first_null), support_first_nulls, kernel_subsamples)
            self.kernels.append(kernel)
            radius = kernel.shape[0]//2
            centered = np.zeros(self.fft_shape, dtype=np.float64)
            # Center at Fourier origin. Negative offsets wrap into the padded
            # region; cropping [:height,:width] then gives the same linear
            # convolution as a conventional full convolution plus origin crop.
            positions_y = np.arange(-radius, radius+1) % self.fft_shape[0]
            positions_x = np.arange(-radius, radius+1) % self.fft_shape[1]
            centered[np.ix_(positions_y, positions_x)] = kernel
            transform = np.fft.rfft2(centered)
            if np.abs(transform.imag).max() > 1e-10:
                raise AssertionError('centrally symmetric PSF has a complex OTF')
            self._transfers.append(transform.real.astype('float32'))
        self.metadata = {
            'model': 'wavelength-resolved Airy intensity PSF of an unobstructed circular aperture',
            'aperture_m': aperture_m,
            'observer_distance_rg': observer_distance_rg,
            'field_horizontal_rg': fov_rg,
            'detector_width_pixels': width,
            'detector_height_pixels': height,
            'wavelength_nm': [float(self.wavelengths_nm[0]), float(self.wavelengths_nm[-1])],
            'wavelength_step_nm': wavelength_step_nm,
            'wavelength_samples': len(self.wavelengths_nm),
            'integration': 'actual-wavelength trapezoidal quadrature, including half endpoint weights',
            'order': 'scalar Planck radiance; detector subpixel average; per-wavelength PSF; CIE XYZ; linear sRGB',
            'first_null_pixels_360_830_nm': [float(self.first_nulls_pixels[0]), float(self.first_nulls_pixels[-1])],
            'kernel_quadrature': 'Gauss-Legendre integration over each detector pixel; more nodes for undersampled kernels',
            'kernel_minimum_quadrature_nodes_per_axis': kernel_subsamples,
            'kernel_support_first_nulls': support_first_nulls,
            'recommended_overscan_pixels_each_edge': self.padding_pixels,
            'boundary': 'zero outside supplied input field; crop an overscan margin to avoid missing external light',
            'finite_support': 'pixel-integrated circular cutoff, normalized to unit sum at each wavelength',
            'kernel_fourier_storage': 'real float32; radiance/color accumulation float64',
            'kernel_cache_bytes': sum(x.nbytes for x in self._transfers),
            'setup_seconds': time.perf_counter()-started,
            'astrophysical_effect': False,
            'cie_doi': '10.25039/CIE.DS.xvudnb9b',
            'cie_md5': CIE_MD5,
            'display_matrix': 'ICC sRGB interpretation, section 7',
            'exposure_or_tone_mapping': 'none; caller applies one fixed exposure and display transform',
        }

    def _convolve(self, radiance, wavelength_index):
        full = np.fft.irfft2(np.fft.rfft2(radiance, s=self.fft_shape) *
                             self._transfers[wavelength_index], s=self.fft_shape)
        image = full[:self.height, :self.width].copy()
        # Physical scalar light and the Airy intensity kernel are nonnegative.
        # Remove negative FFT roundoff before signed color-matching conversion.
        np.maximum(image, 0., out=image)
        return image

    def apply(self, temperature_obs, dilution=1., *, return_ideal=False, diffraction=True):
        """Return absolute linear sRGB radiance, or both camera and ideal views.

        ``return_ideal=True`` returns ``{'camera': rgb, 'ideal': rgb}`` from
        one spectral pass. ``diffraction=False`` returns the ideal spectrum
        alone. No exposure or display clipping is performed; small negative
        out-of-gamut sRGB primaries are intentionally retained.
        """
        temperature, amplitude = np.broadcast_arrays(np.asarray(temperature_obs, dtype=np.float64),
                                                     np.asarray(dilution, dtype=np.float64))
        if temperature.ndim != 2 or not np.isfinite(temperature).all() or (temperature < 0).any():
            raise ValueError('temperature must be a finite nonnegative 2D image')
        if not np.isfinite(amplitude).all() or (amplitude < 0).any():
            raise ValueError('dilution must be finite and nonnegative')
        sh, sw = temperature.shape
        if sh < self.height or sw < self.width:
            raise ValueError('source map must be at least detector resolution')
        if abs(sw/self.width-sh/self.height) > max(1/self.width, 1/self.height)+1e-12:
            raise ValueError('source and detector aspect ratios must agree to within one source pixel')
        needs_resize = temperature.shape != (self.height, self.width)
        camera = np.zeros((self.height, self.width, 3), dtype=np.float64)
        ideal = np.zeros_like(camera) if return_ideal and diffraction else None
        # Sub-kelvin Planck radiance is below float64 underflow throughout this
        # visible band. A 1 K floor keeps zero-temperature sentinel pixels from
        # producing division warnings while still returning exactly zero light.
        evaluation_temperature = np.maximum(temperature, 1.)
        started = time.perf_counter()
        for index, wavelength in enumerate(self.wavelengths_m):
            monochromatic = amplitude * planck_lambda(wavelength, evaluation_temperature)
            if needs_resize:
                monochromatic = np.asarray(Image.fromarray(monochromatic.astype('float32')).resize(
                    (self.width, self.height), Image.Resampling.BOX), dtype=np.float64)
            if ideal is not None:
                for channel in range(3):
                    ideal[..., channel] += monochromatic*self.rgb_weights[index, channel]
            response = self._convolve(monochromatic, index) if diffraction else monochromatic
            for channel in range(3):
                camera[..., channel] += response*self.rgb_weights[index, channel]
        self.metadata['last_apply_seconds'] = time.perf_counter()-started
        self.metadata['last_source_shape_pixels'] = [sh, sw]
        if return_ideal:
            return {'camera': camera, 'ideal': ideal if ideal is not None else camera.copy()}
        return camera

    def apply_fast(self, temperature_obs, dilution=1., *, return_ideal=True):
        """Apply the same wavelength-dependent optics with a thermal SVD basis.

        The only new approximation is a tested scalar-spectrum factorization
        and its temperature coefficient lookup. Linearity permits exchanging
        the basis, detector average, Fourier convolution and spectral sums.
        In particular, absolute radiance is averaged BEFORE diffraction, and
        the same finite-support Airy OTFs and zero-padded origin crop are used.
        Signed basis coefficients must never be clipped as if they were light.

        The basis extends from the direct Planck implementation's exact-zero
        temperature (~24.76 K at 830 nm) through 100000 K. Hotter inputs and
        entirely cold scenes use the unchanged direct implementation.
        """
        temperature, amplitude = np.broadcast_arrays(
            np.asarray(temperature_obs, dtype=np.float64),
            np.asarray(dilution, dtype=np.float64))
        if temperature.ndim != 2 or not np.isfinite(temperature).all() or (temperature < 0).any():
            raise ValueError('temperature must be a finite nonnegative 2D image')
        if not np.isfinite(amplitude).all() or (amplitude < 0).any():
            raise ValueError('dilution must be finite and nonnegative')
        sh, sw = temperature.shape
        if sh < self.height or sw < self.width:
            raise ValueError('source map must be at least detector resolution')
        if abs(sw/self.width-sh/self.height) > max(1/self.width, 1/self.height)+1e-12:
            raise ValueError('source and detector aspect ratios must agree to within one source pixel')
        preparation_started = time.perf_counter()
        # Evaluate the actual temperature here: clipping a faint sub-100 K
        # source upward would invent light and trigger an unnecessary fallback.
        evaluation_temperature = np.maximum(temperature, 1.)
        peak_wavelength = np.clip(WIEN_WAVELENGTH_T/evaluation_temperature,
                                  self.wavelengths_m[0], self.wavelengths_m[-1])
        scale = amplitude*planck_lambda(peak_wavelength, evaluation_temperature)
        emitting = scale > 0.
        zero_temperature = H*C/(self.wavelengths_m[-1]*KB*700.)
        out_of_domain = emitting & ((temperature < zero_temperature) | (temperature > 1e5))
        cold_scene = bool(emitting.any() and temperature[emitting].max() < 2000.)
        if out_of_domain.any() or cold_scene:
            self.metadata['acceleration'] = {
                'method': 'direct fallback',
                'reason': 'temperature outside basis domain' if out_of_domain.any() else 'entirely cold scene',
                'optical_model_changed': False,
            }
            return self.apply(temperature, amplitude, return_ideal=return_ideal)
        preparation_seconds = time.perf_counter()-preparation_started
        started = time.perf_counter()
        if not hasattr(self, '_basis_transfers'):
            basis, log_grid, coefficients, report = _thermal_basis(tuple(self.wavelengths_nm))
            self._thermal_basis_data = basis, log_grid, coefficients, report
            self._basis_transfers = []
            # Sum in float64, then cache real float32 OTFs just like apply().
            # Building one component at a time avoids a large wavelength stack.
            for component in basis:
                channels = []
                for channel in range(3):
                    transfer = np.zeros_like(self._transfers[0], dtype=np.float64)
                    for index, value in enumerate(component*self.rgb_weights[:, channel]):
                        transfer += value*self._transfers[index]
                    channels.append(transfer.astype(np.float32))
                self._basis_transfers.append(channels)
            self._basis_setup_seconds = time.perf_counter()-started
        basis, log_grid, coefficients, report = self._thermal_basis_data
        self.metadata['acceleration'] = dict(report)
        self.metadata['acceleration']['setup_seconds'] = self._basis_setup_seconds
        self.metadata['acceleration']['basis_OTF_cache_bytes'] = sum(
            x.nbytes for channels in self._basis_transfers for x in channels)
        started = time.perf_counter()
        safe_temperature = np.clip(temperature, math.exp(log_grid[0]), 1e5)
        log_temperature = np.log(safe_temperature)
        # Uniform log-T grid permits O(1) coefficient interpolation per pixel.
        location = np.clip((log_temperature-log_grid[0])/(log_grid[1]-log_grid[0]),
                           0., len(log_grid)-1.)
        lower = np.minimum(location.astype(np.int32), len(log_grid)-2)
        fraction = location-lower
        needs_resize = temperature.shape != (self.height, self.width)
        ideal = np.zeros((self.height, self.width, 3), dtype=np.float64) if return_ideal else None
        frequency = [np.zeros(self._transfers[0].shape, dtype=np.complex128) for _ in range(3)]
        ideal_weights = basis @ self.rgb_weights
        for index, transfers in enumerate(self._basis_transfers):
            source = scale*(coefficients[lower, index]*(1.-fraction) +
                            coefficients[lower+1, index]*fraction)
            if needs_resize:
                source = np.asarray(Image.fromarray(source.astype('float32')).resize(
                    (self.width, self.height), Image.Resampling.BOX), dtype=np.float64)
            if ideal is not None:
                ideal += source[..., None]*ideal_weights[index]
            transformed = np.fft.rfft2(source, s=self.fft_shape)
            for channel in range(3):
                frequency[channel] += transformed*transfers[channel]
        camera = np.stack([np.fft.irfft2(channel, s=self.fft_shape)[:self.height, :self.width]
                           for channel in frequency], axis=-1)
        self.metadata['last_apply_seconds'] = time.perf_counter()-started+preparation_seconds
        self.metadata['last_source_shape_pixels'] = [sh, sw]
        self.metadata['acceleration']['last_apply_seconds'] = self.metadata['last_apply_seconds']
        if return_ideal:
            return {'camera': camera, 'ideal': ideal}
        return camera


def fast_diagnostics(include_raymaps=True):
    """Compare accelerated optics with the unmodified wavelength loop.

    Errors use absolute linear radiance BEFORE exposure, clipping or display.
    The 270-pixel Kerr maps include the production overscan, 1.5x source
    sampling, atmosphere, redshift and inclination range.
    """
    started = time.perf_counter()
    records = []
    rgb_to_xyz = np.linalg.inv(XYZ_TO_LINEAR_SRGB)
    def compare(camera, temperature, dilution, name, crop=None):
        begin = time.perf_counter()
        direct = camera.apply(temperature, dilution, return_ideal=True)
        direct_seconds = time.perf_counter()-begin
        accelerated = camera.apply_fast(temperature, dilution, return_ideal=True)
        fast_seconds = camera.metadata['last_apply_seconds']
        result = {'name': name, 'direct_seconds': direct_seconds,
                  'fast_seconds_excluding_once_only_setup': fast_seconds,
                  'speedup': direct_seconds/max(fast_seconds, 1e-12), 'views': {}}
        for key in direct:
            original, fast = direct[key], accelerated[key]
            if crop is not None:
                original, fast = original[crop], fast[crop]
            difference = fast-original
            peak = max(float(np.max(np.abs(original))), np.finfo(float).tiny)
            luminance = original @ rgb_to_xyz[1]
            luminance_difference = difference @ rgb_to_xyz[1]
            peak_y = max(float(np.max(np.abs(luminance))), np.finfo(float).tiny)
            result['views'][key] = {
                'max_RGB_error_relative_to_frame_peak': float(np.max(np.abs(difference))/peak),
                'RMS_RGB_error_relative_to_frame_peak': float(np.sqrt(np.mean(difference**2))/peak),
                'max_CIE_Y_error_relative_to_frame_peak_Y': float(np.max(np.abs(luminance_difference))/peak_y),
                'finite': bool(np.isfinite(fast).all()),
            }
        records.append(result)
        return camera.metadata['acceleration']

    camera = SpectralCamera(128, 160, 24.)
    yy, xx = np.mgrid[:240, :192]
    zero_temperature = H*C/(camera.wavelengths_m[-1]*KB*700.)
    temperature = np.geomspace(zero_temperature*(1.+1e-12), 1e5, yy.size).reshape(yy.shape)
    temperature[(xx+3*yy)%17 < 4] = 0.
    dilution = np.random.default_rng(2718).uniform(.1, 1., yy.shape)
    acceleration = compare(camera, temperature, dilution, 'Planck numerical zero threshold--100000 K; zero sentinels; variable dilution; 1.5x source')
    sub_100 = np.full((240, 192), 30.)
    sub_100[:, :96] = 99.
    sub_100[:, :48] = zero_temperature*.999
    sub_100[120:125, 96:102] = 30000.
    compare(camera, sub_100, .7, 'sub-100 K source with hot pixels; must use basis without whole-frame fallback')
    sub_100_uses_basis = camera.metadata['acceleration']['method'].startswith('normalized scalar')
    checker = np.where((xx+yy)%2, 3000., 30000.)
    compare(camera, checker, .7, 'two-temperature source checkerboard; radiance averaging')
    # Isolated hot pixels also exercise dark fields, chromatic diffraction,
    # and the noncommuting temperature average that must never be performed.
    impulses = np.zeros((240, 192))
    impulses[40, 60], impulses[121, 100], impulses[200, 170] = 2000., 15000., 100000.
    compare(camera, impulses, 1., 'three monochromatic-temperature impulses')
    cold = np.full((160, 128), 120.)
    compare(camera, cold, .5, '120 K all-cold direct fallback')
    outside = np.full((160, 128), 200000.)
    compare(camera, outside, .3, '200000 K out-of-domain direct fallback')
    zero = camera.apply_fast(np.full((160, 128), zero_temperature*.999), return_ideal=True)
    zero_exact = all(np.count_nonzero(value) == 0 for value in zero.values())
    if include_raymaps:
        from kerr_emission import KerrDiskParameters, KerrEmission
        from kerr_raytrace import KerrTracer
        from optics import color_correction, limb_factor
        width, hero_h, fov = 270, 405, 24.
        pad = max(2, math.ceil(FIRST_NULL_FACTOR*830e-9/10./(fov/(1e6*width))*24.))
        detector_w, detector_h = width+2*pad, hero_h+2*pad
        camera_fov = fov*detector_w/width
        camera = SpectralCamera(detector_w, detector_h, camera_fov)
        disk = KerrEmission(KerrDiskParameters(mass_solar=1e10, eddington_ratio=.05, spin_a_star=.9))
        tracer = KerrTracer(round(detector_w*1.5), round(detector_h*1.5),
                            fov=camera_fov, roll_deg=14, spin=.9, rout=500, angle_step=.018)
        for angle in (85., 86.5, 88.):
            ray = tracer.trace(angle)
            physical_t = disk.temperature(ray['r'])
            fcol = color_correction(physical_t)
            observed_t = ray['g']*fcol*physical_t
            dilution = limb_factor(ray['mu'])/fcol**4
            acceleration = compare(camera, observed_t, dilution, f'270px Kerr raymap at {angle} degrees',
                                   (slice(pad, pad+hero_h), slice(pad, pad+width)))
    errors = [view for record in records for view in record['views'].values()]
    passed = bool(zero_exact and sub_100_uses_basis and all(v['finite'] and
                  v['max_RGB_error_relative_to_frame_peak'] < 2e-6 and
                  v['max_CIE_Y_error_relative_to_frame_peak_Y'] < 2e-6 for v in errors))
    report = {'acceleration': acceleration, 'comparisons': records,
              'below_direct_Planck_threshold_exact_zero': zero_exact,
              'sub_100_K_sources_use_acceleration': sub_100_uses_basis,
              'acceptance_max_RGB_and_Y_error_relative_to_frame_peak': 2e-6,
              'elapsed_seconds': time.perf_counter()-started, 'passed': passed}
    if not passed:
        raise AssertionError(report)
    return report


def diagnostics():
    """Independent quadrature, kernel, convolution, and subpixel checks."""
    started = time.perf_counter()
    wavelengths, weights, cmf = spectral_grid(10.)
    temp = np.geomspace(2000., 1e5, 101)
    sampled = planck_lambda(wavelengths[None, :]*1e-9, temp[:, None])
    coarse_xyz = sampled @ (weights[:, None]*cmf)
    coarse_rgb = coarse_xyz @ XYZ_TO_LINEAR_SRGB.T
    reference = ThermalSpectrum()
    exact_rgb = reference._direct(temp)
    exact_xyz = exact_rgb @ np.linalg.inv(XYZ_TO_LINEAR_SRGB).T
    relative_rgb = np.max(np.abs(coarse_rgb-exact_rgb), axis=1)/np.max(np.abs(exact_rgb), axis=1)
    relative_y = np.abs(coarse_xyz[:, 1]/exact_xyz[:, 1]-1)
    xy_error = np.abs(coarse_xyz[:, :2]/coarse_xyz.sum(axis=1, keepdims=True) -
                      exact_xyz[:, :2]/exact_xyz.sum(axis=1, keepdims=True))

    # Independent direct J1 quadrature, rather than the shared interpolation.
    lookup_kernel = spectral_airy_kernel(2.3, 24., 3)
    direct_kernel = spectral_airy_kernel(2.3, 24., 3, _use_lookup=False)
    kernel_l1_error = float(np.abs(lookup_kernel-direct_kernel).sum())
    refined = spectral_airy_kernel(2.3, 24., 9)
    pixel_integration_l1_change = float(np.abs(lookup_kernel-refined).sum())

    test = SpectralCamera(128, 128, 26., aperture_m=5.)
    impulse = np.zeros((128, 128)); impulse[64, 64] = 1.
    blue = test._convolve(impulse, 0)
    red = test._convolve(impulse, len(test.wavelengths_nm)-1)
    flux_error = float(max(abs(blue.sum()-1), abs(red.sum()-1)))
    yy, xx = np.mgrid[:128, :128]
    squared_radius = (yy-64)**2+(xx-64)**2
    width_blue = float(math.sqrt((blue*squared_radius).sum()/blue.sum()))
    width_red = float(math.sqrt((red*squared_radius).sum()/red.sum()))
    linearity = float(np.max(np.abs(test._convolve(impulse*3, 0)-blue*3)))
    # Independent conventional full FFT convolution plus origin crop.
    kernel = test.kernels[0]
    radius = kernel.shape[0]//2
    direct = np.fft.irfft2(np.fft.rfft2(impulse, s=test.fft_shape)*
                           np.fft.rfft2(kernel, s=test.fft_shape), s=test.fft_shape)
    direct = direct[radius:radius+128, radius:radius+128]
    fourier_storage_error = float(np.max(np.abs(blue-direct)))
    first_null_ratio_error = float(abs(test.first_nulls_pixels[-1]/test.first_nulls_pixels[0]-830/360))

    # A two-temperature subpixel checkerboard must average radiance, not T.
    super_temperature = np.where(np.indices((256, 256)).sum(axis=0)%2, 3000., 30000.)
    rendered = test.apply(super_temperature, .7, diffraction=False)
    expected = .7*np.mean(planck_lambda(wavelengths[None, :]*1e-9, np.array([3000., 30000.])[:, None]), axis=0)
    expected_rgb = expected @ test.rgb_weights
    subpixel_error = float(np.max(np.abs(rendered-expected_rgb))/np.max(np.abs(expected_rgb)))

    result = {
        'wavelength_samples': len(wavelengths),
        'wavelength_step_nm': 10,
        'quadrature_weight_sum_m': float(weights.sum()),
        'temperature_test_range_k': [2000, 100000],
        'temperature_test_samples': len(temp),
        'max_RGB_relative_to_peak_channel_error_against_1nm': float(relative_rgb.max()),
        'max_CIE_Y_relative_error_against_1nm': float(relative_y.max()),
        'max_CIE_xy_absolute_error_against_1nm': float(xy_error.max()),
        'J1_lookup_kernel_L1_error_against_direct_quadrature': kernel_l1_error,
        'kernel_3_to_9_Gauss_nodes_L1_change_at_2_3pixel_first_null': pixel_integration_l1_change,
        'Airy_first_null_wavelength_ratio_absolute_error': first_null_ratio_error,
        'nonnegative_scalar_convolution': bool(blue.min() >= 0 and red.min() >= 0),
        'impulse_flux_error': flux_error,
        'impulse_RMS_radius_360nm_pixels': width_blue,
        'impulse_RMS_radius_830nm_pixels': width_red,
        'red_impulse_broader_than_blue': width_red > width_blue,
        'convolution_linearity_absolute_error': linearity,
        'float32_Fourier_kernel_max_error_against_float64': fourier_storage_error,
        'supersampled_radiance_averaging_relative_error': subpixel_error,
        'elapsed_seconds': time.perf_counter()-started,
    }
    result['passed'] = bool(relative_rgb.max() < .003 and relative_y.max() < .003 and
                            xy_error.max() < .001 and kernel_l1_error < 1e-5 and
                            pixel_integration_l1_change < .001 and flux_error < 2e-6 and
                            first_null_ratio_error < 1e-12 and width_red > width_blue and
                            linearity < 1e-12 and fourier_storage_error < 1e-7 and
                            subpixel_error < 2e-7)
    if not result['passed']:
        raise AssertionError(result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--fast-diagnostics', action='store_true',
                        help='validate the optional thermal-basis acceleration including Kerr raymaps')
    args = parser.parse_args()
    report = fast_diagnostics() if args.fast_diagnostics else diagnostics()
    if args.output:
        args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))
