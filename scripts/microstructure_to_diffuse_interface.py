#!/usr/bin/env python3
"""
Convert a grayscale/backscatter micrograph into a binary (or diffuse-interface)
image suitable for use as an inclusion field in Alamo's diffuse-interface
codes. The dark matrix becomes black (0), particle inclusions become white
(1), and an optional Gaussian blur produces a smooth (diffuse) interface
between the two phases.

Usage:
    python3 microstructure_to_diffuse_interface.py input.png output.png \
        --threshold 128 --blur 2.0

    # Auto-threshold (Otsu) instead of a fixed value:
    python3 microstructure_to_diffuse_interface.py input.png output.png \
        --auto-threshold --blur 3.0
"""
import argparse
from typing import Optional

import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter


def otsu_threshold(gray: np.ndarray) -> float:
    hist, bin_edges = np.histogram(gray, bins=256, range=(0, 255))
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    weight1 = np.cumsum(hist)
    weight2 = np.cumsum(hist[::-1])[::-1]

    mean1 = np.cumsum(hist * bin_centers) / np.where(weight1 == 0, 1, weight1)
    mean2 = (np.cumsum((hist * bin_centers)[::-1]) / np.where(weight2[::-1] == 0, 1, weight2[::-1]))[::-1]

    variance12 = weight1[:-1] * weight2[1:] * (mean1[:-1] - mean2[1:]) ** 2
    idx = np.argmax(variance12)
    return bin_centers[idx]


def process(
    input_path: str,
    output_path: str,
    threshold: Optional[float] = None,
    auto_threshold: bool = False,
    blur: float = 0.0,
    blur_x: Optional[float] = None,
    blur_y: Optional[float] = None,
    invert: bool = False,
) -> None:
    img = Image.open(input_path).convert("L")
    gray = np.asarray(img, dtype=np.float64)

    if auto_threshold:
        threshold = otsu_threshold(gray)
        print(f"Auto (Otsu) threshold: {threshold:.1f}")
    elif threshold is None:
        threshold = 128.0

    # Particles (bright) -> 1 (white), matrix (dark) -> 0 (black)
    binary = (gray >= threshold).astype(np.float64)
    if invert:
        binary = 1.0 - binary

    # sigma is (row, col) = (y, x); blur_x/blur_y override the isotropic
    # --blur value, needed when the image is stretched non-uniformly onto
    # the physical domain (fit=coord with different x/y um-per-pixel).
    sigma_x = blur_x if blur_x is not None else blur
    sigma_y = blur_y if blur_y is not None else blur
    if sigma_x > 0 or sigma_y > 0:
        field = gaussian_filter(binary, sigma=(sigma_y, sigma_x))
    else:
        field = binary

    out_img = Image.fromarray((np.clip(field, 0, 1) * 255).astype(np.uint8))
    out_img.save(output_path)
    print(f"Saved {output_path}  (shape={field.shape}, blur sigma={blur})")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", help="Path to the input micrograph image")
    parser.add_argument("output", help="Path to write the processed image")
    parser.add_argument(
        "--threshold", type=float, default=None,
        help="Fixed grayscale threshold (0-255) separating matrix from particles. Default 128.",
    )
    parser.add_argument(
        "--auto-threshold", action="store_true",
        help="Use Otsu's method to pick the threshold automatically (overrides --threshold).",
    )
    parser.add_argument(
        "--blur", type=float, default=0.0,
        help="Gaussian blur sigma (in pixels) applied after thresholding to create a diffuse interface. Default 0 (sharp interface).",
    )
    parser.add_argument(
        "--blur-x", type=float, default=None,
        help="Blur sigma in pixels along x (columns), overriding --blur. Use with --blur-y when the image is stretched anisotropically onto the physical domain.",
    )
    parser.add_argument(
        "--blur-y", type=float, default=None,
        help="Blur sigma in pixels along y (rows), overriding --blur.",
    )
    parser.add_argument(
        "--invert", action="store_true",
        help="Invert the result if particles come out black instead of white.",
    )
    args = parser.parse_args()

    process(
        args.input,
        args.output,
        threshold=args.threshold,
        auto_threshold=args.auto_threshold,
        blur=args.blur,
        blur_x=args.blur_x,
        blur_y=args.blur_y,
        invert=args.invert,
    )


if __name__ == "__main__":
    main()
