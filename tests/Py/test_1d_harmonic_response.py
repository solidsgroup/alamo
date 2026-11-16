import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.analyze_1d_harmonic_response import fit_pair_response


class HarmonicResponseTest(unittest.TestCase):
    def test_recovers_complex_response_and_control_correction(self):
        frequency_hz = 17.0
        omega = 2.0 * np.pi * frequency_hz
        pressure_amplitude = 0.02
        mean_pressure = 2.0
        pressure_fraction = pressure_amplitude / mean_pressure
        forcing_start = 0.13
        phase = -0.7
        trim = 0.37
        time = np.arange(0.0, 3.6, 1.0 / (frequency_hz * 30.0))
        selected = time >= trim
        time = time[selected]
        center = float(np.mean(time))

        def front_history(mean_rate, response):
            rate_phasor = mean_rate * pressure_fraction * response
            front_phasor = 1j * rate_phasor / omega
            theta = omega * (time - forcing_start) + phase
            drift = -mean_rate * (time - center) + 1.0e-8 * (time - center) ** 2
            return 0.01 + drift + np.real(front_phasor * np.exp(1j * theta))

        drive_response = 0.30 + 0.20j
        control_response = 0.04 - 0.03j
        result = fit_pair_response(
            time,
            front_history(4.0e-4, drive_response),
            front_history(4.0e-4, control_response),
            frequency_hz,
            pressure_amplitude,
            mean_pressure,
            forcing_start,
            phase,
        )

        self.assertAlmostEqual(result["driven"]["response"].real, drive_response.real, places=9)
        self.assertAlmostEqual(result["driven"]["response"].imag, drive_response.imag, places=9)
        self.assertAlmostEqual(result["control"]["response"].real, control_response.real, places=9)
        self.assertAlmostEqual(result["control"]["response"].imag, control_response.imag, places=9)
        expected_corrected = drive_response - control_response
        self.assertAlmostEqual(result["corrected"].real, expected_corrected.real, places=9)
        self.assertAlmostEqual(result["corrected"].imag, expected_corrected.imag, places=9)


if __name__ == "__main__":
    unittest.main()
