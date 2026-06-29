#!/usr/bin/env python3
from __future__ import annotations

import os
import stat
import tempfile
import unittest
from pathlib import Path

import testlib_gpu


class FindBinaryTest(unittest.TestCase):
    def test_cuda_arch_glob_selects_highest_matching_executable(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bin_dir = root / "bin"
            bin_dir.mkdir()

            older = bin_dir / "alamo_gpu-2d-cuda80-g++"
            newer = bin_dir / "alamo_gpu-2d-cuda100-g++"
            ignored = bin_dir / "alamo_gpu-2d-cuda90-g++"
            for path in (older, newer, ignored):
                path.write_text("#!/bin/sh\n", encoding="utf-8")
            older.chmod(older.stat().st_mode | stat.S_IXUSR)
            newer.chmod(newer.stat().st_mode | stat.S_IXUSR)
            ignored.chmod(ignored.stat().st_mode | stat.S_IXUSR)

            self.assertEqual(
                testlib_gpu.find_binary(root, "bin/alamo_gpu-2d-cuda*-g++"),
                newer,
            )

    def test_find_binary_ignores_non_executable_matches(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bin_dir = root / "bin"
            bin_dir.mkdir()

            candidate = bin_dir / "alamo_gpu-3d-cuda90-g++"
            candidate.write_text("#!/bin/sh\n", encoding="utf-8")
            candidate.chmod(candidate.stat().st_mode & ~stat.S_IXUSR & ~stat.S_IXGRP & ~stat.S_IXOTH)

            self.assertIsNone(testlib_gpu.find_binary(root, "bin/alamo_gpu-3d-cuda*-g++"))


if __name__ == "__main__":
    unittest.main()
