"""Tests for the Part 3 response analysis."""

from __future__ import annotations

import unittest

from analyze_response import benjamini_hochberg, calculate_statistics
from load_data import POPULATIONS


class ResponseAnalysisTests(unittest.TestCase):
    def test_benjamini_hochberg_adjustment(self) -> None:
        adjusted = benjamini_hochberg([0.01, 0.04, 0.03, 0.002])

        expected = [0.02, 0.04, 0.04, 0.008]
        for actual, target in zip(adjusted, expected):
            self.assertAlmostEqual(actual, target)

    def test_statistics_include_all_populations_and_group_sizes(self) -> None:
        groups = {
            population: {
                "yes": [20.0, 21.0, 22.0, 23.0],
                "no": [5.0, 6.0, 7.0, 8.0],
            }
            for population in POPULATIONS
        }

        results = calculate_statistics(groups)

        self.assertEqual(
            [result["population"] for result in results], list(POPULATIONS)
        )
        self.assertTrue(
            all(result["responder_subjects"] == 4 for result in results)
        )
        self.assertTrue(
            all(result["non_responder_subjects"] == 4 for result in results)
        )
        self.assertTrue(
            all(result["median_difference"] == 15.0 for result in results)
        )

    def test_statistics_require_both_response_groups(self) -> None:
        groups = {
            population: {"yes": [10.0], "no": [9.0]}
            for population in POPULATIONS
        }
        groups["b_cell"]["no"] = []

        with self.assertRaisesRegex(ValueError, "Both response groups"):
            calculate_statistics(groups)


if __name__ == "__main__":
    unittest.main()
