"""Smoke tests for the interactive Streamlit dashboard."""

from __future__ import annotations

import unittest

from streamlit.testing.v1 import AppTest


class DashboardTests(unittest.TestCase):
    def test_dashboard_renders_core_sections_and_metrics(self) -> None:
        app = AppTest.from_file("dashboard.py", default_timeout=30)

        app.run()

        self.assertEqual(len(app.exception), 0)
        self.assertEqual(
            [tab.label for tab in app.tabs],
            [
                "Trial overview",
                "Cell frequencies",
                "Miraclib response",
                "Baseline subset",
            ],
        )
        metrics = {metric.label: metric.value for metric in app.metric}
        self.assertEqual(metrics["Projects"], "3")
        self.assertEqual(metrics["Samples"], "10500")
        self.assertEqual(metrics["Baseline samples"], "656")
        self.assertEqual(
            metrics["Male responder B-cell average"], "10206.15"
        )


if __name__ == "__main__":
    unittest.main()
