"""Tests for enhydra.plotting."""

import re
import pytest

from enhydra.plotting import _rdbu_color


def _rgb(color_str: str) -> tuple[int, int, int]:
    m = re.match(r"rgb\((\d+),(\d+),(\d+)\)", color_str)
    assert m, "not an rgb(...) string: %r" % color_str
    return tuple(int(g) for g in m.groups())


def _distance_from_white(color_str: str) -> float:
    """Euclidean distance from MID (247, 247, 247) — a proxy for
    'how saturated/far from white this color is', used to check that
    saturation increases monotonically with |t| regardless of sign."""
    r, g, b = _rgb(color_str)
    return ((r - 247) ** 2 + (g - 247) ** 2 + (b - 247) ** 2) ** 0.5


class TestRdbuColor:

    def test_zero_is_white(self):
        assert _rgb(_rdbu_color(0.0)) == (247, 247, 247)

    def test_full_positive_is_full_blue(self):
        assert _rgb(_rdbu_color(1.0)) == (33, 102, 172)

    def test_full_negative_is_full_red(self):
        assert _rgb(_rdbu_color(-1.0)) == (214, 96, 77)

    def test_clamped_above_one(self):
        assert _rgb(_rdbu_color(5.0)) == (33, 102, 172)

    def test_clamped_below_negative_one(self):
        assert _rgb(_rdbu_color(-5.0)) == (214, 96, 77)

    def test_red_saturation_increases_with_magnitude(self):
        """Regression test for the reported bug: on the red (negative)
        side, a point further from zero must be MORE saturated (further
        from white), not less. A prior version of _rdbu_color() had this
        backwards — t=0 rendered as full red and t=-1 rendered as white,
        exactly inverted."""
        near_zero = _distance_from_white(_rdbu_color(-0.1))
        far       = _distance_from_white(_rdbu_color(-0.9))
        assert far > near_zero

    def test_blue_saturation_increases_with_magnitude(self):
        near_zero = _distance_from_white(_rdbu_color(0.1))
        far       = _distance_from_white(_rdbu_color(0.9))
        assert far > near_zero

    def test_red_and_blue_saturate_symmetrically(self):
        """Equal-magnitude positive and negative t should be equally far
        from white — the two branches should behave as mirror images of
        each other in terms of saturation (not colour, obviously)."""
        for magnitude in (0.1, 0.3, 0.5, 0.7, 0.9):
            red_dist  = _distance_from_white(_rdbu_color(-magnitude))
            blue_dist = _distance_from_white(_rdbu_color(magnitude))
            assert red_dist == pytest.approx(blue_dist, abs=1.0)

    def test_monotonic_along_negative_range(self):
        """Distance from white must increase monotonically as t moves
        from 0 toward -1 — no fold-back partway through the range."""
        ts = [-0.0, -0.2, -0.4, -0.6, -0.8, -1.0]
        distances = [_distance_from_white(_rdbu_color(t)) for t in ts]
        assert distances == sorted(distances)

    def test_monotonic_along_positive_range(self):
        ts = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
        distances = [_distance_from_white(_rdbu_color(t)) for t in ts]
        assert distances == sorted(distances)
