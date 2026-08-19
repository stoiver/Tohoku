"""Tests for the Okada (1985) dislocation model in okada.py.

These need only numpy/scipy/matplotlib -- no ANUGA -- so they run in CI on a
bare Python environment.
"""
import numpy as np
import pytest

import okada


KM = 1000.0

# A thrust fault roughly like the Tohoku source in the notebooks.
THRUST = dict(depth=20 * KM, length=200 * KM, width=50 * KM,
              strike=195.0, dip=14.0, rake=90.0, slip=10.0, opening=0.0,
              nu=0.25)


def grid(half_width=400 * KM, n=41):
    """Observation grid centred on the fault."""
    a = np.linspace(-half_width, half_width, n)
    return np.meshgrid(a, a)


def test_zero_slip_gives_zero_deformation():
    x, y = grid()
    p = dict(THRUST, slip=0.0, opening=0.0)
    ue, un, uz = okada.forward(x, y, **p)
    assert np.allclose(ue, 0.0)
    assert np.allclose(un, 0.0)
    assert np.allclose(uz, 0.0)


def test_deformation_is_real_and_finite():
    """A complex or NaN displacement field is the failure that crashed a macOS
    run: it propagates out of the KL slip field and into set_quantity."""
    x, y = grid()
    for component in okada.forward(x, y, **THRUST):
        assert np.isrealobj(component), 'displacement came back complex'
        assert np.isfinite(component).all(), 'displacement contains NaN or inf'


def test_displacement_scales_linearly_with_slip():
    """Okada is a linear elastic solution: doubling slip doubles displacement."""
    x, y = grid()
    u1 = okada.forward(x, y, **dict(THRUST, slip=1.0))
    u7 = okada.forward(x, y, **dict(THRUST, slip=7.0))
    for a, b in zip(u1, u7):
        np.testing.assert_allclose(7.0 * a, b, rtol=1e-10, atol=1e-12)


def test_thrust_uplifts_and_subsides():
    """A dipping thrust lifts the hanging wall and drops the footwall, so the
    vertical field must change sign and be bounded by the slip."""
    x, y = grid()
    _, _, uz = okada.forward(x, y, **THRUST)
    assert uz.max() > 0.0, 'no uplift anywhere'
    assert uz.min() < 0.0, 'no subsidence anywhere'
    assert uz.max() < THRUST['slip'], 'uplift exceeds the slip on the fault'


def test_deformation_decays_with_distance():
    """Displacement must fall off away from the fault."""
    near_x, near_y = np.array([[0.0]]), np.array([[0.0]])
    far_x, far_y = np.array([[5000 * KM]]), np.array([[5000 * KM]])
    _, _, uz_near = okada.forward(near_x, near_y, **THRUST)
    _, _, uz_far = okada.forward(far_x, far_y, **THRUST)
    assert abs(uz_far).max() < 1e-3 * abs(uz_near).max()


def test_offset_shifts_the_field():
    """xoff/yoff translate the source rigidly."""
    x, y = grid()
    _, _, uz = okada.forward(x, y, **THRUST)
    _, _, uz_shifted = okada.forward(x + 30 * KM, y - 20 * KM,
                                     xoff=30 * KM, yoff=-20 * KM, **THRUST)
    np.testing.assert_allclose(uz, uz_shifted, rtol=1e-10, atol=1e-12)


@pytest.mark.parametrize('strike', [0.0, 90.0, 180.0, 270.0])
def test_rotating_strike_rotates_the_field(strike):
    """Rotating the fault by its strike rotates the horizontal field with it,
    leaving the vertical field's extremes unchanged."""
    x, y = grid()
    _, _, uz_ref = okada.forward(x, y, **dict(THRUST, strike=0.0))
    _, _, uz = okada.forward(x, y, **dict(THRUST, strike=strike))
    np.testing.assert_allclose(sorted([uz.min(), uz.max()]),
                               sorted([uz_ref.min(), uz_ref.max()]),
                               rtol=1e-6)


def test_subfaults_sum_matches_a_single_fault():
    """okada_subfaults with one subfault must reproduce okada.forward."""
    import okada_subfaults

    x, y = grid(n=21)
    ue, un, uz = okada.forward(x, y, **THRUST)
    ue_s, un_s, uz_s = okada_subfaults.forward(
        x, y, E_subfault=1, N_subfault=1, **THRUST)
    np.testing.assert_allclose(uz_s, uz, rtol=1e-8, atol=1e-9)
    np.testing.assert_allclose(ue_s, ue, rtol=1e-8, atol=1e-9)
    np.testing.assert_allclose(un_s, un, rtol=1e-8, atol=1e-9)
