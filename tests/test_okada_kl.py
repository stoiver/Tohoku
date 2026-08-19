"""Tests for the Karhunen-Loeve slip field in okada_kl_subfaults.py.

Two of these guard bugs that were live in the repository:

* `kl_correlation_matrices` used `np.linalg.eig`, the general non-symmetric
  LAPACK solver, on a symmetric covariance matrix.  That may return
  complex-conjugate eigenpairs, which propagate through the slip field into
  okada() and crash the run.  The trigger is the numpy version, not the
  platform: a CI probe with `eig` restored failed on numpy 2.5.2 under both
  OpenBLAS and Accelerate, and passed on numpy 2.4.6 under both.  So
  `test_eigendecomposition_is_real` needs a matrix spanning numpy versions to
  have teeth -- which is what the Python-version spread in CI buys.
* `sample='sobol'` fed raw Sobol points, uniform on [0, 1), into an expansion
  that wants standard normal deviates -- biasing every mode's coefficient by
  +0.5 and truncating the tails.
"""
import numpy as np
import pytest
from numpy import linalg as LA

import okada_kl_subfaults as okl


KM = 1000.0
FAULT = dict(length=200 * KM, width=50 * KM, slip=60.0)


@pytest.fixture(scope='module')
def geometry():
    return okl.subfaults(E_subfault=10, N_subfault=10, dip=14.0, strike=195.0,
                         length=FAULT['length'], width=FAULT['width'])


@pytest.fixture(scope='module')
def matrices(geometry):
    E, N, D = geometry
    return okl.kl_correlation_matrices(E, N, D, FAULT['length'],
                                       FAULT['width'], FAULT['slip'])


def test_covariance_is_symmetric_and_psd(matrices):
    *_, C = matrices
    assert np.array_equal(C, C.T), 'covariance is not exactly symmetric'
    assert LA.eigvalsh(C).min() > -1e-8 * LA.eigvalsh(C).max(), 'not PSD'


def test_eigendecomposition_is_real(matrices):
    """The macOS crash: complex eigenpairs from a symmetric matrix."""
    _, _, _, D, V, sqrtD, _ = matrices
    for name, a in (('D', D), ('V', V), ('sqrtD', sqrtD)):
        assert np.isrealobj(a), f'{name} came back complex'
        assert np.isfinite(a).all(), f'{name} contains NaN or inf'


def test_eigenvalues_are_non_negative(matrices):
    """sqrtD is NaN if any eigenvalue slips below zero."""
    _, _, _, D, *_ = matrices
    assert np.diag(D).min() >= 0.0


def test_factorisation_reproduces_the_covariance(matrices):
    """V sqrtD (V sqrtD)^T must equal the covariance it came from -- the whole
    point of the decomposition."""
    _, _, _, _, V, sqrtD, C = matrices
    L = V @ sqrtD
    np.testing.assert_allclose(L @ L.T, C, rtol=1e-8, atol=1e-8 * abs(C).max())


def test_eigenvectors_are_orthonormal(matrices):
    _, _, _, _, V, _, _ = matrices
    np.testing.assert_allclose(V.T @ V, np.eye(V.shape[0]), atol=1e-10)


@pytest.mark.parametrize('sample', ['random', 'sobol'])
def test_slip_field_is_real_and_finite(geometry, sample):
    E, N, D = geometry
    s, *_ = okl.kl_slipfield(E, N, D, FAULT['length'], FAULT['width'],
                             FAULT['slip'], sample=sample, iseed=1234)
    assert np.isrealobj(s), 'slip field came back complex'
    assert np.isfinite(s).all()
    assert s.shape == (10, 10)


def test_sobol_deviates_are_standard_normal(geometry):
    """The uniform-deviate bug: Sobol points are U[0,1), so a buggy z has no
    negative values and a mean near 0.5 rather than 0."""
    E, N, D = geometry
    _, _, _, z, _ = okl.kl_slipfield(E, N, D, FAULT['length'], FAULT['width'],
                                     FAULT['slip'], sample='sobol', iseed=1234)
    assert (z < 0).any(), 'no negative deviates -- z looks uniform, not normal'
    assert abs(z.mean()) < 0.3, f'z mean {z.mean():.3f} is not near 0'
    assert 0.7 < z.std() < 1.3, f'z std {z.std():.3f} is not near 1'


def test_slip_field_mean_matches_mu(geometry):
    """Averaged over draws the field must centre on the prescribed slip."""
    E, N, D = geometry
    means = [okl.kl_slipfield(E, N, D, FAULT['length'], FAULT['width'],
                              FAULT['slip'], sample='random', iseed=seed)[0].mean()
             for seed in range(40)]
    assert abs(np.mean(means) - FAULT['slip']) < 0.25 * FAULT['slip']


def test_seeded_draws_are_reproducible(geometry):
    E, N, D = geometry
    kw = dict(sample='random', iseed=7)
    a = okl.kl_slipfield(E, N, D, FAULT['length'], FAULT['width'], FAULT['slip'], **kw)[0]
    b = okl.kl_slipfield(E, N, D, FAULT['length'], FAULT['width'], FAULT['slip'], **kw)[0]
    np.testing.assert_array_equal(a, b)


def test_kl_deformation_end_to_end():
    """The path the notebooks actually call."""
    a = np.linspace(-300 * KM, 300 * KM, 25)
    x, y = np.meshgrid(a, a)
    uE, uN, uZ, slips = okl.kl_deformation(
        x, y, E_subfault=10, N_subfault=10, sample='sobol', iseed=1234,
        depth=20 * KM, length=200 * KM, width=50 * KM,
        strike=195.0, dip=14.0, rake=87.0, slip=60.0, opening=0.0, nu=0.25)
    for name, a_ in (('uE', uE), ('uN', uN), ('uZ', uZ), ('slips', slips)):
        assert np.isrealobj(a_), f'{name} came back complex'
        assert np.isfinite(a_).all(), f'{name} contains NaN or inf'
    assert uZ.max() > 0.0 and uZ.min() < 0.0
