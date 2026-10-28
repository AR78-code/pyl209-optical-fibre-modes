"""Ghatak, Optical Electronics, §§13.5–13.6: scalar step-index modes.

All internal lengths are micrometres, beta is rad/µm. Time convention:
Psi = psi(r,phi) exp[i(omega*t-beta*z)], Eq. (13.15).
Only one fixed polarization is propagated. No fitted/artificial modal weights.
"""
from dataclasses import dataclass, asdict
from functools import cached_property
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import jv, jvp, jn_zeros, kve
from numpy.polynomial.legendre import leggauss


@dataclass(frozen=True)
class Fibre:
    wavelength: float = 0.650
    a: float = 6.25
    clad_radius: float = 62.5
    n1: float = 1.50
    n2: float = 1.47

    @property
    def k0(self):
        return 2 * np.pi / self.wavelength

    @property
    def NA(self):
        return np.sqrt(self.n1**2 - self.n2**2)

    @property
    def V(self):
        return self.k0 * self.a * self.NA


def k_ratio(ell, x, denominator):
    """K_ell(x)/K_ell(denominator), evaluated without exponential underflow."""
    return kve(ell, x) / kve(ell, denominator) * np.exp(denominator - x)


def k_log_derivative(ell, w):
    return -0.5 * (kve(abs(ell - 1), w) + kve(ell + 1, w)) / kve(ell, w)


def characteristic(u, ell, V):
    """Eq. (13.36), cross multiplied to remove J denominator poles."""
    w = np.sqrt((V - u) * (V + u))
    return u * jvp(ell, u) - w * k_log_derivative(ell, w) * jv(ell, u)


@dataclass
class Mode:
    fibre: Fibre
    ell: int
    m: int
    u: float
    cutoff: float

    @property
    def label(self):
        return f"LP({self.ell},{self.m})"

    @property
    def w(self):
        return np.sqrt(self.fibre.V**2 - self.u**2)

    @property
    def beta(self):
        return np.sqrt((self.fibre.n1 * self.fibre.k0)**2 - (self.u / self.fibre.a)**2)

    @property
    def neff(self):
        return self.beta / self.fibre.k0

    @property
    def b(self):
        return self.w**2 / self.fibre.V**2

    def raw_radial(self, r):
        """Equivalent to Eq. (13.35), but core coefficient is 1 before normalization."""
        r = np.asarray(r, dtype=float)
        core = r <= self.fibre.a
        result = np.empty_like(r)
        result[core] = jv(self.ell, self.u * r[core] / self.fibre.a)
        result[~core] = jv(self.ell, self.u) * k_ratio(
            self.ell, self.w * r[~core] / self.fibre.a, self.w)
        return result

    @cached_property
    def radial_integral(self):
        a = self.fibre.a
        inside = quad(lambda r: float(self.raw_radial(r))**2 * r, 0, a,
                      epsabs=1e-11, epsrel=1e-11)[0]
        outside = quad(lambda r: float(self.raw_radial(r))**2 * r, a, np.inf,
                       epsabs=1e-11, epsrel=1e-11)[0]
        return inside + outside

    @property
    def norm(self):
        # Integral of 1 over phi is 2*pi; cos²/sin² gives pi for ell > 0.
        return np.sqrt((2*np.pi if self.ell == 0 else np.pi) * self.radial_integral)

    def radial(self, r):
        return self.raw_radial(r) / self.norm

    def field(self, r, phi, orientation="cos"):
        angular = np.cos(self.ell*np.asarray(phi))
        if orientation == "sin":
            if self.ell == 0:
                raise ValueError("ell=0 has no independent sine mode")
            angular = np.sin(self.ell*np.asarray(phi))
        return self.radial(r) * angular

    def record(self):
        return dict(label=self.label, ell=self.ell, m=self.m, U=self.u, W=self.w,
                    beta_rad_per_um=self.beta, neff=self.neff, b=self.b,
                    cutoff_V=self.cutoff, spatial_degeneracy=1 if self.ell == 0 else 2)


def solve_modes(fibre=Fibre(), xtol=2e-12):
    """Bracket each root using the Bessel cutoffs in Table 13.1.

    LP(0,1): 0 < U < j_(0,1).
    LP(0,m>1): j_(1,m-1) < U < j_(0,m).
    LP(ell>=1,m): j_(ell-1,m) < U < j_(ell,m).
    Upper endpoint is also bounded by V. No sign-change search across poles.
    """
    V = fibre.V
    if V <= 0:
        raise ValueError("Positive V required")
    count = int(np.ceil(V / np.pi)) + 3
    modes = []
    for ell in range(int(np.ceil(V)) + 1):
        upper = jn_zeros(ell, count)
        lower = np.r_[0., jn_zeros(1, count - 1)] if ell == 0 else jn_zeros(ell-1, count)
        for m, (lo, hi) in enumerate(zip(lower, upper), 1):
            if lo >= V:
                break
            left = max(lo, 1e-10)
            right = min(hi, np.nextafter(V, 0))
            root = brentq(characteristic, left, right, args=(ell, V),
                          xtol=xtol, rtol=1e-14)
            modes.append(Mode(fibre, ell, m, root, lo))
    return modes


def quadrature(start, end, n):
    x, w = leggauss(n)
    return start + (x+1)*(end-start)/2, w*(end-start)/2


class Simulation:
    """Centred launch in the supplied geometry; overlap coefficients computed numerically.

    Lens diameter D=20 mm; q=D/f changes f with the entrance kept at focus.
    Gaussian amplitude exp[-rho²/wL²], wL/(D/2)=0.25 by default.
    Lens-passed scalar power is normalized to 1. Scalar Fresnel focusing.
    """
    def __init__(self, fibre=Fibre(), fill=0.25, nrad=240, npupil=200):
        self.fibre = fibre
        self.fill = fill
        self.modes = solve_modes(fibre)
        self.active = [m for m in self.modes if m.ell == 0]
        self.D = 20000.0
        self.nrad = nrad
        self.npupil = npupil
        # Truncate overlap integration only after at least 28 evanescent lengths.
        tail = max(28*fibre.a / m.w for m in self.active)
        self.rmax = max(4*fibre.a, fibre.a + tail)
        r1, wr1 = quadrature(0, fibre.a, nrad)
        r2, wr2 = quadrature(fibre.a, self.rmax, nrad)
        self.r = np.r_[r1, r2]
        self.weights = 2*np.pi*self.r*np.r_[wr1, wr2]
        self.radials = np.array([m.radial(self.r) for m in self.active])
        self.pupil_t, self.pupil_weights = quadrature(0, 1, npupil)

    def focused_field(self, r, ratio):
        if not 0.1 <= ratio <= 1.0:
            raise ValueError("D/f must lie between 0.1 and 1.0 for this demonstration")
        r = np.asarray(r)
        R = self.D / 2
        f = self.D / ratio
        wL = self.fill * R
        norm = np.sqrt(np.pi*wL*wL/2 * (-np.expm1(-2/self.fill**2)))
        pupil = self.pupil_t * np.exp(-(self.pupil_t/self.fill)**2) * self.pupil_weights
        # Angular integration of exp(+i k r rho cos(theta)/f) is 2*pi*J0.
        field = (2*np.pi*R*R/(self.fibre.wavelength*f)/norm) * (
            jv(0, self.fibre.k0 * r.reshape(-1, 1) * R/f * self.pupil_t) @ pupil)
        # Keep the focal-plane quadratic phase; drop only a global prefactor phase.
        phase = np.exp(-1j*self.fibre.k0*r.reshape(-1)**2/(2*f))
        return (field * phase).reshape(r.shape)

    def coefficients(self, ratio):
        return self.radials @ (self.focused_field(self.r, ratio) * self.weights)

    def output(self, r, length_mm, ratio, coefficients=None):
        c = self.coefficients(ratio) if coefficients is None else coefficients
        beta = np.array([m.beta for m in self.active])
        # Remove common carrier phase: intensities and relative phases unchanged.
        phase = np.exp(-1j*(beta-beta[0])*length_mm*1000)
        basis = np.array([m.radial(r) for m in self.active])
        return np.einsum("j,j...->...", c*phase, basis)

    def summary(self):
        spatial = sum(1 if m.ell == 0 else 2 for m in self.modes)
        return dict(parameters=asdict(self.fibre), NA=self.fibre.NA, V=self.fibre.V,
                    LP_families=len(self.modes), spatial_modes_one_polarization=spatial,
                    modes_both_polarizations=2*spatial,
                    centred_launch_families=len(self.active),
                    gaussian_radius_at_lens_um=self.fill*self.D/2)
