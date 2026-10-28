"""Physics/numerical checks, including an independent finite-volume eigenproblem."""
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.special import jv, jvp
from fibre import Fibre, Simulation, solve_modes, quadrature, k_log_derivative


def run_validation(sim=None):
    s = sim or Simulation()
    f = s.fibre
    checks = {}
    def record(name, value, limit):
        checks[name] = {"value": float(value), "limit": float(limit), "passed": bool(value < limit)}

    residuals = []
    for m in s.modes:
        lhs = m.u*jvp(m.ell, m.u)
        rhs = m.w*k_log_derivative(m.ell, m.w)*jv(m.ell, m.u)
        residuals.append(abs(lhs-rhs)/max(1, abs(lhs), abs(rhs)))
    record("maximum_normalized_derivative_boundary_residual", max(residuals), 1e-9)
    record("invalid_guided_effective_indices", sum(not f.n2 < m.neff < f.n1 for m in s.modes), .5)

    gram_error = 0
    tail_boundary = 0
    for ell in sorted(set(m.ell for m in s.modes)):
        modes = [m for m in s.modes if m.ell == ell]
        rmax = f.a + max(30*f.a/m.w for m in modes)
        r1, w1 = quadrature(0, f.a, 260)
        r2, w2 = quadrature(f.a, rmax, 300)
        r = np.r_[r1, r2]
        weights = r*np.r_[w1,w2]*(2*np.pi if ell == 0 else np.pi)
        basis = np.array([m.radial(r) for m in modes])
        gram = (basis*weights) @ basis.T
        gram_error = max(gram_error, np.max(abs(gram-np.eye(len(modes)))))
        for m in modes:
            tail_boundary = max(tail_boundary, abs(float(m.radial(f.clad_radius))) /
                                max(abs(m.radial(np.linspace(0,f.a,1000)))))
    record("orthonormality_max_error", gram_error, 2e-8)
    record("field_at_outer_cladding_over_core_peak", tail_boundary, 1e-6)

    refined = Simulation(nrad=400, npupil=320)
    convergence = max(np.max(abs(s.coefficients(q)-refined.coefficients(q)))
                      for q in [.1,.27,.5,.83,1.])
    record("overlap_quadrature_refinement_max_amplitude_change", convergence, 1e-7)

    gaussian_error=0
    for q in [.1,.3,1.]:
        focal=s.D/q
        waist=f.wavelength*focal/(np.pi*s.fill*s.D/2)
        rcheck=np.linspace(0,4*waist,500)
        analytic=np.sqrt(2/np.pi)/waist*np.exp(-(rcheck/waist)**2)*np.exp(-1j*f.k0*rcheck**2/(2*focal))
        gaussian_error=max(gaussian_error,np.max(abs(s.focused_field(rcheck,q)-analytic))/abs(analytic[0]))
    record("diffraction_vs_untruncated_gaussian_relative_peak_error",gaussian_error,2e-7)

    # The offline demo interpolates 0.01-um radial samples only for plotting.
    # Compare midpoint interpolation against direct Bessel evaluation for all modes.
    grid=np.linspace(0,18,1801);mid=(grid[:-1]+grid[1:])/2
    interp_error=max(np.max(abs((m.radial(grid[:-1])+m.radial(grid[1:]))/2-m.radial(mid))) /
                     np.max(abs(m.radial(grid))) for m in s.modes)
    record("browser_radial_interpolation_relative_peak_error",interp_error,1e-4)

    c = s.coefficients(1.)
    expected = np.sum(abs(c)**2)
    power_errors = []
    for length in [0, .05, .37, 1, 10, 100, 1000]:
        power = np.sum(abs(s.output(s.r, length, 1., c))**2 * s.weights)
        power_errors.append(abs(power-expected))
    record("guided_power_conservation_max_error", max(power_errors), 2e-8)
    record("guided_power_above_incident", max(0., max(np.sum(abs(s.coefficients(q))**2)
           for q in np.linspace(.1,1,13))-1.), 1e-7)

    # Angular overlap of a symmetric launch with all ell>0 is zero.
    phi = np.linspace(0,2*np.pi,256,endpoint=False)
    symmetry = max(abs(np.mean(np.cos(m.ell*phi))) for m in s.modes if m.ell)
    record("centred_launch_forbidden_angular_overlap", symmetry, 1e-12)

    r = np.linspace(0,2*f.a,500)
    one = np.zeros(len(s.active), complex); one[1] = 1.
    record("single_mode_intensity_length_invariance", np.max(abs(
        abs(s.output(r,0,1,one))**2-abs(s.output(r,123.456,1,one))**2)), 1e-12)
    two = np.zeros(len(s.active),complex); two[:2] = 1/np.sqrt(2)
    beat_mm = 2*np.pi/abs(s.active[0].beta-s.active[1].beta)/1000
    record("two_mode_one_beat_recurrence", np.max(abs(
        abs(s.output(r,0,1,two))**2-abs(s.output(r,beat_mm,1,two))**2)), 1e-12)

    # Independent real-space radial finite-volume discretization for ell=0.
    # Transform the r-weighted operator by sqrt(r), then diagonalize tridiagonal H.
    N = 4800
    h = 4*f.a/N
    ri = (np.arange(N)+.5)*h
    left = np.arange(N)*h
    right = (np.arange(N)+1)*h
    potential = np.where(ri<f.a,0.,f.k0**2*(f.n1**2-f.n2**2))
    diag = (left+right)/(ri*h*h)+potential
    off = -right[:-1]/(h*h*np.sqrt(ri[:-1]*ri[1:]))
    values = eigh_tridiagonal(diag,off,select="i",select_range=(0,len(s.active)-1),
                             eigvals_only=True)
    u_fd = f.a*np.sqrt(values)
    u_exact = np.array([m.u for m in s.active])
    record("independent_finite_volume_U_relative_error", np.max(abs(u_fd/u_exact-1)), 3e-5)

    # LP11 cutoff: V=2 is single mode; V=2.5 supports LP11 too.
    def at_V(v):
        return solve_modes(Fibre(wavelength=2*np.pi*f.a*f.NA/v))
    record("single_mode_cutoff_test_errors", int(len(at_V(2.))!=1)+
           int(not any(m.ell==1 for m in at_V(2.5))), .5)
    report = {"all_passed": all(c["passed"] for c in checks.values()),
              "checks": checks, "LP01_LP02_beat_length_mm": beat_mm,
              "summary": s.summary()}
    return report


if __name__ == "__main__":
    report = run_validation()
    Path("results").mkdir(exist_ok=True)
    Path("results/validation.json").write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
    if not report["all_passed"]:
        raise SystemExit(1)
