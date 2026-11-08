"""Independent, synthetic asteroseismic scaling and sampling diagnostics."""
import math

def scaling(nu_max, delta_nu, teff, solar=(3090., 135.1, 5772.)):
    """Return (M/Msun, R/Rsun) from uncorrected solar scaling relations."""
    if min(nu_max, delta_nu, teff, *solar) <= 0:
        raise ValueError('frequencies and temperatures must be positive')
    a, b, t = nu_max / solar[0], delta_nu / solar[1], teff / solar[2]
    return a**3 * b**-4 * t**1.5, a * b**-2 * math.sqrt(t)

def log_uncertainty(nu_max, delta_nu, teff, covariance):
    """Propagate covariance of (nu_max, delta_nu, Teff) to log(M), log(R)."""
    vals = (nu_max, delta_nu, teff)
    if min(vals) <= 0 or len(covariance) != 3 or any(len(row) != 3 for row in covariance):
        raise ValueError('expected positive observables and a 3x3 covariance')
    for i in range(3):
        for j in range(3):
            if not math.isfinite(covariance[i][j]) or abs(covariance[i][j]-covariance[j][i]) > 1e-10:
                raise ValueError('covariance must be finite and symmetric')
    gradients = ((3/nu_max, -4/delta_nu, 1.5/teff),
                 (1/nu_max, -2/delta_nu, .5/teff))
    out = [[sum(gradients[a][i]*covariance[i][j]*gradients[b][j]
                for i in range(3) for j in range(3)) for b in range(2)] for a in range(2)]
    if any(out[i][i] < -1e-12 for i in range(2)):
        raise ValueError('invalid covariance: negative propagated variance')
    return out

def window_response(times_days, frequencies_microhz):
    """Normalized complex spectral-window power; identifies cadence aliases."""
    if not times_days or any(not math.isfinite(t) for t in times_days):
        raise ValueError('finite observation times required')
    n = len(times_days)
    return [sum(math.cos(2*math.pi*f*1e-6*t*86400) for t in times_days)**2 / n**2
            + sum(math.sin(2*math.pi*f*1e-6*t*86400) for t in times_days)**2 / n**2
            for f in frequencies_microhz]
