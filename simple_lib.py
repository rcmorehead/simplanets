"""
Useful classes and functions for SIMPLE.
"""
import numpy as np
import warnings
import math
from scipy import integrate

r_sun_au = 0.004649
r_earth_r_sun = 0.009155
day_hrs = 24.0


#@profile
def impact_parameter(a, e, i, w, r_star):
    """
    Compute the impact parameter at for a transiting planet.

    Parameters
    ----------
    a : int, float or numpy array
        Semimajor axis of planet's orbit in AU
    e : int, float or numpy array
        Eccentricity of planet. WARNING! This function breaks down at
        high eccentricity (>> 0.9), so be careful!
    i : int, float or numpy array
        Inclination of planet in degrees. 90 degrees is edge-on.
    w : int, float or numpy array
        Longitude of ascending node defined with respect to sky-plane.
    r_star : int, float or numpy array
        Radius of star in solar radii.

    Returns
    -------
    b : float or numpy array
        The impact parameter, ie transit latitude in units of stellar radius.

    Examples
    --------
    >>> impact_parameter(1, 0, 90, 0, 1)
    1.3171077641937547e-14
    >>> a = np.linspace(.1, 1.5, 3)
    >>> e = np.linspace(0, .9, 3)
    >>> i = np.linspace(89, 91, 3)
    >>> w = np.linspace(0, 360, 3)
    >>> r_star = np.linspace(0.1, 10, 3)
    >>> impact_parameter(a, e, i, w, r_star)
    array([  3.75401300e+00,   1.66398961e-15,   1.06989371e-01])

    Notes
    -----
    Using Eqn. (7), Chap. 4, Page  56 of Exoplanets, edited by S. Seager.
    Tucson, AZ: University of Arizona Press, 2011, 526 pp.
    ISBN 978-0-8165-2945-2.
    """

    return abs(a/(r_star * r_sun_au) * np.cos(np.radians(i)) *
              (1 - e**2) / (1 + e * np.sin(np.radians(w))))


#@profile
def inclination(fund_plane, mutual_inc, node):
    """
    Compute the inclination of a planet.

    Uses the law a spherical cosines to compute the sky plane of a orbit
    given a reference plane inclination, angle from reference plane (ie mutual
    inclination) and a nodal angle.

    Parameters
    ----------
    fund_plane: int, float or numpy array
        Inclination of of the fundamental plane of the system in degrees with
        respect to the sky plane 90 degrees is edge-on.
    mutual_inc : int, float or numpy array
        Angle in degrees of the orbital plane of the planet with respect to the
        fundamental plane of the system.
    node : int, float or numpy array
        Rotation in degrees of the planet's orbit about the perpendicular of
        the reference plane. I.e. the longitude of the node with respect to the
        reference plane.

    Returns
    -------
    i : float or numpy array
        The inclination of the planet's orbit with respect to the sky plane.


    Examples
    --------
    >>> inclination(90, 3, 0)
    87.0
    >>> fun_i = np.linspace(80, 110, 3)
    >>> mi = np.linspace(0, 10, 3)
    >>> node = np.linspace(30,100,3)
    >>> inclination(fun_i, mi, node)
    array([  80.        ,   92.87347869,  111.41738591])

    Notes
    -----
    See eqn. () in
    """

    fund_plane = np.radians(fund_plane)
    mutual_inc = np.radians(mutual_inc)
    node = np.radians(node)

    return np.degrees(np.arccos(np.cos(fund_plane) * np.cos(mutual_inc) +
                      np.sin(fund_plane) * np.sin(mutual_inc) * np.cos(node)))


#@profile
def semimajor_axis(period, mass):
    """
    Compute the semimajor axis of an object.

    This is a simple implementation of the general form Kepler's Third law.

    Parameters
    ----------
    period : int, float or numpy array
        The orbital period of the orbiting body in units of days.
    mass : int, float or array-like
        The mass of the central body (or mass sum) in units of solar mass.

    Returns
    -------
    a : float or numpy array
        The semimajor axis in AU.

    Examples
    --------
    >>> semimajor_axis(365.256363,1.00)
    0.999985270598628

    >>> semimajor_axis(np.linspace(1, 1000, 5),np.linspace(0.08, 4, 5))
    array([ 0.00843254,  0.7934587 ,  1.56461631,  2.33561574,  3.10657426])

    """
    return (((2.959E-4*mass)/(4*np.pi**2))*period**2.0) ** (1.0/3.0)


#@profile
def transit_depth(r_star, r_planet):
    """
    One-line description

    Full description

    Parameters
    ----------

    Returns
    -------

    Examples
    --------


    """
    return ((r_planet * r_earth_r_sun)/r_star)**2 * 1e6


#@profile
def transit_duration(p, a, e, i, w, b, r_star, r_planet):
    """
    Compute the full (Q1-Q4) transit duration.

    Full description

    Parameters
    ----------
    p : int, float or numpy array
        Period of planet orbit in days
    a : int, float or numpy array
        Semimajor axis of planet's orbit in AU
    e : int, float or numpy array
        Eccentricity of planet. WARNING! This function breaks down at
        high eccentricity (>> 0.9), so be careful!
    i : int, float or numpy array
        Inclination of planet in degrees. 90 degrees is edge-on.
    w : int, float or numpy array
        Longitude of ascending node defined with respect to sky-plane.
    b : int, float or numpy array
        Impact parameter of planet.
    r_star : int, float or numpy array
        Radius of star in solar radii.
    r_planet : int, float or numpy array
        Radius of planet in Earth radii

    Returns
    -------
    T : float or numpy array
        The Q1-Q4 (full) transit duration of the planet in hours.

    Examples
    --------

    Notes
    -----
    Using Eqns. (15) and (16), Chap. 4, Page  58 of Exoplanets, edited by S.
    Seager. Tucson, AZ: University of Arizona Press, 2011, 526 pp.
    ISBN 978-0-8165-2945-2.
    """

    #TODO Make this robust against b > 1
    #warnings.simplefilter("always")
    #print "pars", p, a, e, i, w, b, r_star, r_planet
    #print ""
    #print (1 - (r_planet * r_earth_r_sun) / r_star)**2 - b**2
    #print (1 - e**2)
    #print ""
    duration = np.where(e < 1.0, (p / np.pi *
            np.arcsin((r_star * r_sun_au) / a * 1 / np.sin(np.radians(i)) *
                      np.sqrt((1 - (r_planet * r_earth_r_sun) / r_star)**2
                              - b**2)) *
            1 / (1 + e*np.sin(np.radians(w))) * np.sqrt(1 - e**2)) * day_hrs, 0)

    return duration


#@profile
def snr(catalog):
    """
    Calculate Signal to Noise ratio for a planet transit

    Full description

    Parameters
    ----------

    Returns
    -------

    Examples
    --------


    """


    return catalog['depth']/catalog['cdpp6'] * np.sqrt((catalog['days_obs'] /
                                                        catalog['period']) *
                                                        catalog['T']/6.0)


#@profile
def xi(catalog):
    """
    One-line description

    Full description

    Parameters
    ----------

    Returns
    -------

    Examples
    --------


    """
    catalog.sort(order=['ktc_kepler_id', 'period'])
    p_in = np.roll(catalog['period'], 1)
    t_in = np.roll(catalog['T'], 1)
    kic_id = np.roll(catalog['ktc_kepler_id'], 1)

    idx = np.where(catalog['ktc_kepler_id'] == kic_id)

    P_ratio = catalog['period'][idx]/p_in[idx]
    D_ratio = t_in[idx]/catalog['T'][idx]

    #idx = np.where(P_ratio >= 1.0)
    #print P_ratio
    logxi = np.log10(D_ratio * P_ratio**(1./3.))
    if logxi.size < 1:
        xi_fraction = 0.0
    else:
        xi_fraction = logxi[logxi >= 0.0].size/float(logxi.size)
    return logxi, xi_fraction


#@profile
def multi_count(catalog, stars):
    """
    One-line description

    Full description

    Parameters
    ----------

    Returns
    -------

    Examples
    --------


    """
    count = np.zeros(stars['ktc_kepler_id'].size)

    bincount = np.bincount(catalog['ktc_kepler_id'])
    bincount = bincount[bincount > 0]
    count[:bincount.size] = bincount

    return count


#@profile
def multies_only(catalog):

    unq, unq_idx, unq_cnt = np.unique(catalog['ktc_kepler_id'],
                                      return_inverse=True,
                                      return_counts=True)
    cnt_mask = unq_cnt > 1
    cnt_idx, = np.nonzero(cnt_mask)
    idx_mask = np.in1d(unq_idx, cnt_idx)
    return catalog[idx_mask]

def duration_anomaly(catalog):
    """
    Returns T/T_nu where T is the transit duration and T_nu is the
    duration for a e = 0, b = 0 transit.

    Full description

    Parameters
    ----------

    Returns
    -------

    Examples
    --------

    """
    catalog['T_nu'] = (catalog['T'] /
                        ((catalog['radius'] * r_sun_au * catalog['period'])
                         /(np.pi * catalog['a']) * day_hrs))
    return catalog

#@profile
def normed_duration(catalog):
    """
    One-line description

    Full description

    Parameters
    ----------

    Returns
    -------

    Examples
    --------


    """

    return (catalog['T']/day_hrs)/(catalog['period'])**(1/3.0)

def _anderson_ksamp_midrank(samples, Z, Zstar, k, n, N):
    """
    Compute A2akN equation 7 of Scholz and Stephens.
    Parameters
    ----------
    samples : sequence of 1-D array_like
        Array of sample arrays.
    Z : array_like
        Sorted array of all observations.
    Zstar : array_like
        Sorted array of unique observations.
    k : int
        Number of samples.
    n : array_like
        Number of observations in each sample.
    N : int
        Total number of observations.
    Returns
    -------
    A2aKN : float
        The A2aKN statistics of Scholz and Stephens 1987.
    """

    A2akN = 0.
    Z_ssorted_left = Z.searchsorted(Zstar, 'left')
    if N == Zstar.size:
        lj = 1.
    else:
        lj = Z.searchsorted(Zstar, 'right') - Z_ssorted_left
    Bj = Z_ssorted_left + lj / 2.
    for i in np.arange(0, k):
        s = np.sort(samples[i])
        s_ssorted_right = s.searchsorted(Zstar, side='right')
        Mij = s_ssorted_right.astype(np.float)
        fij = s_ssorted_right - s.searchsorted(Zstar, 'left')
        Mij -= fij / 2.
        inner = lj / float(N) * (N * Mij - Bj * n[i])**2 / \
            (Bj * (N - Bj) - N * lj / 4.)
        A2akN += inner.sum() / n[i]
    A2akN *= (N - 1.) / N
    return A2akN


def _anderson_ksamp_right(samples, Z, Zstar, k, n, N):
    """
    Compute A2akN equation 6 of Scholz & Stephens.
    Parameters
    ----------
    samples : sequence of 1-D array_like
        Array of sample arrays.
    Z : array_like
        Sorted array of all observations.
    Zstar : array_like
        Sorted array of unique observations.
    k : int
        Number of samples.
    n : array_like
        Number of observations in each sample.
    N : int
        Total number of observations.
    Returns
    -------
    A2KN : float
        The A2KN statistics of Scholz and Stephens 1987.
    """

    A2kN = 0.
    lj = Z.searchsorted(Zstar[:-1], 'right') - Z.searchsorted(Zstar[:-1],
                                                              'left')
    Bj = lj.cumsum()
    for i in np.arange(0, k):
        s = np.sort(samples[i])
        Mij = s.searchsorted(Zstar[:-1], side='right')
        inner = lj / float(N) * (N * Mij - Bj * n[i])**2 / (Bj * (N - Bj))
        A2kN += inner.sum() / n[i]
    return A2kN


def anderson_ksamp(samples, midrank=True):
    """The Anderson-Darling test for k-samples.
    The k-sample Anderson-Darling test is a modification of the
    one-sample Anderson-Darling test. It tests the null hypothesis
    that k-samples are drawn from the same population without having
    to specify the distribution function of that population. The
    critical values depend on the number of samples.
    Parameters
    ----------
    samples : sequence of 1-D array_like
        Array of sample data in arrays.
    midrank : bool, optional
        Type of Anderson-Darling test which is computed. Default
        (True) is the midrank test applicable to continuous and
        discrete populations. If False, the right side empirical
        distribution is used.
    Returns
    -------
    A2 : float
        Normalized k-sample Anderson-Darling test statistic.
    critical : array
        The critical values for significance levels 25%, 10%, 5%, 2.5%, 1%.
    logp : float
        The log (ln) of an approximate significance level at which the null hypothesis for the
        provided samples can be rejected.
    Raises
    ------
    ValueError
        If less than 2 samples are provided, a sample is empty, or no
        distinct observations are in the samples.
    See Also
    --------
    ks_2samp : 2 sample Kolmogorov-Smirnov test
    anderson : 1 sample Anderson-Darling test
    Notes
    -----
    [1]_ Defines three versions of the k-sample Anderson-Darling test:
    one for continuous distributions and two for discrete
    distributions, in which ties between samples may occur. The
    default of this routine is to compute the version based on the
    midrank empirical distribution function. This test is applicable
    to continuous and discrete data. If midrank is set to False, the
    right side empirical distribution is used for a test for discrete
    data. According to [1]_, the two discrete test statistics differ
    only slightly if a few collisions due to round-off errors occur in
    the test not adjusted for ties between samples.
    .. versionadded:: 0.14.0
    References
    ----------
    .. [1] Scholz, F. W and Stephens, M. A. (1987), K-Sample
           Anderson-Darling Tests, Journal of the American Statistical
           Association, Vol. 82, pp. 918-924.
    """
    k = len(samples)
    if (k < 2):
        raise ValueError("anderson_ksamp needs at least two samples")

    samples = list(map(np.asarray, samples))
    Z = np.sort(np.hstack(samples))
    N = Z.size
    Zstar = np.unique(Z)
    if Zstar.size < 2:
        raise ValueError("anderson_ksamp needs more than one distinct "
                         "observation")

    n = np.array([sample.size for sample in samples])
    if any(n == 0):
        raise ValueError("anderson_ksamp encountered sample without "
                         "observations")

    if midrank:
        A2kN = _anderson_ksamp_midrank(samples, Z, Zstar, k, n, N)
    else:
        A2kN = _anderson_ksamp_right(samples, Z, Zstar, k, n, N)

    h = (1. / np.arange(1, N)).sum()
    H = (1. / n).sum()
    g = 0
    for l in np.arange(1, N-1):
        inner = np.array([1. / ((N - l) * m) for m in np.arange(l+1, N)])
        g += inner.sum()

    a = (4*g - 6) * (k - 1) + (10 - 6*g)*H
    b = (2*g - 4)*k**2 + 8*h*k + (2*g - 14*h - 4)*H - 8*h + 4*g - 6
    c = (6*h + 2*g - 2)*k**2 + (4*h - 4*g + 6)*k + (2*h - 6)*H + 4*h
    d = (2*h + 6)*k**2 - 4*h*k
    sigmasq = (a*N**3 + b*N**2 + c*N + d) / ((N - 1.) * (N - 2.) * (N - 3.))
    m = k - 1
    A2 = (A2kN - m) / math.sqrt(sigmasq)
    return A2


def hellinger_funct(x,P,Q):
    """
    P,Q should be numpy stats gkde objects
    """
    return np.sqrt(P(x) * Q(x))

def hellinger_cont(P,Q):
    """
    P,Q should be numpy stats gkde objects
    F should be the hellinger_funct method
    """
    return 1 - integrate.quad(hellinger_funct, -np.inf, np.inf, args=(P,Q))[0]

def hellinger_disc(P,Q):
    """
    P,Q should be numpy histogram objects that have density=True
    """
    if P[0].size == Q[0].size:
        pass
    else:
        if P[0].size > Q[0].size:
            Q[0].resize(P[0].size)
        else:
            P[0].resize(Q[0].size)
        
    return  1 - np.sum(np.sqrt(P[0]*Q[0]))

def KL_funct(x, P, Q):
    """
    P,Q should be numpy stats gkde objects
    """
    return (P(x)) * ( np.log(P(x)) - np.log(Q(x)) )

def KL_cont(P,Q, limits):
    """
    P,Q should be numpy stats gkde objects
    """
    return integrate.quad(KL_funct, limits[0], limits[1], args=(P,Q))[0]

def KL_disc(P,Q):
    """
    P,Q should be numpy histogram objects that have density=False
    """
    
    if P[0].size == Q[0].size:
        pass
    else:
        if P[0].size > Q[0].size:
            Q[0].resize(P[0].size)
        else:
            P[0].resize(Q[0].size)
            
    #KL divergence is defined only if Q(i)=0 implies P(i)=0, for all i 
    #(absolute continuity)
    P = np.where(Q[0] > 0, P[0],0)
    #(re)normalize
    P = P/P.sum()
    Q = Q[0]/Q[0].sum()
    #P(i) = 0 -> P(i)log(P(i)/Q(i)) = 0
    ind = P > 0

    return np.sum((P[ind]) * ((np.log(P[ind]) - np.log(Q[ind]))))

def KL_cont_sym(P, Q, P_points, Q_points):
    Plim = (P_points.min(), P_points.max())
    Qlim = (Q_points.min(), Q_points.max())
    
    KL = KL_cont(P, Q, Plim) + KL_cont(Q, P, Qlim)
    

    if np.isfinite(KL) == False:
        KL = 1e6
    
    return KL

def KL_disc_sym(P,Q):
    KL = KL_disc(P, Q) + KL_disc(Q, P)
    
    if np.isfinite(KL) == False:
        KL = 1e6
    
    return KL


