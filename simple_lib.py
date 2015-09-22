"""
Useful classes and functions for SIMPLE.
"""
import numpy as np
from scipy import stats
from scipy import integrate

import warnings
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
    id = np.roll(catalog['ktc_kepler_id'], 1)

    idx = np.where(catalog['ktc_kepler_id'] == id)

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


#def _kull_lei(x, P, Q):
#   return P(x) * (np.log(P(x)) - np.log(Q(x)))

#def kullback_leibler(O,E, discrete=False):
    # '''
    # If discrete, then O,E should be in frequencies/probabilities
    #
    # :param O:
    # :param E:
    # :param discrete:
    # :return:
    # '''
    #
    # if discrete == True:
    #     return sum(O * (np.log(O) - np.log(E)))
    # else:
    #     P = stats.gaussian_kde(O)
    #     Q = stats.gaussian_kde(E)
    #     return integrate.quad(_kull_lei, -np.inf, np.inf, args=(P, Q))
    #

def g_test(O, E):
    E = E[np.where(O > 0)]
    O = O[np.where(O > 0)]
    return 2 * sum(O * (np.log(O) - np.log(E)))