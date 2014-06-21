"""
Useful classes and functions for SIMPLE.
"""
import numpy as np
import warnings
r_sun_au = 0.004649
r_earth_r_sun = 0.009155
day_hrs = 24.0

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

    #duration = (p / np.pi *
    #        np.arcsin((r_star * r_sun_au) / a * 1 / np.sin(np.radians(i)) *
    #                  np.sqrt((1 - (r_planet * r_earth_r_sun) / r_star)**2
    #                          - b**2)) *
    #        1 / (1 + e*np.sin(np.radians(w))) * np.sqrt(1 - e**2)) * day_hrs#

    #return duration


# def xi(catalog):
#     """
#     Compute xi for all planet pairs in catalog.
#
#     """
#     catalog.sort(order=['ktc_kepler_id', 'period'])
#
#     out = []
#     for j in xrange(1, catalog['period'].size):
#         if catalog['period'][j] > catalog['period'][j-1]:
#             if catalog['T'][j] > 0.0 and catalog['T'][j-1] > 0.0:
#                 out.append(np.log10(((catalog['T'][j-1]/24.0) /
#                                      (catalog['T'][j]/24.0)) *
#                            (catalog['period'][j] /
#                             catalog['period'][j-1])**(1/3.)))
#
#
#     return np.array(out)

# def xi(catalog):
#     """
#     Compute xi for all planet pairs in catalog.
#
#     """
#
#     catalog.sort(order=['ktc_kepler_id', 'period'])
#
#     out = []
#     for j in xrange(1, catalog['period'].size):
#         if catalog['ktc_kepler_id'][j] == catalog['ktc_kepler_id'][j-1]:
#             out.append(np.log10(((catalog['T'][j-1]/24.0) /
#                                 (catalog['T'][j]/24.0)) *
#                                 (catalog['period'][j] /
#                                  catalog['T'][j-1])**(1/3.)))
#
#
#     return np.array(out)


def snr(catalog):
    """
    Calculate Signal to Noise ratio for a planet transit
    """

    return catalog['depth']/catalog['cdpp6'] * np.sqrt((catalog['days_obs'] /
                                                        catalog['period']) *
                                                        catalog['T']/6.0)

def xi(catalog):
    catalog.sort(order=['ktc_kepler_id', 'period'])
    p_in = np.roll(catalog['period'], 1)
    t_in = np.roll(catalog['T'], 1)
    id = np.roll(catalog['ktc_kepler_id'], 1)

    idx = np.where(catalog['ktc_kepler_id'] == id)

    P_ratio = catalog['period'][idx]/p_in[idx]
    D_ratio = t_in[idx]/catalog['T'][idx]

    #idx = np.where(P_ratio >= 1.0)
    #print P_ratio
    return np.log10(D_ratio * P_ratio**(1./3.))
