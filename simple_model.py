#import model
from simpleabc.simple_abc import Model
from scipy import stats
import numpy as np
import simple_lib
from  kports.KeplerPORTs_utils import detection_efficiency as detect

class MyModel(Model):

    """
    Model is a follows:
    """
    #@profile
    def __init__(self, stars):
        self.stars = stars
        #self.data = data
        #self.data_sum_stats = self.summary_stats(self.data)

    #@profile
    def draw_theta(self):
        theta = []
        for p in self.prior:
            theta.append(p.rvs())
        return theta

    #@profile
    def generate_data(self, theta):

        #Draw the random model parameters.
        if (theta[0] < 0 or theta[1] < 0 or theta[2] < 0 or
            theta[0] > 90.0 or theta[1] > 1 or theta[2] > 20):

            #planet_numbers = np.ones(1)
            #total_planets = planet_numbers.sum()
            #catalog, star_header, planet_header = self.init_catalog(
                                                        #total_planets)
            return np.array([])

        else:
            planet_numbers = (self.planets_per_system(theta[2],
                          self.stars['ktc_kepler_id'].size))
            total_planets = planet_numbers.sum()
            catalog, star_header, planet_header = self.init_catalog(
                                                        total_planets)




        catalog['period'] = self.planet_period(total_planets)
        catalog['mi'] = self.mutual_inclination(theta[0], total_planets)

        catalog['fund_node'] = self.fundamental_node(total_planets)
        catalog['e'] = self.eccentricity(theta[1], total_planets)
        catalog['w'] = self.longitude_ascending_node(total_planets)
        catalog['planet_radius'] = self.planet_radius(total_planets)
        for h in star_header:
            catalog[h] = np.repeat(self.stars[h], planet_numbers)

        # print catalog.dtype.names

        catalog = catalog[(catalog['period'] >= 10.0) &
                          (catalog['period'] <= 320.0)]

        catalog['fund_plane'] = self.fundamental_plane(catalog)
        #Compute derived parameters.
        catalog['a'] = simple_lib.semimajor_axis(catalog['period'],
                                                 catalog['mass'])

        catalog['i'] = simple_lib.inclination(catalog['fund_plane'],
                                              catalog['mi'],
                                              catalog['fund_node'])

        catalog['b'] = simple_lib.impact_parameter(catalog['a'], catalog['e'],
                                                   catalog['i'], catalog['w'],
                                                   catalog['radius'])

        #Strip non-transiting planets/ unbound
        catalog = np.extract((catalog['b'] < 1.0) & (catalog['e'] <= 1.0),
                             catalog)

        catalog['T'] = simple_lib.transit_duration(catalog['period'],
                                                   catalog['a'], catalog['e'],
                                                   catalog['i'], catalog['w'],
                                                   catalog['b'],
                                                   catalog['radius'],
                                                   catalog['planet_radius'])

        catalog['depth'] = simple_lib.transit_depth(catalog['radius'],
                                                    catalog['planet_radius'])

        catalog['snr'] = simple_lib.snr(catalog)

        # #Strip nans from T  (planets in giant stars)
        catalog = np.extract((~np.isnan(catalog['snr'])
                              == True), catalog)
        rand_detect = stats.uniform.rvs(size=catalog.size)
        catalog = catalog[ detect(catalog['snr'], 7.1, 2) >= rand_detect ]

        return catalog

    #@profile
    def init_catalog(self, total_planets):

        star_header = ['ktc_kepler_id', 'teff', 'teff_err1', 'logg', 'feh',
                   'feh_err1', 'mass', 'mass_err1', 'radius', 'radius_err1',
                   'cdpp3', 'cdpp6', 'cdpp12', 'kepmag', 'days_obs']

        planet_header = ['b', 'i', 'a', 'planet_mass', 'planet_radius', 'T',
                         'period', 'mi', 'fund_plane', 'fund_node', 'e',
                         'w', 'depth', 'snr']

        #Initalize synthetic catalog.
        catalog = np.zeros(total_planets,
                            dtype={'names': star_header + planet_header,
                            'formats': (['i8'] + ['f8'] *
                                        (len(star_header + planet_header)
                                        - 1))})
        return catalog, star_header, planet_header

    #@profile
    def summary_stats(self, data):
        #xi(data)
        #return [0,0,0]
        #return (simple_lib.normed_duration(data), simple_lib.multi_count(data),
        #        simple_lib.xi(data)[1])
        #print data.dtype.names
        if data.size == 0:
            return False
        else:
            return (simple_lib.xi(simple_lib.multies_only(data))[0],
                    simple_lib.multi_count(data, self.stars))

    #@profile
    def distance_function(self, summary_stats, summary_stats_synth):
        N = self.stars['ktc_kepler_id'].size
        if summary_stats == False or summary_stats_synth == False:
            return 1e9
        #KS Distance for xi
        d1 = 1 - stats.ks_2samp(summary_stats[0], summary_stats_synth[0])[1]

        #Histogram distance for count

        max1 = summary_stats[1].max()
        max2 = summary_stats_synth[1].max()

        maxbin = int(max(max1, max2))
        h1 = np.histogram(summary_stats[1], bins=range(0, maxbin+1) )
        h2 = np.histogram(summary_stats_synth[1], bins=range(0, maxbin+1) )

        p = np.zeros(maxbin)
        for i in xrange(maxbin):
            p[i] = 1 - stats.binom_test(h1[0][i], N, h2[0][i]/N)

        d = np.sqrt(np.sum(p**2) + d1**2)

        return d

    #@profile
    def planets_per_system(self, Lambda, size):
        return stats.poisson.rvs(Lambda, size=size)

    #@profile
    def planet_period(self, size):
        return 10**stats.uniform.rvs(0, 3, size=size)

    #@profile
    def fundamental_node(self, size):
        return stats.uniform.rvs(0, 360, size=size)

    #@profile
    def fundamental_plane(self, catalog):
        draws = np.degrees(np.arccos(2*stats.uniform.rvs(0, 1,
                         size=np.arange(0, catalog['ktc_kepler_id'].max() + 1,
                                        1).size) -1 ))
        return draws[catalog['ktc_kepler_id']]

    #@profile
    def mutual_inclination(self, scale, size):
        return stats.rayleigh.rvs(scale, size=size)

    #@profile
    def eccentricity(self, scale, size):
        return stats.rayleigh.rvs(scale=scale, size=size)

    #@profile
    def longitude_ascending_node(self, size):
        return stats.uniform.rvs(0, 360, size)

    #@profile
    def planet_radius(self, size):
        return (10**stats.uniform.rvs(np.log10(1.0), np.log10(19.0),
                size=size))

