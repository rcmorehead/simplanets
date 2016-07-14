#import model
from simpleabc.simple_abc import Model
from scipy import stats
import numpy as np
import simple_lib
from kports.KeplerPORTs_utils import detection_efficiency as detect


class MyModel(Model):

    """
    Model is a follows:
    """
    #@profile
    def __init__(self, stars):
        self.stars = stars
        #self.data = data
        #self.data_sum_stats = self.summary_stats(self.data)

    #functions for pickling
    #@profile
    def __getstate__(self):
        result = self.__dict__.copy()
        result['prior'] = [p.kwds for p in self.prior]
        return result

       #@profile
    def __setstate__(self, state):
        np.random.seed()
        self.__dict__ = state
        new_prior = [stats.uniform(**state['prior'][0]),
                     stats.uniform(**state['prior'][1]),
                     stats.uniform(**state['prior'][2]),
                     stats.uniform(**state['prior'][3])]
        self.__dict__['prior'] = new_prior

    #@profile
    def draw_theta(self):
        theta = []
        for p in self.prior:
            theta.append(p.rvs())
        return theta

    #@profile
    def generate_data(self, theta):

        #Draw the random model parameters.
        if (theta[0] < 0 or theta[1] < 0 or theta[2] < 0 or theta[3] < 0 or
            theta[0] > 90.0 or theta[1] > 1 or theta[2] > 20 or theta[3] > 1):


            #planet_numbers = np.ones(1)
            #total_planets = planet_numbers.sum()
            #catalog, star_header, planet_header = self.init_catalog(
                                                        #total_planets)
            return np.array([])

        else:
            select_stars = np.random.choice(self.stars,
                size=int(np.around(theta[3]*self.stars.size)), replace=False)

            planet_numbers = (self.planets_per_system(theta[2],
                          select_stars['ktc_kepler_id'].size))
            total_planets = planet_numbers.sum()
            catalog, star_header, planet_header = self.init_catalog(
                                                        total_planets)


        fund_plane_draw = self.fundamental_plane(select_stars.size)
        catalog['fund_plane'] = np.repeat(fund_plane_draw, planet_numbers)

        catalog['period'] = self.planet_period(total_planets)
        catalog['mi'] = self.mutual_inclination(theta[0], total_planets)

        catalog['fund_node'] = self.fundamental_node(total_planets)
        catalog['e'] = self.eccentricity(theta[1], total_planets)
        catalog['w'] = self.longitude_ascending_node(total_planets)
        catalog['planet_radius'] = self.planet_radius(total_planets)
        for h in star_header:
            catalog[h] = np.repeat(select_stars[h], planet_numbers)


        # print catalog.dtype.names

        catalog = catalog[(catalog['period'] >= 10.0) &
                          (catalog['period'] <= 320.0)]


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
        catalog = np.extract((catalog['b'] < 1.0), catalog)

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
            multies = simple_lib.multi_count(data, self.stars)
            h = np.histogram(multies, bins=range(0, int(multies.max()) + 1),
                density=True)
            #h = H[0][1::]/float(sum(H[0][1::]))
            #multie_ratio = h[0][2:].sum()/float(h[0][1])
            #return (simple_lib.xi(simple_lib.multies_only(data))[0],
            #        multies, multie_ratio, data.size)

            xi = simple_lib.xi(simple_lib.multies_only(data))[0]
            #xi array must have at lest two unique elements for kde
            if xi.size == 1:
                xi = np.concatenate([xi, xi+0.0005])
            if xi.size == 0:
                xi = np.array([0.0, 0.0005])
            g = stats.gaussian_kde(xi)
            return(multies, xi)
            

    #@profile
    def distance_function(self, summary_stats, summary_stats_synth):

        if summary_stats == False or summary_stats_synth == False:
            return 1e9
        #KS Distance for xi

        d1 = simple_lib.hellinger_cont(summary_stats_synth[0], 
                                          summary_stats[0])

        d2 = simple_lib.hellinger_disc(summary_stats_synth[1], 
                                          summary_stats[1])

        #Thresholds set ahead of time by inspection
        d1_threshold, d2_threshold = 100, 100

        if self.epsilon > d1_threshold and self.epsilon > d2_threshold:
            d = max((d1, d2))

        if d1_threshold > self.epsilon > d2_threshold:
            if d1 < d1_threshold:
                d = d2
            else: 
                d = max((d1, d2))

        if d1_threshold >  d2_threshold > self.epsilon:
            if d1 < d1_threshold and d2 < d2_threshold:         
                d = self.epsilon
            else: 
                d = max((d1, d2))

        return d1, d2

    #@profile
    def planets_per_system(self, Lambda, size):
        return stats.poisson.rvs(Lambda, size=size)

    #@profile
    def planet_period(self, size):
        return 10**stats.uniform.rvs(1, np.log10(320) - 1, size=size)

    #@profile
    def fundamental_node(self, size):
        return stats.uniform.rvs(0, 360, size=size)

    #@profile
    def fundamental_plane(self, size):
        return np.degrees(np.arccos(2*stats.uniform.rvs(0, 1,
                         size=size) -1 ))

    #@profile
    def mutual_inclination(self, scale, size):
        return stats.rayleigh.rvs(scale=scale, size=size)

    #@profile
    def eccentricity(self, scale, size):
        edraw = stats.rayleigh.rvs(scale=scale, size=size)
        return np.where(edraw >= 1.0, 0.99, edraw)

    #@profile
    def longitude_ascending_node(self, size):
        return stats.uniform.rvs(0, 360, size)

    #@profile
    def planet_radius(self, size):
        return (10**stats.uniform.rvs(np.log10(1.0), np.log10(19.0),
                size=size))
