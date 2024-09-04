

import anuga
import numpy as np
import okada_kl_subfaults as okl

def create_domain(sample_id=0):

    import project
    print ('project name: ', project.name_stem)

    domain_name = f'tohoku_source_example_{int(sample_id):03d}'
    print (f'project domain name: {domain_name}') 

    domain = anuga.create_domain_from_regions(project.bounding_polygon,
                                            boundary_tags={'bottom': [0],
                                                        'ocean_east': [1],
                                                        'top': [2],
                                                        'onshore': [3]},
                                            maximum_triangle_area=project.res_whole,
                                            mesh_filename=project.meshname,
                                            interior_regions=project.interior_regions,
                                            use_cache=False,
                                            verbose=False)



    domain.set_hemisphere('northern')
    domain.set_zone(54)
    
    print ('Number of triangles = ', len(domain))
    print ('The extent is ', domain.get_extent())
    print ('Geo reference is ',domain.geo_reference)

    domain.set_quantity('elevation',
                            filename=project.name_stem + '.pts',
                            use_cache=True,
                            verbose=False,
                            alpha=0.1)

    domain.set_name(domain_name)


    #tide = -0.45
    tide = 0.0
    domain.set_quantity('stage', tide)

    Elevation = domain.quantities['elevation'].centroid_values
    Stage     = domain.quantities['stage'].centroid_values

    Stage[:] = np.maximum(Elevation, Stage)

    return domain

def apply_deformation(domain, z, parameters, x0, y0, E_subfault, N_subfault):

    # x,y and x_off, y_off are in UTM coordinates or relative to the domain
    x = domain.centroid_coordinates[:,0] # relative
    y = domain.centroid_coordinates[:,1] # relative

    xll = domain.geo_reference.xllcorner
    yll = domain.geo_reference.yllcorner

    xoff = x0 - xll
    yoff = y0 - yll

    

    uE,uN,uZ, _ = okl.kl_deformation(x, y, xoff=xoff, yoff=yoff, E_subfault=E_subfault, N_subfault=N_subfault, 
                                    sample=z, iseed=1234, **parameters)


    Elevation = domain.quantities['elevation'].centroid_values
    Stage     = domain.quantities['stage'].centroid_values

    print ('Stage extent (orig)',np.max(Stage), np.min(Stage))

    Elevation[:] = Elevation + uZ
    Stage[:]     = Stage + uZ

    print ('uZ extent          ',np.max(uZ), np.min(uZ))

    print ('Stage extent       ',np.max(Stage), np.min(Stage))
    print ('Elevation extent   ',np.max(Elevation), np.min(Elevation))

def evolve_domain(domain, gauges):

    tide = 0.0

    # Setup boundaries
    Br = anuga.Reflective_boundary(domain)
    Bf = anuga.Flather_external_stage_zero_velocity_boundary(domain,lambda t :tide)
    # Boundary conditions
    domain.set_boundary({'onshore': Br,
                            'bottom': Bf,
                            'ocean_east': Bf,
                            'top': Bf})


    from anuga.operators.collect_max_quantities_operator import Collect_max_quantities_operator
    max_min_collector = Collect_max_quantities_operator(domain)

    # Evolve system through time
    import time
    t0 = time.time()
    min = 60.
    hour = 3600.

    
    gauge_series = {}
    gauge_keys = [21418, 0, 1, 2]

    for key in gauge_keys:
        gauge_series[key] = Gauge_series(gauges[key], domain)


    # Initial run without any event
    for t in domain.evolve(yieldstep=2*min, finaltime=2.0*hour):

        for key in gauge_keys:
            gauge_series[key].append_stage()

        domain.print_timestepping_statistics()

    print ('That took %.2f seconds' %(time.time()-t0))

    return gauge_series, max_min_collector



class Gauge:

    def __init__(self, long, lat, t1=0.0, t2= 1.e10):

        self.lat  = lat
        self.long = long

        self.t1 = t1
        self.t2 = t2

        self.time_series = []
        self.stage_series = []

        import utm

        east, north, zone_no, zone_letter =  utm.from_latlon(self.lat, self.long, force_zone_number=54)

        self.east = east
        self.north = north
        self.zone_no = zone_no
        self.zone_letter = zone_letter


class Gauge_series:

    def __init__(self, gauge, domain):

        self.gauge = gauge
        self.domain = domain

        self.time_series = []
        self.stage_series = []

        self.set_triangle()

    def set_triangle(self):

        domain = self.domain
        gauge = self.gauge
        self.triangle = domain.get_triangle_containing_point([gauge.east, gauge.north])

    def append_stage(self):

        tid = self.triangle
        domain = self.domain

        stage = domain.quantities['stage'].centroid_values[tid]
        time = domain.get_time()

        self.time_series.append(time)
        self.stage_series.append(stage)



