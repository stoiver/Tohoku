
import okada
import numpy as np
verbose = False


def kl_deformation(x, y, xoff=0, yoff=0, E_subfault=10, N_subfault=10, sample='random', iseed=None,
                depth=32000.0,
                length=300000,
                width=150000,
                strike=195.0,
                dip=14.0,
                rake=87.0,
                nu=0.25,
                slip = 20.0,
                opening = 0.0,
               ):
    """
    Calculate sea bed deformations due to a KL defined random slip field on a fault. The 
    default values are appropriate for the Tohoku earth quake.
    """
    

    # Calculate subfault coordinates
    epicenters_E, epicenters_N, epicenters_D = subfaults(E_subfault, N_subfault, dip, strike, length, width)

    # Create Karhunen–Loève correlation matrices
    slips, D, V, z, C_hat = kl_slipfield(epicenters_E, epicenters_N, epicenters_D, length, width, slip, sample, iseed)

    uE_sum, uN_sum, uZ_sum = sum_subfault_deformation(
        x, y, slips, xoff=xoff, yoff=yoff, depth=depth, length=length,
        width=width, strike=strike, dip=dip, rake=rake, nu=nu, opening=opening)

    if verbose:
        print(50*'=')
        print(np.max(uZ_sum))
        print(np.min(uZ_sum))

    return uE_sum, uN_sum, uZ_sum, slips


def sum_subfault_deformation(x, y, slips, xoff=0, yoff=0,
                             depth=32000.0,
                             length=300000,
                             width=150000,
                             strike=195.0,
                             dip=14.0,
                             rake=87.0,
                             nu=0.25,
                             opening=0.0,
                             ):
    """
    Sum the Okada surface displacement of every subfault for a given slip field.

    `slips` is indexed [i, j] to match subfaults(): i runs along strike and j
    runs down dip, with j = 0 the shallow (trench) edge and j = E_subfault - 1
    the deep landward edge.  Shared by kl_deformation() and
    deterministic_deformation() so the two differ only in how `slips` is built.
    """
    N_subfault, E_subfault = slips.shape

    epicenters_E, epicenters_N, epicenters_D = subfaults(E_subfault, N_subfault, dip, strike, length, width)

    openings = opening*np.ones_like(epicenters_E)

    # initialise the value of the sum of the displacement of each subfaults
    uE_sum = np.zeros_like(x)
    uN_sum = np.zeros_like(x)
    uZ_sum = np.zeros_like(x)

    # the length and width of subfaults
    length_E = width/E_subfault
    length_N = length/N_subfault

    # calculate the sum of displacements of subfaults
    for i in range(N_subfault):
        for j in range(E_subfault):
            x_convert = x-epicenters_E[i,j]
            y_convert = y-epicenters_N[i,j]
            d_convert = depth-epicenters_D[i,j]

            slipij    = slips[i,j]
            openingij = openings[i,j]

            params = dict(x=x_convert, y=y_convert, xoff=xoff, yoff=yoff,
                          depth=d_convert, length=length_N, width=length_E,
                          slip=slipij, opening=openingij,
                          strike=strike, dip=dip, rake=rake,
                          nu=nu)

            uE,uN,uZ = okada.forward(**params)

            uE_sum = uE_sum+uE
            uN_sum = uN_sum+uN
            uZ_sum = uZ_sum+uZ

    return uE_sum, uN_sum, uZ_sum


def deterministic_slip(N_subfault=10, E_subfault=10,
                       length=400000.0, width=150000.0,
                       M0=4.2e22, rigidity=40e9,
                       u0=0.5, sig_u=0.20,
                       v0=0.5, sig_v=0.30):
    """
    Smooth reproducible slip field: a separable Gaussian taper over the fault
    plane, scaled so the plane carries a prescribed seismic moment.

    `u` is the along-strike coordinate and `v` the down-dip coordinate, both
    normalised to (0, 1) over the plane with v = 0 at the trench.  `u0`/`v0`
    place the peak and `sig_u`/`sig_v` set how tightly slip is concentrated.

    The amplitude is *not* a free parameter -- it follows from
    M0 = rigidity * area * mean(slip).  Changing the taper therefore
    redistributes slip without touching the moment, which is what makes the
    far field and the coastal loading independently tunable: displaced volume
    (and so the DART peak) is set by M0 alone, while the down-dip distribution
    sets how hard the coast is loaded.  See CLAUDE.md, *Deterministic slip*.

    The defaults are the recommended Tohoku source: Mw 9.02, mean slip 17.5 m,
    peak slip 49.5 m, peak sea-bed uplift 13.1 m and ~1.2 m of coastal
    subsidence, with slip above 20 m spanning ~240 km along strike.
    """
    u = (np.arange(N_subfault) + 0.5)/N_subfault
    v = (np.arange(E_subfault) + 0.5)/E_subfault
    U, V = np.meshgrid(u, v, indexing='ij')

    shape = (np.exp(-(U - u0)**2/(2*sig_u**2))
             * np.exp(-(V - v0)**2/(2*sig_v**2)))

    return shape*(M0/(rigidity*length*width*shape.mean()))


def deterministic_deformation(x, y, xoff=0, yoff=0, E_subfault=10, N_subfault=10,
                              M0=4.2e22, rigidity=40e9,
                              u0=0.5, sig_u=0.20,
                              v0=0.5, sig_v=0.30,
                              depth=23000.0,
                              length=400000,
                              width=150000,
                              strike=195.0,
                              dip=14.0,
                              rake=87.0,
                              nu=0.25,
                              opening=0.0,
                              ):
    """
    Sea-bed deformation from the deterministic slip field of deterministic_slip().

    Drop-in alternative to kl_deformation(), returning the same
    (uE, uN, uZ, slips) tuple.  There is no `slip` or `iseed` argument: source
    strength is set by `M0`, so repeated calls give an identical source.
    """
    slips = deterministic_slip(N_subfault=N_subfault, E_subfault=E_subfault,
                               length=length, width=width,
                               M0=M0, rigidity=rigidity,
                               u0=u0, sig_u=sig_u, v0=v0, sig_v=sig_v)

    uE_sum, uN_sum, uZ_sum = sum_subfault_deformation(
        x, y, slips, xoff=xoff, yoff=yoff, depth=depth, length=length,
        width=width, strike=strike, dip=dip, rake=rake, nu=nu, opening=opening)

    if verbose:
        print(50*'=')
        print(np.max(uZ_sum))
        print(np.min(uZ_sum))

    return uE_sum, uN_sum, uZ_sum, slips

def subfaults(E_subfault, N_subfault, dip, strike, length, width):
    """  
    Given the coordinate of the fault is (0,0,0), find out the coordinates of subfaults after rotation.
    E_subfault denote the number of subfaults in a row from east to west, 
    N_subfault denote the number of subfaults in a row from north to south
    L_north: length North
    W_east: width East

    epicenters_E denote the E_axis of the center of these subfaults, 
    epicenters_N denote the N_axis of the center of these subfaults.
    epicenters_H denote the depth of the center of these subfaults.
    """

    from numpy import zeros, pi, sin, cos
    
    epicenters_E=zeros((N_subfault,E_subfault))
    epicenters_N=zeros((N_subfault,E_subfault))
    epicenters_D=zeros((N_subfault,E_subfault))

    #Compute the length and width of each subfault.
    subfault_width=width/E_subfault
    subfault_length=length/N_subfault

    # Convert the Angle system to the radian system
    dip_angle=dip/180*pi
    strike_angle=strike/180*pi

    for i in range(N_subfault):
        for j in range(E_subfault):
        
            # Compute the E,N coordinates of the center of these subfaults before rotation.
            E=-width/2+subfault_width*(j+1)  -subfault_width/2
            N=length/2-subfault_length*(i+1)+subfault_length/2

            # Compute the E,N coordinates of the center of these subfaults after rotation
            epicenters_E[i,j]=E*cos(dip_angle)*cos(strike_angle)+N*sin(strike_angle)
            epicenters_N[i,j]=N*cos(strike_angle)-E*cos(dip_angle)*sin(strike_angle)
            epicenters_D[i,j]=-E*sin(dip_angle)
                  
    return epicenters_E, epicenters_N, epicenters_D

def kl_correlation_matrices(epicenters_E, epicenters_N, epicenters_D, length, width, slip):

    from math import exp, sqrt
    from numpy import linalg as LA

    n,m = epicenters_E.shape

    vector_E = epicenters_E.flatten()
    vector_N = epicenters_N.flatten()
    vector_D = epicenters_D.flatten()

    N=len(vector_E)
    C_hat=np.zeros((N,N),dtype=float)

    #print(C_hat)

    mu=slip

    # parameters to define correlation function.
    #
    # alpha is the coefficient of variation of the slip field (sigma = alpha*mu)
    # and the expansion is Gaussian with no positivity constraint, so a large
    # alpha produces *negative* slip.  It was 0.75, which on a 400 x 100 km
    # plane at 10x10 gave slip from -26.9 m to 100.5 m with 21 of 100 subfaults
    # slipping backwards, and pulled the realised mean well below the nominal
    # `slip`.  0.4 is where negative slip disappears (CoV 0.47), and is close to
    # the scatter of published Tohoku inversions.  Changing it changes the
    # realised mean slip for a given `slip`, so recalibrate after touching it.
    #
    # r0 is the correlation length.  At 0.2*width it is shorter than the
    # along-strike subfault spacing, so the field is nearly uncorrelated in that
    # direction; a longer r0 smooths it but *increases* the spread at fixed
    # alpha, so the two want tuning together.
    alpha=0.4
    sigma=alpha*mu
    r0=0.2*width

    for i in range(N):
        for j in range(N):
            K = sqrt((vector_E[i]-vector_E[j])**2 + (vector_N[i]-vector_N[j])**2 + (vector_D[i]-vector_D[j])**2)
            C_hat[i,j] = sigma**2 * exp(-K/r0)

    #print(C_hat)

    # C_hat is a symmetric positive-semidefinite covariance matrix, so use the
    # symmetric solver.  LA.eig() calls the general LAPACK routine (dgeev),
    # which does not exploit symmetry and may return complex-conjugate
    # eigenpairs; complex D and V then propagate through sqrtD into the slip
    # field and on into okada(), crashing the run.
    #
    # Whether that happens is numpy-version dependent, not platform dependent:
    # CI with LA.eig() restored fails on numpy 2.5.2 under *both* OpenBLAS
    # (Linux) and Accelerate (macOS), and passes on numpy 2.4.6 under both.
    # LA.eigh() cannot return complex on any version, and is more accurate
    # here besides.
    D,V = LA.eigh(C_hat)

    idx = D.argsort()[::-1]
      
    D = D[idx]
    V = V[:,idx]
    # PSD in exact arithmetic, but rounding can leave eigenvalues at ~-1e-15,
    # and np.sqrt of those is NaN.
    D = np.clip(D, 0.0, None)
    D = np.diag(D)
    sqrtD = np.sqrt(D)

    return mu, n, m, D, V, sqrtD, C_hat


def kl_slipfield(epicenters_E, epicenters_N, epicenters_D, length, width, slip, sample='random', iseed=None):

    from math import exp, sqrt
    from numpy import linalg as LA
    from scipy.stats import qmc

    mu, n, m, D, V, sqrtD, C_hat = kl_correlation_matrices(epicenters_E, epicenters_N, epicenters_D, length, width, slip)

    N = len(D)

    try:
        sample.shape == (N,1)
        z = sample
    except:
        if sample == 'sobol':
            sobol_sampler = qmc.Sobol(d=N, scramble=True, seed=iseed)

        if sample is not None:
            if sample == 'random':
                if iseed is not None:
                    np.random.seed(iseed)
                z = np.random.normal(size=(N,1))
            elif sample == 'sobol':
                # The KL expansion needs standard normal deviates.  Sobol
                # points are uniform on [0,1), so map them through the normal
                # quantile function; feeding the raw uniforms in gives every
                # mode a coefficient with mean 0.5 instead of 0, which biases
                # the slip field upward and truncates its tails.
                from scipy.stats import norm
                z = norm.ppf(sobol_sampler.random()).reshape((N,1))
            else:
                msg = 'Unknown sample type %s' % sample
                raise ValueError(msg)
    

    #print(mu)
    #print(z)

    s = mu + np.dot(V,np.dot(sqrtD,z))

    s = np.reshape(s,(n,m))

    return s, D, V, z, C_hat

