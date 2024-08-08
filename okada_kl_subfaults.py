
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

    openings = opening*np.ones_like(epicenters_E)

    # initialise the value of the sum of the displacement of each subfaults
    uE_sum=np.zeros_like(x)
    uN_sum=np.zeros_like(x)
    uZ_sum=np.zeros_like(x)

    # the length and width of subfaults
    length_E= width/E_subfault
    length_N= length/N_subfault

    # calculate the sum of displacements of subfaults

    for i in range(N_subfault):
        for j in range(E_subfault):
            x_convert = x-epicenters_E[i,j]
            y_convert = y-epicenters_N[i,j]
            d_convert = depth-epicenters_D[i,j]

            #slip_local = slip*(1.0 + 0.5*np.random.normal())
            slipij    = slips[i,j]
            openingij = openings[i,j]

            params =dict(x=x_convert, y=y_convert, xoff=xoff, yoff=yoff,
                            depth=d_convert, length=length_N, width=length_E,
                            slip=slipij, opening=openingij,
                            strike=strike, dip=dip, rake=rake,
                            nu=0.25)

            uE,uN,uZ = okada.forward(**params)

            uE_sum=uE_sum+uE
            uN_sum=uN_sum+uN
            uZ_sum=uZ_sum+uZ

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

    # parameters to define correlation function
    alpha=0.75
    sigma=alpha*mu
    r0=0.2*width

    for i in range(N):
        for j in range(N):
            K = sqrt((vector_E[i]-vector_E[j])**2 + (vector_N[i]-vector_N[j])**2 + (vector_D[i]-vector_D[j])**2)
            C_hat[i,j] = sigma**2 * exp(-K/r0)

    #print(C_hat)
    
    D,V = LA.eig(C_hat)

    idx = D.argsort()[::-1]
      
    D = D[idx]
    V = V[:,idx]
    D = np.diag(D)
    sqrtD = np.sqrt(D)

    return mu, n, m, D, V, sqrtD, C_hat


def kl_slipfield(epicenters_E, epicenters_N, epicenters_D, length, width, slip, sample='random', iseed=None):

    from math import exp, sqrt
    from numpy import linalg as LA
    from scipy.stats import qmc

    if sample == 'sobol':
        sobol_sampler = qmc.Sobol(size=(N,1))

    mu, n, m, D, V, sqrtD, C_hat = kl_correlation_matrices(epicenters_E, epicenters_N, epicenters_D, length, width, slip)

    N = len(D)

    if sample is not None:
        if sample == 'random':
            if iseed is not None:
                np.random.seed(iseed)
            z = np.random.normal(size=(N,1))
        elif sample == 'sobol':
            z = sobol_sampler.random()
        elif sample.shape == (N,1):
            z = sample
 

    #print(mu)
    #print(z)

    s = mu + np.dot(V,np.dot(sqrtD,z))

    s = np.reshape(s,(n,m))

    return s, D, V, z, C_hat

