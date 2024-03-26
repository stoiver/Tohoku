
import okada


def forward(x,y, xoff=0, yoff=0, 
                   depth=5e3, length=1e3, width=1e3,
                   slip=0.0, opening=10.0,
                   strike=0.0, dip=0.0, rake=0.0,
                   nu=0.25, 
                   E_subfault=1, N_subfault=1):
    
    from numpy import zeros_like

    epicenters_E, epicenters_N, epicenters_D = subfaults(E_subfault, N_subfault, dip, strike, length, width)


    # initialise the value of the sum of the displacement of each subfaults
    uE_sum=zeros_like(x)
    uN_sum=zeros_like(x)
    uZ_sum=zeros_like(x)

    # the length and width of subfaults
    length_E=width/E_subfault
    length_N=length/N_subfault

    # calculate the sum of displacements of subfaults

    for i in range(N_subfault):
        for j in range(E_subfault):
            E_convert = x-epicenters_E[i,j]
            N_convert = y-epicenters_N[i,j]
            D_convert = depth-epicenters_D[i,j]

            #slip_local = slip*(1.0 + 0.5*np.random.normal())

            params =dict(x=E_convert, y=N_convert, xoff=0.0, yoff=0.0,
                         depth=D_convert, length=length_N, width=length_E,
                         slip=slip, opening=opening,
                         strike=strike, dip=dip, rake=rake,
                         nu=0.25)

            uE,uN,uZ = okada.forward(**params)

            uE_sum=uE_sum+uE;
            uN_sum=uN_sum+uN;
            uZ_sum=uZ_sum+uZ;

    return uE_sum, uN_sum, uZ_sum

def subfaults(E_subfault, N_subfault, dip, strike, length, width):
    """  
    Given the coordinate of the fault is (0,0,0), find out the coordinates of subfaults after rotation.
    E_subfault denote the number of subfaults in a row from east to west, 
    N_subfault denote the number of subfaults in a row from north to south
    length: length North
    width: width East

    epicenters_E denote the E_axis of the center of these subfaults, 
    epicenters_N denote the N_axis of the center of these subfaults.
    epicenters_H denote the depth of the center of these subfaults.
    """

    from numpy import zeros, pi, sin, cos
    
    epicenters_E=zeros((N_subfault,E_subfault));
    epicenters_N=zeros((N_subfault,E_subfault));
    epicenters_D=zeros((N_subfault,E_subfault));

    #Compute the length and width of each subfault.
    subfault_width=width/E_subfault;
    subfault_length=length/N_subfault;

    # Convert the Angle system to the radian system
    dip_angle=dip/180*pi;
    strike_angle=strike/180*pi;

    for i in range(N_subfault):
        for j in range(E_subfault):
        
            # Compute the E,N coordinates of the center of these subfaults before rotation.
            E=-width/2+subfault_width*(j+1)  -subfault_width/2;
            N=length/2-subfault_length*(i+1)+subfault_length/2;

            # Compute the E,N coordinates of the center of these subfaults after rotation
            epicenters_E[i,j]=E*cos(dip_angle)*cos(strike_angle)+N*sin(strike_angle);
            epicenters_N[i,j]=N*cos(strike_angle)-E*cos(dip_angle)*sin(strike_angle);
            epicenters_D[i,j]=-E*sin(dip_angle);
                  
    return epicenters_E, epicenters_N, epicenters_D


