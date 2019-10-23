""" =================================================================
Project: Simulating the solar system planetary orbits (3D)
================================================================= """

"""
-========- Notes -========- 

- Angular momentum is conserved. Gravitational force eild is scentral and since the angular momentum vector is perpendicular
  to the spatial and velocity vector. 
  
- Total energy is conserved in a closed path/conservative field

- Eccentricity vector is a vector constant of integration. It can be used as reference direction on the orbital plane (it is perpendicular to angular momentum vector)

- True anomoly is the angle between the spatial vecotr (r_) and eccentricity vector
    - Recall e_.r_ = ercos(theta) (vector stuff)
    - e < 1 = elipse 
    - e = 1 = parabola
    - e > 1 = hyperbola

- The eccentricity is the ratio between the distance between the foci and the major axis

Key Definitions: 
    Pericentre: 
        The point closest to the orbital focus point (also known as perihelion)
        
    Apocentre:  
        The point farthet to the orbital focus point (also known as the aphelion)
        
    Semi-major axis:
        
        

List of variables: 
    a  = Semi-major axis
    b  = semi-minor axis
    mu = gravitational parameter
    r  = radius
    rp = distance betweeen focal point and pericentre
    ra = distance between focal point and apocentre
    theta = true anomoly (angle between radius and eccentricity vector)
    
    
Equations to possibly implement 
- Distance from focal point to peri and apocentre
    - rp = a(1-e), ra = a(1+e)
- Semi major axis 
    - a = rp + ra / 2
- Semi minor axis 
    - b = a(1-e**2)**(1/2);
"""

# Import external modules
import numpy as np
import matplotlib.pylab as plt
import random as rand

# Import simulation files
import Solar_system_celestial_body_classes as CB
import Solar_system_functions as SSF
# Define simulation parameters 
max_time = 1000;    # Units: (days)
timestep = 1;      # (days)


# Define planet parameters and set up celestial body arrays for processing during simulation 
# Class inputs go as:
#   Celestial_body(         name                mass (kg),      diameter (km), semi_major_axis (km), pericentre (km), apocentre (km), orbital_period (days), orbital_inclination (deg), eccentricity (n/a)):

# Mercury = CB.Celestial_body(                  0.330e24,       4879,           57.9e6,                 46e6,           69.8e6,             88,                     7,                      0.205);    
# Venus   = CB.Celestial_body(                  4.8e24,         12104,          108.2e6,                107.5e6,        108.9e6,            224.7,                  3.4,                    0.007);
# Earth   = CB.Celestial_body(                  5.97e24,        12756,          149.6e6,                147.1e6,        152.1e6,            365.2,                  0,                      0.017);
# Mars    = CB.Celestial_body(                  0.642e24,       6792,           227.9e6,                206.6e6,        249.2e6,            687,                    1.9,                    0.094);
# Jupiter = CB.Celestial_body(                  1898e24,        142984,         778.6e6,                740.5e6,        816.6e6,            4331,                   1.3,                    0.049);
# Saturn  = CB.Celestial_body(                  568e24,         120536,         1433.5e6,               1352.6e6,       1514.5e6,           10747,                  2.5,                    0.057);
# Uranus  = CB.Celestial_body(                  86.8e24,        51118,          2872.5e6,               2741.3e6,       3003.6e6,           30589,                  0.8,                    0.046);
# Neptune = CB.Celestial_body(                  102e24,         49528,          4495.1e6,               4444.5e6,       4545.7e6,           59800,                  1.8,                    0.011);
# Pluto   = CB.Celestial_body(                  0.0146e24,      2370,           5906.4e6,               4436.8e6,       7375.9e6,           90560,                  17.2,                   0.244);

Planets = [  CB.Celestial_body("Mercury",   0.330e24,       4879,           57.9e6,                 46e6,           69.8e6,             88,                     7,                      0.205)   
            ,CB.Celestial_body("Venus",     4.8e24,         12104,          108.2e6,                107.5e6,        108.9e6,            224.7,                  3.4,                    0.007)
            ,CB.Celestial_body("Earth",     5.97e24,        12756,          149.6e6,                147.1e6,        152.1e6,            365.2,                  0,                      0.017)
            ,CB.Celestial_body("Mars",      0.642e24,       6792,           227.9e6,                206.6e6,        249.2e6,            687,                    1.9,                    0.094)
            ,CB.Celestial_body("Jupiter",   1898e24,        142984,         778.6e6,                740.5e6,        816.6e6,            4331,                   1.3,                    0.049)
            ,CB.Celestial_body("Saturn",    568e24,         120536,         1433.5e6,               1352.6e6,       1514.5e6,           10747,                  2.5,                    0.057)
            ,CB.Celestial_body("Uranus",    86.8e24,        51118,          2872.5e6,               2741.3e6,       3003.6e6,           30589,                  0.8,                    0.046)
            ,CB.Celestial_body("Neptune",   102e24,         49528,          4495.1e6,               4444.5e6,       4545.7e6,           59800,                  1.8,                    0.011)
            ,CB.Celestial_body("Pluto",     0.0146e24,      2370,           5906.4e6,               4436.8e6,       7375.9e6,           90560,                  17.2,                   0.244) ]

    
""" --------------------------- Start of simulation --------------------------- """

# Calculate maximum iterations from time step and maximum time
max_iteration = int(max_time/timestep);

# Randomly generate true anomoly positions for start and calculate respective values 
for i in range (0, np.size(Planets)):
    Planets[i].theta = rand.uniform(0,2*np.pi);
    
    Planets[i].r = SSF.oribital_radius_GP(Planets[i].a,Planets[i].e,Planets[i].theta);
    Planets[i].x = SSF.x_coordinate(Planets[i].OF,Planets[i].r,Planets[i].theta);
    Planets[i].y = SSF.y_coordinate(Planets[i].r, Planets[i].theta); 
    Planets[i].z = SSF.z_coordinate(Planets[i].x, Planets[i].inc, Planets[i].theta);

""" Plotting """   
# Plot trajectories 
ax =  plt.figure().add_subplot(111, projection='3d', aspect='equal');

for i in range(0,np.size(Planets)) : 
    ax.plot_wireframe(Planets[i].x_plot,Planets[i].y_plot,Planets[i].z_plot.reshape(1,-1));
    ax.scatter3D(Planets[i].x,Planets[i].y,Planets[i].z);

    plt.draw() 


"""
for it in range (0,max_iteration) :
    
    start = time.time()
    
    for i in range (0, np.size(Planets)) :
        
        # Calculate spatial position of planets
        Planets[i].r = SSF.oribital_radius_GP(Planets[i].a,Planets[i].e,Planets[i].theta);
        Planets[i].x = SSF.x_coordinate(Planets[i].OF,Planets[i].r,Planets[i].theta);
        Planets[i].y = SSF.y_coordinate(Planets[i].r, Planets[i].theta); 
        Planets[i].z = SSF.z_coordinate(Planets[i].x, Planets[i].inc, Planets[i].theta);
    
       
        
        # Update pplots of new position of planets
        ax.scatter3D(Planets[i].x,Planets[i].y,Planets[i].z)
        plt.pause(0.0001)
        
    end  = time.time()
    print(f'Time taken: {end-start:.2f} seconds')

"""
