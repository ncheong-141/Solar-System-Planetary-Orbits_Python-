"""
Classes for Solar system simulation 
"""

# Import external modules
import numpy as np 

# Import simulation files 
import Solar_system_functions as SSF

class Celestial_body :  # i.e. Sun and planets 
    # Anything not in __INIT__ is a static class variable. (all classes have it) 
   
    # Focal point parameters 
    m_sun = 1.989e30;
    
    # Discretization parameters for plotting
    NP_theta = 360;

    # When you define __INIT__ it is the constructor, so will execture the constructor when created
    def __init__(self, name, mass, diameter, semi_major_axis, pericentre, apocentre, orbital_period, orbital_inclination, eccentricity):
        
        # Set input parameter variables 
        self.name       = name; 
        self.m          = mass; 
        self.diameter   = diameter;
        self.T          = orbital_period; 
        self.inc        = orbital_inclination; 
        self.a          = semi_major_axis; 
        self.ra         = apocentre; 
        self.rp         = pericentre;
        self.e          = eccentricity;
        
        # Predefine variables for simulation (Dont need to in python but for C++ yes) 
        self.r = float(0);         # Current radius values
        self.x = float(0);         # "" "" x position 
        self.y = float(0);         # "" "" y position  
        self.z = float(0);         # "" "" z position 
        self.h = float(0);         # "" ""  Angular momentum 
        self.theta = float(0);     # "" "" angular position (True anomoly)

        
        # Convert units to correct data for simulation
        self.inc = self.inc*180/np.pi;
        
        # -- Pre-calculate key characterstics for the orbit of the celestial body --        
        self.b  = self.a*(1-self.e**2)**(1/2);                                  # Semi-minor axis
        self.mu = SSF.gravitational_parameter(self.m_sun, self.m);              # Gravitational parameter
        self.OF = self.a*self.e;                                                # Distance between ellipse origin and focal point.
       
        # Key points on orbit  (3D) self.OF*np.cos(self.inc + np.pi)
        self.orbit_origin   = np.array([ self.OF*np.cos(self.inc + np.pi),      0,      self.OF*np.sin(self.inc + np.pi) ]);   # Orbit origin coordinates 
        self.peri_coord     = np.array([ self.rp*np.cos(self.inc),              0,      self.rp*np.sin(self.inc) ]);           # Pericentre coordinates 
        self.apoc_coord     = np.array([ self.ra*np.cos(self.inc + np.pi),      0,      self.ra*np.sin(self.inc + np.pi) ]);
        
        """ Calculate orbit trajectories for plotting """
        self.r_plot = np.zeros(self.NP_theta);                                       # Predefine arrays
        self.x_plot = np.zeros(self.NP_theta);
        self.y_plot = np.zeros(self.NP_theta); 
        self.z_plot = np.zeros(self.NP_theta); 
        
        self.theta_plot = 0;                                                    # Initiate true anomoly for plotting. 
        self.increment = (360/self.NP_theta)*(np.pi/180);                       # Calculate increment based off discretization parameters 
        
        
        for i in range(0,self.NP_theta) :
            
            # Calculate radius of orbit at anguler point, theta_plot 
            self.r_plot[i] = SSF.oribital_radius_GP(self.a, self.e, self.theta_plot) # Calculate radius
            
            # Calculate cartesian coordinates based off radius and true anomoly 
            self.x_plot[i] = SSF.x_coordinate(self.OF, self.r_plot[i], self.theta_plot); 
            self.y_plot[i] = SSF.y_coordinate(self.r_plot[i], self.theta_plot);
            self.z_plot[i] = SSF.z_coordinate(self.x_plot[i], self.inc, self.theta_plot); 
            
            
            # Increment theta_plot for next iteration 
            self.theta_plot = self.theta_plot + self.increment; 
        
        
        

class Sub_Celestial_body : # i.e Moons 
    pass


