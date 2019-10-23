"""
Solar system simulation functions 
"""

# External modules 
import numpy as np

# Constants
G = 6.6695e-11; # Gravitational constant


""" Gravitational attraction force: Newtons law of gravitational attraction (from III law) """
def gravitational_attraction(m1,m2,r) :
    return G*((m1*m2)/(r**2));

""" Gravitational parameter (From derrivation of gravitational attraction """
def gravitational_parameter(m1,m2):
    return G*(m1 + m2);

""" Polar Orbit equation """
# Oribital radius calculation using geometric parameters (returns radius)
# a     = semi-major axis (aphelion?)
# e     = orbital eccentricity 
# theta =  
def oribital_radius_GP(a,e,theta) :     # (m)
    return (a*(1 - e**2)) / (1 + e*np.cos(theta));


# Orbital radius calculation using dynamic parameters 
def orbital_radius_DP(h,mu,e,theta) :   # (m)
    return ((h**2)/mu) / (1 + e*np.cos(theta));


""" Calculating spatial coordinates of orbit """
def x_coordinate(OF,r,theta) : 
    return OF + (r*np.cos(theta));

def y_coordinate(r, theta) : 
    return r*np.sin(theta); 

def z_coordinate(x_pos, inc, theta) :
    if (theta >= 0 and theta <= np.pi/2):
        return x_pos*np.tan(inc); 
    elif (theta >= (3/2)*np.pi and theta <= 2*np.pi): 
        return x_pos*np.tan(inc); 
    else: 
        return x_pos*np.tan(inc + np.pi); 



""" Orbital energy derived from..
    Note, orbital energy is conserved anywhere on the orbit. 
"""
def orbital_energy(v,mu,r) :    # Total energy per unit mass (j/kg)
    #     Kinetic energy - Potential energy
    return ((v**2)/2) - (mu/r);


def orbital_energy_GP(mu,a) :
    return -mu/2*a; 


""" Flight path angle (angle between velocity vector and theta^ )""" 
def flight_path_angle(e,theta):
    return np.arctan((e*np.sin(theta))/(1+(e*np.cos(theta))));


""" Orbital velocity at arbitrary point, r """
def orbital_velocity(mu, r, a): 
    return ( (2*mu/r) - (mu/a) )**(1/2)


""" Angular momentum at arbitrary  point, r CHANGE THIS"""
def orbital_angular_momentum(r,mu,a,e,theta):
    v       = orbital_velocity(mu,r,a); 
    gamma   = flight_path_angle(e,theta); 
    return r*v*np.cos(gamma);


""" Calculating eccentricity """ 
def eccentricity(ra,rp) :
    return (ra-rp)/(ra+rp);


""" Velocity at pericentre and apocentre """
def velocity_pericentre(mu,a,e) :
    return ( (mu/a) * ((1+e)/(1-e)) )**(1/2);


def velocity_apocentre(mu,a,e) :
    return ( (mu/a) * ((1-e)/(1+e)) )**(1/2);

