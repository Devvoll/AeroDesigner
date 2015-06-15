# -*- coding: utf-8 -*-
"""
Created on Tue Jun 02 15:32:17 2015

@author: Devin Vollmer

    This code visualizes the interaction between wing shape, velocity, 
and subsequent power requirments for small fixed wing RC airplanes. 
Adjust Cl, Propsize, Wing Size, and then see the impact of these variables have
on the power requirments. This is a useful tool for when considering
optimization within a prescribed design space, or to generate a first order 
estimate of wing size and power requirments.

    This code accounts for viscous flows when calculating wing efficiency,
and the presence of a fuselage when, determining total drag.  

    This code will list important design paramenters in the console as 
both MKS and IPS units, but the GUI sliders are in MKS.


Resources:
http://www.fzt.haw-hamburg.de/pers/Scholz/OPerA/OPerA_PUB_DLRK_12-09-10.pdf
http://www.aerostudents.com/files/flightDynamics/longitudinalStabilityAndControl.pdf
http://www.mh-aerotools.de/airfoils/flywing1.htm
http://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node86.html
http://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node97.html

"""
import numpy as np
import matplotlib.pyplot as plt
from numpy import pi
from matplotlib.widgets import Slider
import os

def cls():
    os.system(['clear','cls'][os.name == 'nt'])

# Clear the screen
cls()
plt.close('all')

#conversion factors
m2in = 39.3701
kg2slug = 0.0685217659
N2lbf = 0.224808943
watt2hp = 0.00134102209

fig, ax = plt.subplots(1, sharex=True)
plt.subplots_adjust(left=0.07, right=.55, bottom=0.1)
mngr = plt.get_current_fig_manager()
mngr.window.setGeometry(225,125,900, 545)

# Initial Parameters and Mission Requirments
mass = .6 #Target Mass
ro = 1.125 #Air density
g = 9.81 #Gravity
VelC = 20.0 #Desired Cruise Velocity
v = np.arange(5, VelC*1.75, .2) #Range of velocity to plot
df = .07 #Diameter of Fuselage

weight = mass*g*v/v #Weight
eta_p = .8 #Assumed Propulsive efficieny.  .8 is typical peak efficiency
Root = .2 #Length of Root Chord
Taper = .5  #Amount of taper
Tip = Root*Taper #Length of Tip Chord
Span = 1.0 #Complete Span
S = (0.5 * Span) * (Root+Tip) #Wing Area
AR = (pow(Span,2))/S #Wing Aspect Ratio
Cl = .4 #Lift Coefficient
Cd0 = .02 #Fuselage skin and body drag

u = .99 #Theoretical Oswald Efficiency Factor
Dp = .1 #Propellor 'Disk' Diameter
K = .38 #Emperical value to estimate e

# Initialization Equations
s = 1.0 - 2.0 * ( df/Span )    
Q = 1 / ( u*s )
P = K*Cd0
e = 1 / (Q + P*pi*AR)
L = 0.5*ro*pow(v,2)*S*Cl
Cdi = pow(L,2) * ( 1/( .25*pow(ro,2)*pow(v,4)*pow(S,2)*pi*e*AR)) 
D =  0.5*ro*pow(v,2)*S*Cdi   #Induced Drag/Drag Due to lift
Power = 0.5*ro*pow(v,3)*S*(Cdi+Cd0) + (( pow(mass*9.81,2) )/(0.5*ro*v*S))*(1/(pi*e*AR))
Adisk = pi*pow((0.5*Dp),2)
T = (Power/v)
TW = T / (mass*g)
PropPow = (.5*T*v) * ( pow( (T/(Adisk*pow(v,2)*ro*0.5)) + 1, .5) + 1 )
 
# Create Line2D Objects   
l, = ax.plot(v,Power, lw=2, color='blue', linestyle='--', label='Required Power Curve (Watts)')
m, = ax.plot(v,PropPow, lw=2, color='blue', linestyle='-', label='Prop Power Needed (Watts)')
k, = ax.plot(v,L, lw=2, color='green', label='Lift Force Generated (N)')
j, = ax.plot(v,weight, lw=2, color='k', linestyle='--', label='Vehicle Weight (N)')

ax.axvline(VelC, lw=2, color='red', linestyle='--', label='Mission Cruise Velocity (m/s)')

plt.legend( loc=2, borderaxespad=0.)

ax.set_ylabel("Power (Watts), Force (N)")
ax.set_xlabel("Velocity (m/s)")
ax.axis([0, VelC*1.75, 0, mass*200]) #Plotting Area Scaled to design requirments
ax.grid('on')

axcolor = 'white'

#Control location of sliders on GUI
axDp = plt.axes([0.75, 0.6, 0.2, 0.03], axisbg=axcolor) 
axSpan = plt.axes([0.75, 0.65, 0.2, 0.03], axisbg=axcolor)  
axM  = plt.axes([0.75, 0.7, 0.2, 0.03], axisbg=axcolor)
axCl  = plt.axes([0.75, 0.75, 0.2, 0.03], axisbg=axcolor)
axTaper = plt.axes([0.75, 0.8, 0.2, 0.03], axisbg=axcolor)
axRoot = plt.axes([0.75, 0.85, 0.2, 0.03], axisbg=axcolor) 

#These control the limits on the sliders
sDp = Slider(axDp, 'Propellor Diameter (m)', .05, .4, valinit=Dp)
sSpan = Slider(axSpan, 'Span (m) ', 0.1, 2.0, valinit=Span)  
sM = Slider(axM, 'Mass (kg) ', 0.1, mass*4, valinit=mass)
sCl = Slider(axCl, 'Lift Coefficient ', 0.0, 2.0, valinit=Cl)
sTaper = Slider(axTaper, 'Wing Taper Ratio ', .01, 1.0, valinit=Taper)
sRoot = Slider(axRoot, 'Root Chord Length (m)', .01, 0.5, valinit=Root)

def updatePower(val):
    
    Root = sRoot.val
    mass = sM.val
    Tip = Root*sTaper.val    
    S = (0.5 * sSpan.val) * (Root+Tip)
    AR = pow((sSpan.val),2)/S
    s = 1.0 - 2.0 * ( df/sSpan.val )    
    Q = 1 / ( u*s )
    P = K*Cd0
    e = 1 / (Q + P*pi*AR)
    Cl = sCl.val
    L = 0.5*ro*pow(v,2)*S*Cl
    Cdi = pow(L,2) * ( 1/( .25*pow(ro,2)*pow(v,4)*pow(S,2)*pi*e*AR)) 
    Power = (1/eta_p)*(0.5*ro*pow(v,3)*S*(Cdi+Cd0) + (( pow(mass*g,2) )/(0.5*ro*v*S))*(1/(pi*e*AR))) 
    l.set_ydata(  Power )
    fig.canvas.draw_idle()

sRoot.on_changed(updatePower)
sM.on_changed(updatePower)
sCl.on_changed(updatePower)  
sTaper.on_changed(updatePower)
sSpan.on_changed(updatePower)

def updateThrust(val):
    
    Root = sRoot.val
    mass = sM.val
    Tip = Root*sTaper.val    
    S = (0.5 * sSpan.val) * (Root+Tip)
    AR = pow((sSpan.val),2)/S
    s = 1.0 - 2.0 * ( df/sSpan.val )    
    Q = 1 / ( u*s )
    P = K*Cd0
    e = 1 / (Q + P*pi*AR)
    Cl = sCl.val
    L = 0.5*ro*pow(v,2)*S*Cl
    Cdi = pow(L,2) * ( 1/( .25*pow(ro,2)*pow(v,4)*pow(S,2)*pi*e*AR)) 
    Power = (1/eta_p)*(0.5*ro*pow(v,3)*S*(Cdi+Cd0) + (( pow(mass*g,2) )/(0.5*ro*v*S))*(1/(pi*e*AR))) 
    T = (Power/v)
    Adisk = pi*pow((0.5*sDp.val),2)
    PropPow = (.5*T*v) * ( pow( (T/(Adisk*pow(v,2)*ro*0.5)) + 1, .5) + 1 )
    m.set_ydata( PropPow )

sRoot.on_changed(updateThrust)
sDp.on_changed(updateThrust)
sM.on_changed(updateThrust)
sCl.on_changed(updateThrust)  
sTaper.on_changed(updateThrust)
sSpan.on_changed(updateThrust)

def updateLift(val):

    Root = sRoot.val
    weight = sM.val * g
    Tip = Root*sTaper.val  
    S = (0.5 * sSpan.val) * (Root+Tip)
    Cl = sCl.val
    k.set_ydata(0.5*ro*pow(v,2)*S*Cl)
    j.set_ydata(  weight  )
    
    fig.canvas.draw_idle()

sRoot.on_changed(updateLift)
sM.on_changed(updateLift)
sCl.on_changed(updateLift)
sSpan.on_changed(updateLift)
sTaper.on_changed(updateLift)

def Output(val):
    
    Root = sRoot.val
    mass = sM.val
    Tip = Root*sTaper.val    
    S = (0.5 * sSpan.val) * (Root+Tip)
    Output.Span = sSpan.val
    Output.S = S
    AR = pow((sSpan.val),2)/S
    s = 1.0 - 2.0 * ( df/sSpan.val )    
    Q = 1 / ( u*s )
    P = K*Cd0
    e = 1 / (Q + P*pi*AR)
    Cl = sCl.val
    L = 0.5*ro*pow(v,2)*S*Cl
    Cdi = pow(L,2) * ( 1/( .25*pow(ro,2)*pow(v,4)*pow(S,2)*pi*e*AR)) 
    LD = Cl/(Cdi[0]+Cd0)
    Power = (1/eta_p)*(0.5*ro*pow(v,3)*S*(Cdi+Cd0) + (( pow(mass*g,2) )/(0.5*ro*v*S))*(1/(pi*e*AR)))
    T = (Power/v)
    Adisk = pi*pow((0.5*sDp.val),2)
    PropPow = (.5*T*v) * ( pow( (T/(Adisk*pow(v,2)*ro*0.5)) + 1, .5) + 1 )
    
    for x in range(0,1):
        x = 1
        if x == 1:
            cls()
            print " "
            print "---MKS OUTPUTS---"
            print "Mass                   = " + str(round(mass,2)) + " kg"
            print "Weight                 = " + str(round(mass*g,2)) + " N"
            print "Oswald Efficiency      = " + str(round(e,2))
            print "Root Chord Length      = " + str(round(sRoot.val,2)) + " m"
            print "Tip Chord Length       = " + str(round(Tip,2)) + " m"
            print "Taper Ratio            = " + str(round(sTaper.val,2))
            print "Span                   = " + str(round(sSpan.val,2)) + " m"
            print "Aspect Ratio           = " + str(round(AR,2))
            print "Lift to Drag Ratio     = " + str(round(LD,2))
            print "Prop Diameter          = " + str(round(sDp.val,2)) + " m"
            print "Minimum Required Power = " + str(round(min(PropPow),2)) + " Watts"
            print " "
            
            print "---IPS OUTPUTS---"
            print "Mass                   = " + str(round(kg2slug*mass,2)) + " slug"
            print "Weight                 = " + str(round(N2lbf*mass*g,2)) + " lb"
            print "Oswald Efficiency      = " + str(round(e,2))
            print "Root Chord Length      = " + str(round(m2in*sRoot.val,2)) + " in"
            print "Tip Chord Length       = " + str(round(m2in*Tip,2)) + " in"
            print "Taper Ratio            = " + str(round(sTaper.val,2))
            print "Span                   = " + str(round(m2in*sSpan.val,2)) + " in"
            print "Aspect Ratio           = " + str(round(AR,2))
            print "Lift to Drag Ratio     = " + str(round(LD,2))
            print "Prop Diameter          = " + str(round(m2in*sDp.val,2)) + " in"
            print "Minimum Required Power = " + str(round(watt2hp*min(PropPow),2)) + " hp"
            print " "
            x = 0
        else:
            cls() 

    return Output.Span, Output.S     
    
sRoot.on_changed(Output)
sDp.on_changed(Output)
sM.on_changed(Output)
sCl.on_changed(Output)  
sTaper.on_changed(Output)
sSpan.on_changed(Output)
