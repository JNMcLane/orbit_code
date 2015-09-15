"""
Written by Jacob N. McLane
Last Modified Spet. 14th, 2015

Code for doing basic analysis of a two body system composed of a stellar primary
and a planetary mass secondary. Code could be easily adapted for use in a system
with two stellar components.

The program can also model 100 observations of the system by GAIA over 5 years
taking into account stellar motion, stellar motion+parallax, and stellar motion+
parallax+proper motion. The gaia models must be activated using the 'gaia' keyword.

*********Synatax*********

python orbit.py 'filename' year month day 'gaia'

    INPUTS

    Filename:
    The filename should be a file that list the following parameters of the system
    in one column with the following syntax:

    SEMI-MAJOR AXIS                                  in AU
    ECCENTRICITY
    INCLINATION                                        in degrees
    LONGITUDE OF ASCENDING NODE       in degrees
    ARGUMENT OF PERIAPSIS                      in degrees
    TIME OF PERIAPSIS PASSAGE                  in JD
    PRIMARY MASS                                      in solar masses
    SECONDARY MASS                                in Jupiter masses
    DISTANCE TO TARGET                          in parsecs
    RADIUS OF PRIMARY                             in solar radii
    RADIUS OF SECONDARY                       in Jupiter radii
    *STELLAR MAGNITUDE                           in the V-band
    *RA OF STAR                                          in decimal hours
    *DEC OF STAR                                        in decimal degrees
    *PARALLAX MOTION OF STAR                in mas
    *PROPER MOTION IN RA OF STAR           in mas
    *PROPER MOTION IN DEC OF STAR         in mas
_______________________________________________________________________
    **Parameters marked with a '*' are only needed if doing GAIA modeling.

    YEAR, MONTH, DAY

    What you would think, denotes date you wish to start the simulations on.

    OPTIONAL INPUT

    'gaia'

    Triggers gaia modeling. If not desired, use an empty sting (' ') instead.
===================================================

    OUTPUTS

    orbit_B_wrt_A.txt:

    X, Y, and Z positions of the secondary relative to the primary, along with the postion
    angle and angular seperation of the secondary relative to the primary at each timestep.
    ______________________________________________________________________
    orbit_A_wrt_com.txt:

    X, Y, and Z positions, postion angle and angular seperation of the primary relative to the
    center of mass at each time step. Also reports the RV of target.
    ______________________________________________________________________
    orbit_B_wrt_com.txt:

    X, Y, and Z positions, postion angle and angular seperation of the secondary relative to the
    center of mass at each time step. Also reports the RV of target.
    _____________________________________________________________________
    transit_times.txt:

    Reconds primary and secondary transits of the primary by the secondary.
    _____________________________________________________________________

    ******   CODE CITED  **********************************
    This program makes use of the coords.py written by Erin Sheldon, avaliable at:
    https://code.google.com/p/esutil/source/browse/trunk/esutil/coords.py?spec=svn414&r=414

    coords.py is used to transform from celestial to ecliptic coordinate systems

    This program also uses jdutil.py written by Matt Davis, avaliable at:
    https://gist.github.com/jiffyclub/1294443

    jdutil.py is used to turn calander dates into JD

    *******************************************************************
    *  THIS PROGRAM WAS WRITTEN FOR USE WITH PYTHON 3.4,    *
    *   COMPATIBILITY WITH OTHER VERSIONS NOT GUARANTEED   * 
    ********************************************************************
    
"""



###  Import Packages Used  ###
from math import *
from jdutil import *       #program to get JD from calander date
from coords import *    #program to transform from celestial to ecliptic coordinates
import sys
from random import uniform as rand
import matplotlib.pyplot as plt
##########################

###  Read in arguments from command line   ###
infile=sys.argv[1]             #Read in file with star/planet parameters
jd1=float(sys.argv[2])      #read in desired start date
jd2=float(sys.argv[3])
jd3=float(sys.argv[4])
t=date_to_jd(jd1,jd2,jd3) #calculate starting t in JD

###  Pull out values from parameter file
with open(infile) as file:
     param=file.readlines()

a=float(param[0])		    #semi-major axis in AU
ec=float(param[1])    		#eccentricity
i=float(param[2])     	    #inclination
om=float(param[3])         #longitude of ascending node
w=float(param[4])	    	#argument of periapsis in degrees
t0=float(param[5])   	    #time of periapsis passage in jd
m1=float(param[6])  		#mass primary in solar masses
m2=float(param[7])         #mass secondary in jupiter masses
D=float(param[8])            #distance in parsecs
R1=float(param[9])           #radius of primary in solar radii
R2=float(param[10])         #radius of secondary in Jupiter radii

G=6.67384e-11             #Gravitation constant

au=a
a=a*149597870700      #convert Au to M
i=radians(i)                    #convert i from degrees to radians
w=radians(w)                 #convert w from degrees to radians
om=radians(om)            #conver om from degrees to radians
m1=m1*1.989e30         #convert mass to kg
m2=m2*1.898e27         #convert mass to kg
R1=R1*0.0046491         #convert stellar radius to AU
R2=R2*69911/149597871 #convert jupiter radius to AU
################################################


###  g(E) and dg(E) are used in calculations of the eccentric anomaly

def g(E):
    return E-ec*sin(E)-M

def dg(E):
    return 1-ec*cos(E)
######################################################

#T=2*pi*sqrt(pow(a,3)/(G*(m1+m2)))/86400   #period in days
T=111.43637

### Open output logs ###
log1=open('orbit_B_wrt_A.txt','w')
log2=open('orbit_A_wrt_com.txt','w')
log3=open('orbit_B_wrt_com.txt','w')
log4=open('transit_times.txt','w')

for k in range(15200):
    M=(2*pi*((t-t0)/T))%(2*pi) #Mean Anomaly
    #E=M #Use of M for eccentricities < 0.8
    E=pi  #Use of pi is suggested over M for eccentricities >= 0.8

###  Loop to calculate E  ###
    for j in range(500):
        Eold=E
        E=E-g(E)/dg(E) #Calculate Eccentric anomaly
        if abs(g(E)/dg(E)) < 1e-5:
            break
    
    f=2*atan2(sqrt(1+ec)*sin(E/2),sqrt(1-ec)*cos(E/2)) #true anomaly
    #f=2*atan(sqrt((1+ec)/(1-ec))*tan(E/2))  #true anomaly
    r=au*(1-pow(ec,2))/(1+ec*cos(f))  #radius, will have same units as a

    ###########################

    #Do Relative Astrometry

    x=r*cos(f)
    y=r*sin(f)
    z=0

    #Do Absolute Astrometry of B wrt A
    
    X=r*(cos(om)*cos(w+f)-sin(om)*sin(w+f)*cos(i))
    Y=r*(sin(om)*cos(w+f)+cos(om)*sin(w+f)*cos(i))
    Z=r*sin(w+f)*sin(i)

###  Find position angle of B wrt A
    if X <= 0:
        PA=abs(degrees(atan2(X,Y)))

    else :
        PA =360-degrees(atan2(X,Y))

###  Projected seperation
    PS=tan((r/D))*1000 #in mas

    print(t,X,Y,Z,PA,PS, file=log1)

    #################
    
    rs=m2/(m1+m2)*r #DIstance of star from com
    rp=m1/(m1+m2)*r #DIstance of planet from com
    n=(2*pi)/(T*86400) #DIstance of star from com
    ast=m2/(m1+m2)*a #semi-major axis of star
    ap=m1/(m1+m2)*a #semi-major axis of planet 
    
    #Do absolute astromery and RV of A wrt com
    Xs=-rs*(cos(om)*cos(w+(f))-sin(om)*sin(w+(f))*cos(i))
    Ys=-rs*(sin(om)*cos(w+(f))+cos(om)*sin(w+(f))*cos(i))
    Zs=-rs*sin(w+(-f))*sin(i)

    if PA >= 180:
        PAs=PA-180

    else :
        PAs = PA+180
        
    PSs=tan((rs/D))*1000

    Ks=(m2/(m1+m2))*((n*a*sin(i))/sqrt(1-pow(ec,2)))
    rvs=Ks*(cos(w+(f))+ec*cos(w))

    print(t,Xs,Ys,Zs,PAs,PSs,rvs, file=log2)

    #Do absolute astromery and RV of B wrt com
    Xp=rp*(cos(om)*cos(w+f)-sin(om)*sin(w+f)*cos(i))
    Yp=rp*(sin(om)*cos(w+f)+cos(om)*sin(w+f)*cos(i))
    Zp=rp*sin(w+f)*sin(i)

    PAp=PA

    PSp=tan((rp/D))*1000

    #Kp=(m1/(m1+m2))*((n*a*sin(i))/sqrt(1-pow(ec,2)))
    #rvp=Kp*(cos(w+f)+ec*cos(w))
    rvp=(-rvs*m1)/m2
    print(t,Xp,Yp,Zp,PAp,PSp,rvp, file=log3)

    
    #Check for transits
    sep=(R1+R2)
    dist=sqrt(pow((abs(Xp)+abs(Xs)),2)+pow((abs(Yp)+abs(Ys)),2))
    if dist < sep:
        if Zp < 0:
            print(0,t, file=log4)  #Write in Secondary Transits
        if Zp > 0:
            print(t,0, file=log4)  #Write in Primary Transits
            
    
    t=t+0.01

### Close output logs ###
log1.close()
log2.close()
log3.close()
log4.close()
#####################

### GAIA Astrometry Simulation ###
if sys.argv[5] == 'gaia':  #only run if keyword is found

    log=open('gaia_motion.txt','w') #open output logs


### Create empty vectors to write values to
    jd=[]
    x1=[]
    y1=[]
    x2=[]
    y2=[]
    x3=[]
    y3=[]
    er1=[]
    er2=[]
    er3=[]
#####

    for k in range(100):   #Run for 100 points over 5 years
        tn=rand(t,t+1825)
        M=2*pi*((tn-t0)/T)
        E=M
        mag=float(param[11])  #Get magnitude of star from input file

### Assign GAIA Errors depending upon mag of star  ###
### All errors are in mas  #########################
        if mag <= 13:
            upos=6e-3
            upm=5e-3
            upar=8e-3
        elif mag <=14:
            upos=10e-3
            upm=7e-3
            upar=13e-3
        elif mag <=15:
            upos=16e-3
            upm=11e-3
            upar=21e-3
        elif mag <=16:
            upos=25e-3
            upm=18e-3
            upar=34e-3
        elif mag <=17:
            upos=55e-3
            upm=30e-3
            upar=40e-3
        elif mag <=18:
            upos=70e-3
            upm=50e-3
            upar=90e-3
        elif mag <=19:
            upos=115e-3
            upm=80e-3
            upar=155e-3
        else:
            upos=275e-6
            upm=145e-6
            upar=205e-6
#########################################
            
        for j in range(100):  #fit for E
            Eold=E
            E=E-g(E)/dg(E)

        f=2*atan2(sqrt(1+ec)*sin(E/2),sqrt(1-ec)*cos(E/2))
    
        r=au*(1-pow(ec,2))/(1+ec*cos(f))  #radius, will have same units as a
        rs=m2/(m1+m2)*r
    ###########################
    
    #Do Absolute Astrometry of star wrt com
    
        Xs=rs*(cos(om)*cos(w+(-f))-sin(om)*sin(w+(-f))*cos(i))
        Ys=rs*(sin(om)*cos(w+(-f))+cos(om)*sin(w+(-f))*cos(i))

        ra=float(param[12])              #ra of star in decimal hours
        dec=float(param[13])           #declination of star in decimal degrees
        parallax=float(param[14])    #parallax motion of star in mas
        pmra=float(param[15])        #proper motion in RA of star in mas
        pmdec=float(param[16])      #proper motion in Dec of star in mas
        
        ddec=tan((Ys/D))*1000      #angular distance in dec from com in mas
        dra=tan((Xs/D))*1000         #angular distance in ra from com in mas

        tyr=(tn-t)/365                               #years since start date
        delra=pmra*cos(radians(dec))*tyr  #proper motion in ra
        deldec=pmdec*tyr                         #proper motion in dec
        lon,lat=eq2ec(ra,dec)                     #convert to ecliptic coordinates
        lon,lat=float(lon),float(lat)

###  Apply parallax as a sinusoid wit amplitude +/- 1/2 of parallax over the year
        dparallaxt=(tyr-int(tyr))
        dparallax=sin(2*pi*dparallaxt)*parallax  
################################################################

        lon=lon+dparallax/(3600*1e3)
        ra2,dec2=ec2eq(lon+(5/(3600*1e3)),lat) #Convert back to celestial coordinates
        ra2,dec2=float(ra2),float(dec2)

###  Find shifts in RA and Dec due to parallax  ###
        dparra=(ra2-ra)*3600/15*1e3  
        dpardec=(dec2-dec)*3600*1e3
########################################
            
        print(tn,dra,ddec,dra+dparra,ddec+dpardec,dra+dparra+delra,ddec+dpardec+deldec,upos,upos+upar,upos+upm+upar, file=log)
#Write output file

#Write to vectors for plotting
        jd.append(tn-t)
        x1.append(dra)
        y1.append(ddec)
        x2.append(dra+dparra)
        y2.append(ddec+dpardec)
        x3.append(dra+dparra+delra)
        y3.append(ddec+dpardec+deldec)
        er1.append(5) #upos
        er2.append(upos+upar)
        er3.append(upos+upm+upar)

    log.close() #Close output log

###  Plot up results  ###
    plt.subplot(3,1,1)
    plt.plot(jd,x1,'r.')
    plt.title('Star Motion Only')
    plt.xerr=er1
    plt.yerr=er1
    plt.xlabel('Time (JD-'+str(t)+')')
    plt.ylabel('Change in RA (mas)')

    plt.subplot(3,1,2)
    plt.plot(jd,y1,'r.')
    plt.xerr=er1
    plt.yerr=er1
    plt.xlabel('Time (JD-'+str(t)+')')
    plt.ylabel('Change in Dec (mas)')

    plt.subplot(3,1,3)
    plt.plot(x1,y1,'r.')
    plt.xerr=er1
    plt.yerr=er1
    plt.xlabel('Change in RA (mas)')
    plt.ylabel('Change in Dec (mas)')

    plt.subplots_adjust(hspace=0.7)
    plt.savefig('Gaia_star.png')
    plt.close()

    ###############
    
    plt.subplot(3,1,1)
    plt.plot(jd,x2,'b.')
    plt.title('Star Motion + Parallax')
    plt.xerr=er2
    plt.yerr=er2
    plt.xlabel('Time (JD-'+str(t)+')')
    plt.ylabel('Change in RA (mas)')

    plt.subplot(3,1,2)
    plt.plot(jd,y2,'b.')
    plt.xerr=er2
    plt.yerr=er2
    plt.xlabel('Time (JD-'+str(t)+')')
    plt.ylabel('Change in Dec (mas)')

    plt.subplot(3,1,3)
    plt.plot(x2,y2,'b.')
    plt.xerr=er2
    plt.yerr=er2
    plt.xlabel('Change in RA (mas)')
    plt.ylabel('Change in Dec (mas)')

    plt.subplots_adjust(hspace=0.7)
    plt.savefig('Gaia_star_parallax.png')
    plt.close()

     ###############
    
    plt.subplot(3,1,1)
    plt.plot(jd,x3,'g.')
    plt.title('Star Motion + Parallax + Proper Motion')
    plt.xerr=er3
    plt.yerr=er3
    plt.xlabel('Time (JD-'+str(t)+')')
    plt.ylabel('Change in RA (mas)')

    plt.subplot(3,1,2)
    plt.plot(jd,y3,'g.')
    plt.xerr=er3
    plt.yerr=er3
    plt.xlabel('Time (JD-'+str(t)+')')
    plt.ylabel('Change in Dec (mas)')

    plt.subplot(3,1,3)
    plt.plot(x3,y3,'g.')
    plt.xerr=er3
    plt.yerr=er3
    plt.xlabel('Change in RA (mas)')
    plt.ylabel('Change in Dec (mas)')

    plt.subplots_adjust(hspace=0.7)
    plt.savefig('Gaia_star_parallax_pm.png')
