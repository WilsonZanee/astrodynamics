# MAE 469 - Introduction to Astrodynamics
# Makenzie Nicole Karakula
# Project Attempt
# =======================================
# ---- IMPORT LIBRARIES ----
import math
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt

# ---- FUNCTIONS ----

# ---- UNIT CONVERSIONS ---- 
def daysTOhours(days):
    # function DAYSTOHOURS : Convert days to hours.
    hours = days * (24 / 1)
    return hours

def hoursTOmin(hours):
    # function HOURSTOMIN : Convert hours to minutes.
    min = hours * (60 /1)
    return min

def minutesTOseconds(min):
    # function MINUTESTOSECONDS : Convert minutes to seconds.
    seconds = min * (60 / 1)
    return seconds

# ---- CALCULATE TIME PAST SINCE J2000 EPOCH ----
def TOFofInterest(M,D,Y,HR,MIN,RETURN):
    # function TOFOFINTEREST : Calculate time past to a certain date since J2000 Epoch.
    # ---- J2000 EPOCH PARAMETERS : January 1, 2000 at 11:58 AM ----
    J2000epochY = 2000 # Year 2000
    J2000epochHR = 11 # Military Time Hour 11 of 24
    J2000epochMIN = 58 # Minute 58 of 60
    # ---- NUMBER OF DAYS PER MONTH ----
    Month_standard = [31,28,31,30,31,30,31,31,30,31,30,31]
    Month_leap = [31,29,31,30,31,30,31,31,30,31,30,31] 
    # ---- EXCLUDE FULL CURRENT DAY ----
    D = D - 1 # Do not include the full day in current calculation, the HR and MIN represents this day.
    # ---- MONTH REFERENCE FOR NUMBER OF DAYS ----
    if M == "Jan":
        monthAdd = 0
    elif M == "Feb":
        monthAdd = 1
    elif M == "Mar":
        monthAdd = 2
    elif M == "Apr":
        monthAdd = 3
    elif M == "May":
        monthAdd = 4
    elif M == "Jun":
        monthAdd = 5
    elif M == "Jul":
        monthAdd = 6
    elif M == "Aug":
        monthAdd = 7
    elif M == "Sep":
        monthAdd = 8
    elif M == "Oct":
        monthAdd = 9
    elif M == "Nov":
        monthAdd = 10
    elif M == "Dec":
        monthAdd = 11
    # ---- CALCULATE DAYS ----
    k = J2000epochY # Set iteration to start at year 2000.
    totalDays = 0 # Set the zero start date as of Jan. 1, 2000.
    while k < Y: # Account for all full years leading up to the present year.
        if math.remainder((k-2000),4) == 0 or k == 2000: # Add total days in a full leap year.
            totalDays = totalDays + np.sum(Month_leap)
        else: # Add total days in a standard year.
            totalDays = totalDays + np.sum(Month_standard)
        k += 1 # Iterate to next year.
    while k == Y: # Account for the days so far in the current year.
        if math.remainder((k-2000),4) == 0 or k == 2000: # Add days in leap year thus far.
            if monthAdd == 0 : # January (Do not sum days, they are set.)
                daysThisYear = D 
            else: # Sum days until current month (range excludes the latter boundary) and add the set days of the current month.
                daysThisYear = np.sum(Month_leap[0:monthAdd]) + D
        else:  # Add days in standard year thus far.
            if monthAdd == 0 : # January (Do not sum days, they are set.)
                daysThisYear = D
            else: # Sum days until current month (range excludes the latter boundary) and add the set days of the current month.
                daysThisYear = np.sum(Month_standard[0:monthAdd]) + D
        k += 1 # Exit the final year.
    # ---- SUMMATIONS AND ACCOUNTING FOR J2000 EPOCH ----
    totalDays =  totalDays + daysThisYear # Sum the days from all previous years and the current year thus far.
    totalHours = daysTOhours(totalDays) + HR - J2000epochHR # Sum hours from previous calculation, add the hours in the present day, subtract the hours from the start of epoch (did not start at 00:00).
    totalMinutes = hoursTOmin(totalHours) + MIN - J2000epochMIN # Sum minutes from previous calculation, add the minutes in the present day, subtract the minutes from the start of epoch (did not start at 00:00).
    totalSeconds = minutesTOseconds(totalMinutes) # Convert total sum minutes to seconds.
    if RETURN == "seconds": # User gets return value of total seconds.
        return totalSeconds
    elif RETURN == "TUsun": # User gets return value of total time units TU with respect to heliocentric reference frame.
        TUsun = totalSeconds * (1/3600) * (1/24) * (1 / 58.13) # Convert seconds to hours, hours to days, days to heliocentric TU
        return TUsun

def OrbitalElementsToRV(a,e,i,raan,omega,theta,mu,raddeg,table,distance,time):
    # function ORBITALELEMENTSTORV : Take orbital elements and find position and velocity vectors.
    if raddeg == "degree": # Convert degrees to radians.
        i = math.radians(i)
        raan = math.radians(raan)
        omega = math.radians(omega)
        theta = math.radians(theta)
    elif raddeg == "radian":
        pass # Keep values in radians
    p = a*(1-e**2)
    r_scalar = (a*(1-e**2))/(1+e*math.cos(theta))
    r_vec_PQW = [r_scalar*math.cos(theta), r_scalar*math.sin(theta),0]
    v_vec_PQW = [(math.sqrt(mu/p))*-math.sin(theta), (math.sqrt(mu/p))*(e+math.cos(theta)), 0]
    rotationVec = [[math.cos(raan)*math.cos(omega)-math.sin(raan)*math.sin(omega)*math.cos(i), -1*math.cos(raan)*math.sin(omega)-math.sin(raan)*math.cos(omega)*math.cos(i), math.sin(raan)*math.sin(i)],[math.sin(raan)*math.cos(omega)+math.cos(raan)*math.sin(omega)*math.cos(i), -1*math.sin(raan)*math.sin(omega)+math.cos(raan)*math.cos(omega)*math.cos(i), -1*math.cos(raan)*math.cos(i)],[math.sin(omega)*math.sin(i), math.cos(omega)*math.sin(i), math.cos(i)]]
    r_IJK = [rotationVec[0][0]*r_vec_PQW[0]+rotationVec[0][1]*r_vec_PQW[1]+rotationVec[0][2]*r_vec_PQW[2], rotationVec[1][0]*r_vec_PQW[0]+rotationVec[1][1]*r_vec_PQW[1]+rotationVec[1][2]*r_vec_PQW[2], rotationVec[2][0]*r_vec_PQW[0]+rotationVec[2][1]*r_vec_PQW[1]+rotationVec[2][2]*r_vec_PQW[2]]
    v_IJK = [rotationVec[0][0]*v_vec_PQW[0]+rotationVec[0][1]*v_vec_PQW[1]+rotationVec[0][2]*v_vec_PQW[2], rotationVec[1][0]*v_vec_PQW[0]+rotationVec[1][1]*v_vec_PQW[1]+rotationVec[1][2]*v_vec_PQW[2], rotationVec[2][0]*v_vec_PQW[0]+rotationVec[2][1]*v_vec_PQW[1]+rotationVec[2][2]*v_vec_PQW[2]]
    # Distance units must be AU, DU, or traditional unit [km, m].
    # Time units must be TU, or traditional unit [s].
    if table == "full":
        # All calculated elements will be printed and values are retured.
        print('r_pqw = ',np.around(r_vec_PQW,4),'%s' % (distance))
        print('v_pqw = ',np.around(v_vec_PQW,4),'%s/%s' % (distance,time))
        print('Rotation Matrix =\n',np.around(rotationVec,4))
        print('r_ijk = ',np.around(r_IJK,4),'%s' % (distance))
        print('v_ijk = ',np.around(v_IJK,4),'%s/%s' % (distance,time))
        return r_IJK, v_IJK 
    elif table == "elements":
        # Only values are returned.
        return r_IJK, v_IJK 
    
def universalTOF_SCZ(r0,v0,tof0,mu):
    # function UNIVERSALTOF_SCZ : Calculate r2 and v2 from initial vectors and time of flight (zSC method).
    # ---- MAGNITUDES OF POSITION AND VELOCITY VECTORS ----
    r0_scalar = math.sqrt(r0[0]**2 + r0[1]**2 + r0[2]**2)
    v0_scalar = math.sqrt(v0[0]**2 + v0[1]**2 + v0[2]**2)
    orbitalEnergy = ((v0_scalar**2)/2) - (mu/r0_scalar) # Orbital energy of current orbit
    a = -mu/(2*orbitalEnergy) # Semi-major axis
    rdotv = r0[0]*v0[0] + r0[1]*v0[1] + r0[2]*v0[2]
    tol = 1*10**(-7) # Iteration error tolerance
    if a > 0: # Guess for elliptical case.
        x = (math.sqrt(mu)*tof0)/a
    else: 
        x = 1
        pass # Determine x guess
    z = (x**2) / a
    while True:
        z = (x**2) / a
        if z > 0.00001 : # Ellipse
            S = (math.sqrt(z) - math.sin(math.sqrt(z))) / math.sqrt(z**3)
            C = (1 - math.cos(math.sqrt(z))) / z
        elif z < -0.00001 : # Hyperbola
            S = (math.sinh(math.sqrt(-z))-math.sqrt(-z)) / math.sqrt((-z)**3)
            C = (1 - math.cosh(math.sqrt(-z))) / z
        else:
            print("Orbit is a parabola. This will likely cause overflow errors in calculation.")
            break
        t = (1/math.sqrt(mu))*((x**3)*S + (rdotv/math.sqrt(mu))*(x**2)*C + r0_scalar*x*(1-z*S))
        r = (x**2)*C + (rdotv/math.sqrt(mu))*x*(1-z*S) + r0_scalar*(1 - z*C)
        dtdx  = r / math.sqrt(mu)
        xnew = x + (tof0-t)/dtdx
        if abs(tof0-t) <= tol:
            f = 1 - ((x**2)/r0_scalar)*C
            g = t - ((x**3)/math.sqrt(mu))*S
            fdot = ((math.sqrt(mu)*x)/(r0_scalar*r))*(z*S-1)
            gdot = 1 - ((x**2)/r)*C
            r2 = [f*r0[0]+g*v0[0], f*r0[1]+g*v0[1], f*r0[2]+g*v0[2]]
            v2 = [fdot*r0[0]+gdot*v0[0], fdot*r0[1]+gdot*v0[1], fdot*r0[2]+gdot*v0[2]]
            #print('x = ',np.around(x,4))
            #print('f = ', np.around(f,4))
            #print('g = ', np.around(g,4))
            #print('fdot = ', np.around(fdot,4))
            #print('gdot = ', np.around(gdot,4))            
            #print("r2 = ",np.around(r2,4))
            #print("v2 = ",np.around(v2,4))
            break
        x = xnew
    return r2,v2

def rvToOrbitalElements(r_vec,v_vec,mu,table,distance,time):
    # function RVTOORBITALELEMENTS : Take position and velocity vectors to calculate the orbital elements.
    # ---- Establish position vector r_ijk and |r| ----
    r_scalar = math.sqrt(r_vec[0]**2 + r_vec[1]**2 + r_vec[2]**2) # Magnitude of the position vector
    # ---- Establish velocity vector v_ijk and |v| ----
    v_scalar = math.sqrt(v_vec[0]**2 + v_vec[1]**2 + v_vec[2]**2)
    # ---- Establish angular momentum : h = r x v (or h = rp*vp) ----
    h_vec = [r_vec[1]*v_vec[2]-r_vec[2]*v_vec[1],-1*(r_vec[0]*v_vec[2]-r_vec[2]*v_vec[0]),r_vec[0]*v_vec[1]-r_vec[1]*v_vec[0]]
    h_scalar = math.sqrt(h_vec[0]**2 + h_vec[1]**2 + h_vec[2]**2)
    # ---- Establish the node vector : n = k x h ----
    n_vec = [-1*h_vec[1], -1*-1*h_vec[0], 0]
    n_scalar = math.sqrt(n_vec[0]**2 + n_vec[1]**2 + n_vec[2]**2)
    r_dot_v = r_vec[0]*v_vec[0] + r_vec[1]*v_vec[1] + r_vec[2]*v_vec[2]
    # ---- Establish eccentricity e vector ----
    e_vec1 = (1/mu)*(v_scalar**2 - (mu/r_scalar))
    e_vec1_vec = [e_vec1*r_vec[0], e_vec1*r_vec[1], e_vec1*r_vec[2]]
    e_vec2 = (1/mu)*((r_dot_v))
    e_vec2_vec = [e_vec2*v_vec[0], e_vec2*v_vec[1], e_vec2*v_vec[2]]
    e_vec = [e_vec1_vec[0]-e_vec2_vec[0], e_vec1_vec[1]-e_vec2_vec[1], e_vec1_vec[2]-e_vec2_vec[2]]
    e_scalar = math.sqrt(e_vec[0]**2 + e_vec[1]**2 + e_vec[2]**2)
    p = h_scalar**2 / mu
    a = p / (1 - e_scalar**2)
    e_dot_r = r_vec[0]*e_vec[0] + r_vec[1]*e_vec[1] + r_vec[2]*e_vec[2]
    # ---- Calculate orbital elements -----
    inclination = math.degrees(math.acos(h_vec[2]/h_scalar))
    n_dot_e = n_vec[0]*e_vec[0] + n_vec[1]*e_vec[1] + n_vec[2]*e_vec[2]
    if n_scalar == 0:
        pass
    else:
        longAscendNode = math.degrees(math.acos(n_vec[0]/n_scalar))
        ArgOfPeri = math.degrees(math.acos(n_dot_e/(n_scalar*e_scalar)))
    TrueAnomaly = math.degrees(math.acos(e_dot_r/(e_scalar*r_scalar)))
    # ---- Run sanity checks and adjust for orbital elements ----
    if inclination >= 0 and inclination <180:
        pass
    else:
        print("INCLINATION ERROR")
    if n_vec[1] > 0:
        pass
    else:
        if n_scalar == 0 :
            pass
        else: 
            longAscendNode = 360 - longAscendNode
    if e_vec[2] > 0:
        pass
    else: 
        if n_scalar == 0 :
            pass
        else: 
            ArgOfPeri = 360 - ArgOfPeri
    if r_dot_v > 0:
        pass
    else:
        TrueAnomaly = 360 - TrueAnomaly
    # Distance units must be AU, DU, or traditional unit [km, m].
    # Time units must be TU, or traditional unit [s].
    if table == "full":
        # Print all calculated elements as well as the resultant orbital elements.
        print("Calculations")
        print("------------")
        print(' r_ijk  = ',np.around(r_vec,8), '%s' % (distance))
        print('|r_ijk| = ',np.around(r_scalar,4), '%s' % (distance))
        print(' v_ijk  = ',np.around(v_vec,8),'%s/%s' % (distance,time))
        print('|v_ijk| = ',np.around(v_scalar,4),'%s/%s' % (distance,time))
        print('   h    = ',np.around(h_vec,4),'%s^2/%s' % (distance,time))
        print('  |h|   = ',np.around(h_scalar,4),'%s^2/%s' % (distance,time))
        print('   n    = ',np.around(n_vec,4),'%s^2/%s' % (distance,time))
        print('  |n|   = ',np.around(n_scalar,4),'%s^2/%s' % (distance,time))
        print('r dot v = ',np.around(r_dot_v,4) ,'%s^2/%s' % (distance,time))
        print('   e    = ',np.around(e_vec,4))
        print('e dot r = ',np.around(e_dot_r,4),'%s' % (distance))
        print('n dot e = ',np.around(n_dot_e,4),'%s^2/%s' % (distance,time))
        print("\nOrbital Objects")
        print("---------------")
        print("|e| = %.3f" % e_scalar)
        print("a = %.3f %s" % (a,distance))
        print("i = %.3f degrees" % (inclination))
        print("θ = %.3f degrees" % (TrueAnomaly))
        if n_scalar == 0:
            print("Ω = DNE")
            print("ω = DNE")
            return e_scalar, a, inclination, TrueAnomaly
        else:
            print("Ω = %.3f degrees" % (longAscendNode))
            print("ω = %.3f degrees" % (ArgOfPeri))
            return e_scalar, a, inclination, TrueAnomaly, longAscendNode, ArgOfPeri
    elif table == "elements":
        # Do not print the resultant calculations. Only return the values of interest.
        if n_scalar == 0:
            return e_scalar, a, inclination, TrueAnomaly
        else:
            return e_scalar, a, inclination, TrueAnomaly, longAscendNode, ArgOfPeri
        
# ---- ITERATE TO FIND LAUNCH DATE ----
def gauss(r1,r2,tof,mu,zguess,SHORTLONG):
    # function GAUSS: Calculate v1 and v2 using Gauss' problem methods
    tol = 1*10**(-7)
    # -- Normalize the position vectors -- 
    r1_scalar = math.sqrt(r1[0]**2 + r1[1]**2 + r1[2]**2)
    r2_scalar = math.sqrt(r2[0]**2 + r2[1]**2 + r2[2]**2)
    #print('|r1| = %.3f' % r1_scalar)
    #print('|r2| = %.3f' % r2_scalar)
    # -- Calculate angles of long and short methods -- 
    r1dotr2 = r1[0]*r2[0] + r1[1]*r2[1] + r1[2]*r2[2]
    thetaSHORT = math.acos(r1dotr2/(r1_scalar*r2_scalar))
    thetaLONG = 2*math.pi - thetaSHORT
    print(thetaSHORT)
    print(thetaLONG)
    #print('thetaShort = %.5f rad' % thetaSHORT)
    #print('thetaLong = %.5f rad' % thetaLONG)
    # -- Calculate area for short and long methods -- 
    ASHORT =(math.sqrt(r1_scalar*r2_scalar)*math.sin(thetaSHORT))/(math.sqrt(1-math.cos(thetaSHORT)))
    ALONG = (math.sqrt(r1_scalar*r2_scalar)*math.sin(thetaLONG))/(math.sqrt(1-math.cos(thetaLONG)))
    #print('A Short = %.4f' % ASHORT)
    #print('A Long  = %.4f' % ALONG)
    # -- Take user Z guess to be her new z --
    Zn = zguess
    # -- Approximate S and C from Z -- 
    if SHORTLONG == 'S':
        #print('------------\nShort Way Chosen.\n------------')
        A = ASHORT
    else:
        #print('------------\nLong Way Chosen.\n------------')
        A = ALONG   
    #print('zn         |    S(z)    |    C(z)    |      y       |      x      |     tn     |    TOF-t   |    dS/dZ    |    dC/dZ    |   dt/dZ  |    Zn+1')
    #print('-------------------------------------------------------------------------------------------------------------------------------------------------') 
    while True:
        if Zn > 0 :
            Sz = (math.sqrt(Zn)-math.sin(math.sqrt(Zn))) / (math.sqrt(Zn**3))
            Cz = (1 - math.cos(math.sqrt(Zn))) / Zn
        elif Zn < 0 :
            Sz = (math.sinh(math.sqrt(-Zn))-math.sqrt(-Zn)) / math.sqrt((-Zn)**3)
            Cz = (1-math.cosh(math.sqrt(-Zn))) / Zn
        else:
            Sz = (1/math.factorial(3)) - (Zn/math.factorial(5)) + ((Zn**2)/math.factorial(7)) - ((Zn**3)/math.factorial(9))
            Cz = (1/math.factorial(2)) - (Zn/math.factorial(4)) + ((Zn**2)/math.factorial(6)) - ((Zn**3)/math.factorial(8))
        y = r1_scalar + r2_scalar - A*(1-Zn*Sz)/math.sqrt(Cz)
        x = math.sqrt(y/Cz)
        tn = (1/math.sqrt(mu)) * ((x**3)*Sz + A*math.sqrt(y))
        dSdZ = (Cz-3*Sz) / (2*Zn)
        dCdZ = (1-Zn*Sz-2*Cz) / (2*Zn)
        dtdZ = (1/math.sqrt(mu)) * ((x**3)*(dSdZ-(3*Sz*dCdZ)/(2*Cz)) + (A/8)*((3*Sz*math.sqrt(y)/Cz)+(A/x)))
        Zn1 = Zn + (tof-tn)/dtdZ
        if abs(tof-tn) <= tol: #abs(Zn-Zn1) >= 0:
            #print('%.4f     |   %.4f   |   %.4f   |    %.4f    |    %.4f   |   %.4f   |   %.4f   |   %.4f   |   %.4f   |   %.4f | %.4f' %(Zn,Sz,Cz,y,x,tn,tof-tn,dSdZ,dCdZ,dtdZ,Zn1))
            f = 1 - y/r1_scalar
            g = A *math.sqrt(y/mu)
            fdot = (-1*math.sqrt(mu)*x)*(1-Zn*Sz)/(r1_scalar*r2_scalar)
            gdot = 1 - y/r2_scalar
            v1 = [(1/g)*r2[0]-(f/g)*r1[0], (1/g)*r2[1]-(f/g)*r1[1], (1/g)*r2[2]-(f/g)*r1[2]]
            v2 = [(gdot/g)*r2[0]-(1/g)*r1[0], (gdot/g)*r2[1]-(1/g)*r1[1], (gdot/g)*r2[2]-(1/g)*r1[2]]
            #print('-----\nf = %.5f' % f)
            #print('g = %.5f' % g)
            #print('fdot = %.5f' % fdot)
            #print('gdot = %.5f\n-----' % gdot)
            #print(v1)
            #print(v2)
            break
        #print('%.4f     |   %.4f   |   %.4f   |    %.4f    |    %.4f   |   %.4f   |   %.4f   |   %.4f   |   %.4f   |   %.4f | %.4f' %(Zn,Sz,Cz,y,x,tn,tof-tn,dSdZ,dCdZ,dtdZ,Zn1))
        Zn = Zn1
    return v1, v2

def plotSynodicPeriod():
    # function PLOTSYNODICPERIOD : Plot the synodic period of both Earth and Mars.
    format_data = "%m/%d/%y %H:%M"
    locationOfSun = [0, 0, 0] # [AU] Sun Position Vector in Heliocentric Frame
    mu = 1 # [AU3/TU2] Heliocentric gravitational parameter
    # ---- PLOT LOCATION OF THE SUN ----
    plt.plot(locationOfSun[0],locationOfSun[1],'y',marker=".",markersize=30)
    # ---- PLANETARY ORBITAL OBJECTS ----
    # -- EARTH --
    aE = 1.000000; eE = 0.01671; iE = 0.00005; raanE = -11.26064; omegaE = 114.20783; thetaE = -2.48284
    # -- MARS --
    aM = 1.523662; eM = 0.093412; iM = 1.85061; raanM = 49.57854; omegaM = 286.4623; thetaM = 19.41248
    # -- ENSURE RESULTS ARE FORWARDED FOR ORBITAL ELEMENTS --
    raddeg = "degree" # The 'degree' or 'radian' for parameters.
    table = "elements" # The 'full' prints all calculations and results, 'elements' only forwards results.
    distance = "AU" # The distance unit for table printing. Mainly used for 'full' prints.
    time = "TU" # The time unit for table printing. Mainly used for 'full' prints.
    # -- ESTABLISH POSITIONAL AND VELOCITY VECTORS AT MOMENT OF J2000 EPOCH --
    J2000Epoch = datetime.strptime('1/1/00 11:58',format_data) # Capture the time of J2000 Epoch
    [repochE,vepochE] = OrbitalElementsToRV(aE,eE,iE,raanE,omegaE,thetaE,mu,raddeg,table,distance,time)
    [repochM,vepochM] = OrbitalElementsToRV(aM,eM,iM,raanM,omegaM,thetaM,mu,raddeg,table,distance,time)
    day = 0 # [day] Iteration Reference for Syndoic Period Starting at J2000 Epoch
    SyndoicPeriod = 781 # [days] Approximately equal to the 2.14 year E/M syndoic period
    while True:
        TOF = day / 58.13 # [TUsun] Time of Flight since Epoch
        if day == 0:
            plt.plot(repochE[0],repochE[1],'b',marker=".",markersize=20)
            plt.plot(repochM[0],repochM[1],'r',marker=".",markersize=15)
            plt.legend(["Sun","Synodic Period of Earth","Synodic Period of Mars"],loc="best")
        else:
            [rE,vE] = universalTOF_SCZ(repochE,vepochE,TOF,mu) 
            [rM,vM] = universalTOF_SCZ(repochM,vepochM,TOF,mu)
            plt.plot(rE[0],rE[1],'b',marker=".",markersize=3)
            plt.plot(rM[0],rM[1],'r',marker=".",markersize=3)
        if day == SyndoicPeriod:
            break
        else:
            day += 1
    plt.title("Heliocentric Synodic Period of Earth and Mars")
    plt.xlabel("X-Position [AU]")
    plt.ylabel("Y-Position [AU]")
    plt.show()

def plotTOF(time_data):
    # function PLOTTOF : Plot the time of flight and compare the position of departure and arrival.
    format_data = "%m/%d/%y %H:%M"
    locationOfSun = [0, 0, 0] # [AU] Sun Position Vector in Heliocentric Frame
    mu = 1 # [AU3/TU2] Heliocentric gravitational parameter
    # ---- PLOT LOCATION OF THE SUN ----
    plt.plot(locationOfSun[0],locationOfSun[1],'y',marker=".",markersize=30)
    # ---- PLANETARY ORBITAL OBJECTS ----
    # -- EARTH --
    aE = 1.000000; eE = 0.01671; iE = 0.00005; raanE = -11.26064; omegaE = 114.20783; thetaE = -2.48284
    # -- MARS --
    aM = 1.523662; eM = 0.093412; iM = 1.85061; raanM = 49.57854; omegaM = 286.4623; thetaM = 19.41248
    # -- ENSURE RESULTS ARE FORWARDED FOR ORBITAL ELEMENTS --
    raddeg = "degree" # The 'degree' or 'radian' for parameters.
    table = "elements" # The 'full' prints all calculations and results, 'elements' only forwards results.
    distance = "AU" # The distance unit for table printing. Mainly used for 'full' prints.
    time = "TU" # The time unit for table printing. Mainly used for 'full' prints.
    # -- ESTABLISH POSITIONAL AND VELOCITY VECTORS AT MOMENT OF J2000 EPOCH --
    J2000Epoch = datetime.strptime('1/1/00 11:58',format_data) # Capture the time of J2000 Epoch
    [repochE,vepochE] = OrbitalElementsToRV(aE,eE,iE,raanE,omegaE,thetaE,mu,raddeg,table,distance,time)
    [repochM,vepochM] = OrbitalElementsToRV(aM,eM,iM,raanM,omegaM,thetaM,mu,raddeg,table,distance,time)
    departureTime = datetime.strptime(time_data,format_data)
    
    day = 0 # [day] Iteration Reference for Syndoic Period Starting at J2000 Epoch
    SyndoicPeriod = 781 # [days] Approximately equal to the 2.14 year E/M syndoic period
    while True:
        TOF = day / 58.13 # [TUsun] Time of Flight since Epoch
        if day == 0:
            plt.plot(repochE[0],repochE[1],'b',marker=".",markersize=20)
            plt.plot(repochM[0],repochM[1],'r',marker=".",markersize=15)
            plt.legend(["Sun","Synodic Period of Earth","Synodic Period of Mars"],loc="best")
        else:
            [rE,vE] = universalTOF_SCZ(repochE,vepochE,TOF,mu) 
            [rM,vM] = universalTOF_SCZ(repochM,vepochM,TOF,mu)
            plt.plot(rE[0],rE[1],'b',marker=".",markersize=3)
            plt.plot(rM[0],rM[1],'r',marker=".",markersize=3)
        if day == SyndoicPeriod:
            break
        else:
            day += 1
    plt.title("Heliocentric Synodic Period of Earth and Mars")
    plt.xlabel("X-Position [AU]")
    plt.ylabel("Y-Position [AU]")
    plt.show()

def gaussSHORTWAYcheck(r1,r2):
    # function GAUSSSHORTWAYCHECK : See if the 190 day transfer is reasonable.
    # -- Normalize the position vectors -- 
    r1_scalar = math.sqrt(r1[0]**2 + r1[1]**2 + r1[2]**2)
    r2_scalar = math.sqrt(r2[0]**2 + r2[1]**2 + r2[2]**2)
    #print('|r1| = %.3f' % r1_scalar)
    #print('|r2| = %.3f' % r2_scalar)
    # -- Calculate angles of long and short methods -- 
    r1dotr2 = r1[0]*r2[0] + r1[1]*r2[1] + r1[2]*r2[2]
    thetaSHORT = math.acos(r1dotr2/(r1_scalar*r2_scalar))
    return thetaSHORT
    
def reasonableLaunchSearch(time_data):
    # function LAUNCHDATE : Approximate the departure date from Earth to find 190 day SHORT transfer flight to Mars
    format_data = "%m/%d/%y %H:%M"
    tolerance = 1*10**(-3) # Iterative tolerance
    TUtimeOfFlight = 190/58.13 # [TUsun] 190-Day Transfer Converted
    # ---- PLANETARY ORBITAL OBJECTS ----
    mu = 1 # [AU3/TU2] Heliocentric gravitational parameter
    # -- EARTH --
    aE = 1.000000; eE = 0.01671; iE = 0.00005; raanE = -11.26064; omegaE = 114.20783; thetaE = -2.48284
    # -- MARS --
    aM = 1.523662; eM = 0.093412; iM = 1.85061; raanM = 49.57854; omegaM = 286.4623; thetaM = 19.41248
    # -- ENSURE RESULTS ARE FORWARDED FOR ORBITAL ELEMENTS --
    raddeg = "degree" # The 'degree' or 'radian' for parameters.
    table = "elements" # The 'full' prints all calculations and results, 'elements' only forwards results.
    distance = "AU" # The distance unit for table printing. Mainly used for 'full' prints.
    time = "TU" # The time unit for table printing. Mainly used for 'full' prints.
    # -- ESTABLISH POSITIONAL AND VELOCITY VECTORS AT MOMENT OF J2000 EPOCH --
    J2000Epoch = datetime.strptime('1/1/00 11:58',format_data) # Capture the time of J2000 Epoch
    [repochE,vepochE] = OrbitalElementsToRV(aE,eE,iE,raanE,omegaE,thetaE,mu,raddeg,table,distance,time)
    [repochM,vepochM] = OrbitalElementsToRV(aM,eM,iM,raanM,omegaM,thetaM,mu,raddeg,table,distance,time)
    # -- MAKE INITIAL ESTIMATE FOR DATE OF EARTH DEPARTURE --
    earth_departure_date = datetime.strptime(time_data, format_data) # Establish initial guess of Earth departure at 1/1/21 00:00
    failure_date = date = datetime.strptime("2/21/23 12:00", format_data)
    mars_arrival_date = earth_departure_date + timedelta(days=190) # Establish initial guess of Mars arrival after 190 day TOF
    i = 0 # Time start at midnight 00:00 (minutes passes counter)
    # -- ENTER INTERATIVE LOOP --
    
    while True:
        # -- CALCULATE INITIAL TIMES OF FLIGHT -- 
        earthTIMECHANGE = (earth_departure_date - J2000Epoch) # Amount of time passed since J2000 EPOCH to Earth Departure
        TOFearth = earthTIMECHANGE.total_seconds() / 58.13 # [TUsun] Time of Flight of Earth since Epoch
        marsTIMECHANGE = (mars_arrival_date - J2000Epoch) # Amount of time passed since J2000 EPOCH to Mars Arrival
        TOFmars = marsTIMECHANGE.total_seconds() / 58.13 # [TUsun] Time of Flight of Mars since Epoch
        # -- CALCULATE INITIAL ESTIMATES FOR DEPARTURE AND ARRIVAL VECTORS --
        [rEDepart,vEDepart] = universalTOF_SCZ(repochE,vepochE,TOFearth,mu) 
        [rMArrival,vMArrival] = universalTOF_SCZ(repochM,vepochM,TOFmars,mu)
        rE = math.sqrt(rEDepart[0]**2 + rEDepart[1]**2 + rEDepart[2]**2) # [AU] Magnitude of rEDepart
        rM = math.sqrt(rMArrival[0]**2 + rMArrival[1]**2 + rMArrival[2]**2) # [AU] Magnitude of rMArrival
        vE = math.sqrt(vEDepart[0]**2 + vEDepart[1]**2 + vEDepart[2]**2) # [AU/TU] Magnitude of vEDepart
        vM = math.sqrt(vMArrival[0]**2 + vMArrival[1]**2 + vMArrival[2]**2) # [AU/TU] Magnitude of vMArrival
        # -- IMPLEMENT GAUSS'S PROBLEM --
        thetaSHORT = gaussSHORTWAYcheck(rEDepart,rMArrival)

        if thetaSHORT < 1*math.pi/10:
            if earth_departure_date == failure_date:
                #print("Number of dates valid out of 781 : ",i)
                #print("PROGRAM ENDED. SYNODIC PERIOD PASSED.")
                break
            if i == 0:
                first_guess_departure = earth_departure_date
                first_guess_arrival = mars_arrival_date
            #print("A DATE WAS FOUND", earth_departure_date ," , ", mars_arrival_date," , ", thetaSHORT)
            earth_departure_date = earth_departure_date + timedelta(days=1)
            mars_arrival_date = earth_departure_date + timedelta(days=190)
            i = i+1
        else:
            if earth_departure_date == failure_date:
                #print("Number of dates valid out of 781 : ",i)
                #print("PROGRAM ENDED. SYNODIC PERIOD PASSED.")
                break
            earth_departure_date = earth_departure_date + timedelta(days=1)
            mars_arrival_date = earth_departure_date + timedelta(days=190)
    return first_guess_departure, first_guess_arrival
    
def launchDatePositions(departure_data, arrival_data):
    # function LAUNCHDATE : Approximate the departure date from Earth to find 190 day SHORT transfer flight to Mars
    format_data = "%m/%d/%y %H:%M"
    tolerance = 1*10**(-3) # Iterative tolerance
    TUtimeOfFlight = 190/58.13 # [TUsun] 190-Day Transfer Converted
    # ---- PLANETARY ORBITAL OBJECTS ----
    mu = 1 # [AU3/TU2] Heliocentric gravitational parameter
    # -- EARTH --
    aE = 1.000000; eE = 0.01671; iE = 0.00005; raanE = -11.26064; omegaE = 114.20783; thetaE = -2.48284
    # -- MARS --
    aM = 1.523662; eM = 0.093412; iM = 1.85061; raanM = 49.57854; omegaM = 286.4623; thetaM = 19.41248
    # -- ENSURE RESULTS ARE FORWARDED FOR ORBITAL ELEMENTS --
    raddeg = "degree" # The 'degree' or 'radian' for parameters.
    table = "elements" # The 'full' prints all calculations and results, 'elements' only forwards results.
    distance = "AU" # The distance unit for table printing. Mainly used for 'full' prints.
    time = "TU" # The time unit for table printing. Mainly used for 'full' prints.
    # -- ESTABLISH POSITIONAL AND VELOCITY VECTORS AT MOMENT OF J2000 EPOCH --
    J2000Epoch = datetime.strptime('1/1/00 11:58',format_data) # Capture the time of J2000 Epoch
    [repochE,vepochE] = OrbitalElementsToRV(aE,eE,iE,raanE,omegaE,thetaE,mu,raddeg,table,distance,time)
    [repochM,vepochM] = OrbitalElementsToRV(aM,eM,iM,raanM,omegaM,thetaM,mu,raddeg,table,distance,time)
    # -- MAKE INITIAL ESTIMATE FOR DATE OF EARTH DEPARTURE --
    earth_departure_date = datetime.strptime(departure_data, format_data) # Establish initial guess of Earth departure at 1/1/21 00:00
    mars_arrival_date = datetime.strptime(arrival_data, format_data) # Establish initial guess of Mars arrival after 190 day TOF
    # -- CALCULATE INITIAL TIMES OF FLIGHT -- 
    earthTIMECHANGE = (earth_departure_date - J2000Epoch) # Amount of time passed since J2000 EPOCH to Earth Departure
    TOFearth = earthTIMECHANGE.total_seconds() / 58.13 # [TUsun] Time of Flight of Earth since Epoch
    marsTIMECHANGE = (mars_arrival_date - J2000Epoch) # Amount of time passed since J2000 EPOCH to Mars Arrival
    TOFmars = marsTIMECHANGE.total_seconds() / 58.13 # [TUsun] Time of Flight of Mars since Epoch
    # -- CALCULATE INITIAL ESTIMATES FOR DEPARTURE AND ARRIVAL VECTORS --
    [rEDepart,vEDepart] = universalTOF_SCZ(repochE,vepochE,TOFearth,mu) 
    [rMArrival,vMArrival] = universalTOF_SCZ(repochM,vepochM,TOFmars,mu)
    return rEDepart, rMArrival

def launchDate(departure_data, arrival_data):
    # function LAUNCHDATE : Approximate the departure date from Earth to find 190 day SHORT transfer flight to Mars
    format_data = "%m/%d/%y %H:%M"
    tolerance = 1*10**(-3) # Iterative tolerance
    TUtimeOfFlight = 190/58.13 # [TUsun] 190-Day Transfer Converted
    # ---- PLANETARY ORBITAL OBJECTS ----
    mu = 1 # [AU3/TU2] Heliocentric gravitational parameter
    # -- EARTH --
    aE = 1.000000; eE = 0.01671; iE = 0.00005; raanE = -11.26064; omegaE = 114.20783; thetaE = -2.48284
    # -- MARS --
    aM = 1.523662; eM = 0.093412; iM = 1.85061; raanM = 49.57854; omegaM = 286.4623; thetaM = 19.41248
    # -- ENSURE RESULTS ARE FORWARDED FOR ORBITAL ELEMENTS --
    raddeg = "degree" # The 'degree' or 'radian' for parameters.
    table = "elements" # The 'full' prints all calculations and results, 'elements' only forwards results.
    distance = "AU" # The distance unit for table printing. Mainly used for 'full' prints.
    time = "TU" # The time unit for table printing. Mainly used for 'full' prints.
    # -- ESTABLISH POSITIONAL AND VELOCITY VECTORS AT MOMENT OF J2000 EPOCH --
    J2000Epoch = datetime.strptime('1/1/00 11:58',format_data) # Capture the time of J2000 Epoch
    [repochE,vepochE] = OrbitalElementsToRV(aE,eE,iE,raanE,omegaE,thetaE,mu,raddeg,table,distance,time)
    [repochM,vepochM] = OrbitalElementsToRV(aM,eM,iM,raanM,omegaM,thetaM,mu,raddeg,table,distance,time)
    # -- MAKE INITIAL ESTIMATE FOR DATE OF EARTH DEPARTURE --
    earth_departure_date = datetime.strptime(departure_data, format_data) # Establish initial guess of Earth departure at 1/1/21 00:00
    failure_date = date = datetime.strptime("1/1/24 12:00", format_data)
    mars_arrival_date = datetime.strptime(arrival_data, format_data) # Establish initial guess of Mars arrival after 190 day TOF
    i = 0 # Time start at midnight 00:00 (minutes passes counter)
    # -- ENTER INTERATIVE LOOP --
    
    
    
    
    
    
    
    
    while True:
        # -- CALCULATE INITIAL TIMES OF FLIGHT -- 
        earthTIMECHANGE = (earth_departure_date - J2000Epoch) # Amount of time passed since J2000 EPOCH to Earth Departure
        TOFearth = earthTIMECHANGE.total_seconds() / 58.13 # [TUsun] Time of Flight of Earth since Epoch
        marsTIMECHANGE = (mars_arrival_date - J2000Epoch) # Amount of time passed since J2000 EPOCH to Mars Arrival
        TOFmars = marsTIMECHANGE.total_seconds() / 58.13 # [TUsun] Time of Flight of Mars since Epoch
        # -- CALCULATE INITIAL ESTIMATES FOR DEPARTURE AND ARRIVAL VECTORS --
        [rEDepart,vEDepart] = universalTOF_SCZ(repochE,vepochE,TOFearth,mu) 
        [rMArrival,vMArrival] = universalTOF_SCZ(repochM,vepochM,TOFmars,mu)
        rE = math.sqrt(rEDepart[0]**2 + rEDepart[1]**2 + rEDepart[2]**2) # [AU] Magnitude of rEDepart
        rM = math.sqrt(rMArrival[0]**2 + rMArrival[1]**2 + rMArrival[2]**2) # [AU] Magnitude of rMArrival
        vE = math.sqrt(vEDepart[0]**2 + vEDepart[1]**2 + vEDepart[2]**2) # [AU/TU] Magnitude of vEDepart
        vM = math.sqrt(vMArrival[0]**2 + vMArrival[1]**2 + vMArrival[2]**2) # [AU/TU] Magnitude of vMArrival
        # -- IMPLEMENT GAUSS'S PROBLEM --
        [vE_guess,vM_guess] = gauss(rEDepart,rMArrival,TUtimeOfFlight,mu,1,"S")
        vEguess_mag = math.sqrt(vE_guess[0]**2 + vE_guess[1]**2 + vE_guess[2]**2)
        vMguess_mag = math.sqrt(vM_guess[0]**2 + vM_guess[1]**2 + vM_guess[2]**2)
        vEartherror = abs(vEguess_mag - vE)
        vMarserror = abs(vMguess_mag - vM)
        if (vEartherror <= tolerance and vMarserror <= tolerance) or earth_departure_date == failure_date:
            if earth_departure_date == failure_date:
                print("PROGRAM FAILED. SYNODIC PERIOD PASSED.")
            else:
                print("A DATE WAS FOUND")
                print(earth_departure_date)
                print(mars_arrival_date)
                print(vEartherror)
                print(vMarserror)
            break
        else:
            earth_departure_date = earth_departure_date + timedelta(days=1)
            mars_arrival_date = earth_departure_date + timedelta(days=190)
            i = i+1
            print(earth_departure_date)
    return
