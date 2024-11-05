# MAE 469 - Introduction to Astrodynamics
# Makenzie Nicole Karakula
# Project Attempt
# =======================================
# ---- IMPORT LIBRARIES ----
import math
import numpy as np

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