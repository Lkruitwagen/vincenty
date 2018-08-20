from __future__ import print_function 
import math


# WGS 84
a = 6378137  # meters
f = 1 / 298.257223563
b = 6356752.314245  # meters; b = (1 - f)a

MILES_PER_KILOMETER = 0.62137119

MAX_ITERATIONS = 200
CONVERGENCE_THRESHOLD = 1e-12  # .000,000,000,001


def V_inv(point1, point2, miles=False):


    """
    Vincenty's formula (inverse method) to calculate the distance (in
    kilometers or miles) between two points on the surface of a spheroid

    ### INPUT ###
    point1: tuple/list, (lat,lon) coordinates of departure point
    point2: tuple/list, (lat,lon) coordinates of destination point
    miles: bool, distance output give in miles or km

    ### OUTPUT ###
    dist: float, distance from point 1 to point2
    alpha1: float, departure azimuth angle (clockwise from N=0.0)
    alpha2: float, arrival azimuth angle (clockwise from N=0.0)

    Doctests:
    >>> V_inv((0.0, 0.0), (0.0, 0.0))  # coincident points
    (0.0, 0.0, 0.0)
    >>> V_inv((0.0, 0.0), (0.0, 1.0))
    (111.319491, 90.0, 90.0)
    >>> V_inv((0.0, 0.0), (1.0, 0.0))
    (110.574389, 0.0, 0.0)
    >>> V_inv((0.0, 0.0), (0.5, 179.5))  # slow convergence
    (19936.288579, 25.67187285623632, 154.32708548199773)
    >>> V_inv((0.0, 0.0), (0.5, 179.7))  # failure to converge
    >>> boston = (42.3541165, -71.0693514)
    >>> newyork = (40.7791472, -73.9680804)
    >>> V_inv(boston, newyork)
    (298.396057, 235.0838926194711, 233.1602005520544)
    >>> V_inv(boston, newyork, miles=True)
    (185.414713, 235.0838926194711, 233.1602005520544)
    >>> V_inv(newyork, boston)
    (298.396057, 53.1602005520544, 55.0838926194711)
    """

    # short-circuit coincident points
    if point1[0] == point2[0] and point1[1] == point2[1]:
        return 0.0,0.0,0.0

    U1 = math.atan((1 - f) * math.tan(math.radians(point1[0])))
    U2 = math.atan((1 - f) * math.tan(math.radians(point2[0])))
    L = math.radians(point2[1] - point1[1])
    Lambda = L

    sinU1 = math.sin(U1)
    cosU1 = math.cos(U1)
    sinU2 = math.sin(U2)
    cosU2 = math.cos(U2)

    for iteration in range(MAX_ITERATIONS):
        sinLambda = math.sin(Lambda)
        cosLambda = math.cos(Lambda)
        sinSigma = math.sqrt((cosU2 * sinLambda) ** 2 +
                             (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda) ** 2)
        if sinSigma == 0:
            return 0.0  # coincident points
        cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda
        sigma = math.atan2(sinSigma, cosSigma)
        sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma
        cosSqAlpha = 1 - sinAlpha ** 2
        try:
            cos2SigmaM = cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha
        except ZeroDivisionError:
            cos2SigmaM = 0
        C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha))
        LambdaPrev = Lambda
        Lambda = L + (1 - C) * f * sinAlpha * (sigma + C * sinSigma *
                                               (cos2SigmaM + C * cosSigma *
                                                (-1 + 2 * cos2SigmaM ** 2)))
        if abs(Lambda - LambdaPrev) < CONVERGENCE_THRESHOLD:
            break  # successful convergence
    else:
        return None  # failure to converge

    uSq = cosSqAlpha * (a ** 2 - b ** 2) / (b ** 2)
    A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)))
    B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)))
    deltaSigma = B * sinSigma * (cos2SigmaM + B / 4 * (cosSigma *
                 (-1 + 2 * cos2SigmaM ** 2) - B / 6 * cos2SigmaM *
                 (-3 + 4 * sinSigma ** 2) * (-3 + 4 * cos2SigmaM ** 2)))
    s = b * A * (sigma - deltaSigma)


    num = (math.cos(U2)*math.sin(Lambda))
    den = (math.cos(U1)*math.sin(U2)-math.sin(U1)*math.cos(U2)*math.cos(Lambda))

    #print 'num',num
    #print 'den',den
    alpha1 = math.atan2(num,den)

    if alpha1<0:
        alpha1+=2*math.pi



    num = (math.cos(U1)*math.sin(Lambda))
    den = (-1.0*math.sin(U1)*math.cos(U2)+math.cos(U1)*math.sin(U2)*math.cos(Lambda))
    #print 'num',num
    #print 'den',den
    alpha2 = math.atan2(num,den)

    if alpha2<0:
        alpha2+=2*math.pi


    s /= 1000  # meters to kilometers
    if miles:
        s *= MILES_PER_KILOMETER  # kilometers to miles

    return round(s, 6), math.degrees(alpha1), math.degrees(alpha2)


def V_dir(point1, s, alpha1,miles=False):

    """
    Vincenty's formula (direct method) to calculate the coordinates of a point on a spheroid 
    given an initial point, a geodesic distance, and a departure azimuth.

    ### INPUT ###
    point1: tuple/list, (lat,lon) coordinates of departure point
    s: float, geodesic distance in miles/km
    alpha1: float departure azimuth angle in degrees (clockwise from N=0.0)
    miles: bool, distance input given in miles or km

    ### OUTPUT ###
    point: tuple, (lat,lon) coordinates of arrival point
    alpha2: float, arrival azimuth angle (clockwise from N=0.0)

    Doctests:
    >>> V_dir((0.0, 0.0), 0.0, 0.0)  #no direction
    ((0.0, 0.0), 0.0)
    >>> boston = (42.3541165, -71.0693514)
    >>> newyork = (40.7791472, -73.9680804)
    >>> V_dir(boston, V_inv(boston,newyork)[0], V_inv(boston,newyork)[1])
    ((40.779147202568396, -73.96808039548996), 233.16020055500005)
    >>> V_dir(boston, V_inv(boston,newyork,miles=True)[0], V_inv(boston,newyork)[1],miles=True)
    ((40.77914720282384, -73.96808039504143), 233.16020055529302)
    """

    if miles:
        s *=1.0/MILES_PER_KILOMETER  #put s into km

    s*=1000 #s from km to m

    #alpha1 in degrees
    alpha1=math.radians(alpha1)
    U1 = math.atan((1.0-f)*math.tan(math.radians(point1[0])))
    #print U1
    sigma1 = math.atan2((math.tan(U1)),(math.cos(alpha1)))
    sinAlpha=math.cos(U1)*math.sin(alpha1)
    cosSqAlpha=1.0-(sinAlpha**2)
    uSq = cosSqAlpha*(a**2-b**2)/(b**2)
    A = 1 + uSq/16384.0*(4096.0+uSq*(-768.0+uSq*(320.0-175*uSq)))
    B = uSq/1024*(256+uSq*(-128+uSq*(74-47*uSq)))

    sigma=s/b/A
    #print sigma
    for iteration in range(MAX_ITERATIONS):

        sigma2m = 2*sigma1+sigma
        deltasigma = B*math.sin(sigma)*(math.cos(sigma2m)+1.0/4*B*(math.cos(sigma)*(-1+2*(math.cos(sigma2m)**2))-1.0/6*B*math.cos(sigma2m)*(-3+4*(math.sin(sigma)**2))*(-3+4*(math.cos(sigma2m)**2))))
        sigmaprev = sigma
        sigma = s/b/A+deltasigma
        #print sigma
        if abs(sigma - sigmaprev) < CONVERGENCE_THRESHOLD:
            #print 'converge'
            break  # successful convergence
    else:
        print ('no convergence')
        return None  # failure to converge


    num = math.sin(U1)*math.cos(sigma)+math.cos(U1)*math.sin(sigma)*math.cos(alpha1)
    den = (1.0-f)*math.sqrt(sinAlpha**2+(math.sin(U1)*math.sin(sigma)-math.cos(U1)*math.cos(sigma)*math.cos(alpha1))**2)
    #print num
    #print den
    lat2= math.atan2(num,den)

    num=math.sin(sigma)*math.sin(alpha1)
    den = math.cos(U1)*math.cos(sigma)-math.sin(U1)*math.sin(sigma)*math.cos(alpha1)
    Lambda = math.atan2(num,den)

    C = f/16.0*(cosSqAlpha*(4+f*(4.0-3.0*cosSqAlpha)))
    L = Lambda - (1.0-C)*f*sinAlpha*(sigma+C*math.sin(sigma)*(math.cos(sigma2m)+C*math.cos(sigma)*(-1+2.0*(math.cos(sigma2m)**2))))

    L2 = math.radians(point1[1])+L
    num = sinAlpha
    den = -1*math.sin(U1)*math.sin(sigma)+math.cos(U1)*math.cos(sigma)*math.cos(alpha1)
    #print num
    #print den
    alpha2 = math.atan2(num,den)
    if alpha2<0:
        alpha2+=math.pi*2
    #print alpha2
    # short-circuit coincident points
    return (math.degrees(lat2),math.degrees(L2)),math.degrees(alpha2)



if __name__ == '__main__':
    import doctest
    doctest.testmod()
