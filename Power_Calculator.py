## UR5e


# ***** lib
import numpy as np
from numpy import linalg


import cmath
import math
from math import cos
from math import sin
from math import atan2
from math import acos
from math import asin
from math import sqrt
from math import pi
import pandas as pd
import datetime
from datetime import timedelta
from scipy.signal import resample
from numpy import diff


global mat
mat=np.matrix


# ****** Coefficients ******


global d1, a2, a3, a7, d4, d5, d6
d1 =  0.1625
a2 = -0.425
a3 = -0.3922
a7 = 0.075
d4 =  0.1333
d5 =  0.0997
d6 =  0.0996 + 0.05 # additional 

T = []
c1 = mat([[0],[ -0.02561],[ 0.00193], [1]])
c2 = mat([[0.2125],[0],[0.1136], [1]])
c3 = mat([[0.15], [0.0], [0.0265], [1]])
c4 = mat([[0],[-0.0018],[0.01634], [1]])
c5 = mat([[0],[0.0018],[0.01634], [1]])
c6 = mat([[0],[0],[-0.001159 + 0.025], [1]])

m = [3.761, 8.058, 2.846, 1.37, 1.3, (0.365 + 4.484)] #Kg ur5e

global d, a, alph


d = mat([0.1625, 0, 0, 0.1333, 0.0997, 0.0996])#ur5e mm

a =mat([0 ,-0.425 ,-0.3922 ,0 ,0 ,0])#ur5e mm

alph = mat([pi/2, 0, 0, pi/2, -pi/2, 0 ]) # u5e


# ************************************************** FORWARD KINEMATICS

def AH( n,th,c  ):

  T_a = mat(np.identity(4), copy=False)
  T_a[0,3] = a[0,n-1]
  T_d = mat(np.identity(4), copy=False)
  T_d[2,3] = d[0,n-1]

  Rzt = mat([[cos(th[n-1,c]), -sin(th[n-1,c]), 0 ,0],
	         [sin(th[n-1,c]),  cos(th[n-1,c]), 0, 0],
	         [0,               0,              1, 0],
	         [0,               0,              0, 1]],copy=False)
      

  Rxa = mat([[1, 0,                 0,                  0],
			 [0, cos(alph[0,n-1]), -sin(alph[0,n-1]),   0],
			 [0, sin(alph[0,n-1]),  cos(alph[0,n-1]),   0],
			 [0, 0,                 0,                  1]],copy=False)

  A_i = T_d * Rzt * T_a * Rxa
	    

  return A_i

def HTrans(th,c ):  
  A_1=AH( 1,th,c  )
  A_2=AH( 2,th,c  )
  A_3=AH( 3,th,c  )
  A_4=AH( 4,th,c  )
  A_5=AH( 5,th,c  )
  A_6=AH( 6,th,c  )
  
  T_01=A_1
  T_02=A_1*A_2
  T_03=A_1*A_2*A_3
  T_04=A_1*A_2*A_3*A_4
  T_05=A_1*A_2*A_3*A_4*A_5
  T_06=A_1*A_2*A_3*A_4*A_5*A_6


  T = [T_01,T_02,T_03,T_04,T_05,T_06]
  return T


# ************************************************ Center of Mass

def CM(th, c):
    
    R = HTrans(th,c )
    c_theta1 = R[0]*c1*m[0]
    c_theta2 = R[1]*c2*m[1]
    c_theta3 = R[2]*c3*m[2]
    c_theta4 = R[3]*c4*m[3]
    c_theta5 = R[4]*c5*m[4]
    c_theta6 = R[5]*c6*m[5]
    
    M = m[0] + m[1]+ m[2]+ m[3]+ m[4]+ m[5]
    C_M = (c_theta1 + c_theta2 + c_theta3 + c_theta4 + c_theta5 + c_theta6)/M
    
    return C_M

########################################################################

################################################################
force = pd.read_csv('test_2_ompl-forcePlate.csv', usecols= ['fx','fy','fz','mz'])
#print(force)


velocity = pd.read_csv('test_2_ompl-joint_states.csv', usecols= ['V1','V2','V3','V4','V5','V6'])
#print(velocity)

current = pd.read_csv('test_2_ompl-joint_states.csv', usecols= ['I1','I2','I3','I4','I5','I6'])
#print(current)

position = pd.read_csv('test_2_ompl-joint_states.csv', usecols= ['P1','P2','P3','P4','P5','P6'])
#print(position)

time = pd.read_csv('test_2_ompl-joint_states.csv', usecols= ['time'])
#print(time)

#######################################################

#data_time
datetimeFormat = '%Y/%m/%d/%H:%M:%S.%f'
date2 = str(time.time.iloc[0])
date1 = str(time.time.iloc[-1])
diff = datetime.datetime.strptime(date1, datetimeFormat)\
    - datetime.datetime.strptime(date2, datetimeFormat)

newNumSamples = int((diff.seconds + diff.microseconds/1000000) * 120)
##################################################################

#resampling
newForce = pd.DataFrame(resample(force, newNumSamples), columns = ['fx','fy','fz','mz']).interpolate()
newVelocity = pd.DataFrame(resample(velocity, newNumSamples), columns = ['V3','V2','V1','V4','V5','V6']).interpolate()
newCurrent = pd.DataFrame(resample(current, newNumSamples), columns = ['I3','I2','I1','I4','I5','I6']).interpolate()
newPosition = pd.DataFrame(resample(position, newNumSamples), columns = ['P3','P2','P1','P4','P5','P6']).interpolate()


#print(newForce)
#print(newVelocity['V1'])
#print(newCurrent)
#print(newPosition['P1'])

#########################################

"""Joint Angle (in degrees) reading from encoders.
theta1 = np.radians(302)
theta2 = np.radians(190)
theta3 = np.radians(59)
theta4 = np.radians(199)
theta5 = np.radians(120)
theta6 = np.radians(90)"""
Cm = []
cx = []
cy = []
cz = []
row,column  = newPosition.shape
for x in range(row):
    
    theta1 = newPosition['P1'][x]
    theta2 = newPosition['P2'][x]
    theta3 = newPosition['P3'][x]
    theta4 = newPosition['P4'][x]
    theta5 = newPosition['P5'][x]
    theta6 = newPosition['P6'][x]

    th = np.matrix([[theta1], [theta2], [theta3], [theta4], [theta5], [theta6]])
    c = [0]
    #location = HTrans(th,c )
    #print(location)
    Z = CM(th,c)
    X = Z.tolist()
    Cm.append(X)
#print(len(Cm))


###########################################################################
#print(Cm)

dt = 0.00833

for x in range(len(Cm)):
    cx.append(Cm[x][0][0])
    cy.append(Cm[x][1][0])
    cz.append(Cm[x][2][0])

cx_ = pd.DataFrame (cx, columns = ['cx'])
cy_ = pd.DataFrame (cy, columns = ['cy'])
cz_ = pd.DataFrame (cz, columns = ['cz'])

#print(cz_)

Cx_velocity = cx_.diff()/dt
Cy_velocity = cy_.diff()/dt
Cz_velocity = cz_.diff()/dt


newFZ = pd.DataFrame (newForce, columns = ['fz']).abs()
newFY = pd.DataFrame (newForce, columns = ['fy']).abs()
newFX = pd.DataFrame (newForce, columns = ['fx']).abs()
newMZ = pd.DataFrame (newForce, columns = ['mz']).abs()
newCz_velocity = Cz_velocity.fillna(0).abs()
newCy_velocity = Cy_velocity.fillna(0).abs()
newCx_velocity = Cx_velocity.fillna(0).abs()

newVz_base = pd.DataFrame (newVelocity, columns = ['V1']).abs()


P2z = pd.DataFrame(newFZ.values*newCz_velocity.values, columns = ['P2z'])
P2y = pd.DataFrame(newFY.values*newCy_velocity.values, columns = ['P2y'])
P2x = pd.DataFrame(newFX.values*newCx_velocity.values, columns = ['P2x'])

P2 = pd.DataFrame(P2z.values + P2y.values + P2x.values, columns = ['P2'])

P1 = pd.DataFrame(newMZ.values*newVz_base.values, columns = ['P1'])

#print(newCz_velocity)
#print(newFZ)
#print(newMZ)
print(P2z.describe())
print(P1.describe())


P = pd.DataFrame((newCurrent['I1'].values + newCurrent['I2'].values + newCurrent['I3'].values + newCurrent['I4'].values +newCurrent['I5'].values + newCurrent['I6'].values)*47.5, columns = ['P']).abs()
print(P.describe())
