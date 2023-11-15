#!/usr/bin/env python3
# PROGRAM-NAME: caseCloudLinks.py

import os,sys,shutil, fileinput
import math
import matplotlib.pyplot as plt

import numpy as np




currentDir = os.getcwd()+'/'

#Variables
points = np.zeros((110,3))
nb_points = 105
c=0.6 #chord
L = 0.046*c # Distance between points 0 and 82 to define the layer

#==================== START OF FIRST LAYER PARAMETER REGION ============
#Kinematic viscosity
nu = 1.5e-05

#Reynolds
Re = 1600000

#y+ in first layer
yPlus = 10

# Calculation of first layer height
# Reference : https://www.computationalfluiddynamics.com.au/tips-tricks-cfd-estimate-first-cell-height/

#1) Compute velocity
U = nu*Re/c
print 'U =', U

#2) Compute Cf
#For external flows Cf = 0.058 Re^{-0.20}
C_f = 0.058 *Re**(-0.20)

#3) Compute tau_w / rho
Tau_w_rho = 0.5 * C_f * U**2

#4) Compute U_tau
U_tau = Tau_w_rho**(0.5)

#5) Compute first layer height
y1 = yPlus * nu / U_tau
# print 'y1 =', y1

dl0 = y1 #Heigth of the smallest cell around the profile
# print 'dl0 =', dl0
# spanwise ratio
# dzx = 200*dl0
dzx = 0.5*c
print 'dzx = ', dzx
# total spanwise length
deltaZ = -0.5*c
print 'deltaZ = ', deltaZ
print 'int(deltaZ/dzx) = ', int(deltaZ/dzx)

# sys.exit()

cle = 170 # number of cells on the leading edge of the profile (between points 0 and 1)
cte = 200 # number of cells on the trailing edge of the profile (between points 1 and 3)
smc = 4 # number of cells between the smaller
cd = 40 #number of cells downstream in the static mesh (between points 49 and 38)
cv = 40 #number of cells in the vertical direction in the static part
#compute from the variables
clayer = int(round(L/dl0)) #Number of cells in the layer around the profil
# scyl = c + c/3
scyl = 2*c
mcyl = scyl+c/10
lcyl = mcyl+c/10
#Get the profil coordinates
# spline = currentDir+'spline_NACA0015.txt'
spline = currentDir+'spline_FFA-W3-301.txt'
xu = []
zu = []
xl = []
zl = []
for line in fileinput.input(spline):
	words = line.split()
	if words[0]=='upper':
		xu.append(c*float(words[1]))
		zu.append(c*float(words[3]))
	if words[0]=='lower':
		xl.append(c*float(words[1]))
		zl.append(c*float(words[3]))

#Enlargement of the upper profile
# alphau = []
betau = []
gammau = []
counter = 0
x2u = []
z2u = []
# z3u = []
xu_av = 0.0
zu_av = 0.0
x0 = c/2
xL = c
for i in xu:
	# alphau.append(np.arctan(zu[counter]/(c/2-xu[counter])))
	if zu[counter] < 0:
		gammau.append(np.arctan((zu_av-zu[counter])/(xu[counter]-xu_av)))
	else:
		gammau.append(np.arctan((zu[counter]-zu_av)/(xu[counter]-xu_av)))
	betau.append(1.570796327-gammau[counter])
	x2u.append(xu[counter]-L*np.cos(betau[counter]))
	z2u.append(zu[counter]+L*np.sin(betau[counter]))
	# z3u.append(-z2u[counter])
	xu_av = xu[counter]
	zu_av = zu[counter]

	counter += 1

#Enlargement of the lower profile (now working)
# alphal = []
betal = []
gammal = []
counter = 0
x2l = []
z2l = []
# z3l = []
xl_av = 0.0
zl_av = 0.0
x0 = c/2
xL = c
for i in xl:
	# alphal.append(np.arctan(zl[counter]/(c/2-xl[counter])))
	gammal.append(np.arctan((zl_av-zl[counter])/(xl[counter]-xl_av)))
	betal.append(1.570796327-gammal[counter])
	x2l.append(xl[counter]-L*np.cos(betal[counter]))
	z2l.append(zl[counter]-L*np.sin(betal[counter]))
	# z3l.append(-z2l[counter])
	xl_av = xl[counter]
	zl_av = zl[counter]

	counter += 1

# To picture the blade profile, upper, and lower layer
# plt.plot(xu,zu)
# plt.plot(xl,zl)
# plt.plot(x2u,z2u)
# plt.plot(x2l,z2l)
autoRunDir = os.getcwd()+'/'

# Get the outputDir
if os.path.exists(autoRunDir + 'profileImage'):
	pass
else:
	os.makedirs(autoRunDir + 'profileImage')

outputDir = autoRunDir + 'profileImage/'
outputFile = outputDir+'profileLayers_nonRotating.png'
plt.savefig(outputFile)
plt.close()
# sys.exit()
# end picturing

def Q_N(L,dl0,dlN):
	Q=0.0
	Q=dlN/dl0
	q=0.0
	q=(dl0/L-1)/(Q*dl0/L-1)
	N=0
	N=int(round(np.log(Q)/np.log(q)+1))
	return [L,N,Q,q,dl0,dlN]

def Q(L,dl0,N):
	if N*dl0>L:
		q=0.0+1E-10
	else:
		q=1.0+1E-10
	equation=0.0
	equation = dl0*(1-q**N)/(1-q)-L
	while equation<0:
		q += 1E-7
		equation = dl0*(1-q**N)/(1-q)-L
	dlN=0.0
	dlN=dl0*q**(N-1)
	Q=0.0
	Q=dlN/dl0
	return [L,N,Q,q,dl0,dlN]

def Q_N_bis(L,dl0,q):
	N=int(round(np.log(L/dl0*(q-1)+1)/np.log(q)))
	dlN=dl0*q**(N-1)
	Q=dlN/dl0
	return [L,N,Q,q,dl0,dlN]

#shape
counter = 0
z0 = zu[counter]
points[0,:]=[0.0,0,z0] #profile noze

counter = 0
z1 = 0.0
for i in xu:
	if i < 2*c/5:
		z1= (zu[counter]+zu[counter+1])/2
		x1= (xu[counter]+xu[counter+1])/2
		# z1= zu[counter]
	counter += 1
points[1,:]=[x1,0,z1]

counter = 0
z2 = 0.0
for i in xl:
	if i < 2*c/5:
		z2= (zl[counter]+zl[counter+1])/2
		x2= (xl[counter]+xl[counter+1])/2
		# z2= zl[counter]
	counter += 1
points[2,:]=[x2,0,z2]

counter = 0
z3 = 0.0
for i in xl:
	if i <= c:
		# z61= (z2u[counter]+z2u[counter+1])/2
		z3= zl[counter]
	counter += 1
points[3,:]=[c,0,z3]

points[4,:]=[0.0,deltaZ,z0]
points[5,:]=[x1,deltaZ,z1]
points[6,:]=[x2,deltaZ,z2]
points[7,:]=[c,deltaZ,z3]
#end shape

#smaller cylinder
points[8,:]=[x1, 0 ,scyl]
points[9,:]=[(2*c/5) + scyl, 0 ,z3]
points[10,:]=[x2, 0 ,-scyl]
points[11,:]=[2*c/5 - scyl, 0 ,z0]
points[12,:]=[x1,deltaZ,scyl]
points[13,:]=[2*c/5 + scyl,deltaZ,z3]
points[14,:]=[x2,deltaZ,-scyl]
points[15,:]=[2*c/5 - scyl,deltaZ,z0]

#Inlet
points[16,:]=[-6*c, 0 ,z0]
points[17,:]=[0.0, 0 ,6*c]
points[18,:]=[0.0, 0 ,-6*c]
points[19,:]=[-6*c,deltaZ,z0]
points[20,:]=[0.0,deltaZ,6*c]
points[21,:]=[0.0,deltaZ,-6*c]

#Points downstream
points[22,:]=[2*c/5+lcyl*np.cos(0.7854),0,lcyl*np.sin(0.7854)]
points[23,:]=[2*c/5+lcyl*np.cos(0.7854),0,-lcyl*np.sin(0.7854)]
points[24,:]=[2*c/5+lcyl*np.cos(0.7854),0,6*c]
points[25,:]=[2*c/5+lcyl*np.cos(0.7854),0,-6*c]
points[26,:]=[2*c/5+lcyl*np.cos(0.7854),deltaZ,lcyl*np.sin(0.7854)]
points[27,:]=[2*c/5+lcyl*np.cos(0.7854),deltaZ,-lcyl*np.sin(0.7854)]
points[28,:]=[2*c/5+lcyl*np.cos(0.7854),deltaZ,6*c]
points[29,:]=[2*c/5+lcyl*np.cos(0.7854),deltaZ,-6*c]
points[30,:]=[20*c,0,lcyl*np.sin(0.7854)]
points[31,:]=[20*c,0,-lcyl*np.sin(0.7854)]
points[32,:]=[20*c,deltaZ,lcyl*np.sin(0.7854)]
points[33,:]=[20*c,deltaZ,-lcyl*np.sin(0.7854)]
points[34,:]=[20*c,0,6*c]
points[35,:]=[20*c,0,-6*c]
points[36,:]=[20*c,deltaZ,6*c]
points[37,:]=[20*c,deltaZ,-6*c]
points[38,:]=[20*c,0,z3]
points[39,:]=[20*c,deltaZ,z3]

#medium cylinder
points[40,:]=[x1, 0 ,mcyl]
points[41,:]=[2*c/5 + mcyl, 0 ,z3]
points[42,:]=[x2, 0 ,-mcyl]
points[43,:]=[2*c/5 - mcyl, 0 ,z0]
points[44,:]=[x1,deltaZ,mcyl]
points[45,:]=[2*c/5 + mcyl,deltaZ,z3]
points[46,:]=[x2,deltaZ,-mcyl]
points[47,:]=[2*c/5 - mcyl,deltaZ,z0]

#larger cylinder
points[48,:]=[x1,0,lcyl]
points[49,:]=[2*c/5+lcyl,0,z3]
points[50,:]=[x2,0,-lcyl]
points[51,:]=[2*c/5-lcyl,0,z0]
points[52,:]=[x1,deltaZ,lcyl]
points[53,:]=[2*c/5+lcyl,deltaZ,z3]
points[54,:]=[x2,deltaZ,-lcyl]
points[55,:]=[2*c/5-lcyl,deltaZ,z0]

#Points at pi/4 on medium cylider
points[56,:]=[2*c/5+mcyl*np.cos(0.7854),0,mcyl*np.sin(0.7854)]
points[57,:]=[2*c/5+mcyl*np.cos(0.7854),0,-mcyl*np.sin(0.7854)]
points[58,:]=[2*c/5+mcyl*np.cos(0.7854),deltaZ,mcyl*np.sin(0.7854)]
points[59,:]=[2*c/5+mcyl*np.cos(0.7854),deltaZ,-mcyl*np.sin(0.7854)]

#Points in the rotating part
alpha60 = np.arccos((c - 2*c/5)/scyl)
points[60,:]=[c,0,scyl*np.sin(alpha60)]

counter = 0
z61 = 0.0
for i in x2u:
	if i <= c:
		# z61= (z2u[counter]+z2u[counter+1])/2
		z61= z2u[counter]
	counter += 1
points[61,:]=[c,0,z61]

counter = 0
z62 = 0.0
for i in x2l:
	if i <= c:
		# z62= (z2l[counter]+z2l[counter-1])/2
		z62= z2l[counter]
	counter += 1
points[62,:]=[c,0,z62]

points[63,:]=[c,0,-scyl*np.sin(alpha60)]
points[65,:]=[c+12*L, 0 ,z61*1]

#placement of the point 64 to minimize the non-orthogonality
#angle between (65,66) and (65,64) should be 135 degrees
#point 64 is at the intersection of the circle (z=scyl*sin(alpha64) and x=c/2+scyl*cos(alpha64) and the straight line with an angle of 45 degrees (z=ax+b and a=1)
#This lead to an quadratic equation
b=points[65,2]-points[65,0]
delta =(2*(2*c/5+b)/scyl)**2-4*2*(((2*c/5+b)/scyl)**2-1)
alpha64 = np.arccos((-2*(2*c/5+b)/scyl+np.sqrt(delta))/(2*2)) / 1.5

points[64,:]=[2*c/5+scyl*np.cos(alpha64),0,scyl*np.sin(alpha64)]
points[66,:]=[c+12*L, 0 ,z3]
points[67,:]=[c+12*L, 0 ,z62*1]
points[68,:]=[2*c/5+scyl*np.cos(alpha64),0,-scyl*np.sin(alpha64)]
points[69,:]=[c,deltaZ,scyl*np.sin(alpha60)]
points[70,:]=[c,deltaZ,z61]
points[71,:]=[c,deltaZ,z62]
points[72,:]=[c,deltaZ,-scyl*np.sin(alpha60)]
points[73,:]=[2*c/5+scyl*np.cos(alpha64),deltaZ,scyl*np.sin(alpha64)]
points[74,:]=[c+12*L,deltaZ,z61*1]
points[75,:]=[c+12*L,deltaZ,z3]
points[76,:]=[c+12*L,deltaZ,z62*1]
points[77,:]=[2*c/5+scyl*np.cos(alpha64),deltaZ,-scyl*np.sin(alpha64)]

counter = 0
z78 = 0.0
for i in x2u:
	if i <= 2*c/5:
		z78= (z2u[counter]+z2u[counter+1])/2
	counter += 1
points[78,:]=[x1,0,z78]

counter = 0
z79 = 0.0
for i in x2l:
	if i <= 2*c/5:
		z79= (z2l[counter]+z2l[counter+1])/2
	counter += 1
points[79,:]=[x2,0,z79]

points[80,:]=[x1,deltaZ,z78]
points[81,:]=[x2,deltaZ,z79]
points[82,:]=[points[0,0]-L,0,z0]
points[83,:]=[points[0,0]-L,deltaZ,z0]
points[84,:]=[2*c/5+mcyl*np.cos(alpha60),0,mcyl*np.sin(alpha60)]
points[85,:]=[2*c/5+mcyl*np.cos(alpha64),0,mcyl*np.sin(alpha64)]
points[86,:]=[2*c/5+mcyl*np.cos(alpha64),0,-mcyl*np.sin(alpha64)]
points[87,:]=[2*c/5+mcyl*np.cos(alpha60),0,-mcyl*np.sin(alpha60)]
points[88,:]=[2*c/5+mcyl*np.cos(alpha60),deltaZ,mcyl*np.sin(alpha60)]
points[89,:]=[2*c/5+mcyl*np.cos(alpha64),deltaZ,mcyl*np.sin(alpha64)]
points[90,:]=[2*c/5+mcyl*np.cos(alpha64),deltaZ,-mcyl*np.sin(alpha64)]
points[91,:]=[2*c/5+mcyl*np.cos(alpha60),deltaZ,-mcyl*np.sin(alpha60)]

# new points fo rthe flatBack profile
counter = 0
z92 = 0.0
for i in xu:
	if i <= c:
		# z62= (z2l[counter]+z2l[counter-1])/2
		z92= zu[counter]
	counter += 1
points[92,:]=[c,0,z92]
points[93,:]=[c+12*L, 0 ,z92*1.5]
points[94,:]=[2*c/5 + scyl*np.cos(alpha64/2) , 0 , scyl*np.sin(alpha64/2)]
points[95,:]=[2*c/5 + mcyl*np.cos(alpha64/2) , 0 , mcyl*np.sin(alpha64/2)]
points[96,:]=[c,deltaZ,z92]
points[97,:]=[c+12*L,deltaZ,z92*1.5]
points[98,:]=[2*c/5 + scyl*np.cos(alpha64/2) ,deltaZ, scyl*np.sin(alpha64/2)]
points[99,:]=[2*c/5 + mcyl*np.cos(alpha64/2) ,deltaZ, mcyl*np.sin(alpha64/2)]

points[100,:]=[2*c/5 + mcyl*np.cos(alpha64/2), 0 ,mcyl*np.sin(alpha64/2)]
points[101,:]=[2*c/5 + lcyl*np.cos(alpha64/2), 0 ,lcyl*np.sin(alpha64/2)]
points[102,:]=[20*c, 0 ,lcyl*np.sin(alpha64/2)]
points[103,:]=[2*c/5+mcyl*np.cos(alpha64/2),deltaZ,mcyl*np.sin(alpha64/2)]
points[104,:]=[2*c/5+lcyl*np.cos(alpha64/2),deltaZ,lcyl*np.sin(alpha64/2)]
points[105,:]=[20*c,deltaZ,lcyl*np.sin(alpha64/2)]


# ending new points

def Q_N(L,dl0,dlN):
	Q=0.0
	Q=dlN/dl0
	q=0.0
	q=(dl0/L-1)/(Q*dl0/L-1)
	N=0
	N=int(round(np.log(Q)/np.log(q)+1))
	return [L,N,Q,q,dl0,dlN]

def Q(L,dl0,N):
	if N*dl0>L:
		q=0.0+1E-10
	else:
		q=1.0+1E-10
	equation=0.0
	equation = dl0*(1-q**N)/(1-q)-L
	while equation<0:
		q += 1E-7
		equation = dl0*(1-q**N)/(1-q)-L
	dlN=0.0
	dlN=dl0*q**(N-1)
	Q=0.0
	Q=dlN/dl0
	return [L,N,Q,q,dl0,dlN]

def Q_N_bis(L,dl0,q):
	N=int(round(np.log(L/dl0*(q-1)+1)/np.log(q)))
	dlN=dl0*q**(N-1)
	Q=dlN/dl0
	return [L,N,Q,q,dl0,dlN]


lastCellSizeLayer = 2*3.14159*(c/2)/(4*cle) #cell size of last cell within the layer
#from here until line 452, dl10 was changed by lastCellSizeLayer = 2*3.14159*(c/2)/(4*cle)

vd=Q_N(points[8,2]-points[78,2],lastCellSizeLayer,2*3.14159*scyl/(4*cle))
cvd=vd[1]#Number of cells in the vertical direction (between points 78 and 8)
Qvd=vd[2] #Grading in the vertical direction (between points 78 and 8)
dlNvd=vd[5]

ds=Q_N(points[66,0]-points[3,0],np.sqrt((points[92,0]-points[1,0])**2+(points[1,2]-points[92,2])**2)/cte,lastCellSizeLayer)
cds=ds[1]#number of cells downstream the shape (between points 3 and 66)
Qds=ds[2]  #Grading downstream the shape (between points 3 and 66)

f11_82=Q(points[82,0]-points[11,0],lastCellSizeLayer,cvd)
Q11_82=f11_82[2]
dlN11_82=f11_82[5]

f61_60=Q(points[60,2]-points[61,2],lastCellSizeLayer,cvd)
Q61_60=f61_60[2]
dlN61_60=f61_60[5]

f65_64=Q(np.sqrt((points[64,2]-points[65,2])**2+(points[64,0]-points[65,0])**2),lastCellSizeLayer,cvd)
Q65_64=f65_64[2]
dlN65_64=f65_64[5]

f93_94=Q(np.sqrt((points[94,2]-points[93,2])**2+(points[94,0]-points[93,0])**2),lastCellSizeLayer,cvd)
Q93_94=f93_94[2]
dlN93_94=f93_94[5]

f66_9=Q(points[9,0]-points[66,0],lastCellSizeLayer,cvd)
Q66_9=f66_9[2]
dlN66_9=f66_9[5]

dl0_8_60=(1.5708-alpha60)*scyl/cte
Q8_11=Q(1.5708*scyl,dl0_8_60,cle)[2]

f60_64=Q((alpha60-alpha64)*scyl,dl0_8_60,cds)
Q60_64=f60_64[2]
dlN60_64=f60_64[5]

Q64_9=Q(alpha64*scyl,dlN60_64,clayer)[2]
Q64_94=Q((alpha64/2)*scyl,dlN60_64,clayer)[2]
Q94_9=Q((alpha64/2)*scyl,dlN60_64,clayer)[2]

f11_43=Q(mcyl-scyl,dlN11_82,smc)
Q11_43=f11_43[2]
dlN11_43=f11_43[5]

f8_40=Q(mcyl-scyl,dlNvd,smc)
Q8_40=f8_40[2]
dlN8_40=f8_40[5]

Q60_84=Q(mcyl-scyl,dlN61_60,smc)[2]

Q64_85=Q(mcyl-scyl,dlN65_64,smc)[2]

dlN65_64=Q(mcyl-scyl,dlN65_64,smc)[5]

f9_41=Q(mcyl-scyl,dlN66_9,smc)
Q9_41=f9_41[2]
dlN9_41=f9_41[5]

Q94_95=Q(mcyl-scyl,dlN93_94,smc)[2]

f43_51=Q_N(lcyl-mcyl,2*3.14*mcyl/(4*cle),2*3.14*lcyl/(4*cle))
Q43_51=f43_51[2]
N43_51=f43_51[1]-1


Q49_38=Q(points[38,0]-points[49,0],2*3.14*lcyl/(4*cle),cd)[2]
Q22_30=Q(points[30,0]-points[22,0],2*3.14*lcyl/(4*cle),cd)[2]
Q101_102=Q(points[102,0]-points[101,0],2*3.14*lcyl/(4*cle),cd)[2]


Q48_17=Q(points[17,2]-points[48,2],2*3.14*lcyl/(4*cle),cv)[2]

Q22_24=Q(points[24,2]-points[22,2],2*3.14*lcyl/(4*cle),cv)[2]

Q17_16=Q(1.5708*13*c,(points[24,0]-points[17,0])/(cle/2),cle)[2]
#Write the blockMeshDict File

bmd = open("mesh/nonRotating/blockMeshDict","w")

bmd.write("/*-------------------------------*- C++ -*----------------------------------*\\"+'\n')
bmd.write('| =========                 |                                                |'+'\n')
bmd.write('| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |'+'\n')
bmd.write('|  \\    /   O peration     | Version:  2.2.2                                 |'+'\n')
bmd.write('|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |'+'\n')
bmd.write('|    \\/     M anipulation  |                                                 |'+'\n')
bmd.write("\\*--------------------------------------------------------------------------*/"+'\n')

bmd.write('FoamFile'+'\n')
bmd.write('{'+'\n')
bmd.write('    version     2.0;'+'\n')
bmd.write('    format      ascii;'+'\n')
bmd.write('    class       dictionary;'+'\n')
bmd.write('   object      blockMeshDict;'+'\n')
bmd.write('}'+'\n')
bmd.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //'+'\n')
bmd.write(''+'\n')
bmd.write('convertToMeters 1;'+'\n')
bmd.write(''+'\n')
bmd.write('vertices'+'\n')
bmd.write('('+'\n')
bmd.write(''+'\n')

counter1 = 0
for i in points[:,0]:
	if counter1 <= nb_points:
		bmd.write('	('+str(points[counter1,0])+' '+str(points[counter1,1])+' '+str(points[counter1,2])+') //'+str(counter1)+'\n')
	counter1 += 1

bmd.write(''+'\n')
bmd.write(');'+'\n')
bmd.write(''+'\n')
bmd.write('blocks'+'\n')
bmd.write('('+'\n')
bmd.write(''+'\n')
bmd.write('    hex (43 40 48 51 47 44 52 55) shape ('+str(cle/2)+' '+str(N43_51)+' '+str(int(-1*deltaZ/dzx))+ ') edgeGrading (1 1 1 1 '+str(Q43_51)+' '+str(Q43_51)+' '+str(Q43_51)+' '+str(Q43_51)+' 1 1 1 1)'+'\n')
bmd.write('    hex (40 56 22 48 44 58 26 52) shape ('+str(cle/4)+' '+str(N43_51)+' '+str(int(-1*deltaZ/dzx))+ ') edgeGrading (1 1 1 1 '+str(Q43_51)+' '+str(Q43_51)+' '+str(Q43_51)+' '+str(Q43_51)+' 1 1 1 1)'+'\n')
# bmd.write('    hex (56 41 49 22 58 45 53 26) shape ('+str(cle/2)+' '+str(N43_51)+' '+str(int(-1*deltaZ/dzx))+ ') edgeGrading (1 1 1 1 '+str(Q43_51)+' '+str(Q43_51)+' '+str(Q43_51)+' '+str(Q43_51)+' 1 1 1 1)'+'\n')
#2
bmd.write('    hex (56 100 101 22 58 103 104 26) shape ('+str(cle/2)+' '+str(N43_51)+' '+str(int(-1*deltaZ/dzx))+ ') edgeGrading (1 1 1 1 '+str(Q43_51)+' '+str(Q43_51)+' '+str(Q43_51)+' '+str(Q43_51)+' 1 1 1 1)'+'\n')
#3
bmd.write('    hex (100 41 49 101 103 45 53 104) shape ('+str(cle/8)+' '+str(N43_51)+' '+str(int(-1*deltaZ/dzx))+ ') edgeGrading (1 1 1 1 '+str(Q43_51)+' '+str(Q43_51)+' '+str(Q43_51)+' '+str(Q43_51)+' 1 1 1 1)'+'\n')
#4
bmd.write('    hex (41 57 23 49 45 59 27 53) shape ('+str(cle/2)+' '+str(N43_51)+' '+str(int(-1*deltaZ/dzx))+ ') edgeGrading (1 1 1 1 '+str(Q43_51)+' '+str(Q43_51)+' '+str(Q43_51)+' '+str(Q43_51)+' 1 1 1 1)'+'\n')
#5
bmd.write('    hex (57 42 50 23 59 46 54 27) shape ('+str(cle/4)+' '+str(N43_51)+' '+str(int(-1*deltaZ/dzx))+ ') edgeGrading (1 1 1 1 '+str(Q43_51)+' '+str(Q43_51)+' '+str(Q43_51)+' '+str(Q43_51)+' 1 1 1 1)'+'\n')
#6
bmd.write('    hex (42 43 51 50 46 47 55 54) shape ('+str(cle/2)+' '+str(N43_51)+' '+str(int(-1*deltaZ/dzx))+ ') edgeGrading (1 1 1 1 '+str(Q43_51)+' '+str(Q43_51)+' '+str(Q43_51)+' '+str(Q43_51)+' 1 1 1 1)'+'\n')
#7
bmd.write('    hex (51 48 17 16 55 52 20 19) shape ('+str(cle/2)+' '+str(cv)+' '+str(int(-1*deltaZ/dzx))+ ') edgeGrading (1 '+str(1/Q17_16)+' '+str(1/Q17_16)+' 1 '+str(1)+' '+str(Q48_17)+' '+str(Q48_17)+' '+str(1)+' 1 1 1 1)'+'\n')
#8
bmd.write('    hex (50 51 16 18 54 55 19 21) shape ('+str(cle/2)+' '+str(cv)+' '+str(int(-1*deltaZ/dzx))+ ') edgeGrading (1 '+str(Q17_16)+' '+str(Q17_16)+' 1 '+str(Q48_17)+' '+str(1)+' '+str(1)+' '+str(Q48_17)+' 1 1 1 1)'+'\n')
#9
bmd.write('    hex (48 22 24 17 52 26 28 20) shape ('+str(cle/4)+' '+str(cv)+' '+str(int(-1*deltaZ/dzx))+ ') edgeGrading (1 1 1 1 '+str(Q48_17)+' '+str(Q22_24)+' '+str(Q22_24)+' '+str(Q48_17)+' 1 1 1 1)'+'\n')
#10
bmd.write('    hex (23 50 18 25 27 54 21 29) shape ('+str(cle/4)+' '+str(cv)+' '+str(int(-1*deltaZ/dzx))+ ') edgeGrading (1 1 1 1 '+str(Q22_24)+' '+str(Q48_17)+' '+str(Q48_17)+' '+str(Q22_24)+' 1 1 1 1)'+'\n')
#11
bmd.write('    hex (22 30 34 24 26 32 36 28) shape ('+str(cd*2)+' '+str(cv)+' '+str(int(-1*deltaZ/dzx))+ ') edgeGrading ('+str(Q22_30/2)+' '+str(Q22_30/2)+' '+str(Q22_30/2)+' '+str(Q22_30/2)+' '+str(Q22_24)+' '+str(Q22_24)+' '+str(Q22_24)+' '+str(Q22_24)+' 1 1 1 1)'+'\n')
#12
bmd.write('    hex (23 25 35 31 27 29 37 33) shape ('+str(cv)+' '+str(cd*2)+' '+str(int(-1*deltaZ/dzx))+ ') edgeGrading ('+str(Q22_24)+' '+str(Q22_24)+' '+str(Q22_24)+' '+str(Q22_24)+' '+str(Q22_30/2)+' '+str(Q22_30/2)+' '+str(Q22_30/2)+' '+str(Q22_30/2)+' 1 1 1 1)'+'\n')
# bmd.write('    hex (49 38 30 22 53 39 32 26) shape ('+str(cd)+' '+str(cle/2)+' '+str(int(-1*deltaZ/dzx))+ ') edgeGrading ('+str(Q49_38)+' '+str(Q22_30)+' '+str(Q22_30)+' '+str(Q49_38)+' 1 1 1 1 1 1 1 1)'+'\n')
#13
bmd.write('    hex (49 38 102 101 53 39 105 104) shape ('+str(cd*2)+' '+str(cle/8)+' '+str(int(-1*deltaZ/dzx))+ ') edgeGrading ('+str(Q49_38/2)+' '+str(Q101_102/2)+' '+str(Q101_102/2)+' '+str(Q49_38/2)+' 1 1 1 1 1 1 1 1)'+'\n')
#14
bmd.write('    hex (101 102 30 22 104 105 32 26) shape ('+str(cd*2)+' '+str(cle/2)+' '+str(int(-1*deltaZ/dzx))+ ') edgeGrading ('+str(Q101_102/2)+' '+str(Q22_30/2)+' '+str(Q22_30/2)+' '+str(Q101_102/2)+' 1 1 1 1 1 1 1 1)'+'\n')
#15
bmd.write('    hex (49 23 31 38 53 27 33 39) shape ('+str(cle/2)+' '+str(cd*2)+' '+str(int(-1*deltaZ/dzx))+ ') edgeGrading (1 1 1 1 '+str(Q49_38/2)+' '+str(Q22_30/2)+' '+str(Q22_30/2)+' '+str(Q49_38/2)+' 1 1 1 1)'+'\n')

bmd.write(''+'\n')
bmd.write(');'+'\n')

bmd.write(''+'\n')
bmd.write('edges'+'\n')
bmd.write('('+'\n')
bmd.write(''+'\n')
bmd.write('arc 43 40 ('+str(2*c/5-mcyl*np.cos(0.7853981634))+' 0 '+str(mcyl*np.sin(0.7853981634))+')'+'\n')
bmd.write('arc 40 56 ('+str(2*c/5+mcyl*np.cos(1.570796327-0.2))+' 0 '+str(mcyl*np.sin(1.570796327-0.2))+')'+'\n')
# bmd.write('arc 56 41 ('+str(2*c/5+mcyl*np.cos(0.2))+' 0.5 '+str(mcyl*np.sin(0.2))+')'+'\n')
bmd.write('arc 56 100 ('+str(2*c/5+mcyl*np.cos(3*alpha64/4))+' 0 '+str(mcyl*np.sin(3*alpha64/4))+')'+'\n')
bmd.write('arc 100 41 ('+str(2*c/5+mcyl*np.cos(alpha64/4))+' 0 '+str(mcyl*np.sin(alpha64/4))+')'+'\n')

bmd.write('arc 41 57 ('+str(2*c/5+mcyl*np.cos(0.2))+' 0 '+str(-mcyl*np.sin(0.2))+')'+'\n')
bmd.write('arc 57 42 ('+str(2*c/5+mcyl*np.cos(1.570796327-0.2))+' 0 '+str(-mcyl*np.sin(1.570796327-0.2))+')'+'\n')
bmd.write('arc 42 43 ('+str(2*c/5-mcyl*np.cos(0.7853981634))+' 0 '+str(-mcyl*np.sin(0.7853981634))+')'+'\n')
bmd.write('arc 47 44 ('+str(2*c/5-mcyl*np.cos(0.7853981634))+ ' '+str(deltaZ)  +' '+str(mcyl*np.sin(0.7853981634))+')'+'\n')
bmd.write('arc 44 58 ('+str(2*c/5+mcyl*np.cos(1.570796327-0.2))+ ' '+str(deltaZ)  +' '+str(mcyl*np.sin(1.570796327-0.2))+')'+'\n')
# bmd.write('arc 58 45 ('+str(c/2+mcyl*np.cos(0.2))+ ' '+str(deltaZ)  +' '+str(mcyl*np.sin(0.2))+')'+'\n')
bmd.write('arc 58 103 ('+str(2*c/5+mcyl*np.cos(3*alpha64/4))+ ' '+str(deltaZ)  +' '+str(mcyl*np.sin(3*alpha64/4))+')'+'\n')
bmd.write('arc 103 45 ('+str(2*c/5+mcyl*np.cos(alpha64/4))+ ' '+str(deltaZ)  +' '+str(mcyl*np.sin(alpha64/4))+')'+'\n')

bmd.write('arc 45 59 ('+str(2*c/5+mcyl*np.cos(0.2))+ ' '+str(deltaZ)  +' '+str(-mcyl*np.sin(0.2))+')'+'\n')
bmd.write('arc 59 46 ('+str(2*c/5+mcyl*np.cos(1.570796327-0.2))+ ' '+str(deltaZ)  +' '+str(-mcyl*np.sin(1.570796327-0.2))+')'+'\n')
bmd.write('arc 46 47 ('+str(2*c/5-mcyl*np.cos(0.7853981634))+ ' '+str(deltaZ)  +' '+str(-mcyl*np.sin(0.7853981634))+')'+'\n')
bmd.write('arc 51 48 ('+str(2*c/5-lcyl*np.cos(0.7853981634))+' 0 '+str(lcyl*np.sin(0.7853981634))+')'+'\n')
bmd.write('arc 48 22 ('+str(2*c/5+lcyl*np.cos(1.570796327-0.2))+' 0 '+str(lcyl*np.sin(1.570796327-0.2))+')'+'\n')
# bmd.write('arc 22 49 ('+str(c/2+lcyl*np.cos(0.2))+' 0.5 '+str(lcyl*np.sin(0.2))+')'+'\n')
bmd.write('arc 22 101 ('+str(2*c/5+lcyl*np.cos(3*alpha64/4))+' 0 '+str(lcyl*np.sin(3*alpha64/4))+')'+'\n')
bmd.write('arc 101 49 ('+str(2*c/5+lcyl*np.cos(alpha64/4))+' 0 '+str(lcyl*np.sin(alpha64/4))+')'+'\n')

bmd.write('arc 49 23 ('+str(2*c/5+lcyl*np.cos(0.2))+' 0 '+str(-lcyl*np.sin(0.2))+')'+'\n')
bmd.write('arc 23 50 ('+str(2*c/5+lcyl*np.cos(1.570796327-0.2))+' 0 '+str(-lcyl*np.sin(1.570796327-0.2))+')'+'\n')
bmd.write('arc 50 51 ('+str(2*c/5-lcyl*np.cos(0.7853981634))+' 0 '+str(-lcyl*np.sin(0.7853981634))+')'+'\n')
bmd.write('arc 55 52 ('+str(2*c/5-lcyl*np.cos(0.7853981634))+ ' '+str(deltaZ)  +' '+str(lcyl*np.sin(0.7853981634))+')'+'\n')
bmd.write('arc 52 26 ('+str(2*c/5+lcyl*np.cos(1.570796327-0.2))+ ' '+str(deltaZ)  +' '+str(lcyl*np.sin(1.570796327-0.2))+')'+'\n')
# bmd.write('arc 26 53 ('+str(c/2+lcyl*np.cos(0.2))+ ' '+str(deltaZ)  +' '+str(lcyl*np.sin(0.2))+')'+'\n')
bmd.write('arc 26 104 ('+str(2*c/5+lcyl*np.cos(3*alpha64/4))+ ' '+str(deltaZ)  +' '+str(lcyl*np.sin(3*alpha64/4))+')'+'\n')
bmd.write('arc 104 53 ('+str(2*c/5+lcyl*np.cos(alpha64/4))+ ' '+str(deltaZ)  +' '+str(lcyl*np.sin(alpha64/4))+')'+'\n')

bmd.write('arc 53 27 ('+str(2*c/5+lcyl*np.cos(0.2))+ ' '+str(deltaZ)  +' '+str(-lcyl*np.sin(0.2))+')'+'\n')
bmd.write('arc 27 54 ('+str(2*c/5+lcyl*np.cos(1.570796327-0.2))+ ' '+str(deltaZ)  +' '+str(-lcyl*np.sin(1.570796327-0.2))+')'+'\n')
bmd.write('arc 54 55 ('+str(2*c/5-lcyl*np.cos(0.7853981634))+ ' '+str(deltaZ)  +' '+str(-lcyl*np.sin(0.7853981634))+')'+'\n')
bmd.write('arc 16 17 ('+str(-6*c*np.cos(0.7853981634))+' 0 '+str(+6*c*np.cos(0.7853981634))+')'+'\n')
bmd.write('arc 16 18 ('+str(-6*c*np.cos(0.7853981634))+' 0 '+str(-6*c*np.cos(0.7853981634))+')'+'\n')
bmd.write('arc 19 20 ('+str(-6*c*np.cos(0.7853981634))+ ' '+str(deltaZ)  +' '+str(+6*c*np.cos(0.7853981634))+')'+'\n')
bmd.write('arc 19 21 ('+str(-6*c*np.cos(0.7853981634))+ ' '+str(deltaZ)  +' '+str(-6*c*np.cos(0.7853981634))+')'+'\n')

bmd.write(''+'\n')
bmd.write(');'+'\n')

bmd.write(''+'\n')
bmd.write('boundary'+'\n')
bmd.write('('+'\n')
bmd.write(''+'\n')
bmd.write('	inlet'+'\n')
bmd.write('	{'+'\n')
bmd.write('		type patch;'+'\n')
bmd.write('		faces'+'\n')
bmd.write('		('+'\n')
bmd.write('		(17 20 19 16) '+'\n')
bmd.write('		(16 19 21 18)'+'\n')
bmd.write('		);'+'\n')
bmd.write('	}'+'\n')
bmd.write(''+'\n')
bmd.write('	outlet'+'\n')
bmd.write('	{'+'\n')
bmd.write('		type patch;'+'\n')
bmd.write('		faces'+'\n')
bmd.write('		('+'\n')
bmd.write('		(38 31 33 39)'+'\n')
bmd.write('		(34 30 32 36)'+'\n')
# bmd.write('		(30 38 39 32)'+'\n')
bmd.write('		(30 102 105 32)'+'\n')
bmd.write('		(102 38 39 105)'+'\n')

bmd.write('		(31 35 37 33)'+'\n')
bmd.write('		);'+'\n')
bmd.write('	}'+'\n')
bmd.write(''+'\n')
bmd.write('	topAndBottom'+'\n')
bmd.write('	{'+'\n')
bmd.write('		type patch;'+'\n')
bmd.write('		faces'+'\n')
bmd.write('		('+'\n')
bmd.write('		(24 34 36 28)'+'\n')
bmd.write('		(17 24 28 20)'+'\n')
bmd.write('		(18 21 29 25)'+'\n')
bmd.write('		(25 29 37 35)'+'\n')
bmd.write('		);'+'\n')
bmd.write('	}'+'\n')
bmd.write(''+'\n')
bmd.write('	rotbox'+'\n')
bmd.write('	{'+'\n')
bmd.write('		type patch;'+'\n')
bmd.write('		faces'+'\n')
bmd.write('		('+'\n')
bmd.write('		(43 40 44 47)'+'\n')
bmd.write('		(40 56 58 44)'+'\n')
# bmd.write('		(56 41 45 58)'+'\n')
bmd.write('		(56 100 103 58)'+'\n')
bmd.write('		(100 41 45 103)'+'\n')

bmd.write('		(41 57 59 45)'+'\n')
bmd.write('		(57 42 46 59)'+'\n')
bmd.write('		(42 43 47 46)'+'\n')
bmd.write('		);'+'\n')
bmd.write('	}'+'\n')
bmd.write(''+'\n')
bmd.write('	front'+'\n')
bmd.write('	{'+'\n')
bmd.write('		type empty;'+'\n')
bmd.write('		faces'+'\n')
bmd.write('		('+'\n')
bmd.write('		(16 51 48 17)'+'\n')
bmd.write('		(16 18 50 51)'+'\n')
bmd.write('		(51 43 40 48)'+'\n')
bmd.write('		(48 22 24 17)'+'\n')
bmd.write('		(40 56 22 48)'+'\n')
# bmd.write('		(56 41 49 22)'+'\n')
bmd.write('		(56 100 101 22)'+'\n')
bmd.write('		(100 41 49 101)'+'\n')

bmd.write('		(18 25 23 50)'+'\n')
bmd.write('		(50 23 57 42)'+'\n')
bmd.write('		(23 49 41 57)'+'\n')
bmd.write('		(51 50 42 43)'+'\n')
bmd.write('		(22 30 34 24)'+'\n')
# bmd.write('		(49 38 30 22)'+'\n')
bmd.write('		(49 38 102 101)'+'\n')
bmd.write('		(101 102 30 22)'+'\n')

bmd.write('		(23 31 38 49)'+'\n')
bmd.write('		(25 35 31 23)'+'\n')
bmd.write('		);'+'\n')
bmd.write('	}'+'\n')
bmd.write(''+'\n')
bmd.write('	back'+'\n')
bmd.write('	{'+'\n')
bmd.write('		type empty;'+'\n')
bmd.write('		faces'+'\n')
bmd.write('		('+'\n')
bmd.write('		(19 20 52 55)'+'\n')
bmd.write('		(21 19 55 54)'+'\n')
bmd.write('		(47 55 52 44)'+'\n')
bmd.write('		(54 55 47 46)'+'\n')
bmd.write('		(26 52 20 28)'+'\n')
bmd.write('		(58 44 52 26)'+'\n')
# bmd.write('		(45 58 26 53)'+'\n')
bmd.write('		(45 103 104 53)'+'\n')
bmd.write('		(103 58 26 104)'+'\n')

bmd.write('		(21 54 27 29)'+'\n')
bmd.write('		(27 54 46 59)'+'\n')
bmd.write('		(53 27 59 45)'+'\n')
bmd.write('		(32 26 28 36)'+'\n')
# bmd.write('		(39 53 26 32)'+'\n')
bmd.write('		(39 53 104 105)'+'\n')
bmd.write('		(105 104 26 32)'+'\n')

bmd.write('		(33 27 53 39)'+'\n')
bmd.write('		(37 29 27 33)'+'\n')
bmd.write('		);'+'\n')
bmd.write('	}'+'\n')
bmd.write(''+'\n')
bmd.write(');'+'\n')

bmd.write(''+'\n')
bmd.write('mergePatchPairs'+'\n')
bmd.write('('+'\n')
bmd.write(''+'\n')
bmd.write(''+'\n')
bmd.write(');'+'\n')
bmd.write(''+'\n')
bmd.write('// ************************************************************************* //'+'\n')








bmd.close()
