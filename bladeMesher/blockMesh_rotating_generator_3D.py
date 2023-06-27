#!/usr/bin/env python3 
# -*- coding: utf-8 -*- 
#----------------------------------------------------------------------------
# Created By  : Javier Anez & Julien Reveillon
# Created Date: 2022
# ---------------------------------------------------------------------------
""" generate blockMesh for external aero profiles """ 
# How to use :  blockMesh_rotating_generator_3D.py fileWithXYZprofile.txt
# ---------------------------------------------------------------------------


import os,sys,shutil, fileinput
import math
import matplotlib.pyplot as plt

import numpy as np




currentDir = os.getcwd()+'/'

#Variables
points = np.zeros((110,3))
nb_points = 105
c=0.6 #chord
scyl = c + c/3
mcyl = scyl+c/10
lcyl = mcyl+c/10
L = 0.046*c # Distance between points 0 and 82 to define the layer

#==================== START OF FIRST LAYER PARAMETER REGION ============
#Kinematic viscosity
nu = 1.5e-05

#Reynolds
Re = 1600000

#y+ in first layer
yPlus = 50

# Calculation of first layer height
# Reference : https://www.computationalfluiddynamics.com.au/tips-tricks-cfd-estimate-first-cell-height/

#1) Compute velocity
U = nu*Re/c
print('U =', U)

#2) Compute Cf
#For external flows Cf = 0.058 Re^{-0.20}
C_f = 0.058 *Re**(-0.20)

#3) Compute tau_w / rho
Tau_w_rho = 0.5 * C_f * U**2

#4) Compute U_tau
U_tau = Tau_w_rho**(0.5)

#5) Compute first layer height
y1 = yPlus * nu / U_tau
# print('y1 =', y1)

dl0 = y1 #Heigth of the smallest cell around the profile
# print('dl0 =', dl0)
# sys.exit()

clayer = int(round(L/dl0)) #Number of cells in the layer around the profil
cle = 180 # number of cells on the leading edge of the profile (between points 0 and 1)
cte = 200 # number of cells on the trailing edge of the profile (between points 1 and 3)
smc = 20 # number of cells between the smaller and the medium cylinder

#Get the profil coordinates

spline = spline = sys.argv[1]
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
plt.plot(xu,zu)
plt.plot(xl,zl)
plt.plot(x2u,z2u)
plt.plot(x2l,z2l)
autoRunDir = os.getcwd()+'/'

#SAve the blade profile
outputFile = autoRunDir+'profileLayers.png'
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

#shape
counter = 0
z0 = zu[counter]
points[0,:]=[0.0,1.4,z0] #profile noze

counter = 0
z1 = 0.0
for i in xu:
	if i < c/2:
		z1= (zu[counter]+zu[counter+1])/2
		x1= (xu[counter]+xu[counter+1])/2
		# z1= zu[counter]
	counter += 1
points[1,:]=[x1,1.4,z1]

counter = 0
z2 = 0.0
for i in xl:
	if i < c/2:
		z2= (zl[counter]+zl[counter+1])/2
		x2= (xl[counter]+xl[counter+1])/2
		# z2= zl[counter]
	counter += 1
points[2,:]=[x2,1.4,z2]

counter = 0
z3 = 0.0
for i in xl:
	if i <= c:
		# z61= (z2u[counter]+z2u[counter+1])/2
		z3= zl[counter]
	counter += 1
points[3,:]=[c,1.4,z3]

points[4,:]=[0.0,-0.5,z0]
points[5,:]=[x1,-0.5,z1]
points[6,:]=[x2,-0.5,z2]
points[7,:]=[c,-0.5,z3]
#end shape

#smaller cylinder
points[8,:]=[x1,1.4,scyl]
points[9,:]=[c/2+scyl,1.4,z3]
points[10,:]=[x2,1.4,-scyl]
points[11,:]=[c/2-scyl,1.4,z0]
points[12,:]=[x1,-0.5,scyl]
points[13,:]=[c/2+scyl,-0.5,z3]
points[14,:]=[x2,-0.5,-scyl]
points[15,:]=[c/2-scyl,-0.5,z0]

#Inlet
points[16,:]=[-13*c,1.4,z0]
points[17,:]=[0.0,1.4,13*c]
points[18,:]=[0.0,1.4,-13*c]
points[19,:]=[-13*c,-0.5,z0]
points[20,:]=[0.0,-0.5,13*c]
points[21,:]=[0.0,-0.5,-13*c]

#Points downstream
points[22,:]=[c/2+lcyl*np.cos(0.7854),1.4,lcyl*np.sin(0.7854)]
points[23,:]=[c/2+lcyl*np.cos(0.7854),1.4,-lcyl*np.sin(0.7854)]
points[24,:]=[c/2+lcyl*np.cos(0.7854),1.4,13*c]
points[25,:]=[c/2+lcyl*np.cos(0.7854),1.4,-13*c]
points[26,:]=[c/2+lcyl*np.cos(0.7854),-0.5,lcyl*np.sin(0.7854)]
points[27,:]=[c/2+lcyl*np.cos(0.7854),-0.5,-lcyl*np.sin(0.7854)]
points[28,:]=[c/2+lcyl*np.cos(0.7854),-0.5,13*c]
points[29,:]=[c/2+lcyl*np.cos(0.7854),-0.5,-13*c]
points[30,:]=[27*c,1.4,lcyl*np.sin(0.7854)]
points[31,:]=[27*c,1.4,-lcyl*np.sin(0.7854)]
points[32,:]=[27*c,-0.5,lcyl*np.sin(0.7854)]
points[33,:]=[27*c,-0.5,-lcyl*np.sin(0.7854)]
points[34,:]=[27*c,1.4,13*c]
points[35,:]=[27*c,1.4,-13*c]
points[36,:]=[27*c,-0.5,13*c]
points[37,:]=[27*c,-0.5,-13*c]
points[38,:]=[27*c,1.4,z3]
points[39,:]=[27*c,-0.5,z3]

#medium cylinder
points[40,:]=[x1,1.4,mcyl]
points[41,:]=[c/2+mcyl,1.4,z3]
points[42,:]=[x2,1.4,-mcyl]
points[43,:]=[c/2-mcyl,1.4,z0]
points[44,:]=[x1,-0.5,mcyl]
points[45,:]=[c/2+mcyl,-0.5,z3]
points[46,:]=[x2,-0.5,-mcyl]
points[47,:]=[c/2-mcyl,-0.5,z0]

#larger cylinder
points[48,:]=[c/2,1.4,lcyl]
points[49,:]=[c/2+lcyl,1.4,z3]
points[50,:]=[c/2,1.4,-lcyl]
points[51,:]=[c/2-lcyl,1.4,z0]
points[52,:]=[c/2,-0.5,lcyl]
points[53,:]=[c/2+lcyl,-0.5,z3]
points[54,:]=[c/2,-0.5,-lcyl]
points[55,:]=[c/2-lcyl,-0.5,z0]

#Points at pi/4 on medium cylider
points[56,:]=[c/2+mcyl*np.cos(0.7854),1.4,mcyl*np.sin(0.7854)]
points[57,:]=[c/2+mcyl*np.cos(0.7854),1.4,-mcyl*np.sin(0.7854)]
points[58,:]=[c/2+mcyl*np.cos(0.7854),-0.5,mcyl*np.sin(0.7854)]
points[59,:]=[c/2+mcyl*np.cos(0.7854),-0.5,-mcyl*np.sin(0.7854)]

#Points in the rotating part
alpha60 = np.arccos((c/2)/scyl)
points[60,:]=[c,1.4,scyl*np.sin(alpha60)]

counter = 0
z61 = 0.0
for i in x2u:
	if i <= c:
		# z61= (z2u[counter]+z2u[counter+1])/2
		z61= z2u[counter]
	counter += 1
points[61,:]=[c,1.4,z61]


counter = 0
z62 = 0.0
for i in x2l:
	if i <= c:
		# z62= (z2l[counter]+z2l[counter-1])/2
		z62= z2l[counter]
	counter += 1
points[62,:]=[c,1.4,z62]

points[63,:]=[c,1.4,-scyl*np.sin(alpha60)]
points[65,:]=[c+4*L,1.4,z61]

#placement of the point 64 to minimize the non-orthogonality
#angle between (65,66) and (65,64) should be 135 degrees
#point 64 is at the intersection of the circle (z=scyl*sin(alpha64) and x=c/2+scyl*cos(alpha64) and the straight line with an angle of 45 degrees (z=ax+b and a=1)
#This lead to an quadratic equation
b=points[65,2]-points[65,0]
delta =(2*(c/2+b)/scyl)**2-4*2*(((c/2+b)/scyl)**2-1)
alpha64 = np.arccos((-2*(c/2+b)/scyl+np.sqrt(delta))/(2*2))

# print('z61 =', z61)
# print('c/2 + 4*L =', c/2 + 4*L)
# jap = z61 / (c/2 + 4*L))
# print('alpha64_jap =', np.arctan( jap )* (180/np.pi))
# print('alpha64_jap =', np.arctan( (scyl*np.sin(alpha64) ) / (scyl*np.cos(alpha64)) ))
# print('alpha64 =', alpha64)
# sys.exit()

points[64,:]=[c/2+scyl*np.cos(alpha64),1.4,scyl*np.sin(alpha64)]
points[66,:]=[c+4*L,1.4,z3]
points[67,:]=[c+4*L,1.4,z62]
points[68,:]=[c/2+scyl*np.cos(alpha64),1.4,-scyl*np.sin(alpha64)]
points[69,:]=[c,-0.5,scyl*np.sin(alpha60)]
points[70,:]=[c,-0.5,z61]
points[71,:]=[c,-0.5,z62]
points[72,:]=[c,-0.5,-scyl*np.sin(alpha60)]
points[73,:]=[c/2+scyl*np.cos(alpha64),-0.5,scyl*np.sin(alpha64)]
points[74,:]=[c+4*L,-0.5,z61]
points[75,:]=[c+4*L,-0.5,z3]
points[76,:]=[c+4*L,-0.5,z62]
points[77,:]=[c/2+scyl*np.cos(alpha64),-0.5,-scyl*np.sin(alpha64)]

counter = 0
z78 = 0.0
for i in x2u:
	if i <= c/2:
		z78= (z2u[counter]+z2u[counter+1])/2
	counter += 1
points[78,:]=[x1,1.4,z78]

counter = 0
z79 = 0.0
for i in x2l:
	if i <= c/2:
		z79= (z2l[counter]+z2l[counter+1])/2
	counter += 1
points[79,:]=[x2,1.4,z79]

points[80,:]=[x1,-0.5,z78]
points[81,:]=[x2,-0.5,z79]
points[82,:]=[points[0,0]-L,1.4,z0]
points[83,:]=[points[0,0]-L,-0.5,z0]
points[84,:]=[c/2+mcyl*np.cos(alpha60),1.4,mcyl*np.sin(alpha60)]
points[85,:]=[c/2+mcyl*np.cos(alpha64),1.4,mcyl*np.sin(alpha64)]
points[86,:]=[c/2+mcyl*np.cos(alpha64),1.4,-mcyl*np.sin(alpha64)]
points[87,:]=[c/2+mcyl*np.cos(alpha60),1.4,-mcyl*np.sin(alpha60)]
points[88,:]=[c/2+mcyl*np.cos(alpha60),-0.5,mcyl*np.sin(alpha60)]
points[89,:]=[c/2+mcyl*np.cos(alpha64),-0.5,mcyl*np.sin(alpha64)]
points[90,:]=[c/2+mcyl*np.cos(alpha64),-0.5,-mcyl*np.sin(alpha64)]
points[91,:]=[c/2+mcyl*np.cos(alpha60),-0.5,-mcyl*np.sin(alpha60)]

# new points fo rthe flatBack profile
counter = 0
z92 = 0.0
for i in xu:
	if i <= c:
		# z62= (z2l[counter]+z2l[counter-1])/2
		z92= zu[counter]
	counter += 1
points[92,:]=[c,1.4,z92]
points[93,:]=[c+4*L,1.4,z92]
points[94,:]=[c/2+scyl*np.cos(alpha64/8) ,1.4, scyl*np.sin(alpha64/8)]
points[95,:]=[c/2+mcyl*np.cos(alpha64/8) ,1.4, mcyl*np.sin(alpha64/8)]
points[96,:]=[c,-0.5,z92]
points[97,:]=[c+4*L,-0.5,z92]
points[98,:]=[c/2+scyl*np.cos(alpha64/8) ,-0.5, scyl*np.sin(alpha64/8)]
points[99,:]=[c/2+mcyl*np.cos(alpha64/8) ,-0.5, mcyl*np.sin(alpha64/8)]

points[100,:]=[c/2+mcyl*np.cos(alpha64/2),1.4,mcyl*np.sin(alpha64/2)]
points[101,:]=[c/2+lcyl*np.cos(alpha64/2),1.4,lcyl*np.sin(alpha64/2)]
points[102,:]=[27*c,1.4,lcyl*np.sin(alpha64/2)]
points[103,:]=[c/2+mcyl*np.cos(alpha64/2),-0.5,mcyl*np.sin(alpha64/2)]
points[104,:]=[c/2+lcyl*np.cos(alpha64/2),-0.5,lcyl*np.sin(alpha64/2)]
points[105,:]=[27*c,-0.5,lcyl*np.sin(alpha64/2)]
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

cvd=Q_N(points[8,2]-points[78,2],dl0,2*3.14159*scyl/(4*cle))[1] #Number of cells in the vertical direction (between points 78 and 8)
Qvd=Q_N(points[8,2]-points[78,2],dl0,2*3.14159*scyl/(4*cle))[2] #Grading in the vertical direction (between points 78 and 8)
dlNvd=Q_N(points[8,2]-points[78,2],dl0,2*3.14159*scyl/(4*cle))[5]
cds=Q_N(points[66,0]-points[3,0],np.sqrt((points[92,0]-points[1,0])**2+(points[1,2]-points[92,2])**2)/cte,dl0)[1] #number of cells downstream the shape (between points 3 and 66)

Qds=Q_N(points[66,0]-points[3,0],np.sqrt((points[92,0]-points[1,0])**2+(points[1,2]-points[92,2])**2)/cte,dl0)[2]  #Grading downstream the shape (between points 3 and 66)

f11_82=Q(points[82,0]-points[11,0],dl0,cvd)
Q11_82=f11_82[2]
dlN11_82=f11_82[5]

f61_60=Q(points[60,2]-points[61,2],dl0,cvd)
Q61_60=f61_60[2]
dlN61_60=f61_60[5]

f65_64=Q(np.sqrt((points[64,2]-points[65,2])**2+(points[64,0]-points[65,0])**2),dl0,cvd)
Q65_64=f65_64[2]
dlN65_64=f65_64[5]

f93_94=Q(np.sqrt((points[94,2]-points[93,2])**2+(points[94,0]-points[93,0])**2),dl0,cvd)
Q93_94=f93_94[2]
dlN93_94=f93_94[5]

f66_9=Q(points[9,0]-points[66,0],dl0,cvd)
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

Q11_43=Q(mcyl-scyl,dlN11_82,smc)[2]
Q8_40=Q(mcyl-scyl,dlNvd,smc)[2]
Q60_84=Q(mcyl-scyl,dlN61_60,smc)[2]
Q64_85=Q(mcyl-scyl,dlN65_64,smc)[2]
Q9_41=Q(mcyl-scyl,dlN66_9,smc)[2]
Q94_95=Q(mcyl-scyl,dlN93_94,smc)[2]



#Write the blockMeshDict File

bmd = open("blockMeshDict","w")

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
bmd.write('    hex (0 1 78 82 4 5 80 83) shape ('+str(cle)+' '+str(clayer)+' 10) edgeGrading (1 1 1 1 1 1 1 1 1 1 1 1)'+'\n')
bmd.write('    hex (82 78 8 11 83 80 12 15) shape ('+str(cle)+' '+str(cvd)+' 10) edgeGrading (1 '+str(1/Q8_11)+' '+str(1/Q8_11)+' 1 '+str(Q11_82)+' '+str(Qvd)+' '+str(Qvd)+' '+str(Q11_82)+' 1 1 1 1)'+'\n')
bmd.write('    hex (1 92 61 78 5 96 70 80) shape ('+str(cte)+' '+str(clayer)+' 10) edgeGrading (1 1 1 1 1 1 1 1 1 1 1 1)'+'\n')
bmd.write('    hex (78 61 60 8 80 70 69 12) shape ('+str(cte)+' '+str(cvd)+' 10) edgeGrading (1 1 1 1 '+str(Qvd)+' '+str(Q61_60)+' '+str(Q61_60)+' '+str(Qvd)+' 1 1 1 1)'+'\n')
bmd.write('    hex (92 93 65 61 96 97 74 70) shape ('+str(cds)+' '+str(clayer)+' 10) edgeGrading ('+str(Qds)+' '+str(Qds)+' '+str(Qds)+' '+str(Qds)+' 1 1 1 1 1 1 1 1)'+'\n')
#5
bmd.write('    hex (3 66 93 92 7 75 97 96) shape ('+str(cds)+' '+str(3*clayer/8)+' 10) edgeGrading ('+str(Qds)+' '+str(Qds)+' '+str(Qds)+' '+str(Qds)+' 1 1 1 1 1 1 1 1)'+'\n')
bmd.write('    hex (61 65 64 60 70 74 73 69) shape ('+str(cds)+' '+str(cvd)+' 10) edgeGrading ('+str(Qds)+' '+str(Q60_64)+' '+str(Q60_64)+' '+str(Qds)+' '+str(Q61_60)+' '+str(Q65_64)+' '+str(Q65_64)+' '+str(Q61_60)+' 1 1 1 1)'+'\n')
#7
# bmd.write('    hex (93 94 64 65 97 98 73 74) shape ('+str(cvd)+' '+str(clayer)+' 1) edgeGrading ('+str(Q93_94)+' '+str(Q65_64)+' '+str(Q65_64)+' '+str(Q93_94)+' 1 '+str(1/Q64_94)+' '+str(1/Q64_94)+' 1 1 1 1 1)'+'\n')
bmd.write('    hex (93 94 64 65 97 98 73 74) shape ('+str(cvd)+' '+str(clayer)+' 10) edgeGrading ('+str(Q93_94)+' '+str(Q65_64)+' '+str(Q65_64)+' '+str(Q93_94)+' 1 '+str(1/1)+' '+str(1/1)+' 1 1 1 1 1)'+'\n')

#8
bmd.write('    hex (66 9 94 93 75 13 98 97) shape ('+str(cvd)+' '+str(3*clayer/8)+' 10) edgeGrading ('+str(Q66_9)+' '+str(Q93_94)+' '+str(Q93_94)+' '+str(Q66_9)+' 1 '+str(1/1)+' '+str(1/1)+' 1 1 1 1 1)'+'\n')
bmd.write('    hex (2 0 82 79 6 4 83 81) shape ('+str(cle)+' '+str(clayer)+' 10) edgeGrading (1 1 1 1 1 1 1 1 1 1 1 1)'+'\n')
bmd.write('    hex (79 82 11 10 81 83 15 14) shape ('+str(cle)+' '+str(cvd)+' 10) edgeGrading (1 '+str(Q8_11)+' '+str(Q8_11)+' 1 '+str(Qvd)+' '+str(Q11_82)+' '+str(Q11_82)+' '+str(Qvd)+' 1 1 1 1)'+'\n')
bmd.write('    hex (3 2 79 62 7 6 81 71) shape ('+str(cte)+' '+str(clayer)+' 10) edgeGrading (1 1 1 1 1 1 1 1 1 1 1 1)'+'\n')
bmd.write('    hex (62 79 10 63 71 81 14 72) shape ('+str(cte)+' '+str(cvd)+' 10) edgeGrading (1 1 1 1 '+str(Q61_60)+' '+str(Qvd)+' '+str(Qvd)+' '+str(Q61_60)+' 1 1 1 1)'+'\n')
bmd.write('    hex (66 3 62 67 75 7 71 76) shape ('+str(cds)+' '+str(clayer)+' 10) edgeGrading ('+str(1/Qds)+' '+str(1/Qds)+' '+str(1/Qds)+' '+str(1/Qds)+' 1 1 1 1 1 1 1 1)'+'\n')
bmd.write('    hex (67 62 63 68 76 71 72 77) shape ('+str(cds)+' '+str(cvd)+' 10) edgeGrading ('+str(1/Qds)+' '+str(1/Q60_64)+' '+str(1/Q60_64)+' '+str(1/Qds)+' '+str(Q65_64)+' '+str(Q61_60)+' '+str(Q61_60)+' '+str(Q65_64)+' 1 1 1 1)'+'\n')
#15
bmd.write('    hex (9 66 67 68 13 75 76 77) shape ('+str(cvd)+' '+str(clayer)+' 10) edgeGrading ('+str(1/Q66_9)+' '+str(1/Q65_64)+' '+str(1/Q65_64)+' '+str(1/Q66_9)+' '+str(1/Q64_9)+' 1 1 '+str(1/Q64_9)+' 1 1 1 1)'+'\n')
bmd.write('    hex (11 8 40 43 15 12 44 47) shape ('+str(cle)+' '+str(smc)+' 10) edgeGrading ('+str(1/Q8_11)+' '+str(1/Q8_11)+' '+str(1/Q8_11)+' '+str(1/Q8_11)+' '+str(Q11_43)+' '+str(Q8_40)+' '+str(Q8_40)+' '+str(Q11_43)+' 1 1 1 1)'+'\n')
bmd.write('    hex (8 60 84 40 12 69 88 44) shape ('+str(cte)+' '+str(smc)+' 10) edgeGrading (1 1 1 1 '+str(Q8_40)+' '+str(Q60_84)+' '+str(Q60_84)+' '+str(Q8_40)+' 1 1 1 1)'+'\n')
bmd.write('    hex (60 64 85 84 69 73 89 88) shape ('+str(cds)+' '+str(smc)+' 10) edgeGrading ('+str(Q60_64)+' '+str(Q60_64)+' '+str(Q60_64)+' '+str(Q60_64)+' '+str(Q60_84)+' '+str(Q64_85)+' '+str(Q64_85)+' '+str(Q60_84)+' 1 1 1 1)'+'\n')
#19
bmd.write('    hex (94 9 41 95 98 13 45 99) shape ('+str(3*clayer/8)+' '+str(smc)+' 10) edgeGrading ('+str(1)+' '+str(1)+' '+str(1)+' '+str(1)+' '+str(Q94_95)+' '+str(Q9_41)+' '+str(Q9_41)+' '+str(Q94_95)+' 1 1 1 1)'+'\n')
#20
# bmd.write('    hex (64 94 95 85 73 98 99 89) shape ('+str(clayer)+' '+str(smc)+' 1) edgeGrading ('+str(Q64_94)+' '+str(Q64_94)+' '+str(Q64_94)+' '+str(Q64_94)+' '+str(Q64_85)+' '+str(Q94_95)+' '+str(Q94_95)+' '+str(Q64_85)+' 1 1 1 1)'+'\n')
bmd.write('    hex (64 94 95 85 73 98 99 89) shape ('+str(clayer)+' '+str(smc)+' 10) edgeGrading ('+str(1)+' '+str(1)+' '+str(1)+' '+str(1)+' '+str(Q64_85)+' '+str(Q94_95)+' '+str(Q94_95)+' '+str(Q64_85)+' 1 1 1 1)'+'\n')

bmd.write('    hex (9 68 86 41 13 77 90 45) shape ('+str(clayer)+' '+str(smc)+' 10) edgeGrading ('+str(1/Q64_9)+' '+str(1/Q64_9)+' '+str(1/Q64_9)+' '+str(1/Q64_9)+' '+str(Q9_41)+' '+str(Q64_85)+' '+str(Q64_85)+' '+str(Q9_41)+' 1 1 1 1)'+'\n')
bmd.write('    hex (68 63 87 86 77 72 91 90) shape ('+str(cds)+' '+str(smc)+' 10) edgeGrading ('+str(1/Q60_64)+' '+str(1/Q60_64)+' '+str(1/Q60_64)+' '+str(1/Q60_64)+' '+str(Q64_85)+' '+str(Q60_84)+' '+str(Q60_84)+' '+str(Q64_85)+' 1 1 1 1)'+'\n')
bmd.write('    hex (63 10 42 87 72 14 46 91) shape ('+str(cte)+' '+str(smc)+' 10) edgeGrading (1 1 1 1 '+str(Q60_84)+' '+str(Q8_40)+' '+str(Q8_40)+' '+str(Q60_84)+' 1 1 1 1)'+'\n')
bmd.write('    hex (10 11 43 42 14 15 47 46) shape ('+str(cle)+' '+str(smc)+' 10) edgeGrading ('+str(Q8_11)+' '+str(Q8_11)+' '+str(Q8_11)+' '+str(Q8_11)+' '+str(Q8_40)+' '+str(Q11_43)+' '+str(Q11_43)+' '+str(Q8_40)+' 1 1 1 1)'+'\n')
bmd.write(''+'\n')
bmd.write(');'+'\n')

bmd.write(''+'\n')
bmd.write('edges'+'\n')
bmd.write('('+'\n')
bmd.write(''+'\n')
bmd.write('arc 11 8 ('+str(c/2-scyl*np.cos(0.7853981634))+' 1.4 '+str(scyl*np.sin(0.7853981634))+')'+'\n')
bmd.write('arc 8 60 ('+str(c/2+scyl*np.cos(1.570796327-0.2))+' 1.4 '+str(scyl*np.sin(1.570796327-0.2))+')'+'\n')
bmd.write('arc 60 64 ('+str(c/2+scyl*np.cos(alpha60-0.2))+' 1.4 '+str(scyl*np.sin(alpha60-0.2))+')'+'\n')
# bmd.write('arc 64 9 ('+str(c/2+scyl*np.cos(0.2))+' 0.5 '+str(scyl*np.sin(0.2))+')'+'\n')
bmd.write('arc 64 94 ('+str(c/2+scyl*np.cos(5*alpha64/8))+' 1.4 '+str(scyl*np.sin(5*alpha64/8))+')'+'\n')
bmd.write('arc 94 9 ('+str(c/2+scyl*np.cos(alpha64/16))+' 1.4 '+str(scyl*np.sin(alpha64/16))+')'+'\n')

bmd.write('arc 9 68 ('+str(c/2+scyl*np.cos(0.2))+' 1.4 '+str(-scyl*np.sin(0.2))+')'+'\n')
bmd.write('arc 68 63 ('+str(c/2+scyl*np.cos(alpha60-0.2))+' 1.4 '+str(-scyl*np.sin(alpha60-0.2))+')'+'\n')
bmd.write('arc 63 10 ('+str(c/2+scyl*np.cos(1.570796327-0.2))+' 1.4 '+str(-scyl*np.sin(1.570796327-0.2))+')'+'\n')
bmd.write('arc 10 11 ('+str(c/2-scyl*np.cos(0.7853981634))+' 1.4 '+str(-scyl*np.sin(0.7853981634))+')'+'\n')
bmd.write('arc 15 12 ('+str(c/2-scyl*np.cos(0.7853981634))+' -0.5 '+str(scyl*np.sin(0.7853981634))+')'+'\n')
bmd.write('arc 12 69 ('+str(c/2+scyl*np.cos(1.570796327-0.2))+' -0.5 '+str(scyl*np.sin(1.570796327-0.2))+')'+'\n')
bmd.write('arc 69 73 ('+str(c/2+scyl*np.cos(alpha60-0.2))+' -0.5 '+str(scyl*np.sin(alpha60-0.2))+')'+'\n')
# bmd.write('arc 73 13 ('+str(c/2+scyl*np.cos(0.2))+' -0.5 '+str(scyl*np.sin(0.2))+')'+'\n')
bmd.write('arc 73 98 ('+str(c/2+scyl*np.cos(5*alpha64/8))+' -0.5 '+str(scyl*np.sin(5*alpha64/8))+')'+'\n')
bmd.write('arc 98 13 ('+str(c/2+scyl*np.cos(alpha64/16))+' -0.5 '+str(scyl*np.sin(alpha64/16))+')'+'\n')

bmd.write('arc 13 77 ('+str(c/2+scyl*np.cos(0.2))+' -0.5 '+str(-scyl*np.sin(0.2))+')'+'\n')
bmd.write('arc 77 72 ('+str(c/2+scyl*np.cos(alpha60-0.2))+' -0.5 '+str(-scyl*np.sin(alpha60-0.2))+')'+'\n')
bmd.write('arc 72 14 ('+str(c/2+scyl*np.cos(1.570796327-0.2))+' -0.5 '+str(-scyl*np.sin(1.570796327-0.2))+')'+'\n')
bmd.write('arc 14 15 ('+str(c/2-scyl*np.cos(0.7853981634))+' -0.5 '+str(-scyl*np.sin(0.7853981634))+')'+'\n')

bmd.write('arc 43 40 ('+str(c/2-mcyl*np.cos(0.7853981634))+' 1.4 '+str(mcyl*np.sin(0.7853981634))+')'+'\n')
bmd.write('arc 40 84 ('+str(c/2+mcyl*np.cos(1.570796327-0.2))+' 1.4 '+str(mcyl*np.sin(1.570796327-0.2))+')'+'\n')
bmd.write('arc 84 85 ('+str(c/2+mcyl*np.cos(alpha60-0.2))+' 1.4 '+str(mcyl*np.sin(alpha60-0.2))+')'+'\n')
# bmd.write('arc 85 41 ('+str(c/2+mcyl*np.cos(0.2))+' 0.5 '+str(mcyl*np.sin(0.2))+')'+'\n')
bmd.write('arc 85 95 ('+str(c/2+mcyl*np.cos(5*alpha64/8))+' 1.4 '+str(mcyl*np.sin(5*alpha64/8))+')'+'\n')
bmd.write('arc 95 41 ('+str(c/2+mcyl*np.cos(alpha64/16))+' 1.4 '+str(mcyl*np.sin(alpha64/16))+')'+'\n')

bmd.write('arc 41 86 ('+str(c/2+mcyl*np.cos(0.2))+' 1.4 '+str(-mcyl*np.sin(0.2))+')'+'\n')
bmd.write('arc 86 87 ('+str(c/2+mcyl*np.cos(alpha60-0.2))+' 1.4 '+str(-mcyl*np.sin(alpha60-0.2))+')'+'\n')
bmd.write('arc 87 42 ('+str(c/2+mcyl*np.cos(1.570796327-0.2))+' 1.4 '+str(-mcyl*np.sin(1.570796327-0.2))+')'+'\n')
bmd.write('arc 42 43 ('+str(c/2-mcyl*np.cos(0.7853981634))+' 1.4 '+str(-mcyl*np.sin(0.7853981634))+')'+'\n')
bmd.write('arc 47 44 ('+str(c/2-mcyl*np.cos(0.7853981634))+' -0.5 '+str(mcyl*np.sin(0.7853981634))+')'+'\n')
bmd.write('arc 44 88 ('+str(c/2+mcyl*np.cos(1.570796327-0.2))+' -0.5 '+str(mcyl*np.sin(1.570796327-0.2))+')'+'\n')
bmd.write('arc 88 89 ('+str(c/2+mcyl*np.cos(alpha60-0.2))+' -0.5 '+str(mcyl*np.sin(alpha60-0.2))+')'+'\n')
# bmd.write('arc 89 45 ('+str(c/2+mcyl*np.cos(0.2))+' -0.5 '+str(mcyl*np.sin(0.2))+')'+'\n')
bmd.write('arc 89 99 ('+str(c/2+mcyl*np.cos(5*alpha64/8))+' -0.5 '+str(mcyl*np.sin(5*alpha64/8))+')'+'\n')
bmd.write('arc 99 45 ('+str(c/2+mcyl*np.cos(alpha64/16))+' -0.5 '+str(mcyl*np.sin(alpha64/16))+')'+'\n')

bmd.write('arc 45 90 ('+str(c/2+mcyl*np.cos(0.2))+' -0.5 '+str(-mcyl*np.sin(0.2))+')'+'\n')
bmd.write('arc 90 91 ('+str(c/2+mcyl*np.cos(alpha60-0.2))+' -0.5 '+str(-mcyl*np.sin(alpha60-0.2))+')'+'\n')
bmd.write('arc 91 46 ('+str(c/2+mcyl*np.cos(1.570796327-0.2))+' -0.5 '+str(-mcyl*np.sin(1.570796327-0.2))+')'+'\n')
bmd.write('arc 46 47 ('+str(c/2-mcyl*np.cos(0.7853981634))+' -0.5 '+str(-mcyl*np.sin(0.7853981634))+')'+'\n')




bmd.write('spline 0 1'+'\n')
bmd.write('	('+'\n')
counter = 0
for i in xu:
	if i<c/2:
		bmd.write('        ('+str(xu[counter])+' 1.4 '+str(zu[counter])+')'+'\n')
	counter += 1
bmd.write('	)'+'\n')
bmd.write('spline 1 92'+'\n')
bmd.write('	('+'\n')
counter = 0
for i in xu:
	if i>c/2 and i<c:
		bmd.write('        ('+str(xu[counter])+' 1.4 '+str(zu[counter])+')'+'\n')
	counter += 1
bmd.write('	)'+'\n')
bmd.write('spline 0 2'+'\n')
bmd.write('	('+'\n')
counter = 0
for i in xl:
	if i<c/2:
		bmd.write('        ('+str(xl[counter])+' 1.4 '+str(zl[counter])+')'+'\n')
	counter += 1
bmd.write('	)'+'\n')
bmd.write('spline 2 3'+'\n')
bmd.write('	('+'\n')
counter = 0
for i in xl:
	if i>c/2 and i<c:
		bmd.write('        ('+str(xl[counter])+' 1.4 '+str(zl[counter])+')'+'\n')
	counter += 1
bmd.write('	)'+'\n')
bmd.write('spline 4 5'+'\n')
bmd.write('	('+'\n')
counter = 0
for i in xu:
	if i<c/2:
		bmd.write('        ('+str(xu[counter])+' -0.500000 '+str(zu[counter])+')'+'\n')
	counter += 1
bmd.write('	)'+'\n')
bmd.write('spline 5 96'+'\n')
bmd.write('	('+'\n')
counter = 0
for i in xu:
	if i>c/2 and i<c:
		bmd.write('        ('+str(xu[counter])+' -0.500000 '+str(zu[counter])+')'+'\n')
	counter += 1
bmd.write('	)'+'\n')
bmd.write('spline 4 6'+'\n')
bmd.write('	('+'\n')
counter = 0
for i in xl:
	if i<c/2:
		bmd.write('        ('+str(xl[counter])+' -0.500000 '+str(zl[counter])+')'+'\n')
	counter += 1
bmd.write('	)'+'\n')
bmd.write('spline 6 7'+'\n')
bmd.write('	('+'\n')
counter = 0
for i in xl:
	if i>c/2 and i<c:
		bmd.write('        ('+str(xl[counter])+' -0.500000 '+str(zl[counter])+')'+'\n')
	counter += 1
bmd.write('	)'+'\n')
bmd.write('spline 82 78'+'\n')
bmd.write('	('+'\n')
counter = 0
for i in x2u:
	if i<c/2:
		bmd.write('        ('+str(x2u[counter])+' 1.4 '+str(z2u[counter])+')'+'\n')
	counter += 1
bmd.write('	)'+'\n')
bmd.write('spline 78 61'+'\n')
bmd.write('	('+'\n')
counter = 0
for i in x2u:
	if i>c/2 and i<c:
		bmd.write('        ('+str(x2u[counter])+' 1.4 '+str(z2u[counter])+')'+'\n')
	counter += 1
bmd.write('	)'+'\n')
bmd.write('spline 82 79'+'\n')
bmd.write('	('+'\n')
counter = 0
for i in x2l:
	if i<c/2:
		bmd.write('        ('+str(x2l[counter])+' 1.4 '+str(z2l[counter])+')'+'\n')
	counter += 1
bmd.write('	)'+'\n')
bmd.write('spline 79 62'+'\n')
bmd.write('	('+'\n')
counter = 0
for i in x2l:
	if i>c/2 and i<c:
		bmd.write('        ('+str(x2l[counter])+' 1.4 '+str(z2l[counter])+')'+'\n')
	counter += 1
bmd.write('	)'+'\n')
bmd.write('spline 83 80'+'\n')
bmd.write('	('+'\n')
counter = 0
for i in x2u:
	if i<c/2:
		bmd.write('        ('+str(x2u[counter])+' -0.500000 '+str(z2u[counter])+')'+'\n')
	counter += 1
bmd.write('	)'+'\n')
bmd.write('spline 80 70'+'\n')
bmd.write('	('+'\n')
counter = 0
for i in x2u:
	if i>c/2 and i<c:
		bmd.write('        ('+str(x2u[counter])+' -0.500000 '+str(z2u[counter])+')'+'\n')
	counter += 1
bmd.write('	)'+'\n')
bmd.write('spline 83 81'+'\n')
bmd.write('	('+'\n')
counter = 0
for i in x2l:
	if i<c/2:
		bmd.write('        ('+str(x2l[counter])+' -0.500000 '+str(z2l[counter])+')'+'\n')
	counter += 1
bmd.write('	)'+'\n')
bmd.write('spline 81 71'+'\n')
bmd.write('	('+'\n')
counter = 0
for i in x2l:
	if i>c/2 and i<c:
		bmd.write('        ('+str(x2l[counter])+' -0.500000 '+str(z2l[counter])+')'+'\n')
	counter += 1
bmd.write('	)'+'\n')
bmd.write(''+'\n')
bmd.write(');'+'\n')

bmd.write(''+'\n')
bmd.write('boundary'+'\n')
bmd.write('('+'\n')
bmd.write(''+'\n')
bmd.write('	shape'+'\n')
bmd.write('	{'+'\n')
bmd.write('		type wall;'+'\n')
bmd.write('		faces'+'\n')
bmd.write('		('+'\n')
bmd.write('		(4 0 1 5)'+'\n')
bmd.write('		(5 1 92 96)'+'\n')
bmd.write('		(7 3 92 96)'+'\n')
bmd.write('	        (0 4 6 2)'+'\n')
bmd.write('	        (2 6 7 3)'+'\n')
bmd.write('		);'+'\n')
bmd.write('	}'+'\n')
bmd.write(''+'\n')
bmd.write('	rotbox_slave'+'\n')
bmd.write('	{'+'\n')
bmd.write('		type patch;'+'\n')
bmd.write('		faces'+'\n')
bmd.write('		('+'\n')
bmd.write('	        (40 43 47 44)'+'\n')
bmd.write('	        (84 40 44 88)'+'\n')
bmd.write('	        (85 84 88 89)'+'\n')
# bmd.write('	        (41 85 89 45)'+'\n')
bmd.write('	        (41 95 99 45)'+'\n')
bmd.write('	        (95 85 89 99)'+'\n')

bmd.write('	        (86 41 45 90)'+'\n')
bmd.write('	        (87 86 90 91)'+'\n')
bmd.write('	        (42 87 91 46)'+'\n')
bmd.write('	        (43 42 46 47)'+'\n')
bmd.write('		);'+'\n')
bmd.write('	}'+'\n')
bmd.write(''+'\n')
bmd.write('	front'+'\n')
bmd.write('	{'+'\n')
bmd.write('		type patch;'+'\n')
bmd.write('		faces'+'\n')
bmd.write('		('+'\n')
bmd.write('		(43 11 8 40)'+'\n')
bmd.write('		(11 82 78 8)'+'\n')
bmd.write('		(82 0 1 78)'+'\n')
bmd.write('		(82 79 2 0)'+'\n')
bmd.write('		(11 10 79 82)'+'\n')
bmd.write('		(43 42 10 11)'+'\n')
bmd.write('		(1 92 61 78)'+'\n')
bmd.write('		(78 61 60 8)'+'\n')
bmd.write('		(8 60 84 40)'+'\n')
# bmd.write('		(3 66 65 61)'+'\n')
bmd.write('		(3 66 93 92)'+'\n')
bmd.write('		(92 93 65 61)'+'\n')

bmd.write('		(61 65 64 60)'+'\n')
# bmd.write('		(66 9 64 65)'+'\n')
bmd.write('		(66 9 94 93)'+'\n')
bmd.write('		(93 94 64 65)'+'\n')

bmd.write('		(60 64 85 84)'+'\n')
#bmd.write('		(9 41 85 64)'+'\n')
bmd.write('		(9 41 95 94)'+'\n')
bmd.write('		(94 95 85 64)'+'\n')

bmd.write('		(79 62 3 2)'+'\n')
bmd.write('		(10 63 62 79)'+'\n')
bmd.write('		(42 87 63 10)'+'\n')
bmd.write('		(87 86 68 63)'+'\n')
bmd.write('		(63 68 67 62)'+'\n')
bmd.write('		(62 67 66 3)'+'\n')
bmd.write('		(67 68 9 66)'+'\n')
bmd.write('		(68 86 41 9) '+'\n')
bmd.write('		);'+'\n')
bmd.write('	}'+'\n')
bmd.write(''+'\n')
bmd.write('	back'+'\n')
bmd.write('	{'+'\n')
bmd.write('		type patch;'+'\n')
bmd.write('		faces'+'\n')
bmd.write('		('+'\n')
bmd.write('		(15 47 44 12)'+'\n')
bmd.write('		(83 15 12 80)'+'\n')
bmd.write('		(4 83 80 5)'+'\n')
bmd.write('		(46 47 15 14)'+'\n')
bmd.write('		(15 14 81 83)'+'\n')
bmd.write('		(83 81 6 4)'+'\n')
bmd.write('		(12 69 88 44)'+'\n')
bmd.write('		(80 70 69 12)'+'\n')
bmd.write('		(5 96 70 80)'+'\n')
bmd.write('		(81 71 7 6)'+'\n')
bmd.write('		(14 72 71 81)'+'\n')
bmd.write('		(46 91 72 14)'+'\n')
bmd.write('		(69 73 89 88)'+'\n')
bmd.write('		(70 74 73 69)'+'\n')
# bmd.write('		(7 75 74 70)'+'\n')
bmd.write('		(7 75 97 96)'+'\n')
bmd.write('		(96 97 74 70)'+'\n')

bmd.write('		(71 76 75 7)'+'\n')
bmd.write('		(72 77 76 71)'+'\n')
bmd.write('		(91 90 77 72)'+'\n')
# bmd.write('		(75 13 73 74)'+'\n')
bmd.write('		(75 13 98 97)'+'\n')
bmd.write('		(97 98 73 74)'+'\n')

# bmd.write('		(13 45 89 73)'+'\n')
bmd.write('		(13 45 99 98)'+'\n')
bmd.write('		(98 99 89 73)'+'\n')

bmd.write('		(76 77 13 75)'+'\n')
bmd.write('		(77 90 45 13)'+'\n')
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
