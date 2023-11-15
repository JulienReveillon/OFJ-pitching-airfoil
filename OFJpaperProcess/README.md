# OFJ-pitching-airfoil

Three test cases are available :

- airfoil_OFJ_static:

    * static_RANS (simpleFoam)
    * static_URANS (pimpleFoam)

- airfoil_OFJ_pitching (pimpleFoam)

-------------------
airfoil_OFJ_static:
-------------------
(static_RANS) airflow around a static airfoil oriented with a given AOA
---------------------------------------------------------------------

In the "mesh" folder, you will find:
  - nonRotating : background mesh
  - rotating    : mesh surrounding the airfoil
  - final       : the fusion of the two previous meshes

1. Setup cord and AoA, in the file mesh/rotating/makeMesh :
  - translation of the airfoil so that half of the chord "c" is positioned in the axis of the marker [0,0].
    --> transformPoints -translate '(-c/2 0 0)'   (in this case c=0.6)

  - rotation of an angle of attack AOA with reference to the z axis (0 0 1).
    -->transformPoints -rotate '( (0 0 1) (sin(AOA) 0 cos(AOA)) )'   (in this case : AOA = 6 deg)

  - Finally, return to the initial position of the airfoil
    --> transformPoints -translate '(+c/2 0 0)'  (in this case c=0.6)

2. run simulation :  ./Allrun




---------------------------------------------------------------------
(static_URANS) airflow around a static airfoil oriented with a given AOA
---------------------------------------------------------------------
Computational Strategy: :
---------------------------
1. First, you need to run a static-RANS calculation. 

2. Once the calculation has converged, 
   - we copy the directory from the last iteration (exemple: "2569") and paste it into the "static_URANS" directory, renaming it "0". so that it can be taken as the initial condition. 

to run simulation: ./Allrun (everything is automated)

NOTE!!!: This strategy reduces computation time significantly. 



---------------------
airfoil_OFJ_pitching:
---------------------------------------------------------------------------------------------------------
airflow around a pitching airfoil initially oriented at a given angle of attack AOA.
the airfoil pitches around the Y axis at the center (c/2 0 0) with c: chord.   (in this case c=0.6)
---------------------------------------------------------------------------------------------------------
Computational Strategy: :
-------------------------

1. First, you need to run a static-URANS calculation. 

2. Once the calculation has finished, we copy the directory from the last time (exemple: "1") and paste it into the airfoil_OFJ_pitching directory, renaming it "0". so that it can be taken as the initial condition. 

to run simulation: ./Allrun (everything is automated)

NOTE!!!: This strategy reduces computation time significantly. 
