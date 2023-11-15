# OFJ-pitching-airfoil

Two test cases are present :
- airfoil_OFJ_static
- airfoil_OFJ_pitching

-------------------
airfoil_OFJ_static:
-------------------
(U-RANS) airflow around a static airfoil oriented with a given AOA
---------------------------------------------------------------------

In the "mesh" folder, you will find:
  - nonRotating : background mesh
  - rotating    : mesh surrounding the airfoil
  - final       : the fusion of the two previous meshes


1. Setup cord ans AoA, in the file mesh/rotating/makeMesh :
  - translation of the airfoil so that half of the chord "c" is positioned in the axis of the marker [0,0].
    --> transformPoints -translate '(-c/2 0 0)'   (in this case c=0.6)

  - rotation of an angle of attack AOA with reference to the z axis (0 0 1).
    -->transformPoints -rotate '( (0 0 1) (sin(AOA) 0 cos(AOA)) )'   (in this case : AOA = 10 deg)

  - Finally, return to the initial position of the airfoil
    --> transformPoints -translate '(+c/2 0 0)'  (in this case c=0.6)

2. Type .Allrun


---------------------
airfoil_OFJ_pitching:
---------------------
(U-RANS) airflow around a pitching airfoil initially oriented at a given angle of attack AOA.
the airfoil pitches around the Y axis at the center (c/2 0 0) with c: chord.   (in this case c=0.6)
---------------------------------------------------------------------------------------------------------------
Computational Strategies: :
---------------------------
Strategy #1 :
-------------
- start directly the calculation with the initial conditions from "0.orig" : ./Allrun

  Notes: in this case, the simulation will take several days of calculation to reach the established statistics.

---------------------------
- Strategy #2:
--------------
  1 - cary out a calculation with airfoil_OFJ_static case (the flow needs to be established).
  2 - copy the results folder of the last iteration of the stationary case and rename it to "0" then paste it
      in the folder case_test_pitching/
  4 - Then, ./Allrun

  Notes: the computation time is considerably reduced compared to the strategy #1 to reach the established statistics.
