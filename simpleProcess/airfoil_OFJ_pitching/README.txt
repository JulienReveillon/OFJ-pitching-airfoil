To execute the case  :

1/ make sure that the zipping process did not change the scripts +x execution mode or type :

chmod +x Allclean Allrun makeMesh mesh/rotating/makeMesh mesh/mergeMeshes

2/ edit the file system/decomposeParDict and indicate the number of available processors (numberOfSubdomains)

3/ execute the Allrun script: ./Allrun
