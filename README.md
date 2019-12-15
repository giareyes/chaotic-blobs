# chaotic-blobs

This project is an implementation of Barbic's 2005 paper Real-Time Subspace Integration for St.Venant-Kirchhoff Deformable Models, with the addition of motion.

## Barbic Paper Instructions
To see the results of just Barbic's paper, run:

 ./blob MOTION -b

To see the results of Barbic's paper with the mesh drawn, run:

 ./blob MOTION -m -b

## Motion Instructions 
 To run motion results, you can run either:

 ./blob

 or

 ./blob MOTION -m

 to see the mesh drawn in.

## Quasistatics Instructions 
 To run quasistatics on the blobs, run:

 ./blob [SINGLE, LSHEAR, RSHEAR, SQUASH, STRETCH] [-m]

 choose one test from the first list, and optionally show the mesh. 
