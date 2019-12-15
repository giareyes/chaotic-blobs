# chaotic-blobs

This project is an implementation of Barbic's 2005 paper Real-Time Subspace Integration for St.Venant-Kirchho Deformable Models, with the addition of motion.

To see the results of just Barbic's paper:

Step 1:go into TRIANGLE_MESH.cpp.
Step 2: Comment out the checkCollisions call in stepMotion.
Step 3: Also comment out

intermediate.col(0) = T.col(0);
intermediate.col(1) = T.col(1);
for(int i = 2; i < 11; i++)
{
  intermediate.col(i) = svddiag.col(i - 2) - (T.col(0).transpose()\*svddiag.col(i - 2))\*T.col(0);
  intermediate.col(i) = intermediate.col(i) - (T.col(1).transpose()\*intermediate.col(i))\*T.col(1);
}

\_U = intermediate;

in the function setBasisReduction, and uncomment

\_U = svddiag;


Step 4: go into blobs.cpp
Step 5: in GlutIdle, comment out

if (frame == 0) // if we are on the first frame, insert a force so they jump
{
  triangleMesh.addBodyForce(VEC2(10.0, 200.0));
  blob2.addBodyForce(VEC2(-10.0, 200.0));
}

and

TriangleMesh.addBodyForce(bodyForce);
blob2.addBodyForce(bodyForce);
