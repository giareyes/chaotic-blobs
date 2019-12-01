Name: Gillian (Gia) Reyes
netid: gpr9
HW #2

1. (10 points) Compute the correct F
Everything seems to work fine with the F I computed, so I hope its right

2. (10 points) Implement the correct StVK strain energy psi
I believe this is correct because my finite difference method using psi
works. See the function PK1Difference() which uses it in STVK.cpp

3. (10 points) At startup, verify that your StVK PK1 is correct
I have shown that this work every time you launch my code. If you
look at PK1Difference() in STVK.cpp, you will see that this function
tests a random F force gradient and shows that the derivative of psi
by finite difference method converges to PK1. This function is called
every time the program is launched using ./HW2

4. (10 points) At startup, verify that your StVK Hessian is correct
I have shown that this work every time you launch my code. If you
look at HessianDifference() in STVK.cpp, you will see that this function
tests a random F force gradient and shows that the derivative of PK1
by finite difference method converges to DPDF, aka the hessian.
This function is called every time the program is launched using ./HW2

5. (20 points) Form the K matrix correctly for the quasistatic solver.
For the single triangle case, print this matrix to the command line.
Pretty sure this works properly

6. (20 points) Form the material force vector correctly in the quasistatics solver.
For the single triangle case, print this vector to the command line.
Also pretty sure this works properly

7. (10 points) Verify that the simulation works under the SINGLE, HANG, SHEAR, SQUASH,
and STRETCH scenarios. In the SQUASH scenario, the mesh should initially compress, but
then fail to recover from inversion.
Verified that it looks like its working!

8. (10 points) Follow the formatting instructions correctly.
I think I did this?
