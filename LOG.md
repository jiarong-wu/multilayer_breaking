2024/12/12:

We corrected the `vorticity.c` file. 
Still need to make the same changes to `vorticity.h`, and then push to github, and then tell Rui about it. And it's weird that using central difference or backward difference changes the value of omegax and omegay in a large portion of the domain. Also the sign of omegay still seems incorrect. Right now I'm computing inner layers with central difference and top layer with backward difference for best comparison with the post-processing approach. 

Now the vorticity and more or less the same as to the one computed in post-processing. We can just keep using the post-processing code without recomputing vorticity in Basilisk. There is still some "barotropic" difference but seems ok after averaging.

\<sijsij> is still different from \<omegaiomegai>. Need to compute the boundary contribution term given by Pavel.