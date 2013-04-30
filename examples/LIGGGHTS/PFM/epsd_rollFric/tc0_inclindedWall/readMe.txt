INCLINED WALL LIGGGHTS INPUT SCRIPT

After measuring the mean angle, where the particles start to roll, you may use the in_inclinedWall.liggghts script to determine the rolling friction coefficient of the particles.

The first variable "targetAngle" is the angle, measured at the experiment in degrees.

The script consists of three parts:

1) particle insertion and relaxation
	100 particles are inserted onto an even plate and the script runs until the kinetic energy is lower than 10^-7 J.

2) rotation of the gravity vector to targetAngle
	The gravitation vector is rotated to get an inclined wall.

3) iteration of rolling friction
	The rolling friction coefficient is decreased from 1.5 by 0.5 until the particles start to roll.
	Then velocities will be set to zero and the last not rolling coefficient is reloaded.
	The step size is divided by 5 and the loop starts again.
	The check if the particles are rolling, is done by calculating the mean y-velocity.
	This value is compared with "${vyThresh}".


After the last refinement loop, the coefficient when the particles last started rolling is displayed as "final rf".
