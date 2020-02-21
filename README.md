# SimRabbit

This repo includes some basic Rabbit model simulation.

Run 'sim_rabbit_gait.m' to get everything started.

I've included a 'data' folder that has a library of FROST gaits, 
as well as a function that loads the FROST gaits and gets them into a usable format.

The function 'animateRabbit' draws links, and you can use it to figure out how to go from local to global locations.

Exploring the various parts of 'sim_rabbit_gait.m' and functions called there should
help get a handle on the Rabbit simulation framework.

The folders "SRC", "rabbit_sim_files", and "autogen" contain various functions for simulation that I didn't feel like sorting
So, a lot of them may be useless for this application.
