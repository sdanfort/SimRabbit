# SimRabbit

This repo provides an exmaple for simulating several gaits of a five-link Rabbit model (see https://web.eecs.umich.edu/~grizzle/CDC_2004_Running/Chevallereau_CDC_2004_submission.pdf for model specifics).

Run 'sim_rabbit_gait.m' to get everything started.

I've included a "data" folder that has a library of walking gaits found with FROST (https://ayonga.github.io/frost-dev/), 
as well as a function that loads the FROST gaits and puts them into a usable format.

The function 'animateRabbit' draws the links.

The 'sim_rabbit_gait.m' function should help users understand the Rabbit simulation framework.

The folders "SRC", "rabbit_sim_files", and "autogen" may contain functions that are not useful for this application.
