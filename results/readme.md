# Results obtained with the help of the Runaphys library

Filenames contain:
- the type of generation
- the resolution
- "withI" denotes that the calculations with the REs in magnetic islands were done separately

# A few select images
Dreicer generation with 0 RE density initial condition, avalanche generation with uniform density initial condition.

Dreicer generation 1000x1000 resolution
![Dreicer](https://github.com/leferi99/szakdolgozat/blob/main/results/Dreicer_1000x1000.jpg)

Dreicer generation 1000x1000 resolution with islands treated differently (probably not correct!)
![Dreicer_withI](https://github.com/leferi99/szakdolgozat/blob/main/results/Dreicer_withI_1000x1000_v2.jpg)

Avalanche generation 1000x1000 resolution
![avalanche](https://github.com/leferi99/szakdolgozat/blob/main/results/avalanche_1000x1000.jpg)

Avalanche generation 1000x1000 resolution with islands treated differently - This is the closest to the reference material, although that code did not use generation. The generation is almost non-existent because the collision time is around 0.02s for the used parameters, and I only examine 0.001s in these plots.
![avalanche_withI](https://github.com/leferi99/szakdolgozat/blob/main/results/avalanche_withI_1000x1000.jpg)

0.1s examined time with 0.001 times the original transport coefficients, just to see the avalanche generation.
![avalanche_0.1s](https://github.com/leferi99/szakdolgozat/blob/main/results/avalanche%200.1%20sec%2C%20transport_coeffs%200.001%20times%20original.jpg)
