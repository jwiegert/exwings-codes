# exwings-codes

Source codes for (primarily) translating time-dependent data from radiation-hydrodynamical simulations of evolved stars and surrounding dust from the code CO5BOLD to the Monte-Carlo based radiative transfer simulating code RADMC-3D. This can be used to create synthetic observations.


## Running the translations

Note: Most links are hard-coded. Assumes access to CO5BOLD-data in IDL-formatted .sav-files.

- ipynb-files contains various tests, and code snippets to run the various functions contained in the .py-files.

- To run the translation automatically use either

> scriptbash_creategridstar.sh

to create data containing only the gas portion, i.e., the star. Or

> scriptbash_creategridstardust.sh

to create RADMC-3D data with both gas (star), dust, and dust opacities. Assumes availability of CO5BOLD-data, the code OPTOOL (to create opacity files), and laboratory dust opacity data from the Jena Optical Database.


## Links

https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/index.php

https://github.com/cdominik/optool

https://www.astro.uni-jena.de/Laboratory/Database/databases.html

