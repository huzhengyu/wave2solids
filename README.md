olaFluid
======
# Description
[olaFuid](https://github.com/huzhengyu/olaFluid) is a numerical wave-structure-interaction (WSI) model combining the wave generation toolbox [olaFlow](https://github.com/phicau/olaFlow) and the fluid-strucutre-interaction (FSI) toolbox [solids4Foam](https://bitbucket.org/philip_cardiff/solids4foam-release/src/master/) in the framework of OpenFOAM.

This free and open-source project is committed to enabling the fully-coupled simulation of wave interactions with flexible structures to the foam-extend communities.

The paper for implementing and applying this model has been submitted to Coastal Engineering.

# Download and compilation
## Basic download guide
To get set up, you must first install [foam-extend-4.0](https://github.com/Unofficial-Extend-Project-Mirror/foam-extend-foam-extend-4.0).

To get a copy of [olaFlow](https://github.com/phicau/olaFlow) source and reference materials, run in a terminal:

`git clone https://github.com/phicau/olaFlow.git`

To get a copy of [solids4Foam](https://bitbucket.org/philip_cardiff/solids4foam-release/src/master/) source and reference materials, run in a terminal:

`git clone https://bitbucket.org/philip_cardiff/solids4foam-release.git`

Go to the [solids4Foam](https://bitbucket.org/philip_cardiff/solids4foam-release/src/master/) base folder and change to fluidModels directory:

`cd solids4foam-release/src/solids4FoamModels/fluidModels/`

To get a copy of the combined fluid model [olaFuid](https://github.com/huzhengyu/olaFluid), run in a terminal:

`git clone https://github.com/huzhengyu/olaFluid.git`

Make some modifications in the compilation files, run in a terminal:

`sed -i '35i fluidModels/olaFluid/olaFluid.C' ../Make/files.foamextend`

`sed -i '$s/$/& \\\n\t-lwaveGeneration \\\n\t-lwaveAbsorption/' ../Make/options ../../../applications/solvers/solids4Foam/Make/options`

## Basic compilation guide
The compilation is straightforward, simply run the following script from the [olaFlow](https://olaflow.github.io) base folder:

`./allMake`

Also, simply run the following script from the [solids4Foam](https://bitbucket.org/philip_cardiff/solids4foam-release/src/master/) base folder:

`./Allwmake`


# Tutorials
A tutorial case of a flexible wall in nonlinear periodic waves is provided.

Go to the Tutorials/Wall folder, and run in a terminal:

`./Allrun`

Feel free to modify the phrase and adapt it to your own needs.

# References
Higuera, P., Lara, J.L. and Losada, I.J., 2013. Realistic wave generation and active wave absorption for navier–stokes models: Application to openfoam®. Coastal Engineering, 71: 102-118.

Cardiff, P. et al., 2018. An open-source finite volume toolbox for solid mechanics and fluid-solid interaction simulations. arXiv preprint arXiv:1808.10736.

Tuković, Ž., Karač, A., Cardiff, P., Jasak, H. and Ivanković, A., 2018. Openfoam finite volume solver for fluid-solid interaction. Transactions of FAMENA, 42(3): 1-31.
