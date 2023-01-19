Black hole renderer (BHR)
=========================

BHR is a C++ code for visualizing black holes and neutron stars with simple accretion disks.

Implemented features:

- integration of geodesics (ray tracing) in flat and curved spacetime (rotating and non-rotating black holes) [[1]](#cite_spacetimes)
- computation of polarization of light along the geodesic [[2]](#cite_polarization)
- relativistic Doppler effect [[2]](#cite_polarization)
- bending of light caused by nonlinear interaction with very strong magnetic fields of neutron stars [[3]](#cite_magnetic)
- Shakura-Sunyaev accretion disk (outer layer) [[4, 5, 6]](#cite_disk)

Other aesthetic features:

- "sky" background texture
- quadtree-based antialiasing and optimization


Compilation
-----------

Use the following commands to compile the code:

    mkdir -p build
    cd build
    cmake ..
    make


Usage
-----

**Note:** The BHR code is quite experimental and not entirely user-friendly at the moment.
The settings like the type of the spacetime, the spacetime parameters, the disk type and what quantities are visualized are all currently set at compile time (see [](src/bhr/config.h), [](src/bhr/parameters.h) and [](src/bhr/main.cpp)).
This is because the integration of geodesics is quite slow, so the idea was to disable everything unused and help the compiler optimize the code as much as possible.
Maybe one day I refactor the code and make everything configurable at runtime.
The output image resolution, camera position and similar parameters are currently configurable through command-line arguments (see [](src/bhr/settings.h)).

To run the code, execute the following:

    # Go back from build/ to the repository root.
    cd ..

    # Prepare output folders.
    mkdir -p output/preprocess
    mkdir -p output/raw

    # Run.
    ./build/bhr --width 1024 --aspect 16:9 --output_image output.tga --no-cache

By default, the code renders a rotating black hole of mass 1.25 M<sub>Sub</sub> with an angular moment of 0.998 (relative).
(The mass is so small because it actually represents a neutron star, but its radius is set to 0.)

The expected output is the following, with a different color scheme and no legend:

![](output_annotated.png "Output of bhr + legend.")

The black lines denote the polarization, and the colors the light intensity.
Due to the Doppler effect, the left side (rotating towards the camera) is a few orders of magnitude brighter than the right side.


Gallery
-------

More renders are available [here](https://drive.google.com/drive/folders/0B1mAEaKMwKIVMklmUGx0VUJwNWs?resourcekey=0-gK3q6emMQNOAdDPQ5rEK3A&usp=sharing).


Lite version
------------

A smaller single-file code for rendering a black hole with a disk is available in `src/lite.cpp`.
The code supports only the Kerr metric (rotating black hole) and dummy disk coloring.

Use the following commands to compile and run it:

    mkdir -p build
    cd build
    cmake ..
    make lite
    ./lite

Alternatively, compile manually using:

    cd src
    g++ -O3 -march=native -std=c++17 lite.cpp -o lite

The lite code produces the following image (stored as ``output_lite.tga``).

![](output_lite.png "Output of lite.cpp")

The black hole is marked in red.
The colors on the disk denote the location on the disk.
Multiple images of the same disk can be seen.


References
----------

1. <a id="cite_spacetimes"></a>
   Müller, T., & Grave, F.,
   *Catalogue of spacetimes.*
   arXiv preprint (2009)
   [arXiv:0904.4184](https://arxiv.org/abs/0904.4184v3)

2. <a id="cite_polarization"></a>
  Chen, B., Kantowski, R., Dai, X., Baron, E., & Maddumage, P.,
  *Algorithms and programs for strong gravitational lensing in Kerr space-time including polarization.*
  The Astrophysical Journal Supplement Series (2015)
  [10.1088/0067-0049/218/1/4](https://doi.org/10.1088/0067-0049/218/1/4)
  [arXiv:1505.02714](https://arxiv.org/abs/1505.02714)

3. <a id="cite_magnetic"></a>
  Kim, J. Y., & Lee, T.,
  *Light bending by nonlinear electrodynamics under strong electric and magnetic field.*
  Journal of Cosmology and Astroparticle Physics (2011)
  [10.1088/1475-7516/2011/11/017](https://doi.org/10.1088/1475-7516/2011/11/017)
  [arXiv:1101.3433](https://arxiv.org/abs/1101.3433)

4. <a id="cite_disk"></a>
  Shakura, N. I., & Sunyaev, R. A.,
  *Black holes in binary systems.*
  Observational appearance. Astronomy and Astrophysics (1973)
  [[pdf]](https://adsabs.harvard.edu/pdf/1973A&A....24..337S)

5. <a id="cite_disk2"></a>
  Novikov, I. D., & Thorne, K. S.,
  *Astrophysics of black holes.*
  Black holes. (1973)

6. <a id="cite_disk3"></a>
  Page, D. N., & Thorne, K. S.,
  *Disk-accretion onto a black hole. Time-averaged structure of accretion disk.*
  The Astrophysical Journal (1974)
