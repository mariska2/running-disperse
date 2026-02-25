# Running DisPerSE: A Tutorial and Functions for Running the Persistent Structures Identifier
### Introduction
The Discrete Persistent Structures Extractor, developed by Theirry Sousbie (2012), is an industry standard cosmic web extractor for astrophysics and cosmological research. It identifies persistent topological features, such as voids, filaments, and peaks from a dataset. In astronomy, it is used to determine the large scale structure of the universe, tracing out the filamentary structure of galaxies. 

While it is an essential program, there is a learning curve for getting it started and running, along with understanding and interpreting its outputs. This notebook (a product of my own research using DisPerSE) provides an in-depth tutorial on setting up DisPerSE (including laying out everything required for the program to function), running the program, interpreting and visualizing outputs, and performing scientific analysis with the results. 

Note: the setup of this notebook and the provided functions surround using DisPerSE in the context of galaxies and cosmological structure, and may need adjustment if working with data that is not from observational or simulated astronomy. 

### How to Use
The entire tutorial is contained within the Jupyter notebook in this repository. You can either download the files and ensure they are all contained within the same folder, or fork the repository and clone it from there. The python file, "disperse_function.py", contains the functions used within the notebook, so make sure to run the first cell of the notebook to import it. The .csv and .txt files are example datasets used to show the analysis, but I encourage changing out the data for your own. Below are markdown cells from **within the notebook**, for reference.

### Downloading DisPerSE Files & Necessary Programs
To get DisPerSe running, go to https://github.com/thierry-sousbie/DisPerSE and download the files. Older versions had a few errors that needed to be fixed, but the newest version available here is good to go.

If on a Windows device, a Linux environment is needed. I recommend WSL Ubuntu, which can be installed in Powershell with wsl --install, and started with wsl ~. After this, it will use regular bash syntax. You will have to create a venv, reinstall any packages/programs since it is a new environment, and restart your instance of this notebook through it (i.e. run jupyter lab after switching to WSL and then load the file, switching the kernel to your created venv).

Ensure that CMake is installed on your device, if not it can be downloaded here: https://cmake.org/download/. You will also need to have a C/C++ compiler downloaded somewhere on your device for CMake to work, I used Visual Studio build tools and made sure "C++ Cmake" was installed.

### Running Data through DisPerSE
The important part of running your data through the program is to have everything in a properly formatted tab separated text file. For galaxies, such as the ones loaded below from the SDSS catalog, the coordinates are in right ascension, declination, and redshift, but your points can technically be in any unit as long as they are all consistent with each another. DisPerSE doesn't like headers or IDs included in 3D fields, as (in my experience) it will throw an error saying it's a 2D field instead.

There are three major functions in DisPerSE for analyzing galaxy data:

- delaunay_3D; creates a 3D delaunay tesselation based off an inputted field (ra, dec, and z in .txt file), output of a .NDfield
- mse; generates the morse-smale complex from the .NDfield, outputs another one
- skeleconv; takes the .NDfield from binary back to the readable format

which is the outline I will be using in this notebook.

### Credits
Credit to Theirry Sousbie for the development of DisPerSE and Nicholas Luber for guidance on this project
