=================
**Installation**
=================

-----------------
**Prerequisites**
-----------------

Install the project mandatory dependencies:  

* **ROOT** [https://root.cern.ch/], to be built with FITSIO, PyROOT, RInterface options enabled. Make sure that the FindROOT.cmake is present in $ROOTSYS/etc/cmake directory after installation.
* **OpenCV** [http://opencv.org/]
* **log4cxx** [https://logging.apache.org/log4cxx/]
* **boost** [http://www.boost.org/] 
* **cfitsio** [https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html], to be built with multithread support (e.g. give --enable-reentrant option in configure)
* **protobuf** [https://github.com/google/protobuf]
* **jsoncpp** [https://github.com/open-source-parsers/jsoncpp]
* **python** (>=2.7) [https://www.python.org/], install also these additional modules: pyfits, astropy
* **cmake** (>=2.8) [https://cmake.org]  
  
Optional dependencies are:

* **MPICH** [https://www.mpich.org/] or **OpenMPI** [https://www.open-mpi.org/], needed when the build option ENABLE_MPI=ON    
* **OpenMP** [http://www.openmp.org/], needed when the build option BUILD_WITH_OPENMP=ON
* **R** [https://www.r-project.org/] and additional packages: RInside, Rcpp, rrcovHD, truncnorm, FNN, akima. Needed when the build option ENABLE_R=ON
* **GoogleTest** [https://github.com/google/googletest], needed for unit testing when the build option ENABLE_TEST=ON
* **Doxygen** [www.doxygen.org/], needed to generate the API documentation
* **Sphinx** [http://www.sphinx-doc.org] & **Breathe** [https://pypi.org/project/breathe], needed to generate the Sphinx API & wiki documentation

Make sure you have set the following environment variables to the external library installation dirs:

* ROOTSYS
* OPENCV_DIR
* R_DIR
* BOOST_ROOT
* LOG4CXX_ROOT
* JSONCPP_ROOT

Add also the following paths to the PKG_CONFIG_PATH environment var: 

* $LOG4CXX_ROOT/lib/pkgconfig
* $JSONCPP_ROOT/lib/pkgconfig

CAESAR depends also on the wcstools and linterp libraries which are already provided in the external/ directory. 

.. warning::

   The provided wcslib was slightly modified with respect to the original release to avoid naming conflicts with the R package due to some #define macros used in WCS.

cmake should find all needed include dirs and libraries used to build the project.

---------------------
**Build and install**
---------------------

To build and install the project:

* Clone this repository into your local $SOURCE_DIR:  

	``git clone https://github.com/SKA-INAF/caesar.git $SOURCE_DIR``

* Create the build and install directories: $BUILD_DIR, $INSTALL_DIR  

* Configure build. In the build directory type:

	``cmake -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR -DBUILD_WITH_OPENMP=ON -DENABLE_MPI=ON \``
	         ``-DBUILD_APPS=ON -DENABLE_TEST=ON $SOURCE_DIR``

* Building project. In the build directory type:

	``make``

* Install project in the $INSTALL_DIR directory specified in the configuration step above. in the build directory type:

	``make install``

In the installation directory you should find the following directories:

* `include`: Directory in which project source headers are installed.
* `bin`: Directory in which app executables are installed.
* `lib`: Directory in which Caesar shared library is installed.
* `share`: Directory in which documentation and other misc files are installed.
* `script`: Directory in which python & shell scripts (simulation, submission, etc) are installed.
* `macros`: Directory in which util macros (e.g. ROOT macros, atc) are installed.
* `data`: Directory in which test data images are installed.


----------------------------
**Documentation generation**
----------------------------

To generate and install the API documentation you must have Doxygen installed. Enter the build directory and type:


``make doc``


To generate and install the Sphinx wiki documentation you must have sphinx+breathe installed. Enter the build directory and type:


``make doc-sphinx``


The generated documentation will be installed in the $INSTALL_DIR/doc directory.


----------------------
**Running unit tests**
----------------------

To build the unit tests you must have Google Test installed and the ENABLE_TEST option set to ON when building Caesar. To run the test:   

``make test``    

or alternatively run the script `runUnitTests` installed in the Caesar installation dir.

