=================================
Frequently Asked Questions (FAQs)
=================================

Below some of the questions you might have in mind when using CAESAR:

**1. In which language is implemented CAESAR?**   

CAESAR is implemented in C++.


**2. Does CAESAR support python?**   

I know astronomers love python...CAESAR is developed in C++ and for this version you need to use it as it is. I used C++ for both personal preference, to re-use existing code and for performance reason. It is however possible to use CAESAR classes in python through the PyROOT interface. This feature is at present experimental and to be tested. Please read the `Usage` documentation section for further details.

**3. In which operating systems CAESAR is supported?**   

CAESAR is specifically developed for Linux OS. It was used on Ubuntu and RedHat/Centos distributions. In principle it should be usable also in MacOS as all external dependencies (ROOT, OpenCV, boost, ...) are supported in MacOS. For this you will need to manually modify the CAESAR Makefile. No support for Windows OS is available and never will be. If you are in trouble with the installation and you want to get started in short time, we provide a Singularity container with all software and dependencies already installed, ready to be used. Please read the `Container` documentation section for further details.   

**4. Does CAESAR support distributed or parallel processing?**   

Yes, CAESAR source finding supports two level of parallelism since 2018. Images can be partitioned in adjacent/overlapping tiles and source findind can be carried out on different computing processors using MPI library. Source finding tasks per each tile can be splitted among different threads using OpenMP library. Alternative programming models will be explored
in the future.

**5. What data input formats can be handled?**

CAESAR currently works on 2D images in different formats:

- FITS
- ROOT native Caesar::Image
- Standard image formats (png, jpeg, etc)
- OpenCV Mat
- VTK (to be tested)

Other formats are planned to be added:

- CASA images
- HDF5

Cubes are not supported.

**6. What data outputs are delivered to user?**

CAESAR returns the following outputs to user at the end of processing:

- A ROOT file containing the following information:
    - Computed maps (bkg, rms, significance, residual, etc) stored as Caesar::Image objects
    - Detected sources stored as a ROOT TTree of Caesar::Source objects
    - Run configuration options stored as a ROOT TTree

- Two DS9 file regions containing:
    - Detected sources (also called islands in other finders)
    - Gaussian components fitted to detected sources
 
- Two ascii catalog files containing parameters relative to:
    - Detected sources (also called islands in other finders)
    - Gaussian components fitted to detected sources

- Computed maps (bkg, rms, significance, residual, etc) in FITS format



