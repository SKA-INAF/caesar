=====================================
**Frequently Asked Questions (FAQs)**
=====================================

Below some of the questions you might have in mind when using CAESAR:

**1. What operating systems are supported?**   

CAESAR is specifically developed for Linux OS. It was used on Ubuntu and RedHat distributions. In principle it should be usable also in MacOS as all external dependencies (ROOT, OpenCV, boost, ...) are supported in MacOS. For this you will need to manually modify the CAESAR Makefile. No support for Windows OS is available and never will be. If you are in trouble with the installation and you want to get started in short time, we provide a Singularity container with all software and dependencies already installed, ready to be used. Please read the `Container` documentation section for further details.   

**2. What about python support?**   

I know astronomers love python...CAESAR is developed in C++ and for this version you need to use it as it is. I used C++ for both personal preference, to re-use existing code and for performance reason. It is however possible to use CAESAR classes in python through the PyROOT interface. This feature is at present experimental and to be tested. Please read the `Usage` documentation section for further details.
