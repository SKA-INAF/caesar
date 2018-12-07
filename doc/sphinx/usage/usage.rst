=================
How to use CAESAR
=================

I managed to build the project with success or I got a copy of CAESAR container. Now what should I do?   
CAESAR library is built upon the ROOT framework and can be used in Linux OS (in principle usable also in MacOS but not tested) in different ways for simple tasks as well to build more complex applications:   

* Interactively from the ROOT CLI
* In C++ macros run by the ROOT interpreter
* In C++ standalone programs
* Interactively from the python/ipython CLI
* In python scripts run by the python interpreter


If you have worked with ROOT you will find these guidelines familiar.


------------------------------
Using CAESAR from the ROOT CLI
------------------------------

CAESAR generates dictionaries for all relevant classes, including the `Image` class. This enables to use the CAESAR classes in the ROOT CLI. For this you should put these lines in your ``.rootlogon.C`` file.   

.. code::

    gSystem->Load("[path to CAESAR library libCaesar.so]");
    gROOT->ProcessLine(".include [path to CAESAR include dir]");
    gInterpreter->AddIncludePath("[path to CAESAR include dir]");
    using namespace Caesar;

Additionally add CAESAR library (``libCaesar.so``) and path to dictionary file (``CaesarDict_rdict.pcm``) to your LD_LIBRARY_PATH environment variable. For example for a Bash shell:
  
.. code:: bash

    export LD_LIBRARY= $LD_LIBRARY_PATH:[path-to-CAESAR-lib]:[path-to-CAESAR-dict-file]

The `.rootlogon.C` is loaded every time the ROOT console is started from Linux prompt. To start ROOT CLI type ``root`` at the Linux shell prompt. At this point you should be ready to use CAESAR inside ROOT. The following example shows how to create a Caesar image from a FITS file (say it is named `map.fits`), compute statistics, background and noise maps:

.. code:: bash
    
    Image* img= new Image
    img->ReadFITS("map.fits")
    img->ComputeStats(true)
    ImgBkgData* bkgData= img->ComputeBkg(eMedianBkg,true,100,100,10,10)

This is useful for simple tasks or for drawing purposes. For more complex tasks, you should actually write macros or higher level classes/applications using CAESAR API, like it is briefly discussed below and done in the examples reported in the Tutorial section.

---------------------------
Using CAESAR in ROOT macros
---------------------------

ROOT macros are simply C/C++ code run by the ROOT interpreter. Once you have followed the configuration steps illustrated in the previous paragraph you should be able to run a macro using CAESAR objects in ROOT. Let's prepare a simple macro (file named ``MyMacro.C``) below:   
    
.. code-block:: cpp
    
    //== MyMacro.C ==
    #include <Image.h> //needed when compiling the macro
     
    void MyMacro()
    {     
      Image* img= new Image;//create a new empty image       
      img->ReadFITS("map.fits");//fill image from fits   
      img->ComputeStats(true);//compute standard & robust stats
      ImgBkgData* bkgData= img->ComputeBkg(eMedianBkg,true,100,100,10,10);//compute bkg data	
		
      Image* rmsMap= bkgData->NoiseMap;//get compute rms map
      rmsMap->GetHisto2D("histo")->Draw("COLZ");//draw rms map as histo 2d

      // etc etc...     
    }

To execute this macro in ROOT you can do the following:

* Run from the Linux shell prompt:
  
  .. code:: bash

      root MyMacro.C


* Run from the ROOT prompt:

  .. code::

      .x MyMacro.C


* Compile the macro and run from the ROOT prompt:

  .. code::
    
      .L MyMacro.C+
      MyMacro()


We refer the reader to the ROOT manual for more details on running macros, passing arguments to them, etc.


--------------------------------
Using CAESAR in C++ applications
--------------------------------

To use CAESAR library in your C++ application you just need to add the CAESAR headers (in the _include_ directory) in your compilation and link against CAESAR library (_libCaesar.so_). Let's prepare a simple C++ program (file named ``MyApp.cc``:   

.. code-block:: cpp

    #include <Image.h>
    #include <TApplication.h>//if you need to draw in ROOT canvas

    int main(int argc, char **argv) 
    {      
      //Needed if you want interactivity and draw in ROOT (not needed for batch applications)    
      TApplication* app= new TApplication("Application", 0, 0);
       
      Image* img= new Image;//create a new empty image       
      img->ReadFITS("map.fits");//fill image from fits   
      img->ComputeStats(true);//compute standard & robust stats
      ImgBkgData* bkgData= img->ComputeBkg(eMedianBkg,true,100,100,10,10);//compute bkg data	

      Image* rmsMap= bkgData->NoiseMap;//get compute rms map
      rmsMap->GetHisto2D("histo")->Draw("COLZ");//draw rms map as histo 2d
      
      //This will draw the image and suspend execution (remove in batch apps)
      app->Run();
       
      return 0;
    }

Now compile and execute the program:

.. code:: bash

    g++ -std=c++11 -g -o MyApp MyApp.cc \
        -I[path-to-CAESAR-headers] -I`root-config --incdir` \ 
        `pkg-config $OPENCV_DIR/lib/pkgconfig/opencv.pc --cflags` \
        -L[path-to-CAESAR-lib-dir] -lCaesar `root-config --libs` \
        `pkg-config $OPENCV_DIR/lib/pkgconfig/opencv.pc --libs`

    ./MyApp


----------------------------------
Using CAESAR in python CLI/scripts
----------------------------------

PyROOT interface enables using Caesar classes in python. For example:

.. code:: python

		from ROOT import gSystem
		gSystem.Load('libCaesar')
		from ROOT import Caesar

		img= Caesar.Image()
		img.ReadFITS('recmap.fits')  
		img.ComputeStats(True)    
		bkgData= img.ComputeBkg(Caesar.eMedianBkg,True,100,100,10,10)

python support is currently experimental and not fully tested. 

--------------------
Running CAESAR tasks
--------------------

Long-running tasks and applications such as source finding should be run in batch mode directly on the operating system or inside a Singularity container. A number of applications are available in the `bin` installation directory. Task run can be customized via a configuration file passed as argument. For example to run source finding:

.. code:: bash

    FindSourceMPI --config=config.cfg

For "container" run (assuming to have a caesar container image named `caesar.simg`) source finding can be run as:

.. code:: bash

    singularity run --app sfinder caesar.simg --config=config.cfg

Submission shell scripts, available in the `scripts` installation directory, enable running tasks on batch systems (PBS, SLURM).


That's it! See the Tutorial section for more examples.
