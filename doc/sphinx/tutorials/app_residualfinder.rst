
Remove/subtract compact sources from image 
==========================================

Say you have a FITS image (e.g. `input.fits`) with your radio observations and you want to produce a residual image with compact sources subtracted or removed. This is typically used as a pre-processing step before searching for extended sources.     
You can use the executable ``FindSourceResidual`` installed in the CAESAR `bin` directory to this aim:    

.. code:: bash

    $ ./FindSourceResidual [options]

    =========== USAGE ===========
    Usage: FindSourceResidual [options]

    *** Mandatory options ***    
    -i, --inputfile=[FILENAME] 	 Filename (fits/root) with input image.    
    
 
    *** Optional options ***     
    -h, --help 	 Show help message and exit    
    -s, --sourcefile=[FILENAME] 	 Caesar ROOT file with source list. If provided no sources will be searched.     
    -o, --outputfile=[FILENAME] 	 Filename where to store output residual image (default=resmap.fits)   
    -O, --outputfile_mask=[FILENAME] 	 Filename where to store output source mask image (default=smask.fits)     
    -p, --psSubtractionMethod=[METHOD] - Method used to remove point sources (1=DILATION,2=MODEL SUBTRACTION) (default=1)   
    -r, --resZThr=[NSIGMAS] - Significance threshold (in sigmas) above which sources are removed (if selected for removal) (default=5)   
    -R, --resZHighThr=[NSIGMAS] - Significance threshold (in sigmas) above which sources are always removed (even if they have nested or different type) (default=10)     
    -l, --removedSourceType=[TYPE] - Type of bright sources to be dilated from the input image (-1=ALL,1=COMPACT,2=POINT-LIKE,3=EXTENDED)    
    -a, --removeNestedSources - If a source has nested sources, remove nested rather than mother source (default=no)   
    -k, --dilateKernelSize=[SIZE] - Kernel size in pixel used to dilate image around sources (default=9)    
    -z, --bkgAroundSource - Use bkg computed in a box around source and not from the bkg map (default=use bkg map)    
    -Z, --bkgBoxThickness=[THICKNESS] - Bkg box thickness in pixels (default=20)   
    -j, --randomizeBkg - Randomize bkg in dilated pixels (default=no)   
    -T, --seedthr=[NSIGMAS] - Seed threshold in flood-fill algorithm in nsigmas significance (default=5)   
    -t, --mergethr=[NSIGMAS] - Merge threshold in flood-fill algorithm in nsigmas significance (default=2.6)    
    -m, --minnpixels=[NPIX] - Minimum number of pixels in a blob (default=5)    
    -N, --no-nested - Do not search nested sources (default=search)    
    -b, --boxsize=[SIZE] - Size of sampling box in pixels or expressed as a multiple of the image beam size (if --sizeinbeam option is given) (default=100 pixels)    
    -B, --sizeinbeam - Consider box size option expressed in multiple of beam size (beam info read from image) (default=no)   
    -g, --gridsize=[SIZE] - Size of the interpolation grid expressed as fraction of the sampling box (default=0.25)   
    -e, --estimator=[ESTIMATOR] - Bkg estimator used in the sampling box (1=mean, 2=median, 3=biweight, 4=clipped median) (default=2)    
    -P, --2ndpass - If given, perform a 2nd pass in bkg calculation (default=no)   
    -S, --skipblobs - If given, skip blobs using a flood-fill algorithm (default=no)   
    -A, --no-selection - Do not select and retag input sources (default=apply selection)   
    -d, --no-maxnpixcut - Do not apply max n pixel cut (default=apply)    
    -D, --maxnpix=[MAX_NPIX] - max number of pixel to consider a source as point-like (default=1000)   
    -f, --no-elongcut - Do not apply elongation cut (default=apply)    
    -G, --no-circratiocut - Do not apply circular ratio cut (default=apply)   
    -q, --no-ellipsearearatiocut - Do not apply ellipse area ratio cut (default=apply)    
    -u, --no-nbeamscut - Do not apply nbeams cut (default=apply)    
    -U, --maxnbeams=[MAX_NBEAMS] - Max number of beams in source to consider it as point-like (default=10)    
    -n, --nthreads 	 Number of threads to be used (default=1)   
    -v [LEVEL], --verbosity=[LEVEL] - Log level (<=0=OFF, 1=FATAL, 2=ERROR, 3=WARN, 4=INFO, >=5=DEBUG) (default=INFO)   

    Alternatively you can provide options in a configuration file with: 
    -c, --config=[FILENAME] 	 Config file containing option settings    

    ==============================

`FindSourceResidual` app mainly uses the CAESAR ``Image``, ``BkgFinder`` and ``MorphFilter`` class methods. 
The meaning of each command line option is briefly reported in the program help. The main algorithm steps are:    

1) Extract compact sources using flood-fill algorithm with options provided in the command line (``--seedthr``, ``--mergeThr``, etc) or configuration file. This step is skipped if you provide an input file with the list of sources to be removed (``--sourcefile`` option);    

2) Select sources to be removed from the input image. These are the sources passing the quality criteria, if selection is enabled (``--no-selection``, ``--no-maxnpixcut``, etc.), and sources above a chosen significance threshold. Sources above a high significance threshold (``--resZHighThr``) are always removed no matter what their type flag is. Sources above a lower significance threshold (``--resZThr``) are removed only if their type is equal to desired type (``--removedSourceType``). If ``--removeNestedSources`` option is given, only the nested sources are removed, not the mother source;   

3) Selected sources are "removed" as follows. Sources with available fit information are "subtracted" if ``--psSubtractionMethod``=2 option is given, e.g. the fitted gaussian model is subtracted from the input image. In the other cases, source pixel and neighbouring pixel (depending on the dilation kernel size) values are replaced with the estimated background values. Background values can be obtained either from the computed local background map (see `FindBkg` app) or from the pixel median value (excluding source pixels) computed in a box centred around the source (if ``--bkgAroundSource`` option is given);   

This application produces two outputs:   
 
- `resMap`: a FITS file with residual image;    
- `mask`: a FITS file with a binary source mask;   

