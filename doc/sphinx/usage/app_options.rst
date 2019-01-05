
Configuration Options 
=====================

Most of the example applications provided can be configured from command line arguments, as described in the Tutorial section.
Some applications, like source finding, however, require a large set of configuration options, specified inside a configuration file, passed to the application
as a command line argument ``--config=[FILE]``.

In this section we report a list of the main configuration options defined in CAESAR to customize tasks. The full list of options defined is kept in ``ConfigParser.cc`` class.
To print the full list of defined options use the `ConfigParser::PrintOptions()` method. For example from ROOT prompt type:


``Caesar::ConfigParser::Instance().PrintOptions()``


or from the python CLI:

.. code:: python

    from ROOT import gSystem                     
    gSystem.Load('libCaesar')
    from ROOT import Caesar
    
    Caesar.ConfigParser.Instance().PrintOptions()



-------------
Input Options
-------------

These options enable control of input data to be given to CAESAR applications.  

+------------------------+-------------------------------------------+--------------+------------+
|       Option           |             Description                   |   Default    |   Values   |
+========================+===========================================+==============+============+
| ``inputFile``          | Input image filename (.root/.fits)        |     ""       |            |
+------------------------+-------------------------------------------+--------------+------------+
| ``inputImage``         | Image name to be read in input ROOT file  |     ""       |            |
+------------------------+-------------------------------------------+--------------+------------+
| ``readTileImage``      | | Read sub-image                          |    false     | | true     |
|                        | | If false read the entire image          |              | | false    |
+------------------------+-------------------------------------------+--------------+------------+
| ``tileMinX``           | | Min image x pixel coordinate to be read |      0       |            |
|                        | | Used only when readTileImage is true    |              |            |
+------------------------+-------------------------------------------+--------------+------------+
| ``tileMaxX``           | | Max image x pixel coordinate to be read |      0       |            |
|                        | | Used only when readTileImage is true    |              |            |
+------------------------+-------------------------------------------+--------------+------------+
| ``tileMinY``           | | Min image y pixel coordinate to be read |      0       |            |
|                        | | Used only when readTileImage is true    |              |            |
+------------------------+-------------------------------------------+--------------+------------+
| ``tileMinY``           | | Max image y pixel coordinate to be read |      0       |            |
|                        | | Used only when readTileImage is true    |              |            |
+------------------------+-------------------------------------------+--------------+------------+

--------------
Output Options
--------------

These options enable control of information & data reported in output by CAESAR applications.  

+--------------------------------+----------------------------------+-----------------------+-------------+
|       Option                   |             Description          |      Default          |   Values    |
+================================+==================================+=======================+=============+
| ``saveToFile``                 | | Save results & maps to output  |        true           | | true      |
|                                | | ROOT file                      |                       | | false     |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``saveConfig``                 | | Save config options to output  |        true           | | true      |
|                                | | ROOT file                      |                       | | false     |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``saveSources``                | | Save sources to output ROOT    |        true           | | true      |
|                                | | file                           |                       | | false     |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``saveInputMap``               | | Save input map to output       |        false          | | true      |
|                                | | ROOT file                      |                       | | false     |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``saveBkgMap``                 | | Save computed background map   |        true           | | true      |
|                                | | to output ROOT file            |                       | | false     |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``saveNoiseMap``               | | Save computed rms map          |        true           | | true      |
|                                | | to output ROOT file            |                       | | false     |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``saveResidualMap``            | | Save computed residual map     |        true           | | true      |
|                                | | to output ROOT file            |                       | | false     |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``saveSignificanceMap``        | | Save computed significance map |        true           | | true      |
|                                | | to output ROOT file            |                       | | false     |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``saveSaliencyMap``            | | Save computed saliency map     |        true           | | true      |
|                                | | to output ROOT file            |                       | | false     |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``saveSegmentedMap``           | | Save computed segmented map    |        true           | | true      |
|                                | | to output ROOT file            |                       | | false     |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``outputFile``                 | | Name of ROOT file where to     |    `output.root`      |             |
|                                | | store output data (images,     |                       |             |
|                                | | run config, sources, etc)      |                       |             |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``saveToCatalogFile``          | | Save island and fitted         |        true           | | true      |
|                                | | components to ascii files      |                       | | false     |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``outputCatalogFile``          | | Name of ascii file where to    |    `catalog.dat`      |             |
|                                | | store source catalog           |                       |             |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``outputComponentCatalogFile`` | | Name of ascii file where to    | `catalog_fitcomp.dat` |             |
|                                | | store fitted source component  |                       |             |
|                                | | catalog                        |                       |             |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``saveDS9Region``              | | Save sources & fit components  |        true           | | true      |
|                                | | to DS9 region files            |                       | | false     |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``ds9RegionFile``              | | Name of DS9 region file where  |      `ds9.reg`        |             |
|                                | | to store source catalog        |                       |             |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``ds9FitRegionFile``           | | Name of ascii file where to    |   `ds9_fitcomp.reg`   |             |
|                                | | store source fitted components |                       |             |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``convertDS9RegionsToWCS``     | | Store DS9 regions in WCS       |        false          | | true      |
|                                | | coordinates                    |                       | | false     |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``ds9WCSType``                 | | DS9 region WCS type to be used |         0             | | 0=J2000   |
|                                | | if convertDS9RegionsToWCS=true |                       | | 1=B1950   |
|                                |                                  |                       | | 2=GAL     |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``ds9RegionFormat``            | | Shape to be used to store      |         2             | | 1=ellipse |
|                                | | source islands in DS9 region   |                       | | 2=polygon |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``saveToFITSFile``             | | Save output data images to     |        false          | | true      |
|                                | | FITS files                     |                       | | false     |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``inputMapFITSFile``           | | Name of FITS file where        |     `input.fits`      |             |
|                                | | to store input map read        |                       |             |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``bkgMapFITSFile``             | | Name of FITS file where        |      `bkg.fits`       |             |
|                                | | to store computed bkg map      |                       |             |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``noiseMapFITSFile``           | | Name of FITS file where        |      `rms.fits`       |             |
|                                | | to store computed rms map      |                       |             |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``significanceMapFITSFile``    | | Name of FITS file where to     |  `significance.fits`  |             |
|                                | | store computed significance map|                       |             |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``residualMapFITSFile``        | | Name of FITS file where        |    `residual.fits`    |             |
|                                | | to store computed residual map |                       |             |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``saliencyMapFITSFile``        | | Name of FITS file where        |    `saliency.fits`    |             |
|                                | | to store computed saliency map |                       |             |
+--------------------------------+----------------------------------+-----------------------+-------------+

------------------------------------
Run & Distributed Processing Options
------------------------------------

These options enable control of application run (e.g. logging levels) and distributed processing (e.g. number of threads). 

+--------------------------------+----------------------------------+-----------------------+-------------+
|       Option                   |             Description          |      Default          |   Values    |
+================================+==================================+=======================+=============+
| ``logLevel``                   | Log level threshold              |        INFO           | | DEBUG     |
|                                |                                  |                       | | INFO      |
|                                |                                  |                       | | WARN      |
|                                |                                  |                       | | ERROR     |
|                                |                                  |                       | | FATAL     |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``nThreads``                   | | Number of threads used if      |        -1             |             |
|                                | | OPENMP is enabled. If set to   |                       |             |
|                                | | -1 a number of threads equal   |                       |             |
|                                | | to the available cores is used |                       |             |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``splitInTiles``               | | Split input image in tiles     |       false           | | true      |
|                                | | for parallel processing        |                       | | false     |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``tileSizeX``                  | | Size of tile in pixels along X |        1000           |             |
|                                | | coordinate used for partition  |                       |             |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``tileSizeY``                  | | Size of tile in pixels along Y |        1000           |             |
|                                | | coordinate used for partition  |                       |             |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``useTileOverlap``             | | Enable tile overlap in image   |        false          | | true      |
|                                | | partition for parallel         |                       | | false     |
|                                | | processing                     |                       |             |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``tileStepSizeX``              | | Tile overlap fraction along    |          1            |             |
|                                | | X coordinate to partition the  |                       |             |
|                                | | input image for parallel       |                       |             |
|                                | | processing (1=no overlap,      |                       |             |
|                                | | 0.5=half overlap)              |                       |             |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``tileStepSizeY``              | | Tile overlap fraction along    |          1            |             |
|                                | | Y coordinate to partition the  |                       |             |
|                                | | input image for parallel       |                       |             |
|                                | | processing (1=no overlap,      |                       |             |
|                                | | 0.5=half overlap)              |                       |             |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``mergeSourcesAtEdge``         | | Merge overlapping sources      |         true          | | true      |
|                                | | found at tile edge by each     |                       | | false     |
|                                | | worker when aggregating the    |                       |             |
|                                | | final catalog                  |                       |             |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``mergeSources``               | | Merge overlapping sources      |         false         | | true      |
|                                | | found in each tile. If true    |                       | | false     |
|                                | | compact and extended sources   |                       |             |
|                                | | found by different algorithms  |                       |             |
|                                | | in a tile are merged if        |                       |             |
|                                | | overlapping. If you want to    |                       |             |
|                                | | keep sources distinct set      |                       |             |
|                                | | option to false                |                       |             |
+--------------------------------+----------------------------------+-----------------------+-------------+


----------------------------------
Stats & Background Compute Options
----------------------------------

These options enable control of image background calculation. Background can be either computed globally or locally.
Local background maps (bkg, rms) are obtained by interpolating background estimator values computed on a grid of sampling image rectangular boxes.

+--------------------------------+----------------------------------+-----------+------------------------+
|       Option                   |             Description          |  Default  |   Values               |
+================================+==================================+===========+========================+
| ``bkgEstimator``               | | Stat estimator used to compute |    2      | | 1=Mean/RMS           |
|                                | | image background & noise       |           | | 2=Median/MAD         |
|                                | | image background & noise       |           | | 3=BiWeight           |
|                                | | image background & noise       |           | | 4=Clipped Median/RMS |
+--------------------------------+----------------------------------+-----------+------------------------+
| ``useParallelMedianAlgo``      | | Use C++ parallel algorithm     |   true    | | true                 |
|                                | | to compute median estimator    |           | | false                |
+--------------------------------+----------------------------------+-----------+------------------------+
| ``useLocalBkg``                | | Compute local background       |   true    | | true                 |
|                                | | and noise maps and use them    |           | | false                |
|                                | | instead of global bkg info     |           |                        |
+--------------------------------+----------------------------------+-----------+------------------------+
| ``use2ndPassInLocalBkg``       | | Use 2nd pass to refine local   |   true    | | true                 |
|                                | | noise map                      |           | | false                |
+--------------------------------+----------------------------------+-----------+------------------------+
| ``skipOutliersInLocalBkg``     | | Exclude pixels belonging to    |   false   | | true                 |
|                                | | detected bright blobs when     |           | | false                |
|                                | | computing local background     |           |                        |
|                                | | estimators. Blob find seed thr |           |                        |
|                                | | parameters are reported in     |           |                        |
|                                | | source finding option table    |           |                        |
|                                | | below                          |           |                        |
+--------------------------------+----------------------------------+-----------+------------------------+
| ``boxSizeX``                   | | Size of sampling box along x   |    20     |                        |
|                                | | coordinate for local bkg       |           |                        |
|                                | | calculation in pixels. Size is |           |                        |
|                                | | instead assumed as multiple of |           |                        |
|                                | | beam size if                   |           |                        |
|                                | | ``useBeamInfoInBkg`` is true   |           |                        |
+--------------------------------+----------------------------------+-----------+------------------------+
| ``boxSizeY``                   | | Size of sampling box along y   |    20     |                        |
|                                | | coordinate for local bkg       |           |                        |
|                                | | calculation in pixels. Size is |           |                        |
|                                | | instead assumed as multiple of |           |                        |
|                                | | beam size if                   |           |                        |
|                                | | ``useBeamInfoInBkg`` is true   |           |                        |
+--------------------------------+----------------------------------+-----------+------------------------+
| ``gridSizeX``                  | | Size of grid along x           |    0.2    |                        |
|                                | | coordinate used for local bkg  |           |                        |
|                                | | interpolation expressed as     |           |                        |
|                                | | fraction of sampling box x     |           |                        |
|                                | | size                           |           |                        |
+--------------------------------+----------------------------------+-----------+------------------------+
| ``gridSizeY``                  | | Size of grid along y           |    0.2    |                        |
|                                | | coordinate used for local bkg  |           |                        |
|                                | | interpolation expressed as     |           |                        |
|                                | | fraction of sampling box y     |           |                        |
|                                | | size                           |           |                        |
+--------------------------------+----------------------------------+-----------+------------------------+
| ``sourceBkgBoxBorderSize``     | | Border pad size in pixels of   |    20     |                        |
|                                | | box around source bounding box |           |                        |
|                                | | used to estimate bkg for       |           |                        |
|                                | | fitting                        |           |                        |
+--------------------------------+----------------------------------+-----------+------------------------+
| ``useBeamInfoInBkg``           | | Use beam information in bkg    |   true    | | true                 |
|                                | | sampling box size definition.  |           | | false                |
|                                | | Beam info are taken from image |           |                        |
|                                | | when available, otherwise from |           |                        |
|                                | | user beam parameter below.     |           |                        |
+--------------------------------+----------------------------------+-----------+------------------------+
| ``pixSize``                    | | User-supplied map pixel area   |     1     |                        |
|                                | | in arcsec. Used when CDELT     |           |                        |
|                                | | info is not available in       |           |                        |
|                                | | image metadata                 |           |                        |
+--------------------------------+----------------------------------+-----------+------------------------+
| ``beamFWHM``                   | | User-supplied circular beam    |    6.5    |                        |
|                                | | FWHM in arcsec (BMAJ=BMIN).    |           |                        |
|                                | | Used when beam info is not     |           |                        |
|                                | | available in image metadata    |           |                        |
+--------------------------------+----------------------------------+-----------+------------------------+
| ``beamBmaj``                   | | User-supplied beam ellipse     |    10     |                        |
|                                | | major axis in arcsec.          |           |                        |
|                                | | Used when beam info is not     |           |                        |
|                                | | available in image metadata    |           |                        |
+--------------------------------+----------------------------------+-----------+------------------------+
| ``beamBmin``                   | | User-supplied beam ellipse     |     5     |                        |
|                                | | minor axis in arcsec.          |           |                        |
|                                | | Used when beam info is not     |           |                        |
|                                | | available in image metadata    |           |                        |
+--------------------------------+----------------------------------+-----------+------------------------+
| ``beamTheta``                  | | User-supplied beam position    |     0     |                        |
|                                | | angle in degrees and measured  |           |                        |
|                                | | CCW from North (pa=0 North).   |           |                        |
|                                | | Used when beam info is not     |           |                        |
|                                | | available in image metadata    |           |                        |
+--------------------------------+----------------------------------+-----------+------------------------+


----------------------
Source Finding Options
----------------------

These options enable control of source detection. This is performed using a flood-fill algorithm
aggregating pixels around significant seeds if above a given merge threshold. Detected blobs form a collection
of candidate sources.

+--------------------------------+----------------------------------+-----------+------------------------+
|       Option                   |             Description          |  Default  |   Values               |
+================================+==================================+===========+========================+
| ``searchCompactSources``       | | Enable/disable search of       |   true    | | true                 |
|                                | | compact sources                |           | | false                |
+--------------------------------+----------------------------------+-----------+------------------------+
| ``minNPix``                    | | Minimum number of pixels       |    5      |                        |
|                                | | to consider a blob as source   |           |                        |
|                                | | candidate                      |           |                        |
+--------------------------------+----------------------------------+-----------+------------------------+
| ``seedThr``                    | | Seed threshold in blob finding |    5      |                        |
|                                | | given as number of sigmas      |           |                        |
|                                | | above background               |           |                        |
+--------------------------------+----------------------------------+-----------+------------------------+
| ``mergeThr``                   | | Merge/aggregation threshold    |   2.6     |                        |
|                                | | in blob finding given as       |           |                        |
|                                | | number of sigmas above         |           |                        |
|                                | | background. Pixels above this  |           |                        |
|                                | | threshold are added to the blob|           |                        |
+--------------------------------+----------------------------------+-----------+------------------------+
| ``mergeBelowSeed``             | | Add to blob only pixels above  |   false   | | true                 |
|                                | | merge threshold but below seed |           | | false                |
|                                | | threshold                      |           |                        |
+--------------------------------+----------------------------------+-----------+------------------------+
| ``searchNegativeExcess``       | | Search for holes (i.e. blobs   |   false   | | true                 |
|                                | | with negative significance)    |           | | false                |
|                                | | along with "positive" blobs    |           |                        |
+--------------------------------+----------------------------------+-----------+------------------------+
| ``compactSourceSearchNIters``  | | Number of iterations to be     |     1     |                        |
|                                | | performed in compact source    |           |                        |
|                                | | search. At each iteration the  |           |                        |
|                                | | seed threshold is decreased by |           |                        |
|                                | | ``seedThrStep``                |           |                        |
+--------------------------------+----------------------------------+-----------+------------------------+
| ``seedThrStep``                | | Seed threshold decrease step   |    0.5    |                        |
|                                | | size between iterations.       |           |                        |
|                                | | Effective only when            |           |                        |
|                                | | ``compactSourceSearchNIters``>1|           |                        |
+--------------------------------+----------------------------------+-----------+------------------------+

		
-----------------------------
Nested Source Finding Options
-----------------------------

These options enable control of nested source detection. Nested sources are blobs inside another mother blobs.
Detection of nested blob uses a blob detection algorithm, based on the thresholding of a filter blob map (LoG or Gaus2D smoothed),
which increases the computation time, particularly if blob search is done at multiple spatial scales. In presence of extended/diffuse object you can consider turning off
this calculation. If however you have extended and bright object and you turn off nested source search you may see that 
compact/point-source located inside the extended one will be included in the mother and not fitted.

+---------------------------------------+----------------------------------+-----------+------------------------+
|       Option                          |             Description          |  Default  |   Values               |
+=======================================+==================================+===========+========================+
| ``searchNestedSources``               | | Enable/disable search of       |   true    | | true                 |
|                                       | | compact nested sources         |           | | false                |
+---------------------------------------+----------------------------------+-----------+------------------------+
| ``blobMaskMethod``                    | | Filter map used in nested      |    2      | | 1=gaus smoothed Lapl |
|                                       | | blob finder to search blobs    |           | | 2=multi-scale LoG    |
+---------------------------------------+----------------------------------+-----------+------------------------+
| ``nestedBlobKernFactor``              | | Filter kernel size factor par  |    6      |                        |
|                                       | | so that kern size=             |           |                        |
|                                       | | factor x sigma (sigma is the   |           |                        |
|                                       | | filter scale par in pixels)    |           |                        |
+---------------------------------------+----------------------------------+-----------+------------------------+
| ``sourceToBeamAreaThrToSearchNested`` | | Mother source area/beam thr to |    10     |                        |
|                                       | | add nested sources. If         |           |                        |
|                                       | | npix<=thr*beamArea no nested   |           |                        |
|                                       | | sources are added to the       |           |                        |
|                                       | | mother source even if detected.|           |                        |
|                                       | | If thr=0 nested sources are    |           |                        |
|                                       | | always added if                |           |                        |
|                                       | | ``searchNestedSources`` is     |           |                        |
|                                       | | enabled                        |           |                        |
+---------------------------------------+----------------------------------+-----------+------------------------+
| ``nestedBlobThrFactor``               | | Threshold factor param used in |    0      |                        |
|                                       | | blob filter map to create mask |           |                        |
|                                       | | (thr=thrFactor*<img>).         |           |                        |
+---------------------------------------+----------------------------------+-----------+------------------------+
| ``minNestedMotherDist``               | | Minimum distance in pixels     |    2      |                        |
|                                       | | (in x or y) between nested and |           |                        |
|                                       | | parent blob centroids below    |           |                        |
|                                       | | which nested source is skipped |           |                        |
|                                       | | as most probably equal to the  |           |                        |
|                                       | | parent (avoid duplicates)      |           |                        |
+---------------------------------------+----------------------------------+-----------+------------------------+
| ``maxMatchingPixFraction``            | | Maximum fraction of matching   |   0.5     |                        |
|                                       | | pixels between nested and      |           |                        |
|                                       | | parent blob above which nested |           |                        |
|                                       | | is skipped as most probably    |           |                        |
|                                       | | equal to the parent (avoid     |           |                        |
|                                       | | duplicates)                    |           |                        |
+---------------------------------------+----------------------------------+-----------+------------------------+
| ``nestedBlobPeakZThr``                | | Nested blob significance       |    5      |                        |
|                                       | | seed thr in sigmas (in filter  |           |                        |
|                                       | | blob map) below which nested   |           |                        |
|                                       | | blob is skipped                |           |                        |
+---------------------------------------+----------------------------------+-----------+------------------------+
| ``nestedBlobPeakZMergeThr``           | | Nested blob peak significance  |   2.5     |                        |
|                                       | | merge thr in sigmas (in filter |           |                        |
|                                       | | blob map) below which nested   |           |                        |
|                                       | | blob is skipped                |           |                        |
+---------------------------------------+----------------------------------+-----------+------------------------+
| ``nestedBlobMinScale``                | | Nested blob min search scale   |    1      |                        |
|                                       | | factor parameter so that blob  |           |                        |
|                                       | | filter scale in pixels is      |           |                        |
|                                       | | = scaleFactor x beam width     |           |                        |
+---------------------------------------+----------------------------------+-----------+------------------------+
| ``nestedBlobMaxScale``                | | Nested blob max search scale   |    3      |                        |
|                                       | | factor parameter so that blob  |           |                        |
|                                       | | filter scale in pixels is      |           |                        |
|                                       | | = scaleFactor x beam width     |           |                        |
+---------------------------------------+----------------------------------+-----------+------------------------+
| ``nestedBlobScaleStep``               | | Nested blob scale factor step  |    1      |                        |
|                                       | | so that scaleFactor=           |           |                        |
|                                       | | minScaleFactor + step          |           |                        |
+---------------------------------------+----------------------------------+-----------+------------------------+


------------------------
Source Selection Options
------------------------

These options enable control of quality selection cuts applied to detected blobs to select good source candidates and tag point-source candidates
(used later in source residual map and fitting stage). Options are also provided to select sources to be stored in the final catalog. 

+---------------------------------------+----------------------------------+-----------+------------------------+
|       Option                          |             Description          |  Default  |   Values               |
+=======================================+==================================+===========+========================+
| ``applySourceSelection``              | | Enable/disable source          |   true    | | true                 |
|                                       | | selection                      |           | | false                |
+---------------------------------------+----------------------------------+-----------+------------------------+
| ``useMinBoundingBoxCut``              | | Apply minimum bounding box cut |   false   | | true                 |
|                                       | | to detected blobs              |           | | false                |
+---------------------------------------+----------------------------------+-----------+------------------------+
| ``sourceMinBoundingBox``              | | Minimum bounding box cut value |    2      |                        |
|                                       | | in pixel. Blobs with minimum   |           |                        |
|                                       | | bounding box size below the    |           |                        |
|                                       | | threshold are tagged as bad    |           |                        |
+---------------------------------------+----------------------------------+-----------+------------------------+
| ``useCircRatioCut``                   | | Apply cut on blob circular     |   false   | | true                 |
|                                       | | ratio param to detected blobs  |           | | false                |
+---------------------------------------+----------------------------------+-----------+------------------------+
| ``psCircRatioThr``                    | | Circular ratio cut value.      |    0.4    | 0                      |
|                                       | | in pixel. Blobs with circ      |           | 1                      |
|                                       | | ratio above this threshold     |           |                        |
|                                       | | passed the point-like cut      |           |                        |
|                                       | | (1=circle, 0=line)             |           |                        |
+---------------------------------------+----------------------------------+-----------+------------------------+
| ``useElongCut``                       | | Apply cut on blob elongation   |   false   | | true                 |
|                                       | | param to detected blobs        |           | | false                |
+---------------------------------------+----------------------------------+-----------+------------------------+
| ``psElongThr``                        | | Elongation cut value.          |    0.7    | 0                      |
|                                       | | Blobs with elongation param    |           | 1                      |
|                                       | | below this threshold           |           |                        |
|                                       | | passed the point-like cut      |           |                        |
+---------------------------------------+----------------------------------+-----------+------------------------+
| ``useMaxNPixCut``                     | | Apply cut on blob maximum      |   false   | | true                 |
|                                       | | number of pixels.              |           | | false                |
+---------------------------------------+----------------------------------+-----------+------------------------+
| ``psMaxNPix``                         | | Max number of pixels cut value.|   1000    |                        |
|                                       | | Blobs with a number of pixels  |           |                        |
|                                       | | below this threshold           |           |                        |
|                                       | | passed the point-like cut      |           |                        |
+---------------------------------------+----------------------------------+-----------+------------------------+
| ``useEllipseAreaRatioCut``            | | Apply cut on ratio between     |   false   | | true                 |
|                                       | | blob area and blob ellipse     |           | | false                |
|                                       | | bounding box area.             |           |                        |
+---------------------------------------+----------------------------------+-----------+------------------------+
| | ``psEllipseAreaRatioMinThr``        | | Area/EllipseArea ratio min and |   0.6     |                        |
| | ``psEllipseAreaRatioMaxThr``        | | max cut values.                |   1.4     |                        |
|                                       | | Blobs in cut range passes the  |           |                        |
|                                       | | point-like cut                 |           |                        |
+---------------------------------------+----------------------------------+-----------+------------------------+
| ``useNBeamsCut``                      | | Apply cut on number of beams   |   false   | | true                 |
|                                       | | found in detected blob         |           | | false                |
|                                       | | (NBeams=blob npix/beam npix)   |           |                        |
+---------------------------------------+----------------------------------+-----------+------------------------+
| ``psNBeamsThr``                       | | Max number of beams cut value. |    10     |                        |
|                                       | | Blobs with a number of beams   |           |                        |
|                                       | | below this threshold           |           |                        |
|                                       | | passed the point-like cut      |           |                        |
+---------------------------------------+----------------------------------+-----------+------------------------+



----------------------
Source Fitting Options
----------------------

These options enable control of source fitting stage: minimization algorithm and relative parameters, starting parameters and limits, etc.

+---------------------------------+----------------------------------+-----------+------------------------+
|       Option                    |             Description          |  Default  |   Values               |
+=================================+==================================+===========+========================+
| ``fitSources``                  | | Enable/disable source          |   false   | | true                 |
|                                 | | fitting stage                  |           | | false                |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fitUseThreads``               | | Split source fitting among     |   false   | | true                 |
|                                 | | multiple threads. Multithread  |           | | false                |
|                                 | | is not supported by Minuit     |           |                        |
|                                 | | minimizer                      |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fitScaleDataToMax``           | | Scale source flux data to max  |   false   | | true                 |
|                                 | | peak flux if true, otherwise   |           | | false                |
|                                 | | scale to mJy units             |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fitMinimizer``                | | Minimizer used in source       |  Minuit2  | | Minuit               |
|                                 | | fitting                        |           | | Minuit2              |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fitMinimizerAlgo``            | | Minimization algorithm used in |  minimize | | migrad               |
|                                 | | source fitting                 |           | | simplex              |
|                                 |                                  |           | | scan                 |
|                                 |                                  |           | | minimize             |
|                                 |                                  |           | | fumili               |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fitPrintLevel``               | | Minimizer printout level       |     1     |                        |
|                                 | |                                |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fitStrategy``                 | | Minimizer strategy parameter   |     2     |                        |
|                                 | | (larger means more accurate    |           |                        |
|                                 | | minimization but more fcn      |           |                        |
|                                 | | calls)                         |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fitFcnTolerance``             | | Fit function minimization      |   1.e-2   |                        |
|                                 | | tolerance (smaller means more  |           |                        |
|                                 | | accurate minimization but more |           |                        |
|                                 | | fcn calls)                     |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fitMaxIters``                 | | Maximum number of iterations   |   100000  |                        |
|                                 | | that can be done by minimizer  |           |                        |
|                                 | | before giving up and returning |           |                        |
|                                 | | not converged fit              |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fitImproveConvergence``       | | Try to improve convergence by  |   true    | | true                 |
|                                 | | iterating fit if not converged |           | | false                |
|                                 | | or converged with pars at      |           |                        |
|                                 | | limits                         |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fitNRetries``                 | | Number of times fit is         |   1000    |                        |
|                                 | | repeated (with enlarged        |           |                        |
|                                 | | limits) if improve convergence |           |                        |
|                                 | | flag is enabled                |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fitParBoundIncreaseStepSize`` | | Par bound rel increase step    |   0.1     |                        |
|                                 | | size set when trying to improve|           |                        |
|                                 | | convergence:                   |           |                        |
|                                 | | parmax= parmax_old+(1+nretry)* |           |                        |
|                                 | |  *fitParBoundIncreaseStepSize* |           |                        |
|                                 | |  *0.5*|max-min|                |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fitDoFinalMinimizerStep``     | | If enabled run HESSE minimizer |   true    | | true                 |
|                                 | | at convergence to improve      |           | | false                |
|                                 | | minimum and par error estimate |           |                        |
|                                 | | limits                         |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fitRetryWithLessComponents``  | | If fit fails to converge,      |   true    | | true                 |
|                                 | | repeat it iteratively with one |           | | false                |
|                                 | | component less at each cycle   |           |                        |
|                                 | | until convergence or until no  |           |                        |
|                                 | | more components are available  |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``nBeamsMaxToFit``              | | Maximum number of beams        |    20     |                        |
|                                 | | for a compact source to be     |           |                        |
|                                 | | fitted (if above this threshold|           |                        |
|                                 | | the fit is not performed)      |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fitUseNestedAsComponents``    | | If true use nested sources     |   false   | | true                 |
|                                 | | (if any) as starting fit       |           | | false                |
|                                 | | components, otherwise estimate |           |                        |
|                                 | | blended components in blob     |           |                        |
|                                 | | using a peak finding +         |           |                        |
|                                 | | segmentation algorithm         |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fitMaxNComponents``           | | Maximum number of components   |     5     |                        |
|                                 | | fitted in a blob               |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``peakMinKernelSize``           | | Minimum dilation kernel size   |     3     |                        |
|                                 | | in pixels used to detect start |           |                        |
|                                 | | fit components                 |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``peakMaxKernelSize``           | | Maximum dilation kernel size   |     7     |                        |
|                                 | | in pixels used to detect start |           |                        |
|                                 | | fit components                 |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``peakKernelMultiplicityThr``   | | Requested peak multiplicity    |     1     |                        |
|                                 | | across different dilation      |           |                        |
|                                 | | kernels. A multiplicity=-1     |           |                        |
|                                 | | imposes that a peak must be    |           |                        |
|                                 | | found in all given dilation    |           |                        |
|                                 | | kernels (within a tolerance)   |           |                        |
|                                 | | to be considered a component   |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``peakShiftTolerance``          | | Peak max position offset in    |     2     |                        |
|                                 | | pixels above which two peaks   |           |                        |
|                                 | | are considered distincs.       |           |                        |
|                                 | | Used to compare peaks found    |           |                        |
|                                 | | in different dilation kernels  |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``peakZThrMin``                 | | Minimum peak flux significance |     1     |                        |
|                                 | | (in nsigmas wrt source avg     |           |                        |
|                                 | | bkg and rms) below which peak  |           |                        |
|                                 | | is skipped and not considered  |           |                        |
|                                 | | as a fit component             |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fitWithCentroidLimits``       | | Apply limits to source         |   true    | | true                 |
|                                 | | centroid pars in fit           |           | | false                |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fixCentroidInPreFit``         | | Fix source centroid pars       |   false   | | true                 |
|                                 | | in pre-fit                     |           | | false                |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fitCentroidLimit``            | | Source centroid par limits     |     3     |                        |
|                                 | | given as max offset in pixel   |           |                        |
|                                 | | with respect to starting fit   |           |                        |
|                                 | | centroid pars                  |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fitWithFixedBkg``             | | Fix bkg level par in fit       |   true    | | true                 |
|                                 |                                  |           | | false                |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fitWithBkgLimits``            | | Apply limits to bkg level par  |   true    | | true                 |
|                                 | | in fit                         |           | | false                |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fitUseEstimatedBkgLevel``     | | Use estimated (avg bkg) as     |   true    | | true                 |
|                                 | | starting bkg level par in fit  |           | | false                |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fitUseBkgBoxEstimate``        | | Use bkg estimated in a box     |   true    | | true                 |
|                                 | | around source (if available)   |           | | false                |
|                                 | | as bkg level par in fit        |           | | false                |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fitBkgLevel``                 | | Starting bkg level par in fit  |     0     |                        |
|                                 | | (used when option              |           |                        |
|                                 | | fitParBoundIncreaseStepSize is |           |                        |
|                                 | | false                          |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fitWithAmplLimits``           | | Apply limits to amplitude par  |   true    | | true                 |
|                                 | | in fit                         |           | | false                |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fixAmplInPreFit``             | | Fix amplitude par in pre-fit   |   true    | | true                 |
|                                 | |                                |           | | false                |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fitAmplLimit``                | | Amplitude par limit given as   |   0.3     |                        |
|                                 | | max relative offset with       |           |                        |
|                                 | | respect to starting source     |           |                        |
|                                 | | component peak                 |           |                        |
|                                 | | Speak*(1+-fitAmplLimit))       |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fitWithSigmaLimits``          | | Apply limits to sigma pars     |   true    | | true                 |
|                                 | | in fit                         |           | | false                |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fixSigmaInPreFit``            | | Fix sigma pars in pre-fit      |   false   | | true                 |
|                                 | |                                |           | | false                |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fitSigmaLimit``               | | Sigma par limit given as max   |   0.3     |                        |
|                                 | | relative offset with respect   |           |                        |
|                                 | | to starting component sigma    |           |                        |
|                                 | | pars                           |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fitWithFixedSigma``           | | Fix sigma pars in fit          |   false   | | true                 |
|                                 | |                                |           | | false                |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fitWithThetaLimits``          | | Apply limits to theta par      |   true    | | true                 |
|                                 | | in fit                         |           | | false                |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fixThetaInPreFit``            | | Fix theta par in pre-fit       |   false   | | true                 |
|                                 | |                                |           | | false                |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fitWithFixedTheta``           | | Fix theta par in fit           |   false   | | true                 |
|                                 | |                                |           | | false                |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fitThetaLimit``               | | Theta par limit given as max   |     5     |                        |
|                                 | | offset in degrees with respect |           |                        |
|                                 | | to starting component theta    |           |                        |
|                                 | | par                            |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``useFluxZCutInFit``            | | If enabled only blob pixels    |   false   | | true                 |
|                                 | | above a significance threshold |           | | false                |
|                                 | | are included in chi2. Pixels   |           |                        |
|                                 | | below threshold are included   |           |                        |
|                                 | | in a regularization chi2 term  |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fitZCutMin``                  | | Blob significance              |    2.5    |                        |
|                                 | | threshold below which pixels   |           |                        |
|                                 | | are included in the            |           |                        |
|                                 | | regularization chi2 term but   |           |                        |
|                                 | | not in the chi2                |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``fitChi2RegPar``               | | Fit chi2 regularization par    |     0     |                        |
|                                 | | so that total chi2 is given by |           |                        |
|                                 | | chi2(Z>thr)+regPar*chi2(Z<thr) |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+


--------------------------
Source Fit Selection Cuts
--------------------------

These options enable control of source fit selection cuts. These cuts are used to assign flag to source fitted components.

+----------------------------------------+----------------------------------+-----------+--------------+
|       Option                           |             Description          |  Default  |   Values     |
+========================================+==================================+===========+==============+
| ``fitApplyRedChi2Cut``                 | | Apply fit Chi2/NDF cut.        |   true    | | true       |
|                                        | | Used to set fit quality flag.  |           | | false      |
|                                        | | If Chi2/NDF>cut the good fit   |           |              |
|                                        | | cut is not passed              |           |              |
+----------------------------------------+----------------------------------+-----------+--------------+
| ``fitRedChi2Cut``                      | | Chi2/NDF cut value             |     5     |              |
+----------------------------------------+----------------------------------+-----------+--------------+
| ``fitApplyFitEllipseCuts``             | | Apply fit ellipse selection    |   false   | | true       |
|                                        | | cuts. Used to set component    |           | | false      |
|                                        | | flags. If not passed, fit      |           |              |
|                                        | | component is tagged as fake    |           |              |
+----------------------------------------+----------------------------------+-----------+--------------+
| ``fitEllipseEccentricityRatioMinCut``  | | Ellipse eccentricity ratio     |    0.5    |              |
|                                        | | (fit/beam) min cut value       |           |              |
+----------------------------------------+----------------------------------+-----------+--------------+
| ``fitEllipseEccentricityRatioMaxCut``  | | Ellipse eccentricity ratio     |    1.5    |              |
|                                        | | (fit/beam) max cut value       |           |              |
+----------------------------------------+----------------------------------+-----------+--------------+
| ``fitEllipseAreaRatioMinCut``          | | Ellipse area ratio             |    0.01   |              |
|                                        | | (fit/beam) min cut value       |           |              |
+----------------------------------------+----------------------------------+-----------+--------------+
| ``fitEllipseAreaRatioMaxCut``          | | Ellipse area ratio             |    10     |              |
|                                        | | (fit/beam) max cut value       |           |              |
+----------------------------------------+----------------------------------+-----------+--------------+
| ``fitEllipseRotAngleCut``              | | Ellipse rot angle diff         |    45     |              |
|                                        | | (|fit-beam|) cut value         |           |              |
|                                        | | in degrees                     |           |              |
+----------------------------------------+----------------------------------+-----------+--------------+

	
-----------------------
Source Residual Options
-----------------------

These options enable control of source residual map. Residual map is made by removing and/or subtracting detected sources from the input map. 
Source removal is done by replacing source pixel flux values (along with surrounding pixel around them, controlled by a dilation filter) with a residual model value, chosen among: average estimated background, median of source pixels. Residual model value can be randomized if desired. Source removal is controlled by two significance thresholds. 
Sources with fluxes above the higher threshold are removed regardless of any other conditions (e.g. on source type, etc). Sources with fluxes above the lower threshold (but below the higher threshold) are removed conditionally on chosen source type assigned in the finding process (e.g. point-like, compact, extended). Sources tagged as point-like can be removed with two different algorithms. The first one is described above and consists of replacing source pixel values by model values. The second method uses source fit model (if available) and subtract flux model from the input image. Removal of sources with nested components is controlled by the ``removeNestedSources`` flag. If enabled, the removal/subtraction process is done
on nested sources and not on parent source pixels. On the contrary, sources are removed as described above and nested sources are removed, being part of the parent.

+----------------------------+----------------------------------+-----------+------------------------+
|       Option               |             Description          |  Default  |   Values               |
+============================+==================================+===========+========================+
| ``residualZHighThr``       | | High source significance       |    10     |                        |
|                            | | threshold (in nsigmas wrt bkg) |           |                        |
|                            | | used to remove sources         |           |                        |
+----------------------------+----------------------------------+-----------+------------------------+
| ``residualZThr``           | | Source significance            |     5     |                        |
|                            | | threshold (in nsigmas wrt bkg) |           |                        |
|                            | | used to remove sources         |           |                        |
+----------------------------+----------------------------------+-----------+------------------------+
| ``removeNestedSources``    | | Remove nested sources instead  |   true    | | true                 |
|                            | | of parent source               |           | | false                |
|                            | | is not supported by Minuit     |           |                        |
|                            | | minimizer                      |           |                        |
+----------------------------+----------------------------------+-----------+------------------------+
| ``dilateKernelSize``       | | Dilation filter kernel size in |     9     |                        |
|                            | | pixels used to remove sources. |           |                        |
|                            | | NB: Must be an odd number >1   |           |                        |
|                            | | This option controls the halo  |           |                        |
|                            | | size around source to be       |           |                        |
|                            | | removed                        |           |                        |
+----------------------------+----------------------------------+-----------+------------------------+
| ``removedSourceType``      | | Type of sources to be removed  |     2     | | -1=all types         |
|                            | | threshold (in nsigmas wrt bkg) |           | | 1=compact            |
|                            | | used to remove sources         |           | | 2=point-like         |
|                            |                                  |           | | 3=extended           |
+----------------------------+----------------------------------+-----------+------------------------+
| ``residualModel``          | | Residual model used to replace |     1     | | 1=bkg                |
|                            | | source pixel values            |           | | 2=source median      |
+----------------------------+----------------------------------+-----------+------------------------+
| ``residualModelRandomize`` | | Randomize residual model pixel |   false   | | true                 |
|                            | | values                         |           | | false                |
+----------------------------+----------------------------------+-----------+------------------------+
| ``psSubtractionMethod``    | | Method used to subtract point  |     1     | | 1=model removal      |
|                            | | sources                        |           | | 2=fit model subtract |
+----------------------------+----------------------------------+-----------+------------------------+


-------------------------------
Extended Source Finding Options
-------------------------------

These options enable control of extended source search. Specific options for the available algorithms are 
reported in the Tables below. Superpixel Hierarchical Clustering algorithm is not currently available (not ported yet 
from CAESAR old repository).

+---------------------------------+----------------------------------+-----------+------------------------+
|       Option                    |             Description          |  Default  |   Values               |
+=================================+==================================+===========+========================+
| ``searchExtendedSources``       | | Enable/disable search of       |   false   | | true                 |
|                                 | | extended sources after compact |           | | false                |
|                                 | | source finding                 |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``extendedSearchMethod``        | | Extended source search method  |     4     | | 1=Wavelet Transform  |
|                                 | |                                |           | | 2=SP Hier Clustering |
|                                 | |                                |           | | 3=Active Contour     |
|                                 | |                                |           | | 4=Saliency Filter    |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``useResidualInExtendedSearch`` | | Use residual map as input for  |   true    | | true                 |
|                                 | | extended source search         |           | | false                |
|                                 | | source finding                 |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``usePreSmoothing``             | | Apply smoothing to residual    |   true    | | true                 |
|                                 | | map before performing extended |           | | false                |
|                                 | | source finding                 |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``smoothFilter``                | | Filter used to smooth residual |    2      | | 1=gaus               |
|                                 | | map                            |           | | 2=guided             |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``gausFilterKernSize``          | | Gaussian filter kernel size    |    5      |                        |
|                                 | | in pixels. NB: Must be an odd  |           |                        |
|                                 | | value                          |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``gausFilterSigma``             | | Gaussian filter sigma par      |    1      |                        |
|                                 | | in pixels                      |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``guidedFilterRadius``          | | Guided filter radius par       |    12     |                        |
|                                 | | in pixels                      |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``guidedFilterColorEps``        | | Guided filter epsilon par      |   0.04    |                        |
|                                 | | (regularization parameter)     |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+

		
		
		
-----------------------------------
Wavelet Transform Algorithm Options
-----------------------------------

These options enable control of extended source search with the Wavelet Transform method. 

+-----------------------+----------------------------------+-----------+------------------+
|       Option          |           Description            |  Default  |   Values         |
+=======================+==================================+===========+==================+
| ``wtScaleSearchMin``  | | Minimum Wavelet scale to be    |    3      |                  |
|                       | | used for extended source       |           |                  |
|                       | | search                         |           |                  |
+-----------------------+----------------------------------+-----------+------------------+
| ``wtScaleSearchMax``  | | Maximum Wavelet scale to be    |    6      |                  |
|                       | | used for extended source       |           |                  |
|                       | | search                         |           |                  |
+-----------------------+----------------------------------+-----------+------------------+

	
	
--------------------------------
Active Contour Algorithm Options
--------------------------------

These options enable control of extended source search with the Active Contour method. Two algorithms are 
provided: Chan-Vese, Linear Region-based Active Contour (LRAC).

+---------------------------------+----------------------------------+-----------+------------------------+
|       Option                    |             Description          |  Default  |   Values               |
+=================================+==================================+===========+========================+
| ``acMethod``                    | | Active contour method          |    1      | | 1=Chan-Vese          |
|                                 | |                                |           | | 2=LRAC               |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``acNIters``                    | | Maximum number of iterations   |   1000    |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``acInitLevelSetMethod``        | | Level set initialization       |     1     | | 1=circle             |
|                                 | | method                         |           | | 2=checkerboard       |
|                                 |                                  |           | | 3=saliency           |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``acInitLevelSetSizePar``       | | Level set size fraction wrt    |    0.1    |                        |
|                                 | | to minimum image size (e.g.    |           |                        |
|                                 | | circle radius=fraction x image |           |                        |
|                                 | | size)                          |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``acTolerance``                 | | Tolerance parameter to stop    |    0.1    | | 0                    |
|                                 | | main iteration loop            |           | | 1                    |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``cvNItersInner``               | | Number of iteration done in    |     5     |                        |
|                                 | | inner cycle in Chan-Vese algo  |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``cvNItersReInit``              | | Number of iteration done in    |     5     |                        |
|                                 | | re-initialization step in      |           |                        |
|                                 | | Chan-Vese algo                 |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``cvTimeStepPar``               | | Chan-Vese time step par        |   0.007   |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``cvWindowSizePar``             | | Chan-Vese window size par      |     1     |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``cvLambda1Par``                | | Chan-Vese lambda1 par          |     1     |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``cvLambda2Par``                | | Chan-Vese lambda2 par          |     2     |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``cvMuPar``                     | | Chan-Vese mu par               |    0.5    |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``cvNuPar``                     | | Chan-Vese nu par               |     0     |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``cvPPar``                      | | Chan-Vese p par                |     1     |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``lracLambdaPar``               | | LRAC regularization par        |    0.1    |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``lracRadiusPar``               | | LRAC radius of locatization    |    10     |                        |
|                                 | | ball par                       |           |                        |
+---------------------------------+----------------------------------+-----------+------------------------+
| ``lracEpsPar``                  | | LRAC convergence par           |   0.01    |                        |
+---------------------------------+----------------------------------+-----------+------------------------+

	

------------------------------------
Saliency Filtering Algorithm Options
------------------------------------

These options enable control of extended source search with the Saliency Filtering method. 

+------------------------------------+-------------------------------+-----------+------------------------+
|            Option                  |         Description           |  Default  |   Values               |
+====================================+===============================+===========+========================+
| ``spBeta``                         | | Superpixel regularization   |    1      |                        |
|                                    | | parameter                   |           |                        |
+------------------------------------+-------------------------------+-----------+------------------------+
| ``spMinArea``                      | | Superpixel min area         |    10     |                        |
|                                    | | parameter in pixels         |           |                        |
+------------------------------------+-------------------------------+-----------+------------------------+
| ``saliencyResoMin``                | | Superpixel min scale par in |    20     |                        |
|                                    | | pixels used in multi-scale  |           |                        |
|                                    | | saliency calculation        |           |                        |
+------------------------------------+-------------------------------+-----------+------------------------+
| ``saliencyResoMax``                | | Superpixel max scale par in |    60     |                        |
|                                    | | pixels used in multi-scale  |           |                        |
|                                    | | saliency calculation        |           |                        |
+------------------------------------+-------------------------------+-----------+------------------------+
| ``saliencyResoStep``               | | Superpixel scale step par in|    10     |                        |
|                                    | | pixels used in multi-scale  |           |                        |
|                                    | | saliency calculation        |           |                        |
+------------------------------------+-------------------------------+-----------+------------------------+
| ``saliencyNNFactor``               | | Fraction of most similar    |    0.2    | | 0                    |
|                                    | | superpixel neighbors used   |           | | 1                    |
|                                    | | in saliency map computation |           |                        |
+------------------------------------+-------------------------------+-----------+------------------------+
| ``saliencyUseRobustPars``          | | Use robust stats pars in    |   false   | | true                 |
|                                    | | saliency map computation    |           | | false                |
+------------------------------------+-------------------------------+-----------+------------------------+
| ``saliencyDissExpFalloffPar``      | | Superpixel dissimilarity    |    100    |                        |
|                                    | | exponential decay parameter |           |                        |
|                                    | | used in saliency map        |           |                        |
|                                    | | computation                 |           |                        |
+------------------------------------+-------------------------------+-----------+------------------------+
| ``saliencySpatialDistRegPar``      | | Regularization parameter    |     1     |                        |
|                                    | | controlling superpixel      |           |                        |
|                                    | | spatial-intensity balance in|           |                        |
|                                    | | in distance measure used for|           |                        |
|                                    | | saliency map computation    |           |                        |
|                                    | | (1 means equal weights)     |           |                        |
+------------------------------------+-------------------------------+-----------+------------------------+
| ``saliencyMultiResoCombThrFactor`` | | Fraction of resolution      |    0.7    | | 0                    |
|                                    | | scales required             |           | | 1                    |
|                                    | | above threshold to          |           |                        |
|                                    | | consider a pixel salient.   |           |                        |
|                                    | | If set to 1 a pixel is      |           |                        |
|                                    | | considered salient if its   |           |                        |
|                                    | | saliency value at all       |           |                        |
|                                    | | scales is above threshold   |           |                        |
+------------------------------------+-------------------------------+-----------+------------------------+
| ``saliencyUseBkgMap``              | | Add background map to       |   false   | | true                 |
|                                    | | total saliency map          |           | | false                |
+------------------------------------+-------------------------------+-----------+------------------------+
| ``saliencyUseNoiseMap``            | | Add noise map to            |   false   | | true                 |
|                                    | | total saliency map          |           | | false                |
+------------------------------------+-------------------------------+-----------+------------------------+
| ``saliencyThrFactor``              | | Saliency threshold factor   |    2.8    |                        |
|                                    | | parameter. Threshold is     |           |                        |
|                                    | | computed as:                |           |                        |
|                                    | | thr=<saliency>*factor       |           |                        |
|                                    | | (<saliency> is the median)  |           |                        |
|                                    | | if ``saliencyUseOptimalThr``|           |                        |
|                                    | | disabled                    |           |                        |
+------------------------------------+-------------------------------+-----------+------------------------+
| ``saliencyUseOptimalThr``          | | Use optimal threshold in    |   true    | | true                 |
|                                    | | multiscale saliency         |           | | false                |
|                                    | | thresholding. If true the   |           |                        |
|                                    | | threshold is computed as    |           |                        |
|                                    | | max(min(otsuThr,valleyThr), |           |                        |
|                                    | | medianThr)                  |           |                        |
+------------------------------------+-------------------------------+-----------+------------------------+
| ``saliencyImgThrFactor``           | | Threshold factor on input   |    1      |                        |
|                                    | | map to consider a pixel as  |           |                        |
|                                    | | salient. Threshold is set as|           |                        |
|                                    | | thr=<img>*factor (<img> is  |           |                        |
|                                    | | the median). Pixel below    |           |                        |
|                                    | | threshold are not set as    |           |                        |
|                                    | | salient even if saliency is |           |                        |
|                                    | | above saliency threshold    |           |                        |
+------------------------------------+-------------------------------+-----------+------------------------+


