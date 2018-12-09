
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
| ``saveInputMap``               | | Save input map to output       |        true           | | true      |
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
| ``nestedBlobKernFactor``              | | Filter kernel size factor par  |    1      |                        |
|                                       | | so that kern size=             |           |                        |
|                                       | | factor x sigma (sigma is the   |           |                        |
|                                       | | filter scale par in pixels)    |           |                        |
+---------------------------------------+----------------------------------+-----------+------------------------+
| ``sourceToBeamAreaThrToSearchNested`` | | Mother source area/beam thr to |    0      |                        |
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
|                                       | | (NBeams=blob npix/beam npix)   |           | | false                |
+---------------------------------------+----------------------------------+-----------+------------------------+
| ``psNBeamsThr``                       | | Max number of beams cut value. |   5       |                        |
|                                       | | Blobs with a number of beams   |           |                        |
|                                       | | below this threshold           |           |                        |
|                                       | | passed the point-like cut      |           |                        |
+---------------------------------------+----------------------------------+-----------+------------------------+

		

