
Configuration Options 
=====================

Most of the example applications provided can be configured from command line arguments, as described in the Tutorial section.
Some applications, like source finding, however, require a large set of configuration options, specified inside a configuration file, passed to the application
as a command line argument ``--config=[FILE]``.

In this section we report a list of the configuration options defined in CAESAR to customize tasks.

-------------
Input Options
-------------

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

+--------------------------------+----------------------------------+-----------------------+-------------+
|       Option                   |             Description          |      Default          |   Values    |
+================================+==================================+=======================+=============+
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
| ``mergeSourcesAtEdge``         | | Merge sources found at tile    |         true          | | true      |
|                                | | edge by each worker            |                       | | false     |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``mergeSources``               | | Merge overlapping sources      |         false         | | true      |
|                                | | found by workers.              |                       | | false     |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``mergeCompactSources``        | | Merge overlapping compact      |         false         | | true      |
|                                | | sources found by workers       |                       | | false     |
+--------------------------------+----------------------------------+-----------------------+-------------+
| ``mergeExtendedSources``       | | Merge overlapping compact-     |         false         | | true      |
|                                | | extended sources found by      |                       | | false     |
|                                | | workers                        |                       |             |
+--------------------------------+----------------------------------+-----------------------+-------------+



