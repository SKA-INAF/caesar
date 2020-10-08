
Data Products 
=============

Several data outputs are produced by \caesar{} at the end of the processing if corresponding options are activated in the configuration file:

-----------
ROOT Output
-----------

If `saveToFile` option is enabled, a ROOT file is produced with different stored objects (depending on the activated options):

* ROOT TTree (named `SourceInfo`) with caesar Source objects, containing summary parameters plus detailed information at pixel level (see Source API section);   
* ROOT TTree (named `ConfigInfo`) with a list of configuration options used in the run;   
* ROOT TTree (named `PerformanceInfo`) with a list of run parameters (runtimes at different stages, used memory, etc);   
     
    - Input map, stored as a caesar Image object (see Image API section)   
    - Background, noise and significance maps, stored as caesar Image objects (see Image API section)  
    - Residual map, stored as a caesar Image object (see Image API section)   
    - Segmentation and saliency maps, stored as a caesar Image object (see Image API section)   

----------
DS9 Output
----------

If `saveDS9Region` option is enabled, two ds9 region files are produced with the list of catalogued source islands and fitted components, respectively reported as labeled polygon or ellipse regions, each with a series of tags assigned (e.g. compact/point-like vs extended, real vs spurious, etc).

------------
Ascii Output
------------

If `saveToCatalogFile` option is enabled, two tabular ascii files are produced with a series of summary parameters for each catalogued source islands and fitted components, respectively.  

Table format for source islands is described below:    

+---------+-------------------------+---------+--------------------------------------------------+
| Col No. |         Name            |   Unit  |             Description                          |
+=========+=========================+=========+==================================================+
|  1      | ``name``                |         | Source name assigned by finder                   |
+---------+-------------------------+---------+--------------------------------------------------+
|  2      | ``iauName``             |         | Source name in IAU notation                      |
+---------+-------------------------+---------+--------------------------------------------------+
|  3      | ``nPix``                |         | Number of pixels in island                       |
+---------+-------------------------+---------+--------------------------------------------------+
|  4      | ``nComp``               |         | Number of fitted components                      |
|         |                         |         | (=0 if fit not performed or failed)              |
+---------+-------------------------+---------+--------------------------------------------------+
|  5      | ``nNested``             |         | Number of nested sources found in island         |
+---------+-------------------------+---------+--------------------------------------------------+
|  6-7    | ``(x,y)``               |         | Source x,y position                              |
+---------+-------------------------+---------+--------------------------------------------------+
|  8-9    | ``(xw,yw)``             |         | Source x,y position weighted by pixel fluxes     |
+---------+-------------------------+---------+--------------------------------------------------+
|  10-11  | ``(x_wcs,y_wcs)``       | deg     | Source x,y sky position                          |
+---------+-------------------------+---------+--------------------------------------------------+
|  12-13  | ``(xw_wcs,yw_wcs)``     | deg     | Source x,y sky position weighted by pixel fluxes |
+---------+-------------------------+---------+--------------------------------------------------+
|  14-15  | ``(xmin,xmax)``         |         | Source pixel min, max x coord.                   |
+---------+-------------------------+---------+--------------------------------------------------+
|  16-17  | ``(ymin,ymax)``         |         | Source pixel min, max y coord.                   |
+---------+-------------------------+---------+--------------------------------------------------+
|  18-19  | ``(xmin_wcs,xmax_wcs)`` | deg     | Source pixel min, max x sky coord.               |
+---------+-------------------------+---------+--------------------------------------------------+
|  20-21  | ``(ymin_wcs,ymax_wcs)`` | deg     | Source pixel min, max y sky coord.               |
+---------+-------------------------+---------+--------------------------------------------------+
|  22     | ``nu``                  | GHz     | Frequency value extracted from image header      |
+---------+-------------------------+---------+--------------------------------------------------+
|  23     | ``Ssum``                | Jy/beam | Sum of island pixel brightness                   |
+---------+-------------------------+---------+--------------------------------------------------+
|  24     | ``Smax``                | Jy/beam | Max pixel brightness in island                   |
+---------+-------------------------+---------+--------------------------------------------------+
|  25-26  | ``(S,Serr)``            | Jy/beam | Source fitted flux brightness (not               |
|         |                         |         | corrected by beam area) and its error            |
+---------+-------------------------+---------+--------------------------------------------------+
|  27     | ``npix_beam``           |         | Number of pixels in beam                         |
+---------+-------------------------+---------+--------------------------------------------------+
|  28     | ``bkg_sum``             | Jy/beam | Background summed over island pixels             |
+---------+-------------------------+---------+--------------------------------------------------+
|  29     | ``rms_sum``             | Jy/beam | Background noise summed over island pixels       |
+---------+-------------------------+---------+--------------------------------------------------+
|  30     | ``F1``                  |         | Source morph. flag {1=COMPACT,2=POINT-LIKE,      |
|         |                         |         | 3=EXTENDED,4=COMPACT-EXTENDED}                   |
+---------+-------------------------+---------+--------------------------------------------------+
|  31     | ``F2``                  |         | Sourceness flag{1=REAL,2=CANDIDATE,3=FAKE}       |
+---------+-------------------------+---------+--------------------------------------------------+
|  32     | ``F3``                  |         | Boolean flag indicating if source was tagged as  |
|         |                         |         | good (true) or bad (false) in finding process    |
+---------+-------------------------+---------+--------------------------------------------------+
|  33     | ``F4``                  |         | Source depth level flag                          |
|         |                         |         | (nested source if >0)                            |
+---------+-------------------------+---------+--------------------------------------------------+


Table format for source components is described below:    

+---------+---------------------------------------------------+----------+--------------------------------------+
| Col No. |         Name                                      |   Unit   |             Description              |
+=========+===================================================+==========+======================================+
|  1      | ``name``                                          |          | Source name assigned by finder       |
+---------+---------------------------------------------------+----------+--------------------------------------+
|  2      | ``nPix``                                          |          | Number of pixels in island           |
+---------+---------------------------------------------------+----------+--------------------------------------+
|  3      | ``compId``                                        |          | Component id                         |
+---------+---------------------------------------------------+----------+--------------------------------------+
|  4      | ``iauName``                                       |          | Source name in IAU notation          |
+---------+---------------------------------------------------+----------+--------------------------------------+
|  5-6    | ``(x,y)``                                         |          | Component x,y position               |
+---------+---------------------------------------------------+----------+--------------------------------------+
|  7-8    | ``(xerr,yerr)``                                   |          | Component x,y position errors        |
+---------+---------------------------------------------------+----------+--------------------------------------+
|  9-10   | ``(x_wcs,y_wcs)``                                 | deg      | Component x,y sky position           |
+---------+---------------------------------------------------+----------+--------------------------------------+
|  11-12  | ``(xerr_wcs,yerr_wcs)``                           | deg      | Component x,y sky position errors    |
+---------+---------------------------------------------------+----------+--------------------------------------+
|  13     |  ``nu``                                           | GHz      | Frequency value extracted from       |
|         |                                                   |          | image header                         |
+---------+---------------------------------------------------+----------+--------------------------------------+
|  14-15  | ``(Speak,Speak_err)``                             | Jy/beam  | Component peak brightness and error  |
+---------+---------------------------------------------------+----------+--------------------------------------+
|  16-17  | ``(S,Serr)``                                      | Jy/beam  | Component brightness (not corrected  |
|         |                                                   |          | by beam area) and error              | 
+---------+---------------------------------------------------+----------+--------------------------------------+
|  18-19  | ``(S_i,Serr_i)``                                  | Jy/beam  | Island brightness (not corrected by  |
|         |                                                   |          | beam area) and error                 |  
+---------+---------------------------------------------------+----------+--------------------------------------+
|  20     | ``npix_beam``                                     |          | Number of pixels in beam             |
+---------+---------------------------------------------------+----------+--------------------------------------+
|  21-23  | ``(a,b,theta)``                                   |          | Component ellipse pars               |
+---------+---------------------------------------------------+----------+--------------------------------------+
|  24-26  | ``(a_err,b_err,theta_err)``                       |          | Component ellipse par errors         |
+---------+---------------------------------------------------+----------+--------------------------------------+
|  27-29  | ``(a_wcs,b_wcs,theta_wcs)``                       | ",",deg  | Component ellipse pars in sky coords |
+---------+---------------------------------------------------+----------+--------------------------------------+
|  30-32  | ``(a_err_wcs,b_err_wcs,theta_err_wcs)``           | ",",deg  | Component ellipse par                |
|         |                                                   |          | errors in sky coords                 |
+---------+---------------------------------------------------+----------+--------------------------------------+
|  33-35  | ``(bmaj,bmin,pa)``                                | ",",deg  | Beam ellipse pars in sky coords      |
+---------+---------------------------------------------------+----------+--------------------------------------+
|  36-38  | ``(a_deconv_wcs,b_deconv_wcs,theta_deconv_wcs)``  | ",",deg  | Component ellipse par                |
|         |                                                   |          | in sky coords, deconvolved by beam   |
+---------+---------------------------------------------------+----------+--------------------------------------+
|  39     | ``E/E_beam``                                      |          | Ratio of fitted and beam ellipse     |
|         |                                                   |          | eccentricities                       |
+---------+---------------------------------------------------+----------+--------------------------------------+
|  40     | ``A/A_beam``                                      |          | Ratio of fitted and beam ellipse     |
|         |                                                   |          | areas                                |
+---------+---------------------------------------------------+----------+--------------------------------------+
|  41     | ``psi``                                           | deg      | Rotation angle between fit and beam  |
|         |                                                   |          | ellipse                              |
+---------+---------------------------------------------------+----------+--------------------------------------+
|  42     | ``bkg_sum``                                       | Jy/beam  | Background summed over island pixels |
+---------+---------------------------------------------------+----------+--------------------------------------+
|  43     | ``rms_sum``                                       | Jy/beam  | Background noise summed over island  |
|         |                                                   |          | pixels                               |
+---------+---------------------------------------------------+----------+--------------------------------------+
|  44-45  | ``(chi2,ndf)``                                    |          | Source island fit chisquare and      |
|         |                                                   |          | degrees of freedom                   |
+---------+---------------------------------------------------+----------+--------------------------------------+
|  46     | ``F1``                                            |          | Source fit quality flag              |
|         |                                                   |          | {0=BAD,1=LOW,2=MEDIUM,3=HIGH}        |
+---------+---------------------------------------------------+----------+--------------------------------------+
|  47     | ``F2``                                            |          | Component sourceness flag            |
|         |                                                   |          | {1=REAL,2=CANDIDATE,3=FAKE}          |
+---------+---------------------------------------------------+----------+--------------------------------------+
|  48     | ``F3``                                            |          | Component morph. flag                |
|         |                                                   |          | {1=COMPACT,2=POINT-LIKE,3=EXTENDED,  |
|         |                                                   |          | 4=COMPACT-EXTENDED}                  |
+---------+---------------------------------------------------+----------+--------------------------------------+  


-----------
FITS Output
-----------

If `saveToFITSFile` option is enabled, a series of FITS image files are produced (if corresponding options are activated):   

* Input map     
* Background, noise and significance maps   
* Residual map   
* Saliency map




