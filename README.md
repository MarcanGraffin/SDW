![GIF_low](https://github.com/MarcanGraffin/SDW/assets/148250755/8e7ae46e-00e3-403b-afce-110d8269253f)


Here you will find documentation and codes relative to satellite-dervied waterlines (SDW) extraction and analysis.

# SDW

Many methods enable the extraction of waterlines from optical satellite imagery, usually provided by the cloud-platform Google Earth Engine (GEE, https://earthengine.google.com/). 
Here are listed some of them :
   - CoastSat : https://www.sciencedirect.com/science/article/pii/S1364815219300490
   - CASSIE : https://www.sciencedirect.com/science/article/abs/pii/S1364815221000761
   - SHOREX : https://www.sciencedirect.com/science/article/pii/S0378383918306070
   - ShorelineMonitor (see https://www.nature.com/articles/s41598-018-24630-6)
   - Shoreliner (by _Bergsma and colleagues_, in progress)

These algorithms share common characteristics and generally follow a similar operational framework :
![image](https://github.com/MarcanGraffin/SDW/assets/148250755/aaa40f5f-3954-4ae2-b751-a46e6789ac2b)

   **I. Pre-processing** <br />
   Optical satellite imagery quality, when being used for land and ocean observation, is greatly impacted by cloud coverage. So, a basic first step when generating a collection of satellite images is to apply a cloud mask. Some methods also allow to artificially improve the pixel resolution of some images (e.g. pansharpening : using the 15 m-resolution panchromatic band of Landsat (7, 8 and 9) to refine the other bands).

   **II. Band combination** <br />
   Key step for all satellite-based applications using optical images. Band combination involves the mathematical combination of different spectral bands in multispectral imagery (input) to create a new composite band (output). It is used for enhancing specific features or properties on an image, such as highlighting land cover, water bodies, geological features, and more.
   For coastal applications, usual indexes used are the (Modified) Normalized Difference Water Index ((M)NDWI), the Automated Water Extraction Index (AWEI) and the Subtractive Coastal Water Index (SCoWI). <br />
   <br />
   Here is shown the band combination for a Sentinel-2 multispectral image of Saint-Louis, Senegal. On the left, the RGB image, on the right, the SCoWI image.

   ![image](https://github.com/MarcanGraffin/SDW/assets/148250755/712187a2-414b-4ea1-827f-12429958f413)

   $SCoWI = B + 2(G - NIR) - 0.75 \times SWIR1 - 0.5 \times SWIR2$, <br />
   with B, G, NIR, SWIR1 and SWIR2 being respectively the blue, green, near-infrared, short-wave infrared 1 and 2 bands.

   **III. Contouring** <br />
   Contouring involves the extraction of an interface line from the composite image generated at the previous step. Two kind of approaches can be used : a thresholding method (generally an optimized Otsu's method) followed by a Marching Square calculation, or a edge-based segmentation (e.g. Canny Filter). This step is essential to extract the waterline from the composite image, but shows similar results regardless of the contouring process followed as it mostly relies on the computation of the composite image, meaning that if the composite image dissociates clearly the 2 classes of pixel (land/sea in our case) the contouring method will roughly extract the same interface.
   Additionnal (optional) contouring steps can be performed to refine the extracted interface, such as the Minimal Shoreline Variability (MSV, Almar et al., 2012).

   **IV. Post-processing** <br />
   Waterlines are usually projected over transects perpendicular to the coast in order to obtain cross-shore waterline position data over a range of spatial points.
   Once this projection is done, the data can be cleaned and formated into a desired format to make them more presentable. Outlier correction is part of this process, indeed raw data usually contain some outliers (due to unmasked clouds, bad weather, issues in the georeferencing, etc...). 
   As outliers can greatly bias any analysis performed on the data, it is recommanded to remove them. Note that the outlier correction used depend on the need (i.e. on the raw data).
   - The most common outlier-correction is the IQR (InterQuartile Range) correction, which defines two outlier thresholds (min and max values) based on the first and third quartile values computed from the raw data.
   - For high numbers of systematic outliers, a modal approach can be considered.
   
   Smoothing/interpolation can also be used regarding some specific need (generate monthly aggregates for example).
