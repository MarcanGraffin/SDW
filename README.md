Here you will find documentation and codes relative to satellite-dervied waterlines (SDW) extraction and analysis.

# SDW

Many methods enable the extraction of waterlines from optical satellite imagery, usually provided by the cloud-platform Google Earth Engine (GEE, https://earthengine.google.com/). 
Here are listed some of them :
   - CoastSat : https://www.sciencedirect.com/science/article/pii/S1364815219300490
   - CASSIE : https://www.sciencedirect.com/science/article/abs/pii/S1364815221000761
   - SHOREX : https://www.sciencedirect.com/science/article/pii/S0378383918306070
   - ShorelineMonitor (see https://www.nature.com/articles/s41598-018-24630-6)
   - Shoreliner (in review for _Remote Sensing_)

These algorithms share common characteristics and generally follow a similar operational framework :
![image](https://github.com/MarcanGraffin/SDW/assets/148250755/aaa40f5f-3954-4ae2-b751-a46e6789ac2b)

   **I. Pre-processing** n/
   Optical satellite imagery quality, when being used for land and ocean observation, is greatly impacted by cloud coverage. So, a basic first step when generating a collection of satellite images is to apply a cloud mask.

   II. Band combination
   Key step for all satellite-based applications using optical images. Band combination involves the mathematical combination of different spectral bands in multispectral imagery (input) to create a new composite band (output). It is used for enhancing specific features or properties on an image, such as highlighting land cover, water bodies, geological features, and more.
   For coastal applications, usual indexes used are the (Modified) Normalized Difference Water Index ((M)NDWI), the Automated Water Extraction Index (AWEI) and the Subtractive Coastal Water Index (SCowI).

   III. Contouring
   Contouring involves the extraction of an interface line from the composite image generated at the previous step. Two kind of approaches can be used : a thresholding method (generally an optimized Otsu's method) followed by a Marching Square calculation, or a edge-based segmentation (e.g. Canny Filter). This step is essential to extract the waterline from the composite image, but shows similar results regardless of the contouring process followed as it mostly relies on the computation of the composite image, meaning that if the composite image dissociates clearly the 2 classes of pixel (land/sea in our case) the contouring method will roughly extract the same interface.
   Additionnal (optional) contouring steps can be performed to refine the extracted interface, such as the Minimal Shoreline Variability (MSV, Almar et al., 2012).

   IV. Post-processing
   Waterlines are usually projected over transects perpendicular to the coast in order to obtain cross-shore waterline position data over a range of spatial points.
   Once this projection is done, the data can be cleaned and formated into ss desired to make them more presentable. Outlier correction is part of this process, indeed raw data usually contain some outliers (due to unmasked clouds, bad weather, issues in the georeferencing, etc...). 
   As outliers can greatly biased any analysis performed on the data, it is recommanded them. Note that the outlier correction used depend on the need (i.e. on the raw data).
   - The most common outlier-correction is the IQR (InterQuartile Range) correction, which defines two outlier thresholds (min and max values) based on the first and third quartile values.
   - For high numbers of systematic outliers, a modal approach can be considered.
   
   Smoothing/interpolation can also be used regarding some specific need (generate monthly data for example).

# Data extraction

For shoreline extraction we use a Python-GEE module inspired by Shoreliner algorithm (Bergsma et al., 2023).

In a few words, the method uses a water index (can be defined by the user, SCoWI index by default) combined with a thresholding method (refined Otsu) to extract waterlines, between 1984-now, from publicly available satellite imagery provided by Sentinel-2 and Landsat (5, 7, 8 and 9) missions. 

A Python-GEE module has been developed (largely inspired by the CoastSat algorithm, see Vos et al., 2019) to generate and save water-indexed images from the Google Earth Engine servers.
