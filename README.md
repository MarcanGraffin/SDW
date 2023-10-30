# SDW
Codes relative to satellite-derived shoreline/waterline data extraction and analysis.

![image](https://github.com/MarcanGraffin/SDS/assets/148250755/80954324-8c85-4fe7-9f01-917ecf9f8988)

# Data extraction

For shoreline extraction we use a Python-GEE module inspired by Shoreliner algorithm (Bergsma et al., 2023).

In a few words, the method uses a water index (can be defined by the user, SCoWI index by default) combined with a thresholding method (refined Otsu) to extract waterlines, between 1984-now, from publicly available satellite imagery provided by Sentinel-2 and Landsat (5, 7, 8 and 9) missions. 

A Python-GEE module has been developed (largely inspired by the CoastSat algorithm, see Vos et al., 2019) to generate and save water-indexed images from the Google Earth Engine servers.
