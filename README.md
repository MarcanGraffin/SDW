Here you will find documentation and codes relative to satellite-dervied waterlines (SDW) extraction and analysis.

# SDW

Many methods enable the extraction of waterline from optical satellite imagery, usually provided by the cloud-platform Google Earth Engine (GEE). 
Here are listed some of them :
   - CoastSat : https://www.sciencedirect.com/science/article/pii/S1364815219300490
   - CASSIE : https://www.sciencedirect.com/science/article/abs/pii/S1364815221000761
   - SHOREX : https://www.sciencedirect.com/science/article/pii/S0378383918306070
   - ShorelineMonitor (see https://www.nature.com/articles/s41598-018-24630-6)
   - Shoreliner (in review for _Remote Sensing_)

These algorithms share common characteristics and generally follow a similar operational framework :
![image](https://github.com/MarcanGraffin/SDW/assets/148250755/743c96af-0a41-4ab4-8ed2-f4ac9aba5b72)


# Data extraction

For shoreline extraction we use a Python-GEE module inspired by Shoreliner algorithm (Bergsma et al., 2023).

In a few words, the method uses a water index (can be defined by the user, SCoWI index by default) combined with a thresholding method (refined Otsu) to extract waterlines, between 1984-now, from publicly available satellite imagery provided by Sentinel-2 and Landsat (5, 7, 8 and 9) missions. 

A Python-GEE module has been developed (largely inspired by the CoastSat algorithm, see Vos et al., 2019) to generate and save water-indexed images from the Google Earth Engine servers.
