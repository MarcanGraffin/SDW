Shoreliner : Satellite-derived waterline extractor by Bergsma et al., 2024

Inputs : contains the config file, open and edit it whenever you want to extract waterlines on a new site or using a new method.*
### launcher.py needs files in pickle format, with a specific format (see example) to launch.

Codes : contains scripts and functions to extract the waterline.

Projects : contains all your current or past waterline extraction projects, with codes, config files and outputs.

### You need a certified google earth engine account to use this tool. Check on https://earthengine.google.com/ to create an account, and for details.

PREWORK : 
1 - Create a conda environment with basic packages : "conda create -n shoreliner spyder numpy scipy pandas matplotlib sympy cython" (shoreliner is a proposed name for the env)
2 - Activate the environment (eg, conda activate shoreliner)
3 - Download all the packages using first "pip install gdal" then "pip install -r PATH/requirements.txt"
4 - After the account has been validated, you would need to authenticate using the command "earthengine authenticate"
  OR launch spyder (type "spyder" in the shell) and type in the terminal type "import ee" then "ee.Authenticate()".

WATERLINE EXTRACTION : 
1 - Edit the config.yaml file in ShorelinerTool/Inputs to design the extraction process as you intend it to be (site, date range, methods, ...) and save.
2 - Launch launcher.py in the ShorelinerTool folder.
3 - Go on ShorelinerTool/Projects/NAMEOFYOURPROJECT folder, you will find poly_0 to poly_n (n being the number of sub ROIs) folders, for each of them you'll have the same 3 scripts in : download.py (download images from Google Earth Engine), main_process.py (extract waterline from images and save it) and ploting.py (if you want to observe the outputs). Launch first download.py, then main_process.py (and optionally ploting.py).
 
