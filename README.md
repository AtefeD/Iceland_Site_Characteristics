# Iceland_Site_Characteristics

This repository provides Python code for computing and mapping the site proxies and seismic site amplification predictions for Iceland based on both Icelandic and large-scale, non-Icelandic, models. 
Full detail characteristics of the Icelandic recording stations in SISZ including geological units, site amplification factors, the δS2S_s statistic along with all associated site effect proxies at each station are given as well. 

## Reference:
Darzi Atefe, Halldorsson Benedikt, Cotton Fabrice, Rahpeyma Sahar (2024) Nationwide frequency-dependent seismic site amplification models for Iceland, Soil Dynamics and Earthquake Engineering, https://doi.org/10.1016/j.soildyn.2024.108798
<br/>
If you use the codes, please cite the relevant article appropriately.

The SiteProxyIceland.py code specifically designed for Iceland but can be easily adopted to other regions. 
In case of any questions or comments, do not hesitate to contact me via  atefe@hi.is
<br/>
<br/>

Before running the code, download the following data: 
<br/>1- The 30-arcsecond gridded thickness of geomorphological sedimentary deposit layer estimated by Pelletier et al. (2016) can be downloaded from https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1304 
<br/>2- The slope-based inferred V_S30 raster file developed by Wald and Allen (2007) can be downloaded from the U.S. Geological Survey (USGS) via https://earthquake.usgs.gov/static/lfs/data/vs30/vs30.zip 
<br/>3- The slope and geology-based site amplification map used in the European Seismic Risk Model (ESRM20) are constructed using the open source ‘Exposure To Site’ tool available here:  https://gitlab.seismo.ethz.ch/efehr/esrm20_sitemodel. 
<br/>4- The datasets of slope, geological era, and inferred VS30 from slope that were used in ESRM20 are provided through web-service of EFEHR European Site Response Model Datasets Viewer (https://maps.eu-risk.eucentre.it/map/european-site-response-model-datasets/#4/53.98/4.53) as well as from https://nextcloud.gfz-potsdam.de/s/93ZR4ky8D4mDXb9
<br/>5- The proxy-based site amplification models proposed by Loviknes et al. (2024) are available at: https://doi.org/10.5281/zenodo.8072116

Change the corresponding path in the script


The section-by-section description of code is given within SiteProxyIceland.py script. 
