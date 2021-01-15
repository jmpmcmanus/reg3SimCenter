# reg3SimCenter
Code for Ajay at the SimCenter

## displayXXXXXXNetCDF Jupyter Notebooks
There are three Jupyter Notebook files, one for each ADCIRC data type, for displaying data the the netCDF files created 
from the code described below:

    displayFort63NetCDF.ipynb  
    displayFort64NetCDF.ipynb  
    displaySwan63NetCDF.ipynb

This notebooks are designed to run localy on your own computer, after you have downloaded the netCDF file from 
adcirc-db.edc.renci.org. To download data from adcirc-db.edc.renci.org using the following commmand:

scp -r bhajay@adcirc-db.edc.renci.org:/projects/regionthree/ingestProcessing/ajay/data data/

## Code Directory
This code, in the code directory, is set up to run on adcirc-db.edc.renci.org at RENCI. It is designed to run inside of a 
docker container on adcirc-db.edc.renci.org named region3db_container. This container contains an account called ajay, 
which is owned by Dr. Ajay B Harish at UC Berkley SimCenter.

The home directory of the ajay account, in the region3db_container container, containes directory work (/home/ajay/work), 
which is a symbolic link to the directory /home/data/ingestProcessing/ajay, which in turn is a volume mount of the 
/projects/regionthree/ingestProcessing/ajay directory on adcirc-db.edc.renci.org. 

In the work directory there are two more directories, /home/ajay/work/code and /home/ajay/work/data. The code directory
contains the code to extract data from the database. It is set up to output data to the /home/ajay/work/data directory.

The contents of the code directory are:

    extractStormSubset.py (code to extract data from database)  
    roi.geojson (region of Interest file used in extractStormSubset.py  to extract data from database)  
    timestamps.csv (timestamps used in extractStormSubset.py  to extract data from database)  
    runExtractStormSubset.sh (shell script with example command for running extractStormSubset.py. Can be used to run cronjob)  
    shell_region3db_ajay.sh (shell script to connect to ajay account on region3db_container)  
    psql_ajay_adcirc-db.sh (shell script to connect to reg3sim database, in the region3db_container, from your bhajay account on adcirc-db.edc.renci.org)  
    psql_ajay_region3db.sh (shell script to connect to reg3sim database, in the region3db_container, from your ajay account in the region3db_container)  
    testquery.py (test program you have already used)  

To run extractStormSubset.py in a cronjob, from the bhajay account on adcirc-db.edc.renci.org, put the following line 
in crontab, using crontab -e:

43 14 15 01 * docker exec -t --user ajay region3db_container bash -c "source /home/ajay/.bashrc; conda activate adcirc; sg postgres -c '/home/ajay/work/code/runExtractStormSubset.sh" >> ~/cron.log 2>&1

Change the start date and time (43 14 15 01 ) to the time you want to run the data. The shell file runExtractStormSubset.sh 
is currently set up to process the tables:

    var_dp3r3b1c1h1l1_fort63  
    var_dp3r3b1c1h1l1_swan63  
    var_dp3r3b1c1h1l1_fort64  

The output netCDF files are written to /home/ajay/work/data.

