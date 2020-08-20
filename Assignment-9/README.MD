A flow cytometer detects the fluorescence signal from cells that are flowing past a laser-detector. The data is arriving at regular intervals and the task is to identify 

the location (record number) of the peak of the fluorescence and the width (full-width-half-max) of each cell
 the mean arrival time between cells (how does it match with the median arrival time)
the mean width of the fluorescence peaks
Since the data file is very large, and the objective is to learn to do realtime analysis, do not read the full file. Instead, read one point at at time, process it, and make a measurement. Sample code attached. Your output file should only be the locationof the peak, and the FWHM associated with that peak. Here is the data

https://owncloud.iitm.ac.in/index.php/s/OMo0zj6k583s4Uy

The mean and median times between cells, and the mean FWHM, can also be estimated in real time. These numbers will be inaccurate to begin with, but will improve (and stabilise) as you read more data from the file. You can log the mean and median times in a separate output file and observe them converge.
