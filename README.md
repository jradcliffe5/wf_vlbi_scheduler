# Wide-field VLBI scheduler

This code generates the phase centres for wide-field VLBI correlation. The final version shall have capabilities to cover an entire area but currently only takes an input catalogue and outputs a VEX suitable format for correlators. 

## Prerequisites
The code uses python 3.7 and requires the following dependencies that can be installed using `pip`
* traceback
* logging
* astropy
* numpy

## Catalogue pre-processing

Before using this code, it is advisable to do some catalogue pre-processing. In particular:

* Ensure that the prior catalogue can be readable using the astropy `Table.read()` class. The supported formats are at the bottom of this README and should be used in the `table_format` input parameter.
* The RA and Dec columns should be in <b>degrees</b>. The code doesn't have the capabilities to convert just yet.

## Usage 
1. Clone the repository using `git clone https://github.com/jradcliffe5/wf_vlbi_scheduler.git`.
2. To run the code copy the input file (`WFVLBI_inputs.txt`) to your current working directory.
3. Edit the input script.
4. Run the script using `python <path to repository>/wf_vlbi_scheduler.py` or `python3` if you have both 2.7 & 3 installed.
5. Check the output is as you expected (see the `*.vex` file).

### Supported catalogue formats


|           Format           |
|:----------------------------|
|                      ascii |        
|               ascii.aastex |          
|                ascii.basic |          
|                  ascii.cds |          
|    ascii.commented_header |          
|                  ascii.csv |                           
|              ascii.daophot |                            
|                 ascii.ecsv |                          
|          ascii.fast_basic |                           
|ascii.fast_commented_header |                           
|             ascii.fast_csv |                           
|       ascii.fast__header   |                           
|             ascii.fast_rdb |                           
|             ascii.fast_tab |                          
|          ascii.fixed_width |                          
|ascii.fixed_width__header   |                        
| ascii.fixed_width_two_line |                          
|                 ascii.html |                         
|                 ascii.ipac |                          
|                ascii.latex |                         
|            ascii._header   |                        
|                  ascii.rdb |                         
|                  ascii.rst |                          
|           ascii.sextractor |                           
|                  ascii.tab |                          
|                       fits |                         
|                       hdf5 |                         
|                    votable |                         
|                     aastex |                       
|                        cds |                        
|                        csv |                      
|                    daophot |                        
|                       html |                       
|                      ipac |                       
|                      latex |                       
|                        rdb |                       

