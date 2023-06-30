# Wide-field VLBI scheduler

[![DOI](https://zenodo.org/badge/201935735.svg)](https://zenodo.org/badge/latestdoi/201935735)

This code generates the phase centres for wide-field VLBI correlation. The code that can generate expected simulated data sets and sensitivity maps using realistic EVN primary beams has been moved to another repo. The code has been tested for the SPARCS (Njeri et al., 2023) and COSMOS projects (Radcliffe+in prep.). The final version shall have capabilities to cover an entire area but currently only takes an input catalogue and outputs a VEX format that's suitable format for correlators.



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

## Outputs
The code should produce an output `vex` file format that should be readable by most VLBI correlators. The format is of the `$SOURCE` format required. It should look similar to the following (otherwise raise an issue).

```
def ER047001;
    source_name = ER047001;
    ra = 11h56m10.6548s; dec =  62d13'54.7893"; ref_coord_frame = J2000;
enddef;
```

If you set `do_plot=True` then you can generate a plot of all the phase centres and their FoVs. The initial source positions along with reduced phase centre locations (as some could be overlapping) are also plotted. An example is shown below. 

<img src="https://raw.githubusercontent.com/jradcliffe5/wf_vlbi_scheduler/master/testing/random_catalogue_correlation_params.png" width="500">

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

