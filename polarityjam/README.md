# Feature extraction pipeline

The feature extraction pipeline to study vascular cell dynamics and more. 
This instruction is for installing the feature extraction algorithm only. How to 
install the app to visualize results is explained [here](https://github.com/wgiese/vascu-ec-app) or
if you cloned the whole project and now look at this instructions visit [this file](../README.md).  

## Install the feature extraction pipeline:

- Install via pip:
    - clone https://github.com/wgiese/vascu-ec-app via git or download via browser
    - change to the base directory
    - run `pip install -e polarityjam`
    
## Execute the feature extraction pipeline

- from the commandline:
    - run `polarityjam` to obtain further help on how to start a feature extraction process
    - run `polarityjam_test` to run the provided test suit and check your installation 
    
### Extract features from single tif

run `polarityjam run paramfile.yml input.tif outputpath`. 
Optionally define a `--filename_prefix` to rename the input file.

### Extract features from a folder containing tif files

run `polarityjam run-stack paramfile.yml inputpath outputpath`. 

### Extract features from a given list of conditions (folders) provided within a csv file

run `polarityjam run-key paramfile.yml inputpath inputkey.csv outputpath`. 
csv file might look like this:

|folder_name  |short_name |
|-------------|-----------|
|/set_1/ | cond_1 |
|/set_2/ | cond_2 |

Folder name contains all folders relative to the `inputpath`.
Folder structure will also be created in the `outputpath`. 
Specify a `short_name` different to the `folder_name` to rename each folder.
(e.g. folder `set_1` will be named `cond_1` in the `putputpath`)