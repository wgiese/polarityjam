# Polarityjam-app

Polarityjam is a web app to study endothelial cell behaviour and more.
First, most important features from at least one input image are extracted via the [feature extraction pipeline](polarityjam/README.md) 
also provided within this repository.
The result of the extraction process is a csv file which can then be analysed in the browser with the help
of the polarityjam-app. The app is based on R-shiny and can be used without installation by
 visiting the following link: [www.polarityjam.com](www.polarityjam.com). If you rather want to work on your own system follow the 
 installation steps below.

For instructions how to install and use the feature extraction pipeline see [this file](polarityjam/README.md).


## Polarityjam-app installation

Make sure you have [anaconda](https://www.anaconda.com/) installed. 
Clone the repository and then execute the following steps on the commandline:

```commandline
cd polarityjam
conda env create -f polarityjam.yml
conda activate polarityjam
Rscript app/app.R
```

Open the browser in the URL given in the output of the R-shiny call. (Usually smth. like `http://127.0.0.1:6660`)
