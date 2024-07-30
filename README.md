Manual for SRNN-SZ

## Introduction

Usage of SRNN-SZ need to deploy 2 code repositories:

* This one (the compression framework): https://github.com/Meso272/SRNN-SZ
* HAT (my customization):https://github.com/Meso272/HAT

## Installation of compression framework

### Dependencies

* cmake>=3.13
* gcc>=6.0

The following ones are not used but need to get installed before compilation (will be removed later):

* Python >= 3.6
* numpy 
* PyWavelets
* pybind 11

The following is not mandatory:

* Zstandard (https://facebook.github.io/zstd/). Not mandatory to be manually installed as Zstandard v1.4.5 is included and will be used if libzstd can not be found by
  pkg-config.

### Installation 

* mkdir build && cd build
* cmake -DCMAKE_INSTALL_PREFIX:PATH=[INSTALL_DIR] ..
* make
* make install

Then, you'll find all the executables in [INSTALL_DIR]/bin and header files in [INSTALL_DIR]/include.

## Installation and deployment of HAT

Please follow the readme here to install HAT: https://github.com/Meso272/HAT

Important: please clone the HAT repository to the following path: your_home_path/lossycompression/HAT (check https://github.com/Meso272/SRNN-SZ/blob/sr/include/QoZ/preprocessor/SRNet.hpp 

## Run SRNN-SZ

The executable name of SRNN-SZ is srnz. It is derived from QoZ so its usage is mostly similar to the QoZ command (run srnz with no arguments to check the help information) (QoZ: https://github.com/szcompressor/QoZ)

Some important points:

* use -q to specify interpolation optimization level. The default is 1, 3, or 4 recommended to test.
* use -k to pass the trained model checkpoint path.


## Train the HAT network

(This section is under construction)

The training is done within the customized HAT code (so please make it installed before you perform the training).

First, please collect the training data, and preprocess them to 2D-slice data files.

Second, write an option file for the dataset path, network configuration, and training parameters (example: https://github.com/JLiu-1/HAT/blob/main/options/train/train_HAT_SRx2_2D_assorted.yml). 

Then run the training script provided in the HAT code for generating the model checkpoint to be used in SRNN-SZ.










