# MindSpark
### Muse Device-based Matlab application for stress relief.
![MindSpark Logo](https://github.com/macdvdup/MindSpark/blob/main/interface/logo2.png?raw=true)

With this work, we porpused a MatLab application for meditation subject using the Muse Wearable device. 

## How to test the app

In order to test the application, there are some preliminar steps to be take, specifically to install the uvicMuse app and dependencies, as well as the MatLab toolbox that comes with. A detail tutorial how to these tools can be found in https://github.com/krigolson/uvicMuse.
It may also be necessary to install aditional MatLab add-ons.

Having the environment ready, it is necessary to first connect the Muse to the uvicMuse app. Then, our app can be launched by running the 'inicio.mlapp' file inside the 'Interface' folder.

## Other contents in the repository

The repository contains other files used along the project development.
* GAMEEMO: dataset used to train the Machine Learing model
* functions: EEGLab functions that are used in the MindSpark program
* helper: scripts used to test the the functions and EEG data acquire to visualize those functions
* creating dataset: scripts, features data and KNN model used
* loose files: the functions used to determine the necessary features from the EEG signal
