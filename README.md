HyperSeed is a MATLAB software to segment seeds from hyperspectral images.

# Introduction

High-throughput, nondestructive, and precise measurement of seeds is critical for the evaluation of seed quality and the improvement of agricultural productions. To this end, we have developed a novel end-to-end platform named HyperSeed to provide hyperspectral information for seeds. As a test case, the hyperspectral images of rice seeds are obtained from a high-performance line-scan image spectrograph covering the spectral range from 600 to 1700 nm. The acquired images are processed via a graphical user interface (GUI)-based open-source software for background removal and seed segmentation. The output is generated in the form of a hyperspectral cube and curve for each seed. In our experiment, we presented the visual results of seed segmentation on different seed species. Moreover, we conducted a classification of seeds raised in heat stress and control environments using both traditional machine learning models and neural network models. The results show that the proposed 3D convolutional neural network (3D CNN) model has the highest accuracy, which is 97.5% in seed-based classification and 94.21% in pixel-based classification, compared to 80.0% in seed-based classification and 85.67% in seed-based classification from the support vector machine (SVM) model. Moreover, our pipeline enables systematic analysis of spectral curves and identification of wavelengths of biological interest.

## workflow

![](https://raw.githubusercontent.com/tgaochn/HyperSeed/main/fig/workflow.png)

## GUI

![](https://raw.githubusercontent.com/tgaochn/HyperSeed/main/fig/GUI.png)

# Software and Data

## Open-source Software

run "hyperSeed.m"

## Dataset

https://uofnelincoln-my.sharepoint.com/:u:/g/personal/tgao6_unl_edu/EWFoM_MYgTFNr6pIOQmskJEBppRnPvUEd-dZrO-lBkDGzQ?e=MOephm

## Stand-alone Software

https://uofnelincoln-my.sharepoint.com/:u:/g/personal/tgao6_unl_edu/ERtAA4p3NW9Gksl_tpxlM2QBIUvhrXxXa1SPgSZbwKLfbA?e=xuf4Z3

# Reference

If you find our software helpful, please cite our paper:
Gao, T.; Chandran, A.K.N.; Paul, P.; Walia, H.; Yu, H. HyperSeed: An End-to-End Method to Process Hyperspectral Images of Seeds. Sensors 2021, 21, 8184. https://doi.org/10.3390/s21248184
