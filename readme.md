# Lookback Option Pricer with Monte Carlo based on Monte Carlo approach:

## Description:

The pricer returns the values of the price of the options and the greeks with the Price and Delta curves as a function of the underlying at 0. 

This repository contains the visual studio implementation of a lookback pricer in C++ and the user interface written in VBA. 

A dll library is created using C++ which is then used by the VBA code. 

Please note that this project works only on 64bits windows devices. 

## Steps: 

Before testing the program: 

Make sure to update the paths in the VBA code of the user interface.

The user interface is accessible by opening the excel in: C++\lookback\x64

Then in Excel, select Developer => Visual Basic => Right click on ModernUserForm and select code

The Lib (dll C++) directory must be changed by adding the correct path to the C++ folder

"path_to_C++_folder\C++\lookback\x64\Debug\lookback.dll"

Do not remove the images from the folder. 

## Theory: 

The report explains in details the implementation of the Price.  

## screenshots: 

The user has to specify the characteristics of the option for the price and the greeks to be returned. The use of the interface is intuitive and don't require prior knowledge of how to code in C++ or in VBA. 

Here you'll find a screenshot of user interface.

The user interface detects the compatibility of the values entrered by the user: 

When the values entered are correct, the price is deduced based on the code of the C++ library. 


