# K2pi

## Introduction

This is the report for the K2pi assignement for the Datenanalyse project.
In this report we will document our results and how we got to said results.

The project is a collaboration effort between Nora Salgo, Tino Heuberger, Manuel Sommerhalder, Stefan Hochi and Anar Bold

### The experiment

The experiment, which is called K2Pi, analyses the decay of a K<sup>+</sup> particle into a π<sup>0</sup> and a π<sup>+</sup>
Where the lifetime of a kaon can be determined by analasying the data from anohter experiment called LifeK.

[image Experiment PDF]

As we can see, we have a collimator which accelerates the K<sup>+</sup> and generates a stream of particles. Then we have two detectors placed at a certain distance from each other along the flight path of the kaon. The first detector is called the "Up stream detectors" and measures the direction and the momentum of the Kaons from the stream: The second detector, called "Downstream detector", which is composed of a tracking detector to detect the π<sup>+</sup> and a calorimeter to measure the π<sup>0</sup>.

However the pions also have a decay time which is given by τ<sub>π</sub> = 2.6 * 10<sup>-8</sup>s, which gives us an average travel distance of Beta\*Gamma\*c\*τ<sub>π</sub> = 4.188km before a pion decays.

For the purpose of the project we were told to assume that the π<sup>0</sup> is stable and does not decay.
We now have to create a simulation to simulate this decay, and find the optimal distance between the two sensors to capture as many events as possible. In the first part of the simulation we will assume a stream that has no spread and moves exactly on the z-axis. For the second part we will add a spread to the stream coming from the colimator.

## Results

### Determining τ<sub>K<sup>+</sup></sub> from the given data
First part of our assignement was to determine the average lifetime of a K<sup>+</sup> particle. For this we have the file dec_lengths.dat that contains 100'000 measurements of K<sup>+</sup> particles and their travel time before they decay. Each row in that documents represents one measurements with the unit [m] (meters).
To create a best fit of the data, we created tauestimator.py which plots different functions against the data to find the best fit.
![Fitting the measurement from the experiment KLife](https://raw.githubusercontent.com/abold/DAProject/master/K2Pi/plot2.jpg)
