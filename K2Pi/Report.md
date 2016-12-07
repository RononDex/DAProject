# K2pi

## Introduction

This is the report for the K2pi assignement for the Datenanalyse project.
In this report we will document our results and how we got to said results.

The project is a collaboration effort between Nora Salgo, Tino Heuberger, Manuel Sommerhalder, Stefan Hochrein and Anar Bold

### The experiment

The experiment, which is called K2Pi, analyses the decay of a K<sup>+</sup> particle into a π<sup>0</sup> and a π<sup>+</sup>
Where the lifetime of a kaon can be determined by analasying the data from anohter experiment called LifeK.

![Experiment setup](https://raw.githubusercontent.com/RononDex/DAProject/master/K2Pi/ExperimentSetup.jpg)

As we can see, we have a collimator which accelerates the K<sup>+</sup> and generates a stream of particles. Then we have two detectors placed at a certain distance from each other along the flight path of the kaon. The first detector is called the "Up stream detector" and measures the direction and the momentum of the Kaons from the stream: The second detector, called "Downstream detector", which is composed of a tracking detector to detect the π<sup>+</sup> and a calorimeter to measure the π<sup>0</sup>.

However the pions also have a decay time which is given by τ<sub>π</sub> = 2.6 * 10<sup>-8</sup>s, this gives us an average travel distance of Beta\*Gamma\*c\*τ<sub>π</sub> = 4.188km before a pion decays.

For the purpose of the project we were told to assume that the π<sup>0</sup> is stable and does not decay.
We now have to create a simulation to simulate this decay, and find the optimal distance between the two sensors to capture as many events as possible. In the first part of the simulation we will assume a stream that has no spread and moves exactly on the z-axis. For the second part we will add a spread to the stream coming from the colimator.

## Results

### Determining τ<sub>K<sup>+</sup></sub> from the given data
First part of our assignement was to determine the average lifetime of a K<sup>+</sup> particle. For this we have the file dec_lengths.dat that contains 100'000 measurements of K<sup>+</sup> particles and their travel time before they decay. Each row in that documents represents one measurements with the unit [meters].
To create a best fit of the data, we created tauestimator.py which plots different functions against the data to find the best fit.
![Fitting the measurement from the experiment KLife](https://raw.githubusercontent.com/abold/DAProject/master/K2Pi/plot2.jpg)
As we can see here, the function called "f1" is the best fit for the dataset. We can then use this function to determine the decay time τ<sub>K<sup>+</sup></sub>. 

### Simulating the K<sup>+</sup> decay
To simulate the decay, we wrote a numerical simulation in python. Bascially we simulate the collimator shooting one K<sup>+</sup> particle at a time and with τ<sub>K<sup>+</sup></sub> that we got from the first part of the assignement we can calculate when this particle should decay. When the kaon decays we generate the π<sup>0</sup> and the π<sup>+</sup> particle with random angles θ and φ in the K<sup>+</sup> restframe and the calculate the repective angles in the lab frame (using a Lorentz boost, matrix multiplication). Since the detectors are given in the instructions to be circular with 4m in diameter, we can now easily calculate where the particles will be at z=[location of sensor]. If this location is within the 4m of the sensors diameter we have a positive event registration.

This is how the decay looks like from the two different lab frames:
![Decay in K<sup>+</sup> and lab frame](https://raw.githubusercontent.com/RononDex/DAProject/master/K2Pi/Frames.png)
