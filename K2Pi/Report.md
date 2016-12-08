# K2Pi

## Introduction

This is the report for the K2Pi assignement for the data analysis project.
In this report we will document our results and how we got to said results.

The project is a collaboration effort between Nora Salgo, Tino Heuberger, Manuel Sommerhalder, Stefan Hochrein and Anar Bold

### The experiment

The experiment, which is called K2Pi, analyses the decay of positively charged kaons (K<sup>+</sup>) into one positively charged and one neutral pion (π<sup>+</sup> and a π<sup>0</sup>), where the lifetime of a K<sup>+</sup> can be determined by analasying the data from an earlier experiment called LifeK.

![Experiment setup](https://raw.githubusercontent.com/RononDex/DAProject/master/K2Pi/ExperimentSetup.jpg)

As we can see, we have a collimator which generates a beam of K<sup>+</sup> particles. Furthermore, we have two detectors placed at a certain distance from each other along the flight path of the K<sup>+</sup>. The first detector is called the "Upstream detector" and measures the direction and the momentum of the Kaons in the beam: The second detector, called "Downstream detector", is composed of a tracking detector to detect the π<sup>+</sup> and a calorimeter to measure the π<sup>0</sup>.

However, the pions also have an average lifetime which is given by τ<sub>π</sub> = 2.6 * 10<sup>-8</sup>s. This gives us an average travel distance of Beta\*Gamma\*c\*τ<sub>π</sub> = 4.188km before a pion decays on average.

For the purpose of the project we were told to assume that the π<sup>0</sup> is stable and does not decay during the process.
We now have to create a simulation to simulate this decay, and find the optimum distance between the two detectors to capture as many events as possible, where a succesful event is considered the detection of both pions corresponding to the same kaon decay. In the first part of the simulation, we will assume a stream that has no spread and moves exactly on the z-axis. For the second part we will add a spread to the stream coming from the collimator following a Gaussian distribution.

## Results

### Determining τ<sub>K<sup>+</sup></sub> from the given data
The first part of the assignement is to determine the average lifetime of a K<sup>+</sup> particle. For this purpose, we are given the data from the earlier LifeK experiment listed in the file dec_lengths.dat that containing 100'000 measurements of K<sup>+</sup> particles and their decay length. Each row in that documents represents one measurements with the unit [meters].
To create a best fit of the data, we created tauestimator.py which plots different functions against the data to find the best fit.
![Fitting the measurement from the experiment KLife](https://raw.githubusercontent.com/abold/DAProject/master/K2Pi/plot2.jpg)
As we can see here, the function called "f1" is the best fit for the dataset. We can then use this function to determine the decay time τ<sub>K<sup>+</sup></sub> = 560m.

### Simulating the K<sup>+</sup> decay
In order to simulate the decay, we created a Monte Carlo simulation in python. Basically, we simulate the collimator shooting one K<sup>+</sup> particle at a time and with τ<sub>K<sup>+</sup></sub> that we got from the first part of the assignement we can calculate when this particle should decay. When the kaon decays we generate the π<sup>0</sup> and the π<sup>+</sup> particle with random angles θ and φ in the K<sup>+</sup> restframe and calculate the repective angles in the lab frame (using a Lorentz boost, matrix multiplication). Since the detectors are given in the instructions to be circular with 4m in diameter, we can now easily calculate where the particles will be at z=[location of sensor]. If this location is within the 4m of the sensors diameter we have a positive event registration.

This is how the decay looks like from the two different lab frames:
![Decay in K<sup>+</sup> and lab frame](https://raw.githubusercontent.com/RononDex/DAProject/master/K2Pi/Frames.png)
Where the angles θ and φ are uniformly spread between 0 to π for θ, and 0 to 2π for φ for all the decays.
The time that the kaon needs to decay is calculated by creating a contineous random number which is distributed exponentially using the following code: 
```python
vlen= np.array(stats.expon.rvs(loc=0, scale=tau, size=1))
```
Where tau is τ<sub>K<sup>+</sup></sub> and stats is the library scipystats.
After the decay happend, we just let the π<sup>+</sup> and π<sup>0</sup> flying in their respective direction as seen from the lab frame and check if they hit the sensors.

The simulation script "K2Pi.py" produces a plot where on the z axis those particles can be detected by the sensor and plots the amount of positive events for a given distance between the sensors, given that the first one is right in front of the collimator to pick up all the K<sup>+</sup>. This is the simulation with no spread applied at the collimator:
![Simulation with 300 K<sup>+</sup>, no spread](https://raw.githubusercontent.com/RononDex/DAProject/master/K2Pi/Simulation50NoSpread.png)

The point where the most events get registered is marked in the plot, which in this case is at **250m** distance with a maximum efficency of **~30%**.
Note that every run of the simulation will produce a different result, since this process has (as explained above) an element of randomness in it, however the bigger the number of simulated K<sup>+</sup> the less the spread of the found optimal distances between the simulations get.

Simulation with a spreading beam of K<sup>+</sup> at the collimator of 0.001:
![Simulation with 300 K<sup>+</sup>, with spread](https://raw.githubusercontent.com/RononDex/DAProject/master/K2Pi/Simulation50.png)

As we can see, the results are a bit different. The optimal distance is determined to be **364m** with an effiency of **~35%**

For the simulation we used the following parameters σ<sub>x</sub>=σ<sub>y</sub> = 1mrad:
 - E<sub>π<sup>0</sup></sub> = 245.563588 MeV     Energy of π<sup>0</sup> in K<sup>+</sup> rest frame
 - E<sub>π<sup>+</sup></sub> = 248.118174 MeV     Energy of π<sup>+</sup> in K<sup>+</sup> rest frame
 - p = 205.14091 Mev / c                          Impuls of the pions (same for both), in  K<sup>+</sup> rest frame
 - β = 0.99997833784995                           Beta factor of a K<sup>+</sup> rest frame
 - γ = 151.92756392754                            Corresponding gamma factor to the beta factor
 - τ = 560 m                                      τ as calculated in the first part of the experiment
