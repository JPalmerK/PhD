Work flow
1) PropModel.m Calls ESCBellhop_transmission.m function to create the transmission loss grids
for a given frequency

2) MakeMeanTL grids combines TLgrids from multiple frequencies to make a mean transmission loss
grid which (generally) corrilates with the average transmission loss within an octave or 1/3octave
band

3a) Area vs 40khz NL loads the mean TL grid for the 40.1kHz octave band and calculates
the total area monitored as a function of ambient noise.

	This calls ENRMonitored which calculates the maximum range
	and the total area monitored greater than the RL threshold
	
	ENRMonitored now considers the integration time tao over which the detector is looking
	
	To plot the relationship between area monitored, tao and SNR threshold use 
	ENRMonitored_demoMode.m

3b) Pdet as a funciton of NL table (Pdet 41 kHz oct.txt) was created using the 
UnconditionalPdet as fx of noise.m which is housed under the Code file (as opposed to prop model)

I think (3a) is better way to proceed is to use 3a and just divide by the pi*r^2. 
Yes, this was done and 3 text files were produced:

-AreaVsNL_41KhzBand             Total area monitored
-MaximumRange41kHzBand	        Maximum range detected
-Pdet_41kHzBand_5km             Probability of detecting a click within 5km




4) Dist of Pdet Each SM2M.m then calculates the detection probability for
the times and locations for which there is SM coverage. Reads in AreaVsNL_41KhzBand from previous
and applies noise levels. Also includes the matlab date, noise level and area monitored
 
Outputs _NLandPdetVals.csv files which are output to 
W:\KJP PHD\4-Bayesian Habitat Use\Pdet at Time\2013




_NLandPdetVals used by R code to incorporate into Occupancy data

5) OcctableRWrite.r (in GitHub) creates the hourly occupancy table

6) Occ with and WO TL and NL.R  reads in the occupancy table (whch has the area monitored) and creates a Bayesian model to
estimate occupancy.

7) Plotting JAGS Model results.R plots the results from step 6
