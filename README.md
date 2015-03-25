# PermVsBoot
## Exploring differences between resampling methods: Bootstrap vs Permutation

### The goal of this project is to exlpore differences in the distributions of p-values from two different resampling methods, bootstrap and permution under the null hypothesis for EEG data. 

  The statisic under exploration is peak to peak.  Peak to peak (p2p) is a method commonly used in signal processing and is calculated by finding the distance between a peak and a trough. Usually instead of identifying a single point to describe the peak/trough a mean amplitude window is used. For example one will look for the maximum 100 ms window within a range of the signal and then look for the minimum 100 ms window and calculate the difference between them. The p2p function can be found in the permVsBoot.R script.
```
peak2peak_stat<-function(signal,start,max_end,end,size){    
  ...
  ...
  return (max_avg_peak - min_avg_peak)
}
```
  The test creates noise EEG data sets using the createnoiseeeg() function, which generates noise at the spectrum of human EEG. The noiseR function is adopted from the: [Generation of simulated EEG data project](http://www.cs.bris.ac.uk/~rafal/phasereset/). 

  The **bootstrap** method is applied as following: At each resampling a new surrogate ERP is generated for each condition, by selecting with replacement from the available single trials (of each condition). Then the p2p value is calculated for each of the two surrogate ERPs and the difference of p2p values is then stored. The p-value is calculated as the proportion of the values (i.e. differences) that are greater than 0.
  
  For the **permutation**, at each resampling a new surrogate ERP is generated for each condition by randomly assigning single trials to conditions. Then the p2p value is calculated for each of the surrogate ERPs and the difference between them is stored. The p-value is calculated as the proportion of values (in this case p2p differences) that are greater than the true observerd value (value of the original set)

## Output
  The script saves the results (p-values) to an .Rda file, a .csv file and also creates a .pdf image that compares the two distributions.

## Parametrisation
  There are several variables that can be parametrised for the tests. These are:
* trials: number of trials per condition
* len: number of time points per trial. It is set to 2200 based on the time file (time.txt)
* ns: scaling factor of the noise (not important for this test)
* p2p_start: point where to start looking for the peak. Selecting realistic time periods to search for p2p based on the time file availabe here, in order to keep consistent for tests with signal. It is set to 801 which corresponds to 300ms
* p2p_max_end: point where to stop looking for the peak. It is set to 2001 which corresponds to 1000ms
* p2p_end: point where to stop looking for the miminum. It is set to the end of the signal 
* p2p_win: size of the p2p internal window. It is set to 200 timepoints which corresponds to a 100ms window as is commonly used
* numberoftests: how many tests to perform. This is the number of p-values produced by the test
* numberofreps: how many resamplings to do at each step. Usually set to a minimum of 1000 repetitions

## Data
  There are three data files available:
* meanpower.txt: used in the noiseR() function to generate noise at the spectrum of human EEG
* signal.txt: an ERP signal. It can be used to generate EEG data set with signal using the createnoiseeeg_withsignal() function. It is not currently used for this test.
* time.txt: a file with time recording. 

## Pseudocode
**START TEST**
* **FOR EACH TEST** (tests = 10000)
  * Generate 2 EEG sets (first condition vs second condition)
    * first condition i.e. (fc<-createnoiseeeg(trials,len,ns))
    * sc<-createnoiseeeg(trials,len,ns)
    * No difference between them
  * Generate 2 ERPs
    * First condition ERP 
    * Second Condition ERP
    * Calculate p2p measurement for each one of these
    * Save difference to be used as the true observed value for the permutation
  * **START RESAMPLING**
    * **FOR EACH RESAMPLING**
      * Generate bootstraped first condition by selecting with replacement from the set
		of single trials of the first condition.Generate boostraped second condition by selecting with replacement from the set 
		of single trials of the second condition. Generate two ERPs  from the bootstraped EEGs. Calculate p2p on each ERP. Save the difference.
      * Create two permuted conditions by randomly choosing from the set of all available trials. Generate two ERPs, calculate p2p and save difference.
      * **END RESAMPLING**
  * **END RESAMPLINGS/TEST**
  * Calculate bootstrap p-value = 1- (number of bootstrapped differences greater 
	  									  than 0/ number of resamplings)
  * Calculate permutation p-value = number of permuted differences greater than true
	  								    observed/ number of resampling
* **END TESTS**

Plot comparative histogram of bootstrap p-values and permutation p-values. Save .Rda, .csv with the results.

## Published work

The results of this test have been published in: [Zoumpoulaki, Alexia, Abdulmajeed Alsufyani, and Howard Bowman. "Resampling the peak, some dos and don'ts." *Psychophysiology* (2014)] (http://onlinelibrary.wiley.com/doi/10.1111/psyp.12363/abstract?systemMessage=Wiley+Online+Library+will+be+disrupted+on+7th+March+from+10%3A00-13%3A00+GMT+%2805%3A00-08%3A00+EST%29+for+essential+maintenance.++Apologies+for+the+inconvenience.&userIsAuthenticated=true&deniedAccessCustomisedMessage=). 

Comparison of p values obtained from 10,000 bootstrap tests vs. 10,000 permutation tests. p Values were obtained for p2p measurement on simulated noise EEG data. Using the chi-squared test to check for uniformity, we can reject the hypothesis that the bootstrap distribution is uniform (p value < 22 × 10−17 while for permutation p value = 0.3915).


![Comparision of distributions generated from the test. Results from 10000 tests.] [image]
[image]: https://github.com/CogSys-Kent/PermVsBoot/blob/master/figure2.jpg
