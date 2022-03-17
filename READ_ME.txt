The source code is broken up into different codes that can handle a different amount of data production and filtering.

Star_main.c has the ability to return one document called STAR_RESULT.csv and STAR_RESULT_2.csv. It will not create a large amount of result files and will return any results as it will not filter for certain types of graphs.

Star_main_data.c has the abilities that the previous contains but also can filter out  non oscillating states. It also will produce files for every one of these states within the range N=1 to N=2.

Star_main_freq.c has the abilities that the previous contains but also can filter out  non oscillating states. It also will produce files for every one of these states within the range N=1 to N=2. Finally it can find the frequency and amplitude of each oscillating state in the Frequency_and_Amplitude_vs_n.csv file.

Star_main_data_opp.c prints out the stationary state of each stationary state system in the range of n values 1 to 3.5. These numbers get printed out in a file called convergent.csv