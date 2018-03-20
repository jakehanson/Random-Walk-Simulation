# QUORUM LOADING (ENC RATE)

This is really the heart and soul of the quorum loading codes. I have got rid of brood for now and am just looking at probability to transport, assuming quorum is sensed relative to a threshold encounter rate.

I vary different parameters to try and get a feel for how they affect the nest site selection. Ideally, I will want to see how transport probability affects speed and accuracy, but for now, I am simply studing how the random walk parameters affect transport probability for a fixed quorum threshold.

I have some computationally costly (relatively speaking) dataframes created for the following parameters:

filename / a / inverse discovery rate / number of recruiters / velocity / r_enc 

dataframe1 / 0.5 / 300 / 50 / 0.1 / 0.1
dataframe2 / 0.5 / 300 / 100 / 0.1 / 0.1
dataframe3 / 0.5 / 150 / 50 / 0.1 / 0.1
dataframe4 / 0.5 / 300 / 50 / 0.05 / 0.1
dataframe5 / 1.0 / 300 / 50 / 0.05 / 0.1
dataframe6 / 2.0 / 300 / 50 / 0.05 / 0.1
dataframe7 / 1.0 / 300 / 50 / 0.1 / 0.1
dataframe8 / 2.0 / 300 / 50 / 0.1 / 0.1

Note, within the file, most info is there (in addition to collision and time inside info) but it is missing velocity and encounter radius
