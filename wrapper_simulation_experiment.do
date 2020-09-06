cls
clear all
set seed 10101

glo numobs = 1000
glo numrep = 100

code_simulation_experiment
return list

simulate ATEr = r(ATE) ATTr = r(ATT) ATUTr = r(ATUT) LATEr = r(LATE) naiver = r(naive) olsxr = r(olsx) heckr = r(heck) sharprdd1r = r(sharprdd1) sharprdd2r = r(sharprdd2) ivr = r(iv) fuzzyrdr = r(fuzzyrd) psattr=r(psatt) psater=r(psate), ///
reps($numrep) noleg nodots saving(selfselect, replace): selfselect
mean ATEr ATTr ATUTr LATEr naiver olsxr heckr sharprdd1r sharprdd2r ivr fuzzyrdr psattr psater

