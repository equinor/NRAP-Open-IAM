# file: sum.W31-0.2.TDS.th.v6.merged_var_mass_cb2.log_dz.psuade.th00100.tab
# dir: /erdfilespace/userspace/mansoor1/work/co2/kimberlina/nuftresults/prod07/thresholds/pyearth
# variables: 
# x[0]   time
# x[1]   brnFlux
# x[2]   log_brnMass
# x[3]   co2Flux
# x[4]   log_co2Mass
# x[5]   logK_sand1
# x[6]   logK_sand2
# x[7]   logK_sand3
# x[8]   logK_caprock
# x[9]   corrlationLengthX
# x[10]  corrlationLengthZ
# x[11]  sandFraction
# x[12]  groundwater_gradient
# x[13]  leak_depth
# output: log_dz
# MSE: 0.0218, GCV: 0.0219, RSQ: 0.6708, GRSQ: 0.6700

def model(example_iterator):
    accessors = [lambda x: 2.1451701745728085,
		lambda x: 0.3008349050534548 * max(0, x[4] - 7.05083367049),
		lambda x: -0.24904001543763954 * max(0, 7.05083367049 - x[4]),
		lambda x: 0.0002948665559537632 * x[13],
		lambda x: 0.30807566263046143 * max(0, x[3] - 0.0180977812),
		lambda x: -70.68956724396432 * max(0, 0.0180977812 - x[3]),
		lambda x: 0.061105260533145855 * max(0, x[2] - 3.87795780554) * max(0, 7.05083367049 - x[4]),
		lambda x: -0.07805467157571436 * max(0, 3.87795780554 - x[2]) * max(0, 7.05083367049 - x[4]),
		lambda x: 8.823896411191967e-05 * max(0, 7.96221573551 - x[4]) * x[13],
		lambda x: -2.551539175499329e-07 * x[13] * x[13],
		lambda x: -4.829913310568507 * x[6] * max(0, 0.0180977812 - x[3]),
		lambda x: 10.048455350271944 * max(0, x[2] - 4.64560328174) * max(0, 0.0180977812 - x[3]),
		lambda x: 5.284024112202553 * max(0, 4.64560328174 - x[2]) * max(0, 0.0180977812 - x[3]),
		lambda x: 661.1085455784487 * max(0, x[3] - 2.45883139e-05) * max(0, 0.0180977812 - x[3]),
		lambda x: 971387.5896588733 * max(0, 2.45883139e-05 - x[3]) * max(0, 0.0180977812 - x[3]),
		lambda x: -2.4103754276083578e-05 * x[9] * max(0, x[4] - 7.05083367049),
		lambda x: 0.1962228378911135 * x[11]]
    for x in example_iterator:
        yield sum(accessor(x) for accessor in accessors)
    
