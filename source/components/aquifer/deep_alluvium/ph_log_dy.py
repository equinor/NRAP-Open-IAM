# file: sum.W31-0.2.pH.th.v6.merged_var_mass_cb2.log_dy.psuade.th00675.tab
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
# output: log_dy
# MSE: 0.1143, GCV: 0.1145, RSQ: 0.6470, GRSQ: 0.6464

def model(example_iterator):
    accessors = [lambda x: 4.906431315925561,
		lambda x: -0.5458854995991794 * max(0, x[4] - 6.61070052675),
		lambda x: 0.03585625258438087 * max(0, x[2] - 4.0496099444),
		lambda x: -0.23163824136357627 * max(0, 4.0496099444 - x[2]),
		lambda x: 0.03196922118532711 * x[6],
		lambda x: -0.021744518640556626 * x[7] * x[6],
		lambda x: -1.1298222456126328 * max(0, x[4] - 5.08027277978) * max(0, 6.61070052675 - x[4]),
		lambda x: -0.11350520531221436 * max(0, 5.08027277978 - x[4]) * max(0, 6.61070052675 - x[4]),
		lambda x: 0.007422444632870437 * max(0, x[4] - 6.58863492394) * x[7] * x[6],
		lambda x: 0.0053923999802051314 * max(0, 6.58863492394 - x[4]) * x[7] * x[6],
		lambda x: -19.63318026670267 * x[12] * x[6],
		lambda x: 0.02492798342265612 * max(0, x[3] - 2.61517514e-07) * x[6],
		lambda x: -254053.5011882582 * max(0, 2.61517514e-07 - x[3]) * x[6]]
    for x in example_iterator:
        yield sum(accessor(x) for accessor in accessors)
    
