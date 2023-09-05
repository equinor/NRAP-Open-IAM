# file: sum.W31-0.2.TDS.th.v6.merged_var_mass_cb2.log_dy.psuade.th00100.tab
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
# MSE: 0.0331, GCV: 0.0331, RSQ: 0.7626, GRSQ: 0.7621

def model(example_iterator):
    accessors = [lambda x: 2.150143637852732,
		lambda x: 0.7389040007103429 * max(0, x[4] - 6.48181194674),
		lambda x: -0.20603514788653995 * max(0, 6.48181194674 - x[4]),
		lambda x: 12.287155263146344 * max(0, x[1] - 8.22559475e-05),
		lambda x: -1778.3830941682072 * max(0, 8.22559475e-05 - x[1]),
		lambda x: -0.00037325398639580926 * x[13],
		lambda x: 0.8160584937133333 * max(0, x[3] - 0.000145693584),
		lambda x: 8097.92617295 * max(0, 0.000145693584 - x[3]),
		lambda x: 806.1367744340608 * max(0, x[4] - 6.6092090464) * max(0, 8.22559475e-05 - x[1]),
		lambda x: -3900.7385633353415 * max(0, 6.6092090464 - x[4]) * max(0, 8.22559475e-05 - x[1]),
		lambda x: 1.233497641806025e-06 * x[13] * x[13],
		lambda x: 7.35240064386744e-08 * x[7] * x[13] * x[13],
		lambda x: -2.7773188779519273 * max(0, x[3] - 0.00047641042) * max(0, x[3] - 0.000145693584),
		lambda x: 0.026589640101065015 * x[6] * max(0, x[4] - 6.48181194674),
		lambda x: 1.2107330803701188e-05 * x[13] * x[6] * max(0, x[4] - 6.48181194674),
		lambda x: -537.2237478865609 * max(0, x[1] - 7.97520842e-07) * max(0, x[3] - 0.00047641042) * max(0, x[3] - 0.000145693584)]
    for x in example_iterator:
        yield sum(accessor(x) for accessor in accessors)
    
