# file: sum.W31-0.2.P.deltabg.th.v6.merged_var_mass_cb2.log_th_vol.psuade.th00500.tab
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
# output: log_th_vol
# MSE: 0.1097, GCV: 0.1101, RSQ: 0.8401, GRSQ: 0.8395

def model(example_iterator):
    accessors = [lambda x: 6.508524484277225,
		lambda x: 0.19444823189623392 * max(0, x[4] - 6.92157813213),
		lambda x: -0.17420591042403347 * max(0, 6.92157813213 - x[4]),
		lambda x: -0.10172396187755231 * x[6] * max(0, x[4] - 6.92157813213),
		lambda x: -51.081340371660886 * max(0, x[3] - 7.63819188e-08),
		lambda x: 11341918.676514866 * max(0, 7.63819188e-08 - x[3]),
		lambda x: -0.0009185706055635373 * x[11] * x[13],
		lambda x: -0.0013551208748508875 * max(0, x[4] - 7.50653954789) * x[13],
		lambda x: 0.001845092101625142 * max(0, 7.50653954789 - x[4]) * x[13],
		lambda x: 5.462586705223657 * x[8] * max(0, x[3] - 7.63819188e-08),
		lambda x: -1.3870967095717788 * max(0, x[4] - 2.42078479894) * x[8] * max(0, x[3] - 7.63819188e-08),
		lambda x: -7569.94799003337 * max(0, 2.42078479894 - x[4]) * x[8] * max(0, x[3] - 7.63819188e-08),
		lambda x: 8.148327469825745e-05 * x[9],
		lambda x: -5.1860231906175613e-05 * max(0, x[2] - 3.6498546069) * x[13],
		lambda x: -0.00027544843032956123 * max(0, 3.6498546069 - x[2]) * x[13],
		lambda x: 0.01475366041995585 * max(0, x[3] - 0.0105054201) * max(0, x[4] - 7.50653954789) * x[13],
		lambda x: 0.026837811808945844 * max(0, x[3] - 0.0220667229) * x[13],
		lambda x: -0.027340178057784215 * max(0, 0.0220667229 - x[3]) * x[13],
		lambda x: 0.29002419945027214 * max(0, x[4] - 8.64355222755) * max(0, x[4] - 2.42078479894) * x[8] * max(0, x[3] - 7.63819188e-08),
		lambda x: -0.36131829337682575 * max(0, 8.64355222755 - x[4]) * max(0, x[4] - 2.42078479894) * x[8] * max(0, x[3] - 7.63819188e-08),
		lambda x: 0.00042243837378919125 * x[13] * max(0, x[4] - 2.42078479894) * x[8] * max(0, x[3] - 7.63819188e-08),
		lambda x: 0.0012739412486553192 * max(0, x[4] - 6.91047396298) * x[13],
		lambda x: -0.001739646540954709 * max(0, 6.91047396298 - x[4]) * x[13],
		lambda x: -28653078.419741187 * max(0, 4.72469861e-07 - x[1]) * max(0, x[4] - 2.42078479894) * x[8] * max(0, x[3] - 7.63819188e-08)]
    for x in example_iterator:
        yield sum(accessor(x) for accessor in accessors)
    
