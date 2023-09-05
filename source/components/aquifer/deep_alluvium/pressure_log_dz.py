# file: sum.W31-0.2.P.deltabg.th.v6.merged_var_mass_cb2.log_dz.psuade.th00500.tab
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
# MSE: 0.0428, GCV: 0.0429, RSQ: 0.6329, GRSQ: 0.6316

def model(example_iterator):
    accessors = [lambda x: -0.9618991576879029,
		lambda x: 0.15809779083300654 * max(0, x[4] - 5.17513821994),
		lambda x: 0.06742255864978584 * max(0, 5.17513821994 - x[4]),
		lambda x: 0.0025453471294845115 * x[13],
		lambda x: -0.0021996756209102294 * max(0, x[3] - 0.015892439) * x[13],
		lambda x: -0.11256619897150606 * x[6],
		lambda x: -0.002944116174703887 * max(0, x[4] - 6.37591487044) * max(0, 0.015892439 - x[3]) * x[13],
		lambda x: 0.029684709092167263 * max(0, x[2] - 4.77308819634) * max(0, x[4] - 6.37591487044) * max(0, 0.015892439 - x[3]) * x[13],
		lambda x: 0.012533235779843721 * max(0, 4.77308819634 - x[2]) * max(0, x[4] - 6.37591487044) * max(0, 0.015892439 - x[3]) * x[13],
		lambda x: 7.67674738714249e-05 * x[9] * max(0, 5.17513821994 - x[4]),
		lambda x: 0.07159807738829141 * max(0, x[4] - 7.88005928038) * x[6],
		lambda x: 0.02098284454534656 * max(0, 7.88005928038 - x[4]) * x[6],
		lambda x: 0.0015873677100408448 * max(0, x[4] - 7.84773482564) * x[13],
		lambda x: 9.595629720024305e-05 * max(0, 7.84773482564 - x[4]) * x[13],
		lambda x: -7.019685525355701e-07 * x[13] * max(0, x[4] - 7.84773482564) * x[13],
		lambda x: -0.5598243834681581 * max(0, 0.0758275379 - x[3]) * max(0, x[4] - 5.17513821994),
		lambda x: -2.202828021038123 * max(0, x[2] - 5.90220405517) * max(0, 0.0758275379 - x[3]) * max(0, x[4] - 5.17513821994),
		lambda x: -0.16920148335941068 * max(0, 5.90220405517 - x[2]) * max(0, 0.0758275379 - x[3]) * max(0, x[4] - 5.17513821994),
		lambda x: 0.00024332108739748282 * x[5] * x[13],
		lambda x: -0.1901291991266642 * x[5],
		lambda x: 11.392341206463431 * max(0, x[3] - 0.00804998469) * max(0, 5.90220405517 - x[2]) * max(0, 0.0758275379 - x[3]) * max(0, x[4] - 5.17513821994),
		lambda x: -75.28193442078607 * max(0, 0.00804998469 - x[3]) * max(0, 5.90220405517 - x[2]) * max(0, 0.0758275379 - x[3]) * max(0, x[4] - 5.17513821994),
		lambda x: -0.00021073497304101296 * max(0, x[3] - 1.3816757e-07) * x[5] * x[13],
		lambda x: -109.4698204820532 * max(0, 1.3816757e-07 - x[3]) * x[5] * x[13]]
    for x in example_iterator:
        yield sum(accessor(x) for accessor in accessors)
    
