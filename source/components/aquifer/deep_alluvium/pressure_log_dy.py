# file: sum.W31-0.2.P.deltabg.th.v6.merged_var_mass_cb2.log_dy.psuade.th00500.tab
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
# MSE: 0.0910, GCV: 0.0914, RSQ: 0.6219, GRSQ: 0.6205

def model(example_iterator):
    accessors = [lambda x: 2.8323242414056833,
		lambda x: -1.7838820731304588 * max(0, x[4] - 6.66270992527),
		lambda x: -0.04914003109069534 * max(0, 6.66270992527 - x[4]),
		lambda x: 0.0020050879588844016 * x[13] * max(0, x[4] - 6.66270992527),
		lambda x: -5.335240621244437 * max(0, x[3] - 0.0140519626),
		lambda x: -24.372130646994155 * max(0, 0.0140519626 - x[3]),
		lambda x: -756.582670386443 * max(0, x[6] - -10.722) * max(0, x[4] - 6.66270992527),
		lambda x: 0.13584937320717796 * max(0, -10.722 - x[6]) * max(0, x[4] - 6.66270992527),
		lambda x: -0.8008797192494356 * x[11],
		lambda x: 3.9773438299098958 * max(0, x[3] - 1.63616083e-07) * x[11],
		lambda x: 0.0024884261984610134 * x[13] * x[11],
		lambda x: -0.0002272579855002732 * max(0, x[4] - 7.03638221095) * x[13] * x[11],
		lambda x: 0.00015679506448584846 * max(0, 7.03638221095 - x[4]) * x[13] * x[11],
		lambda x: -2.8919591841258807e-06 * x[13] * x[13] * x[11],
		lambda x: 6.084101755732263e-08 * x[5] * x[13] * x[13] * x[11],
		lambda x: -0.1372334317569015 * x[8] * max(0, x[4] - 6.66270992527),
		lambda x: 0.00011763645250084664 * x[8] * x[13] * max(0, x[4] - 6.66270992527),
		lambda x: -5.211177509067966e-06 * max(0, x[3] - 0.0147711665) * x[13] * x[13] * x[11],
		lambda x: 2.2703664640744137e-05 * max(0, 0.0147711665 - x[3]) * x[13] * x[13] * x[11],
		lambda x: 1.6003725011870529e-09 * x[13] * x[13] * x[13] * x[11],
		lambda x: -0.06117600320708419 * max(0, x[4] - 7.15839256957) * max(0, -10.722 - x[6]) * max(0, x[4] - 6.66270992527),
		lambda x: 0.501745440364882 * max(0, 7.15839256957 - x[4]) * max(0, -10.722 - x[6]) * max(0, x[4] - 6.66270992527),
		lambda x: 0.00914398996075777 * max(0, x[3] - 0.0845492753) * x[13] * x[11],
		lambda x: -0.010432719164409576 * max(0, 0.0845492753 - x[3]) * x[13] * x[11],
		lambda x: -0.00023401765261610308 * x[9],
		lambda x: 1.2417553563182082e-07 * x[9] * x[9]]
    for x in example_iterator:
        yield sum(accessor(x) for accessor in accessors)
    
