# file: sum.W31-0.2.pH.th.v6.merged_var_mass_cb2.log_th_vol.psuade.th00675.tab
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
# MSE: 0.3789, GCV: 0.3797, RSQ: 0.5984, GRSQ: 0.5975

def model(example_iterator):
    accessors = [lambda x: 10.1658576990974,
		lambda x: 2.219678468003535 * max(0, x[4] - 6.79785387005),
		lambda x: -0.6589910341249411 * max(0, 6.79785387005 - x[4]),
		lambda x: 538.452037110101 * x[6],
		lambda x: -0.020936024047622057 * x[7] * x[6],
		lambda x: 307237.72974987247 * max(0, x[3] - 0.0017521741) * x[6],
		lambda x: -307253.1457442677 * max(0, 0.0017521741 - x[3]) * x[6],
		lambda x: -1901.564124722666 * x[12] * x[6],
		lambda x: -0.00024488011011647437 * x[13] * max(0, x[3] - 0.0017521741) * x[6],
		lambda x: -305.65101748545567 * x[5] * x[12] * x[6],
		lambda x: 0.1035993630066514 * max(0, x[4] - 7.8707083207) * x[6],
		lambda x: -0.12780277954880148 * max(0, 7.8707083207 - x[4]) * x[6],
		lambda x: -12.488332869092119 * x[5] * x[5] * x[12] * x[6],
		lambda x: -0.8908438824801124 * max(0, x[4] - 4.91260925765) * max(0, 6.79785387005 - x[4]),
		lambda x: -0.1425999758648686 * max(0, 4.91260925765 - x[4]) * max(0, 6.79785387005 - x[4]),
		lambda x: -307237.35733849934 * max(0, x[3] - 4.00110558e-07) * x[6]]
    for x in example_iterator:
        yield sum(accessor(x) for accessor in accessors)
    
