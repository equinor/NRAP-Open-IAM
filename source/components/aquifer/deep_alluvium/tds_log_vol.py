# file: sum.W31-0.2.TDS.th.v6.merged_var_mass_cb2.log_th_vol.psuade.th00100.tab
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
# MSE: 0.0719, GCV: 0.0722, RSQ: 0.8713, GRSQ: 0.8707

def model(example_iterator):
    accessors = [lambda x: 31.55352744989693,
		lambda x: -4.616185250343749 * max(0, x[4] - 7.12034869592),
		lambda x: -0.5702666597242113 * max(0, 7.12034869592 - x[4]),
		lambda x: 4.179339340712326 * x[6],
		lambda x: -0.013064514307434778 * x[13] * max(0, x[4] - 7.12034869592),
		lambda x: 19.344174760484762 * x[11],
		lambda x: 0.02198703904091611 * max(0, x[4] - 5.92828372946) * x[6],
		lambda x: -0.04852058002655423 * max(0, 5.92828372946 - x[4]) * x[6],
		lambda x: -1.7485025414594144 * max(0, x[3] - 0.0532711818) * max(0, x[4] - 7.12034869592),
		lambda x: 18.59246904941476 * max(0, 0.0532711818 - x[3]) * max(0, x[4] - 7.12034869592),
		lambda x: -0.0068004346073278765 * max(0, 5.1907280921 - x[2]) * x[6],
		lambda x: 1.9289513765951034 * x[7] * max(0, 0.0532711818 - x[3]) * max(0, x[4] - 7.12034869592),
		lambda x: 2.9808632546197678 * x[5] * x[11],
		lambda x: 0.17129145323426886 * x[6] * x[6],
		lambda x: -0.0016831799621908061 * x[7] * x[13] * max(0, x[4] - 7.12034869592),
		lambda x: -0.00019068258412513472 * x[6] * x[7] * x[13] * max(0, x[4] - 7.12034869592),
		lambda x: -0.001713956825598828 * x[6] * x[13] * max(0, x[4] - 7.12034869592),
		lambda x: -0.35552551210589983 * x[7] * max(0, x[4] - 7.12034869592),
		lambda x: -8.142143857653135 * x[12] * max(0, x[4] - 5.92828372946) * x[6],
		lambda x: -0.0066204721969937985 * max(0, x[3] - 0.00467543483) * x[13] * max(0, x[4] - 7.12034869592),
		lambda x: 0.8548090171843457 * max(0, 0.00467543483 - x[3]) * x[13] * max(0, x[4] - 7.12034869592),
		lambda x: -0.019524647379064716 * max(0, x[3] - 0.0857157005) * max(0, x[2] - 5.1907280921) * x[6],
		lambda x: -3.1046866555851778 * max(0, 0.0857157005 - x[3]) * max(0, x[2] - 5.1907280921) * x[6],
		lambda x: 0.010648050334052028 * max(0, x[4] - 7.19070084294) * x[6] * x[6],
		lambda x: -1144.9748872008258 * max(0, x[1] - 6.21416514e-05) * x[7] * max(0, 0.0532711818 - x[3]) * max(0, x[4] - 7.12034869592),
		lambda x: -9524.946620049328 * max(0, 6.21416514e-05 - x[1]) * x[7] * max(0, 0.0532711818 - x[3]) * max(0, x[4] - 7.12034869592),
		lambda x: 0.11896278669769589 * x[5] * x[5] * x[11],
		lambda x: -0.0008785166554119428 * max(0, x[3] - 0.0581129332) * x[6] * x[13] * max(0, x[4] - 7.12034869592),
		lambda x: 0.0003741723071666314 * max(0, 0.0581129332 - x[3]) * x[6] * x[13] * max(0, x[4] - 7.12034869592),
		lambda x: 0.0002989497507916994 * x[9] * x[7] * max(0, 0.0532711818 - x[3]) * max(0, x[4] - 7.12034869592),
		lambda x: -0.24365408725713658 * x[7] * max(0, 0.0857157005 - x[3]) * max(0, x[2] - 5.1907280921) * x[6],
		lambda x: -2.680464873172639e-06 * x[13] * max(0, x[3] - 0.00467543483) * x[13] * max(0, x[4] - 7.12034869592)]
    for x in example_iterator:
        yield sum(accessor(x) for accessor in accessors)
    
