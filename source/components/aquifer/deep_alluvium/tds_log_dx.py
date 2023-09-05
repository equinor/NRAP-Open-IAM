# file: sum.W31-0.2.TDS.th.v6.merged_var_mass_cb2.log_dx.psuade.th00100.tab
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
# output: log_dx
# MSE: 0.0361, GCV: 0.0363, RSQ: 0.8300, GRSQ: 0.8293

def model(example_iterator):
    accessors = [lambda x: 25.715068096710024,
		lambda x: 0.3155440799548896 * max(0, x[4] - 6.20557856835),
		lambda x: -0.12581896287812747 * max(0, 6.20557856835 - x[4]),
		lambda x: 3.6707199503728773 * x[6],
		lambda x: -0.23362857767781955 * max(0, x[3] - 0.00969850686) * max(0, x[4] - 6.20557856835),
		lambda x: 7.090704362254155 * max(0, 0.00969850686 - x[3]) * max(0, x[4] - 6.20557856835),
		lambda x: -0.007242026622033131 * x[13],
		lambda x: -0.000407711463874999 * x[7] * x[13],
		lambda x: 0.1779780283441383 * x[6] * x[6],
		lambda x: 0.16956397949073182 * x[7],
		lambda x: 0.26881216070499353 * x[11] * max(0, x[4] - 6.20557856835),
		lambda x: -0.00034883636494909 * x[5] * x[6] * x[6],
		lambda x: -4.9623819591504414e-05 * x[5] * x[7] * x[13],
		lambda x: -0.03506394662043026 * x[12] * x[7] * x[13],
		lambda x: -0.28861020073279015 * max(0, x[3] - 9.19812492e-05),
		lambda x: 13945.122050364427 * max(0, 9.19812492e-05 - x[3]),
		lambda x: -0.01318138618129924 * x[6] * x[12] * x[7] * x[13],
		lambda x: 37.39850208934062 * max(0, x[2] - 4.39785883353) * max(0, 0.00969850686 - x[3]) * max(0, x[4] - 6.20557856835),
		lambda x: 25.40520116246589 * max(0, 4.39785883353 - x[2]) * max(0, 0.00969850686 - x[3]) * max(0, x[4] - 6.20557856835),
		lambda x: 74.0274567524958 * max(0, 0.000466385203 - x[3]) * x[12] * x[7] * x[13],
		lambda x: -50.407488035179476 * x[12] * x[12] * x[7] * x[13],
		lambda x: -6.042682884510597e-06 * x[7] * x[5] * x[7] * x[13],
		lambda x: -0.00019758739013298054 * x[7] * x[5] * x[6] * x[6],
		lambda x: 1.1969696241731498e-05 * x[13] * x[13],
		lambda x: 1.0107693434457943e-06 * x[13] * x[7] * x[13],
		lambda x: -2.394305467223745e-07 * max(0, x[3] - 0.139022421) * x[13] * x[13],
		lambda x: 1.24037964255308e-06 * max(0, 0.139022421 - x[3]) * x[13] * x[13],
		lambda x: 90.64482035606054 * max(0, x[1] - 4.51317652e-07) * max(0, x[3] - 9.19812492e-05)]
    for x in example_iterator:
        yield sum(accessor(x) for accessor in accessors)
    
