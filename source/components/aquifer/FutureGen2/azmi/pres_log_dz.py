'''
output: pres_log_dz
prediction score: 0.8450292115285494
Forward Pass
---------------------------------------------------------------
iter  parent  var  knot  mse       terms  gcv    rsq    grsq
---------------------------------------------------------------
0     -       -    -     3.709868  1      3.711  0.000  0.000
1     0       10   1596  2.217294  3      2.222  0.402  0.401
2     0       9    3216  1.359155  5      1.365  0.634  0.632
3     3       10   4374  1.064131  7      1.070  0.713  0.712
4     0       3    -1    0.902533  8      0.909  0.757  0.755
5     7       10   145   0.829518  10     0.837  0.776  0.775
6     7       8    -1    0.780063  11     0.788  0.790  0.788
7     7       9    1950  0.741474  13     0.750  0.800  0.798
8     0       9    1950  0.699875  15     0.709  0.811  0.809
9     7       9    2100  0.666295  17     0.677  0.820  0.818
10    0       10   145   0.639918  19     0.651  0.828  0.825
11    1       0    -1    0.625352  20     0.637  0.831  0.828
12    14      8    -1    0.616864  21     0.629  0.834  0.831
13    14      10   4112  0.608916  23     0.622  0.836  0.832
14    18      10   3801  0.601435  25     0.615  0.838  0.834
15    2       10   1716  0.596078  27     0.611  0.839  0.835
16    13      8    -1    0.591245  28     0.606  0.841  0.837
17    13      9    440   0.585822  30     0.602  0.842  0.838
18    4       8    -1    0.580958  31     0.598  0.843  0.839
19    18      9    2840  0.576982  33     0.595  0.844  0.840
20    2       9    4672  0.570849  35     0.589  0.846  0.841
21    7       5    -1    0.566787  36     0.586  0.847  0.842
22    17      8    -1    0.563507  37     0.583  0.848  0.843
---------------------------------------------------------------
Stopping Condition 2: Improvement below threshold

Pruning Pass
--------------------------------------------
iter  bf  terms  mse   gcv    rsq    grsq
--------------------------------------------
0     -   37     0.56  0.583  0.848  0.843
1     3   36     0.56  0.582  0.848  0.843
2     16  35     0.56  0.582  0.848  0.843
3     2   34     0.56  0.581  0.848  0.843
4     12  33     0.56  0.581  0.848  0.843
5     6   32     0.56  0.580  0.848  0.844
6     1   31     0.56  0.580  0.848  0.844
7     34  30     0.57  0.582  0.847  0.843
8     32  29     0.57  0.582  0.847  0.843
9     29  28     0.57  0.584  0.846  0.843
10    36  27     0.57  0.586  0.846  0.842
11    21  26     0.57  0.589  0.845  0.841
12    28  25     0.58  0.591  0.844  0.841
13    22  24     0.58  0.594  0.843  0.840
14    35  23     0.58  0.597  0.842  0.839
15    23  22     0.59  0.601  0.841  0.838
16    30  21     0.59  0.605  0.840  0.837
17    27  20     0.60  0.610  0.839  0.836
18    25  19     0.61  0.618  0.836  0.833
19    18  18     0.61  0.623  0.835  0.832
20    33  17     0.62  0.628  0.833  0.831
21    31  16     0.62  0.630  0.833  0.830
22    20  15     0.63  0.638  0.830  0.828
23    19  14     0.66  0.666  0.823  0.821
24    10  13     0.72  0.728  0.806  0.804
25    4   12     0.78  0.793  0.789  0.786
26    14  11     0.79  0.793  0.788  0.786
27    11  10     0.79  0.799  0.786  0.785
28    13  9      0.81  0.812  0.783  0.781
29    26  8      0.88  0.884  0.763  0.762
30    8   7      0.94  0.941  0.748  0.746
31    17  6      1.10  1.107  0.703  0.702
32    5   5      1.30  1.309  0.649  0.647
33    7   4      1.62  1.625  0.563  0.562
34    24  3      1.82  1.821  0.510  0.509
35    9   2      2.51  2.508  0.325  0.324
36    15  1      3.71  3.711  0.000  0.000
--------------------------------------------
Selected iteration: 6

Earth Model
--------------------------------------------------------------------------
Basis Function                                       Pruned  Coefficient
--------------------------------------------------------------------------
(Intercept)                                          No      -17.9631
h(log_brine_mass-14.2749)                            Yes     None
h(14.2749-log_brine_mass)                            Yes     None
h(log_co2_mass-10.2625)                              Yes     None
h(10.2625-log_co2_mass)                              No      -2.69095
h(log_brine_mass-11.582)*h(log_co2_mass-10.2625)     No      -0.016525
h(11.582-log_brine_mass)*h(log_co2_mass-10.2625)     Yes     None
log_permh                                            No      -0.645299
h(log_brine_mass-18.7198)*log_permh                  No      0.189814
h(18.7198-log_brine_mass)*log_permh                  No      0.0311623
time*log_permh                                       No      0.0021865
h(log_co2_mass-14.9869)*log_permh                    No      0.396221
h(14.9869-log_co2_mass)*log_permh                    Yes     None
h(log_co2_mass-14.9869)                              No      2.10704
h(14.9869-log_co2_mass)                              No      2.71575
h(log_co2_mass-10.2979)*log_permh                    No      -0.240337
h(10.2979-log_co2_mass)*log_permh                    Yes     None
h(log_brine_mass-18.7198)                            No      2.40441
h(18.7198-log_brine_mass)                            No      -0.328042
thick*h(log_brine_mass-14.2749)                      No      0.00199012
time*h(14.9869-log_co2_mass)                         No      0.0040078
h(log_brine_mass-20.7015)*h(14.9869-log_co2_mass)    No      -0.0195462
h(20.7015-log_brine_mass)*h(14.9869-log_co2_mass)    No      -0.0168903
h(log_brine_mass-17.4138)*h(18.7198-log_brine_mass)  No      1.23066
h(17.4138-log_brine_mass)*h(18.7198-log_brine_mass)  No      0.0498319
h(log_brine_mass-6.93715)*h(14.2749-log_brine_mass)  No      0.0345807
h(6.93715-log_brine_mass)*h(14.2749-log_brine_mass)  No      -0.0486245
time*h(log_co2_mass-14.9869)                         No      0.00300549
h(log_co2_mass-21.0521)*h(log_co2_mass-14.9869)      No      -0.0301984
h(21.0521-log_co2_mass)*h(log_co2_mass-14.9869)      No      -0.0316196
time*h(10.2625-log_co2_mass)                         No      -0.00395007
h(log_co2_mass-13.3837)*h(18.7198-log_brine_mass)    No      0.051253
h(13.3837-log_co2_mass)*h(18.7198-log_brine_mass)    No      0.0251944
h(log_co2_mass-13.008)*h(14.2749-log_brine_mass)     No      -0.0546435
h(13.008-log_co2_mass)*h(14.2749-log_brine_mass)     No      -0.0127379
rel_vol_frac_calcite*log_permh                       No      -0.0177712
time*h(log_brine_mass-18.7198)                       No      0.00233833
--------------------------------------------------------------------------
MSE: 0.5640, GCV: 0.5801, RSQ: 0.8480, GRSQ: 0.8437
parameters:
X[ 0 ] =  thick
X[ 1 ] =  depth
X[ 2 ] =  por
X[ 3 ] =  log_permh
X[ 4 ] =  log_aniso
X[ 5 ] =  rel_vol_frac_calcite
X[ 6 ] =  log_co2_rate
X[ 7 ] =  log_brine_rate
X[ 8 ] =  time
X[ 9 ] =  log_co2_mass
X[ 10 ] =  log_brine_mass
'''
def model(example_iterator):
    accessors = [
        lambda x: -17.963051223317887,
        lambda x: -2.690954842503831 * max(0, 10.262462567761771 - x[9]),
        lambda x: -0.016524997455952622 * max(0, x[10] - 11.582000641026282) *\
            max(0, x[9] - 10.262462567761771),
        lambda x: -0.6452990536390911 * x[3],
        lambda x: 0.18981411530323722 * max(0, x[10] - 18.719795392664363) * x[3],
        lambda x: 0.031162319038241182 * max(0, 18.719795392664363 - x[10]) * x[3],
        lambda x: 0.002186504347829299 * x[8] * x[3],
        lambda x: 0.39622146637743627 * max(0, x[9] - 14.986948904468676) * x[3],
        lambda x: 2.1070437459016067 * max(0, x[9] - 14.986948904468676),
        lambda x: 2.7157528219800056 * max(0, 14.986948904468676 - x[9]),
        lambda x: -0.24033727419701936 * max(0, x[9] - 10.297915923330548) * x[3],
        lambda x: 2.404407355986845 * max(0, x[10] - 18.719795392664363),
        lambda x: -0.3280424150343512 * max(0, 18.719795392664363 - x[10]),
        lambda x: 0.001990115232980094 * x[0] * max(0, x[10] - 14.274869313706027),
        lambda x: 0.0040077954621718925 * x[8] * max(0, 14.986948904468676 - x[9]),
        lambda x: -0.019546225554649254 * max(0, x[10] - 20.701468467610617) *\
            max(0, 14.986948904468676 - x[9]),
        lambda x: -0.01689029884767701 * max(0, 20.701468467610617 - x[10]) *\
            max(0, 14.986948904468676 - x[9]),
        lambda x: 1.2306553028277991 * max(0, x[10] - 17.413838661845336) *\
            max(0, 18.719795392664363 - x[10]),
        lambda x: 0.04983192783293622 * max(0, 17.413838661845336 - x[10]) *\
            max(0, 18.719795392664363 - x[10]),
        lambda x: 0.0345806713665546 * max(0, x[10] - 6.9371532496084) *\
            max(0, 14.274869313706027 - x[10]),
        lambda x: -0.04862450332185715 * max(0, 6.9371532496084 - x[10]) *\
            max(0, 14.274869313706027 - x[10]),
        lambda x: 0.0030054859473457873 * x[8] * max(0, x[9] - 14.986948904468676),
        lambda x: -0.030198391224332814 * max(0, x[9] - 21.05214080475075) *\
            max(0, x[9] - 14.986948904468676),
        lambda x: -0.03161963455128847 * max(0, 21.05214080475075 - x[9]) *\
            max(0, x[9] - 14.986948904468676),
        lambda x: -0.003950065814491255 * x[8] * max(0, 10.262462567761771 - x[9]),
        lambda x: 0.05125295101998732 * max(0, x[9] - 13.383697091420652) *\
            max(0, 18.719795392664363 - x[10]),
        lambda x: 0.02519437505745703 * max(0, 13.383697091420652 - x[9]) *\
            max(0, 18.719795392664363 - x[10]),
        lambda x: -0.05464351911722548 * max(0, x[9] - 13.007992634603008) *\
            max(0, 14.274869313706027 - x[10]),
        lambda x: -0.012737939138882104 * max(0, 13.007992634603008 - x[9]) *\
            max(0, 14.274869313706027 - x[10]),
        lambda x: -0.017771196828043796 * x[5] * x[3],
        lambda x: 0.002338330704485886 * x[8] * max(0, x[10] - 18.719795392664363)]
    for x in example_iterator:
        yield sum(accessor(x) for accessor in accessors)
