#include "common.h"

scalar *cyclotomic_factorisation_tree[LOG_R + 1];

scalar *bezout_coefficients_tree[LOG_R + 1];

#if (PARAM_Q == 1073741441)
scalar cyclotomic_factorisation_array[2*PARAM_R - 1] = {1, 67419063, 1006322378, 150088098, 923653343, 254756945, 818984496, 458419187, 615322254, 410770664, 662970777, 114739670, 959001771, 322382370, 751359071, 257349718, 816391723, 290773124, 782968317, 68809976, 1004931465, 83833776, 989907665, 304307791, 769433650, 157773404, 915968037, 238377628, 835363813, 292773971, 780967470, 38317637, 1035423804, 46040912, 1027700529, 54004485, 1019736956, 177815152, 895926289, 158070168, 915671273, 242891671, 830849770, 490914661, 582826780, 280198643, 793542798, 407600992, 666140449, 163749671, 909991770, 191657608, 882083833, 235535230, 838206211, 11928080, 1061813361, 251010649, 822730792, 71932388, 1001809053, 223226193, 850515248, 378408708, 695332733, 251403611, 822337830, 261404352, 812337089, 377362640, 696378801, 187213060, 886528381, 453044675, 620696766, 201106414, 872635027, 89337068, 984404373, 165994870, 907746571, 253065452, 820675989, 237940497, 835800944, 311676845, 762064596, 86155504, 987585937, 158416699, 915324742, 200145781, 873595660, 24913014, 1048828427, 146389750, 927351691, 33329164, 1040412277, 359173175, 714568266, 69633870, 1004107571, 293845527, 779895914, 408538962, 665202479, 108009024, 965732417, 298706058, 775035383, 400902748, 672838693, 363927126, 709814315, 314993962, 758747479, 37787724, 1035953717, 470658221, 603083220, 339662003, 734079438, 222916136, 850825305, 82110815, 991630626};

scalar bezout_coefficients_array[2*PARAM_R - 1] = {0, 503161189, 570580252, 664249193, 409492248, 75044049, 998697392, 57369835, 1016371606, 912550256, 161191185, 766080314, 307661127, 868356109, 205385332, 384716825, 689024616, 994854739, 78886702, 390483735, 683257706, 119188814, 954552627, 945066582, 128674859, 928354879, 145386562, 41916888, 1031824553, 1039336453, 34404988, 618745556, 454995885, 869940945, 203800496, 117767615, 955973826, 977912637, 95828804, 648483817, 425257624, 35966194, 1037775247, 5964040, 1067777401, 411365396, 662376045, 1050720985, 23020456, 556029539, 517711902, 984833865, 88907576, 563872963, 509868478, 782328051, 291413390, 396771399, 676970042, 658316556, 415424885, 79035084, 994706357, 357284133, 716457308, 1038924506, 34816935, 16664582, 1057076859, 1000546566, 73194875, 1019736929, 54004512, 924388412, 149353029, 869471960, 204269481, 683793484, 389947957, 962283373, 111458068, 495815313, 577926128, 772199831, 301541610, 367039719, 706701722, 873290067, 200451374, 181963563, 891777878, 1054847579, 18893862, 916244460, 157496981, 885060121, 188681320, 130702176, 943039265, 884537087, 189204354, 411168915, 662572526, 44668534, 1029072907, 973188234, 100553207, 980134911, 93606530, 310348383, 763393058, 1030663689, 43077752, 616079070, 457662371, 1061284934, 12456507, 436797830, 636943611, 655840969, 417900472, 381032298, 692709143, 990744006, 82997435, 947208715, 126532726};

#elif (PARAM_Q == 1073740609)
scalar cyclotomic_factorisation_array[2*PARAM_R - 1] = {1, 462868367, 610872242, 109389081, 964351528, 67138279, 1006602330, 210412476, 863328133, 147527534, 926213075, 154235627, 919504982, 205880521, 867860088, 436348658, 637391951, 314037703, 759702906, 315980002, 757760607, 350459901, 723280708, 114039180, 959701429, 119964082, 953776527, 101692609, 972048000, 12155428, 1061585181, 13210439, 1060530170, 330514882, 743225727, 290151346, 783589263, 89303044, 984437565, 173862666, 899877943, 13750973, 1059989636, 141044114, 932696495, 527549350, 546191259, 399263929, 674476680, 164505006, 909235603, 112561345, 961179264, 441703667, 632036942, 155994955, 917745654, 278495011, 795245598, 197155324, 876585285, 238371586, 835369023};

scalar bezout_coefficients_array[2*PARAM_R - 1] = {0, 305436121, 768304488, 570439444, 503301165, 591564845, 482175764, 613988118, 459752491, 433930044, 639810565, 105206238, 968534371, 999976842, 73763767, 1016721019, 57019590, 1013758568, 59982041, 6077714, 1067662895, 587716609, 486024000, 855566280, 218174329, 379851453, 693889156, 712100255, 361640354, 157990001, 915750608, 82252503, 991488106, 736502269, 337238340, 316018471, 757722138, 593150977, 480589632, 975162947, 98577662, 954554816, 119185793, 614867782, 458872827, 397622799, 676117810, 165257441, 908483168, 543475524, 530265085, 44651522, 1029089087, 928664936, 145075673, 70522057, 1003218552, 809965934, 263774675, 986809276, 86931333, 529994818, 543745791};

#elif (PARAM_Q == 1073739937)
scalar cyclotomic_factorisation_array[2*PARAM_R - 1] = {1, 319779557, 753960380, 133260203, 940479734, 287332869, 786407068, 355562456, 718177481, 313594583, 760145354, 311636094, 762103843, 512358091, 561381846, 460376567, 613363370, 311135302, 762604635, 159161129, 914578808, 279289587, 794450350, 223936993, 849802944, 271496178, 802243759, 135444559, 938295378, 244752693, 828987244};

scalar bezout_coefficients_array[2*PARAM_R - 1] = {0, 376980190, 696759747, 680536403, 393203534, 603500070, 470239867, 155818047, 917921890, 793049014, 280690923, 177781228, 895958709, 693667260, 380072677, 648838465, 424901472, 937991848, 135748089, 604592248, 469147689, 414493622, 659246315, 767058252, 306681685, 918172286, 155567651, 616450533, 457289404, 397225175, 676514762};

#elif(PARAM_Q == 134221313)
scalar cyclotomic_factorisation_array[2*PARAM_R - 1] = {0};

scalar bezout_coefficients_array[2*PARAM_R - 1] = {0};

#elif (PARAM_N == 512) && (PARAM_Q == 8000033)
scalar cyclotomic_factorisation_array[2*PARAM_R - 1] = {1, 228281, 7771752, 770814, 7229219, 1464899, 6535134, 1273894, 6726139, 3403369, 4596664, 1651039, 6348994, 3279263, 4720770};

scalar bezout_coefficients_array[2*PARAM_R - 1] = {0, 3885876, 4114157, 4732466, 3267567, 385407, 7614626, 4825536, 3174497, 5639648, 2360385, 636947, 7363086, 5701701, 2298332};

#elif (PARAM_Q == 1073741969)
scalar cyclotomic_factorisation_array[2*PARAM_R - 1] = {1, 443087604, 630654365, 264924626, 808817343, 56119300, 1017622669, 299249113, 774492856, 405878247, 667863722, 112073183, 961668786, 65537676, 1008204293};

scalar bezout_coefficients_array[2*PARAM_R - 1] = {0, 852198167, 221543802, 1045682319, 28059650, 941279656, 132462313, 1040973131, 32768838, 480834393, 592907576, 333931861, 739810108, 387246428, 686495541};

#elif (PARAM_Q == 134218289)
scalar cyclotomic_factorisation_array[2*PARAM_R - 1] = {1, 63431645, 70786644, 59190834, 75027455, 19132594, 115085695, 27887844, 106330445, 4441712, 129776577, 25881438, 108336851, 45284023, 88934266};

scalar bezout_coefficients_array[2*PARAM_R - 1] = {0, 35393322, 98824967, 9566297, 124651992, 29595417, 104622872, 12940719, 121277570, 44467133, 89751156, 13943922, 120274367, 131997433, 2220856};

#elif (PARAM_Q == 33554641)
scalar cyclotomic_factorisation_array[2*PARAM_R - 1] = {1, 14986550, 18568091, 15249757, 18304884, 15488453, 18066188, 14864901, 18689740, 8066139, 25488502, 9867322, 23687319, 318514, 33236127};

scalar bezout_coefficients_array[2*PARAM_R - 1] = {0, 26061366, 7493275, 24521547, 9033094, 24402199, 9152442, 4933661, 28620980, 33395384, 159257, 24209771, 9344870, 12744251, 20810390};

#elif (PARAM_Q == 8388593)
scalar cyclotomic_factorisation_array[2*PARAM_R - 1] = {1, 325432, 8063161, 1663925, 6724668, 2373857, 6014736, 4044329, 4344264, 610614, 7777979, 2095617, 6292976, 3390763, 4997830};

scalar bezout_coefficients_array[2*PARAM_R - 1] = {0, 8225877, 162716, 5381225, 3007368, 5026259, 3362334, 5242105, 3146488, 5889678, 2498915, 6216461, 2172132, 305307, 8083286};


#elif (PARAM_Q == 8381777)
scalar cyclotomic_factorisation_array[2*PARAM_R - 1] = {1, 3357658, 5024119, 1884573, 6497204, 3881900, 4499877, 1769246, 6612531, 2790443, 5591334, 2351358, 6030419, 408177, 7973600};

scalar bezout_coefficients_array[2*PARAM_R - 1] = {0, 6702948, 1678829, 6440827, 1940950, 3248602, 5133175, 3986800, 4394977, 1175679, 7206098, 5586110, 2795667, 7497154, 884623};

#elif (PARAM_Q == 17)
scalar cyclotomic_factorisation_array[2*PARAM_R - 1] = {1, 4, 13, 8, 9, 2, 15, 3, 14, 5, 12, 7, 10, 6, 11};

scalar bezout_coefficients_array[2*PARAM_R - 1] = {0, 15, 2, 16, 1, 13, 4, 3, 14, 12, 5, 11, 6, 10, 7};

#elif ((PARAM_N == 8) && (PARAM_Q == 8000053))
scalar cyclotomic_factorisation_array[2*PARAM_R - 1] = {1, 3193400, 4806653};

scalar bezout_coefficients_array[2*PARAM_R - 1] = {0, 6403353, 1596700};

#elif (PARAM_Q == 13)
scalar cyclotomic_factorisation_array[2*PARAM_R - 1] = {1, 5, 8};

scalar bezout_coefficients_array[2*PARAM_R - 1] = {0, 4, 9};

#elif (PARAM_Q == 7681)
scalar cyclotomic_factorisation_array[2*PARAM_R - 1] = {1, 3383, 4298, 1213, 6468, 1925, 5756, 849, 6832, 527, 7154, 1728, 5953, 583, 7098, 2648, 5033, 2138, 5543, 1366, 6315, 2784, 4897, 2132, 5549, 97, 7584, 2446, 5235, 2381, 5300, 1846, 5835, 365, 7316, 2753, 4928, 3654, 4027, 2645, 5036, 330, 7351, 2273, 5408, 878, 6803, 3000, 4681, 2399, 5282, 1112, 6569, 1794, 5887, 2268, 5413, 675, 7006, 3092, 4589, 1286, 6395, 732, 6949, 3074, 4607, 3477, 4204, 3080, 4601, 2469, 5212, 3380, 4301, 693, 6988, 1714, 5967, 1382, 6299, 2423, 5258, 1908, 5773, 2724, 4957, 1875, 5806, 1381, 6300, 695, 6986, 799, 6882, 2881, 4800, 766, 6915, 243, 7438, 202, 7479, 1080, 6601, 2516, 5165, 3411, 4270, 2551, 5130, 3449, 4232, 528, 7153, 2508, 5173, 2941, 4740, 1655, 6026, 584, 7097, 2774, 4907, 1740, 5941, 1402, 6279, 3789, 3892, 2819, 4862, 3125, 4556, 3180, 4501, 3141, 4540, 257, 7424, 1478, 6203, 1155, 6526, 2264, 5417, 3566, 4115, 3073, 4608, 2573, 5108, 1886, 5795, 1220, 6461, 2563, 5118, 1800, 5881, 1633, 6048, 1996, 5685, 869, 6812, 3837, 3844, 319, 7362, 2897, 4784, 405, 7276, 3501, 4180, 219, 7462, 880, 6801, 3188, 4493, 2900, 4781, 2063, 5618, 1587, 6094, 198, 7483, 2722, 4959, 993, 6688, 1044, 6637, 1408, 6273, 3041, 4640, 2844, 4837, 1003, 6678, 1853, 5828, 2880, 4801, 3532, 4149, 1415, 6266, 1682, 5999, 2562, 5119, 3078, 4603, 648, 7033, 3099, 4582, 1591, 6090, 2028, 5653, 2044, 5637, 1952, 5729, 1097, 6584, 1228, 6453, 1848, 5833, 550, 7131, 2681, 5000, 1438, 6243, 2990, 4691, 707, 6974, 417, 7264, 2593, 5088, 1125, 6556, 3780, 3901, 3139, 4542, 3586, 4095, 2372, 5309, 2169, 5512, 1959, 5722, 1406, 6275, 2838, 4843, 296, 7385, 346, 7335, 3006, 4675, 2757, 4924, 2197, 5484, 1230, 6451, 2012, 5669, 2002, 5679, 1876, 5805, 3081, 4600, 94, 7587, 3394, 4287, 1193, 6488, 1131, 6550, 1035, 6646, 2996, 4685, 3452, 4229, 1266, 6415, 3120, 4561, 2173, 5508, 542, 7139, 1437, 6244, 702, 6979, 1065, 6616, 506, 7175, 1683, 5998, 1968, 5713, 1667, 6014, 1607, 6074, 3626, 4055, 201, 7480, 1979, 5702, 2875, 4806, 62, 7619, 2359, 5322, 1604, 6077, 3546, 4135, 398, 7283, 2259, 5422, 1129, 6552, 1950, 5731, 3694, 3987, 185, 7496, 2799, 4882, 1656, 6025, 2358, 5323, 3445, 4236, 2922, 4759, 321, 7360, 2583, 5098, 2689, 4992, 2668, 5013, 669, 7012, 763, 6918, 413, 7268, 3799, 3882, 1704, 5977, 1994, 5687, 1784, 5897, 793, 6888, 2050, 5631, 3086, 4595, 1459, 6222, 2671, 5010, 3137, 4544, 111, 7570, 856, 6825, 1393, 6288, 3615, 4066, 3265, 4416, 217, 7464, 2951, 4730, 2067, 5614, 2110, 5571, 2481, 5200, 1499, 6182, 1657, 6024, 1717, 5964, 1775, 5906, 2395, 5286, 1170, 6511, 2433, 5248, 3193, 4488, 1885, 5796, 1725, 5956, 2546, 5135, 2717, 4964, 572, 7109, 536, 7145, 3751, 3930, 621, 7060, 2811, 4870, 535, 7146, 1036, 6645, 2252, 5429, 2760, 4921, 3016, 4665, 1115, 6566, 674, 7007, 3376, 4305, 639, 7042, 3832, 3849, 1872, 5809, 1211, 6470, 2840, 4841, 3280, 4401, 2805, 4876, 118, 7563, 218, 7463, 335, 7346, 3483, 4198, 738, 6943, 329, 7352, 2457, 5224, 1189, 6492, 113, 7568, 1771, 5910, 3239, 4442, 3250, 4431, 1897, 5784, 3765, 3916};

scalar bezout_coefficients_array[2*PARAM_R - 1] = {0, 2149, 5532, 4803, 2878, 4447, 3234, 864, 6817, 4132, 3549, 4265, 3416, 4104, 3577, 1066, 6615, 3792, 3889, 1223, 6458, 5031, 2650, 1324, 6357, 6612, 1069, 683, 6998, 1392, 6289, 6181, 1500, 5040, 2641, 6784, 897, 7125, 556, 6547, 1134, 4178, 3503, 1546, 6135, 643, 7038, 6758, 923, 4023, 3658, 5854, 1827, 2464, 5217, 2518, 5163, 165, 7516, 4977, 2704, 439, 7242, 7298, 383, 5281, 2400, 3962, 3719, 7580, 101, 2565, 5116, 2135, 5546, 1258, 6423, 7141, 540, 264, 7417, 2116, 5565, 1254, 6427, 5311, 2370, 4668, 3013, 7389, 292, 6294, 1387, 6811, 870, 1537, 6144, 7315, 366, 5579, 2102, 6141, 1540, 6824, 857, 4187, 3494, 5991, 1690, 2606, 5075, 2629, 5052, 691, 6990, 954, 6727, 1362, 6319, 4778, 2903, 3150, 4531, 3493, 4188, 3441, 4240, 6977, 704, 522, 7159, 6320, 1361, 4337, 3344, 5361, 2320, 1422, 6259, 2914, 4767, 4342, 3339, 2291, 5390, 324, 7357, 6142, 1539, 1281, 6400, 4548, 3133, 6840, 841, 5915, 1766, 6241, 1440, 1022, 6659, 976, 6705, 1014, 6667, 3045, 4636, 4389, 3292, 614, 7067, 6757, 924, 275, 7406, 2500, 5181, 719, 6962, 3487, 4194, 6186, 1495, 5137, 2544, 3632, 4049, 5791, 1890, 4403, 3278, 2431, 5250, 5403, 2278, 5735, 1946, 6980, 701, 1590, 6091, 5411, 2270, 739, 6942, 3712, 3969, 2559, 5122, 7071, 610, 5127, 2554, 6738, 943, 5377, 2304, 5898, 1783, 1132, 6549, 3263, 4418, 3406, 4275, 998, 6683, 900, 6781, 4657, 3024, 5759, 1922, 4000, 3681, 2392, 5289, 4043, 3638, 2090, 5591, 3950, 3731, 6087, 1594, 7241, 440, 2809, 4872, 1450, 6231, 99, 7582, 3047, 4634, 5409, 2272, 2505, 5176, 6138, 1543, 3111, 4570, 892, 6789, 6684, 997, 4237, 3444, 1025, 6656, 3896, 3785, 428, 7253, 4537, 3144, 5648, 2033, 2807, 4874, 2365, 5316, 2208, 5473, 3732, 3949, 268, 7413, 286, 7395, 6408, 1273, 2482, 5199, 4703, 2978, 4783, 2898, 2624, 5057, 5437, 2244, 2982, 4699, 2953, 4728, 7096, 585, 5038, 2643, 3012, 4669, 4590, 3091, 5081, 2600, 1055, 6626, 518, 7163, 6555, 1126, 1380, 6301, 1508, 6173, 2435, 5246, 4108, 3573, 3530, 4151, 5716, 1965, 4398, 3283, 7344, 337, 5993, 1688, 3521, 4160, 6745, 936, 1916, 5765, 4446, 3235, 6261, 1420, 5243, 2438, 1640, 6041, 7622, 59, 109, 7572, 4005, 3676, 7312, 369, 5582, 2099, 4008, 3673, 3897, 3784, 2955, 4726, 3246, 4435, 2612, 5069, 5723, 1958, 2892, 4789, 2221, 5460, 6056, 1625, 6978, 703, 4820, 2861, 1419, 6262, 148, 7533, 6495, 1186, 2756, 4925, 5888, 1793, 5410, 2271, 173, 7508, 1503, 6178, 5219, 2462, 4939, 2742, 6680, 1001, 6743, 938, 6675, 1006, 7066, 615, 253, 7428, 4373, 3308, 351, 7330, 3122, 4559, 7048, 633, 6121, 1560, 271, 7410, 2754, 4927, 6183, 1498, 1726, 5955, 4358, 3323, 4406, 3275, 5984, 1697, 3244, 4437, 47, 7634, 5381, 2300, 5278, 2403, 2851, 4830, 5868, 1813, 3941, 3740, 4682, 2999, 6697, 984, 4674, 3007, 4644, 3037, 31, 7650, 2661, 5020, 6879, 802, 5908, 1773, 4970, 2711, 7482, 199, 4405, 3276, 6706, 975, 3933, 3748, 1847, 5834, 2441, 5240, 828, 6853, 4001, 3680, 1461, 6220, 2118, 5563, 1179, 6502, 3506, 4175, 6347, 1334, 5132, 2549, 2496, 5185, 1941, 5740, 6829, 852, 3634, 4047, 4222, 3459};

#elif (PARAM_Q == 1032193)
scalar cyclotomic_factorisation_array[2*PARAM_R - 1] = {1, 130048, 902145, 355573, 676620, 343297, 688896, 354660, 677533, 311668, 720525, 356512, 675681, 372598, 659595, 37421, 994772, 263787, 768406, 157041, 875152, 102730, 929463, 126714, 905479, 58973, 973220, 179521, 852672, 205734, 826459, 243133, 789060, 207785, 824408, 205494, 826699, 425251, 606942, 392728, 639465, 450889, 581304, 146560, 885633, 391135, 641058, 380351, 651842, 166095, 866098, 490491, 541702, 89446, 942747, 46876, 985317, 1810, 1030383, 501302, 530891, 12616, 1019577, 298302, 733891, 363216, 668977, 15634, 1016559, 249778, 782415, 500160, 532033, 133592, 898601, 365641, 666552, 186356, 845837, 98547, 933646, 131968, 900225, 235533, 796660, 268309, 763884, 11488, 1020705, 408153, 624040, 213517, 818676, 434923, 597270, 428329, 603864, 2354, 1029839, 89681, 942512, 85981, 946212, 174497, 857696, 222751, 809442, 268358, 763835, 56339, 975854, 352629, 679564, 425588, 606605, 293258, 738935, 149420, 882773, 19021, 1013172, 508580, 523613, 199319, 832874, 425497, 606696, 161988, 870205, 188487, 843706, 125338, 906855, 435632, 596561, 166107, 866086, 148032, 884161, 48658, 983535, 499699, 532494, 361023, 671170, 11694, 1020499, 397258, 634935, 316541, 715652, 27891, 1004302, 42566, 989627, 23801, 1008392, 274359, 757834, 285507, 746686, 432260, 599933, 254122, 778071, 334575, 697618, 220680, 811513, 101532, 930661, 55468, 976725, 494413, 537780, 370015, 662178, 94747, 937446, 272296, 759897, 104957, 927236, 509913, 522280, 73461, 958732, 72095, 960098, 401541, 630652, 404938, 627255, 77643, 954550, 311732, 720461, 289132, 743061, 353303, 678890, 341535, 690658, 21526, 1010667, 105832, 926361, 510079, 522114, 161546, 870647, 144592, 887601, 440535, 591658, 32491, 999702, 408574, 623619, 413906, 618287, 185269, 846924, 26291, 1005902, 468752, 563441, 125835, 906358, 202258, 829935, 465502, 566691, 515354, 516839, 287641, 744552, 462448, 569745, 431978, 600215, 261274, 770919, 481230, 550963, 105257, 926936, 155582, 876611, 80750, 951443, 374651, 657542, 7069, 1025124, 310574, 721619, 184538, 847655, 221264, 810929, 496411, 535782, 14471, 1017722, 236769, 795424, 14778, 1017415, 94022, 938171, 503526, 528667, 225328, 806865, 332102, 700091, 181390, 850803, 382034, 650159, 211963, 820230, 247910, 784283, 348675, 683518, 127445, 904748, 44359, 987834, 368294, 663899, 78526, 953667, 431827, 600366, 286855, 745338, 322766, 709427, 87770, 944423, 387685, 644508, 191795, 840398, 132025, 900168, 88838, 943355, 3044, 1029149, 496000, 536193, 406245, 625948, 416752, 615441, 307150, 725043, 438486, 593707, 29994, 1002199, 2365, 1029828, 402009, 630184, 109018, 923175, 150801, 881392, 298552, 733641, 391095, 641098, 187515, 844678, 411683, 620510, 267933, 764260, 398007, 634186, 335842, 696351, 342640, 689553, 125090, 907103, 171980, 860213, 97116, 935077, 202448, 829745, 189347, 842846, 169593, 862600, 362633, 669560, 78044, 954149, 87657, 944536, 417089, 615104, 151878, 880315, 96757, 935436, 410527, 621666, 490012, 542181, 450858, 581335, 59775, 972418, 173717, 858476, 496985, 535208, 108392, 923801, 90774, 941419, 214189, 818004, 192830, 839363, 26905, 1005288, 342821, 689372, 326841, 705352, 169833, 862360, 423830, 608363, 262204, 769989, 422156, 610037, 187924, 844269, 93309, 938884, 381458, 650735, 377789, 654404, 486324, 545869, 98137, 934056, 390362, 641831, 481250, 550943, 195784, 836409, 212901, 819292, 219540, 812653, 279540, 752653, 141066, 891127, 184979, 847214, 158017, 874176, 135621, 896572, 72850, 959343, 502747, 529446, 422478, 609715, 182253, 849940, 501974, 530219, 331533, 700660, 514735, 517458, 476844, 555349, 104659, 927534, 196734, 835459, 485414, 546779, 260378, 771815, 511353, 520840, 368726, 663467, 458933, 573260, 144862, 887331, 304197, 727996, 382538, 649655, 144880, 887313, 296782, 735411, 370597, 661596, 243100, 789093, 399929, 632264, 174292, 857901, 490722, 541471, 18045, 1014148, 203097, 829096, 428021, 604172, 435060, 597133, 55778, 976415, 437701, 594492, 207723, 824470, 246635, 785558, 23198, 1008995, 328191, 704002, 434811, 597382};

scalar bezout_coefficients_array[2*PARAM_R - 1] = {0, 967169, 65024, 687745, 344448, 693883, 338310, 853937, 178256, 845894, 186299, 854863, 177330, 876359, 155834, 486610, 545583, 968836, 63357, 102867, 929326, 426336, 605857, 384203, 647990, 497386, 534807, 980828, 51365, 594617, 437576, 44723, 987470, 761342, 270851, 433049, 599144, 706272, 325921, 250651, 781542, 6308, 1025885, 1031288, 905, 23438, 1008755, 728722, 303471, 929446, 102747, 619989, 412204, 637663, 394530, 711664, 320529, 958913, 73280, 196364, 835829, 741541, 290652, 134179, 898014, 487927, 544266, 428848, 603345, 404721, 627472, 559087, 473106, 560937, 471256, 301932, 730261, 1031016, 1177, 506586, 525607, 254290, 777903, 615756, 416437, 728845, 303348, 74710, 957483, 146629, 885564, 692411, 339782, 819399, 212794, 333276, 698917, 939015, 93178, 66796, 965397, 250080, 782113, 1024376, 7817, 907304, 124889, 149151, 883042, 850585, 181608, 622855, 409338, 298635, 733558, 720173, 312020, 5744, 1026449, 466823, 565370, 65984, 966209, 633863, 398330, 650251, 381942, 532342, 499851, 204287, 827906, 608731, 423462, 206953, 825240, 80773, 951420, 261057, 771136, 295829, 736364, 72296, 959897, 10763, 1021430, 52916, 979277, 339445, 692748, 686864, 345329, 477275, 554918, 202469, 829724, 144566, 887627, 155866, 876327, 130637, 901556, 215989, 816204, 791578, 240615, 568725, 463468, 77791, 954402, 40375, 991818, 328771, 703422, 512562, 519631, 372276, 659917, 800969, 231224, 799442, 232751, 774516, 257677, 502951, 529242, 797817, 234376, 931064, 101129, 453179, 579014, 21283, 1010910, 502151, 530042, 653276, 378917, 527997, 504196, 833564, 198629, 674367, 357826, 696608, 335585, 5847, 1026346, 958177, 74016, 599150, 433043, 765946, 266247, 1007864, 24329, 80994, 951199, 610340, 421853, 217816, 814377, 62669, 969524, 261140, 771053, 479366, 552827, 315326, 716867, 480049, 552144, 896045, 136148, 463618, 568575, 331089, 701104, 468723, 563470, 216130, 816063, 658850, 373343, 905132, 127061, 683384, 348809, 110340, 921853, 50766, 981427, 1004459, 27734, 268890, 763303, 961660, 70533, 423607, 608586, 437088, 595105, 448286, 583907, 820954, 211239, 424970, 607223, 36425, 995768, 264723, 767470, 195181, 837012, 791568, 240625, 467028, 565165, 243162, 789031, 139770, 892423, 109770, 922423, 97892, 934301, 622547, 409646, 93962, 938231, 562751, 469442, 190729, 841464, 704991, 327202, 211915, 820278, 601013, 431180, 131102, 901091, 821115, 211078, 623191, 409002, 986806, 45387, 764589, 267604, 977997, 54196, 687507, 344686, 352676, 679517, 96415, 935778, 529549, 502644, 316132, 716061, 945047, 87146, 701395, 330798, 910643, 121550, 525119, 507074, 786832, 245361, 617645, 414548, 302086, 730107, 217530, 814663, 27889, 1004304, 734947, 297246, 412235, 619958, 11599, 1020594, 392779, 639414, 298691, 733502, 680192, 352001, 840924, 191269, 363998, 668195, 883802, 148391, 959753, 72440, 847830, 184363, 260420, 771773, 959762, 72431, 286630, 745563, 350330, 681863, 781206, 250987, 238422, 793771, 258729, 773464, 902004, 130189, 242707, 789486, 98367, 933826, 568426, 463767, 582109, 450084, 987774, 44419, 420199, 611994, 709939, 322254, 719219, 312974, 208376, 823817, 1522, 1030671, 784193, 248000, 659524, 372669, 732010, 300183, 161383, 870810, 988308, 43885, 579819, 452374, 538276, 493917, 184147, 848046, 39263, 992930, 90695, 941498, 866142, 166051, 251763, 780430, 919529, 112664, 123955, 908238, 690434, 341759, 622078, 410115, 191017, 841176, 876906, 155287, 939924, 92269, 921561, 110632, 267891, 764302, 7389, 1024804, 985182, 47011, 508861, 523332, 397712, 634481, 956254, 75939, 307552, 724641, 310833, 721360, 564475, 467718, 602955, 429238, 545984, 486209, 225429, 806764, 787187, 245006, 334780, 697413, 431300, 600893, 472268, 559925, 993171, 39022, 983635, 48558, 946203, 85990, 421423, 610770, 930969, 101224, 14997, 1017196, 514914, 517279, 878618, 153575, 812950, 219243, 977684, 54509, 717101, 315092, 591497, 440696, 882917, 149276, 711644, 320549, 609854, 422339, 721938, 310255, 382130, 650063, 864272, 167921, 715100, 317093, 62545, 969648, 860873, 171320};

#else
#error "The CRT trees are not defined for this value of PARAM_Q"

#endif
