# import sys
# sys.path.append('D:\git\Andre\HGS-Tools\Python\pks_curve')
# import read_mprops

import read_mprops

param={"grok_dirc":r"./test_data",
		"grok_name":"abdul",
		"mprop_dirc":r".\test_data\mprop",
		"mprop_name":"gen.mprops",
		}

# read_mprops.read_mprops(grok_dirc = param["grok_dirc"], grok_name = param["grok_name"], mprop_dirc = param["mprop_dirc"], mprop_name = param["mprop_name"], 
#                 lgentab = True, linctab = False, ldebug = False,
#                 sr = 0.13, kr_min = 1.072e-12, tsf = 1e-4, p_min = -1.0e3)

read_mprops.read_mprops(grok_dirc = param["grok_dirc"], grok_name = param["grok_name"], mprop_dirc = param["mprop_dirc"], mprop_name = param["mprop_name"], 
                lgentab = False, linctab = True, ldebug = True,
                sr = 0.13, kr_min = 1.072e-12, tsf = 1e-4, p_min = -1.0e3 )