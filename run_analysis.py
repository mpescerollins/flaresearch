import os, sys

import numpy as np
from threeML import *
from threeML.io.package_data import get_path_of_data_file
from threeML.utils.statistics.stats_tools import Significance

import likelihhod_analyis as la




if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--goes_t0", type=str, help="GOES flare trigger time ")
    parser.add_argument("--stix_t0", type=str, help="STIX flare trigger time ")
    parser.add_argument("--t0_list", type=str, help="list of trigger times ")
    args = parser.parse_args()
    goes_t0 = args.goes_t0
    stix_t0 = args.stix_t0
    t0_list = args.t0_list  
    if t0_list is not None:
        outputfile = 'flare_search_results.csv'
        outstr = 'FlareName, TS, inFoV\n'
        with open(t0_list,'r') as f:
            lines = f.readlines()
            for line in lines:
             
                t0 = line.split(',')[1].strip()
                
                dt_t0 = la.str2datetime(t0)
                flare_name = la.build_source_name(dt_t0)
                print("Analyzing flare at GOES time %s" % t0)
                ts = la.doAnalysis(t0)
                if ts==-1.0:
                    outstr +='%s, %.2f, False\n' % (flare_name,ts)
                else:
                    outstr +='%s, %.2f, True\n' % (flare_name,ts)
        with open(outputfile,'w') as f:
            f.write(outstr)
        print("Results written to %s" % outputfile)
    else:
        if goes_t0 is not None:
            print("Analyzing flare at GOES time %s" % goes_t0)
            la.doAnalysis(goes_t0)
        if stix_t0 is not None:
            print("Analyzing flare at STIX time %s" % stix_t0)
            la.doAnalysis(stix_t0)