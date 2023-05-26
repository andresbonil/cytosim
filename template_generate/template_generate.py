#!/usr/bin/env python

from string import Template 
import random 
import sys, getopt


# TODO discuss/add randomized parameters/ranges


'''
parameters:

-k --kinesinbindingrate
-b --blobquantity
-t --thresholdbreaking
-s --stiffnessgrafted

'''

argumentList = sys.argv[1:]
print(argumentList)

def main(argv):
    replacements = {
        "kinesin_binding_rate": 10,
        "blob_motor_quantity": 5,
        "breaking_threshold": 30,
        "stiffness_grafted": 40
    }
   
    opts, args = getopt.getopt(argv,"k:b:t:s:h", ["kinesinbindingrate", "blobmotorquantity", "thresholdbreaking", "stiffnessgrafted, help" ])

    for opt, arg in opts:

        if opt == "-h":
            print("Please select which parameter you would like to modify in the config file and specify a value:\n-k <value> --kinesinbindingrate <value> \n-b <value>--blobquantity <value>\n-t <value> --thresholdbreaking <value>\n-s <value> --stiffnessgrafted <value>\n-h --help")
            sys.exit()
        elif opt in ("-k", "--kinesinbindingrate"):
            replacements["kinesin_binding_rate"] = arg
        elif opt in ("-b", "--blobmotorquantity"):
            replacements["blob_motor_quantity"] = arg
        elif opt in ("-t", "--thresholdbreaking"):
            replacements["breaking_threshold"] = arg
        elif opt in ("-s", "--stiffnessgrafted"):
            replacements["stiffnessgrafted"] = arg


    with open('configtemplate.txt', 'r') as file:
        src = Template(file.read())
        print(src)
        result = src.substitute(replacements)
        print(result)
        config = open('config.cym', 'w')
        config.write(result)
        config.close()
    file.close()


if __name__ == "__main__":
    main(sys.argv[1:])