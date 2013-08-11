#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Last-Updated : <2013/08/11 14:42:04 by ymnk>
import sys,os

# 1D map plot
def mkdir(folder):
    folder = os.path.abspath(folder)
    if not os.path.isdir(folder):
        os.system("mkdir {0}".format(folder))


if __name__  == "__main__":
    data_folder =  os.path.abspath(sys.argv[1])
    map_folder = os.path.abspath(os.path.join(sys.argv[1],"../map"))
    mkdir(map_folder)
    
    dump_file = os.path.join(map_folder,"dump.txt")
    png_file =  os.path.join(map_folder,"-".join([len for len in (sys.argv[1].rstrip().split("/")) if not len == ""][-2:])+".png")

    dump = open(dump_file,"w")
    for i in range(0,18000,10):
        data_file = os.path.join(data_folder,"data{0}.txt".format(i))
        data = open(data_file,"r")
        for j,line in enumerate(data.readlines()):
            line = line.rstrip().split(" ")
            dump.write("{0} {1} {2}\n".format(i,line[0],line[2]))
        dump.write("\n")
        data.close()
    dump.close()


    cmd_file = os.path.join(data_folder,"gnuplot.txt")
    gnuplot_file = open(cmd_file,"wt")
    gnuplot_file.write("set term png\n")
    gnuplot_file.write("set pm3d map\n")
    gnuplot_file.write("set output \"{0}\"\n".format(png_file))
    gnuplot_file.write("splot \"{0}\" u 1:2:3 \n".format(dump_file))
    gnuplot_file.write("set output\n")
    gnuplot_file.close()
    
    os.system("gnuplot {0}".format(cmd_file))
    os.system("rm {0}".format(cmd_file))
