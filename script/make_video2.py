#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Last-Updated : <2013/08/10 04:00:44 by samui>
import sys,os

# 1D plot

def mkdir(folder):
    folder = os.path.abspath(folder)
    if not os.path.isdir(folder):
        os.system("mkdir {0}".format(folder))


if __name__  == "__main__":
    data_folder =  os.path.abspath(sys.argv[1])
    movie_folder = os.path.abspath(os.path.join(sys.argv[1],"../movie"))
    
    mkdir(movie_folder)
    mkdir(os.path.join(data_folder,"png"))

    cmd_file = os.path.join(data_folder,"gnuplot.txt")
    
    movie_file = os.path.join(movie_folder,"-".join((sys.argv[1].split("/"))[1:])+".avi")

    gnuplot_file = open(cmd_file,"wt")
    gnuplot_file.write("set term png\n")
    gnuplot_file.write("set yrange[0:1.0]\n")
    gnuplot_file.write("set view 35,45\n")
    for fname in os.listdir(data_folder):
        root, ext = os.path.splitext(fname)
        if not os.path.isfile(os.path.join(data_folder,fname)) or fname in ["gnuplot.txt",".DS_Store","parameter.txt"]:
            continue
        gnuplot_file.write("set output \"{0}\"\n".format(os.path.join(data_folder,"png","data{0:0>5}.png".format(int(root[4:])))))

        gnuplot_file.write("splot \"./{0}\" u 1:2:3 w pm3d\n".format(fname))
        gnuplot_file.write("set output \n")    
    gnuplot_file.close()

    os.system("cd {0}; gnuplot {1}".format(sys.argv[1],cmd_file))
 
    os.system("ffmpeg -r 15 -i {0}/data00%3d.png -vcodec mjpeg {1}".format(os.path.join(data_folder,"png"),movie_file))
    os.system("rm -r {0}".format(os.path.join(data_folder,"png")))
    
    #os.system("rm {0}".format(cmd_file))
