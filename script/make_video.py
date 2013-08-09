#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Last-Updated : <2013/08/10 03:54:24 by samui>
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

    files = open(os.path.join(data_folder,"./data0.txt"))
    fact = len(files.readline().rstrip().split(" "))-1
    files.close()
    plot_label = ["u","v","w"]

    cmd_file = os.path.join(data_folder,"gnuplot.txt")
    
    movie_file = os.path.join(movie_folder,"-".join((sys.argv[1].split("/"))[1:])+".avi")

    gnuplot_file = open(cmd_file,"wt")
    gnuplot_file.write("set term png\n")
    gnuplot_file.write("set yrange[0:1.0]\n")

    for fname in os.listdir(data_folder):
        root, ext = os.path.splitext(fname)
        if not os.path.isfile(os.path.join(data_folder,fname)) or fname in ["gnuplot.txt",".DS_Store","parameter.txt"]:
            continue
        gnuplot_file.write("set output \"{0}\"\n".format(os.path.join(data_folder,"png","data{0:0>5}.png".format(int(root[4:])))))
        plot = ""
        for i in range(fact):
            plot += " \"./{0}\" u 1:{1} title \"{2}\",".format(fname,i+2,plot_label[i])
        gnuplot_file.write("plot {0}\n".format(plot[:-1]))
        gnuplot_file.write("set output \n")    
    gnuplot_file.close()

    os.system("cd {0}; gnuplot {1}".format(sys.argv[1],cmd_file))
 
    os.system("ffmpeg -r 15 -i {0}/data00%3d.png -vcodec mjpeg {1}".format(os.path.join(data_folder,"png"),movie_file))
    os.system("rm -r {0}".format(os.path.join(data_folder,"png")))
    
    #os.system("rm {0}".format(cmd_file))
