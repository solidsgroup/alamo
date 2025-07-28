#!/usr/bin/env python3
import argparse
import glob
import yt
import testlib
import pylab
import os
yt.set_log_level(50)


parser = argparse.ArgumentParser(description='Generate thumbnails');
parser.add_argument('outputs', default=None, nargs='*', help='Output directories')
parser.add_argument('--vars',default=None,nargs='*',help='Specific variables to plot')
args = parser.parse_args()

for output in args.outputs:
    try:
        path = sorted(glob.glob(output+"/*cell"))[-1]
        info = testlib.readHeader(path)
        len_x = info["geom_hi"][0]-info["geom_lo"][0]
        len_y = info["geom_hi"][1]-info["geom_lo"][1]
        resolution=(1000,1000*len_y/len_x)
        ds = yt.load(path)
        slice2d = ds.slice(2,0.0)

        print("Plotting: ", output)

        header_time = os.path.getmtime(f"{path}/Header")
        
        for var in info['vars']:
            if len(args.vars):
                if var not in args.vars:
                    continue
            if os.path.isfile(f"{output}/{var}.png"):
                png_time = os.path.getmtime(f"{output}/{var}.png")
                if header_time < png_time:
                    print("     Skipping {output}/{var}.png (most current)")
                    continue

            data2d = slice2d.to_frb(width=len_x,height=len_y,resolution=(1000,1000*len_y/len_x))[var]

            pylab.clf()
            im = pylab.imshow(data2d,origin='lower',cmap='jet',
                              extent=[info["geom_lo"][0],info["geom_hi"][0],info["geom_lo"][1],info["geom_hi"][1]])
            pylab.colorbar(im)
            pylab.tight_layout()
            pylab.savefig(f"{output}/{var}.png")
            print(f"     {output}/{var}.png")


    except Exception as e:
        print("Skipping ", output, ": ", e)

    
