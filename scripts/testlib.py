import numpy, yt, pylab, os, pandas, sys, math

def integrate(x,y):
    ret = 0
    for i in range(len(x)-1): ret = ret + 0.5*(y[i+1]+y[i]) * (x[i+1]-x[i])
    return ret

def readHeader(path):
    ret = dict()
    f = open(path+"/Header")
    f.readline().replace('\n','') # HyperClaw
    nvars = int(f.readline().replace('\n','')) # number of variables
    ret["vars"] = list()
    for i in range(nvars):
        ret["vars"].append(f.readline().replace('\n',''))
    ret["dim"] = int(f.readline().replace('\n',''))
    ret["time"] = float(f.readline().replace('\n',''))
    ret["amrlevels"] = int(f.readline().replace('\n',''))
    ret["geom_lo"] = [float(x) for x in f.readline().replace('\n','').split()]
    ret["geom_hi"] = [float(x) for x in f.readline().replace('\n','').split()]
    return ret


def validate(path, 
             outdir,
             vars = [],
             start = [0,0,0],
             end = [1,1,1],
             axis = 2,
             intercept = 0,
             generate_ref_data = False,
             reference=None,
             tolerance=1E-8,
             coord = 'x',
             tight_layout = True,
             filename_suffix=None):
            
    info = readHeader(path)
            
    ds = yt.load(path)
    dim = int(ds.domain_dimensions[0] > 1) + int(ds.domain_dimensions[1] > 1) + int(ds.domain_dimensions[2] > 1)

    if not type(start) == list: raise Exception("Invalid start, ",start)
    if not type(end)   == list: raise Exception("Invalid end, ",end)
    
    if len(start) == 2: start.append(0.0)
    if len(end) == 2:   end.append(0.0)
    if dim == 2: 
        start[2] = 0
        end[2] = 0

    prof = ds.ray(start,end)
    slice2d = ds.slice(axis,intercept)

    new_df = prof.to_dataframe([("gas","x"),("gas","y"),("gas","z"),*vars])
    new_df = new_df.sort_values(by=coord)
    
    if not reference:
        reference = "reference/reference-{}d.csv".format(dim)

    if generate_ref_data:
        new_df.to_csv(reference)
        print("Reference data generated, now exiting noisily so we don't accidentally pass the test")
        exit(-1)
    

    all_ok = True
    if not type(tolerance)==list:
        tolerance = [tolerance] * len(vars)
    elif len(tolerance) != len(vars):
        raise Exception("Wrong number of tolerance values")
    for i, (var, tol) in enumerate(zip(vars,tolerance)):
        if isinstance(reference,list):
            ref_df = pandas.read_csv(reference[i])
        else:
            ref_df = pandas.read_csv(reference)
        ref_df = ref_df.sort_values(by=coord)
        
        len_x = info["geom_hi"][0]-info["geom_lo"][0]
        len_y = info["geom_hi"][1]-info["geom_lo"][1]
        data2d = slice2d.to_frb(width=len_x,height=len_y,resolution=(1000,1000*len_y/len_x))[var]
        pylab.clf()
        im = pylab.imshow(data2d,origin='lower',cmap='jet',
                     extent=[info["geom_lo"][0],info["geom_hi"][0],info["geom_lo"][1],info["geom_hi"][1]])
        pylab.colorbar(im)
        pylab.plot([start[0],end[0]],[start[1],end[1]],linestyle='-',color='white')
        if tight_layout: pylab.tight_layout()
        pylab.savefig("{}/2d_{}.png".format(outdir,var))

        new_x,new_y,new_var = new_df["x"], new_df["y"], new_df[var]
        ref_x,ref_y,ref_var = ref_df["x"], ref_df["y"], ref_df[var]


        if coord == 'x':
            new_coord = new_x
            ref_coord = ref_x
        elif coord == 'y':
            new_coord = new_y
            ref_coord = ref_y
            
        #print(type(ref_coord),axis=None)
        final_coord = sorted(numpy.concatenate((ref_coord, new_coord),axis=None))
        ref_final_var   = numpy.interp(final_coord, ref_coord, ref_var)
        new_final_var   = numpy.interp(final_coord, new_coord, new_var)
            
        pylab.clf()
        pylab.plot(final_coord,ref_final_var,color='C0',label='ref')
        pylab.plot(final_coord,new_final_var,color='C1',label='new',linestyle='--',marker='o',markerfacecolor='None')
        pylab.xlabel(coord)
        pylab.ylabel(var)
        pylab.legend()
        pylab.grid()
        pylab.savefig(outdir+"/{}.png".format(var))
        err = numpy.sqrt(integrate(final_coord, (    new_final_var  -     ref_final_var )**2))
        mag = numpy.sqrt(integrate(final_coord, (abs(new_final_var) + abs(ref_final_var))**2))



        relerr = err/mag

        print("{} abs error [tolerance={}]".format(var,tol),err)
        print("{} norm".format(var),mag)
        print("{} rel error [tolerance={}]".format(var,tol),relerr)
        
        if relerr > tol: all_ok = False
        if math.isnan(relerr): all_ok = False

    if not all_ok:
        print("ERROR: Regression test failed")
        raise(Exception("One or more errors out of tolerance or nan"))
