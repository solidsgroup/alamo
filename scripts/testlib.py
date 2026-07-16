import numpy, yt, pylab, os, pandas, sys, math

yt.set_log_level(50) # disable logging

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

def warning(msg):
    print(msg, file=sys.stderr)

def rayCoordinate(df,start,end):
    start = list(start)
    end = list(end)
    if len(start) == 2: start.append(0.0)
    if len(end) == 2:   end.append(0.0)
    start = numpy.array(start[:3], dtype=float)
    end = numpy.array(end[:3], dtype=float)
    direction = end - start
    length = numpy.linalg.norm(direction)
    if length <= 0.0:
        raise Exception("Cannot construct ray coordinate from identical start and end points")

    points = numpy.column_stack([
        df[c].to_numpy(dtype=float) if c in df else numpy.zeros(len(df))
        for c in ["x","y","z"]
    ])
    return numpy.dot(points - start, direction) / length

def hasDuplicates(values):
    return pandas.Series(values).duplicated().any()

def sortForComparison(df,start,end,coord,use_ray_coord):
    ret = df.copy()
    if use_ray_coord:
        compare_coord = "_ray_coord"
        ret[compare_coord] = rayCoordinate(ret,start,end)
    else:
        compare_coord = coord

    if compare_coord not in ret:
        raise Exception("{} is not available for comparison".format(compare_coord))
    if hasDuplicates(ret[compare_coord]):
        raise Exception("{} contains duplicate values and cannot be used for interpolation".format(compare_coord))

    return ret.sort_values(by=compare_coord).reset_index(drop=True), compare_coord
 
def readContours(path,
                 start,
                 end,
                 axis = 2,
                 intercept = 0,
                 coord = 'x',
                 vars=[],
                 returnAll = False):
                 
    ds = yt.load(path)
    dim = int(ds.domain_dimensions[0] > 1) + int(ds.domain_dimensions[1] > 1) + int(ds.domain_dimensions[2] > 1)

    if len(start) == 2: start.append(0.0)
    if len(end) == 2:   end.append(0.0)
    if dim == 2: 
        start[2] = 0
        end[2] = 0

    prof = ds.ray(start,end)
    slice2d = ds.slice(axis,intercept)

    new_df = prof.to_dataframe([("gas","x"),("gas","y"),("gas","z"),*vars])
    new_df = new_df.sort_values(by=coord)

    if returnAll:
        return new_df, slice2d, ds, dim
    else: 
        return new_df

def validate(path, 
             outdir,
             vars = [],
             start = None, #[0,0,0],
             end = None, #[1,1,1],
             axis = 2,
             intercept = 0,
             generate_ref_data = False,
             reference=None,
             tolerance=1E-8,
             abs_tolerance=None,
             coord = 'x',
             tight_layout = True,
             filename_suffix=None):
            
    info = readHeader(path)
            
    if not start: start = info["geom_lo"]
    if not end:   end = info["geom_hi"]
                 
    new_df, slice2d, ds, dim = readContours(path=path,start=start,end=end,axis=axis,
                                       intercept=intercept,coord=coord,vars=vars,returnAll=True)

    if not reference:
        reference = "reference/reference-{}d.csv".format(dim)

    if generate_ref_data:
        use_ray_coord = hasDuplicates(new_df[coord])
        ref_df, compare_coord = sortForComparison(new_df,start,end,coord,use_ray_coord)
        if compare_coord in ref_df and compare_coord.startswith("_"):
            ref_df = ref_df.drop(columns=[compare_coord])
        ref_df.to_csv(reference)
        print("Reference data generated, now exiting noisily so we don't accidentally pass the test")
        exit(-1)
    

    if not type(tolerance)==list:
        tolerance = [tolerance] * len(vars)
    elif len(tolerance) != len(vars):
        raise Exception("Wrong number of tolerance values")

    if not abs_tolerance:
        abs_tolerance = tolerance
    elif not type(abs_tolerance)==list:
        abs_tolerance = [abs_tolerance] * len(vars)
    elif len(abs_tolerance) != len(vars):
        raise Exception("Wrong number of rel_tolerance values")

    all_ok = True
    for i, (var, tol, abs_tol) in enumerate(zip(vars,tolerance, abs_tolerance)):
        if isinstance(reference,list):
            ref_df = pandas.read_csv(reference[i])
        else:
            ref_df = pandas.read_csv(reference)
        use_ray_coord = hasDuplicates(new_df[coord]) or hasDuplicates(ref_df[coord])
        new_compare_df, compare_coord = sortForComparison(new_df,start,end,coord,use_ray_coord)
        ref_df, compare_coord = sortForComparison(ref_df,start,end,coord,use_ray_coord)
        compare_label = "ray coordinate" if use_ray_coord else coord
        
        len_x = info["geom_hi"][0]-info["geom_lo"][0]
        len_y = info["geom_hi"][1]-info["geom_lo"][1]
        data2d = slice2d.to_frb(width=len_x,height=len_y,resolution=(1000,1000*len_y/len_x))[var]
        pylab.clf()
        im = pylab.imshow(data2d,origin='lower',cmap='jet',
                     extent=[info["geom_lo"][0],info["geom_hi"][0],info["geom_lo"][1],info["geom_hi"][1]])
        pylab.colorbar(im)
        pylab.plot([start[0],end[0]],[start[1],end[1]],linestyle='-',color='white')
        pylab.title(var)
        if tight_layout: pylab.tight_layout()
        pylab.savefig("{}/2d_{}.png".format(outdir,var))

        new_x,new_y,new_var = new_compare_df["x"], new_compare_df["y"], new_compare_df[var]
        ref_x,ref_y,ref_var = ref_df["x"], ref_df["y"], ref_df[var]

        new_coord = new_compare_df[compare_coord]
        ref_coord = ref_df[compare_coord]
            
        #print(type(ref_coord),axis=None)
        final_coord = sorted(numpy.concatenate((ref_coord, new_coord),axis=None))
        ref_final_var   = numpy.interp(final_coord, ref_coord, ref_var)
        new_final_var   = numpy.interp(final_coord, new_coord, new_var)
            
        pylab.clf()
        pylab.plot(final_coord,ref_final_var,color='C0',label='ref')
        pylab.plot(final_coord,new_final_var,color='C1',label='new',linestyle='--',marker='o',markerfacecolor='None')
        pylab.xlabel(compare_label)
        pylab.ylabel(var)
        pylab.legend()
        pylab.grid()
        pylab.title(var)
        pylab.savefig(outdir+"/{}.png".format(var))
        err = numpy.sqrt(integrate(final_coord, (    new_final_var  -     ref_final_var )**2))
        mag = numpy.sqrt(integrate(final_coord, (abs(new_final_var) + abs(ref_final_var))**2))


        if mag > 1e-8:
            relerr = err/mag
        else:
            relerr = err/(mag + 1e-12)

        print("{} abs error [tolerance={}]: {}".format(var,abs_tol,err))
        print("{} norm".format(var),mag)
        print("{} rel error [tolerance={}]: {}".format(var,tol,relerr))
        
        if relerr > tol:
            all_ok = False
            print("    {} relative error is too high".format(var))
        if err > abs_tol:
            all_ok = False
            print("    {} absolute error is too high".format(var))
        if math.isnan(relerr):
            all_ok = False
            print("    {} relative error is NAN".format(var))

    if all_ok:
        print("All errors are within tolerance")

    if not all_ok:
        print("ERROR: Regression test failed")
        raise(Exception("One or more errors out of tolerance or nan"))
