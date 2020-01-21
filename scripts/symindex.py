#!/usr/bin/env python3
dim = 2


#sym = "full"
sym = "majorminor"
#sym = "major"

them = []
for i in range(0,dim):
    for j in range(0,dim):
        for k in range(0,dim):
            for l in range(0,dim):
                if sym=="full":
                    if not sorted([i,j,k,l]) in them:
                        them.append(sorted([i,j,k,l]))
                elif sym=="majorminor":
                    if not ((sorted([i,j])+sorted([k,l]) in them) or
                            (sorted([k,l])+sorted([i,j]) in them)):
                        them.append(sorted([i,j])+sorted([k,l]))
                elif sym=="major":
                    if not (([i,j]+[k,l] in them) or
                            ([k,l]+[i,j] in them)):
                        them.append([i,j]+[k,l])
 
print(len(them))
ctr=0
for t in them:
    print("// "+str(t))
    indices = []
    for i in range(0,dim):
        for j in range(0,dim):
            for k in range(0,dim):
                for l in range(0,dim):
                    if sym=="full":
                        srt = sorted([i,j,k,l])
                        if srt == t:
                            indices.append(i + dim*j + dim*dim*k + dim*dim*dim*l)
                    elif sym=="majorminor":
                        if ((sorted([i,j])+sorted([k,l])) == t or
                            (sorted([k,l])+sorted([i,j])) == t):
                            indices.append(i + dim*j + dim*dim*k + dim*dim*dim*l)
                    elif sym=="major":
                        if ([i,j]+[k,l]) == t or ([k,l]+[i,j]) == t:
                            indices.append(i + dim*j + dim*dim*k + dim*dim*dim*l)
    print (("if" if ctr==0 else "else if") + " (0 ",end="")
    for n in indices:
        print("|| uid=="+str(n)+" ",end="")
    print(") return data["+str(ctr)+"];")
    ctr += 1

print("Total of " + str(len(them)) + " unique constants")
