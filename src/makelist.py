# takes the arguments of max angular momentum for unrolling
# and file prefix, then generates a list of all the
# CPP files that will be generated: the unrolled kernels Q(LA, LB, lam)
# with LA <= LB <= max_unrol and LA + LB + lam <= 3 max_unrol (as in generate.cpp)

import sys
max_unrol = int(sys.argv[1])
prefix = str(sys.argv[2])

names = [prefix + "/generated/ecpint_gen.cpp"]
for j in range(max_unrol+1):
    for i in range(j+1):
        for k in range(3*max_unrol - i - j + 1):
            names.append(prefix + "/generated/Q" + str(i) + str(j) + str(k) + ".cpp")
file = open(prefix + "/qlist.txt", "w")
file.write("\n".join(names))
file.close()
