# takes the arguments of max angular momentum
# and file prefix, then generates a list of all the
# CPP files that will be generateds

import sys
max_am = int(sys.argv[1])
prefix = str(sys.argv[2])

file = open(prefix + "/qlist.txt", "w")
file.write(prefix + "/generated/ecpint_gen.cpp\n")
for j in range(max_am+1):
    for i in range(j+1):
        for k in range(max_am+1):
            if j == i == k == max_am:
                file.write(prefix + "/generated/Q" + str(i) + str(j) + str(k) + ".cpp")
            else:
                file.write(prefix + "/generated/Q" + str(i) + str(j) + str(k) + ".cpp\n")
file.close()
