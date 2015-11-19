import glob
import sys
import pickle
import numpy as np

filepath = str(sys.argv[1])
glob_target = filepath[:-5] + '*'
files = glob.glob(glob_target)

for i, f in enumerate(files):
    print "Combining {}".format(f)
    f = pickle.load(file(f))
    if i == 0:
        combo = np.copy(f)
        continue
    else:
        combo[-1][0] = np.hstack((combo[-1][0], f[-1][0]))
        combo[-1][1] = np.hstack((combo[-1][1], f[-1][1]))
        combo[-1][2] = combo[-1][2] + f[-1][2]
        combo[-1][3] = combo[-1][3] + f[-1][3]

out_file_name = glob_target[:-2] + "-combined.pkl"
print "Writing {}".format(out_file_name)
f_out = file(out_file_name, 'w')
pickle.dump(combo, f_out)
f_out.close()
print "Done! Enjoy your burrito"
