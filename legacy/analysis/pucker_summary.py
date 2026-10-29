import numpy as np
with open('pucker.dat') as f:
    hdr = f.readline().replace('#', '').split()
d = np.loadtxt('pucker.dat')
print(f'{len(d)} frames')
print(f'{"ring":6s} {"mean":>7s} {"min":>7s} {"max":>7s}  frames >40 deg off chair')
for k in range(18):
    name = hdr[1 + 3 * k]
    th = d[:, 3 + 3 * k]
    off = np.minimum(th, 180 - th)
    pole = '1C4' if np.median(th) > 90 else '4C1'
    print(f'{name:6s} {th.mean():7.1f} {th.min():7.1f} {th.max():7.1f}  '
          f'{(off > 40).sum():5d} ({100 * (off > 40).mean():4.1f}%)  {pole}')
