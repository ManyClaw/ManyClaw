#!/usr/bin/env python

import matplotlib.pyplot as plt
import clawpack.pyclaw.solution as solution

num_frames = 2
solutions = []
fig = plt.figure()
for frame in xrange(num_frames):
    solutions.append(solution.Solution(frame))
    axes = fig.add_subplot(1,int(num_frames / 2) + 1,frame+1)
    x = solutions[frame].domain.grid.centers[0]
    y = solutions[frame].domain.grid.centers[1]
    plot = axes.pcolor(solutions[frame].q[0,:,:])
    axes.set_title("Solution at t = %s" % solutions[frame].t)
    axes.set_aspect('equal')

plt.show()