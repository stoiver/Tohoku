
import matplotlib.pyplot as plt

import anuga

figsize=(10, 6)
dpi = 80
swwfile = 'Tohoku_Fujii.sww'

splotter = anuga.SWW_plotter(swwfile)


for frame, time in enumerate(splotter.time):
    print frame, time
    splotter.save_stage_frame(frame=frame, vmin=-10.0, vmax=10.0)
    plt.close('all')
