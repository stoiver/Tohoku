
# coding: utf-8

# In[1]:


#get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt

# Allow inline jshtml animations
from matplotlib import rc
rc('animation', html='jshtml')


# In[2]:


#cd _output_Fujii/


# In[3]:


import anuga

figsize=(10,6)
dpi = 80
swwfile = 'Tohoku_Fujii.sww'

splotter = anuga.SWW_plotter(swwfile)


# In[ ]:


for frame, time in enumerate(splotter.time):
    print frame, time
    splotter.save_stage_frame(frame=frame, vmin=-10.0, vmax=10.0)
    plt.close('all')


# In[ ]:


#splotter.make_stage_animation()

