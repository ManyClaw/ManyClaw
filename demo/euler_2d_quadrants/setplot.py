
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of clawpack.visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data
    

    # Figure for density - pcolor
    plotfigure = plotdata.new_plotfigure(name='Density', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Density'
    plotaxes.afteraxes = "pylab.axis('scaled')" 

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 0
    #plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.pcolor_cmin = 0.
    plotitem.pcolor_cmax = 2.0
    # plotitem.pcolor_cmax = 3.0
    plotitem.add_colorbar = True
    plotitem.amr_patchedges_show = [0]
    plotitem.amr_celledges_show = [0]


    # Figure for density - Schlieren
    plotfigure = plotdata.new_plotfigure(name='Schlieren', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Density'
    plotaxes.scaled = True      # so aspect ratio is 1

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_schlieren')
    plotitem.schlieren_cmin = 0.0
    plotitem.schlieren_cmax = 1.0
    plotitem.plot_var = 0
    plotitem.add_colorbar = False


    # Figure for grid cells
    plotfigure = plotdata.new_plotfigure(name='cells', figno=2)

    # This plot is not really relevant to ManyClaw right now
    plotfigure.show = False 

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0,1]
    plotaxes.ylimits = [0,1]
    plotaxes.title = 'Grid patches'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_patch')
    plotitem.amr_patch_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee']
    plotitem.amr_celledges_show = [1,0]
    plotitem.amr_patchedges_show = [1]

    # Profile plots
    plotfigure = plotdata.new_plotfigure(name='x_profile', figno=3)
    plotfigure.show = False

    def x_profile(cd, plot_var=0):
        slice_index = 0
        if plot_var == 0:
            q = cd.q[plot_var, :, slice_index]
        elif plot_var == 1:
            q = cd.q[1, :, slice_index] / cd.q[0, :, slice_index]
        elif plot_var == 2:
            gamma = 1.4
            q = (cd.q[3, : ,slice_index] - 0.5 * (cd.q[1,:,slice_index]**2) / cd.q[0,:,slice_index]) * (gamma - 1.0)

        return cd.x[:,0], q

    # Density
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0,1]
    plotaxes.ylimits = [0.9,3.1]
    plotaxes.title = 'X-Profile of Density'
    plotaxes.scaled = False
    plotaxes.axescmd = 'subplot(3,1,1)'
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.map_2d_to_1d = lambda cd: x_profile(cd, plot_var=0)
    plotitem.plotstyle = 'ob'

    # Velocity
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0,1]
    plotaxes.ylimits = [-0.1,0.5]
    plotaxes.title = 'X-Profile of Velocity'
    plotaxes.scaled = False
    plotaxes.axescmd = 'subplot(3,1,2)'
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.map_2d_to_1d = lambda cd: x_profile(cd, plot_var=1)
    plotitem.plotstyle = 'ob'

    # Pressure
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0,1]
    plotaxes.ylimits = [0.9,3.1]
    plotaxes.title = 'X-Profile of Pressure'
    plotaxes.scaled = False
    plotaxes.axescmd = 'subplot(3,1,3)'
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.map_2d_to_1d = lambda cd: x_profile(cd, plot_var=2)
    plotitem.plotstyle = 'ob'

    # Profile plots
    plotfigure = plotdata.new_plotfigure(name='y_profile', figno=4)
    plotfigure.show = False

    def y_profile(cd, plot_var=0):
        slice_index = 0
        if plot_var == 0:
            q = cd.q[plot_var,slice_index,:]
        elif plot_var == 1:
            q = cd.q[2, slice_index, :] / cd.q[0, slice_index, :]
        elif plot_var == 2:
            gamma = 1.4
            q = (cd.q[3,slice_index,:] - 0.5 * (cd.q[1,slice_index,:]**2) / cd.q[0,slice_index,:]) * (gamma - 1.0)

        return cd.y[0,:], q


    # Density
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0,1]
    plotaxes.ylimits = [0.9,3.1]
    plotaxes.title = 'Y-Profile of Density'
    plotaxes.scaled = False
    plotaxes.axescmd = 'subplot(3,1,1)'
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.map_2d_to_1d = lambda cd: y_profile(cd, plot_var=0)
    plotitem.plotstyle = 'ob'

    # Velocity
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0,1]
    plotaxes.ylimits = [-0.1,0.5]
    plotaxes.title = 'Y-Profile of Velocity'
    plotaxes.scaled = False
    plotaxes.axescmd = 'subplot(3,1,2)'
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.map_2d_to_1d = lambda cd: y_profile(cd, plot_var=1)
    plotitem.plotstyle = 'ob'

    # Pressure
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0,1]
    plotaxes.ylimits = [0.9,3.1]
    plotaxes.title = 'Y-Profile of Pressure'
    plotaxes.scaled = False
    plotaxes.axescmd = 'subplot(3,1,3)'
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.map_2d_to_1d = lambda cd: y_profile(cd, plot_var=2)
    plotitem.plotstyle = 'ob'


    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via clawpack.visclaw.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.html_movie = 'JSAnimation'      # new style, or "4.x" for old style
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

    
