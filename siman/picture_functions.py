# -*- coding: utf-8 -*- 
from __future__ import division, unicode_literals, absolute_import 
import sys, os
import copy

import numpy as np

try:
    import scipy
    from scipy import interpolate
    # print (scipy.__version__)
    # print (dir(interpolate))
except:
    print('picture_functions.py: scipy is not avail')
# from scipy.interpolate import spline 
try:
    ''
    from scipy.interpolate import  CubicSpline
except:
    print('scipy.interpolate.CubicSpline is not avail')
try:
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
except:
    print('mpl_toolkits or matplotlib are not avail')

from matplotlib import gridspec 
try:
    from adjustText import adjust_text
    adjustText_installed = True
except:
    adjustText_installed = False



from siman import header
from siman.header import calc, printlog, printlog
from siman.inout import write_xyz
from siman.small_functions import makedir, is_list_like
from siman.geo import replic

# from siman.chg.chg_func import chg_at_point, cal_chg_diff
# from dos.functions import plot_dos
# from ase.utils.eos import EquationOfState

def plot_mep(atom_pos, mep_energies, image_name = None, filename = None, show = None, plot = 1, fitplot_args = None, style_dic = None):
    """
    Used for NEB method
    atom_pos (list) - xcart positions of diffusing atom along the path or just coordinates along one line (for polarons)
    mep_energies (list) - full energies of the system corresponding to atom_pos

    image_name - deprecated, use filename
    style_dic - dictionary with styles
        'p' - style of points
        'l' - style of labels
        'label' - label of points

    plot - if plot or not

    """


    from siman.analysis import determine_barrier

    if filename is None:
        filename = image_name

    #Create
    if not style_dic:
        style_dic = {'p':'ro', 'l':'b-', 'label':None}

    if 'p' not in style_dic:
        style_dic['p']='ro'

    if not fitplot_args:
        fitplot_args = {}

    # print
    if is_list_like(atom_pos[0]):
        atom_pos = np.array(atom_pos)
        data = atom_pos.T #
        tck, u= interpolate.splprep(data) #now we get all the knots and info about the interpolated spline
        path = interpolate.splev(np.linspace(0,1,500), tck) #increase the resolution by increasing the spacing, 500 in this example
        path = np.array(path)


        diffs = np.diff(path.T, axis = 0)
        path_length =  np.linalg.norm( diffs, axis = 1).sum()
        mep_pos =  np.array([p*path_length for p in u])
    else:
        mep_pos = atom_pos
        path_length  = atom_pos[-1]

    if 0: #plot the path in 3d
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.plot(data[0], data[1], data[2], label='originalpoints', lw =2, c='Dodgerblue')
        ax.plot(path[0], path[1], path[2], label='fit', lw =2, c='red')
        ax.legend()
        plt.show()





    # if '_mep' not in calc:
    calc['_mep'] = [atom_pos, mep_energies] # just save in temp list to use the results in neb_wrapper

    if hasattr(header, 'plot_mep_invert') and header.plot_mep_invert: # for vacancy
        mep_energies = list(reversed(mep_energies) )

    mine = min(mep_energies)
    eners = np.array(mep_energies)-mine

    
    
    xnew = np.linspace(0, path_length, 1000)

    # ynew = spline(mep_pos, eners, xnew )
    # spl = CubicSpline(mep_pos, eners, bc_type = 'natural' ) # second-derivative zero
    # spl = CubicSpline(mep_pos, eners,) #
    # spl = CubicSpline(mep_pos, eners, bc_type = 'periodic') 
    # spl = CubicSpline(mep_pos, eners, bc_type = 'clamped' ) #first derivative zero



    spl = scipy.interpolate.PchipInterpolator(mep_pos, eners)
    ynew = spl(xnew)

    diff_barrier = determine_barrier(mep_pos, eners)


    printlog('plot_mep(): Diffusion barrier =',round(diff_barrier, 2),' eV', imp = 'y')
    # sys.exit()
    # print()

    if 'fig_format' not in fitplot_args:
        fitplot_args['fig_format'] = 'pdf'

    if 'xlim' not in fitplot_args:
        fitplot_args['xlim'] = (-0.05, None  )

    if 'xlabel' not in fitplot_args:
        fitplot_args['xlabel'] = 'Reaction coordinate (${\\rm \AA}$)'


    if 'ylabel' not in fitplot_args:
        fitplot_args['ylabel'] = 'Energy (eV)'

    path2saved = None
    if plot:
        # print(style_dic.get('color'))
        path2saved = fit_and_plot(orig = {'x':mep_pos, 'y':eners, 'fmt':style_dic['p'], 'label':style_dic['label'], 'color':style_dic.get('color')}, 
            spline = {'x':xnew, 'y':ynew, 'fmt':style_dic['l'], 'label':None, 'color':style_dic.get('color')}, 
        image_name =  image_name, filename = filename, show = show, 
        **fitplot_args)

        # print(image_name, filename)
        if 0:
            with open(filename+'.txt', 'w') as f:
                f.write('DFT points:\n')
                for m, e in zip(mep_pos, eners):
                    f.write('{:10.5f}, {:10.5f} \n'.format(m, e))
                f.write('Spline:\n')
                for m, e in zip(xnew, ynew):
                    f.write('{:10.5f}, {:10.5f} \n'.format(m, e))





    return path2saved, diff_barrier


def process_fig_filename(image_name, fig_format):

    makedir(image_name)

    if fig_format in image_name:
        path2saved = str(image_name)

    elif str(image_name).split('.')[-1] in ['eps', 'png', 'pdf']:
        path2saved = str(image_name)
        fig_format = str(image_name).split('.')[-1]

    else:
        path2saved = str(image_name)+'.'+fig_format
    
    dirname = os.path.dirname(image_name)
    if not dirname:
        dirname+='.'

    path2saved_png = dirname+'/png/'+os.path.basename(image_name)+'.png'
    makedir(path2saved_png)

    return path2saved, path2saved_png

def fit_and_plot(ax = None, power = None, xlabel = None, ylabel = None,
    image_name = None, filename = None,
    show = None, pad = None,
    xlim = None, ylim = None, title = None, figsize = None,
    xlog = False,ylog = False, scatter = False, 
    legend = False, ncol = 1, 
    fontsize = None, legend_fontsize=None, markersize = None,  
    linewidth = None, hor = False, ver = False, fig_format = 'pdf', dpi = 300,
    ver_lines = None, hor_lines = None, xy_line = None, x_nbins = None,
    alpha = 0.8, fill = False,
    first = True, last = True, 
    convex = None, dashes = None,
    corner_letter = None, corner_letter_pos = None, hide_ylabels = None, hide_xlabels= None, annotate = None,
    params = None,
    **data):
    """
    Produce complex plots using arbirtary number of axes and arbirtary number of curves on each axis. 
    See tutorial with examples (to be created). 
    
    INPUT: 

        - ax (axes) - matplotlib axes object - to create multiple axes plots. If None than single axes is used
        - data (tuple or dict) - each entry should be 
            - (X, Y, fmt) 
            - (X, Y, fmt, label) 
            - (X, Y, R, fmt) - for scatter = 1, R - size of spots; 
            - {'x':,'y':, 'fmt':, ..., *any argument valid for pyplot.plot()* }  not implemented for powers and scatter yet
                - 'label' (str) - 
                - 'xticks' (list) - 
                - 'annotate_fontsize' (str)
                - 'annotate_arrowprops' (dict)
        - first, last (int) - allows to call this function multiple times to put several axes on one impage. Use first = 1, last = 0 for the first axes, 0, 0 for intermidiate, and 0, 1 for last axes.
        - power (int) - the power of polynom, turns on fitting
        - scatter (bool) - plot scatter points - the data format is slightly different - see *data*
        - convex (bool) - plot convex hull around points like in ATAT
        - fill (bool) - fill under the curves
        - filename (str) - name of file with figure, image_name - deprecated
        - fig_format (str) - format of saved file.
        - dpi    - resolution of saved file
        - ver_lines - list of dic args for  vertical lines {'x':, 'c', 'lw':, 'ls':}
        - hor_lines
        - ver - vertical line at 0
        - hor - horizontal line at 0
        - hide_ylabels - just hide numbers
        - ncol - number of legend columns
        - corner_letter - letter in the corner of the plot
        - corner_letter_pos (list*2 float) - list with [x,y] corner position, default left upper corner is set
        - pad - additional padding, if dict than the same keys as in plt.subplots_adjust() are used
        - annotate - annotate each point, 'annotates' list should be in data dic!
        - linewidth - was 3 !
        - markersize - was 10
        - x_nbins - number of ticks
        - params (dict) - dictionary with parameters, should be used instead of new arguments
            - 'xlim_power' - xlim for power
            - 'y0' - move plot to have y = 0
            - 'xnbins' - number of bins x

    TODO:

        - remove some arguments that can be provided in data dict
        - move all rare arguments to params


    RETURN: 

        filename of the saved image
    
    AUTHOR:

        Aksyonov D.A.


    """


    if image_name == None:
        image_name  = filename

    # print(ver_lines)
    # sys.exit()
    # fontsize = 1
    # print(fontsize)
    # sys.exit()
    if fontsize:
        # header.mpl.rcParams.update({'font.size': fontsize+4})
        # fontsize = 2
        SMALL_SIZE = fontsize
        MEDIUM_SIZE = fontsize
        BIGGER_SIZE = fontsize

        header.mpl.rc('font', size=SMALL_SIZE)          # controls default text sizes
        header.mpl.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
        header.mpl.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
        header.mpl.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        header.mpl.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        header.mpl.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
        header.mpl.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


        # font = {'family' : 'normal',
        #         'weight' : 'bold',
        #         'size'   : fontsize}

        # header.mpl.rc('font', **font)


        if legend_fontsize is None:
            legend_fontsize = fontsize
    



    if legend_fontsize:
        ''
        header.mpl.rc('legend', fontsize= legend_fontsize) 

    if corner_letter_pos is None:
        corner_letter_pos = [0.05,0.8]
    # print('fontsize', fontsize, legend_fontsize)


    if hasattr(header, 'first'):
        first = header.first

    if hasattr(header, 'last'):
        last  = header.last
    # print('fit_and_plot, first and last', first, last)



    # print('ax is', ax)
    
    if ax is None:
        if first:
            # fig, ax = plt.subplots(1,1,figsize=figsize)
            plt.figure(figsize=figsize)

        ax = plt.gca() # get current axes )))
        # ax  = fig.axes
    # print('ax is', ax)

    if title: 
        ax.title(title)
    
    # print(dir(plt))
    # print(ax)

    # plt.ylabel(ylabel, axes = ax)
    # print(ylabel)
    if ylabel is not None:

        ax.set_ylabel(ylabel)
    

    if xlabel is not None:
        ''
        # plt.xlabel(xlabel, axes = ax)
        ax.set_xlabel(xlabel)

    if params is None:
        params = {}


    if corner_letter:
        # print(corner_letter)
        ''
        sz = header.mpl.rcParams['font.size']
        ax.text(corner_letter_pos[0],corner_letter_pos[1], corner_letter, size = sz*1.5, transform=ax.transAxes) # transform = None - by default in data coordinates!

        # text(x, y, s, bbox=dict(facecolor='red', alpha=0.5))


    if convex:
        from scipy.spatial import ConvexHull


    keys = []
    shift = 0
    for key in sorted(data):
        keys.append(key)
        if scatter:
            
            ax.scatter(data[key][0], data[key][1],  s = data[key][2], c = data[key][-1], alpha = alpha, label = key)
        
        else:

            con = data[key]
            # print('keys', keys)
            # print('con type', type(con))
            # sys.exit()
            if type(con) == list or type(con) == tuple:
                try:
                    label = con[3]
                except:
                    label = key

                try:
                    fmt = con[2]
                except:
                    fmt = ''


                xyf = [con[0], con[1], fmt]
                con = {'label':label} #fmt -color style
                # print('con1', con)
                # del con['fmt']

            elif type(con) == dict:
                if 'fmt' not in con:
                    con['fmt'] = ''
                # print(con)

                if 'x' not in con:
                    l = len(con['y'])
                    con['x'] = list(range(l))

                if 'xticks' in con:
                    # print(con['xticks'])
                    ax.set_xticklabels(con['xticks'])
                    ax.set_xticks(con['x'])
                    del con['xticks']

                xyf = [con['x'], con['y'], con['fmt']]
                del con['fmt']

                # if 'lw' in 
            if linewidth:
                con['lw'] = linewidth

            # print('con2', con)
            # sys.exit()
                
            if markersize:
                con['ms'] = markersize

            # print('key is ', key)
            # print('x ', xyf[0])






            con_other_args = copy.deepcopy(con) # given to plot
            # print('con_copy', con_other_args)
            # sys.exit()
            for k in ['x', 'x2', 'x2label', 'y', 'fmt', 'annotates', 'x2_func', 'x2_func_inv', 
            'annotate_fontsize', 'annotate_arrowprops']:
                if k in con_other_args:
                    del con_other_args[k]
            if 'color' in con_other_args:
                if con_other_args['color'] is None:
                    del con_other_args['color']
            # print('con', con_other_args)
            # sys.exit()
            
            if params.get('y0'):
                if key == keys[0]:
                    shift = min(xyf[1])

            if shift:
                xyf[1] = list(np.array(xyf[1])-shift)

            # print(con_other_args)
            if 'ls' in con_other_args:
                ''
                del xyf[-1]
            # print(xyf)
            ax.plot(*xyf, alpha = alpha, **con_other_args)


            #second x axis
            if con.get('x2_func' ):
                # ax2 = ax.twiny()
                ax2 = ax.secondary_xaxis("top", functions=(con['x2_func'],con['x2_func_inv']))
                # ax2.plot(con['x2'], con['x2'])
                ax2.set_xlabel(con['x2label'])
                # ax2.cla()

            if power:
                coeffs1 = np.polyfit(xyf[0], xyf[1], power)        
                
                fit_func1 = np.poly1d(coeffs1)

                if params.get('xlim_power'):
                    x_range = np.linspace(params['xlim_power'][0], params['xlim_power'][1])

                else:
                    x_range = np.linspace(min(xyf[0]), max(xyf[0]))
                
                fit_y1 = fit_func1(x_range); 
         

                # ax.plot(x_range, fit_y1, xyf[2][0]+'--', )
                ax.plot(x_range, fit_y1, '--', )

                # x_min  = fit_func2.deriv().r[power-2] #derivative of function and the second cooffecient is minimum value of x.
                # y_min  = fit_func2(x_min)
                from scipy import stats
                slope, intercept, r_value, p_value, std_err = stats.linregress(xyf[0], xyf[1])
                print ('R^2 = {:5.2f} for {:s}'.format(r_value**2, key))






            adjustText_installed =0 
            if annotate:
                fs = con.get('annotate_fontsize') or fontsize
                arrowprops = con.get('annotate_arrowprops')
                if 'annotate_arrowprops' not in con:
                    arrowprops = dict(arrowstyle='->', connectionstyle='arc3,rad=0.5', color='black')


                if adjustText_installed:
                    # sys.exit()
                    ts = []
                    for t, x, y in zip(con['annotates'], con['x'], con['y']):
                        ''
                        # print(x,y,t)
                        ts.append(ax.text(x, y, t, size = fs, alpha = 1, color = con['fmt'][0]))
                    # ax.text(6, 6, 'h', size = 10, alpha = 0.5, color = 'k')
                    if 1:
                        adjust_text(ts, ax = ax,
                            # force_points=10, force_text=10, force_objects = 0.5, 
                            expand_text=(2, 2), 
                            expand_points=(2, 2), 
                            # lim = 150,
                            expand_align=(2, 2), 
                            # expand_objects=(0, 0),
                            text_from_text=1, text_from_points=0.5,
                            # arrowprops=dict(arrowstyle='->', color='black')
                            )


                else:
                    for name, x, y in zip(con['annotates'], con['x'], con['y']):
                        ax.annotate(name, xy=(x, y),
                        xytext = (-fs//2, fs//2), fontsize = fs, color = 'r',
                        textcoords = 'offset points', ha='center', va='bottom',
                        # bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.3),
                        arrowprops = arrowprops
                        )            





            # print(key)
            # print(con)
            if fill:
                ''
                # print('fill', xyf[0], xyf[1])
                # ax.fill(xyf[0], xyf[1], facecolor = con['c'], alpha = 0.6)
                ax.fill_between(xyf[0], xyf[1], y2=0, facecolor = con['c'], alpha = 0.6)


            if convex:
                points = np.asarray(list(zip(xyf[0], xyf[1])))
                hull = ConvexHull(points)
                for simplex in hull.simplices:
                    if max(points[simplex, 1]) > 0:
                        continue
                    ax.plot(points[simplex, 0], points[simplex, 1], 'k-')



    if not linewidth:
        linewidth = 1
    if hor: 
        ax.axhline(color = 'k', lw = linewidth, alpha = 0.6, ls = '-') #horizontal line

    if ver:

        ax.axvline(color='k', lw = linewidth, alpha = 0.6, ls = '-') # vertical line at 0 always 

    if ver_lines:
        for line in ver_lines:
            ax.axvline(**line)

    if hor_lines:
        for line in hor_lines:
            ax.axhline(**line)

    if xlim: 
        ax.set_xlim(xlim)

    if ylim:
        ax.set_ylim(ymin=ylim[0])
        if ylim[1]: 
            ax.set_ylim(ymax=ylim[1])



    if xlog: 
        ax.set_xscale('log')

    if ylog: 
        if "sym" in str(ylog):
            ax.set_yscale('symlog', linthreshx=0.1)
        else:
            ax.set_yscale('log')

    if xy_line:
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        # print(xlim, ylim)
        # sys.exit()
        x = np.linspace(*xlim)
        y = np.linspace(*ylim)

        # print(x)
        if abs(x[0]-x[-1])>abs(y[0]-y[-1]):
            da = x
        else:
            da = y
        ax.plot(da,da, color = 'k', alpha = 0.6, ls = '-', lw = linewidth)

    if x_nbins:
        ax.locator_params(nbins=x_nbins, axis='x')




    if hide_ylabels:
        ax.yaxis.set_major_formatter(plt.NullFormatter())
        # ax.yaxis.set_ticklabels([])
    if hide_xlabels:
        ax.xaxis.set_major_formatter(plt.NullFormatter())


    if legend: 
        scatterpoints = 1 # for legend

        ax.legend(loc = legend, scatterpoints = scatterpoints, ncol = ncol)
        # plt.legend()

    # plt.tight_layout(pad = 2, h_pad = 0.5)


    if params.get('xnbins'):
        # print(params.get('xnbins'))
        ax.locator_params(tight=True, axis='x', nbins=params['xnbins'])
        # plt.locator_params(axis='x', numticks=params['xnbins'])

    if params.get('xticks_step'):

        start, end = ax.get_xlim()
        ax.xaxis.set_ticks(np.arange(start+params.get('step_shift'), end+params.get('step_shift'), params.get('xticks_step')))


    plt.tight_layout()

    if pad:
        # pad =0.13
        if type(pad) == dict:
            plt.subplots_adjust(**pad)
        else:
            plt.subplots_adjust(left=0.13, bottom=None, right=None, top=None,
                        wspace=0.07, hspace=None)


    path2saved = ''
    
    if last:

        if image_name:
            # plt.subplots_adjust(hspace=0.1)

            path2saved, path2saved_png = process_fig_filename(image_name, fig_format)

            plt.savefig(path2saved, dpi = dpi, format=fig_format, bbox_inches = "tight")
            plt.savefig(path2saved_png, dpi = 300)
            
            printlog("Image saved to ", path2saved, imp = 'y')


        elif show is None:
            show = True
        # printlog(show)
        if show:
            plt.show()
        plt.clf()
        plt.close('all')
    else:
        printlog('Attention! last = False, no figure is saved')

    return path2saved

fitplot = fit_and_plot

def plot_bar(xlabel = "xlabel", ylabel = "ylabel",
    xlim = None, ylim = None, legend = 0,
    image_name = None, title = None, bottom = 0.18, hspace = 0.15, barwidth = 0.2,
    data1 = [],data2 = [],data3 = [],data4 = [], hor_lines = None,
    **data):

    width = barwidth      # the width of the bars

    if data: 
        N = len(data.values()[0][0])
        key = data.keys()[0]
        xlabels = data[key][0]
        # print N
        ind = np.arange(N)  # the x locations for the groups
        shift = 0
        fig, ax = plt.subplots()
        for key in sorted(data):
            # print 'color', data[key][2]
            ax.bar(ind+shift, data[key][1], width, color = data[key][2], label = data[key][-1])# yerr=menStd)
            # print ind
            shift+=width

    elif data1 and data4:
        fig = plt.figure(figsize=(10,5))   #5:7 ratio for A4, 
        gs = gridspec.GridSpec(2, 2,
                               width_ratios =[5,1],
                               height_ratios=[1,1]
                               )
        gs.update(top=0.98, bottom=bottom, left=0.1, right=0.98, wspace=0.15, hspace=hspace)

        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])
        ax3 = plt.subplot(gs[2])
        ax4 = plt.subplot(gs[3])
        # fig, ax = plt.subplots()
        # fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col')#, sharey='row') equal

        for ax, data in (ax1, data1), (ax2,data2), (ax3,data3), (ax4, data4):
            N = len(data[0][0])
            xlabels = data[0][0]
            ind = np.arange(N)  # the x locations for the groups
            shift = 0

            for d in data:
                ax.bar(ind+shift, d[1], width, color = d[2], label = d[-1])# yerr=menStd)
                # print ind
                shift+=width

            ax.axhline(y=0, color='black')
            # ax.set_xticklabels(xlabels , rotation=70 )

            ax.set_xticks(ind+width)
    
        ax3.set_xticklabels(data3[0][0] , rotation=80 )
        ax4.set_xticklabels(data4[0][0] , rotation=80 )

        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.setp(ax2.get_xticklabels(), visible=False)
            
        ax3.set_ylabel(ylabel)
        ax3.yaxis.set_label_coords(-0.1, 1.1)
        # plt.ylabel(ylabel)

        ax1.legend(loc=2, )
        ax3.legend(loc=2, )
        # ax1.axis('tight')
        # ax2.axis('tight')
        # ax3.axis('tight')
        # ax4.axis('tight')
        ax1.margins(0.0, 0.2)
        ax2.margins(0.0, 0.2)
        ax3.margins(0.0, 0.2)
        ax4.margins(0.0, 0.2)





    elif data1 and data2 and not data4:
        fig = plt.figure(figsize=(10,5))   #5:7 ratio for A4, 
        gs = gridspec.GridSpec(1, 2,
                               width_ratios =[5,1],
                               height_ratios=[1,0]
                               )
        gs.update(top=0.95, bottom=bottom, left=0.1, right=0.98, wspace=0.15, hspace=hspace)

        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])

        for ax, data in (ax1, data1), (ax2,data2):
            N = len(data[0][0])
            xlabels = data[0][0]
            ind = np.arange(N)  # the x locations for the groups
            # print ind+width
            # print data[0][0]
            shift = 0.2

            for d in data:
                ax.bar(ind+shift, d[1], width, color = d[2], label = d[-1])# yerr=menStd)
                # print ind
                shift+=width

            ax.axhline(y=0, color='black')
            # ax.set_xticklabels(xlabels , rotation=70 )

            ax.set_xticks(ind+width+len(data)*width/2)
            
        names1 = [ n1  for n1, n2 in  zip( data1[0][0], data1[1][0] ) ] # 
        names2 = [ n1  for n1, n2 in  zip( data2[0][0], data2[1][0] ) ]

        ax1.set_xticklabels( names1, rotation = 80 ) # Names of configurations on x axis
        ax2.set_xticklabels( names2, rotation = 80 ) 

        ax1.set_ylabel(ylabel)

        ax1.legend(loc=2, )
        ax1.axis('tight')
        ax2.axis('tight')

    elif data1 and not data2:
        fig = plt.figure(figsize=(10,5))   #5:7 ratio for A4, 
        # gs = gridspec.GridSpec(1,2, width_ratios =[9,1])#, height_ratios=[1,1])                               
        # gs.update(top=0.95, bottom=bottom, left=0.1, right=0.98, wspace=0.15, hspace=hspace)

        # ax1 = plt.subplot(gs[0])
        ax1 = plt.subplot()
        # ax2 = plt.subplot(gs[1])

        for ax, data in (ax1, data1),:
            N = len(data[0][0])
            xlabels = data[0][0]
            ind = np.arange(N)  # the x locations for the groups
            # print ind+width
            # print data[0][0]
            shift = 0.2

            for d in data:
                ax.bar(ind+shift, d[1], width, color = d[2], label = d[-1])# yerr=menStd)
                # print ind
                shift+=width

            ax.axhline(y=0, color='black')
            # ax.set_xticklabels(xlabels , rotation=70 )

            ax.set_xticks(ind+width+len(data)*width/2)
        
        if len(data1) == 2:
            # names1 = [ n1 + '; ' + n2 for n1, n2 in  zip( data1[0][0], data1[1][0] ) ] # 
            names1 = [ n1  for n1 in  data1[0][0] ] # 

        else:
            names1 = [ n1  for n1 in  data1[0][0] ] # 

        ax1.set_xticklabels( names1, rotation = 80 ) # Names of configurations on x axis

        ax1.set_ylabel(ylabel)
        if legend:
            # ax1.legend(loc=2, )
            ax1.legend( )
        # ax1.axis('tight')
        # ax2.axis('tight')
        if ylim:
            ax1.set_ylim(ylim)
    # ax.set_yscale('log')
    # plt.yscale('symlog', linthreshx=0.1)

    # ax.set_title('Scores by group and gender')
    if hor_lines:
        for line in hor_lines:
            ax1.axhline(**line)

    def autolabel(rects):
        # attach some text labels
        for rect in rects:
            height = rect.get_height()
            ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
                    ha='center', va='bottom')

    # autolabel(rects1)
    # autolabel(rects2)


    # plt.axis('tight')
    # plt.margins(0.05, 0)

    # plt.tight_layout()
    # elif data1: gs.tight_layout(fig)
    plt.tight_layout()
    if image_name:
        printlog( "Saving image ...", str(image_name), imp = 'y')
        plt.savefig(str(image_name)+'.png', dpi = 200, format='png')
    else:
        plt.show()

    return





def plot_and_annotate(power = 2, xlabel = "xlabel", ylabel = "ylabel", filename = None,
    xlim = None, ylim = None, title = None, fit = None,
    legend = None, 
    **data):
    """Should be used in two below sections!
    Creates one plot with two dependecies and fit them;
    return minimum fitted value of x2 and corresponding valume of y2; 
    if name == "" image will not be plotted
    power - the power of polynom

    data - each entry should be (X, Y, 'r-')
    """

    # print data
    # coeffs1 = np.polyfit(x1, y1, power)        
    # coeffs2 = np.polyfit(x2, y2, power)
    
    # fit_func1 = np.poly1d(coeffs1)
    # fit_func2 = np.poly1d(coeffs2)
    
    #x_min  = fit_func2.deriv().r[power-2] #derivative of function and the second cooffecient is minimum value of x.
    #y_min  = fit_func2(x_min)
    
    image_name = filename
    if 1:
        # x_range = np.linspace(min(x2), max(x2))
        # fit_y1 = fit_func1(x_range); 
        # fit_y2 = fit_func2(x_range); 
        
        plt.figure()
        if title: plt.title(title)
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        




        for key in data:
            plt.plot(data[key][0], data[key][1], data[key][-1], markersize = 15, label = key)
            
            for x, y, name in zip(data[key][0], data[key][1], data[key][2]): 
                xytext = (-20,20)
                if 'T1m' in name: xytext = (20,20)
                plt.annotate(name, xy=(x, y), xytext=xytext, 
                textcoords='offset points', ha='center', va='bottom',
                bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.3),
                arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5', 
                                color='red'))










        if fit:
            for key in data:

                f1 = interp1d(data[key][0], data[key][1], kind='cubic')
                x = np.linspace(data[key][0][0], data[key][0][-1], 100) 
                plt.plot(x, f1(x), '-', label = key+'fit')


        plt.axvline(color='k')
        if xlim: 
            plt.xlim(xlim)
            # axes = plt.gca()
            # axes.set_xlim([xmin,xmax])
            # axes.set_ylim([ymin,ymax])
        if ylim:
            plt.ylim(ymin=ylim[0])
            if ylim[1]: plt.ylim(ymax=ylim[1])


        # plt.plot(x2, y2, 'bo', label = 'r'   )
        # plt.plot(x_range, fit_y1, 'r-', label = 'init_fit')
        # plt.plot(x_range, fit_y2, 'b-', label = 'r_fit'   )



        plt.tight_layout()

        if legend: plt.legend(loc = 2)

        if image_name:
            # print "Saving image ..."
            if not os.path.exists('images/'):
                os.makedirs('images/')
            plt.savefig(str(image_name)+'.png', dpi = 300, format='png')
            plt.close()
        else:
            plt.show()


    return 


def plot_bar_simple(xlabel = "xlabel", ylabel = "ylabel",
    xlim = None, ylim = None,
    image_name = None, title = None, 
    data = []):

    width = 0.6       # the width of the bars



    plt.figure(figsize=(10,5))   #5:7 ratio for A4, 
    fig, ax = plt.subplots()


    N = len(data[0])
    xlabels = data[0]
    # print xlabels
    ind = np.arange(N)  # the x locations for the groups
    shift = 0
    # for d in data:
    d = data
    # print d[2]

    rects = ax.bar(ind+shift, d[1], width, color = d[2], label = d[-1],align="center" )# yerr=menStd)
    rects[0].set_color('g')
    rects[-1].set_color('g')
    # rects[1].set_color('b')

    # print ind
    shift+=width

    ax.axhline(y=0, color='black')

    ax.set_xticks(ind)#+width)
    ax.set_xticklabels(xlabels , rotation=50 )
        
    ax.set_ylabel(ylabel)
    
    handles, labels = ax.get_legend_handles_labels()
    import matplotlib.patches as mpatches
    red_patch = mpatches.Patch(color='red', label='Substitutional')
    ax.legend(handles+[red_patch], labels+['Substitutional'], loc = 4)
    # ax.legend(handles=[red_patch], loc = 8)
    

    # ax.legend(loc=2, )

    # ax.set_yscale('log')
    # plt.yscale('symlog', linthreshx=0.1)

    # ax.set_title('Scores by group and gender')
    if xlim: 
        plt.xlim(xlim)
        # axes = plt.gca()
        # axes.set_xlim([xmin,xmax])
        # axes.set_ylim([ymin,ymax])
    if ylim:
        plt.ylim(ymin=ylim[0])
        if ylim[1]: plt.ylim(ymax=ylim[1])

    def autolabel(rects):
        # attach some text labels
        for rect in rects:
            height = rect.get_height()
            ax.text(rect.get_x()+rect.get_width()/2., -1.05*height, '%.0f'%float(height),
                    ha='center', va='top')
    autolabel(rects)
    # autolabel(rects2)
    plt.tight_layout()

    if image_name:
        # print "Saving image ..."
        if not os.path.exists('images/'):
            os.makedirs('images/')
        plt.savefig('images/'+str(image_name)+'.png', dpi = 200, format='png')
        plt.close()
    else:
        plt.show()
    return








def plot_conv(list_of_calculations = None, calc = None, 
    type_of_plot = None, conv_ext = [], labelnames = None, cl = None,
    plot = 1, filename = None):
    """
    Allows to fit and plot different properties;
    Input:
    'type_of_plot' - ("fit_gb_volume"-fits gb energies and volume and plot dependencies without relaxation and after it,
     'dimer'

    cl - calculation to use - new interface, please rewrite the old one

    """

            
    def fit_and_plot(x1, y1, x2, y2, power, name = "", xlabel = "", ylabel = "", image_name = "test", lines = None):
        """Should be used in two below sections!
        Creates one plot with two dependecies and fit them;
        return minimum fitted value of x2 and corresponding valume of y2; 
        if name == "" image will not be plotted
        power - the power of polynom

        lines - add lines at x = 0 and y = 0

        """
        coeffs1 = np.polyfit(x1, y1, power)        
        coeffs2 = np.polyfit(x2, y2, power)
        
        fit_func1 = np.poly1d(coeffs1)
        fit_func2 = np.poly1d(coeffs2)
        
        #x_min  = fit_func2.deriv().r[power-2] #derivative of function and the second cooffecient is minimum value of x.
        #y_min  = fit_func2(x_min)
        
        if name:

            x_range = np.linspace(min(x2), max(x2))
            fit_y1 = fit_func1(x_range); 
            fit_y2 = fit_func2(x_range); 
            
            plt.figure(figsize=(8,6.1))
            # plt.title(name)
            plt.ylabel(ylabel)
            plt.xlabel(xlabel)
            plt.xlim(min(x2)-0.1*abs(min(x2) ), max(x2)+0.1*abs(min(x2)))

            plt.plot(x1, y1, 'ro', label = 'initial')
            plt.plot(x2, y2, 'bo', label = 'relaxed'   )
            plt.plot(x_range, fit_y1, 'r-',) #label = 'init_fit')
            plt.plot(x_range, fit_y2, 'b-',) #label = 'r_fit'   )
            plt.legend(loc =9)
            
            if lines == 'xy':
                plt.axvline(color='k')
                plt.axhline(color='k')



            plt.tight_layout()
            #plt.savefig('images/'+image_name)
            file = header.path_to_images+'/'+str(image_name)+'.png'
            makedir(file)
            printlog( 'Saving file ...',file, imp = 'y' )
            plt.savefig(file,format='png', dpi = 300)
        return fit_func2  



    if list_of_calculations:
        conv = list_of_calculations
        n = conv[0]

        name = []; 
       
        name.append( n[0] )
        
        image_name = n[0]+'_'+n[1]+'_'+str(n[2])

    energies = []; init_energies = []
    volumes = []
    gb_volumes = []
    pressures = []
    pressures_init = []
    sigma_xx = []
    sigma_yy = []
    sigma_zz = []
    e_gbs = [] 
    e_gbs_init = []

    if type_of_plot == "e_imp":
        e_imps = []
        v_imps = []
        lengths = []
        for id in conv:        
            e_imps.append(calc[id].e_imp*1000)
            v_imps.append(calc[id].v_imp)
            l = calc[id].vlength
            lengths.append( "%s\n%.1f\n%.1f\n%.1f"%(id[0],l[0], l[1], l[2]) )
        #l = lengths[0]
        #print str(l[0])+'\n'+str(l[1])+'\n'+str(l[2])


        xlabel = "Sizes, $\AA$"
        ylabel = "Impurity energy, meV"
        ylabel2 = "Impurity volume, $\AA^3$"
        plt.figure()
        plt.title(str(name)+' other cells')
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        x = range( len(e_imps) )
        plt.xticks(x, lengths)
        plt.plot(x, e_imps, 'ro-', label = 'energy')
        plt.legend()
        plt.twinx()
        plt.ylabel(ylabel2)
        plt.plot(x, v_imps, 'bo-', label = 'volume')

        plt.subplots_adjust(left=None, bottom=0.2, right=None, top=None,
                wspace=None, hspace=None)
        #plt.ticker.formatter.set_scientific(True)
        plt.legend(loc =9)
        plt.savefig('images/e_imp_'+str(image_name)+'.png',format='png')#+str(image_name))#+'e_imp')


    if type_of_plot == "e_2imp":

        def dist_between_imp(cl):
            """Only for two impurities"""

            return np.linalg.norm(cl.end.xcart[-1] - cl.end.xcart[-2]) #assuming that impurities are at the end of xcart list.

        e_imps = [] # binding energy
        dist = [] #dist between impurities
        
        e_imps_ex = []
        dist_ex = []
        name_ex = []

        for id in conv:        
            cl = calc[id]
            e_imps.append(cl.e_imp*1000)
            #dist.append( "%s\n%.1f"%(id[0],dist_between_imp(cl) ) )
            dist.append( dist_between_imp(cl)  )







        xlabel = "Distance between atoms, $\AA$"
        ylabel = "Interaction energy, meV"
        plt.figure()
        
        # plt.title(str(name)+' v1-15')

        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        #x = range( len(e_imps) )
        #plt.xticks(x, dist)
        # plt.yscale('log')
        # plt.yscale('semilog')
        if labelnames:
            label = labelnames
        else:
            label = []
            label[0] = str(name)
            label[1] = name_ex[0]
            label[2] = name_ex[1]


        plt.plot(dist, e_imps, 'ro-', label = label[0], linewidth = 2 )
        

        if conv_ext: #manually add 
            for conv in conv_ext:
                e_imps_ex.append([])
                dist_ex.append([])
                for id in conv:        
                    cl = calc[id]
                    e_imps_ex[-1].append(cl.e_imp*1000)
                    #dist.append( "%s\n%.1f"%(id[0],dist_between_imp(cl) ) )
                    dist_ex[-1].append( dist_between_imp(cl)  )
                name_ex.append(id[0])
            plt.plot(dist_ex[0], e_imps_ex[0], 'go-', label = label[1], linewidth = 2)
            plt.plot(dist_ex[1], e_imps_ex[1], 'bo-', label = label[2], linewidth = 2)






        plt.axhline(color = 'k') #horizontal line

        plt.tight_layout()

        # plt.subplots_adjust(left=None, bottom=0.2, right=None, top=None,
        #         wspace=None, hspace=None)
        # #plt.ticker.formatter.set_scientific(True)
        plt.legend(loc =9)
        plt.savefig(path_to_images+'e_2imp_'+str(image_name)+'.png',format='png', dpi = 300)#+str(image_name))#+'e_imp')


    if type_of_plot == "fit_gb_volume_pressure":

        for id in conv:
            #energies.append(calc[id].energy_sigma0)
            #init_energies.append( calc[id].list_e_sigma0[0] ) 
            gb_volumes.append(calc[id].v_gb)
            #volumes.append(calc[id].end.vol)
            pressures.append(calc[id].extpress/1000. )
            pressures_init.append(calc[id].extpress_init/1000. )
            sigma_xx.append( calc[id].stress[0]  )
            sigma_yy.append( calc[id].stress[1]  )
            sigma_zz.append( calc[id].stress[2]  )
            #pressures_init = pressures
            e_gbs.append(calc[id].e_gb)
            e_gbs_init.append(calc[id].e_gb_init )           
            # print calc[id].bulk_extpress

        power = 3

        fit_ve = fit_and_plot(gb_volumes, e_gbs_init,  gb_volumes, e_gbs, power, 
            name, "Grain boundary expansion (m$\AA$)", "Grain boundary energy (mJ/m$^2$)", 
            image_name+"_fit_ve")


        fit = fit_and_plot(pressures_init, e_gbs_init,  pressures, e_gbs, power, 
            name, "External pressure (GPa)", "Grain boundary  energy (mJ/m$^2$)", 
            image_name+"_pe")
        #print fit
        ext_p_min  = fit.deriv().r[power-2] #external pressure in the minimum; derivative of function and the value of x in minimum

        fit_sxe = fit_and_plot(sigma_xx, e_gbs_init,  sigma_xx, e_gbs, power, 
            name, "Sigma xx (MPa)", "Grain boundary energy (mJ/m$^2$)", 
            image_name+"_sxe")
        sxe_min = fit_sxe.deriv().r[power-2] #sigma xx at the minimum of energy
        printlog( "sigma xx at the minimum of energy is", sxe_min," MPa")


        fit1 = fit_and_plot(pressures_init, gb_volumes,  pressures, gb_volumes, 1,
            name, "External pressure (GPa)", "Grain boundary expansion (m$\AA$)", 
            image_name+"_pv", lines = 'xy')
        #print fit1
        pulay = - calc[id].bulk_extpress
        #print " At external pressure of %.0f MPa; Pulay correction is %.0f MPa." % (ext_p_min+pulay, pulay)       
        #print " Egb = %.1f mJ m-2; Vgb = %.0f mA;"%(fit(ext_p_min), fit1(ext_p_min)  )
        printlog ("%s.fit.pe_pv & %.0f & %.0f & %0.f & %0.f \\\\" %
            (n[0]+'.'+n[1], fit(ext_p_min), fit1(ext_p_min),ext_p_min, ext_p_min+pulay   ))


        #print "\n At zero pressure with Pullay correction:"
        #print " Egb = %.1f mJ m-2; Vgb = %.0f mA; " % (fit(-pulay), fit1(-pulay))
        outstring =  ("%s.fit.pe_pv & %.0f & %.0f & %0.f & %0.f\\\\" %(n[0]+'.'+n[1], fit(-pulay), fit1(-pulay),-pulay,0    ))
        # print outstring
        calc[conv[0]].egb = fit(-pulay)
        calc[conv[0]].vgb = fit1(-pulay)

        return outstring #fit(ext_p_min), fit1(ext_p_min) 


    if type_of_plot == "fit_gb_volume":
        """
        should be rewritten using fit_and_plot() function
        """

        for id in conv:
            #energies.append(calc[id].energy_sigma0)
            #init_energies.append( calc[id].list_e_sigma0[0] ) 
            gb_volumes.append(calc[id].v_gb)
            e_gbs.append(calc[id].e_gb)
            e_gbs_init.append(calc[id].e_gb_init )           


        power = 3
        fit_ve = fit_and_plot(gb_volumes, e_gbs_init,  gb_volumes, e_gbs, power, 
            name, "Excess volume ($m\AA$)", "Twin energy ($mJ/m^2$)", 
            image_name+"_fit_ve")

        vgb_min  = fit_ve.deriv().r[power-2]


        #print "Fit of excess volume against energy. Pressure is uknown:"
        #print "Test Egb_min = %.1f mJ m-2; v_min = %.0f mA;"%(fit_ve(vgb_min), vgb_min)
        print ("%s.fit.ve & %.0f & %.0f & - & - \\\\" %
            (n[0]+'.'+n[1], fit_ve(vgb_min), vgb_min,   ))

    if type_of_plot == "fit_gb_volume2":

        for id in conv:
            energies.append(calc[id].energy_sigma0)
            init_energies.append( calc[id].list_e_sigma0[0] ) 
            volumes.append(calc[id].end.vol)
            pressures.append(calc[id].extpress )
            pressures_init.append(calc[id].extpress_init )

        power = 3
        pulay = 500

        fit_ve = fit_and_plot(volumes, init_energies,  volumes, energies, power, 
            name, "Volume ($\AA^3$)", "Energy  sigma->0 ($eV$)", 
            image_name+"_fit_ve")
        
        Vmin  = fit_ve.deriv().r[power-2] # minimum volume at the minimum energy
        Emin  = fit_ve(Vmin)

        fit_pe = fit_and_plot(pressures_init, init_energies,  pressures, energies, power, 
            name, "External pressure ($MPa$)", "Energy  sigma->0 ($eV$)", 
            image_name+"_fit_pe")

        ext_p_min  = fit_pe.deriv().r[power-2] #external pressure in the minimum; derivative of function and the value of x in minimum
        

        fit_pv = fit_and_plot(pressures_init, volumes,  pressures, volumes, 1,
            name, "External pressure ($MPa$)", "Volume of cell ($\AA^3$)", 
            image_name+"_fit_pv")


             
        atP = (" Emin = %.3f meV;  Vmin = %.0f A^3; "%( fit_pe(ext_p_min), fit_pv(ext_p_min)  )  ) + \
              (" for the minimum of energy relative to external pressure. The value of pressure is %.0f MPa; Pulay correction is %.0f MPa." % (ext_p_min+pulay, pulay) )
        
        at_zeroP = (" Emin = %.3f meV;  Vmin = %.0f A^3; " % (fit_pe(-pulay), fit_pv(-pulay) )  ) + \
                   (" the value of energy and volume at zero pressure with Pullay correction" )
        
        #print " Emin = %.3f meV;  Vmin = %.0f A^3;  for the minimum of energy relative to volume at some external pressure"%(fit_ve(Vmin), Vmin )
        #print atP
        #print at_zeroP

        printlog( "Compare V at -pulay and V for energy minimum", fit_pv(-pulay), Vmin, imp = 'y')

        return fit_pe(-pulay), fit_pv(-pulay), Emin, Vmin








    if type_of_plot == "kpoint_conv":
        energies = []
        kpoints = []
        times = []

        for id in list_of_calculations:
            if "4" not in calc[id].state:
                continue
            energies.append(calc[id].potenergy)
            kpoints.append(calc[id].kspacing[2])
            times.append(calc[id].time)

            name.append( id[1] )

        plt.figure()
        plt.title(name)
        plt.plot(kpoints, energies,'bo-')
        plt.ylabel("Total energy (eV)")
        plt.xlabel("KSPACING along 3rd recip. vector ($\AA ^{-1}$)")
        plt.twinx()
        plt.plot(kpoints,times,'ro-')
        plt.ylabel("Elapsed time (min)")
        plt.savefig('images/'+str(conv[0])+'kconv')


    if type_of_plot == "contour":
        alist = [] ;        clist = []
        nn = str(calc[conv[0]].id[0])+"."+str(calc[conv[0]].id[1])
        f = open("a_c_convergence/"+nn+"/"+nn+".out","w")
        f.write("END DATASET(S)\n")
        k = 1
        for id in conv: #Find lattice parameters and corresponding energies
            a = calc[id].a
            c = calc[id].c
            if a not in alist: alist.append(a); 
            if c not in clist: clist.append(c);
            f.write( "acell%i %f %f %f Bohr\n"%(k, calc[id].a/to_ang, calc[id].a/to_ang, calc[id].c/to_ang )   )
            #print "etotal%i %f\n"%(k, calc[id].energy_sigma0/to_eV ),
            k+=1;
        X,Y = np.meshgrid(alist, clist)
        Z = np.zeros(X.shape)
        Zinv = np.zeros(X.shape)

        
        k=1
        for i in range(len(alist)):
            for j in range(len(clist)):
                for id in conv:
                    if calc[id].a == alist[i] and calc[id].c == clist[j]:
                        Z[i][j] = calc[id].energy_sigma0
                        Zinv[j][i] = calc[id].energy_sigma0
                        f.write( "etotal%i %f\n"%(k, calc[id].energy_sigma0/to_eV )   )
                        k+=1
        f.write("+Overall time at end (sec) : cpu=     976300.2  wall=     976512.8")
        f.close

        #Make two plots for different a and c
        plt.figure()
        plt.title(name)
        for i in range(len(alist)):
            plt.plot(clist, Z[i],'o-',label='a='+str(alist[i]))
        plt.legend()
        plt.ylabel("Total energy (eV)")
        plt.xlabel("c parameter ($\AA$)")
        plt.savefig('images/'+str(conv[0])+'c')

        plt.figure()
        plt.title(name)
        for j in range(len(clist)):
            plt.plot(alist, Zinv[j],'o-',label='c='+str(clist[j]))
        plt.legend()
        plt.ylabel("Total energy (eV)")
        plt.xlabel("a parameter ($\AA$)")
        plt.savefig('images/'+str(conv[0])+'a')

        #Make contour
        plt.figure()
        cf = plt.contourf(X, Y, Z, 20,cmap=plt.cm.jet)
        cbar = plt.colorbar(cf)
        cbar.ax.set_ylabel('Energy (eV)')

        plt.xlabel('$a$ ($\AA$)')
        plt.ylabel('$c/a$')

        plt.legend()
        plt.savefig('images/ru-contourf.png')
        #plt.show()


        #Make equation of state
        eos = EquationOfState(clist,Z[2])
        v0, e0, B = eos.fit()
        #print "a = ", alist[2]
        printlog( '''
        v0 = {0} A^3
        E0 = {1} eV
        B  = {2} eV/A^3'''.format(v0, e0, B) )
        eos.plot('images/a[2]-eos.png') 

        eos = EquationOfState(alist,Zinv[2])
        v0, e0, B = eos.fit()
        #print "c = ", clist[2]
        printlog( '''
        v0 = {0} A^3
        E0 = {1} eV
        B  = {2} eV/A^3'''.format(v0, e0, B) )
        eos.plot('images/c[2]-eos.png') 


    if type_of_plot == "dimer":

        if not cl:
            cl =  calc[list_of_calculations[0]]

        x1 = [] #list of distances


        
        if cl.end.natom > 2:
            raise RuntimeError


        # print (cl.end.list_xcart)
        for xcart in cl.end.list_xcart:
            # x = xcart[1]
            # d = (x[0]**2 + x[1]**2 + x[2]**2)**0.5
            d = np.linalg.norm(xcart[1]-xcart[0]) #assuming there are only two atoms
            print(xcart[0], xcart[1], d)

            x1.append(d)

        y1 = cl.list_e_without_entr
        power = 4
        name = 'dimer'
        xlabel = 'Bond length'
        ylabel = 'Full energy'
        # print(x1, y1)
        coeffs1 = np.polyfit(x1, y1, power)        
      
        fit_func1 = np.poly1d(coeffs1)

        x_range = np.linspace(min(x1), max(x1))
        fit_y1 = fit_func1(x_range); 
        f = fit_func1.deriv()
        min_e = fit_func1(f.r[2]).real
        printlog("The minimum energy per atom and optimal length of dimer are {:.3f} eV and {:.3f} A".format( min_e/2., f.r[2].real), imp = 'Y' )
        try:
            printlog("The atomization energy for dimer is {:.3f} eV ; The energy of atom in box is taken from the provided b_id".format(min_e - 2*cl.e_ref), imp = 'Y' )
        except:
            print('Reference energy was not found')
        plt.figure()
        plt.title(name)
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        plt.plot(x1, y1, 'ro', label = 'init')
        plt.plot(x_range, fit_y1, 'r-', label = 'init_fit')
        
        if filename:
            path2saved, path2saved_png = process_fig_filename(filename, fig_format)


        if plot:
            plt.show()

    return











