# -*- coding: utf-8 -*- 

from header import *
from functions import write_xyz, replic
# import header
# from operator import itemgetter
from classes import res_loop , add_loop
# from pairs import 
# from functions import image_distance, local_surrounding
# from chargeden.functions import chg_at_point, cal_chg_diff
# from dos.functions import plot_dos
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = "stix"
import matplotlib.gridspec as gridspec
import scipy

def fit_and_plot(power = None, xlabel = "xlabel", ylabel = "ylabel", image_name = None,
    xlim = None, ylim = None, title = None, figsize = None,
    xlog = False,ylog = False, scatter = False, legend = False, markersize = 10,
    **data):
    """Should be used in two below sections!
    Creates one plot with two dependecies and fit them;
    return minimum fitted value of x2 and corresponding value of y2; 
    if name == "" image will not be plotted
    power - the power of polynom

    data - each entry should be (X, Y, 'r-')
    """

    # print data

    if 1:

        
        plt.figure(figsize=figsize)
        if title: plt.title(title)
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)



        scatterpoints = 1
        for key in sorted(data):



            if scatter:
                
                plt.scatter(data[key][0], data[key][1],  s = data[key][2], c = data[key][-1], alpha = 0.8, label = label)
            else:
                try:
                    label = data[key][3]
                except:
                    label = key
                print 'label is ', label


                plt.plot(data[key][0], data[key][1], data[key][2], linewidth = 2, label = label, markersize = markersize, alpha = 0.8)






        plt.axvline(color='k')
        if xlim: 
            plt.xlim(xlim)
            # axes = plt.gca()
            # axes.set_xlim([xmin,xmax])
            # axes.set_ylim([ymin,ymax])
        if ylim:
            plt.ylim(ymin=ylim[0])
            if ylim[1]: plt.ylim(ymax=ylim[1])


        if power:
            for key in data:
                coeffs1 = np.polyfit(data[key][0], data[key][1], power)        
                
                fit_func1 = np.poly1d(coeffs1)
                x_range = np.linspace(min(data[key][0]), max(data[key][0]))
                fit_y1 = fit_func1(x_range); 
         
                plt.plot(x_range, fit_y1, data[key][-1][0], )

                # x_min  = fit_func2.deriv().r[power-2] #derivative of function and the second cooffecient is minimum value of x.
                # y_min  = fit_func2(x_min)
                slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(data[key][0], data[key][1])
                print 'R^2 = ', r_value**2, key


        if xlog: plt.xscale('log')
        if ylog: 
            if "sim" in str(ylog):
                plt.yscale('symlog', linthreshx=0.1)
            else:
                plt.yscale('log')



        if legend: 
            plt.legend(loc = legend, scatterpoints = scatterpoints)


        plt.tight_layout()
        if image_name:
            try:
                path_to_images
            except:
                path_to_images = ''
            print "Saving image ... to ", path_to_images+str(image_name)+'.png'
            plt.savefig(path_to_images+str(image_name)+'.png', dpi = 300, format='png')
            plt.close()
        else:
            plt.show()


    return 


def plot_bar(xlabel = "xlabel", ylabel = "ylabel",
    xlim = None, ylim = None,
    image_name = None, title = None, bottom = 0.18, hspace = 0.15, barwidth = 0.2,
    data1 = [],data2 = [],data3 = [],data4 = [],
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
            print 'color', data[key][2]
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





    elif data1 and not data4:
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
            print ind+width
            print data[0][0]
            shift = 0.2

            for d in data:
                ax.bar(ind+shift, d[1], width, color = d[2], label = d[-1])# yerr=menStd)
                # print ind
                shift+=width

            ax.axhline(y=0, color='black')
            # ax.set_xticklabels(xlabels , rotation=70 )

            ax.set_xticks(ind+width+len(data)*width/2)
            

        ax1.set_xticklabels(data1[0][0] , rotation=80 )
        ax2.set_xticklabels(data2[0][0] , rotation=80 )

        ax1.set_ylabel(ylabel)

        ax1.legend(loc=2, )
        ax1.axis('tight')
        ax2.axis('tight')


    # ax.set_yscale('log')
    # plt.yscale('symlog', linthreshx=0.1)

    # ax.set_title('Scores by group and gender')


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

    if image_name:
        # print "Saving image ..."
        plt.savefig(path_to_images+str(image_name)+'.png', dpi = 200, format='png')
    else:
        plt.show()

    return





def plot_and_annotate(power = 2, xlabel = "xlabel", ylabel = "ylabel", image_name = None,
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
            plt.savefig(path_to_images+str(image_name)+'.png', dpi = 300, format='png')
            plt.close()
        else:
            plt.show()


    return 


def plot_bar_simple(xlabel = "xlabel", ylabel = "ylabel",
    xlim = None, ylim = None,
    image_name = None, title = None, 
    data = []):

    width = 0.5       # the width of the bars



    plt.figure(figsize=(10,5))   #5:7 ratio for A4, 
    fig, ax = plt.subplots()


    N = len(data[0])
    xlabels = data[0]
    # print xlabels
    ind = np.arange(N)  # the x locations for the groups
    shift = 0
    # for d in data:
    d = data
    print d[2]

    rects = ax.bar(ind+shift, d[1], width, color = d[2], label = d[-1])# yerr=menStd)
    rects[0].set_color('g')
    rects[1].set_color('b')

    # print ind
    shift+=width

    ax.axhline(y=0, color='black')

    ax.set_xticks(ind)#+width)
    ax.set_xticklabels(xlabels , rotation=20 )
        
    ax.set_ylabel(ylabel)

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
            ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
                    ha='center', va='bottom')
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














def plot_stress_gb_v():


    res_loop('t111g','9292',range(1,11), calc, conv, varset, 'gbep',('t111b_r','8302',1), readfiles = False) #Fig 2 and 3

    return




def plot_gb_v():
    """1. Grain boundary energy and volume"""
    #Gb_e - gb_v
    name = ["T1m", "T1g",  "T2",  "C1", "S7",]# "S7" ]
    gb_e = [351, 359, 288, 407, 730, ]#1084] 
    gb_v = [43,  39,  23,  70,  59 , ]#171 ]
      

    if 1:
        xlabel = "Grain boundary expansion (m${\AA}$)"; ylabel = "Grain boundary energy (mJ/m$^2$)"
    else:
        xlabel = u"Зернограничный объём (м${\AA}$)"; ylabel = u"Зернограничная энергия (мДж/м$^2$)"

    if 1:
        plot_and_annotate(image_name = 'gb_e_v',
            xlabel = xlabel, ylabel =ylabel,
                        title = None ,
                        # fit = 'cubic',          
                        xlim = (10, 72), 
                        ylim = (200, 820),#0.4
                        Ti   = (gb_v, gb_e, name, 'ro'), 
                        )
    return








def produce_structure_pictures(calc):
    """Includes coordinates of impurities to show segregation positions"""
    posT1m = [
    [20.72385349,   4.1269786 ,   5.89584523 , 'Be' ],
    [ 18.93552982,   2.48843802,   5.8957184 , 'Be'],
    [ 20.71714221,   0.68616663,   1.47209765, 'Be'],
    [ 19.13746512,   5.88883189,   1.47480645, 'Be'],
    ]
    posT1g = [
    [ 20.9384873     , 0.67306254       ,3.58141074, 'Be' ],
    [  2.10061489e+01,   6.68137381e-01 ,  1.03523281e-02, 'Be' ],
    [ 20.6806174     , 4.12415131       ,2.93928801, 'Be' ],
    [ 19.14006886    , 2.42215779       ,1.47405535, 'Be' ],
    [ 22.53515212    , 2.44143107       ,2.94816354, 'Be' ],
    [  1.91505721e+01,   5.90231582e+00 ,  1.13720244e-04, 'Be' ],
    [  2.06783544e+01,   4.12136016e+00 , -2.63942403e-03, 'Be' ],
    ]
    posT2 = [
    [ 15.78118505 ,  4.78527643 ,  5.11735029, 'Be'],
    [ 14.74191762 ,  2.50200925 ,  4.41093292, 'Be'],
    [ 14.06713264 ,  0.2067005  ,  4.48418856, 'Be'],
    [ 13.40773994 ,  2.86328286 ,  1.93928418, 'Be'],

    ]
    posC1 = [
    [ 23.6121777  ,  2.59936859 ,  1.58603275, 'Be'],
    [ 21.76417242 ,  0.61791592 ,  1.75450068, 'Be'],
    ]

    posS7 = [
    [15.44238904,   6.7381076 ,   2.03964221, 'Be'],
    [20.38037862,   6.73610186,   3.38817988, 'Be'],
    [16.44283507,   3.86499065,   1.76279741, 'Be'],
    [16.44873277,   3.92762434,   4.13582109, 'Be'],
    [15.40052276,   6.73977648,   4.33879376, 'Be'],
    [17.56441714,   1.22090821,   1.95250884, 'Be'],
    ]

    posC1bulk = [ [15.68533491,   1.77127397 ,  1.69868483, 'Be'], ] 

    """Final picture"""
    # write_xyz( replic(calc[('c1g','929',1)].end, (1,2,1)), repeat = 2 , shift = 0.85, gbpos2 = 23.7, gbwidth = 2.7, imp_positions = posC1,  rotate = 90,
    #     full_cell = 1,     include_boundary = (0,2))

    # write_xyz( replic(calc[('t111g','9292',1)].end, (1,2,1)), repeat = 2, shift = 0.8, gbpos2 = 20.9, gbwidth = 4.1, imp_positions = posT1m,  rotate = 90,
    #     full_cell = 1,     include_boundary = (0,2))

    # write_xyz( replic(calc[('t111sg','93kp9',1)].end, (1,2,1)), repeat = 2, shift = 0.8, gbpos2 = 20.9, gbwidth = 4.1, imp_positions = posT1g,  rotate = 90, 
    #     specialcommand = 'select Be2/1, Be7/1\nlabel off\n', full_cell = 1, include_boundary = (2,0)  ) #


    write_xyz( replic(calc[('t21g','93',1)].end, (1,1,1)), repeat = 2, shift = 0.89, gbpos2 = 15.726, gbwidth = 3.3 , imp_positions = posT2, rotate = 90,
        full_cell = 1, include_boundary = (0,2))

    
    # write_xyz( replic(calc[('csl71sg10','93',1)].end, (1,1,1)), repeat = 2, shift = 1.1 , gbpos2 = 17.4865, gbwidth = 3.9, rotate = 90, imp_positions = posS7,
    # full_cell = 1, include_boundary = (0,2.5), specialcommand = 'select Be1/1, Be3/1\nlabel off\n'     )


    #For presentation: 
    # write_xyz( replic(calc[('c1g','929',1)].end, (1,2,1), inv = -1), 
    #     repeat = 2 , shift = 0.9, gbpos2 = 23.7, gbwidth = 2.7, 
    #     imp_positions = posC1bulk, topview = 0, file_name = 'C1_bulk',
    #     specialcommand = 'select Be*\ncpk 400\ncolor red\n')
    # write_xyz( replic(calc[('c1g','929',1)].end, (1,2,1), inv = -1), 
    #     repeat = 2 , shift = 0.9, gbpos2 = 23.7, gbwidth = 2.7, 
    #     imp_positions = posC1[0:1], topview = 0, file_name = 'C1_1',
    #     specialcommand = 'select Be*\ncpk 400\ncolor red\n')
    


    # # write_xyz( replic(calc[('csl71g','93',3)].end, (1,2,1)), repeat = 2, shift = 0.70 )
    # write_xyz( replic(calc[('t112g','93kp9',1)].end, (1,1,1)), repeat = 2, shift = 1.1 )
    return