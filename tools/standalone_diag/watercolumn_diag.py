# Python script to generate plots of BFM STANDALONE results
# with HTML index file in the experiment directory
#
# T.Lovato for BFM System Team (2020)
#
# Prerequisite: install "bfm_standalone_diag" python conda environment
#   `conda env create -f environment.yml`
#
# Usage: python standalone_diag.py <bfm-output.nc>
#

from matplotlib import use
#use("Agg")
from ncPyTools.netCDFView import netCDFView
from sys import argv, exit
#from matplotlib.pyplot import plot, contourf, figure, close, legend, ylabel, axes
import matplotlib.pyplot as plt
from pdb import set_trace
from netCDF4 import num2date
import numpy as np
import os

# elemental mass
CMass = 12.0

# plot resolution
res = 200



def diagnostics():

    # script path
    curpath = os.path.dirname(os.path.abspath(argv[0]))

    # check input arguments
    if len(argv) < 2:
        print('Need to provide BFM output filename as argument for the script')
        exit(1)

    print('Begin processing file: ' + argv[1])

    # input data path
    datapath = os.path.dirname(os.path.abspath(argv[1]))

    # create output folder in experiment
    outpath = datapath + '/diagnostics_html'
    if not os.path.exists(outpath):
        os.makedirs(outpath + '/images')

    # constituents dict
    constituents = {
        'c': 'carbon',
        'n': 'nitrogen',
        'p': 'phosphorus',
    }

    # compartments dict
    compartments = {'pel': 'Pelagic', 'ben': 'Benthic', 'ice': 'Seaice'}

    # initialize component dicts
    layout = {}
    layout['pel'] = {
        'inorganic': ['N1','N3'],
        'phytoplankton': ['Chla'],
        'zooplankton': ['Z4', 'Z5', 'Z6'],
    }

    layout['ben'] = {
    }

    layout['ice'] = {
    }

    # initialize HTML dict
    html = {}
    for cons in constituents.keys():
        html[cons] = {}
        for comp in compartments.keys():
            html[cons][comp] = {'name': [], 'file': []}

    # Process input file
    with netCDFView(argv[1], Quiet=True) as nc:
        # State Variables
        for cons in constituents.keys():
            for comp in compartments.keys():
                for grp in layout[comp].keys():
                    for v in layout[comp][grp]:
                        vcon = v + cons
                        if v in ['Chla']:
                            vcon = v
                        if vcon in nc.variables:
                            filename = "images/{}-{}-{}.png".format(
                                cons, comp, vcon).replace(' ','-')
                            plot_clim(nc, vcon, outpath + '/' + filename)
                            #html[cons][comp]['name'].append(grp)
                            #html[cons][comp]['file'].append(filename)

        # HTML index
        #create_index(curpath, datapath, constituents, compartments, html,
        #             outpath)

    print(
        'Processing completed. \n\nType \'open diagnostics_html/index.html\' or inspect content of diagnostics_html/images.\n'
    )

    return


#
#==========================================================================
# Functions
#==========================================================================


def create_index(curpath, datapath, constituents, compartments, html, outpath):
    ''' create HTML index from template '''
    template = curpath + '/index.template'

    # Read in the template file
    with open(template, 'r') as file:
        filedata = file.read()

    # add experiment path
    filedata = filedata.replace('_EXP_PATH_', datapath)

    # fill in navigation bar
    navs = []
    for cons in constituents.keys():
        navs.append(constituents[cons].upper())
    navs.append('QUOTAS')
    txt = ''
    for nn in navs:
        txt += '\t\t<a href="#' + nn + '">' + nn + '</a>\n'
    filedata = filedata.replace('_NAV_ITEMS_', txt)

    # fill in plots by constituent
    txt = ''
    for cons in html.keys():
        cname = constituents[cons].upper() if cons != 'quotas' else 'QUOTAS'
        txt += '<div id="' + cname + '" class="anchor"> </div>\n'
        txt += '<div class="constituent">\n<h2>' + cname + '</h2>\n'
        # by compartment
        for comp in html[cons].keys():
            if html[cons][comp]['name']:
                cc = html[cons][comp]
                txt += '<div class="compartment">\n'
                txt += '\t<h4>' + compartments[comp] + '</h4>\n'
                for pname, fname in zip(cc['name'], cc['file']):
                    txt += '\t<div class="plot">\n'
                    txt += '\t\t<a href="' + fname + '" ><img class="image" src="' + fname + '" alt="" /></a>\n'
                    txt += '\t\t<h3>' + pname + '</h3>\n\t</div>\n'
                txt += '</div>\n'
        txt += '</div>\n<div style="clear: left; padding-bottom: 2em;"> </div>\n\n'
    filedata = filedata.replace('_PLOT_ITEMS_', txt)

    # Write HTML file to data folder
    with open(outpath + '/index.html', 'w') as file:
        file.write(filedata)

    return


def plot_clim(nc, var, filename):
    ''' Plot state variables of groups '''
    months=['J','F','M','A','M','J','J','A','S','O','N','D']
    #dates = num2date(nc("time"),
    #                 "seconds since 01-01-2000",
    #                 calendar='360_day',
    #                 only_use_cftime_datetimes=False)
    # axes
    dates = range(14)
    depths = nc("z")
    #depths = depths * -1.
    xx, yy = np.meshgrid(dates, depths)
    # data clim
    clim = np.zeros([len(depths), len(dates)])
    data = nc(var)
    name = nc[var].long_name
    units = nc[var].units
    #if var == 'Z5c':
    #    data = nc(var) + nc('Z6c')
    #    name = 'Microzoo. + Het. nanop.'
    # get last 5 yrs
    aaa = data[int(len(data)*0.5):]
    for ii in range(12):
        sel = slice(ii,len(aaa),12)
        clim[:,ii+1] = np.mean(aaa[sel,:], axis=0)

    # add bound replica
    clim[:,0] = clim[:,12]
    clim[:,13] = clim[:,1]

    # fixes
    units = nc[var].units
    if var == 'O3c':
        data = data * (1. / CMass)
        units = 'umol/kg'


    # plot 
    f, ax = plt.subplots()
    cmap = 'jet' 
    # set ranges
    if var in ['Z4c']:
      vmin = 2.
      vmax = 6.
    elif var in ['Z5c']:
      vmin = 10.
      vmax = 30.
    if var in ['N1p']:
      vmin = 0.
      vmax = 0.15
    else:
      vmin = 0.
      vmax = np.ceil(np.max(clim))

    levels = np.linspace(vmin, vmax, 100+1)
    cc = plt.contourf(xx, yy, clim, levels, cmap=cmap)
    cbar = plt.colorbar(cc, ticks=levels[::25])
  
    ax.set_yticks([0., -5., -10., -15.])
    # xaxis
    ax.set_xlim([0.5, 12.5])
    ax.set_xticks(range(1,13))
    ax.set_xticklabels(months, fontweight='bold')
    
    # 
    ax.set_title(name+' ['+units+']',fontsize=14, fontweight='bold')
    ax.set_ylabel('Depth [m]', fontsize=12)

    f.savefig(filename, dpi=res)
    plt.close()
    
    return



#
#==========================================================================
# Main sentinel
#==========================================================================

if __name__ == "__main__":
    diagnostics()
