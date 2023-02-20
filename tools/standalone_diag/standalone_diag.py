#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# MODEL  BFM - Biogeochemical Flux Model
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#
# DESCRIPTION
#   Python script to generate plots of BFM STANDALONE results
#   with HTML index file in the experiment directory
#
# USAGE
#   Install "bfm_standalone_diag" python conda environment
#      `conda env create -f environment.yml`
#   Execution: python standalone_diag.py <bfm-output.nc>
#
# COPYING
#
#   Copyright (C) 2023 BFM System Team (bfm_st@cmcc.it)
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation.
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#   See the GNU General Public License for more details.
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

from matplotlib import use
use("Agg")
from ncPyTools.netCDFView import netCDFView
from sys import argv, exit
from matplotlib.pyplot import plot, figure, close, legend, ylabel
from pdb import set_trace
from netCDF4 import num2date
from numpy import ones
import os

# elemental mass
CMass = 12.0


def diagnostics():
    # plot resolution
    res = 200

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
        's': 'silica',
        'o': 'oxygen',
        'l': 'chlorophyll'
    }

    # compartments dict
    compartments = {'pel': 'Pelagic', 'ben': 'Benthic', 'ice': 'Seaice'}

    # initialize component dicts
    layout = {}
    layout['pel'] = {
        'inorganic': ['O3', 'O2', 'N1', 'N3', 'N5'],
        'phytoplankton': ['P'+str(x) for x in range(1,10)],
        'zooplankton': ['Z'+str(x) for x in range(1,10)],
        'bacteria': ['B1'],
        'detritus': ['R1', 'R6'],
        'gross production': ['ruPPY', 'ruZOO'],
        'respiration': ['resPPY', 'resPBA', 'resZOO'],
        'air-sea flux': ['jsurO3', 'jsurO2'],
        'calcite': ['O5'],
        'quotas': ['phytoplankton', 'zooplankton', 'bacteria', 'detritus']
    }

    layout['ben'] = {
        'inorganic': ['K1', 'K11', 'K4', 'K14', 'G2'],
        'living': ['Y1', 'Y2', 'Y3', 'Y4', 'Y5', 'H1', 'H2'],
        'detritus': ['Q1', 'Q6'],
        'quotas': ['detritus']
    }

    layout['ice'] = {
        'inorganic': ['F2', 'F3', 'I1', 'I3', 'I4', 'I5'],
        'phytoplankton': ['S1', 'S2'],
        'zooplankton': ['X1'],
        'bacteria': ['T1'],
        'detritus': ['U1', 'U6'],
        'quotas': ['phytoplankton', 'zooplankton', 'bacteria', 'detritus']
    }

    # initialize HTML dict
    html = {}
    for cons in constituents.keys():
        html[cons] = {}
        for comp in compartments.keys():
            html[cons][comp] = {'name': [], 'file': []}
    html['quotas'] = {}
    for comp in compartments.keys():
        html['quotas'][comp] = {'name': [], 'file': []}

    # Process input file
    with netCDFView(argv[1], Quiet=True) as nc:
        # State Variables
        for cons in constituents.keys():
            for comp in compartments.keys():
                for grp in layout[comp].keys():
                    if grp == 'quotas': continue
                    f = figure()
                    leg = []
                    for v in layout[comp][grp]:
                        vcon = v + cons
                        pl(nc, vcon, grp, leg)
                    # save plot
                    if leg:
                        filename = "images/{}-{}-{}.png".format(
                            cons, comp, grp).replace(' ','-')
                        f.savefig(outpath + '/' + filename, dpi=res)
                        html[cons][comp]['name'].append(grp)
                        html[cons][comp]['file'].append(filename)
                    close(f)

        # Quotas
        for comp in compartments.keys():
            for grp in layout[comp]['quotas']:
                for v in layout[comp][grp]:
                    f = figure()
                    leg = []
                    plq(nc, v, grp, leg)
                    if leg:
                        filename = "images/quotas-{}-{}.png".format(v, comp).replace(' ','-')
                        f.savefig(outpath + '/' + filename, dpi=res)
                        html['quotas'][comp]['name'].append(' '.join(
                            [v, '-', grp]))
                        html['quotas'][comp]['file'].append(filename)
                    close(f)

        # HTML index
        create_index(curpath, datapath, constituents, compartments, html,
                     outpath)

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


def pl(nc, var, group, leg):
    ''' Plot state variables of groups '''
    dates = num2date(nc("time"),
                     nc.variables["time"].units,
                     only_use_cftime_datetimes=False)
    # plot
    if var in nc.variables:
        data = nc(var)
        leg.append(var)
        units = nc[var].units
        if var == 'O3c':
            data = data * (1. / CMass)
            units = 'umol/kg'
        plot(dates, data)
        if var == 'O3c' and 'O3h' in nc.variables:
            plot(dates, nc('O3h'))
            leg.append('O3h')
        elif var == 'N3n' and 'N4n' in nc.variables:
            plot(dates, nc('N4n'))
            leg.append('N4n')
        if leg:
            legend(leg)
            ylabel(group.upper() + ' ['+ units + ']')
    return


def plq(nc, var, group, leg):
    ''' Plot quota values of constituents '''
    dates = num2date(nc("time"),
                     nc.variables["time"].units,
                     only_use_cftime_datetimes=False)
    N = False
    P = False
    sc = var + 'c'
    if sc in nc.variables:
        sn = var + 'n'
        sp = var + 'p'
        if sn in nc.variables:
            N = True
            plot(dates, nc(sc) / nc(sn) / CMass)
            leg.append("C/N")
        else:
            print("Skipping {} C/N".format(var))
        if sp in nc.variables:
            P = True
            plot(dates, nc(sc) / nc(sp) / CMass)
            leg.append("C/P")
        else:
            print("Skipping {} C/P".format(var))
        if N and P:
            N2P = nc(sn) / nc(sp)
            plot(dates, nc(sn) / nc(sp))
            leg.append("N/P")
        if leg:
            # plot also Redfield ratios
            plot(dates, 106. * ones(dates.shape), "k:")
            leg.append("Redfield")
            legend(leg, loc=1)
            plot(dates, 16. * ones(dates.shape), "k:")
            plot(dates, 106. / 16. * ones(dates.shape), "k:")
            ylabel(' '.join([var.upper(), '-', group.upper(), '[ratio]']))
    return


#
#==========================================================================
# Main sentinel
#==========================================================================

if __name__ == "__main__":
    diagnostics()
