import xarray as xr
import numpy as np
import json
from astropy.utils.misc import JsonCustomEncoder
import astropy.units as u
import os
import pandas as pd
import glob
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotting as niko


def get_individual_planet_data(current_dir,planet_name,planet_letter,event):
    """
    This function loads all the data for an individual planet and event type. 

    Parameters
    ----------
    current_dir : str 
        This defines the current directory you are in and assumes that there is a path structure 
        current_dir/../Transformed Data
    planet_name : str 
        This is the name of the host star e.g. "WASP-39" for WASP-39 b
    planet_letter : str 
        This is the letter e.g. "b" for WASP-39 b
    event : str 
        Transit or Eclipse 

    Returns
    -------
    dictionary with [x,y,e1, and e2] for all the data available for that planet event 
    """
    current_dir = os.getcwd()
    formatted_data_dir = '../Transformed Data/'    
    full_path = os.path.join(current_dir,formatted_data_dir, planet_name, planet_letter,event)
    
    if event == "Transit":
        to_plot='transit_depth'
    elif event == 'Eclipse':
        to_plot='eclipse_depth'
    
    
    #what xarrays are in this directory
    all_files = glob.glob(os.path.join(full_path,'*nc'))
    
    return_dictionary = {}
    #let's loop through and load all the xarrays for one planet
    for ifile in all_files:
        
        ds = xr.open_dataset(ifile)
        
        #lets label them by observing mode for now
        ys = [i for i in ds.data_vars.keys() if 'error' not in i]
        xs = [ds.data_vars[i].dims[0] for i in ds.data_vars.keys() if 'error' not in i]
        es1 = [to_plot+'_'+i.split(to_plot)[0]+'error1'+i.split(to_plot)[-1] for i in ys]
        es2 = [to_plot+'_'+i.split(to_plot)[0]+'error2'+i.split(to_plot)[-1] for i in ys]
    
        #CONSIDER CHANGING TO if "error1" and error2 always come after y name
        #es1 = [i+'_error1' for i in ys]
        #es2 = [i+'_error2' for i in ys]
        
        #read out the data 
        for i,key in enumerate(ys):
            return_dictionary[key] = [ds.coords[xs[i]].values,
                                     ds.data_vars[key].values,
                                     ds.data_vars[es1[i]].values,
                                     ds.data_vars[es2[i]].values]
        
    return return_dictionary

def plot_single_parameter_planet_spectra(parameter, lower_limit, upper_limit, data_base, planet_name, planet_letter, event):
    """This function plots the transmission spectra of all planets within the parameters
        Parameters
        ----------
        parameter : str
            This is the attribute of a planet with a float value e.g. 'planet_radius: 0.9' planet_radius is the string that has a float 0.9. This will not work for attributes that are strings with strings e.g. 'molecules_significant: H2O'.
        lower_limit : float
            This is the lower end of the range for your parameter. The range looks like [lower_limit, upper_limit] e.g. 'if lower_limit = 0.5 then the range looks like [0.5, upper_limit]'. If you don't want a lower limit, let lower_limit = 0
        upper_limit : float
            This is the upper end of the range for your parameter. The range looks like [lower_limit, upper_limit] e.g. 'if upper_limit = 2.5 then the range looks like [lower_limit, 2.5]'. If you don't want an upper limit, let upper_limit = np.inf
        data_base : str
            This is all the planet files which has the transmission spectra and attributes for each planet. These should be xarray files that all end in 'nc'.
        planet_name : str 
            This is the name of the host star e.g. "WASP-39" for WASP-39 b. If you want all planets/host stars, let planet_name = '*'.
        planet_letter : str 
            This is the letter e.g. "b" for WASP-39 b. If you want all planets, let planet_letter = '*'.
        event : str 
            Transit or Eclipse
        Returns
        -------
        a single plot of all trasnsmission spectra within the range and parameter .
        """
    plt.xlabel('Wavelength')
    plt.ylabel('Event Depth')
    plt.title(f'Transmission Spectra for {parameter} [{lower_limit},{upper_limit}]')

    return_dictionary = {}

    for i in data_base:
        ds = xr.open_dataset(i)
    
        #event type and it's corresponding variable
        if ds.attrs['Event'] == 'Transit':
            to_plot = 'transit_depth'
        elif ds.attrs['Event'] == 'Eclipse':
            to_plot = 'eclipse_depth'
        elif ds.attrs['Event'] == 'PhaseC':
            to_plot = 'eclipse_depth'

        #check if parameter is within the desired range
        if float(ds.attrs[parameter]) > lower_limit-0.1 and float(ds.attrs[parameter]) < upper_limit+0.01:
        
            #files with single coord
            if len(ds.coords) == 1:
                ys = [var for var in ds.data_vars.keys() if 'error' not in var]
                xs = [ds.data_vars[var].dims[0] for var in ys]  # Assuming all variables share the same coordinate
                es1 = [to_plot + '_' + var.split(to_plot)[0] + 'error1' + var.split(to_plot)[-1] for var in ys]
                es2 = [to_plot + '_' + var.split(to_plot)[0] + 'error2' + var.split(to_plot)[-1] for var in ys]
            
                #dictionary for each file with 1 coord
                for i, key in enumerate(ys):
                    return_dictionary[key] = [
                        ds.coords[xs[i]].values,
                        ds.data_vars[key].values,
                        ds.data_vars[es1[i]].values,
                        ds.data_vars[es2[i]].values
                    ]

            #files with multiple coords (modes)
            else: 
                
                for coord in ds.coords:
                    ys = [var for var in ds.data_vars.keys() if 'error' not in var]
                    xs = [ds.data_vars[var].dims[0] for var in ys]  # Assuming all variables share the same coordinate
                    es1 = [to_plot + '_' + var.split(to_plot)[0] + 'error1' + var.split(to_plot)[-1] for var in ys]
                    es2 = [to_plot + '_' + var.split(to_plot)[0] + 'error2' + var.split(to_plot)[-1] for var in ys]

                    #dictionary for each mode in the file
                    for i, key in enumerate(ys):
                        return_dictionary[key] = [
                            ds.coords[xs[i]].values,
                            ds.data_vars[key].values,
                            ds.data_vars[es1[i]].values,
                            ds.data_vars[es2[i]].values
                        ]

            #plot all data on a single plot
            for key in ys:
                x, y, e1, e2 = return_dictionary[key]
                plt.errorbar(x, y, yerr=[e1, np.abs(e2)], linestyle='', marker='o', label=key)

    
    plt.legend()
    plt.show()


def interactive_plot_yscaled_planet_data(parameter, lower_limit, upper_limit, data_base, planet_name, planet_letter, event):
    """This function plots the transmission spectra (yscaled) of all planets within the parameters 
        Parameters
        ----------
        parameter : str
            This is the attribute of a planet with a float value e.g. 'planet_radius: 0.9' planet_radius is the string that has a float 0.9. This will not work for attributes that are strings with strings e.g. 'molecules_significant: H2O'.
        lower_limit : float
            This is the lower end of the range for your parameter. The range looks like [lower_limit, upper_limit] e.g. 'if lower_limit = 0.5 then the range looks like [0.5, upper_limit]'. If you don't want a lower limit, let lower_limit = 0
        upper_limit : float
            This is the upper end of the range for your parameter. The range looks like [lower_limit, upper_limit] e.g. 'if upper_limit = 2.5 then the range looks like [lower_limit, 2.5]'. If you don't want an upper limit, let upper_limit = np.inf
        data_base : str
            This is all the planet files which has the transmission spectra and attributes for each planet. These should be xarray files that all end in 'nc'.
        planet_name : str 
            This is the name of the host star e.g. "WASP-39" for WASP-39 b. If you want all planets/host stars, let planet_name = '*'.
        planet_letter : str 
            This is the letter e.g. "b" for WASP-39 b. If you want all planets, let planet_letter = '*'.
        event : str 
            Transit or Eclipse
        Returns
        -------
        a single plot of all trasnsmission spectra within the range and parameter .
        
    """
        
    database = glob.glob(f'../Transformed Data/{planet_name}/{planet_letter}/{event}/*.nc')
    
    
    return_dictionary = {}
    
    for i in database:
        ds = xr.open_dataset(i)
        
        # Determine the type of event and the corresponding variable
        if ds.attrs['Event'] == 'Transit':
            to_plot = 'transit_depth'
        elif ds.attrs['Event'] == 'Eclipse':
            to_plot = 'eclipse_depth'
        elif ds.attrs['Event'] == 'PhaseC':
            to_plot = 'eclipse_depth'
        
        # Check if the planet's radius is within the desired range
        if lower_limit < float(ds.attrs[parameter]) < upper_limit:
    
            # Extract the planet name from the file path
            path = i.split(os.sep)
            planet_full_name = f'{path[-4]} {path[-3]}'
            
            # Files with single coordinate
            if len(ds.coords) == 1:
                ys = [var for var in ds.data_vars.keys() if 'error' not in var]
                xs = [ds.data_vars[var].dims[0] for var in ys]
                es1 = [to_plot + '_' + var.split(to_plot)[0] + 'error1' + var.split(to_plot)[-1] for var in ys]
                es2 = [to_plot + '_' + var.split(to_plot)[0] + 'error2' + var.split(to_plot)[-1] for var in ys]
                
                # Store data for each variable
                for j, key in enumerate(ys):
                    return_dictionary[planet_full_name] = [
                        ds.coords[xs[j]].values,
                        ds.data_vars[key].values,
                        ds.data_vars[es1[j]].values,
                        ds.data_vars[es2[j]].values
                    ]
            
            # Files with multiple coordinates (modes)
            else:
                for coord in ds.coords:
                    ys = [var for var in ds.data_vars.keys() if 'error' not in var]
                    xs = [ds.data_vars[var].dims[0] for var in ys]
                    es1 = [to_plot + '_' + var.split(to_plot)[0] + 'error1' + var.split(to_plot)[-1] for var in ys]
                    es2 = [to_plot + '_' + var.split(to_plot)[0] + 'error2' + var.split(to_plot)[-1] for var in ys]
    
                    # Store data for each variable in different modes
                    for j, key in enumerate(ys):
                        return_dictionary[planet_full_name] = [
                            ds.coords[xs[j]].values,
                            ds.data_vars[key].values,
                            ds.data_vars[es1[j]].values,
                            ds.data_vars[es2[j]].values
                        ]
    
    # Create the figure using Plotly
    fig = go.Figure()
    
    # Plot all data on a single plot
    for i, (planet_name, (x, y, e1, e2)) in enumerate(return_dictionary.items()):
        fig.add_trace(go.Scatter(
            x=x, 
            y= y - np.mean(y) + i/20,  # stack each planet spectra above the previous
            # Plotly's way of doing two-sided error bars
            error_y=dict(type='data', array=e1, arrayminus=np.abs(e2)),
            mode='markers', 
            name=planet_name
        ))
    
    # Label graph based on event type and parameter
    if event == 'Transit':
        spectra = 'Transmission Spectra'
    elif event == 'Eclipse':
        spectra = 'Emission Spectra'
    if parameter == 'planet_radius':
        par = 'all planet radii'
    elif parameter == 'planet_mass':
        par = 'all planet mass'
    
    # Update figure layout
    fig.update_layout(
        title=dict(
            text=f'{spectra} for {par} [{lower_limit}, {upper_limit}]',
            font=dict(
                size=24
            )
        ),
        xaxis=dict(
            title=dict(
                text='Wavelength (μm)',
                font=dict(
                    size=20
                )
            )
        ),
        yaxis=dict(
            title=dict(
                text=f'{event} Depth',
                font=dict(
                    size=20
                )
            )
        ),
        height=1250,
        width=2000,
        legend=dict(
            x=0.01,
            y=0.99,
            xanchor='left',
            yanchor='top',
            bordercolor='Black',
            borderwidth=1
        )
    )
    
    # Show the figure
    fig.show()
    
        