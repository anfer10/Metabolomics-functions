import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from pathlib import Path
import pandas as pd
import numpy.matlib
from scipy.signal import find_peaks

# STOCSY function 
def stocsy(folder_path):
    # Read csv 
    f2,data= csv_read_1H(folder_path) # Insert path 

    # Figure stack spectra 
    plt.figure
    plt.plot(f2, data)
    plt.xlabel("Chemical shift (ppm)")
    plt.ylabel("Intensity")
    plt.xlim(0, 10)
    plt.gca().invert_xaxis()
    plt.show()

    # Find peaks in mean spectrum 
    sum_spectrum= data.sum(axis=1) # Sum all data to obtain mean spectra

    # Figure mean spectra
    plt.figure()
    plt.plot(f2, sum_spectrum)
    plt.xlabel("Chemical shift (ppm)")
    plt.ylabel("Intensity")
    plt.xlim(0,10)
    plt.gca().invert_xaxis()
    plt.show()

    # Threshold 
    threshold= float(input("Enter threshold:"))  # Enter threshold to find peak 
    index_peaks, properties= find_peaks(sum_spectrum, height=threshold) # peaks list, peak index

    # Extract peak position and intensities
    peaks_shift= f2[index_peaks]  # Index position
    peaks_intensities= sum_spectrum[index_peaks] # Peaks intensities selected

    # Plot Peak picking
    plt.figure()
    plt.plot(f2, sum_spectrum)
    plt.plot(peaks_shift, peaks_intensities, "rx", label= "peaks")

    # Add text label with peak numbers
    for i, (x,y) in enumerate(zip(peaks_shift, peaks_intensities)):
        plt.text(x,y+0.5, str(i), color="red", fontsize=9, ha="center")

    plt.xlabel("Chemical shift (ppm)")
    plt.ylabel("Intensity")
    plt.xlim(0,10)
    plt.gca().invert_xaxis()
    plt.show()

    # Calculate covariance matrix and correlation 
    cov_matrix,corr_matrix,f2_list= cov_corr_matrix(data, index_peaks,f2) # Input: data matrix, index of peaks and f2

    # Index position 
    index_selected= int(input("Enter peak number to obtain STOCSY trace:")) 

    # STOCSY plot
    multi_color_plot(f2, cov_matrix[:, index_selected], corr_matrix[:, index_selected], cmap='jet', linewidth=2)
    plt.show()


    # Draft
    # If you want to do cov and correlation of the full matrix 
    # cov_matrix = np.cov(data_invert, rowvar=0)
    # corr_matrix = np.corrcoef(data_invert, rowvar=False)

# MultiColor plot function
def multi_color_plot(x, y, c, cmap='jet', linewidth=1): # Multicolorplot function to built STOCSY plot
    """
    Plots a line with color varying according to values in c using linear interpolation.

    Parameters:
    -----------
    x, y : 1D numpy arrays
        Coordinates of the line.
    c : 1D numpy array
        Color values for interpolation.
    cmap : str
        Colormap name (default: 'jet').
    linewidth : float
        Line width (default: 1).

    Returns:
    --------
    lc : LineCollection
        The LineCollection object for the plot.
    """

    x = np.asarray(x)
    y = np.asarray(y)
    c = np.asarray(c)

    # Ensure 1D row arrays
    if x.ndim > 1: x = x.flatten()
    if y.ndim > 1: y = y.flatten()
    if c.ndim > 1: c = c.flatten()

    # Create segments between points
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    # Create a LineCollection
    lc = LineCollection(segments, cmap=plt.get_cmap(cmap), norm=plt.Normalize(0, 1))
    lc.set_array(c)  # color values
    lc.set_linewidth(linewidth)

    # Add collection to current axis
    fig, ax = plt.gcf(), plt.gca()
    ax.add_collection(lc)
    ax.autoscale()
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(y.min(), y.max())
    
    plt.colorbar(lc, ax=ax, label='Correlation coefficient')  # optional colorbar
    plt.ylabel("Covariance")
    plt.xlabel("Chemical shift (ppm)")
    plt.xlim(0, 10)
    plt.gca().invert_xaxis()
    

    return lc

# Csv read 1H functionn 
def csv_read_1H(folder_path):
# Andres Charris
    all_spectra= []  # List empty to insert spectra

    for file in folder_path.glob("*.csv"):  # Find .csv files in the same directory 
        df= pd.read_csv(file, delimiter="\t")   # Read .csv. Output: pandas dataframe. It is a two-dimensional, size-mutable.
        intensity= df.iloc[:,1] # iloc to select rows or columns based  on their numerical position
        all_spectra.append(intensity)  # Insert each spectrum in the all_spectra list. 

    f2= df.iloc[:,0] # Select f2 
    spectra_df= pd.DataFrame(all_spectra).transpose() # matrix in nx len(f2)
    f2= pd.DataFrame(f2)  # f2 to dataframe
    f2= np.array(f2)  # f2 to Numpy array
    data= spectra_df.values  # spectra_df to Numpy array 
    return f2, data

# Product of column to calculate covariance
def cov_corr_matrix(data, index_peaks,f2):
    # Andres Charris
    """ Input:
            data: data size vxn where n is the sample number and v intensities
            index_peaks: index of selected peaks in f2
            f2: chemical shift vector
    """
    data_invert= data.transpose()  # data is a row matrix vxn, where n is the sample number and v is the intensity
    variables_selected= data_invert[:,index_peaks] # Variable selected. These variables are peaks in the mean spectrum
    n= len(data_invert[:,1])
    mean_data= np.mean(data_invert, axis=0) # Mean of data
    repmat_mean= numpy.matlib.repmat(mean_data, n,1) # Repmat values of mean to subtract
    mean_centered_cada= data_invert-repmat_mean # Mean center data
    variables_mean_centered= variables_selected-np.mean(variables_selected, axis=0) # Mean center variable selected
    mean_centered_invert= mean_centered_cada.transpose()   # Mean centered data invert 
    cov_matrix= (1/(n-1))*np.dot(mean_centered_invert, variables_mean_centered) # Covariance XXT
    std_data= np.std(data_invert,axis=0, ddof=1) # Standar deviation of all variable 
    std_data_peaks= np.std(variables_selected,axis=0, ddof=1) #  Standar deviation of variable selected
    iner_sdt= np.outer(std_data, std_data_peaks) # Pint by point product stdx*stdy. Product of std from all data and variables selected
    corr_matrix= np.divide(cov_matrix,iner_sdt) # Correlation amtrix calculation
    f2_list= f2.flatten().tolist()   # f2 as a list of values
    return cov_matrix, corr_matrix,f2_list    # Return covariance correlation and matrix

