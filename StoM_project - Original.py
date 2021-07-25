import numpy as np
import scipy.stats as stats
import scipy.signal as signal
import matplotlib.pyplot as plt
import STOM_higgs_tools as higgs
import sys
import os
from tqdm import tqdm

#plt.rcParams['agg.path.chunksize'] = 20000

#Define equation for exponential PDF
def f(x, l):
    return np.exp((-1/l) * x)

def g(x, A, l, mu, sigma, signal_amp):
    return A * np.exp(-x/l) + higgs.signal_gaus(x, mu, sigma, signal_amp)

# Estimate the gradient of decay
def estimateGradient(bin_heights, bin_edges):
    l_estimate = (bin_edges[0] - bin_edges[len(bin_edges)-1]) / np.log(bin_heights[len(bin_heights)-1] / bin_heights[0])
    return l_estimate

# Estimate the normalisation factor
def estimateFactor(lambda_estimate, bin_edges, bin_heights):
    A = bin_heights[0] / np.exp((-1 / lambda_estimate) * bin_edges[0])
    return A

def estimateMean(bin_heights, bin_centres):
    peaks = signal.argrelmax(bin_heights)
    prominence = signal.peak_prominences(bin_heights, peaks[0])[0]
    i = np.where(prominence == np.max(prominence))
    mean = bin_centres[peaks[0][i]]
    return mean

def estimateSignalAmp(mean, bin_heights, bin_centres, A, l, sigma):
    no_signal = A * f(mean, l)
    peak_height = bin_heights[np.where(bin_centres == mean)]
    amp = (peak_height - no_signal) * sigma * np.sqrt(2 * np.pi)
    return amp

def estimateSigma(bin_heights, bin_centres, mean):
    index = np.where(bin_centres == mean)[0][0]
    heights = bin_heights[0:index+1]
    for i in range(len(heights) - 1):
        if heights[index - i] < heights[index - i - 1]:
            mark = bin_centres[index - i - 1]
            break
    sigma = (mean - mark) / 2
    return sigma

def getBackgroundChiSquare(A_estimates, l_estimates):
    # Store chi squared results in a matrix because it is a 2D iteration
    chi_squares = np.empty([len(A_estimates), len(l_estimates)])
    # Loop through both arrays of A and lambda and find chi squared for each iteration
    for i in range(len(A_estimates)):
        for j in range(len(l_estimates)):
            chi_squares[i][j] = higgs.get_B_chi(values, [104,155], 30, A_estimates[i], l_estimates[j])
    # Find location of smallest chi square
    result = np.where(chi_squares == np.amin(chi_squares))
    A = A_estimates[result[0]]
    l = l_estimates[result[1]]
    return np.amin(chi_squares), A, l

def getChiSquare(A_estimates, l_estimates, mean, standard_deviation, signal_amp):
    # Store chi squared results in a matrix because it is a 2D iteration
    chi_squares = np.empty([len(A_estimates), len(l_estimates)])
    # Loop through both arrays of A and lambda and find chi squared for each iteration
    for i in range(len(A_estimates)):
        for j in range(len(l_estimates)):
            chi_squares[i][j] = higgs.get_SB_chi(values, [104,155], 30, A_estimates[i], l_estimates[j], mean, standard_deviation, signal_amp)
    # Find location of smallest chi square
    result = np.where(chi_squares == np.amin(chi_squares))
    A = A_estimates[result[0]]
    l = l_estimates[result[1]]
    return np.amin(chi_squares), A, l

for count in tqdm(range(1)):
    # Generate data and plot histogram
    n_signals = 400
    values = higgs.generate_data(n_signals)
    bin_heights, bin_edges = np.histogram(values, range = [104.0, 155.0], bins = 30)
    bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
    err = (bin_edges[1] - bin_edges[0])/2
    plt.errorbar(bin_centres, bin_heights, yerr = None, xerr = err, fmt = "d")

    # Process data in step 2
    data = [x for x in values if (104 <=  x <= 155)]
    l_estimate = estimateGradient(bin_heights, bin_edges)
    A_estimate = estimateFactor(l_estimate, bin_edges, bin_heights)
    mean_estimate = estimateMean(bin_heights, bin_centres)
    sigma_estimate = estimateSigma(bin_heights, bin_centres, mean_estimate)
    signal_amp_estimate = estimateSignalAmp(mean_estimate, bin_heights, bin_centres, A_estimate, l_estimate, sigma_estimate)

    A_lower = A_estimate - 0.5*A_estimate
    A_upper = A_estimate + 0.5*A_estimate
    l_lower = l_estimate - 0.5*l_estimate
    l_upper = l_estimate + 0.5*l_estimate

    # Refine result from step 2 in step 3
    A_estimates = np.arange(A_lower,A_upper,0.1*A_estimate)
    l_estimates = np.arange(l_lower, l_upper, 0.1*l_estimate)
    

    min_chi_square, A, l = getChiSquare(A_estimates, l_estimates, mean_estimate, sigma_estimate, 700)
    print(min_chi_square)
    print(A)
    print(l)
    print(mean_estimate)
    print(sigma_estimate)
    print(signal_amp_estimate)
    #chi_squares.append(min_chi_square)
#min_chi_square_s = getChiSquare(A_estimates, l_estimates, signal = True)
#signal_chi_squares.append(min_chi_square_s)
#print(signal_chi_squares)
#np.savetxt("signal_amp_400.csv", signal_chi_squares, delimiter=',')


#min_chi, A, l = getSignalChiSquare(A_estimates, l_estimates, 1000)
'''
path = "C:/Viraj/University/First Year/Statistics of Measurement and Uncertainty/StoM Project/ChiSquareDistributions"
chi_distributions = []
for filename in os.listdir(path):
    distribution = np.loadtxt("ChiSquareDistributions/" + filename, dtype = float, delimiter=',')
    chi_distributions.append(distribution)
amp = [0, 500, 700, 1000]
for i in tqdm(range(len(chi_distributions))):
    if i == 0:
        plt.hist(chi_distributions[i], bins = 40, histtype = "stepfilled", density = True, lw = 3, fc = (0,0,1,0.5), label = str(amp[i]))
    else:
        plt.hist(chi_distributions[i], bins = 40, histtype = "step", density = True, label = str(amp[i]))
'''
no_signal = [g(x, A, l, mean_estimate, sigma_estimate, signal_amp_estimate) for x in np.sort(data)]
#signal = [A * f(x, l) for x in np.sort(data)]


plt.plot(np.sort(data), no_signal)
#plt.plot(np.sort(data), signal)
plt.xlabel("Rest mass")
plt.ylabel("Number")
'''
plt.hist(chi_squares, 10, histtype = "step")
#plt.plot(n, p_values, "o")
plt.xlabel(r"$\chi^2$")
plt.ylabel('Number')
'''
#plt.legend()
plt.savefig("second_fit_background_and_signal.png")
plt.show()

