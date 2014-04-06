#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <ialvarez.das.dcc@gmail.com> wrote this file. As long as you retain this
# notice you can do whatever you want with this stuff. If we meet some day,
# and you think this stuff is worth it, you can buy me a beer in return
# isma_comi
# ----------------------------------------------------------------------------
#

import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
import pyfits
# from astropy.nddata import NDData # no usado en esta version
# from pylab import *
# import scipy


def centroide(data):
    """
    funcion que calcula el centroide de data:
    se le da un arreglo data y retorna el indice, no entero, en donde se
    encuentra el centroide
    """
    Mr = 0
    M = 0
    for i in range(0, len(data)):
        Mr = Mr + i*data[i]
        M = M + data[i]
    if M == 0:
        return 0
    return Mr/M


def flujo(arreglo, ubicacion):
    """
    retorna el flujo en un punto dado un arreglo:
    como ubicacion no es entero hace una interpolacion lineal de los
    datos que rodean la ubicacion
    """
    a = arreglo[int(ubicacion)]
    b = arreglo[int(ubicacion+1)]
    return (((b-a)*(ubicacion-int(ubicacion)))+a)


def flujosky(arreglo, s1, s2):
    """
    retorna el flujo del arreglo del cielo:
    retorna la mediana de arreglo[s1:s2]
    """
    return np.median(arreglo[int(s1):int(s2)])


def apert_error(gain, apert, a):
    """
    funcion de error dada una apertura
    """
    # N_e_peak = gain * Fpeak
    # N_e_sky = gain * Fsky
    # SNR = N_e_peak /
    # sp.sqrt(
    #     N_e_peak +
    #     apert * (1.0 + float(apert)/apert_sky) * (N_e_sky + ron))
    SNR = np.sqrt(a)

    return SNR
