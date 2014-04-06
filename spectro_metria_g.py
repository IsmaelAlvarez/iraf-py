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
from astrocalc_s import *
from astroplot_s import *
import spectro_metria
# import scipy #libreria interesante para continuar
# from astropy.nddata import NDData # no usado en esta version
# from pylab import * # interesante para continuar


class spect(spectro_metria.spect):

    def set_disp_axis(self, axis=None):
        if axis is None:
            if self.loadimg():
                self.disp_axis = setdispaxis(self.img)
                self.freeimg()
                d.limpiar()
        else:
            self.disp_axis = axis

    def plot_sky_pendent(self):
        plt.plot(self.pendiente)
        plt.ylabel('intensidad luminica')
        plt.xlabel('columnas del CCD')
        plt.title(
            'pendiente de la diferencia entre el cielo derecho' +
            ' y el izquerdo')
        plt.show()

    def plot_final_and_err(self):
        plt.plot(self.final)
        plt.errorbar(range(0, len(self.final)), self.final, yerr=self.errfinal)
        plt.ylabel('intensidad luminica')
        plt.xlabel('columnas del CCD')
        plt.title('espectro del cuerpo con su error')
        plt.show()

    def plot_final(self):
        plt.plot(self.final)
        plt.ylabel('intensidad luminica')
        plt.xlabel('columnas del CCD')
        plt.title('espectro del cuerpo')
        plt.show()

    def plot_sky(self):
        plt.plot(self.skyR, label='cielo derecho')
        plt.plot(self.skyL, label='cielo izquerdo')
        plt.ylabel('intensidad luminica')
        plt.xlabel('columnas del CCD')
        plt.title('espectro del cielo')
        plt.legend()
        plt.show()

    def plot_sky_and_err(self):
        plt.plot(self.skyR, label='cielo derecho')
        plt.plot(self.skyL, label='cielo izquerdo')
        plt.errorbar(
            range(0, len(self.skyR)),
            self.skyR, yerr=self.errskyR, label='error cielo derecho')
        plt.errorbar(
            range(0, len(self.skyL)),
            self.skyL, yerr=self.errskyL, label='error cielo izquerdo')
        plt.ylabel('intensidad luminica')
        plt.xlabel('columnas del CCD')
        plt.title('espectro del cielo con su error')
        plt.legend()
        plt.show()

    def plot_spec(self):
        if self.spectro1 is None:
            print '\033[93m' + 'Warning:' + '\033[0m' +
            ' spec_stract not done yet'
        plt.plot(self.spectro1)
        plt.ylabel('intensidad luminica')
        plt.xlabel('columnas del CCD')
        plt.title('espectro extraido')
        plt.show()

    def plot_spec_and_err(self):
        if self.spectro1 is None:
            print '\033[93m' + 'Warning:' + '\033[0m' +
            ' spec_stract not done yet'
        plt.plot(self.spectro1)
        plt.errorbar(
            range(0, len(self.spectro1)),
            self.spectro1, yerr=self.errspect, label='error espectro')
        plt.ylabel('intensidad luminica')
        plt.xlabel('columnas del CCD')
        plt.title('espectro extraido con su error')
        plt.show()

    def plot_traza_fit_and_img(self):
        if self.fit_traza is None:
            print '\033[93m' + 'Warning:' + '\033[0m' +
            ' fitt not done yet'
            return
        self.loadimg()
        imgplot = plt.imshow(self.img[0].data)
        plt.colorbar()
        plt.plot(self.traza, 'k-', label='traza')
        plt.plot(
            range(0, len(self.traza)), self.fit_traza, 'r-',
            label='fit')
        plt.ylabel('filas del CCD')
        plt.xlabel('columnas del CCD')
        plt.legend()
        plt.title('imagen con traza y fit obtenidos')
        plt.show()
        self.freeimg

    def plot_traza_and_fit(self):
        if self.fit_traza is None:
            print '\033[93m' + 'Warning:' + '\033[0m' +
            ' fitt not done yet'
            return
        plt.plot(self.traza, 'ko', label='traza')
        plt.plot(
            range(0, len(self.traza)), self.fit_traza, 'r-',
            label='fit')
        plt.ylabel('filas del CCD')
        plt.xlabel('columnas del CCD')
        plt.title('fit de la traza obtenida')
        plt.legend()
        plt.show()

    def plot_traza_and_img(self):
        if self.traza is None:
            self.get_traza()
        self.loadimg()
        imgplot = plt.imshow(self.img[0].data)
        plt.colorbar()
        plt.plot(self.traza, 'r-', label='traza')
        plt.ylabel('filas del CCD')
        plt.xlabel('columnas del CCD')
        plt.legend()
        plt.title('imagen con traza obtenida')
        plt.show()
        self.freeimg()

    def set_peak_aperture_and_sky(self):
        self.loadimg()
        setpeakerrsky(self.img, self.disp_axis)
        self.set_peak(sps.xpeak)
        self.set_aperture(sps.apert)
        self.set_skyL(sps.skyizq1, sps.skyizq2)
        self.set_skyR(sps.skyder1, sps.skyder2)
        sps.limpiar()
        self.freeimg()
