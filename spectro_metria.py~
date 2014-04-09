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
# import scipy # libreria interesante para continuar robusteciendo
# la aplicacion
# from astropy.nddata import NDData # no usado en esta version
# from pylab import * # libreria interesante para continuar
# robusteciendo la aplicacion


class spect:
    """
    clase principal para el manejo de espectros obtenidos por un fits
    """
    archive_path = None  # direccion del archivo de origen del spectro
    disp_axis = 2  # horisontal es 2 y vertical es 1
    exp_time = None  # usado para almacenar el tiempo de expocicion
    gain = None
    traza = None  # arreglo de ubicacion del peak en el CCD
    fit_traza = None  # fiteo de la traza
    peak = None
    apertura = None  # numero fijo
    deltaR = None  # num fijo de error derecho del peak al sky derecho
    deltaL = None  # num fijo de error izquerdo del peak al sky izquerdo
    skyDR = None
    skyDL = None
    spectro1 = None  # arreglo de spectro del cuerpo sin ningun retoque
    errspect = None
    skyR = None  # arreglo de spectro de sky derecho
    errskyR = None
    skyL = None  # arreglo de spectro de sky izquerdo
    errskyL = None
    final = None
    errfinal = None
    pendiente = None

    def __init__(self, path=None):
        """
        init apartir de el archivo origen
        """
        self.archive_path = path

    '''pendiente:
    +atributo de tiempo en que fue tomado con:
        -from astropy.time import Time
        -ver http://docs.astropy.org/en/stable/time/index.html
    +atributo de observatorio en que se tomo
    +atributo de telescopio usado
    +coordenadas del objeto con:
    -from astropy.coordinates import ICRS, Galactic
        -ver http://docs.astropy.org/en/stable/coordinates/index.html
    +usar NDData para el cielo y peak etc
    '''

    def get_exp_time(self):
        """
        obtencion del tiempo de expocicion de la imagen
        """
        if self.exp_time is not None:
            return self.exp_time
        else:
            if self.archive_path is None:
                print 'archive path not especified'
            else:
                self.loadimg()
                self.exp_time = self.img[0].header['EXPTIME']
                del self.img
                return self.exp_time

    def get_gain(self):
        """
        geter de la ganancia
        """
        if self.gain is not None:
            return self.gain

    def set_gain(self, gain):
        """
        seter de la ganancia
        """
        self.gain = gain

    def set_path(self, path=None):
        """
        seter del archivo de origen
        """
        self.archive_path = path

    def get_path(self):
        """
        geter del archivo de origen
        """
        return self.archive_path

    def display_path(self):
        """
        display del path de origen, a futuro remplasar print por un ret
        """
        print "path: %s" % self.archive_path

    def loadimg(self):
        """
        carga en ram la imagen de origen
        """
        if self.archive_path is None:
            print "path not especified"
            return False
        else:
            self.img = pyfits.open(self.archive_path)
            return True

    def freeimg(self):
        """
        livera de la ram la imegen de origen
        """
        del self.img

    def set_peak(self, peak):
        """
        seter del peak
        """
        self.peak = peak

    def set_aperture(self, aperture):
        """
        seter de la apertura
        """
        self.apertura = aperture

    def set_skyL(self, s1, s2):
        """
        seter del Lsky
        """
        self.deltaL = self.peak - self.apertura - s2
        self.skyDL = s2 - s1

    def set_skyR(self, s1, s2):
        """
        seter del Rsky
        """
        self.deltaR = s2 - self.peak + self.apertura
        self.skyDR = s2 - s1

    def get_traza(self):
        """
        funcion para obtener la traza suponiendo que ya tenemos el peak
        y la apertura
        """
        if (self.peak is None or self.apertura is None):
            print '\033[93m' + 'Warning:' + '\033[0m' +
            ' peak or aperture not set'
            # colores: '\033[xxm'
            return
        self.loadimg()
        data = self.img[0].data
        if self.disp_axis == 2:
            data = zip(*data)
        self.traza = np.zeros(len(data))
        auxpeak = self.peak

        for i in range((int(len(data) / 2) + 1), len(data)):
            auxpeak = centroide(
                data[i][int(auxpeak - self.apertura):
                        int(auxpeak + self.apertura + 1)]) +
            int(auxpeak - self.apertura)
            self.traza[i] = auxpeak

        auxpeak = self.peak

        for i in range(int(len(data)/2), -1, -1):
            auxpeak = centroide(
                data[i][int(auxpeak - self.apertura):
                        int(auxpeak + self.apertura + 1)]) +
            int(auxpeak - self.apertura)
            self.traza[i] = auxpeak

        del data
        del auxpeak
        self.freeimg()
        return self.traza

    def do_fit_traza(self, s, grado=None):
        """
        funcion para generar un fiteo suponiendo que ya tenemos la traza
        grado es solo usado para el caso de que queramos fitear un
        polinomio
        s=g para gauseana
        s=l para lineal
        s=p para polinomial
        """
        if s == 't':
            f_init = models.Trapezoid1D(
                amplitude=1., x_0=0., width=1., slope=0.5)
        if s == 'g':
            f_init = models.Gaussian1D(amplitude=1., mean=0, stddev=1.)
        if s == 'l':
            f_init = models.Linear1D(slope=1., intercept=1.)
        if s == 'p':
            f_init = models.Polynomial1D(grado)
        f = fitting.NonLinearLSQFitter()
        g = f(f_init, range(0, len(self.traza)), self.traza)
        self.fit_traza = g(range(0, len(self.traza)))

    def get_fit_traza(self):
        """
        geter de la traza
        """
        if self.fit_traza is None:
            print '\033[93m'+'Warning:'+'\033[0m'+' fitt not done yet'
            return
        return self.fit_traza

    def do_spec_stract(self):
        """
        funcion que realiza la extraccion del espectro ya teniendo un
        fiteo
        """
        if self.fit_traza is None:
            print '\033[93m'+'Warning:'+'\033[0m'+' fitt not done yet'
            return
        self.loadimg()
        data = self.img[0].data
        if self.disp_axis == 2:
            data = zip(*data)
        self.errspect = np.zeros(len(data))
        self.spectro1 = np.zeros(len(data))

        for i in range(0, len(data)):
            self.spectro1[i] = flujo(data[i], self.fit_traza[i])
            self.errspect[i] = apert_error(
                self.gain, self.apertura, self.spectro1[i])

        self.freeimg
        self.spectro1 = map(lambda x: x/self.exp_time, self.spectro1)
        self.errspect = map(lambda x: x/self.exp_time, self.errspect)

    def get_spec(self):
        """
        geter del espectro
        """
        if self.spectro1 is None:
            print '\033[93m'+'Warning:'+'\033[0m'+' spec_stract not done yet'
            return
        return self.spectro1

    def get_spec_error(self):
        """
        geter del error del espectro
        """
        return self.errspect

    def do_sky_stract(self):
        """
        funcion que realiza la extraccion del espectro del cielo, tanto
        R como L
        """
        self.loadimg()
        data = self.img[0].data
        if self.disp_axis == 2:
            data = zip(*data)
        self.skyR = np.zeros(len(data))
        self.skyL = np.zeros(len(data))
        self.errskyR = np.zeros(len(data))
        self.errskyL = np.zeros(len(data))

        for i in range(0, len(data)):
            self.skyL[i] = flujosky(
                data[i], self.fit_traza[i] - self.apertura -
                self.deltaL - self.skyDL, self.fit_traza[i] -
                self.apertura - self.deltaL)
            self.skyR[i] = flujosky(
                data[i], self.fit_traza[i] + self.apertura +
                self.deltaR, self.fit_traza[i] + self.apertura +
                self.deltaR + self.skyDR)
            self.errskyL[i] = apert_error(
                self.gain, self.apertura, self.skyL[i])
            self.errskyR[i] = apert_error(
                self.gain, self.apertura, self.skyR[i])

        del self.img
        self.skyL = map(lambda x: x/self.exp_time, self.skyL)
        self.skyR = map(lambda x: x/self.exp_time, self.skyR)
        self.errskyL = map(lambda x: x/self.exp_time, self.errskyL)
        self.errskyR = map(lambda x: x/self.exp_time, self.errskyR)

    def get_sky_L(self):
        """
        geter del espectro del skyL
        """
        return self.skyL

    def get_sky_R(self):
        """
        geter del espectro del skyR
        """
        return self.skyR

    def get_sky_L_error(self):
        """
        geter del error del espectro del skyL
        """
        return self.errskyL

    def get_sky_R_error(self):
        """
        geter del error del espectro del skyR
        """
        return self.errskyR

    def do_final_spect(self):
        """
        restamos el cielo del espectro
        """
        self.final = np.zeros(len(self.spectro1))
        self.errfinal = np.zeros(len(self.spectro1))

        for i in range(0, len(self.spectro1)):
            self.final[i] = self.spectro1[i] -
            (self.skyR[i] - self.skyL[i]) *
            ((self.apertura + self.deltaL + (self.skyDL / 2)) /
                (self.deltaR + (self.skyDR / 2) - self.deltaL -
                    (self.skyDL / 2))) +
            self.skyL[i]
            self.errfinal[i] = np.sqrt(
                np.power(self.errspect[i], 2) +
                np.sqrt(
                    self.errskyR[i] * self.errskyR[i] +
                    self.errskyL[i] * self.errskyL[i]) *
                np.power(
                    ((self.apertura + self.deltaL + (self.skyDL / 2)) /
                        (self.deltaR + (self.skyDR / 2) - self.deltaL -
                            (self.skyDL / 2))), 2) +
                np.power(self.errskyL[i], 2))

    def get_final_spect(self):
        """
        geter del espectro final
        """
        return self.final

    def get_final_spect_error(self):
        """
        geter del error del espectro final
        """
        return self.errfinal

    def get_sky_pendent(self):
        """
        geter de la pendiente entre los skyL y skyR
        """
        self.pendiente = np.zeros(len(self.skyR))

        for i in range(0, len(self.skyR)):
            self.pendiente[i] = (self.skyL[i] - self.skyR[i]) /
            (self.deltaL + (self.skyDL / 2) - self.deltaR +
                (self.skyDR / 2))

        return self.pendiente

    def recalcularpeak(self, spe=None):
        """
        funcion para recalcular el peak
        si spe=none recalcula usando la funcion centroide en el area
        local de la apertura
        spe puede ser un objeto spect y se puede recalcular usando el
        peak y la apertura de este objeto
        """
        if spe is None:
            self.loadimg()
            data = self.img[0].data
            self.peak = centroide(
                data[(int(len(data) / 2))]
                [int(self.peak - self.apertura):
                    int(self.peak + self.apertura + 1)]) +
            self.peak-self.apertura
            self.freeimg()
        else:
            self.loadimg()
            data = self.img[0].data
            self.peak = spe.peak
            self.apertura = spe.apertura
            self.deltaR = spe.deltaR
            self.deltaL = spe.deltaL
            self.skyDR = spe.skyDR
            self.skyDL = spe.skyDR
            self.freeimg()
            self.recalcularpeak()

    def set_peak_aperture_skyL_skyR(
            self, peak, aperture, skyl1, skyl2, skyr1, skyr2):
        """
        funcion para setear rapidamente peak apertura y los sky
        """
        self.set_peak(peak)
        self.set_aperture(aperture)
        self.set_skyL(skyl1, skyl2)
        self.set_skyR(skyr1, skyr2)

    def set_disp_axis(self, axis):
        """
        seter del eje de dispecion
        """
        self.disp_axis = axis
