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


class dispes:
    'usado para definir el eje de dispercion'
    x1 = 0
    y1 = 0
    x2 = 0
    y2 = 0
    confirmacion = 0

    def __init__(self, x=0, y=0):
        self.x1 = x
        self.y1 = y

    def setx1(self, x=0, y=0):
        self.x1 = x
        self.y1 = y

    def setx2(self, x=0, y=0):
        self.x2 = x
        self.y2 = y

    def disper(self):
        if self.x1 - self.x2 > self.y1 - self.y2:
            return 1  # dispersion en el primer indice, x
        else:
            return 2  # dispersion en el segundo indice, y

    def limpiar(self):
        self.setx1()
        self.setx2()
        self.confirmacion = 0

d = dispes(0, 0)


class skypeaksky:
    ypeak = None
    ysky = None
    xpeak = None
    apert = None
    skyder1 = None
    skyder2 = None
    skyizq1 = None
    skyizq2 = None
    confir = 0

    def __init__(self, yp=None):
        ypeak = yp

    def limpiar(self):
        self.ypeak = None
        self.ysky = None
        self.xpeak = None
        self.apert = None
        self.skyder1 = None
        self.skyder2 = None
        self.skyizq1 = None
        self.skyizq2 = None
        self.confir = 0

sps = skypeaksky()


def on_keysky1(event):

    if event.key == 'd':  # restet
        plt.close()

    if event.key == 'b':  # confirm
        sps.confir = 1
        plt.close()

    if event.key == 'enter':
        plt.plot(
            [sps.skyizq1, sps.skyizq2], [sps.ysky, sps.ysky], 'r-',
            linewidth=3.0)
        plt.title('to continue press \'b\' to reset \'d\'')
        plt.show()

    if event.key == 'q':
        sps.ysky = event.ydata
        sps.skyizq1 = event.xdata
        plt.plot(sps.skyizq1, sps.ysky, 'ro')
        plt.show()

    if event.key == 'e':
        sps.skyizq2 = event.xdata
        plt.plot(sps.skyizq2, sps.ysky, 'ro')
        plt.show()


def on_keysky2(event):

    if event.key == 'd':  # restet
        plt.close()

    if event.key == 'b':  # confirm
        sps.confir = 1
        plt.close()

    if event.key == 'enter':
        plt.plot(
            [sps.skyder1, sps.skyder2], [sps.ysky, sps.ysky], 'r-',
            linewidth=3.0)
        plt.title('to continue press \'b\' to reset \'d\'')
        plt.show()

    if event.key == 'q':
        sps.skyder1 = event.xdata
        plt.plot(sps.skyder1, sps.ysky, 'ro')
        plt.show()

    if event.key == 'e':
        sps.skyder2 = event.xdata
        plt.plot(sps.skyder2, sps.ysky, 'ro')
        plt.show()


def on_keypeak(event):

    if event.key == 'enter':
        plt.plot(
            [sps.xpeak-sps.apert, sps.xpeak+sps.apert],
            [sps.ypeak, sps.ypeak], 'r-', linewidth=3.0)
        plt.title('to continue press \'b\' to reset \'d\'')
        plt.show()

    if event.key == 'd':
        # restet
        plt.close()

    if event.key == 'b':
        # confirm
        sps.confir = 1
        plt.close()

    if event.key == 'w':
        # centro
        sps.ypeak = event.ydata
        sps.xpeak = event.xdata
        plt.plot(sps.xpeak, sps.ypeak, 'ro')
        plt.show()

    if event.key == 'q':
        # error izquerdo
        sps.apert = sps.xpeak - event.xdata
        plt.plot(sps.xpeak - sps.apert, sps.ypeak, 'ro')
        plt.plot(sps.xpeak + sps.apert, sps.ypeak, 'ro')
        plt.show()

    if event.key == 'e':
        # error derecho
        sps.apert = event.xdata - sps.xpeak
        plt.plot(sps.xpeak - sps.apert, sps.ypeak, 'ro')
        plt.plot(sps.xpeak + sps.apert, sps.ypeak, 'ro')
        plt.show()


def on_keydisp(event):
    if event.key == 'enter':
        plt.plot([d.x1, d.x2], [d.y1, d.y2], 'r-', linewidth=3.0)
        plt.title('to continue press \'b\' to reset \'d\'')
        plt.show()

    if event.key == 'd':
        plt.close()

    if event.key == 'b':
        d.confirmacion = 1
        plt.close()

    if event.key == 'q':
        d.setx1(event.xdata, event.ydata)
        plt.plot(d.x1, d.y1, 'ro')
        plt.show()

    if event.key == 'e':
        d.setx2(event.xdata, event.ydata)
        plt.plot(d.x2, d.y2, 'ro')
        plt.show()


def setdispaxis(img):
    # en caso de no saber cual es el eje de dispersion

    data = img[0].data
    while d.confirmacion == 0:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        imgplot = plt.imshow(data)
        plt.colorbar()
        cid = fig.canvas.mpl_connect('key_press_event', on_keydisp)
        plt.title('set disp axis using \'q\' and \'e\' then \'enter\'')
        plt.ylabel('filas del CCD')
        plt.xlabel('columnas del CCD')
        plt.show()
    del data
    return d.disper()


def setpeakerrsky(img, eje_de_disp):  # remplaza apall de IRAF
    data = img[0].data
    if eje_de_disp == 2:
        data = zip(*data)  # trasponemos los datos
    # en otro caso es disp==1 y no nesecitamos trasponer los datos
    aux = np.zeros((len(data)))
    a = len(data)

    for i in range(int(a/2), int(a/2) + 10):
        # suma de las columnas 0 a la 20
        aux = aux + data[i]

    while sps.confir == 0:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.plot(aux)
        s = 'suma de las columnas ' + str(int(a/2)) + ' a ' +
        str(int(a / 2) + 10) + '\n presione \'w\' para setear el' +
        ' peak y \'q\' y \'e\' para la apertura luego \'enter\''
        ax.set_title(s)
        cid = fig.canvas.mpl_connect('key_press_event', on_keypeak)
        plt.ylabel('intensidad luminica')
        plt.xlabel('columnas del CCD')
        plt.show()

    sps.confir = 0

    while sps.confir == 0:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.plot(aux)
        plt.plot(sps.xpeak, sps.ypeak, 'r+')
        plt.plot(sps.xpeak - sps.apert, sps.ypeak, 'r+')
        plt.plot(sps.xpeak + sps.apert, sps.ypeak, 'r+')
        plt.plot(
            [sps.xpeak - sps.apert, sps.xpeak + sps.apert],
            [sps.ypeak, sps.ypeak], 'r-', linewidth=1.0)
        s = 'seleccion de sky izquerdo\n suma de las columnas ' +
        str(int(a / 2)) + ' a ' + str(int(a / 2) + 10) +
        '\n presione \'q\' y \'e\' para setear sky luego \'enter\''
        ax.set_title(s)
        cid = fig.canvas.mpl_connect('key_press_event', on_keysky1)
        plt.ylabel('intensidad luminica')
        plt.xlabel('columnas del CCD')
        plt.show()

    sps.confir = 0

    while sps.confir == 0:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.plot(aux)
        plt.plot(sps.xpeak, sps.ypeak, 'r+')
        plt.plot(sps.xpeak - sps.apert, sps.ypeak, 'r+')
        plt.plot(sps.xpeak + sps.apert, sps.ypeak, 'r+')
        plt.plot(
            [sps.xpeak - sps.apert, sps.xpeak + sps.apert],
            [sps.ypeak, sps.ypeak], 'r-', linewidth=1.0)
        plt.plot(sps.skyizq1, sps.ysky, 'r+')
        plt.plot(sps.skyizq2, sps.ysky, 'r+')
        plt.plot(
            [sps.skyizq1, sps.skyizq2], [sps.ysky, sps.ysky],
            'r-', linewidth=1.0)
        s = 'seleccion de sky derecho\n suma de las columnas ' +
        str(int(a / 2)) + ' a ' + str(int(a / 2) + 10) +
        '\n presione \'q\' y \'e\' para setear sky luego \'enter\''
        ax.set_title(s)
        cid = fig.canvas.mpl_connect('key_press_event', on_keysky2)
        plt.ylabel('intensidad luminica')
        plt.xlabel('columnas del CCD')
        plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.plot(aux)
    plt.plot(sps.xpeak, sps.ypeak, 'r+')
    plt.plot(sps.xpeak - sps.apert, sps.ypeak, 'r+')
    plt.plot(sps.xpeak + sps.apert, sps.ypeak, 'r+')
    plt.plot(
        [sps.xpeak - sps.apert, sps.xpeak + sps.apert],
        [sps.ypeak, sps.ypeak], 'r-', linewidth=1.0)
    plt.plot(sps.skyizq1, sps.ysky, 'r+')
    plt.plot(sps.skyizq2, sps.ysky, 'r+')
    plt.plot(
        [sps.skyizq1, sps.skyizq2], [sps.ysky, sps.ysky], 'r-',
        linewidth=1.0)
    plt.plot(sps.skyder1, sps.ysky, 'r+')
    plt.plot(sps.skyder2, sps.ysky, 'r+')
    plt.plot(
        [sps.skyder1, sps.skyder2], [sps.ysky, sps.ysky], 'r-',
        linewidth=1.0)
    s = 'suma de las columnas ' + str(int(a / 2)) + ' a ' +
    str(int(a / 2) + 10) + '\n presione \'b\' para continuar'
    ax.set_title(s)
    cid = fig.canvas.mpl_connect('key_press_event', on_keysky2)
    plt.ylabel('intensidad luminica')
    plt.xlabel('columnas del CCD')
    plt.show()
    del data
