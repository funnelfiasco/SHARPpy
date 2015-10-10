#!/usr/bin/env python
# encoding: utf-8
import os, sys
import numpy as np
import Tkinter as tk
import sharppy as sp
import sharppy.sharptab as tab
from sharppy.sharptab import profile, params, thermo

ver = sp.__version__


class MainApp(tk.Tk):

    def __init__(self, **kwargs):

        # Create a Profile Object if needed
        self.prof = kwargs.get('prof', None)

        # Check to see if profile object is a list, if not, make it one
        if not isinstance(self.prof, list):
            self.prof = [self.prof]
        kwargs['prof'] = self.prof

        # Update Keyword Arguments
        self.__dict__.update(kwargs)

        # Set up the main window
        tk.Tk.__init__(self)
        tk.Tk.configure(self, background=self.bgcolor)
        self.title("SHARPpy: Sounding/Hodograph Analysis and Research" \
        " Program in Python \t [Build: %s]" % ver)

        # Set Geometry Variables
        self.wid = self.skwtwid + self.hodowid
        self.hgt = max(self.skwthgt, self.hodohgt)
        self.screenx = self.winfo_screenwidth()
        self.screeny = self.winfo_screenheight()
        kwargs['screenx'] = self.screenx
        kwargs['screeny'] = self.screeny

        # Set Geometry
        x = self.screenx / 2 - self.wid / 2
        y = self.screeny / 2 - self.hgt / 2
        self.geometry('%sx%s+%s+%s' % (self.wid, self.hgt, x, y))

        # Place the SkewT
        self.skwt_frame = tk.Frame(width=self.skwtwid, height=self.skwthgt)
        self.skwt_frame.pack(side=tk.LEFT, anchor=tk.NW)
        self.skwt = SkewtApp(parent=self.skwt_frame, **kwargs)

        # Place the Hodograph
        self.hodo_frame = tk.Frame(width=self.hodowid, height=self.hodohgt)
        self.hodo_frame.pack(anchor=tk.NE)
        self.hodo = HodoApp(parent=self.hodo_frame, **kwargs)




class SkewtApp(tk.Frame):

    def __init__(self, **kwargs):
        # Setup Defaults
        self.skwtwid = 800
        self.skwthgt = 800
        self.parent = None
        self.prof = None
        self.initial = True
        self.font = ("Helvetica", 8)

        # Update All Keyword Arguments
        self.__dict__.update(kwargs)

        # Create Parent Frame
        self.createWidgets()


    def createWidgets(self):
        # Create Canvas
        self.canvas = tk.Canvas(self.parent, bg="black", width=self.skwtwid,
            height=self.skwthgt)
        self.canvas.bind("<Motion>", self.callback)
        self.canvas.pack(anchor=tk.S)

        # Create SkewT and Profile Classes
        self.sk = sp.SkewT(self.canvas, width=self.skwtwid,
            height=self.skwthgt)

        # Draw SkewT
        self.drawSkewt()


    def drawSkewt(self):
        self.sk.drawSkewT()

        self.sbpcls = []; self.mlpcls = []
        self.mupcls = []; self.epcls = []
        for i in range(len(self.prof)):
            # Create Surface Based Parcel
            sblplvals = params.DefineParcel(self.prof[i], 1)
            self.sbpcls.append(params.parcelx(-1, -1, sblplvals.pres,
                sblplvals.temp, sblplvals.dwpt, self.prof[i],
                lplvals=sblplvals))

            # Create Mixed Layer Parcel
            mllplvals = params.DefineParcel(self.prof[i], 4)
            self.mlpcls.append(params.parcelx(-1, -1, mllplvals.pres,
                mllplvals.temp, mllplvals.dwpt, self.prof[i],
                lplvals=mllplvals))

            # Create Most Unstable Parcel
            mulplvals = params.DefineParcel(self.prof[i], 3)
            self.mupcls.append(params.parcelx(-1, -1, mulplvals.pres,
                mulplvals.temp, mulplvals.dwpt, self.prof[i],
                lplvals=mulplvals))

            # Create Effective Parcel
            elplvals = params.DefineParcel(self.prof[i], 6)
            self.epcls.append(params.parcelx(-1, -1, elplvals.pres,
                elplvals.temp, elplvals.dwpt, self.prof[i],
                lplvals=elplvals))


            # Plot the Mixed Layer Parcels
            pcls = [self.mlpcls]
            cols = ['#FFFF00']
            for pcl, col in zip(pcls, cols):
                init = True
                for i in range(len(pcl)):
                    if pcl[i].bplus > 0:
                        if init: width = 3
                        else: width = 1
                        self.sk.drawVirtualParcelTrace(pcl[i], color=col,
                            width=width)
                        init = False

            # Plot the actual Parcels
            init = True
            for i in range(len(self.prof)):
                if init:
                    width = 4
                    plottxt = True
                else:
                    width = 1
                    plottxt = False
                self.sk.drawProfile(prof=self.prof[i], tracewidth=width,
                    plottxt=plottxt)
                init = False
            self.sk.drawBarbs(size=3.33, prof=self.prof[0])


    def callback(self, event, **kwargs):
        p = self.sk.pix2Pres(event.y)
        if p > self.prof[0].gSndg[self.prof[0].sfc][self.prof[0].pind]:
            p = self.prof[0].gSndg[self.prof[0].sfc][self.prof[0].pind]

        v0 = tab.interp.interp_from_pres(p, self.prof[0], 0)
        v1 = tab.interp.interp_from_pres(p, self.prof[0], 1)
        v2 = tab.interp.interp_from_pres(p, self.prof[0], 2)
        v3 = tab.interp.interp_from_pres(p, self.prof[0], 3)
        v4 = tab.interp.interp_from_pres(p, self.prof[0], 4)
        v5 = tab.interp.interp_from_pres(p, self.prof[0], 5)
        v4, v5 = tab.vector.comp2vec(v4, v5)
        vals = [v0, v1, v2, v3, v4, v5]
        if vals[1] < 0: vals[1] = -9999.0
        if vals[2] < -110: vals[2] = -9999.0
        if vals[3] < -110: vals[3] = -9999.0
        self.drawReadout(event, vals)


    def drawReadout(self, event, vals, **kwargs):
        p = self.sk.pix2Pres(event.y)
        if p > self.sk.pmax or \
            p < self.sk.pmin: return
        if self.initial:
            # Set the readout line
            self.readout = self.canvas.create_line(self.sk.tlx+1, event.y,
                self.sk.wid-1, event.y, fill='#333333', width=0.5)

            # Set the pressure readout
            self.preadoutbox = self.canvas.create_rectangle((self.sk.tlx+1,
                event.y-5, self.sk.tlx+51, event.y+5), fill=self.sk.framebg)
            self.preadout = self.canvas.create_text(self.sk.tlx+6, event.y,
                fill=self.sk.framefg, text='%-4i hPa' % (p), anchor=tk.W,
                font=self.font)

            # Set the height readout
            self.zreadoutbox = self.canvas.create_rectangle((self.sk.wid-51,
                event.y-5, self.sk.wid-1, event.y+5), fill=self.sk.framebg)
            self.zreadout = self.canvas.create_text(self.sk.wid-46, event.y,
                fill=self.sk.framefg, text='%-5i m' % (vals[1]), anchor=tk.W,
                font=self.font)

            # Set the temperature readout
            self.treadoutbox = self.canvas.create_rectangle((self.sk.wid-101,
                event.y-5, self.sk.wid-51, event.y+5), fill=self.sk.framebg)
            self.treadout = self.canvas.create_text(self.sk.wid-76, event.y,
                fill=self.sk.tcolor, text='%-3.1f C' % (vals[2]),
                anchor=tk.CENTER, font=self.font)

            # Set the dew point readout
            self.tdreadoutbox = self.canvas.create_rectangle((self.sk.tlx+51,
                event.y-5, self.sk.tlx+101, event.y+5), fill=self.sk.framebg)
            self.tdreadout = self.canvas.create_text(self.sk.tlx+76, event.y,
                fill=self.sk.tdcolor, text='%-3.1f C' % (vals[3]),
                anchor=tk.CENTER, font=self.font)

            # Set the wet-bulb temperature readout
            self.twreadoutbox = self.canvas.create_rectangle((self.sk.tlx+101,
                event.y-5, self.sk.tlx+151, event.y+5), fill=self.sk.framebg)
            self.twreadout = self.canvas.create_text(self.sk.tlx+126, event.y,
                fill=self.sk.twcolor, text='%-3.1f C' % (thermo.wetbulb(vals[0],
                vals[2], vals[3])), anchor=tk.CENTER, font=self.font)

            self.initial = False

        else:
            # Update the readout line
            x1, y1, x2, y2 = self.canvas.coords(self.readout)
            self.canvas.coords(self.readout, x1, event.y, x2, event.y)

            # Update the pressure readout
            xc, yc = self.canvas.coords(self.preadout)
            x1, y1, x2, y2 = self.canvas.coords(self.preadoutbox)
            self.canvas.coords(self.preadoutbox, (x1, event.y-5, x2,
                event.y+5))
            self.canvas.coords(self.preadout, xc, event.y)
            self.canvas.itemconfig(self.preadout, text='%-4i hPa' % (p))

            # Update the height readout
            xc, yc = self.canvas.coords(self.zreadout)
            x1, y1, x2, y2 = self.canvas.coords(self.zreadoutbox)
            self.canvas.coords(self.zreadoutbox, (x1, event.y-5, x2,
                event.y+5))
            self.canvas.coords(self.zreadout, xc, event.y)
            self.canvas.itemconfig(self.zreadout, text='%-5i m' % vals[1])

            # Update the wet-bulb temperature readout
            xc, yc = self.canvas.coords(self.treadout)
            x1, y1, x2, y2 = self.canvas.coords(self.treadoutbox)
            self.canvas.coords(self.treadoutbox, (x1, event.y-5, x2,
                event.y+5))
            self.canvas.coords(self.treadout, xc, event.y)
            self.canvas.itemconfig(self.treadout, text='%-3.1f C' % \
                (vals[2]))

            # Update the dew point readout
            xc, yc = self.canvas.coords(self.tdreadout)
            x1, y1, x2, y2 = self.canvas.coords(self.tdreadoutbox)
            self.canvas.coords(self.tdreadoutbox, (x1, event.y-5, x2,
                event.y+5))
            self.canvas.coords(self.tdreadout, xc, event.y)
            self.canvas.itemconfig(self.tdreadout, text='%-3.1f C' % \
                (vals[3]))

            # Update the wet-bulb temperature readout
            xc, yc = self.canvas.coords(self.twreadout)
            x1, y1, x2, y2 = self.canvas.coords(self.twreadoutbox)
            self.canvas.coords(self.twreadoutbox, (x1, event.y-5, x2,
                event.y+5))
            self.canvas.coords(self.twreadout, xc, event.y)
            self.canvas.itemconfig(self.twreadout, text='%-3.1f C' % \
                thermo.wetbulb(vals[0], vals[2], vals[3]))





class HodoApp(tk.Frame):

    def __init__(self, **kwargs):
        # Setup Defaults
        self.stn = 'OUN'
        self.hodowid = 600
        self.hodohgt = 600
        self.parent = None
        self.prof = None

        # Update All Keyword Arguments
        self.__dict__.update(kwargs)

        # Create Parent Frame
        self.createWidgets()


    def createWidgets(self):
        self.canvas = tk.Canvas(self.parent, bg="black", width=self.hodowid,
            height=self.hodohgt)
        self.canvas.pack(anchor=tk.S)

        self.ho = sp.Hodo(self.canvas, width=self.hodowid,
            height=self.hodohgt, prof=self.prof[0])
        self.drawHodo()


    def drawHodo(self, **kwargs):
        self.ho.drawHodo()
        for i in range(len(self.prof)):
            if i == 0: width = 4
            else: width = 1
            self.ho.drawProfile(self.prof[i], width)





def string_to_profile(sndgtxt, **kwargs):
    data = np.genfromtxt(cStringIO.StringIO(sndgtxt), names=True)

    sndg = profile.Profile(pres=data['PRES'], hght=data['HGHT'],
        temp=data['TMPC'], dwpt=data['DWPC'], wdir=data['DRCT'],
        wspd=data['SPED'], **kwargs)

    return sndg


def read_file(flines):
    sites = {}
    startlines = []
    for i, data in enumerate(flines):
        if 'STID' in data:
            startlines.append(i);
    for a,z in zip(startlines[:-1], startlines[1:]):
        data = "".join(flines[a+3:z-2])
        id = flines[a].split('=')[1].split()[0]
        if id not in sites:
            sites[id] = [data]
        else:
            sites[id] += [data]

    return sites


def parse_fname(fname):
    parts = fname.split('_')
    name = parts[1] + '_' + parts[2]
    return name, parts[3].split('.')[0]


if __name__ == '__main__':
    stn = sys.argv[1].upper()
    url = r'http://www.spc.noaa.gov/exper/soundings/LATEST/%s.txt' % (stn)
    profs = [profile.Profile(url=url)]
    app = MainApp(prof=profs, skwtwid=750, skwthgt=750, hodowid=400,
        hodohgt=400, ptxtwid=700, ptxthgt=200, plevs=None, mps=True, stn=stn,
        bgcolor='#000000', url=url)
    app.mainloop()



