from distutils.cmd import Command
from logging import root
from operator import le
from pickle import FRAME
from pydoc import describe 
from tkinter import *
from tkinter import ttk
from tkinter import messagebox
from traceback import FrameSummary
from turtle import color, fd

import math 
import numpy as np
import matplotlib.image as mpimg
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt

import lab2_1 as fibra_dibujo

class principal_Atenuaciones():
    def __init__(self,root):
        self.root = root
        self.lambdafo = StringVar()
        self.parv = StringVar()
        self.num_uniones = StringVar()
        self.atxval = []
        self.tipodeatx = []
        self.fibras = []
        self.ftotales = []
        self.fuentepara = []
        self.seleccion = 1
        self.selmultlat = IntVar()
        self.tamaFO = 0
        self.tamfo1 = StringVar()
        self.atxconect = StringVar()
        self.atxco = ""
        self.atxunion = StringVar()
        self.atxuni = ""
        self.atenuaciones = 0
        self.perfil = ""
        self.vali = 0
        self.parav = 0
        self.tipo = ""
        self.nucleoFO = 0
        self.val = 1
        self.uniones = 0
        self.unionnum = 1
        self.Pein = StringVar()
        self.Pout = StringVar()
        self.snrmin = StringVar()
        self.pnoise = []
        self.vecesfinal = 0
        self.vo = StringVar()
        self.rl = StringVar()
        self.res = StringVar()
        self.n = StringVar()
        self.ganancia = StringVar()
        pass
    
    #Formulas

    def numero_modos(self,perfil,parav):
        if perfil == "Lineal":
            modo = (parav**2)/6
        elif perfil == "Parabolico":
            modo = (parav**2)/4
        elif perfil == "Escalonado":
            modo = (parav**2)/2
        modo = round(modo,3)
        return modo
    
    def angulo_critico(self,indice_nucleo,indice_revestimiento):
        angc = math.asin(indice_revestimiento/indice_nucleo)
        angc = math.degrees(angc)
        return angc
    
    def angulo_apertura(self,indice_nucleo,angc):
        anga = math.degrees(math.asin(indice_nucleo*math.cos(math.radians(angc))))
        return anga

    def apertura_numerica(self,anga,indice_nucleo,indice_revestimiento):
        an_f1 = math.sin(math.radians(anga))
        an_f2 = math.sqrt((indice_nucleo**2)-(indice_revestimiento**2))
        na = an_f1
        return na

    def diff(self,na,indice_nucleo):
        dif = (na**2)/(2*(indice_nucleo**2))
        return dif

    def parametro_v(self,nucleo,indices,vlambda,perfil,tam):
        vlambda = float(vlambda)
        angc = self.angulo_critico(indices[0],indices[1])
        anga = self.angulo_apertura(indices[0],angc)
        na = self.apertura_numerica(anga,indices[0],indices[1])           
        
        para_v = round((2*math.pi*nucleo*na)/(vlambda*math.pow(10,-9)),3)
        self.seleccion = 2
        self.parv.set(para_v)
        self.vali = 1
        self.fibras.append(angc)
        self.fibras.append(anga)
        self.fibras.append(na)
        self.parav = para_v
        self.nucleoFO = nucleo
        if self.parav<=2.4:
            self.tipo = "Monomodo"
        else:
            self.tipo = "Multimodo"
        
        self.perfil = perfil
        self.fibraoptica_parv(indices,vlambda,nucleo,perfil,tam)

    def nucleo(self,par_v,indices,vlambda,num_uniones):
        if par_v>0:
            i = []
            i.append(float(indices[1:5]))
            i.append(float(indices[6:10])) 
            indices = i
            angc = self.angulo_critico(indices[0],indices[1])
            anga = self.angulo_apertura(indices[0],angc)
            na = self.apertura_numerica(anga,indices[0],indices[1])
            self.nucleoFO = (par_v*vlambda)/(na*2*math.pi)
            self.parav = par_v
            self.fibras.append(angc)
            self.fibras.append(anga)
            self.fibras.append(na)
            if par_v<2.4:
                self.tipo = "Monomodo"
            else:
                self.tipo = "Multimodo"
            self.guardar_fo(self.tipo,indices,float(self.parv.get()),self.num_uniones.get())
        else:
            messagebox.showwarning("Fibra Optica","Parametro v mayor a 0")
            self.parv.set("")
            self.num_uniones.set("")
    
    def perdidasdistancia(self,distancia,Ps,Pe):
        alpha = (-10/distancia)*math.log10(Ps/Pe)
        perdida = distancia*alpha
        return perdida

    #desalineaciones---------------------------------------------------------------------------------------------------------
    def lateral_mono(self,y,v,nucleo):
        w = nucleo*(0.65+(1.62*math.pow((v),(-3/2)))+(2.8*math.pow(v,-6)))/math.sqrt(2)
        print(w)
        loss_latmono = 2.17*math.pow((y/w),2)
        return loss_latmono
    
    def angular_mono(self,omg,v,nucleo,indice_nucleo,na):
        omg = math.radians(omg)
        w = nucleo*(0.65+(1.62*math.pow((v),(-3/2)))+(2.8*math.pow(v,-6)))/math.sqrt(2)
        print(w)
        loss_angmono = 2.17*math.pow(((omg*indice_nucleo*v*w)/(nucleo*na)),2)
        return loss_angmono
    
    def lateralgap_multi(self,y,indice_nucleo,nucleo):
        n1 = (16*(indice_nucleo**2))/(((1+indice_nucleo)**4)*math.pi)
        n2 = 2*math.acos(y/(2*nucleo))
        n3 = -(y/nucleo)*math.sqrt(1-math.pow((y/(2*nucleo)),2))
        n2 = n2+n3
        loss_latgapmulti = n1*n2
        return loss_latgapmulti
    
    def lateralsingap_multi(self,y,nucleo):
        n1 = 1/math.pi
        n2 = 2*math.acos(y/(2*nucleo))-(y/nucleo)*math.sqrt(1-math.pow((y/(2*nucleo)),2))
        loss_latsingapmulti = n1*n2
        return loss_latsingapmulti
    
    def indicegradual(self,indicegr,y,nucleo):
        if indicegr == "Lineal":
            lt = (3/math.pi)*(y/nucleo)
        elif indicegr == "Parabolico":
            lt = (8/(3*math.pi))*(y/nucleo)
        elif indicegr == "Escalonado":
            lt = (2/math.pi)*(y/nucleo)
        return lt

    def angular_multi(self,nucleo,ang,na):
        ang = math.radians(ang)
        n1 = (16*(nucleo**2))/(((1+nucleo)**4)*math.pi)
        n2 = 1-((1*ang)/(math.pi*na))
        loss_angmulti = n1*n2
        return loss_angmulti
    
    def longitudinal_multi(self,n0,z,na,nucleo):
        loss_longmulti = -10*math.log10(1-((z*na)/(4*nucleo*n0)))
        return loss_longmulti
    
    def pasaradb(self,veces):
        db = -10*math.log10(veces)
        return db

    def pasardbm(self,p):
        dbm = 10*math.log10((p/math.pow(10,-3)))
        return dbm
    
    def pasarapot(self,dbm):
        pot = (10**(dbm/10))*(math.pow(10,-3))
        return pot
    
    def calculomono_lat(self,y):
        try:
            y = float(y)
            lim = self.fibras[3]*math.pow(10,6)
            if (y>0 and y<=lim):
                self.frmmono_lat.destroy()
                y = y*math.pow(10,-6)
                nlat = self.lateral_mono(y,self.fibras[2],self.nucleoFO)
                self.atxval.append(nlat)
                self.val = 0
                
                print(self.num_uniones.get())
                x = int(self.num_uniones.get())-1
                self.num_uniones.set(str(x))

                self.guardar_fo(self.fibras[0],self.fibras[1],self.fibras[2],self.num_uniones.get(),self.atxco,self.atxuni,self.atxco,self.atxuni)
            else:
                messagebox.showwarning("Fibra Optica","y no valido")
        except ValueError:
            messagebox.showerror("Fibra Optica","y no valido")

    def calculomono_ang(self,ang):
        try:
            ang = float(ang)
            if (ang>0 and ang<=10):
                self.frmmono_ang.destroy()
                nang = self.angular_mono(ang,self.fibras[2],self.nucleoFO,self.fibras[1][0],self.fibras[6])
                self.atxval.append(nang)
                self.val = 0
                
                print(self.num_uniones.get())
                x = int(self.num_uniones.get())-1
                self.num_uniones.set(str(x))

                self.guardar_fo(self.fibras[0],self.fibras[1],self.fibras[2],self.num_uniones.get(),self.atxco,self.atxuni)
            else:
                messagebox.showwarning("Fibra Optica","𝜃 no valido")
        except ValueError:
            messagebox.showerror("Fibra Optica","𝜃 no valido")

    def calculomult_latsingap(self,y,indicegr):
        try:
            y = float(y)
            lim = self.fibras[3]*math.pow(10,6)
            if (y>0 and y<=lim):
                self.frmmut_lat.destroy()
                y = y*math.pow(10,-6)
                if self.selmultlat.get() == 1:
                    lt = self.indicegradual(indicegr,y,self.nucleoFO)
                    nlat = 1-lt
                    db = self.pasaradb(nlat)
                    self.atxval.append(db)
                    self.val = 0
                    
                    print(self.num_uniones.get())
                    x = int(self.num_uniones.get())-1
                    self.num_uniones.set(str(x))

                    self.guardar_fo(self.fibras[0],self.fibras[1],self.fibras[2],self.num_uniones.get(),self.atxco,self.atxuni)
                elif self.selmultlat.get() == 2:
                    nlat = self.lateralsingap_multi(y,self.nucleoFO)
                    db = self.pasaradb(nlat)
                    self.atxval.append(db)
                    self.val = 0
                    
                    print(self.num_uniones.get())
                    x = int(self.num_uniones.get())-1
                    self.num_uniones.set(str(x))

                    self.guardar_fo(self.fibras[0],self.fibras[1],self.fibras[2],self.num_uniones.get(),self.atxco,self.atxuni)
                else:
                    messagebox.showerror("Fibra Optica","Escoga si desea indice gradual")
            else:
                messagebox.showwarning("Fibra Optica","y no valido")
        except ValueError:
            messagebox.showerror("Fibra Optica","y no valido")

    def calculomult_latgap(self,y,indice_nucleo,indicegr):
        try:
            lim = self.fibras[3]*math.pow(10,6)
            y = float(y)
            if (y>0 and y<=lim):            
                self.frmmut_latg.destroy()
                y = y*math.pow(10,-6)
                if self.selmultlat.get() == 1:
                    lt = self.indicegradual(indicegr,y,self.nucleoFO)
                    nlat = 1-lt
                    db = self.pasaradb(nlat)
                    self.atxval.append(db)
                    self.val = 0
                    
                    print(self.num_uniones.get())
                    x = int(self.num_uniones.get())-1
                    self.num_uniones.set(str(x))

                    self.guardar_fo(self.fibras[0],self.fibras[1],self.fibras[2],self.num_uniones.get(),self.atxco,self.atxuni)
                elif self.selmultlat.get() == 2:
                    nlat = self.lateralgap_multi(y,indice_nucleo,self.nucleoFO)
                    db = self.pasaradb(nlat)
                    self.atxval.append(db)
                    self.val = 0
                    
                    print(self.num_uniones.get())
                    x = int(self.num_uniones.get())-1
                    self.num_uniones.set(str(x))

                    self.guardar_fo(self.fibras[0],self.fibras[1],self.fibras[2],self.num_uniones.get(),self.atxco,self.atxuni)
                else:
                    messagebox.showerror("Fibra Optica","Escoga si desea indice gradual")
            else:
                messagebox.showwarning("Fibra Optica","y no valido")
        except ValueError:
            messagebox.showerror("Fibra Optica","y no valido")

    def calculomult_ang(self,ang):
        try:
            ang = float(ang)
            if (ang>0 and ang<=10):
                self.frmmut_ang.destroy()
                nang = self.angular_multi(self.fibras[1][0],ang,self.fibras[6])
                db = self.pasaradb(nang)
                self.atxval.append(db)
                self.val = 0
                
                print(self.num_uniones.get())
                x = int(self.num_uniones.get())-1
                self.num_uniones.set(str(x))

                self.guardar_fo(self.fibras[0],self.fibras[1],self.fibras[2],self.num_uniones.get(),self.atxco,self.atxuni)
            else:
                messagebox.showwarning("Fibra Optica","𝜃 no valido")
        except ValueError:
            messagebox.showerror("Fibra Optica","𝜃 no valido")
    
    def calculomult_long(self,z):
        try:
            z = float(z)
            if (z>0 and z<=10):
                self.frmmut_long.destroy()
                z = z*math.pow(10,-6)
                nlong = self.longitudinal_multi(1,z,self.fibras[6],self.nucleoFO)
                self.atxval.append(nlong)
                self.val = 0
                
                print(self.num_uniones.get())
                x = int(self.num_uniones.get())-1
                self.num_uniones.set(str(x))

                self.guardar_fo(self.fibras[0],self.fibras[1],self.fibras[2],self.num_uniones.get(),self.atxco,self.atxuni)
            else:
                messagebox.showwarning("Fibra Optica","z no valido")
        except ValueError:
            messagebox.showerror("Fibra Optica","y no valido")
    

    def mono_lateralgui(self):
        self.frmatenuaciones.destroy()
        self.frmmono_lat = Frame()
        self.frmmono_lat .pack()
        self.frmmono_lat .config(width="400",height="230",background="mediumturquoise",bd=10,relief="raised")

        lbl = ttk.Label(self.frmmono_lat,text="Lateral Monomodo",font=(52),foreground="DarkOrange2")
        lbl.config(background="mediumturquoise")
        lbl.grid(row=0,column=1,columnspan=2)

        img = mpimg.imread('mon_lat.png')
        fig = Figure(figsize=(2,2),dpi=50)
        fig.add_subplot(111).imshow(img)

        canvas = FigureCanvasTkAgg(fig,master=self.frmmono_lat)
        canvas.draw()
        canvas.get_tk_widget().grid(row=1,column=0,columnspan=2)

        lbliy = ttk.Label(self.frmmono_lat,text="y en um",font=(18),foreground="DarkOrange2")
        lbliy.config(background="mediumturquoise")
        lbliy.grid(row=2,column=0)

        y = StringVar()
        ytxt= Entry(self.frmmono_lat,font=("consolas",10),textvariable=y,justify="center",background="lightsteelblue4",fg="lightcyan")
        ytxt.grid(row=2,column=1) 

        btnpv = ttk.Button(self.frmmono_lat,text="Perdidas",style="MyButton.TButton"
            ,command=lambda:self.calculomono_lat(ytxt.get()))
        btnpv.grid(row=3,column=1,columnspan=2)

    def mono_angulargui(self):
        self.frmatenuaciones.destroy()
        self.frmmono_ang = Frame()
        self.frmmono_ang .pack()
        self.frmmono_ang .config(width="400",height="230",background="mediumturquoise",bd=10,relief="raised")

        lbl = ttk.Label(self.frmmono_ang,text="Angular Monomodo",font=(52),foreground="DarkOrange2")
        lbl.config(background="mediumturquoise")
        lbl.grid(row=0,column=0,columnspan=2)

        img = mpimg.imread('mon_ang.png')
        fig = Figure(figsize=(2,2),dpi=50)
        fig.add_subplot(111).imshow(img)

        canvas = FigureCanvasTkAgg(fig,master=self.frmmono_ang)
        canvas.draw()
        canvas.get_tk_widget().grid(row=1,column=0,columnspan=2)

        lbliy = ttk.Label(self.frmmono_ang,text="𝜃° =",font=(18),foreground="DarkOrange2")
        lbliy.config(background="mediumturquoise")
        lbliy.grid(row=2,column=0)

        ang = StringVar()
        angtxt= Entry(self.frmmono_ang,font=("consolas",10),textvariable=ang,justify="center",background="lightsteelblue4",fg="lightcyan")
        angtxt.grid(row=2,column=1) 

        btnpv = ttk.Button(self.frmmono_ang,text="Perdidas",style="MyButton.TButton"
            ,command=lambda:self.calculomono_ang(angtxt.get()))
        btnpv.grid(row=3,column=0,columnspan=2)

    def mult_Lateralsingapgui(self):
        self.frmatenuaciones.destroy()
        self.frmmut_lat = Frame()
        self.frmmut_lat .pack()
        self.frmmut_lat .config(width="400",height="230",background="mediumturquoise",bd=10,relief="raised")

        lbl = ttk.Label(self.frmmut_lat,text="Lateral Multimodo sin gap",font=(52),foreground="DarkOrange2")
        lbl.config(background="mediumturquoise")
        lbl.grid(row=0,column=0,columnspan=4)

        img = mpimg.imread('mult_lat.png')
        fig = Figure(figsize=(2,2),dpi=50)
        fig.add_subplot(111).imshow(img)

        canvas = FigureCanvasTkAgg(fig,master=self.frmmut_lat)
        canvas.draw()
        canvas.get_tk_widget().grid(row=1,column=1,columnspan=2)

        lblingradual = ttk.Label(self.frmmut_lat,text="fibra indice gradual",font=(18),foreground="DarkOrange2")
        lblingradual.config(background="mediumturquoise")
        lblingradual.grid(row=2,column=0,columnspan=2)

        rbsiaper = Radiobutton(self.frmmut_lat,background="mediumturquoise",foreground="DarkOrange2",text="SI",variable=self.selmultlat,value=1)
        rbsiaper.grid(row=3,column=0)

        rbsiaper = Radiobutton(self.frmmut_lat,background="mediumturquoise",foreground="DarkOrange2",text="NO",variable=self.selmultlat,value=2)
        rbsiaper.grid(row=3,column=1)

        lblingradual = ttk.Label(self.frmmut_lat,text="𝛼 =",font=(18),foreground="DarkOrange2")
        lblingradual.config(background="mediumturquoise")
        lblingradual.grid(row=2,column=2)

        x = StringVar()
        x.set(self.fibras[8])
        perfiltxtx = Entry(self.frmmut_lat,font=("consolas",10),textvariable=x,justify="center",background="lightsteelblue4",fg="lightcyan")
        perfiltxtx.grid(row=2,column=3)

        lbliy = ttk.Label(self.frmmut_lat,text="y en um",font=(18),foreground="DarkOrange2")
        lbliy.config(background="mediumturquoise")
        lbliy.grid(row=3,column=2)

        y = StringVar()
        ytxt= Entry(self.frmmut_lat,font=("consolas",10),textvariable=y,justify="center",background="lightsteelblue4",fg="lightcyan")
        ytxt.grid(row=3,column=3) 

        btnpv = ttk.Button(self.frmmut_lat,text="Perdidas",style="MyButton.TButton"
            ,command=lambda:self.calculomult_latsingap(ytxt.get(),self.fibras[8]))
        btnpv.grid(row=4,column=3)

    def mult_Lateralgapgui(self):
        self.frmatenuaciones.destroy()
        self.frmmut_latg = Frame()
        self.frmmut_latg .pack()
        self.frmmut_latg .config(width="400",height="230",background="mediumturquoise",bd=10,relief="raised")

        lbl = ttk.Label(self.frmmut_latg,text="Lateral Multimodo gap",font=(52),foreground="DarkOrange2")
        lbl.config(background="mediumturquoise")
        lbl.grid(row=0,column=0,columnspan=4)

        img = mpimg.imread('mult_lat.png')
        fig = Figure(figsize=(2,2),dpi=60)
        fig.add_subplot(111).imshow(img)

        canvas = FigureCanvasTkAgg(fig,master=self.frmmut_latg)
        canvas.draw()
        canvas.get_tk_widget().grid(row=1,column=1,columnspan=2)

        lblingradual = ttk.Label(self.frmmut_latg,text="fibra indice gradual",font=(18),foreground="DarkOrange2")
        lblingradual.config(background="mediumturquoise")
        lblingradual.grid(row=2,column=0,columnspan=2)

        rbsiaper = Radiobutton(self.frmmut_latg,background="mediumturquoise",foreground="DarkOrange2",text="SI",variable=self.selmultlat,value=1)
        rbsiaper.grid(row=3,column=0)

        rbsiaper = Radiobutton(self.frmmut_latg,background="mediumturquoise",foreground="DarkOrange2",text="NO",variable=self.selmultlat,value=2)
        rbsiaper.grid(row=3,column=1)

        lblingradual = ttk.Label(self.frmmut_latg,text="𝛼 =",font=(18),foreground="DarkOrange2")
        lblingradual.config(background="mediumturquoise")
        lblingradual.grid(row=2,column=2)

        x = StringVar()
        x.set(self.fibras[8])
        perfiltxtx = Entry(self.frmmut_latg,font=("consolas",10),textvariable=x,justify="center",background="lightsteelblue4",fg="lightcyan")
        perfiltxtx.grid(row=2,column=3)

        lbliy = ttk.Label(self.frmmut_latg,text="y en um",font=(18),foreground="DarkOrange2")
        lbliy.config(background="mediumturquoise")
        lbliy.grid(row=3,column=2)

        y = StringVar()
        ytxt= Entry(self.frmmut_latg,font=("consolas",10),textvariable=y,justify="center",background="lightsteelblue4",fg="lightcyan")
        ytxt.grid(row=3,column=3) 

        btnpv = ttk.Button(self.frmmut_latg,text="Perdidas",style="MyButton.TButton"
            ,command=lambda:self.calculomult_latgap(ytxt.get(),self.fibras[1][0],self.fibras[8]))
        btnpv.grid(row=4,column=3)

        # btnpv = ttk.Button(self.frmmut_latg,text="Volver",style="MyButton.TButton"
        #     ,command=lambda:self.calculomult_latgap(ytxt.get(),self.fibras[1][0],self.fibras[4]))
        # btnpv.grid(row=4,column=0)


    def mult_angulargui(self):
        self.frmatenuaciones.destroy()
        self.frmmut_ang = Frame()
        self.frmmut_ang .pack()
        self.frmmut_ang .config(width="400",height="230",background="mediumturquoise",bd=10,relief="raised")

        lbl = ttk.Label(self.frmmut_ang,text="Angular Multimodo",font=(52),foreground="DarkOrange2")
        lbl.config(background="mediumturquoise")
        lbl.grid(row=0,column=0,columnspan=2)

        img = mpimg.imread('mult_angu.png')
        fig = Figure(figsize=(2,2),dpi=60)
        fig.add_subplot(111).imshow(img)

        canvas = FigureCanvasTkAgg(fig,master=self.frmmut_ang)
        canvas.draw()
        canvas.get_tk_widget().grid(row=1,column=0,columnspan=2)

        lbliy = ttk.Label(self.frmmut_ang,text="𝜃° =",font=(18),foreground="DarkOrange2")
        lbliy.config(background="mediumturquoise")
        lbliy.grid(row=2,column=0)

        ang = StringVar()
        angtxt= Entry(self.frmmut_ang,font=("consolas",10),textvariable=ang,justify="center",background="lightsteelblue4",fg="lightcyan")
        angtxt.grid(row=2,column=1) 

        btnpv = ttk.Button(self.frmmut_ang,text="Perdidas",style="MyButton.TButton"
            ,command=lambda:self.calculomult_ang(angtxt.get()))
        btnpv.grid(row=3,column=1)

    def mult_longitudinalgui(self):
        self.frmatenuaciones.destroy()
        self.frmmut_long = Frame()
        self.frmmut_long.pack()
        self.frmmut_long.config(width="400",height="230",background="mediumturquoise",bd=10,relief="raised")

        lbl = ttk.Label(self.frmmut_long,text="Longitudinal Multimodo",font=(52),foreground="DarkOrange2")
        lbl.config(background="mediumturquoise")
        lbl.grid(row=0,column=0,columnspan=2)
        
        img = mpimg.imread('mult_long.png')
        fig = Figure(figsize=(2,2),dpi=60)
        fig.add_subplot(111).imshow(img)

        canvas = FigureCanvasTkAgg(fig,master=self.frmmut_long)
        canvas.draw()
        canvas.get_tk_widget().grid(row=1,column=0,columnspan=2)

        lbl = ttk.Label(self.frmmut_long,text="Z(um) =",font=(18),foreground="DarkOrange2")
        lbl.config(background="mediumturquoise")
        lbl.grid(row=2,column=0)

        z = StringVar()
        ztxt= Entry(self.frmmut_long,font=("consolas",10),textvariable=z,justify="center",background="lightsteelblue4",fg="lightcyan")
        ztxt.grid(row=2,column=1)

        btnpv = ttk.Button(self.frmmut_long,text="Perdidas",style="MyButton.TButton"
            ,command=lambda:self.calculomult_long(ztxt.get()))
        btnpv.grid(row=3,column=1)
    def hallarruido(self,corrs,id,rl,t,bw):
        try:
            bw = float(bw)
            if bw >0:
                bw = bw *math.pow(10,6)
                corrs = float(corrs)
                id = float(id)*math.pow(10,-6)
                rl = self.fuentepara[4]
                t = float(t)
                e = 1.6021*math.pow(10,-19)
                k = 1.38*math.pow(10,-23)
                pns = 2*e*bw*(corrs+id)*rl
                pnt = 4*k*t*bw
                pruido = pns+pnt
                self.pnoise.append(pns)
                self.pnoise.append(pnt)
                self.pnoise.append(pruido)
                self.frmruido.pack_forget()
                self.frmfinal.pack()
                self.vecesfinal = 2
                
                self.fuentepara.append(id)

                Pruido = "Pnoise = "+ str(pruido) + "W"
                lblsnrmin = ttk.Label(self.frmfinal,text=Pruido,font=(18),foreground="DarkOrange2")
                lblsnrmin.config(background="mediumturquoise")
                lblsnrmin.grid(row=4,column=1)

            else:
                messagebox.showerror("Fibra Optica","Ingrese otro ancho de banda")
        except ValueError:
            messagebox.showerror("Fibra Optica","Ingrese valor acorde")
    def hallarfuente(self,Vo,responsividad,rl,pe,n,ganancia):
        try:
            #si ahi cmbio aca
            m = float(ganancia)
            Vo = float(Vo)
            if n == "x":
                responsividad = float(responsividad)
                e = 1.6021*math.pow(10,-19)
                c = 3*math.pow(10,8)
                h =6.626*math.pow(10,-34)
                l = self.fibras[7]*math.pow(10,-9)
                f = (c/l)
                if float(responsividad) >=0 or float(responsividad) < 0:
                        nefi = (responsividad*h*c)/(e*l)

                        Pomax = (responsividad*Vo)/float(rl)
                        corriente_s=(m*nefi*e*Pomax)/(h*f)
                        self.fuentepara.append(Pomax)
                        self.fuentepara.append(corriente_s)
                        self.fuentepara.append(Vo)
                        self.fuentepara.append(float(responsividad))
                        if pe == "x":
                            pele = (math.pow(corriente_s,2))*float(rl)
                            self.fuentepara.append(pele)
                            self.fuentepara.append(float(rl))

                            self.Pein.set(str(pele))
                            print(self.fuentepara)

                            self.frmfuente.pack_forget()
                            self.frmfinal.pack()
                            lblsnrmin = ttk.Label(self.frmfinal,text="SNR min = ",font=(18),foreground="DarkOrange2")
                            lblsnrmin.config(background="mediumturquoise")
                            lblsnrmin.grid(row=3,column=1)

                            snrmintxt= Entry(self.frmfinal,font=("consolas",10),textvariable=self.snrmin,justify="center",background="lightsteelblue4",fg="lightcyan")
                            snrmintxt.grid(row=3,column=2)

                            btnsnr = ttk.Button(self.frmfinal,text="Perdidas",style="MyButton.TButton"
                                ,command=lambda:[self.ruido()])
                            btnsnr.grid(row=4,column=2,pady=2)
                            self.vecesfinal = 1

                        elif float(pe)>0 or float(pe)<100:
                            pele = round((math.pow(corriente_s,2))*float(rl),4)
                            self.fuentepara.append(pele)
                            self.fuentepara.append(float(rl))

                            self.Pein.set(str(pele))
                            print(self.fuentepara)

                            self.frmfuente.pack_forget()
                            self.frmfinal.pack()
                            lblsnrmin = ttk.Label(self.frmfinal,text="SNR min = ",font=(18),foreground="DarkOrange2")
                            lblsnrmin.config(background="mediumturquoise")
                            lblsnrmin.grid(row=3,column=1)

                            snrmintxt= Entry(self.frmfinal,font=("consolas",10),textvariable=self.snrmin,justify="center",background="lightsteelblue4",fg="lightcyan")
                            snrmintxt.grid(row=3,column=2)

                            btnsnr = ttk.Button(self.frmfinal,text="HallarRuido",style="MyButton.TButton"
                                ,command=lambda:[self.ruido()])
                            btnsnr.grid(row=4,column=2,pady=2)
                            self.vecesfinal = 1
                else:
                    messagebox.showerror("Fibra Optica","Ingrese eficiencia Optica")
                
            elif float(n)>0 or float(n<100):
                if responsividad == "x":
                    n = float(n)*(1/100)
                    e = 1.6021*math.pow(10,-19)
                    c = 3*math.pow(10,8)
                    h =6.626*math.pow(10,-34)
                    l = self.fibras[7]*math.pow(10,-9)
                    f = (c/l)
                    res = (n*e*l)/(h*c)
                    
                    Pomax = (res*Vo)/float(rl)
                    corriente_s=(m*n*e*Pomax)/(h*f)

                    self.fuentepara.append(Pomax)
                    self.fuentepara.append(corriente_s)
                    self.fuentepara.append(Vo)

                    if pe == "x" :
                        pele = round((math.pow(corriente_s,2))*float(rl),4)
                        self.fuentepara.append(pele)
                        self.fuentepara.append(float(rl))

                        self.Pein.set(str(pele))
                        print(self.fuentepara)

                        self.frmfuente.pack_forget()
                        self.frmfinal.pack()
                        lblsnrmin = ttk.Label(self.frmfinal,text="SNR min = ",font=(18),foreground="DarkOrange2")
                        lblsnrmin.config(background="mediumturquoise")
                        lblsnrmin.grid(row=3,column=1)

                        snrmintxt= Entry(self.frmfinal,font=("consolas",10),textvariable=self.snrmin,justify="center",background="lightsteelblue4",fg="lightcyan")
                        snrmintxt.grid(row=3,column=2)

                        btnsnr = ttk.Button(self.frmfinal,text="HallarRuido",style="MyButton.TButton"
                            ,command=lambda:[self.ruido()])
                        btnsnr.grid(row=4,column=2,pady=2)
                        self.vecesfinal = 1
                    elif float(pe)>0 or float(pe)<100:
                        pele = round((math.pow(corriente_s,2))*float(rl),4)
                        self.fuentepara.append(pele)
                        self.fuentepara.append(float(rl))

                        self.Pein.set(str(pele))
                        print(self.fuentepara)

                        self.frmfuente.pack_forget()
                        self.frmfinal.pack()
                        lblsnrmin = ttk.Label(self.frmfinal,text="SNR min = ",font=(18),foreground="DarkOrange2")
                        lblsnrmin.config(background="mediumturquoise")
                        lblsnrmin.grid(row=3,column=1)

                        snrmintxt= Entry(self.frmfinal,font=("consolas",10),textvariable=self.snrmin,justify="center",background="lightsteelblue4",fg="lightcyan")
                        snrmintxt.grid(row=3,column=2)

                        btnsnr = ttk.Button(self.frmfinal,text="HallarRuido",style="MyButton.TButton"
                            ,command=lambda:[self.ruido()])
                        btnsnr.grid(row=4,column=2,pady=2)
                        self.vecesfinal = 1

                else:
                    messagebox.showerror("Fibra Optica","Ingrese responsividad Optica")
            else: 
                messagebox.showerror("Fibra Optica","Valor n entre 0 a 100")
        except ValueError:
            messagebox.showerror("Fibra Optica","Ingrese valor acorde")
    

    def hallarpot(self,pin,pout):
        try:
            snrmin = self.snrmin.get()
            tsnrmin = "SNR min = " + snrmin
            lblsnrmin = ttk.Label(self.frmfinal,text=tsnrmin,font=(18),foreground="DarkOrange2")
            lblsnrmin.config(background="mediumturquoise")
            lblsnrmin.grid(row=3,column=1)
            atx = 0

            for i in range(len(self.atxval)):
                atx = atx + self.atxval[i]
            atx = atx+ self.atenuaciones
            if pout == "x":
                pin = float(pin)
                if pin >= 0 or pin<= 0:
                    messagebox.showinfo("Fibra Optica","Hallo Pout")
                    dbmpin = self.pasardbm(pin) 
                    p = dbmpin - atx
                    r = "Pout = " +str(p) +"dBm"
                    p1 = self.pasarapot(p)
                    r1 =  str(p1) +"W"
                    

                    snr = (math.pow(10,(p/10))*math.pow(10,-3))/self.pnoise[2]
                    snrmi = float(self.snrmin.get())

                    txt = "SNR = "+ str(snr)
                    lblpin = ttk.Label(self.frmfinal,text=txt,font=(18),foreground="DarkOrange2")
                    lblpin.config(background="mediumturquoise")
                    lblpin.grid(row=5,column=1,columnspan=2)
                    self.Pout.set(r1)
                    if self.fuentepara[1]> self.fuentepara[6]:
                        if snr > snrmi:
                            messagebox.showinfo("Fibra Optica","funciona el enlace")
                        else:
                            men = "Los parametros no cumplen con el snr min " + self.snrmin.get()+"\n"
                            men = men + "Cambie la fuente Optica "
                            self.Pein.set("Cambiar fuente Optica")
                            messagebox.showinfo("Fibra Optica",men)
                            self.fuenteoptica()
                    else:
                        men = "Los parametros no cumplen con la Id min " + str(self.fuentepara[5])+"\n"
                        men = men + "Cambie la fuente Optica "
                        self.Pein.set("Cambiar fuente Optica")
                        messagebox.showinfo("Fibra Optica",men)
                        self.fuenteoptica()
                    self.fuentepara = []
                    
                    # lblpout = ttk.Label(self.frmfinal,text=r,font=(18),foreground="DarkOrange2")
                    # lblpout.config(background="mediumturquoise")
                    # lblpout.grid(row=3,column=1)

                    # lblpoutW = ttk.Label(self.frmfinal,text=r1,font=(18),foreground="DarkOrange2")
                    # lblpoutW.config(background="mediumturquoise")
                    # lblpoutW.grid(row=4,column=1)
            elif float(pout)>0 or float(pout)<=0:
                pin = float(pin)
                if pin >= 0 or pin<= 0:
                    messagebox.showinfo("Fibra Optica","Hallo Pout")
                    dbmpin = self.pasardbm(pin) 
                    p = dbmpin - atx
                    r = "Pout = " +str(p) +"dBm"
                    p1 = self.pasarapot(p)
                    r1 =  str(p1) +"W"
                    

                    snr = (math.pow(10,(p/10))*math.pow(10,-3))/self.pnoise[2]
                    snrmi = float(self.snrmin.get())
                    if snr > snrmi:
                        self.Pout.set(r1)
                        txt = "SNR ="+ str(snr)
                        lblpin = ttk.Label(self.frmfinal,text=txt,font=(18),foreground="DarkOrange2")
                        lblpin.config(background="mediumturquoise")
                        lblpin.grid(row=5,column=1,columnspan=2)
                    else:
                        self.Pout.set(r1)
                        men = "Los parametros no cumplen con el snr min " + self.snrmin.get()+"\n"
                        men = men + "Cambie la RL o POptica "
                        self.Pein.set("Cambiar RL o Potencia Optica")
                        messagebox.showinfo("Fibra Optica",men)
                        self.fuenteoptica()
                    self.fuentepara = []
            elif pin == "x":
                pout = float(pout)
                if pout >= 0 or pout<= 0:
                    messagebox.showinfo("Fibra Optica","Hallo Pin")
                    dbmpout = self.pasardbm(pout) 
                    p = dbmpout+atx
                    r = "Pin = " +str(p)+" dBm"

                    p1 = self.pasarapot(p)
                    r1 = "Pin = " +str(p1) +" W"

                    # lblpin = ttk.Label(self.frmfinal,text=r,font=(18),foreground="DarkOrange2")
                    # lblpin.config(background="mediumturquoise")
                    # lblpin.grid(row=3,column=1)

                    # lblpinW = ttk.Label(self.frmfinal,text=r1,font=(18),foreground="DarkOrange2")
                    # lblpinW.config(background="mediumturquoise")
                    # lblpinW.grid(row=4,column=1)
                    self.fuentepara = []


            else:
                pout = float(pout)
                pin = float(pin)
                if pin >= 0 and pin<= 0:
                    if pout >= 0 and pout <= 0:
                        messagebox.showwarning("Fibra Optica","Debe ingresar solo una potencia")
                elif pout >= 0 and pout <= 0:
                    if pin >= 0 and pin<= 0:
                        messagebox.showwarning("Fibra Optica","Debe ingresar solo una potencia")
                else:
                    messagebox.showwarning("Fibra Optica","Ingrese una potencia")
        except ValueError:
            messagebox.showerror("Fibra Optica","Ingrese valor de potencias")

    # opciones de atenuaciones
    def perdidas(self,tipo,perdidaM):
        if tipo == "Monomodo":
            if perdidaM == "Lateral":
                self.mono_lateralgui()
                self.tipodeatx.append("SMlat")
            elif perdidaM  == "Angular":
                self.mono_angulargui()
                self.tipodeatx.append("SMang")
        elif tipo == "Multimodo":
            if perdidaM == "Lateral sin gap":
                self.mult_Lateralsingapgui()
                self.tipodeatx.append("MMlat")
            elif perdidaM  == "Lateral con gap":
                self.mult_Lateralgapgui()
                self.tipodeatx.append("MMlatgap")
            elif perdidaM  == "Angular":
                self.mult_angulargui()
                self.tipodeatx.append("MMang")
            elif perdidaM == "Longitudinal":
                self.mult_longitudinalgui()
                self.tipodeatx.append("MMlong")
    
    def ruido(self):
        self.frmfinal.pack_forget()
        self.frmruido = Frame()
        self.frmruido.pack()
        self.frmruido.config(width="400",height="230",background="mediumturquoise",bd=10,relief="raised")


        lbl = ttk.Label(self.frmruido,text="Ruido",font=(52),foreground="DarkOrange2")
        lbl.config(background="mediumturquoise")
        lbl.grid(row=0,column=0,columnspan=4)

        lblIs = ttk.Label(self.frmruido,text="Is = ",font=(18),foreground="DarkOrange2")
        lblIs.config(background="mediumturquoise")
        lblIs.grid(row=1,column=0)

        is1 = StringVar()
        is1.set(str(self.fuentepara[1]))
        corrstxt= Entry(self.frmruido,font=("consolas",10),textvariable=is1,justify="center",background="lightsteelblue4",fg="lightcyan")
        corrstxt.grid(row=1,column=1) 

        lblId = ttk.Label(self.frmruido,text="Id(µA) = ",font=(18),foreground="DarkOrange2")
        lblId.config(background="mediumturquoise")
        lblId.grid(row=2,column=0)

        isd = StringVar()
        idtxt= Entry(self.frmruido,font=("consolas",10),textvariable=isd,justify="center",background="lightsteelblue4",fg="lightcyan")
        idtxt.grid(row=2,column=1) 

        lblrl = ttk.Label(self.frmruido,text="RL = ",font=(18),foreground="DarkOrange2")
        lblrl.config(background="mediumturquoise")
        lblrl.grid(row=3,column=0)

        rl = StringVar()
        rl.set(str(self.fuentepara[4])+"Ω")
        rltxt= Entry(self.frmruido,font=("consolas",10),textvariable=rl,justify="center",background="lightsteelblue4",fg="lightcyan")
        rltxt.grid(row=3,column=1) 

        lblt = ttk.Label(self.frmruido,text="T(°k) = ",font=(18),foreground="DarkOrange2")
        lblt.config(background="mediumturquoise")
        lblt.grid(row=1,column=2)

        temp = StringVar()
        temptxt= Entry(self.frmruido,font=("consolas",10),textvariable=temp,justify="center",background="lightsteelblue4",fg="lightcyan")
        temptxt.grid(row=1,column=3) 

        lblbw = ttk.Label(self.frmruido,text="Ancho banda Mhz = ",font=(18),foreground="DarkOrange2")
        lblbw.config(background="mediumturquoise")
        lblbw.grid(row=2,column=2)

        bw = StringVar()
        bwtxt= Entry(self.frmruido,font=("consolas",10),textvariable=bw,justify="center",background="lightsteelblue4",fg="lightcyan")
        bwtxt.grid(row=2,column=3) 

        btnsnr = ttk.Button(self.frmruido,text="Perdidas",style="MyButton.TButton"
            ,command=lambda:[self.hallarruido(corrstxt.get(),idtxt.get(),rltxt.get(),temptxt.get(),bwtxt.get())])
        btnsnr.grid(row=3,column=2,columnspan=2,pady=2)

        snrmin = self.snrmin.get()
        tsnrmin = "SNR min = " + snrmin
        lblsnrmin = ttk.Label(self.frmfinal,text=tsnrmin,font=(18),foreground="DarkOrange2")
        lblsnrmin.config(background="mediumturquoise")
        lblsnrmin.grid(row=3,column=1)

    def fuenteoptica(self):

        self.frmfinal.pack_forget()
        self.frmfuente = Frame()
        self.frmfuente.pack()
        self.frmfuente.config(width="400",height="230",background="mediumturquoise",bd=10,relief="raised")

        lbl = ttk.Label(self.frmfuente,text="Fuente Optica",font=(52),foreground="DarkOrange2")
        lbl.config(background="mediumturquoise")
        lbl.grid(row=0,column=0,columnspan=5)

        lblpoled = ttk.Label(self.frmfuente,text="Voltaje",font=(18),foreground="DarkOrange2")
        lblpoled.config(background="mediumturquoise")
        lblpoled.grid(row=1,column=0)

        
        votxt= Entry(self.frmfuente,font=("consolas",10),textvariable=self.vo,justify="center",background="lightsteelblue4",fg="lightcyan")
        votxt.grid(row=2,column=0) 

        img = mpimg.imread('FuenteOptica.png')
        fig = Figure(figsize=(2,2),dpi=70)
        fig.add_subplot(111).imshow(img)

        canvas = FigureCanvasTkAgg(fig,master=self.frmfuente)
        canvas.draw()
        canvas.get_tk_widget().grid(row=1,column=1,rowspan=3)

        lblRL = ttk.Label(self.frmfuente,text="RL = ",font=(18),foreground="DarkOrange2")
        lblRL.config(background="mediumturquoise")
        lblRL.grid(row=2,column=2)

        
        rltxt= Entry(self.frmfuente,font=("consolas",10),textvariable=self.rl,justify="center",background="lightsteelblue4",fg="lightcyan")
        rltxt.grid(row=2,column=3)

        lblres = ttk.Label(self.frmfuente,text=" responsividad(ρ) = ",font=(18),foreground="DarkOrange2")
        lblres.config(background="mediumturquoise")
        lblres.grid(row=4,column=0)

        
        restxt= Entry(self.frmfuente,font=("consolas",10),textvariable=self.res,justify="center",background="lightsteelblue4",fg="lightcyan")
        restxt.grid(row=4,column=1,pady=2)

        lblpe = ttk.Label(self.frmfuente,text=" Pe = ",font=(18),foreground="DarkOrange2")
        lblpe.config(background="mediumturquoise")
        lblpe.grid(row=4,column=2)

        pe = StringVar()
        pe.set("x")
        petxt= Entry(self.frmfuente,font=("consolas",10),textvariable=pe,justify="center",background="lightsteelblue4",fg="lightcyan")
        petxt.grid(row=4,column=3,pady=2)

        lbln = ttk.Label(self.frmfuente,text=" efectividad (0-100)% = ",font=(18),foreground="DarkOrange2")
        lbln.config(background="mediumturquoise")
        lbln.grid(row=5,column=0)

        
        ntxt = Entry(self.frmfuente,font=("consolas",10),textvariable=self.n,justify="center",background="lightsteelblue4",fg="lightcyan")
        ntxt.grid(row=5,column=1,pady=2)

        lblnumdinodos = ttk.Label(self.frmfuente,text=" Ganancia M = ",font=(18),foreground="DarkOrange2")
        lblnumdinodos.config(background="mediumturquoise")
        lblnumdinodos.grid(row=6,column=0)

        
        ndinodostxt = Entry(self.frmfuente,font=("consolas",10),textvariable=self.ganancia,justify="center",background="lightsteelblue4",fg="lightcyan")
        ndinodostxt.grid(row=6,column=1,pady=2)

        btnhallarpe = ttk.Button(self.frmfuente,text="Calcular Pe",style="MyButton.TButton"
            ,command=lambda:[self.hallarfuente(votxt.get(),restxt.get(),rltxt.get(),petxt.get(),ntxt.get(),ndinodostxt.get())])
        btnhallarpe.grid(row=5,column=3,columnspan=2,rowspan=3)

        messagebox.showinfo("Fibra Optica","Ingrese una x si quiere hallar responsividad sino ingrese valor")
        

    def fibraopfinal(self):

        self.frmatenuaciones.destroy()
        self.frmfinal = Frame()
        self.frmfinal.pack()
        self.frmfinal.config(width="400",height="230",background="mediumturquoise",bd=10,relief="raised")

        lbl = ttk.Label(self.frmfinal,text="Sistema de Fibra Optica",font=(52),foreground="DarkOrange2")
        lbl.config(background="mediumturquoise")
        lbl.grid(row=0,column=0,columnspan=4)
        max = 500
        w = Canvas(self.frmfinal, width=max, height=35)
        w.grid(row=1,column=1,rowspan=2,columnspan=2)
        ini = 0
        final = 0
        p = ((1/self.uniones)*max)
        for i in range(self.uniones):
            if i == self.uniones-1:
                ini = final
                final = max
            else:
                ini = final
                final = final+p
                
            if self.fibras[0] == "Monomodo":
                relleno = "yellow"
            elif self.fibras[0] == "Multimodo":
                relleno = "orange"
            w.create_rectangle(ini, 0, final, 21, fill=relleno)

            tipoatx = self.tipodeatx[i]
            valor = str(round(self.atxval[i],2))+"dB"
            centro = round(((final-ini)/2))+ini
            w.create_text(centro, 13, text=tipoatx, fill='black')
            w.create_text(centro, 30, text=valor, fill='DarkOrange4')

        lblpin = ttk.Label(self.frmfinal,text="Pin = W",font=(18),foreground="DarkOrange2")
        lblpin.config(background="mediumturquoise")
        lblpin.grid(row=1,column=0)

        Pintxt= Entry(self.frmfinal,font=("consolas",10),textvariable=self.Pein,justify="center",background="lightsteelblue4",fg="lightcyan")
        Pintxt.grid(row=2,column=0) 

        btnpv = ttk.Button(self.frmfinal,text="Perdidas",style="MyButton.TButton"
            ,command=lambda:[self.hallarpot(Pintxt.get(),Pouttxt.get())])
        
        btnpein = ttk.Button(self.frmfinal,text="Calcular Pe-entrada",style="MyButton.TButton"
            ,command=lambda:[btnpv.grid(row=3,column=3),self.fuenteoptica()])
        btnpein.grid(row=3,column=0,pady=2)

        lblpout = ttk.Label(self.frmfinal,text="Pout = W",font=(18),foreground="DarkOrange2")
        lblpout.config(background="mediumturquoise")
        lblpout.grid(row=1,column=3)

        self.Pout.set("x")
        Pouttxt= Entry(self.frmfinal,font=("consolas",10),textvariable=self.Pout,justify="center",background="lightsteelblue4",fg="lightcyan")
        Pouttxt.grid(row=2,column=3) 
        
        btnvolver = ttk.Button(self.frmfinal,text="Volver",style="MyButton.TButton"
            ,command=lambda:[self.frmfinal.destroy(),self.fibraoptica_parv(self.fibras[1],self.fibras[7],self.fibras[3],self.fibras[8],str(self.tamaFO))])
        btnvolver.grid(row=4,column=3,pady=5)
        

    #------------Interfaz Grafica----------------------------------------------------------------------------------------------------------------------------
    
    def guardar_fo(self,tipo,indices,param_v,num_uniones,atxc,atx):
        try:
            num_uniones = int(num_uniones)
            
            if num_uniones>0:
                if num_uniones > 0 and num_uniones<=5:
                    atxconectores = float(atxc)*2
                    self.atxco = str(atxconectores)
                    atxunion = float(atx)*self.tamaFO
                    self.atxuni = str(atxunion)
                    if self.val == 1:
                        self.unionnum = 1 
                        self.uniones = num_uniones
                        self.fibras.insert(0,tipo)
                        self.fibras.insert(1,indices)
                        self.fibras.insert(2,param_v)
                        self.fibras.insert(3,self.nucleoFO)
                        self.fibras.append(float(self.lambdafo.get()))
                        self.fibras.append(self.perfil)
                        modos = self.numero_modos(self.perfil,param_v)
                        self.fibras.append(modos)
                        self.atenuaciones = atxconectores+atxunion
                        print(self.fibras)


                    self.frmuniones.destroy()
                    self.frmatenuaciones = Frame()
                    self.frmatenuaciones.pack()
                    self.frmatenuaciones.config(width="400",height="230",background="mediumturquoise",bd=10,relief="raised")

                    lbl = ttk.Label(self.frmatenuaciones,text="Fibra Optica Características",font=(52),foreground="DarkOrange2")
                    lbl.config(background="mediumturquoise")
                    lbl.grid(row=0,column=0,columnspan=4)
                    
                    percon = "Perdidas por conectores = "+ str(atxconectores)
                    lblPerdidasunion = ttk.Label(self.frmatenuaciones,text=percon,font=(18),foreground="DarkOrange2")
                    lblPerdidasunion.config(background="mediumturquoise")
                    lblPerdidasunion.grid(row=3,column=0)
                    
                    peruni = "Perdidas por uniones = "+ str(atxunion)
                    lblPerdidasunion = ttk.Label(self.frmatenuaciones,text=peruni,font=(18),foreground="DarkOrange2")
                    lblPerdidasunion.config(background="mediumturquoise")
                    lblPerdidasunion.grid(row=4,column=0)
                    if tipo == "Monomodo":
                        
                        s = "Desalineacion #"+str(self.unionnum) +" en Monomodo"
                        lbltipo_at = ttk.Label(self.frmatenuaciones,text=s,font=(18),foreground="DarkOrange2")
                        lbltipo_at.config(background="mediumturquoise")
                        lbltipo_at.grid(row=1,column=0)


                        self.cbmmon= ttk.Combobox(self.frmatenuaciones,state="readondly",values=["Lateral","Angular"])
                        self.cbmmon.grid(row=2,column=0)
                        self.cbmmon.current(0)

                        btnpv = ttk.Button(self.frmatenuaciones,text="Union",style="MyButton.TButton"
                            ,command=lambda:self.perdidas(tipo,self.cbmmon.get()))
                        btnpv.grid(row=5,column=0,columnspan=2)
                        
                    elif tipo == "Multimodo":
                        s = "Desalineacion #"+str(self.unionnum) +" en Multimodo"
                        lbltipo_at = ttk.Label(self.frmatenuaciones,text=s,font=(18),foreground="DarkOrange2")
                        lbltipo_at.config(background="mediumturquoise")
                        lbltipo_at.grid(row=1,column=0)

                        self.cbmmult= ttk.Combobox(self.frmatenuaciones,state="readondly",values=["Lateral sin gap","Lateral con gap","Angular","Longitudinal"])
                        self.cbmmult.grid(row=2,column=0)
                        self.cbmmult.current(0)

                        btnpv = ttk.Button(self.frmatenuaciones,text="Union",style="MyButton.TButton"
                            ,command=lambda:self.perdidas(tipo,self.cbmmult.get()))
                        btnpv.grid(row=5,column=0,columnspan=2)

                    
                    

                    self.unionnum += 1 
                else:
                    messagebox.showwarning("Fibra Optica","Se permiten entre 1 a 5 uniones")
                    self.parv.set("")
                    self.num_uniones.set("")
            else:
                self.fibraopfinal()
                self.vecesfinal = 0
                # messagebox.showinfo("Fibra Optica","Ingrese la potencia que necesita con una X")
        except ValueError:
            messagebox.showerror("Fibra Optica","Ingrese el número de uniones")
            self.num_uniones.set("")
        

    def fibraoptica_parv(self,indices,lambda_fo,nucleo,perfil,tamFO):
        try:
            if float(lambda_fo) >= 800 and float(lambda_fo) <= 1600:
                self.pnoise = []
                self.fuentepara = []
                self.atxval = []
                self.tipodeatx = []
                self.atenuaciones = 0
                if float(tamFO)>0:
                    if self.seleccion >= 1:
                        self.tamaFO = float(tamFO)
                        self.frminicio.destroy()
                        if self.seleccion==1:
                            i = []
                            i.append(float(indices[1:5]))
                            i.append(float(indices[6:10]))
                            self.parametro_v(((float(nucleo)/2)*math.pow(10,-6)),i,lambda_fo,perfil,tamFO)
                        else:
                            self.frmuniones = Frame()
                            self.frmuniones.grid(row=0,column=0)
                            self.frmuniones.config(width="400",height="230",background="mediumturquoise",bd=10,relief="raised")

                            lbl = ttk.Label(self.frmuniones,text="Fibra Optica Uniones",font=(52),foreground="DarkOrange2")
                            lbl.config(background="mediumturquoise")
                            lbl.grid(row=0,column=0,columnspan=4)

                            lblPV = ttk.Label(self.frmuniones,text="Parametro V = ",font=(18),foreground="DarkOrange2")
                            lblPV.config(background="mediumturquoise")
                            lblPV.grid(row=1,column=0)


                            pvtxt = Entry(self.frmuniones,font=("consolas",10),textvariable=self.parv,justify="center",background="lightsteelblue4",fg="lightcyan")
                            pvtxt.grid(row=1,column=1)

                            lblnum_uniones = ttk.Label(self.frmuniones,text="Número de uniones",font=(18),foreground="DarkOrange2")
                            lblnum_uniones.config(background="mediumturquoise")
                            lblnum_uniones.grid(row=3,column=0)

                            numatetxt = Entry(self.frmuniones,font=("consolas",10),textvariable=self.num_uniones,justify="center",background="lightsteelblue4",fg="lightcyan")
                            numatetxt.grid(row=3,column=1) 

                            
                            lblperdidas = ttk.Label(self.frmuniones,text="Fibra Optica Perdidas",font=(52),foreground="DarkOrange2")
                            lblperdidas.config(background="mediumturquoise")
                            lblperdidas.grid(row=4,column=0,columnspan=3)

                            lblPerdidasc = ttk.Label(self.frmuniones,text="Perdidas por conectores = ",font=(18),foreground="DarkOrange2")
                            lblPerdidasc.config(background="mediumturquoise")
                            lblPerdidasc.grid(row=5,column=0)

                            perctxt = Entry(self.frmuniones,font=("consolas",10),textvariable=self.atxconect,justify="center",background="lightsteelblue4",fg="lightcyan")
                            perctxt.grid(row=5,column=1)

                            lblPerdidascdB = ttk.Label(self.frmuniones,text="dB",font=(18),foreground="DarkOrange2")
                            lblPerdidascdB.config(background="mediumturquoise")
                            lblPerdidascdB.grid(row=5,column=2)

                            lblPerdidasunion = ttk.Label(self.frmuniones,text="Coeficiente Atenuacion (𝛼) = ",font=(18),foreground="DarkOrange2")
                            lblPerdidasunion.config(background="mediumturquoise")
                            lblPerdidasunion.grid(row=6,column=0)

                            perutxt = Entry(self.frmuniones,font=("consolas",10),textvariable=self.atxunion,justify="center",background="lightsteelblue4",fg="lightcyan")
                            perutxt.grid(row=6,column=1)


                            if self.vali == 1:
                            #     try:
                            #         btnpV = ttk.Button(self.frmuniones,text="Calcular",style="MyButton.TButton"
                            #             ,command=lambda:self.nucleo(float(pvtxt.get()),indices,float(lambda_fo)*math.pow(10,-9),float(self.num_uniones.get())))
                            #         btnpV.grid(row=7,column=0,columnspan=3)
                            #     except ValueError:
                            #         messagebox.showerror("Fibra Optica","No es valido ese el parametro")
                            #         self.parv.set("")
                            #         self.num_uniones.set("")
                            # else:
                                try:
                                    btnpV = ttk.Button(self.frmuniones,text="Calcular",style="MyButton.TButton"
                                        ,command=lambda:self.guardar_fo(self.tipo,indices,self.parav,self.num_uniones.get(),perctxt.get(),perutxt.get()))
                                    btnpV.grid(row=7,column=0,columnspan=3)                            
                                except ValueError:
                                    messagebox.showerror("Fibra Optica","No es valido ese el parametro")
                                    self.parv.set("")
                                    self.num_uniones.set("")

                    else:
                        messagebox.showwarning("Fibra Optica","Seleccione si va a ingresar el parametro V o no")
                else:
                    messagebox.showwarning("Fibra Optica","Tamaño de fibra optica no valido")
                    self.tamfo1.set("")
            else:
                messagebox.showwarning("Fibra Optica","No es valido ese valor en lambda o tamaño de FO")
                self.lambdafo.set("")
                self.tamfo1.set("")
        except ValueError:
            messagebox.showerror("Fibra Optica","No es valido ese lambda")
        


    def pantalla(self):
        
        self.lambdafo.set("")
        self.parv.set("")
        self.num_uniones.set("")
        self.atxval = []
        self.tipodeatx = []
        self.fibras = []
        self.ftotales = []
        self.seleccion = 1
        self.selmultlat.set(0)
        self.vali = 0
        self.parav = 0
        self.tipo = ""
        self.nucleoFO = 0
        self.val = 1
        self.uniones = 0
        self.unionnum = 1
        
        s = ttk.Style()
        s.configure(
            "MyButtonL.TButton",
            foreground="indianred1",
        )

        s1 = ttk.Style()
        s1.configure(
            "MyButton.TButton",
            foreground="DarkOrange4",
            background="mediumpurple1"
        )
        self.frminicio = Frame()
        self.frminicio.pack()
        self.frminicio.config(width="400",height="230",background="mediumturquoise",bd=10,relief="raised")

        lbl = ttk.Label(self.frminicio,text="Simulador Fibra Optica",font=(52),foreground="DarkOrange2")
        lbl.config(background="mediumturquoise")
        lbl.grid(row=0,column=1,columnspan=4,padx=5,pady=8)

        # lbltipo = ttk.Label(self.frminicio,text="Tipo Fibra",font=(18),foreground="DarkOrange2")
        # lbltipo.config(background="mediumturquoise")
        # lbltipo.grid(row=1,column=0)

        # self.cbmtipo= ttk.Combobox(self.frminicio,state="readondly",values=["Multimodo","Monomodo"])
        # self.cbmtipo.grid(row=1,column=1)
        # self.cbmtipo.current(0)

        lblindices = ttk.Label(self.frminicio,text="Indices",font=(18),foreground="DarkOrange2")
        lblindices.config(background="mediumturquoise")
        lblindices.grid(row=2,column=3)

        self.cbmindices= ttk.Combobox(self.frminicio,state="readondly",
            values=["[1.550 1.520]",
                    "[1.500 1.490]",
                    "[1.485 1.482]",
                    "[1.480 1.460]",
                    "[1.460 1.450]",
                    "[1.533 1.529]",
                    "[1.501 1.499]",
                    "[1.485 1.482]",
                    "[1.483 1.482]",
                    "[1.464 1.450]"
                    ])
        self.cbmindices.grid(row=2,column=4)
        self.cbmindices.current(0)

        lbllambda = ttk.Label(self.frminicio,text="Lambda",font=(18),foreground="DarkOrange2")
        lbllambda.config(background="mediumturquoise")
        lbllambda.grid(row=1,column=0)

        lambdatxt = Entry(self.frminicio,font=("consolas",10),textvariable=self.lambdafo,justify="center",background="lightsteelblue4",fg="lightcyan")
        lambdatxt.grid(row=1,column=1)

        lbllambda = ttk.Label(self.frminicio,text="nm",font=(18),foreground="DarkOrange2")
        lbllambda.config(background="mediumturquoise")
        lbllambda.grid(row=1,column=2)

        lblingradual = ttk.Label(self.frminicio,text="Perfil(𝛼) =",font=(18),foreground="DarkOrange2")
        lblingradual.config(background="mediumturquoise")
        lblingradual.grid(row=1,column=3)

        self.cbmindgr= ttk.Combobox(self.frminicio,state="readondly",values=["Lineal","Parabolico","Escalonado"])
        self.cbmindgr.grid(row=1,column=4)
        self.cbmindgr.current(0)
        
        lblnucleo = ttk.Label(self.frminicio,text="Diametro del Nucleo",font=(18),foreground="DarkOrange2")
        lblnucleo.config(background="mediumturquoise")
        lblnucleo.grid(row=2,column=0)

        self.cbmnucleo = ttk.Combobox(self.frminicio,state="readondly",values=["8","9","50","62.5"])
        self.cbmnucleo.grid(row=2,column=1)
        self.cbmnucleo.current(0)
        
        
        lbllambda = ttk.Label(self.frminicio,text="um",font=(18),foreground="DarkOrange2")
        lbllambda.config(background="mediumturquoise")
        lbllambda.grid(row=2,column=2)

        
        lbltamFO = ttk.Label(self.frminicio,text="longitud Fibra Optica",font=(18),foreground="DarkOrange2")
        lbltamFO.config(background="mediumturquoise")
        lbltamFO.grid(row=3,column=0)

        tamFOtxt = Entry(self.frminicio,font=("consolas",10),textvariable=self.tamfo1,justify="center",background="lightsteelblue4",fg="lightcyan")
        tamFOtxt.grid(row=3,column=1)

        lbltamfokm = ttk.Label(self.frminicio,text="Km",font=(18),foreground="DarkOrange2")
        lbltamfokm.config(background="mediumturquoise")
        lbltamfokm.grid(row=3,column=2)

        btnCalcunion = ttk.Button(self.frminicio,text="Detector",style="MyButton.TButton"
            ,command=lambda:self.fibraoptica_parv(self.cbmindices.get(),lambdatxt.get(),self.cbmnucleo.get(),self.cbmindgr.get(),tamFOtxt.get()))
        btnCalcunion.grid(row=3,column=3,columnspan=2,padx=15)



    def creador(self):
        p = "Curso: Opticas \n"+"Killiam Puentes \n"+"Adriana Alvarez \n"+"Juan Hernandez \n"
        messagebox.showinfo("Autores",p)


class principal_Lineal():
    
    def iniciar(self):
        root = Tk()
        root.title("Simulador Fibra optica")

        barra_menu = Menu(root)
        root.config(menu=barra_menu)
        ayudatr=Menu(barra_menu,tearoff=0)
        ayudatr.add_command(label="Acerca de...",command=lambda:principal_Atenuaciones.creador(self))


        fograf = Menu(barra_menu,tearoff=0)
        fograf.add_command(label="Grafica de fibra",command=lambda:[root.destroy(),fibra_dibujo.principal_Lineal.iniciar(self)])

        fo = Menu(barra_menu,tearoff=0)
        fo.add_command(label="Simulador Fibra Optica",command=lambda:[root.destroy(),principal_Lineal.iniciar(self)])

        barra_menu.add_cascade(label="Simulador Fibra Optica",menu=fo)
        barra_menu.add_cascade(label="Fibra Optica Grafica",menu=fograf)
        barra_menu.add_cascade(label="Ayuda",menu=ayudatr)


        s1 = principal_Atenuaciones(root)
        s1.pantalla()
        root.mainloop()

p = principal_Lineal()
p.iniciar()