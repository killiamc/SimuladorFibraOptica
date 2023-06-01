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
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt

import Atenuaciones_Lab6 as uniones

class principal_optica():
    def __init__(self,root):
        self.root = root
        self.angulo = StringVar()
        self.indices = []
        self.seleccion = IntVar()
        pass

    def limpiar(self):
        l = ""
        self.angulo.set(l)
        self.cbmnucleo.current(0)
        self.cbmindices.current(0)
    
    def reflejado(self,n1,n2,angi):
        try:
            angr = math.asin((n1/n2)*math.sin(math.radians(angi)))
            angr = math.degrees(angr)
            return angr
        except ValueError:
            angr = 180-angi
            return angr

    def ecuacion_recta(self,angr,px,py,x):
        m = math.tan(math.radians(angr))
        v = ((m*-px)+py)
        ec = m*x+v
        return ec
    
    def calculos(self,indices,angc):
        anga = math.degrees(math.asin(self.indices[0]*math.cos(math.radians(angc))))
        an = math.sin(math.radians(anga))
        dif = (an**2)/2*(self.indices[0]**2)
        return anga,an,dif

    # def critico(self,indices):
    #     angc = math.asin(float(indices[1])/float(indices[0]))
    #     return angc

    def Fibra_final(self,angi0,nucleo,indices):
        
        try:
            if self.seleccion.get() == 1 or self.seleccion.get() == 2:
                if float(angi0) >= -90 and float(angi0)<=90:
                    self.frminicio.destroy()
                    self.frmresultado = Frame()
                    self.frmresultado.pack()
                    self.frmresultado.config(width="400",height="230",background="mediumturquoise",bd=10,relief="raised")
                    self.indices=[]
                    self.indices.append(float(indices[1:5]))
                    self.indices.append(float(indices[6:10]))
                    
                    limites = []
                    angi0 = float(angi0)
                        
                    

                    x_l1 = [5,150]
                    x_y1 =[5,5]
                    y_y1 = [0,20]
                    y_r1 = [0,0]
                    y_r2 = [20,20]
                    x_a = [0,5]
                    cor = "orange"
                    if nucleo == "9":
                        limites = [9,11]
                        y_l1 = [9,9]
                        y_l2 = [11,11]
                        co = "violet"
                    elif nucleo == "50":
                        limites = [5,15]
                        y_l1 = [5,5]
                        y_l2 = [15,15]
                        co = "skyblue"
                    elif nucleo == "62.5":
                        limites = [3,17]
                        y_l1 = [3,3]
                        y_l2 = [17,17]
                        co = "darkgreen"
                    
                    #parte 1
                    m = math.tan(math.radians(angi0))
                    b = ((m*-5)+10)
                    
                    x_1 = np.linspace(0,5,6)  
                    # x_1 = list(x_1)
                    y_1 = []
                    for j in range(len(x_1)):  
                        p_1 = m*x_1[j]+b
                        y_1.append(p_1)

                    if angi0<0:
                        angi0 = 180+angi0
                    else:
                        angi0 = angi0

                    x = np.linspace(x_l1[0],x_l1[1],x_l1[1]-x_l1[0]+1)
                    y = []
                    angr0 = self.reflejado(1,self.indices[0],float(angi0))
                    angi = 180-angr0
                    angr = self.reflejado(self.indices[0],self.indices[1],float(angi))
                    print(angr)

                    if angi0>=90:
                        angc1 = 90-angr0
                        angr0=(180-angr0)
                    else:
                        angc1 = 90-angr0
                    # x = list(x)
                    #parte 
                    angr_r = angr0
                    angc = math.asin(float(self.indices[1])/float(self.indices[0]))
                    angc = math.degrees(angc)
                    valor = True
                    x_re = []
                    y_re = []
                    con = 0
                    if angi0 != 90:
                        px = 5
                        py = 10
                        for i in range(len(x)):
                            m = math.tan(math.radians(angr0))
                            b = ((m*-px)+py)
                            punto = m*x[i]+b
                            if punto <= limites[0] or punto >= limites[1]:
                                if punto<=limites[0]:
                                    y_l2 = [punto, punto]
                                elif punto >= limites[1]:
                                    y_l1 = [punto,punto]
                                angr0 = (180-angr0)
                                # angr = math.degrees(angr)
                                px = x[i] 
                                py = punto 
                                if angc1<angc:#refraccion perdida
                                    m1 = math.tan(math.radians(angr_r))
                                    b1 = ((m1*-5)+10)
                                    punto1 = m1*x[i]+b1
                                    if punto1>=y_r2[0]:
                                        if con < 1:    
                                            x_re.append(x[i])
                                            y_re.append(y_r2[0])
                                        con+=1
                                    elif punto1 <= y_r1[0]:
                                        if con < 1:
                                            x_re.append(x[i])
                                            y_re.append(y_r1[0])
                                        con+=1
                                    else:    
                                        x_re.append(x[i])
                                        y_re.append(punto1)
                                    valor = False
                            y.append(punto)
    
                        x = np.concatenate((x_1,x),axis=0)
                        x = list(x)
                        y = np.concatenate((y_1,y),axis=0)
                        y = list(y)
                    else:
                        x = [5,35]
                        y = [10,10]

                    lbl = ttk.Label(self.frmresultado,text="Fibra Optica Dibujo",font=(52),foreground="DarkOrange2")
                    lbl.config(background="mediumturquoise")
                    lbl.grid(row=0,column=0,columnspan=4,padx=5,pady=8)

                    anga,an,dif = self.calculos(indices,angc)
                    fig = Figure(figsize=(3,3),dpi=100)
                    if self.seleccion.get() == 1:
                        m_a = math.tan(math.radians(anga))
                        ba = ((m_a*-5)+10)
                        p = m_a*x_a[0]+ba
                        p2 = m_a*x_a[1]+ba
                        y_a1 =[p,p2]
                        m_a=math.tan(math.radians((180-anga)))
                        ba = ((m_a*-5)+10)
                        p = m_a*x_a[0]+ba
                        p2 = m_a*x_a[1]+ba
                        y_a2 =[p,p2]
                        if valor:   
                            fig.add_subplot(111).plot(x,y,"lightgreen",x_l1,y_l1,co,x_l1,y_l2,co,x_l1,y_r1,cor,x_l1,y_r2,cor,x_y1,y_y1,cor,x_a,y_a1,'darkblue',x_a,y_a2,'darkblue')
                        else:
                            fig.add_subplot(111).plot(x,y,"lightgreen",x_l1,y_l1,co,x_l1,y_l2,co,x_l1,y_r1,cor,x_l1,y_r2,cor,x_y1,y_y1,cor,x_re,y_re,'r',x_a,y_a1,'darkblue',x_a,y_a2,'darkblue')
                    elif self.seleccion.get() == 2:
                        if valor:   
                            fig.add_subplot(111).plot(x,y,"lightgreen",x_l1,y_l1,co,x_l1,y_l2,co,x_l1,y_r1,cor,x_l1,y_r2,cor,x_y1,y_y1,cor)
                        else:
                            fig.add_subplot(111).plot(x,y,"lightgreen",x_l1,y_l1,co,x_l1,y_l2,co,x_l1,y_r1,cor,x_l1,y_r2,cor,x_y1,y_y1,cor,x_re,y_re,'r')

                    canvas = FigureCanvasTkAgg(fig,master=self.frmresultado)
                    canvas.draw()
                    canvas.get_tk_widget().grid(row=1,column=0,columnspan=4)

                    v = "N.A. = "+str(round(an,2))
                    lbl = ttk.Label(self.frmresultado,text=v,font=(10),foreground="DarkOrange2")
                    lbl.config(background="mediumturquoise")
                    lbl.grid(row=2,column=1)

                    v = "angulo a = "+str(round(anga,2))+"°"
                    lbla = ttk.Label(self.frmresultado,text=v,font=(10),foreground="DarkOrange2")
                    lbla.config(background="mediumturquoise")
                    lbla.grid(row=2,column=0)

                    v = "diferencia = "+str(round(dif,2))
                    lbld = ttk.Label(self.frmresultado,text=v,font=(10),foreground="DarkOrange2")
                    lbld.config(background="mediumturquoise")
                    lbld.grid(row=3,column=0)

                    v = "angulo c = "+str(round(angc,2))+"°"
                    lbla = ttk.Label(self.frmresultado,text=v,font=(10),foreground="DarkOrange2")
                    lbla.config(background="mediumturquoise")
                    lbla.grid(row=3,column=1)

                    btnmodifi = ttk.Button(self.frmresultado,text="Volver",style="MyButton.TButton",
                        command=lambda:[self.frmresultado.destroy(),self.pantalla()])
                    btnmodifi.grid(row=2,column=2,rowspan=2)

                    if valor:
                        messagebox.showinfo("Fibra Optica","No se presenta perdidas por refracción")
                    else:
                        messagebox.showinfo("Fibra Optica","se presenta perdidas refracción")
                else:
                    messagebox.showwarning("Fibra Optica","Solo angulos entre -90 a 90")
                    self.angulo.set(" ")
            else:
                messagebox.showwarning("Fibra Optica","Seleccione si quiere apertura numerica o no")
        except ValueError:
            messagebox.showwarning("Fibra Optica","Ingrese solo numero")
            self.angulo.set(" ")



    def pantalla(self):
        
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

        lbl = ttk.Label(self.frminicio,text="Fibra Optica Dibujo",font=(52),foreground="DarkOrange2")
        lbl.config(background="mediumturquoise")
        lbl.grid(row=0,column=0,columnspan=4,padx=5,pady=8)

        lblang = ttk.Label(self.frminicio,text="angulo",font=(18),foreground="DarkOrange2")
        lblang.config(background="mediumturquoise")
        lblang.grid(row=1,column=0,columnspan=2)

        angitxt = Entry(self.frminicio,font=("consolas",10),textvariable=self.angulo,justify="center",background="lightsteelblue4",fg="lightcyan")
        angitxt.grid(row=1,column=2,columnspan=3)

        lblnucleo = ttk.Label(self.frminicio,text="nucleo",font=(18),foreground="DarkOrange2")
        lblnucleo.config(background="mediumturquoise")
        lblnucleo.grid(row=2,column=0)


        self.cbmnucleo = ttk.Combobox(self.frminicio,state="readondly",values=["9","50","62.5"])
        self.cbmnucleo.grid(row=2,column=1)
        self.cbmnucleo.current(0)

        rbsiaper = Radiobutton(self.frminicio,background="mediumturquoise",foreground="DarkOrange2",text="SI Apertura numerica",variable=self.seleccion,value=1)
        rbsiaper.grid(row=2,column=2)

        lblindices = ttk.Label(self.frminicio,text="indices",font=(18),foreground="DarkOrange2")
        lblindices.config(background="mediumturquoise")
        lblindices.grid(row=3,column=0)
        
        self.cbmindices= ttk.Combobox(self.frminicio,state="readondly",values=["[1.50 1.41]","[1.49 1.46]","[1.55 1.45]","[1.46 1.35]","[1.50 1.43]"])
        self.cbmindices.grid(row=3,column=1)
        self.cbmindices.current(0)

        rbsiaper = Radiobutton(self.frminicio,background="mediumturquoise",foreground="DarkOrange2",text="NO Apertura numerica",variable=self.seleccion,value=2)
        rbsiaper.grid(row=3,column=2)

        btningresaprob = ttk.Button(self.frminicio,text="Calcular",style="MyButton.TButton"
            ,command=lambda:self.Fibra_final(angitxt.get(),self.cbmnucleo.get(),self.cbmindices.get()))
        btningresaprob.grid(row=2,column=3,columnspan=2)

        btncallimpiar = ttk.Button(self.frminicio,text="Limpiar",style="MyButton.TButton",
            command=lambda:self.limpiar())
        btncallimpiar.grid(row=3,column=3,columnspan=2)

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
        ayudatr.add_command(label="Autores",command=lambda:principal_optica.creador(self))

        fo = Menu(barra_menu,tearoff=0)
        fo.add_command(label="Uniones de Fibra Optica",command=lambda:[root.destroy(),uniones.principal_Lineal.iniciar(self)])

        fograf = Menu(barra_menu,tearoff=0)
        fograf.add_command(label="Grafica de fibra",command=lambda:[root.destroy(),principal_Lineal.iniciar(self)])

        barra_menu.add_cascade(label="Simulador Fibra Optica",menu=fo)
        barra_menu.add_cascade(label="Fibra Optica Grafica",menu=fograf)
        barra_menu.add_cascade(label="Ayuda",menu=ayudatr)

        s1 = principal_optica(root)
        s1.pantalla()
        root.mainloop()

# p = principal_Lineal()
# p.iniciar()