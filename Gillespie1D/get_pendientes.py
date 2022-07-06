#!/usr/bin/python

import os
import numpy as np

DIRECTORY = "/home/azzorini/Fisica/SistemasComplejos/Trabajo/Pendientes"
OUTPUT = "landa_vs_pendientes_new.txt"

landas_pendientes = {}

# Hacemos un ajuste para todos los ficheros txt en DIRECTORY
for file in os.listdir(DIRECTORY):
    if ".txt" in file:
        landa, n, lado = file.split("_")
        lado = lado.replace(".txt", "")

        x, y = np.loadtxt(os.path.join(DIRECTORY, file)).T
        m, n = np.polyfit(x, y, 1)

        pendiente = m if lado == "R" else -m 

        if landa in landas_pendientes:
            landas_pendientes[landa].append(pendiente)
        else:
            landas_pendientes[landa] = [pendiente]

# Guardamos todos los datos en un fichero de texto
with open(OUTPUT, "w") as f:
    for (landa, pendientes) in landas_pendientes.items():
        f.write(landa)
        for m in pendientes:
            f.write(f"\t{m}")
        f.write('\n')
