En este directorio se muestran códigos que aplican el algoritmo de Gillespie para obtener diferentes parámetros de simulación.

Los archivos de código que hay se corresponden a las siguientes simulaciones:

1. *Gillespie.cpp:* Simulación en la que simplemente se saca una imagen de la cascada.
2. *Gillespie_beta:* (Hilos) Calcula la densidad estacionaria en función del parámetro de control para poder obtener beta.
3. *GillespieCorrLengthsSub.cpp:* Calcula las longitudes de correlación en función del parámetro de control por debajo del punto crítico partiendo de un estado con una única celda infectada en el centro.
4. *GillespieCorrLengthsSubHom.cpp:* Mismo que el programa anterior pero solo para la longitud de correlación temporal y partiendo del estado con todas las celdas activas.
5. *PBC_Gillespie_metodo3_1Isla.cpp:* (Hilos) Intenta calcular las longitudes de correlación por encima del punto crítico usando condiciones de contorno periódicas. Para ello se parte del estado estacionario y se obtiene la anchura y la altura de la isla inactiva más grande de la simulación. No da los coeficientes esperados.
6. *GillespieRelCorrLengthsSup.cpp:* Calcula la relación entre ambas longitudes de correlación por encima del punto crítico. Para ello se parte del estado con una única celda activa en el centro y se deja avanzar la simulación. La relación entre ambas longitudes vendrá de la pendiente. Se almacenan varios ficheros en el directorio Figuras para después ajustarlos usando "get_pendientes.py".