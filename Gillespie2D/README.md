En este directorio se muestran códigos que aplican el algoritmo de Gillespie para obtener diferentes parámetros de simulación en dos dimensiones.

Los archivos de código que hay se corresponden a las siguientes simulaciones:

1. **Gillespie2D.cpp:** Simulación en la que simplemente se saca una imagen de la cascada.
2. **Gillespie2D_beta:** (*Hilos*) Calcula la densidad estacionaria en función del parámetro de control para poder obtener beta.
3. **Gillespie2D_metodo1.cpp:** (*Hilos*) Calcula las longitudes de correlación en función del parámetro de control por debajo del punto crítico partiendo de un estado con una única celda infectada en el centro.
4. **Gillespie2D_metodo2.cpp:** (*Hilos*) Mismo que el programa anterior pero solo para la longitud de correlación temporal y partiendo del estado con todas las celdas activas.
5. **Gillespie2D_metodo4.cpp:** (*Hilos*) Calcula la relación entre ambas longitudes de correlación por encima del punto crítico. Para ello se parte del estado con una única celda activa en el centro y se deja avanzar la simulación. La relación entre ambas longitudes vendrá de la pendiente. Se almacenan varios ficheros en el directorio Figuras para después ajustarlos usando "*get_pendientes.py*".
6. **get_pendientes:** Script en Python en el que tenemos que ajustar la variable DIRECTORY (directorio con todos los archivos generados por "*GillespieRelCorrLengthsSup.cpp*") y podemos ajustar también la variable OUTPUT para indicar donde queremos que se guarden los resultados.

Los archivos en los que viene indicado entre paréntesis la palabra *Hilos* hacen uso de hilos para agilizar los cálculos y tendrán que ser compilados con las etiqueta -pthread y con el estándar de C++ de 2020 (-std=c++20).
