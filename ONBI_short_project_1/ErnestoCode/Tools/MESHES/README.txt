################################################################
#INSTALACION
################################################################

################################
#Sistemas Linux:
################################
#Paso 1:
#En sistemas Linux, la instalacion de un Matlab implica que en los directorios de instalacion de Matlab se copian una serie de librerias "standard" de c que, en el entorno de Matlab, se utilizan con preferencia sobre las del sistema. 
#La libreria mas critica en este momento para nosotros es la libreria libstdc++.so.6. Matlab instala su version de esta libreria en los directorios &MATLAB_INS/sys/os/glnxa64 y &MATLAB_INS/bin/glnxa64 (habitualmente suelen ser equivalente a /usr/local/MATLAB/R2012a/sys/os/glnxa64/ y /usr/local/MATLAB/R2012a/bin/glnxa64/). 
#Si hemos instalado una version de compilador mas nueva o diferente de la que Matlab da soporte, probablemente, a la hora de compilar y/o ejecutar compilaciones de archivos mex nos de algun error similar a este al intentar ejecutar el mex vtkPolyDataReader.mexa64:
#Invalid MEX-file '/home/mia/mmonla/libs/MESHES_VTK/vtkPolyDataReader.mexa64': /usr/local/MATLAB/R2012a/bin/glnxa64/libstdc++.so.6: version `GLIBCXX_3.4.15' not found (required by /usr/local/lib/vtk-5.10/libvtkFiltering.so.5.10).
#Esto se debe a que la version de compilador utilizada para generar el mex utiliza una version diferente de libstdc++.so.6 que la que Matlab encuentra instalada en sus directorios.
#Para solucionar este problema, hacemos una copia de los ficheros libstdc++.so.6 que encontramos en &MATLAB_INS/sys/os/glnxa64 y &MATLAB_INS/bin/glnxa64 y en ambos lugares hacemos un link dinamico a la misma libreria del sistema (ln -s /usr/lib/libstdc++.so.6 libstdc++.so.6)

#Paso 2:
#muchos de nuestros mex van a utilizar librerias dinamicas durante su ejecucion, es decir, que nuestros archivos mex, dependen de librerias de terceros del tipo .so. Estas libreria, Matlab las busca y las carga en tiempo de ejecucion desde los lugares indicados en la variable LD_LIBRARY_PATH.

#Por lo tanto, antes de poder utilizar las bibliotecas .mex como funciones de Matlab, debe modificarse esa variable de entorno LD_LIBRARY_PATH. Esta variable de entorno indica a Matlab en tiempo de ejecucion los directorios donde se pueden localizar todas las librerias dinamicas que nuestro mex utilizara. Para modificar variables dentro de Matlab, se realiza con la instruccion setenv. Si no se utiliza una forma permanente de realizar esta modificacion, unicamente funcionara en la sesion actual de Matlab. 
#Para evitar este inconveniente, deberia ser suficiente con realizar los siguientes pasos (pero no funciona):
#1. Editar el fichero de Matlab startup.m que habitualmente se encuentra en /home/mia/%user%/matlab/ o en /home/mia/%user%/.matlab/%matlab_version%/ o algo similar.
#2. Agregar la instruccion: setenv('LD_LIBRARY_PATH', [getenv('LD_LIBRARY_PATH'),':dir_lib']), donde dir_lib es el directorio donde se guardaran las librerias dinamicas que utilizen nuestros mex.

#Por lo tanto lo que debe hacerse es modificar esta variable de entorno en el .profile de cada usuario:
#1. Editar el fichero .bashrc o .bash_profile y agregar la ruta del dir_lib a la variable LD_LIBRARY_PATH. Debe resultar una linea similar a esta : 'LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/dir_lib/'
#En ese mismo fichero, a continuacion, agregar si no existe la instruccion: export LD_LIBRARY_PATH.
#El directorio que agregaremos tanto para compilaciones como para ejecuciones es el directorio llamado 'vtk_libs'. Se encuentra en un sitio diferente si estamos en una maquina local o en el cluster hermes. En local la ruta completa a agregar es '/extra/disco1/miaTools/MESHES_VTK/vtk_libs' mientras que para un nodo de hermes la ruta debe ser '/home/mia/olmos/miaTools/MESHES_VTK/vtk_libs'. 

#Para poder saber las librerias dinamicas de las que depende cada mex, se puede usar desde dentro de Matlab la instruccion: !ldd name_file.mexa64


################################################################
#COMPILACION de NUEVOS MEX
################################################################

#Habitualmente la forma de trabajar para conseguir nuevos compilados sera utilizar nuestro Matlab local para compilar mediante el wrapper de Matlab 'mex'.
#Si todo esta bien configurado, tal como se indica en el apartado de configuracion anterior, esta instruccion 'mex' llamara a nuestro compilador 'gcc' (actualmente version 4.6).

#Durante una compilacion de un '.mex' que utiliza librerias VTK, el compilador necesitara de los archivos header en compilacion y las librerias en linkado. En el apartado anterior se ha indicado como y donde residen las librerias VTK compiladas en el cluster que seran utilizadas para la compilacion y la ejecucion. Pero vamos a repasar los pasos de configuracion necesarios para compilar nuevos programas.

#Paso 1: Configuracion de variables de entorno de librerias dinamicas.
#1. Editar el fichero .bashrc o .bash_profile.
#2. Modificar la variable de entorno de carga de librerias con la ruta correspondiente. Pueden darse varios casos:
#a. Cluster Hermes. En hermes la variable a modificar es la 'DYLD_LIBRARY_PATH'
#b. Equipo local. En local la variable a modificar es la 'LD_LIBRARY_PATH'
#c. Si existe ya la variable en el fichero, solo la modificaremos agregando una linea como esta: 
#	Cluster:	DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/home/mia/olmos/miaTools/MESHES_VTK/vtk_libs
#	Local:		LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/extra/disco1/miaTools/MESHES_VTK/vtk_libs
#d. Si no existe las instrucciones a escribir seran:
#	Cluster:	DYLD_LIBRARY_PATH=/home/mia/olmos/miaTools/MESHES_VTK/vtk_libs
#	Local:		LD_LIBRARY_PATH=/extra/disco1/miaTools/MESHES_VTK/vtk_libs
#e. nos aseguramos de exportar esta variable revisando que exista o agregando la instruccion:
#	Cluster:	export DYLD_LIBRARY_PATH
#	Local:		export LD_LIBRARY_PATH

#Paso 2: Configuracion de nuestro entorno Matlab para la compilacion
#1. Editar el fichero de configuracion de Matlab 'mexopts.sh' donde se encuentran las opciones que el compilador 'mex' utiliza. Para saber donde se encuentra este fichero, ejecutamos en Matlab 'prefdir'. Si este fichero no existe, se reconstruye ejecutando en Matlab 'mex -setup' y se elige el compilador a utilizar, en nuestro caso gcc y g++ ambos en version 4.6.
#2. Modificar las variables de rutas de inclusion, linkado y librerias que usaran nuestros nuevos '.mex'. Para indicarle las donde debe buscar los ficheros header, tenemos que modificar las variables CFLAGS, CXXFLAGS, CLIBS y CXXLIBS de nuestro sistema, que en este caso es el glna64, con algo similar a esto:
# CXXFLAGS="$CXXFLAGS $INCLUDE_DIRS_VTK", CFLAGS="$CFLAGS $INCLUDE_DIRS_VTK", CXXLIBS="$CXXLIBS $LIBS_VTK", CLIBS="$CLIBS $LIBS_VTK", donde las variables INCLUDE_DIRS_VTK y LIBS_VTK tienen que estar definidas antes y reflejar las rutas de los headers y libs de la VTK. 
#En nuestro caso son asi para el equipo local:
#	INCLUDE_DIRS_VTK=" -I/usr/local/include/vtk-5.10 "
#	LIBS_VTK=" -L/extra/disco1/miaTools/MESHES_VTK/vtk_libs -lvtkCommon -lvtkDICOMParser -lvtkexoIIc -lvtkexpat -lvtkFiltering -lvtkfreetype -lvtkftgl -lvtkGenericFiltering -lvtkGraphics -lvtkHybrid -lvtkImaging -lvtkIO -lvtkjpeg -lvtkNetCDF -lvtkpng -lvtkRendering -lvtksys -lvtktiff -lvtkVolumeRendering -lvtkWidgets -lvtkzlib "

#y asi para el cluster hermes:
#	INCLUDE_DIRS_VTK=" -I/extra/database/mia/soft/vtk_5.10.1/include/vtk-5.10/ "
#	LIBS_VTK=" -L/home/mia/olmos/miaTools/MESHES_VTK/vtk_libs -lvtkCommon -lvtkDICOMParser -lvtkexoIIc -lvtkexpat -lvtkFiltering -lvtkfreetype -lvtkftgl -lvtkGenericFiltering -lvtkGraphics -lvtkHybrid -lvtkImaging -lvtkIO -lvtkjpeg -lvtkNetCDF -lvtkpng -lvtkRendering -lvtksys -lvtktiff -lvtkVolumeRendering -lvtkWidgets -lvtkzlib "

#Observar que el directorio de los header si que se sigue manteniendo direccionado a los originales instalados al compilar la libreria VTK.



