# Proyecto-topicos-minimizers-syncmers



Integrantes:

* Vicente Cuello
* Nicolas Pino



Los comandos para ejecutar los códigos son los siguientes:



* minimizer\_sketch.cpp (Minimizer con Countmin):

1. g++ -std=c++17 -O3 -o minimizer\_sketch minimizer\_sketch.cpp
2. .\\minimizer\_sketch



* minimizer\_sketch\_CU.cpp (Minimizer con Countmin CU):

1. g++ -std=c++17 -O3 -o minimizer\_sketch\_CU minimizer\_sketch\_CU.cpp
2. .\\minimizer\_sketch\_CU



* exact\_counter.cpp (Minimizer sin Countmin):

1. g++ -std=c++17 -O3 -o exact\_counter exact\_counter.cpp
2. .\\exact\_counter



* all\_kmers\_sketch.cpp (Método ingenuo con Countmin):

1. g++ -std=c++17 -O3 -o all\_kmers\_sketch all\_kmers\_sketch.cpp
2. .\\all\_kmers\_sketch



* all\_kmers\_sketch\_CU.cpp (Método ingenuo con Countmin CU):

1. g++ -std=c++17 -O3 -o all\_kmers\_sketch\_CU all\_kmers\_sketch\_CU.cpp
2. .\\all\_kmers\_sketch\_CU



* exact\_all\_kmers\_counter.cpp (Método ingenuo sin Countmin):

1. g++ -std=c++17 -O3 -o exact\_all\_kmers\_counter exact\_all\_kmers\_counter.cpp
2. .\\exact\_all\_kmers\_counter



* syncmer\_sketch.cpp (Syncmer con Countmin):

1. g++ -std=c++17 -O3 -o syncmer\_sketch syncmer\_sketch.cpp
2. .\\syncmer\_sketch



* syncmer\_sketch\_CU.cpp (Syncmer con Countmin CU):

1. g++ -std=c++17 -O3 -o syncmer\_sketch\_CU syncmer\_sketch\_CU.cpp
2. .\\syncmer\_sketch\_CU



* exact\_syncmer\_counter.cpp (Syncmer sin Countmin):

1. g++ -std=c++17 -O3 -o exact\_syncmer\_counter exact\_syncmer\_counter.cpp
2. .\\exact\_syncmer\_counter



* compare\_accuracy (código de comparación del top 100 kmers):

1. g++ -std=c++17 -O3 -o compare\_accuracy compare\_accuracy.cpp
2. Aquí depende de que archivos quieres comparar, por ejemplo:

   * .\\compare\_accuracy ../minimizer/minimizer\_sketch\_CU.bin ../minimizer/exact\_countsT.txt



* compare\_accuracy\_global (código de comparación de todos los kmers):

1. g++ -std=c++17 -O3 -o compare\_global compare\_accuracy\_global.cpp
2. Aquí depende de que archivos quieres comparar, por ejemplo:

   * .\\compare\_global ../minimizer/minimizer\_sketch\_CU.bin ../minimizer/exact\_countsT.txt

* query_sketch.cpp (Codigo para ejecutar archivos .bin y buscar kmers que se añadieron al sketch):

1. g++ -std=c++17 -O3 -o query_sketch query_sketch.cpp
2. .\query_sketch   
   

