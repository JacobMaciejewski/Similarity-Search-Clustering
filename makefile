#some name types for readable make file
OB = lib/mainCluster.o lib/mainHyperCube.o lib/mainLSH.o lib/mainAux.o lib/hashTable.o lib/image.o lib/bucket.o lib/LSH.o lib/imageDuplet.o lib/algoAux.o lib/hyperCube.o lib/cluster.o lib/clustering.o lib/clusterEntry.o
OBLSH = lib/mainLSH.o lib/mainAux.o lib/hashTable.o lib/image.o lib/bucket.o lib/LSH.o lib/imageDuplet.o lib/algoAux.o lib/hyperCube.o lib/cluster.o lib/clustering.o lib/clusterEntry.o
OBHC = lib/mainHyperCube.o lib/mainAux.o lib/hashTable.o lib/image.o lib/bucket.o lib/LSH.o lib/imageDuplet.o lib/algoAux.o lib/hyperCube.o lib/cluster.o lib/clustering.o lib/clusterEntry.o
OBCLUST = lib/mainCluster.o lib/mainAux.o lib/hashTable.o lib/image.o lib/bucket.o lib/LSH.o lib/imageDuplet.o lib/algoAux.o lib/hyperCube.o lib/cluster.o lib/clustering.o lib/clusterEntry.o
HEADERS = headers/mainAux.h headers/hashTable.h headers/image.h headers/bucket.h headers/LSH.h headers/imageDuplet.h headers/algoAux.h headers/hyperCube.h headers/cluster.h headers/clustering.h headers/clusterEntry.h
EXECLSH = lsh
EXECHC = cube
EXECCLUST = cluster
COMP = g++
FLAGS = -g -Wall -c -ggdb3
#executable
all: lsh cube cluster

lib/hashTable.o: source/hashTable.cpp
	$(COMP) $(FLAGS) source/hashTable.cpp
	mv hashTable.o lib/hashTable.o

lib/mainHyperCube.o: source/mainHyperCube.cpp $(HEADERS)
	$(COMP) $(FLAGS) source/mainHyperCube.cpp 
	mv mainHyperCube.o lib/mainHyperCube.o

lib/mainLSH.o: source/mainLSH.cpp $(HEADERS)
	$(COMP) $(FLAGS) source/mainLSH.cpp 
	mv mainLSH.o lib/mainLSH.o

lib/mainCluster.o: source/mainCluster.cpp $(HEADERS)
	$(COMP) $(FLAGS) source/mainCluster.cpp 
	mv mainCluster.o lib/mainCluster.o

lib/mainAux.o: source/mainAux.cpp headers/mainAux.h 
	$(COMP) $(FLAGS) source/mainAux.cpp
	mv mainAux.o lib/mainAux.o

lib/image.o: source/image.cpp headers/image.h 
	$(COMP) $(FLAGS) source/image.cpp
	mv image.o lib/image.o

lib/bucket.o: source/bucket.cpp headers/bucket.h 
	$(COMP) $(FLAGS) source/bucket.cpp
	mv bucket.o lib/bucket.o

lib/LSH.o: source/LSH.cpp headers/LSH.h 
	$(COMP) $(FLAGS) source/LSH.cpp
	mv LSH.o lib/LSH.o

lib/imageDuplet.o: source/imageDuplet.cpp headers/imageDuplet.h 
	$(COMP) $(FLAGS) source/imageDuplet.cpp
	mv imageDuplet.o lib/imageDuplet.o

lib/algoAux.o: source/algoAux.cpp headers/algoAux.h 
	$(COMP) $(FLAGS) source/algoAux.cpp
	mv algoAux.o lib/algoAux.o

lib/hyperCube.o: source/hyperCube.cpp headers/hyperCube.h 
	$(COMP) $(FLAGS) source/hyperCube.cpp
	mv hyperCube.o lib/hyperCube.o

lib/cluster.o: source/cluster.cpp headers/cluster.h 
	$(COMP) $(FLAGS) source/cluster.cpp
	mv cluster.o lib/cluster.o

lib/clustering.o: source/clustering.cpp headers/clustering.h 
	$(COMP) $(FLAGS) source/clustering.cpp
	mv clustering.o lib/clustering.o

lib/clusterEntry.o: source/clusterEntry.cpp headers/clusterEntry.h 
	$(COMP) $(FLAGS) source/clusterEntry.cpp
	mv clusterEntry.o lib/clusterEntry.o

lsh: $(OBLSH)
	$(COMP) -g $(OBLSH) -o $(EXECLSH)

cube: $(OBHC)
	$(COMP) -g $(OBHC) -o $(EXECHC)

cluster: $(OBCLUST)
	$(COMP) -g $(OBCLUST) -o $(EXECCLUST)
#  cleaning command
clean :
	rm -f $(OB) $(EXECLSH) $(EXECHC) $(EXECCLUST)