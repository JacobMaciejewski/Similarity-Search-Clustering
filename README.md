# Similarity-Search-Clustering
Efficient top-K Vector Similarity Search Queries on MNIST Images Dataset using custom Locality Sensitive Hashing 🔢

Suppose we have two sets of entities represented by vectors (images in our case), namely a `query-set` and a `train-set`. We define our problem
as the task of finding for each image of the query-set, the $k$ most similar images in the train-set. One can easily infer that such an endeavour becomes
prohibitively slow, when we increase the number of query or train entities. For example, if the total number of images in both sets is $N$ and both sets
are equal in size, we have to deal with $O(N^2)$ complexity.

## Algorithms:

In order to produce results in a time efficient manner, we have to employ probabilistic methods that group similar items (`clustering`) and allow for
comparisons within smaller groups of similar items with high probability. We implement the following novel methods: </br>

*   [Locality-Sensitive-Hashing](https://en.wikipedia.org/wiki/Locality-sensitive_hashing) - random vectors following normal distribution are chosen
from the hyperspace defined by the images. These constitute borders for subspaces of the initial space. Those areas are populated by the train images.
For each query image, its target subspace is infered and it is compared to its neighbours within the area, resulting in similarity search with a smaller
part of the initial train set.

*   [Hypercube-Random-Projection](https://eclass.uoa.gr/modules/document/file.php/DI460/Lectures/2c.rproj.pdf) - a novel method using an LSH family of
functions for [Randomized-Projections](https://arxiv.org/pdf/1612.07405.pdf) of vectors' points to a [Hamming Cube](https://arxiv.org/abs/2005.09205)
of $log(n)$ dimensionality, where $n$ is the number of pixels of the images. Examines points which are assigned to the same or nearby vertices on the Hamming cube.

*   [Clustering](https://en.wikipedia.org/wiki/Cluster_analysis) - in the simple version the `centroids` are chosen in the images' hyperplane
using a probabilistic algorithm that maximizes the chance for their even distribution. Images are distributed to the corresponding clusters using
[Lloyd's Algorithm](https://en.wikipedia.org/wiki/Lloyd%27s_algorithm), then the centroids are updated using [K-Medians](https://en.wikipedia.org/wiki/K-medians_clustering).
We develop two novel *clustering* methods: </br></br>
&emsp;&emsp; 1. [LSH Range Searching](https://en.wikipedia.org/wiki/Range_searching) - emulates the iterative increase of the range within we search
for neighbours of the centroids, in each 
</br>&emsp;&emsp; iteration points are assigned to the nearest centroids</br></br>
&emsp;&emsp; 2. [Hypercube Range Searching](https://en.wikipedia.org/wiki/Range_searching) - similar to the upper approach but with the usage of Hamming
Hypercubes</br>

The clusters are being updated till centroids' deviation from previous point gets into a specified range. The quality of the produced 
clusters is examined through [Silhouette Evaluation](https://en.wikipedia.org/wiki/Silhouette_(clustering)).


## Dataset:

We use the famous [MNIST Dataset](https://en.wikipedia.org/wiki/MNIST_database). It consists of 60.000 images, consisting of 784 Pixels represented with a number in range [0,255].
*   Training Dataset : **train-images-idx3-ubyte** - contains all the images of the *MNIST Dataset*
*   Query Dataset : **t10k-images-idx3-ubyte** - a subset of the initial dataset with 10.000 images in total
 
## Similarity:

Once we have identified the candidates in the training set, we have to infer how similar they are to the query image. We do that by comparing the two
vectors using a specific similarity metric. In the case of *LSH* we use the [Manhattan Distance](https://en.wikipedia.org/wiki/Taxicab_geometry) and 
[Hamming Distance](https://en.wikipedia.org/wiki/Hamming_distance) for Hypercube. When we execute a *Range Searching* approach, we still use those two
metrics to infer similarity within the given range.

## Compilation & Execution:

In order to compile the code, you simply open the command line and enter: `make`. The compilation produces 3 different executable names with distinct
names, each one for the methods we described above. You can execute the implemented algorithms, as follows:

* *LSH*: 

`./lsh -d <input file> -q <query file> -k <number of hastables> -L <number of projection functions> -ο <output file>
-Ν <number of nearest neighbours> -R <radius>` 

* *Hypercube*:

`./cube -d <input file> -q <query file> -k <hypercube dimensions> -M <maximum check points> -probes <number of probes> -ο <output file> -Ν <number of nearest neighbours> -R <radius>`

* *Clustering*:

`./cluster -i  <input file> -c <configuration file> -o <output file> -complete <optional> -m <clustering method: Classic / LSH / HyperCube>` 

## Further informations

We include a **readME** PDF file that contains an extensive description of the project set up, the inner workings of the respective algorithms and
the design choices of the creators, namely us. Also, we include a description of the contents of each source file. This PDF is written in *Greek*!

Big thanks to my fellow collaborator: [Anastasios Melidonis](https://github.com/Anastasios084)

*Built as part of the course: Software Development for Algorithmic Problems , Winter of 2020. University of Athens, DiT.*
