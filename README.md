DBSCAN: Density-based Spatial Clustering of Applications with Noise

- A cluster is defined as a maximal set of densely connected points, which indicates that clusters are the dense regions of a data space, separated by regions of lower density.

- DBScan makes use of two important parameters:
•	Eps: radius for the neighborhood of a point P.
•	MinPts: minimum number of points in the given neighborhood.

- Optimum values of Eps and MinPts will yield good clustering results. These values can be changed in DBScan.java 

- Also, ValidateCoefficients.java is used to calculate the Jaccard and Rand Coefficient. Silhouette.java is used to calculate the Silhouette coefficient.
