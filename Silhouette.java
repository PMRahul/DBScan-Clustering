import java.util.*;

public class Silhouette
{
	double[][] matrix;
	double silhouetteCoefficient = 0;
	
	public Silhouette(double[][] distanceMatrix)
	{
		matrix = distanceMatrix;
	}
	
	public void calculate(List<List<Integer>> clusterList)
	{
		int clusterCount = 0;
		double silh_eachcluster = 0.0;
		double silh_clustering = 0.0;
		
		for(List<Integer> cluster : clusterList)
		{
			for(int geneId : cluster)
			{
				double inClusterAvg = average_Incluster(geneId, cluster);
				double outClusterAvg = average_Outcluster(geneId, clusterList, cluster);

				if(inClusterAvg!=0 && outClusterAvg!=0)
				{
					silhouetteCoefficient = 1 - (inClusterAvg/outClusterAvg);
				}
				silh_eachcluster += silhouetteCoefficient;
				
			}
			silh_eachcluster = silh_eachcluster/cluster.size();
			clusterCount++;
			System.out.println("Silhouette coefficient for cluster " + clusterCount + ": " + silh_eachcluster);
			silh_clustering += silh_eachcluster;
		}
		System.out.println("Silhouette coefficient of clustering: " + silh_clustering/clusterList.size());
	}
	
	public double average_Incluster(int geneId, List<Integer> cluster)
	{
		double inClusterAvg = 0;
		
		for(int point : cluster)
		{
			inClusterAvg += matrix[geneId-1][point-1];
		}
		
		inClusterAvg = inClusterAvg/(cluster.size());
		return inClusterAvg;
	}
	
	public double average_Outcluster(int geneId,List<List<Integer>> clusterList,List<Integer> cluster)
	{
		List<List<Integer>> otherClusters=new ArrayList<>();
		otherClusters.addAll(clusterList);
		otherClusters.remove(cluster);
		
		double minAvg=100;
	
		for(List<Integer> currentCluster: otherClusters)
		{
			double clusterAvg = 0;
			for(int point: currentCluster)
			{
				clusterAvg += matrix[geneId-1][point-1];
				
			}
			clusterAvg = clusterAvg/currentCluster.size();
			if(clusterAvg<minAvg)
				minAvg = clusterAvg;
		}
		return minAvg;
	}

}
