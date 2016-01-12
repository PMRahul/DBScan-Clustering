import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.List;
import java.util.ArrayList;

public class ValidateCoefficients 
{
	public double JaccardCoefficient(Map<Integer, Integer> clusterMap, HashMap<Integer, Integer> groundTruth) throws IOException
	{
		double externalIndex = 0;
		
		int clusterSize = clusterMap.size();
		int groundTruthSize = groundTruth.size();
		
		int[][] clusterArray = new int[clusterSize][clusterSize];
		int[][] groundTruthArray = new int[groundTruthSize][groundTruthSize];
		
		for(int i = 0; i < clusterSize; i++)
		{
			for(int j = 0; j < clusterSize; j++)
			{
				if(i == j)
				{
					clusterArray[i][j] = 1;
					groundTruthArray[i][j] = 1;
				}
				else
				{
					clusterArray[i][j] = 0;
					groundTruthArray[i][j] = 0;
				}
				
				int key1 = i+1;
				int key2 = j+1;
				
				if(clusterMap.get(key1) == clusterMap.get(key2))
				{
					clusterArray[i][j] = 1;
					clusterArray[j][i] = 1;
				}
				
				if(groundTruth.get(key1) == groundTruth.get(key2))
				{
					groundTruthArray[i][j] = 1;
					groundTruthArray[j][i] = 1;
				}
			}
		}
		
		FileWriter fstream = new FileWriter("Matrix.txt");
        BufferedWriter out = new BufferedWriter(fstream);
		//PRINT BOTH ARRAYS AND CHECK
		for(int i = 0; i < clusterSize; i++)
		{
			for(int j = 0; j < clusterSize; j++)
			{
				out.write(clusterArray[i][j]+"\t");
			}
			out.write("\n");
		}
		
		out.close();
		
		int sameClusterCount = 0;
		int same_diff_ClusterCount = 0;
		int diff_same_ClusterCount = 0;
		int diffClusterCount = 0;
		
		for(int i = 0; i < clusterSize; i++)
		{
			for(int j = 0; j < clusterSize; j++)
			{
				if(clusterArray[i][j] == 1 && groundTruthArray[i][j] == 1)
					sameClusterCount += 1;
				else if(clusterArray[i][j] == 1 && groundTruthArray[i][j] == 0)
					same_diff_ClusterCount += 1;
				else if(clusterArray[i][j] == 0 && groundTruthArray[i][j] == 1)
					diff_same_ClusterCount += 1;
				else
					diffClusterCount += 1;
			}
		}
		
		//Jaccard Coefficient
		externalIndex = (double)sameClusterCount/(double)(sameClusterCount + same_diff_ClusterCount + diff_same_ClusterCount);
		
		//Rand Coefficient: 
		//externalIndex = (double)(sameClusterCount + diffClusterCount)/(double)(sameClusterCount + same_diff_ClusterCount + diff_same_ClusterCount + diffClusterCount);
		return externalIndex;
	}
	
	public void silhouetteCoefficient(Map<Integer, Integer> clusterMap, double[][] distanceMatrix, int numClusters)
	{
		List<List<Integer>> clusterList = new ArrayList<>();
		for(int id = 1; id<numClusters; id++)
		{
			List<Integer> cluster = new ArrayList<>();
			for(Map.Entry<Integer,Integer> entry:clusterMap.entrySet())
			{
				int geneIndex = entry.getKey();
				int clusterId = entry.getValue();

				if(clusterId == id)
				{
					cluster.add(geneIndex);
				}
			}
			clusterList.add(cluster);
		}
		Silhouette silhouetteCoefficient = new Silhouette(distanceMatrix);
		silhouetteCoefficient.calculate(clusterList);
	}
}

