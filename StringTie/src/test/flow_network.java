package test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Vector;

public class flow_network {
	int prefix[];//记录每个节点的前驱
	int flow[];//记录到当前节点的最大流
	int node_size;
	int[][] edges;
	int flow_max=0;
	int [][]f;
	int bias[];
	int increament=Integer.MAX_VALUE;
public void max_flow(int src,int des) {
	// TODO Auto-generated method stub
	for(int i=0;i<node_size;i++){
		for(int j=0;j<node_size;j++){
			System.out.print(edges[i][j]+"    ");
		}
		System.out.println();
	}
	//设置每条边的流量为0
	for(int i=0;i<node_size;i++){
		for(int j=0;j<node_size;j++){
			f[i][j]=0;
		}
	}
	while(augmenting_path(src,des)!=-1){
		//只要找到增广路径，就一直循环
		int u=des;
		bias[u]=1;
		if(prefix[u]!=-1){
			int v=prefix[u];
			increament=Math.min(increament, edges[v][u]-bias[u]*f[v][u]);
		}
		
	}
	
}
public int augmenting_path(int src,int des){
	Queue<Integer> queue=new LinkedList<Integer>();
	for(int i=0;i<node_size;i++){
		prefix[i]=-1;//开始时设置前驱为-1
	}
	prefix[src]=-2;
	queue.add(src);
	while(!queue.isEmpty()){
		int index=queue.poll();
		if(index==des){
			break;
		}
		int max_cov=0;
		//int adjacent[];
		Map<Integer,Integer> adjacent=new HashMap<Integer,Integer>();
		for(int i=0;i<node_size;i++){
			//选出与当前顶点相连的最大cov的边
			if(i!=src&&edges[index][i]>0&&prefix[i]==-1){
				adjacent.put(i, edges[index][i]);
			}
			
		}
		//对于可选边的权值按照从大到小排序
		List<Map.Entry<Integer, Integer>> adjacent_list1=sort_adjacent(adjacent);
		 for (Map.Entry<Integer, Integer> mapping : adjacent_list1) {
               int node=mapping.getKey();
               prefix[node]=index;
               queue.add(node);
               flow[node]=Math.min(mapping.getValue(), flow[index]);
				 }
	
	}
	if(prefix[des]==-1){
		//找不到增广路径
		return -1;
	}
	else{
		return flow[des];
	}
}
public List sort_adjacent(Map<Integer,Integer> adjacent){
	List<Map.Entry<Integer, Integer>> adjacent_list = new ArrayList<Map.Entry<Integer,Integer>>(adjacent.entrySet());
	// 通过比较器实现比较排序
	Collections.sort(adjacent_list, new Comparator<Map.Entry<Integer, Integer>>() {
		public int compare(Map.Entry<Integer, Integer> o1, Map.Entry<Integer, Integer> o2) {
			return o2.getValue().compareTo(o1.getValue()); // 倒序
		}
	});
	return adjacent_list;

}


}
