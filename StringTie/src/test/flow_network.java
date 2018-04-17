package test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Stack;
import java.util.Vector;

public class flow_network {
	int node_size;
	//int[][] edges;
	float flow_max=0;
	
	
	float bv=0.5f;
	
public void max_flow(int src,int des,int[][] edges,int[] node) {
	// TODO Auto-generated method stub
	for(int i=0;i<node_size;i++){
		for(int j=0;j<node_size;j++){
			System.out.print(edges[i][j]+"    ");
		}
		System.out.println();
	}
	int[][] new_edge=new int[edges.length*2][edges.length*2];
	rewrite_graph(edges, node, new_edge);
	for(int i=0;i<edges.length*2;i++){
		for(int j=0;j<edges.length*2;j++){
			System.out.print(new_edge[i][j]+"   ");
		}
		System.out.println();
	}
	//设置每条边的流量为0
	float [][]f=new float[new_edge.length][new_edge.length];
	for(int i=0;i<new_edge.length;i++){
		for(int j=0;j<new_edge.length;j++){
			if(new_edge[i][j]!=0){
				f[i][j]=0;
			}
		}
	}
	int prefix[]=new int[new_edge.length];//记录每个节点的前驱
	while(augmenting_path(src,des, new_edge,prefix)!=-1){
		//只要找到增广路径，就一直循环
		int u=des;
		float bias[]=new float[new_edge.length];
		bias[u]=1;
		float increament=Float.MAX_VALUE;
		while(prefix[u]!=-1){
			int v=prefix[u];
			increament=Math.min(increament, new_edge[v][u]-bias[u]*f[v][u]);
			if(prefix[v]!=-1){
				int w=prefix[v];
				if(v+edges.length!=u){
					//v和u不是同一个节点   v comes before u
					if(w+edges.length!=v){
						//w comes brfore v
						bias[v]=bias[u]*bv;
					}
					else{
						bias[v]=bias[u];
					}
				}
				else{
					if(w+edges.length!=v){
						//w comes before v
						bias[v]=bias[u];
					}
					else{
						bias[v]=bias[u]/bv;
					}
					
				}
			}
			u=v;
		}
		
		u=des;
		while(prefix[u]!=-1){
			int v=prefix[u];
			//f[u][v]+=increament/bias[u];
			new_edge[v][u]-=increament/bias[u];
			//f[v][u]-=increament/bias[u];
			new_edge[u][v]+=increament/bias[u];
			u=v;
		}
		flow_max+=increament;
	}
	System.out.println(flow_max);
	
}
public float augmenting_path(int src,int des,int[][] new_edge,int[] prefix){
	Stack<Integer> stack=new Stack<Integer>();
	float flow[]=new float[new_edge.length];//记录到当前节点的最大流
	for(int i=0;i<new_edge.length;i++){
		prefix[i]=-1;//开始时设置前驱为-1
	}
	prefix[src]=-1;
	flow[src]=Float.MAX_VALUE;
//	queue.add(src);
	stack.add(src);
	while(!stack.isEmpty()){
		int index=stack.pop();
		if(index==des){
			break;
		}
		int max_cov=0;
		//int adjacent[];
		Map<Integer,Integer> adjacent=new HashMap<Integer,Integer>();
		for(int i=0;i<new_edge.length;i++){
			//选出与当前顶点相连的最大cov的边
			if(i!=src&&new_edge[index][i]>0&&prefix[i]==-1){
				adjacent.put(i, new_edge[index][i]);
			}
			
		}
		//对于可选边的权值按照从小到大排序
		List<Map.Entry<Integer, Integer>> adjacent_list1=sort_adjacent(adjacent);
		 for (Map.Entry<Integer, Integer> mapping : adjacent_list1) {
               int node=mapping.getKey();
               prefix[node]=index;
               stack.add(node);
               flow[node]=Math.min(mapping.getValue(), flow[index]);
             //  break;
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
			return o1.getValue().compareTo(o2.getValue()); // 倒序
		}
	});
	return adjacent_list;

}
public void rewrite_graph(int[][] edges,int[] node,int[][] new_edge){
	int size=edges.length;
	System.out.println(size);
	size=2*size;
	
	for(int i=0;i<edges.length;i++){
		//处理每个顶点
		new_edge[i][i+edges.length]=node[i];
		for(int j=0;j<edges.length;j++){
			if(edges[i][j]!=0){
				//处理i的子节点
				new_edge[i+edges.length][j]=edges[i][j];
			}
		}
	}
	
}

	public static void main(String args[]){
		int[][] edges=new int[6][6];
		edges[0][4]=3;
		edges[0][5]=5;
		edges[5][4]=3;
		edges[0][4]=3;
		edges[4][2]=7;
		edges[4][3]=2;
		edges[2][1]=3;
		edges[3][1]=3;
		int[] node=new int[6];
		node[0]=10;
		node[4]=6;
		node[5]=4;
		node[2]=6;
		node[3]=2;
		node[1]=6;
		for(int i=0;i<6;i++){
			for(int j=0;j<6;j++){
				System.out.print(edges[i][j]+"   ");
			}
			System.out.println();
		}
		System.out.println("重构图！");
		test ts=new test();
		//ts.max_flow(0, 7, edges, node);
		//ts.rewrite_graph(edges,node);
 }
}
