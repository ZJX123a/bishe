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
	int prefix[];//��¼ÿ���ڵ��ǰ��
	int flow[];//��¼����ǰ�ڵ�������
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
	//����ÿ���ߵ�����Ϊ0
	for(int i=0;i<node_size;i++){
		for(int j=0;j<node_size;j++){
			f[i][j]=0;
		}
	}
	while(augmenting_path(src,des)!=-1){
		//ֻҪ�ҵ�����·������һֱѭ��
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
		prefix[i]=-1;//��ʼʱ����ǰ��Ϊ-1
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
			//ѡ���뵱ǰ�������������cov�ı�
			if(i!=src&&edges[index][i]>0&&prefix[i]==-1){
				adjacent.put(i, edges[index][i]);
			}
			
		}
		//���ڿ�ѡ�ߵ�Ȩֵ���մӴ�С����
		List<Map.Entry<Integer, Integer>> adjacent_list1=sort_adjacent(adjacent);
		 for (Map.Entry<Integer, Integer> mapping : adjacent_list1) {
               int node=mapping.getKey();
               prefix[node]=index;
               queue.add(node);
               flow[node]=Math.min(mapping.getValue(), flow[index]);
				 }
	
	}
	if(prefix[des]==-1){
		//�Ҳ�������·��
		return -1;
	}
	else{
		return flow[des];
	}
}
public List sort_adjacent(Map<Integer,Integer> adjacent){
	List<Map.Entry<Integer, Integer>> adjacent_list = new ArrayList<Map.Entry<Integer,Integer>>(adjacent.entrySet());
	// ͨ���Ƚ���ʵ�ֱȽ�����
	Collections.sort(adjacent_list, new Comparator<Map.Entry<Integer, Integer>>() {
		public int compare(Map.Entry<Integer, Integer> o1, Map.Entry<Integer, Integer> o2) {
			return o2.getValue().compareTo(o1.getValue()); // ����
		}
	});
	return adjacent_list;

}


}
