package test;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
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
	// int[][] edges;
	float flow_max = 0;

	// float bv=0.5f;
	Vector<Integer> path = new Vector<Integer>();
	Vector<Integer> path_ter = new Vector<Integer>();

	/**
	 * 执行最大流算法，进行预测转录本的定量分析
	 * 
	 * @param src
	 *            源点
	 * @param des
	 *            汇点
	 * @param edges
	 *            边的权重
	 * @param node
	 *            顶点信息
	 * @param kh
	 * @param bv
	 *            偏差因数
	 * @param sg
	 * @param path_highest
	 *            目前运行最大流的路径
	 * @return
	 * @throws FileNotFoundException
	 */
	public float max_flow(int src, int des, int[][] edges, int[] node, kmerHash kh, float[] bv, SplicingGraph sg,
			Vector<Integer> path_highest) throws FileNotFoundException {
		// TODO Auto-generated method stub
		float[][] new_edge = new float[edges.length * 2][edges.length * 2];
		rewrite_graph(edges, node, new_edge, bv);
		// 设置每条边的流量为0

		float[][] fl = new float[new_edge.length][new_edge.length];
		for (int i = 0; i < new_edge.length; i++) {
			for (int j = 0; j < new_edge.length; j++) {
				if (new_edge[i][j] != 0) {
					fl[i][j] = 0;
				}
			}
		}
		des = des + edges.length;
		int temp = des;
		int prefix[] = new int[new_edge.length];// 记录每个节点的前驱
		while (augmenting_path(src, des, new_edge, prefix) != -1) {
			// 只要找到增广路径，就一直循环
			// for(int i=0;i<path.size();i++){
			// if(i%2==0)
			// }
			for (int i = 0; i < edges.length * 2; i++) {
				if (prefix[temp] != -1) {
					// System.out.print(temp+"->");
					// temp=prefix[temp];
					path.add(0, temp);
					temp = prefix[temp];
				}
			}
			temp = des;
			path.add(0, src);

			for (int i = 0; i < path.size(); i++) {
				if (i % 2 == 0) {
					path_ter.add(path.get(i));
				}
			}
			// System.out.print(path_ter);
			// System.out.println();
			int start = 0;
			int end = 0;
			for (int i = 0; i < path_highest.size(); i++) {
				// 寻找起点
				if (path_highest.get(i) == path_ter.get(0)) {
					start = i;
				}
				// 寻找终点
				if (path_highest.get(i) == path_ter.lastElement()) {
					end = i;
				}
			}
			String path_str = "";
			for (int i = start; i <= end; i++) {
				path_str += sg.node_set.get(path_highest.get(i)).getSequence();
			}
			int u = des;
			float bias[] = new float[new_edge.length];
			bias[u] = 1;
			float increament = Float.MAX_VALUE;
			while (prefix[u] != -1) {
				int v = prefix[u];
				increament = Math.min(increament, new_edge[v][u]);
				/*
				 * if(v+edges.length==u){ bias[v]=bias[u]*bv[v]; }else{
				 * bias[v]=1; }
				 */
				/*
				 * if(prefix[v]!=-1){ int w=prefix[v]; if(v+edges.length!=u){
				 * //v和u不是同一个节点 v comes before u if(w+edges.length!=v){ //w
				 * comes brfore v bias[v]=bias[u]*bv[v]; } else{
				 * bias[v]=bias[u]; } } else{ if(w+edges.length!=v){ //w comes
				 * before v bias[v]=bias[u]; } else{ bias[v]=bias[u]/bv[v]; }
				 * 
				 * } }
				 */
				u = v;
			}

			u = des;
			while (prefix[u] != -1) {
				int v = prefix[u];
				// fl[v][u]+=increament/bias[u];
				new_edge[v][u] -= increament / bias[u];
				// fl[u][v]-=increament/bias[u];
				new_edge[u][v] += increament / bias[u];
				u = v;
			}
			flow_max += increament;
			// sg.reset_kmerMap(kh, increament, path_str);
			path.clear();
			path_ter.clear();
		}
		// System.out.println("最大流为："+flow_max);
		return flow_max;
	}

	/**
	 * 寻找流网络中的增广路径
	 * 
	 * @param src
	 *            源点
	 * @param des
	 *            汇点
	 * @param new_edge
	 *            边的权重
	 * @param prefix
	 *            前驱
	 * @return 如果找到了增广路径返回true，否则返回false
	 */
	public float augmenting_path(int src, int des, float[][] new_edge, int[] prefix) {
		Stack<Integer> stack = new Stack<Integer>();
		// path.clear();
		float flow[] = new float[new_edge.length];// 记录到当前节点的最大流
		for (int i = 0; i < new_edge.length; i++) {
			prefix[i] = -1;// 开始时设置前驱为-1
		}
		prefix[src] = -1;
		flow[src] = Float.MAX_VALUE;
		// queue.add(src);
		stack.add(src);
		while (!stack.isEmpty()) {
			int index = stack.pop();
			// path.add(index);
			if (index == des) {
				break;
			}
			int max_cov = 0;
			// int adjacent[];
			Map<Integer, Float> adjacent = new HashMap<Integer, Float>();
			for (int i = 0; i < new_edge.length; i++) {
				// 选出与当前顶点相连的最大cov的边
				if (i != src && new_edge[index][i] > 0 && prefix[i] == -1) {
					adjacent.put(i, new_edge[index][i]);
				}

			}

			// 对于可选边的权值按照从小到大排序
			List<Map.Entry<Integer, Float>> adjacent_list1 = sort_adjacent(adjacent);
			for (Map.Entry<Integer, Float> mapping : adjacent_list1) {
				int node = mapping.getKey();
				prefix[node] = index;
				stack.add(node);
				flow[node] = Math.min(mapping.getValue(), flow[index]);
				// break;
			}

		}
		if (prefix[des] == -1) {
			// 找不到增广路径
			return -1;
		} else {

			return flow[des];
		}
	}

	/**
	 * 对从某一个顶点开始的所有边按照权重降序排序，运行最大流的过程中优先处理较重的边
	 * 
	 * @param adjacent
	 *            邻接顶点和边上的权重
	 * @return 将排序的结果保存在list中，返回
	 */
	public List sort_adjacent(Map<Integer, Float> adjacent) {
		List<Map.Entry<Integer, Float>> adjacent_list = new ArrayList<Map.Entry<Integer, Float>>(adjacent.entrySet());
		// 通过比较器实现比较排序
		Collections.sort(adjacent_list, new Comparator<Map.Entry<Integer, Float>>() {
			public int compare(Map.Entry<Integer, Float> o1, Map.Entry<Integer, Float> o2) {
				return o1.getValue().compareTo(o2.getValue()); // 倒序
			}
		});
		return adjacent_list;

	}

	/**
	 * 重新构图，在运行最大流的过程中，我们把顶点一分为二，分为Vin和Vout，并将这两个顶点之间的边的权重设置为原来的V点的顶点权重。
	 * 
	 * @param edges
	 *            原来的边信息
	 * @param node
	 *            顶点信息
	 * @param new_edge
	 *            重新构图之后边的信息
	 * @param bv
	 *            偏差因数
	 */
	public void rewrite_graph(int[][] edges, int[] node, float[][] new_edge, float[] bv) {
		int size = edges.length;
		// System.out.println(size);
		size = 2 * size;
		for (int i = 0; i < edges.length; i++) {
			// 处理每个顶点
			new_edge[i][i + edges.length] = node[i];
			for (int j = 0; j < edges.length; j++) {
				if (edges[i][j] != 0) {
					// 处理i的子节点
					new_edge[i + edges.length][j] = edges[i][j];
				}
			}
		}
		// for(int i=0;i<size;i++){
		// for(int j=0;j<size;j++){
		// System.out.print(new_edge[i][j]+" ");
		// }
		// System.out.println();
		// }

	}
	/*
	 * public static void main(String args[]){ int[][] edges=new int[6][6];
	 * edges[0][4]=3; edges[0][5]=5; edges[5][4]=3; edges[0][4]=3;
	 * edges[4][2]=7; edges[4][3]=2; edges[2][1]=3; edges[3][1]=3; int[]
	 * node=new int[6]; node[0]=10; node[4]=6; node[5]=4; node[2]=6; node[3]=2;
	 * node[1]=6; for(int i=0;i<6;i++){ for(int j=0;j<6;j++){
	 * System.out.print(edges[i][j]+"   "); } System.out.println(); }
	 * System.out.println("重构图！"); // test ts=new test(); int [][] new_edge=new
	 * int[edges.length*2][edges.length*2]; flow_network fn=new flow_network();
	 * //fn.rewrite_graph(edges,node, edges); fn.max_flow(0, 7, edges, node);
	 * 
	 * }
	 */
}
