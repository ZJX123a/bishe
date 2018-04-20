package test;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;
import java.util.Vector;

public class enter {
	public static void print_message(kmerHash kh){
		System.out.println("基于流分解的转录组组装程序开始运行..........");
		System.out.println("请输入相关参数..........");
		System.out.println("请输入kmer长度（kmer长度需在15-31之间，默认25）.........");
		Scanner scan = new Scanner(System.in);
		int read = scan.nextInt();
		while(read>31||read<15){
			System.out.println("您输入的数据不合法，请重新输入！");
			read=scan.nextInt();
		}
		System.out.println("您输入的kmer长度是："+read);
		kh.kmer_length=read;
		System.out.println("请输入reads测序方式，1：(2)-><-(1)   2:(1)-><-(1)");
		read=scan.nextInt();
		while(read!=1&&read!=2){
			System.out.println("您输入的数据不合法，请重新输入！");
			read=scan.nextInt();
		}
		kh.fr_strand=read;
	}
	
	
	
	public static void main(String args[]) throws IOException {
		kmerHash kh = new kmerHash();
		print_message(kh);
		System.out.println("开始构造kmerHash！");
		long startTime = System.currentTimeMillis();    //获取开始时间
		load_Read.load_reads();
		kh.readsToKmer(kh.kmer_hash);
		kh.sort_kmer(kh.kmer_hash);
		long endTime = System.currentTimeMillis();    //获取开始时间
		System.out.println("kmerHash构造完成！运行时间为："+(endTime-startTime)+"ms.");
		List listK = kh.sort_kmer(kh.kmer_hash);
		if (listK.size() == 0) {
			System.out.println("没有数据！");
			return;
		} 
		int count = 0;
		Set node_jihe = new HashSet();
		Map<Long, Integer> used_kmers_plus = new HashMap<Long, Integer>();
		 FileOutputStream out=new FileOutputStream("D:/test.txt");
         PrintStream p=new PrintStream(out);
         System.out.println("开始构图并计算最大流！");
         startTime = System.currentTimeMillis();
		for (int i = 0; i < listK.size(); i++) {
			if (!used_kmers_plus.containsKey(kh.list.get(i).getKey())) {
				SplicingGraph sg = new SplicingGraph();
				if (sg.init_trunk(kh, kh.list.get(i).getKey(), node_jihe, sg)) {
					sg.forward_check_and_extend(kh, 0);
					sg.reverse_check_and_extend(kh);
					sg.init_parents();
					while (sg.if_can_extend(kh)) {
						for (int k = 0; k < sg.node_set.size(); k++) {
							if (sg.forward_branches.contains(sg.node_set.get(k))
									|| sg.reverse_branches.contains(sg.node_set.get(k))) {
								sg.forward_check_and_extend(kh, k);
								sg.reverse_check_and_extend(kh);
							}
						}
					}

					for (int j = 0; j < sg.node_set.size(); j++) {
						System.out.println("顶点编号：" + j + "    cov:" + sg.node_set.get(j).getcov() + "     顶点序列:"
								+ sg.node_set.get(j).getSequence());
						System.out.println("父节点：" + sg.node_set.get(j).getParents());
						System.out.println("子节点：" + sg.node_set.get(j).getChildren());
					}
					
//			         for(int j=0;j<sg.node_set.size();j++){
//			        	 String node=sg.node_set.get(j).getSequence();
//			        	 for(int k=0;k<=node.length()-kh.kmer_length;k++){
//			        		 String kmer=node.substring(k,k+kh.kmer_length);
//			        		 p.append(kmer+":"+kh.kmer_map.get(baseOptions.kmerToIntval(kmer))+"        ");
//			        	 }
//			         }
					Vector<Integer> src = new Vector<Integer>();
					Vector<Integer> des = new Vector<Integer>();
					int[][] f = new int[sg.node_set.size()][sg.node_set.size()];
					for (int k = 0; k < sg.node_set.size(); k++) {
						if (sg.node_set.get(k).getChildren().size() == 0) {
							des.add(k);
						}
						if (sg.node_set.get(k).getParents().size() == 0) {
							src.add(k);
						}
					}

					float flow_sum = 0;
					int[][] edges = new int[sg.node_set.size()][sg.node_set.size()];
					for(int k=0;k<sg.node_set.size();k++){
						for(int m=0;m<sg.node_set.get(k).getChildren().size();m++){
							int child=(int) sg.node_set.get(k).getChildren().get(m);
							edges[k][child]=1;
						}
					}
					Vector<Vector<Integer>> path_sum = new Vector<Vector<Integer>>();
					for (int k = 0; k < src.size(); k++) {
						for (int h = 0; h < des.size(); h++) {
							sg.dfs(edges, src.get(k), des.get(h), kh, f, path_sum);
							sg.path.clear();
						}
					}
					// 根据read_per_base_cov找出最重的一条路径
					List<Map.Entry<Vector<Integer>, Float>> list_path = (List) sg.find_the_highest_path(kh, path_sum);
					System.out.println("路径与其cov:");
					Map<Integer,Float> nodes_cov=new HashMap<Integer,Float>();
					for (Map.Entry<Vector<Integer>, Float> mapping : list_path) {
						Vector<Integer> path_highest=mapping.getKey();
						int[] nodes = new int[sg.node_set.size()];
						int[][] edges_ter=new int[sg.node_set.size()][sg.node_set.size()];
						String path_str = "";
						path_str=sg.build_graph(path_highest, kh,nodes_cov,nodes,edges_ter);
						for(int k=0;k<sg.node_set.size();k++){
							 System.out.print(nodes[k]+" ");
						 }
						System.out.println();
						for (int k = 0; k < sg.node_set.size(); k++) {
							for (int j = 0; j < sg.node_set.size(); j++) {
								System.out.print(edges_ter[k][j]+"   ");
							}
							System.out.println();
						}
						//设置bv
						float[] bv=new float[sg.node_set.size()];
						float fenzi=0;
						float fenmu=0;
						int self=0;
						for(int k=0;k<sg.node_set.size();k++){
							//处理子节点
							fenzi=0;
							fenmu=0;
							if(nodes[k]!=0){
							self=nodes[k];
							for(int m=0;m<sg.node_set.size();m++){
								if(edges_ter[k][m]!=0){
									fenmu+=(float)1/(edges_ter[k][m]);
									self-=edges_ter[k][m];
								}
							}
							//处理父节点
							for(int m=0;m<sg.node_set.size();m++){
								if(edges_ter[m][k]!=0){
									fenzi+=(float)1/(edges_ter[m][k]);
									self-=edges_ter[m][k];
								}
							}
							bv[k]=(float)(((float)(1.0/self)+fenzi)/((float)(1.0/self)+fenmu));
							}
						}
						
						
//						 for(int k=0;k<sg.node_set.size();k++){
//							 System.out.print(bv[k]+" ");
//						 }
//						System.out.println();
						flow_network fn=new flow_network();
						float max_flow=fn.max_flow(path_highest.get(0), path_highest.lastElement(), edges_ter, nodes, kh, bv, sg);
						p.println("转录本："+path_str);
						p.println("其转录丰度为："+max_flow);
					}
					System.out.println();
					used_kmers_plus.putAll(sg.used_kmers);
				}
				count++;
			}
		}
		endTime = System.currentTimeMillis();
		System.out.println( "运行结束！运行时间为："+(endTime-startTime)+"ms.");
		
		p.flush();
	}
}
