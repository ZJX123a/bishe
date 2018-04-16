package test;

import java.io.IOException;
import java.util.*;
import test.baseOptions;
import test.load_Read;
import test.flow_network;

public class kmerHash {
	final int kmer_length = 25;
	final int min_exon_length = 20;
	final int fr_strand = 2;
	final int min_seed_cov = 2;
	final float min_seed_entry = 1.3f;
	final float g_min_ratio_non_error = 0.5f;
	final int pair_gap_length = 200;
	final int max_pair_gap_length = 500;
	final int min_anchor_length = 20;
	final int min_function_count = 2;
	final boolean is_paired_end = true;
	final float min_ratio_welds = 0.04f;
	final int min_reads_span_junction = 2;
	Map<Long, Vector<Integer>> kmer_hash = new HashMap<Long, Vector<Integer>>();

	Map<Long, Integer> kmer_map = new HashMap<Long, Integer>();
	List<Map.Entry<Long, Integer>> list;

	public void readsToKmer(Map<Long, Vector<Integer>> kmer_hash) {
		// TODO Auto-generated method stub
		String temp_kmer;
		long intval = 0;
		for (int i = 0; i < load_Read.read_vector.size(); i++) {
			String read = load_Read.read_vector.get(i);
			for (int j = 0; j < read.length() - kmer_length + 1; j++) {
				temp_kmer = read.substring(j, j + kmer_length);
				intval = baseOptions.kmerToIntval(temp_kmer);
				Vector<Integer> readset = new Vector<Integer>();
				if (kmer_hash.containsKey(intval)) {
					readset = get_readset(intval, kmer_hash);
					readset.add(i);
				} else {
					readset.add(i);
				}
				kmer_hash.put(intval, readset);
			}

		}
	}

	public Vector<Integer> get_readset(long intval, Map<Long, Vector<Integer>> kmer_hash) {
		return (Vector<Integer>) kmer_hash.get(intval);
	}

	private boolean ifexist(Map<Long, Vector<Integer>> kmer_hash, long intval) {
		if (kmer_hash.containsKey(intval)) {
			return true;
		} else {
			return false;
		}
	}

	public int get_readset_count(Map<Long, Vector<Integer>> kmer_hash, long intval) {
		if (ifexist(kmer_hash, intval)) {
			Vector readset = (Vector) kmer_hash.get(intval);
			int count = readset.size();
			return count;
		} else {
			return 0;
		}
	}

	private void delete_bad_kmers(Map<Long, Vector<Integer>> kmer_hash) {
		Vector delete_list = new Vector();
		for (Iterator<Long> iterator = kmer_hash.keySet().iterator(); iterator.hasNext();) {
			List<Map.Entry<Long, Long>> candidates = new ArrayList<Map.Entry<Long, Long>>();
			Long key = iterator.next();
			Long candidate;
			candidates = get_forward_candidates(key, kmer_hash);
			long dominate_count = 0;
			if (candidates.size() > 0) {
				for (Map.Entry<Long, Long> mapping : candidates) {
					candidate = mapping.getKey();
					if (mapping.getValue() > 0) {
						long candidate_count = mapping.getValue();
						if (dominate_count == 0) {
							dominate_count = candidate_count;
						} else if ((float) candidate_count / dominate_count < g_min_ratio_non_error) {
							delete_list.add(candidate);
						}
					}
				}
			}
		}
		System.out.println("之前的hash");
		System.out.println(kmer_hash.size());
		if (!delete_list.isEmpty()) {
			for (int j = 0; j < delete_list.size(); j++) {
				long delete_kmer = (long) delete_list.get(j);
				kmer_hash.remove(delete_kmer);
			}
		}
		System.out.println(kmer_hash.size());
	}

	private boolean if_find_kmer(long intval, Map<Long, Vector<Integer>> kmer_hash) {
		if (kmer_hash.containsKey(intval)) {
			return true;
		} else {
			return false;
		}
	}

	public List get_forward_candidates(long seed_kmer, Map<Long, Vector<Integer>> kmer_hash) {
		if (ifexist(kmer_hash, seed_kmer)) {
			Map forward_candidates = new HashMap();
			List<Map.Entry<Long, Integer>> list_forward = new ArrayList<Map.Entry<Long, Integer>>();
			long temp_intval = seed_kmer << (65 - kmer_length * 2) >> (63 - kmer_length * 2);
			long a[] = { 0l, 1l, 2l, 3l };
			for (int i = 0; i < 4; i++) {
				long temp2 = temp_intval;
				temp2 |= a[i];
				temp2 = baseOptions.kmerToIntval(baseOptions.intvalToKmer(temp2, kmer_length));
				if (ifexist(kmer_hash, temp2)) {
					int read_count = get_readset_count(kmer_hash, temp2);
					forward_candidates.put(temp2, read_count);
				}
				list_forward = new ArrayList<Map.Entry<Long, Integer>>(forward_candidates.entrySet());

				// 通过比较器实现比较排序
				Collections.sort(list_forward, new Comparator<Map.Entry<Long, Integer>>() {
					public int compare(Map.Entry<Long, Integer> o1, Map.Entry<Long, Integer> o2) {
						return o2.getValue().compareTo(o1.getValue()); // 倒序
					}
				});
			}
			// System.out.println(forward_candidates);
			// System.out.println(list_forward);
			// System.out.println();
			return list_forward;
		}
		return null;

	}

	public List get_reverse_candidates(long seed_kmer, Map<Long, Vector<Integer>> kmer_hash) {
		if (ifexist(kmer_hash, seed_kmer)) {
			Map reverse_candidates = new HashMap();
			List<Map.Entry<Long, Integer>> list_reverse = new ArrayList<Map.Entry<Long, Integer>>();
			long temp_intval = seed_kmer >> 2;
			long a[] = { 0l, 1l, 2l, 3l };
			for (int i = 0; i < 4; i++) {
				long temp2 = (a[i] << (2 * kmer_length - 2)) | temp_intval;
				temp2 = baseOptions.kmerToIntval(baseOptions.intvalToKmer(temp2, kmer_length));
				if (ifexist(kmer_hash, temp2)) {
					int read_count = get_readset_count(kmer_hash, temp2);
					// System.out.println(temp2+" "+read_count);
					reverse_candidates.put(temp2, read_count);
				}

			}
			// if(reverse_candidates.size()==0){
			// return null;
			// }
			list_reverse = new ArrayList<Map.Entry<Long, Integer>>(reverse_candidates.entrySet());

			// 通过比较器实现比较排序
			Collections.sort(list_reverse, new Comparator<Map.Entry<Long, Integer>>() {
				public int compare(Map.Entry<Long, Integer> o1, Map.Entry<Long, Integer> o2) {
					return o2.getValue().compareTo(o1.getValue()); // 倒序
				}
			});
			return list_reverse;
		} else {
			return null;
		}

	}

	private void print_kmerhash(Map<Long, Vector<Integer>> kmer_hash) {
		System.out.println(kmer_hash.size());
		for (Iterator<Long> iterator = kmer_hash.keySet().iterator(); iterator.hasNext();) {
			Long key = iterator.next();
			System.out.println("key-----" + baseOptions.intvalToKmer(key, kmer_length) + "     " + key);
			System.out.println("value--------" + kmer_hash.get(key));
		}
	}

	public List sort_kmer(Map<Long, Vector<Integer>> kmer_hash) {

		Long key;
		int count;
		for (Iterator<Long> iterator = kmer_hash.keySet().iterator(); iterator.hasNext();) {
			key = iterator.next();
			count = kmer_hash.get(key).size();
			kmer_map.put(key, count);
		}
		list = new ArrayList<Map.Entry<Long, Integer>>(kmer_map.entrySet());

		// 通过比较器实现比较排序
		Collections.sort(list, new Comparator<Map.Entry<Long, Integer>>() {
			public int compare(Map.Entry<Long, Integer> o1, Map.Entry<Long, Integer> o2) {
				return o2.getValue().compareTo(o1.getValue()); // 倒序
			}
		});

		// for (Map.Entry<Long, Long> mapping : list) {
				// System.out.println(baseOptions.intvalToKmer(mapping.getKey(),
				// kmer_length) + ":" + mapping.getValue()
				// + " intval:" + mapping.getKey());
				// }
		return list;
	}

	public static void main(String args[]) throws IOException {
		load_Read.load_reads();
		kmerHash kh = new kmerHash();

		kh.readsToKmer(kh.kmer_hash);
		kh.sort_kmer(kh.kmer_hash);
		List listK = kh.sort_kmer(kh.kmer_hash);
		if (listK.size() == 0) {
			System.out.println("没有数据！");
			return;
		} else {
			System.out.println("有" + listK.size() + "个可用kmer");
		}
		int count = 0;
		Set node_jihe = new HashSet();
		Map<Long, Integer> used_kmers_plus = new HashMap<Long, Integer>();
		// System.out.println("初始used_kmers:"+sg.used_kmers.size());
		for (int i = 0; i < listK.size(); i++) {
			if (!used_kmers_plus.containsKey(kh.list.get(i).getKey())) {
				SplicingGraph sg = new SplicingGraph();
				if (sg.init_trunk(kh, kh.list.get(i).getKey(), node_jihe, sg)) {
					//System.out.println(sg.node_set.get(0).getSequence());
					sg.forward_check_and_extend(kh, 0);
					sg.reverse_check_and_extend(kh);
					sg.init_parents();
					// sg.check_again();

					while (sg.if_can_extend(kh)) {
						for (int k = 0; k < sg.node_set.size(); k++) {
							if (sg.forward_branches.contains(sg.node_set.get(k))
									|| sg.reverse_branches.contains(sg.node_set.get(k))) {
								sg.forward_check_and_extend(kh, k);
								sg.reverse_check_and_extend(kh);
							}
						}
					}
					sg.compute_node_cov(kh);
					int[][] edges = new int[sg.node_set.size()][sg.node_set.size()];
					sg.compute_edge_cov(kh,edges);
					float[] bv=new float[sg.node_set.size()];
					//sg.compute_bv(kh, edges, bv);
					flow_network fn=new flow_network();
				//	fn.max_flow(edges, sg.node_set.size());
					
					for (int j = 0; j < sg.node_set.size(); j++) {
					//	System.out.println("顶点编号：" + j + "    cov:"+sg.node_set.get(j).getcov()+"     顶点序列:" + sg.node_set.get(j).getSequence());
					//	System.out.println("父节点：" + sg.node_set.get(j).getParents());
						//System.out.println("子节点：" + sg.node_set.get(j).getChildren());
					}
//					for(int m=0;m<sg.node_set.size();m++){
//						for(int n=0;n<sg.node_set.size();n++){
//							System.out.print(edges[m][n]+"     ");
//						}
//						System.out.println();
//					}
					//System.out.println("运行结束！当前节点个数为：" + sg.node_set.size());
					used_kmers_plus.putAll(sg.used_kmers);
				}
				count++;
			}
		}
		SplicingGraph sg3=new SplicingGraph();
		sg3.compute_node("TCATCTATAAAATACTGAAAATATCATTTTAAG", kh);
		System.out.println("运行结束！");
		System.out.println("开始构图！");
	
		
		//flow_network fn = new flow_network();
//		 for (Iterator i = kh.kmer_map.keySet().iterator(); i.hasNext();) {  
//			   Object obj = i.next();  
//			   //System.out.print(obj+"   ");// 循环输出key  
//			   System.out.println("key=" + obj + " value=" + kh.kmer_map.get(obj));  
//			  }  
//		System.out.println(baseOptions.intvalToKmer(255593747292392l, kh.kmer_length));

	}
}
