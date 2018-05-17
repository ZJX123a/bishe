package test;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;
import test.baseOptions;
import test.load_Read;
import test.flow_network;

public class kmerHash {
	int kmer_length;
	final int min_exon_length = 20;
	int fr_strand;
	final int min_seed_cov = 2;
	final float min_seed_entry = 1.3f;
	final float g_min_ratio_non_error = 0.5f;
	final int pair_gap_length = 200;
	final int max_pair_gap_length = 500;
	final int min_anchor_length = 23;
	final int min_function_count = 2;
	final boolean is_paired_end = true;
	final float min_ratio_welds = 0.04f;
	final int min_reads_span_junction = 2;
	Map<Long, Vector<Integer>> kmer_hash = new HashMap<Long, Vector<Integer>>();

	Map<Long, Integer> kmer_map = new HashMap<Long, Integer>();
	List<Map.Entry<Long, Integer>> list;

	/**
	 * 利用解析好的reads数据构建k-mer字典
	 * 对于每一个read，按照k-mer大小截取所有k-mer，并统计目前k-mer出现的reads的标号，将k-
	 * mer和所覆盖的标号一同存入kmer_hash
	 * 
	 * @param kmer_hash
	 *            将构建好的k-mer字典的结果保存在kmer_hash中
	 * @param read_vector
	 *            输入
	 */
	public void readsToKmer(Map<Long, Vector<Integer>> kmer_hash, Vector<String> read_vector) {
		// TODO Auto-generated method stub
		String temp_kmer;
		long intval = 0;
		for (int i = 0; i < read_vector.size(); i++) {
			String read = read_vector.get(i);
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

	/**
	 * 得到某一k-mer覆盖的所有reads的标号的集合
	 * 
	 * @param intval
	 *            k-mer的long型值
	 * @param kmer_hash
	 * @return 返回Vector，Vector中含有k-mer覆盖的所有reads标号
	 */
	public Vector<Integer> get_readset(long intval, Map<Long, Vector<Integer>> kmer_hash) {
		return (Vector<Integer>) kmer_hash.get(intval);
	}

	/**
	 * 判断某一k-mer在kmer_hash中是否存在
	 * 
	 * @param kmer_hash
	 * @param intval
	 * @return 如果存在返回true，否则返回false
	 */
	private boolean ifexist(Map<Long, Vector<Integer>> kmer_hash, long intval) {
		if (kmer_hash.containsKey(intval)) {
			return true;
		} else {
			return false;
		}
	}

	/**
	 * 得到k-mer覆盖的reads数目
	 * 
	 * @param kmer_hash
	 * @param intval
	 * @return
	 */
	public int get_readset_count(Map<Long, Vector<Integer>> kmer_hash, long intval) {
		if (ifexist(kmer_hash, intval)) {
			Vector readset = (Vector) kmer_hash.get(intval);
			int count = readset.size();
			return count;
		} else {
			return 0;
		}
	}

	/**
	 * 刪除一些不符合要求的k-mer
	 * 
	 * @param kmer_hash
	 */
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

	/**
	 * 得到正向的candidates，用于后续的正向延申 利用种子k-mer（seed_kmer）通过位移得到forward_candidates，
	 * 并根据每个forward_candidate的reads覆盖度进行排序，将排序的结果存放在list中
	 *
	 * @param seed_kmer
	 *            种子k-mer
	 * @param kmer_hash
	 * @return 返回排好序的list
	 */
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

			}
			list_forward = new ArrayList<Map.Entry<Long, Integer>>(forward_candidates.entrySet());

			// 通过比较器实现比较排序
			Collections.sort(list_forward, new Comparator<Map.Entry<Long, Integer>>() {
				public int compare(Map.Entry<Long, Integer> o1, Map.Entry<Long, Integer> o2) {
					return o2.getValue().compareTo(o1.getValue()); // 倒序
				}
			});
			return list_forward;
		}
		return null;

	}

	/**
	 * 得到反向的candidates，用于后续的反向延伸， 利用种子k-mer（seed_kmer）通过位移得到reverse_candidates，
	 * 并根据每个reverse_candidate的reads覆盖度进行排序，将排序的结果存放在list中
	 *
	 * @param seed_kmer
	 *            种子k-mer
	 * @param kmer_hash
	 * @return 返回排好序的list
	 */
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

	/**
	 * 打印k-mer字典
	 * 
	 * @param kmer_hash
	 */
	private void print_kmerhash(Map<Long, Vector<Integer>> kmer_hash) {
		System.out.println(kmer_hash.size());
		for (Iterator<Long> iterator = kmer_hash.keySet().iterator(); iterator.hasNext();) {
			Long key = iterator.next();
			System.out.println("key-----" + baseOptions.intvalToKmer(key, kmer_length) + "     " + key);
			System.out.println("value--------" + kmer_hash.get(key));
		}
	}

	/**
	 * 对于k-mer字典中的所有k-mer，按照其reads覆盖度进行排序
	 * 
	 * @param kmer_hash
	 * @return 将排序的结果放入list中，返回list
	 */
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

}
