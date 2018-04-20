package test;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;
import test.baseOptions;
import test.load_Read;
import test.flow_network;

public class kmerHash {
	int kmer_length ;
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

	
}
