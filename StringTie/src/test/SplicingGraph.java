package test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;
import java.util.Vector;

import test.SplicingGraph.Node;

public class SplicingGraph {
	Vector<Node> node_set = new Vector<Node>();
	Set<Long> restore_kmers = new HashSet<Long>();
	Set<Integer> forward_branches = new HashSet<Integer>();
	Set<Integer> reverse_branches = new HashSet<Integer>();
	Map<Long, Integer> used_kmers = new HashMap<Long, Integer>();

	public class Node {
		private Vector<Integer> parents = new Vector<Integer>();
		private Vector<Integer> children = new Vector<Integer>();
		private String sequence;
		private int cov;

		public boolean addParents(int parent) {
			if (parent < 0) {
				return false;
			} else {
				if (parents.size() != 0) {
					for (int i = 0; i < parents.size(); ++i) {
						if (parents.get(i) == parent) // if exist already
							return false;
					}
				}
				this.parents.add(parent);
				return true;
			}
		}

		public List getParents() {
			return parents;
		}

		public boolean addChildren(int child) {
			if (child < 0) {
				return false;
			} else {
				if (children.size() != 0) {
					for (int i = 0; i < children.size(); ++i) {
						if (children.get(i) == child) // if exist already
							return false;
					}
				}
				this.children.add(child);
				return true;
			}
		}

		public Vector getChildren() {
			return children;
		}

		public void setSequence(String seq) {
			sequence = seq;
		}

		public String getSequence() {
			return sequence;
		}

		public boolean isChild(int child) {
			if (child < 0) {
				return false;
			} else {
				for (int i = 0; i < children.size(); i++) {
					if (children.get(i) == child) {
						return true;
					}
				}
				return false;
			}
		}

		public boolean isParent(int parent) {
			if (parent < 0) {
				return false;
			} else {
				for (int i = 0; i < parents.size(); i++) {
					if (parents.get(i) == parent) {
						return true;
					}
				}
				return false;
			}
		}

		public boolean deleteChild(long child) {
			// TODO Auto-generated method stub
			if (child < 0) {
				return false;
			} else {
				for (int i = 0; i < children.size(); i++) {
					if (children.get(i) == child) {
						children.remove(i);
						return true;
					}
				}
				return false;
			}
		}

		public boolean deleteParent(long parent) {
			// TODO Auto-generated method stub
			if (parent < 0) {
				return false;
			} else {
				for (int i = 0; i < parents.size(); i++) {
					if (parents.get(i) == parent) {
						parents.remove(i);
						return true;
					}
				}
				return false;
			}
		}

		public void clearChildren() {
			// TODO Auto-generated method stub
			children.clear();
		}

		public void clearParents() {
			parents.clear();
		}

		public void clearAll() {
			// TODO Auto-generated method stub
			children.clear();
			parents.clear();
			sequence = "";
		}

		public void setcov(int cove) {
			cov = cove;
		}

		public int getcov() {
			// TODO Auto-generated method stub
			return cov;
		}
	}

	private String forward_extend(long kmer_int, Vector<Long> bifurcation, kmerHash kh) {
		List<Map.Entry<Long, Integer>> candidates = new ArrayList<Map.Entry<Long, Integer>>();

		long sum = 0l;
		String kmer_string = baseOptions.intvalToKmer(kmer_int, kh.kmer_length);
		while (true) {
			boolean flag = false;
			candidates = kh.get_forward_candidates(kmer_int, kh.kmer_hash);
			// System.out.println(candidates);
			if (candidates.size() == 0) {
				break;
			}
			int count = 0;
			long candidate = 0l;
			int cov = 0;
			// System.out.println("其candidates为：");
			for (Map.Entry<Long, Integer> mapping : candidates) {
				candidate = mapping.getKey();
				if (!has_been_used(candidate)) {// 如果有未被使用的candidate
					// System.out.println("未被使用");
					flag = true;
					cov = mapping.getValue();
					// System.out.println("candidate:"+candidate+" cov:"+cov);
					// 如果还有其他未被使用的candidate
					if (count < candidates.size() - 1) {
						bifurcation.add(kmer_int);
					}
					break;
				}
				count++;
				// System.out.println();
			}
			if (!flag) {
				break;
			}
			if (flag == true) {
				used_kmers.put(candidate, cov);
				if (cov == 0) {
					System.out.println(baseOptions.intvalToKmer(candidate, kh.kmer_length) + "000000000");
				}
				sum += cov;
				int base_last_int = (int) (candidate & 3l);
				char base_last_char = baseOptions.intToBase(base_last_int);
				kmer_string += base_last_char;
				kmer_int = candidate;
			}
		}
		return kmer_string;
	}

	// 反向扩展
	private String reverse_extend(long kmer_int, Vector<Long> bifurcation, kmerHash kh) {
		List<Map.Entry<Long, Integer>> candidates = new ArrayList<Map.Entry<Long, Integer>>();

		long sum = 0l;
		String kmer_string = baseOptions.intvalToKmer(kmer_int, kh.kmer_length);
		while (true) {
			candidates.clear();
			candidates = kh.get_reverse_candidates(kmer_int, kh.kmer_hash);
			// System.out.println(candidates);
			if (candidates.size() == 0) {
				break;
			}
			// int count = 0;
			long candidate = 0l;
			int cov = 0;
			// System.out.println("其candidates为：");
			boolean flag = false;
			for (Map.Entry<Long, Integer> mapping : candidates) {
				candidate = mapping.getKey();
				if (!has_been_used(candidate)) {// 如果有未被使用的candidate
					// System.out.println("未被使用");
					flag = true;
					cov = mapping.getValue();
					// if (count < candidates.size() - 1) {
					// bifurcation.add(kmer_int);
					// }
					break;
				}
				// count++;
				// System.out.println();
			}
			if (!flag) {
				break;
			}
			if (flag == true) {
				used_kmers.put(candidate, cov);
				// if(cov==0){
				// System.out.println("234:"+"000000000");
				// }
				sum += cov;
				int base_first_int = (int) (candidate >> (2 * kh.kmer_length - 2));
				char base_first_char = baseOptions.intToBase(base_first_int);
				kmer_string = base_first_char + kmer_string;
				kmer_int = candidate;
			}
		}
		return kmer_string;
	}

	public boolean has_been_used(Long candidate) {
		if (used_kmers.containsKey(candidate)) {
			return true;
		} else {
			return false;
		}
	}

	public void extend_again(kmerHash kh, Vector<Long> bifurcation) {
		List<Map.Entry<Long, Long>> candidates = new ArrayList<Map.Entry<Long, Long>>();
		if (node_set.get(0).sequence.length() > 2 * kh.kmer_length + kh.kmer_length) {
			for (int i = 0; i < 2 * kh.kmer_length; i++) {
				String kmer = node_set.get(0).sequence.substring(i, kh.kmer_length + i);
				// System.out.println(kmer + " &");
				long intval = baseOptions.kmerToIntval(kmer);
				candidates = kh.get_reverse_candidates(intval, kh.kmer_hash);
				// System.out.println("candidates.size():" + candidates.size());
				if (candidates.size() <= 1) {
					continue;
				}
				boolean flag = false;
				long candidate;
				for (Map.Entry<Long, Long> mapping : candidates) {
					candidate = mapping.getKey();
					// System.out.println("candidate:" + candidate);
					if (!has_been_used(candidate)) {
						flag = true;
						break;
					}
				}
				if (flag == true) {
					String sequence = reverse_extend(intval, bifurcation, kh);
					if (sequence.length() > kh.kmer_length + kh.min_exon_length) {
						node_set.get(0).sequence = sequence + node_set.get(0).sequence.substring(i + kh.kmer_length);
						break;
					}
				}
			}
		}
		if (node_set.get(0).sequence.length() > 3 * kh.kmer_length + kh.kmer_length) {
			for (int i = node_set.get(0).sequence.length() - 3 * kh.kmer_length; i < node_set.get(0).sequence.length()
					- kh.kmer_length; i++) {
				String kmer = node_set.get(0).sequence.substring(i, kh.kmer_length + i);
				long intval = baseOptions.kmerToIntval(kmer);
				candidates = kh.get_forward_candidates(intval, kh.kmer_hash);
				// System.out.println("forward_candidates:"+candidates.size());
				if (candidates == null) {
					continue;
				}
				boolean flag = false;
				long candidate;
				for (Map.Entry<Long, Long> mapping : candidates) {
					candidate = mapping.getKey();
					if (!has_been_used(candidate)) {
						flag = true;
						break;
					}
				}
				if (flag == true) {
					String sequence = forward_extend(intval, bifurcation, kh);
					if (sequence.length() > kh.kmer_length + kh.min_exon_length) {
						node_set.get(0).sequence = node_set.get(0).sequence.substring(0, i) + sequence;
						break;
					}
				}
			}
		}
	}

	public int add_node(Node node) {
		node_set.add(node);
		// System.out.println("节点个数：" + node_set.get(0));
		return node_set.size() - 1;

	}

	public void forward_extend_use_pairInfo(kmerHash kh, Vector<String> data, int node_index,
			Vector<Long> bifurcation) {
		// 利用paired_end_reads进行扩展 node_index为顶点标号
		int middle_read_id = data.size() / 2;
		long extend_val = 0l;
		while (true) {
			int start_pos;
			int length = node_set.get(node_index).sequence.length();
			if (length > 2 * kh.kmer_length) {
				start_pos = length - 2 * kh.kmer_length;
			} else {
				start_pos = 0;
			}
			String check = node_set.get(node_index).sequence.substring(start_pos);
			boolean extend_flag = false;
			Set<Integer> reads = new HashSet<Integer>();
			set_reads(kh, check, reads);
			if (reads.size() == 0) {
				return;
			}
			Iterator<Integer> it = reads.iterator();
			int paired_end_read_id = 0;
			String extend_str = "";
			while (it.hasNext()) {
				int read_id = it.next();
				if (kh.fr_strand == 1) { // 2-> <-1 此时应该把1反过来 根据2找1
					if (read_id < middle_read_id || (!compatible(check, data.get(read_id), kh))) {
						continue;
					} else {
						paired_end_read_id = read_id - middle_read_id;
					}
				} else if (kh.fr_strand == 2) {
					if (read_id >= middle_read_id || (!compatible(check, data.get(read_id), kh))) {
						continue;
					} else {
						paired_end_read_id = read_id + middle_read_id;
					}
				}
				String paired_end_read = load_Read.read_vector.get(paired_end_read_id);
				long max = 0l;
				long max_read_kmer = 0l;
				for (int i = 0; i < paired_end_read.length() - kh.kmer_length; i++) {
					String kmer = paired_end_read.substring(i, i + kh.kmer_length);
					long intval = baseOptions.kmerToIntval(kmer);
					if (has_been_used(intval)) {
						max = 0;
						break;
					} else if (kh.get_readset_count(kh.kmer_hash, intval) > max) {
						max = kh.get_readset_count(kh.kmer_hash, intval);
						max_read_kmer = intval;
					}
				}
				// System.out.println("entry:"+baseOptions
				// .computeEntropy(baseOptions.intvalToKmer(max_read_kmer,
				// kh.kmer_length)));
				if (max >= kh.min_seed_cov && baseOptions
						.computeEntropy(baseOptions.intvalToKmer(max_read_kmer, kh.kmer_length)) > kh.min_seed_entry) {
					long stop_kmer = baseOptions
							.kmerToIntval(node_set.get(node_index).sequence.substring(length - kh.kmer_length));
					// String str = "";
					StringBuffer str = new StringBuffer("");
					extend_val = find_head_kmer(kh, max_read_kmer, stop_kmer, str);
					// System.out.println("str:" + str);
					String extend_kmer = baseOptions.intvalToKmer(extend_val, kh.kmer_length);
					// String anchor = extend_kmer.substring(0, 5);
					String anchor = extend_kmer;
					int start = check.indexOf(anchor);
					if (start != -1) { // 如果找到了
						// if (is_similar(check.substring(start), extend_kmer,
						// 'F')) {
						// System.out.println("check.substring(start):"+check.substring(start));
						// System.out.println("extend_kmer:"+extend_kmer);
						// System.out.println("str:"+str);
						// System.out.println("check:"+check);
						// node_set.get(node_index).sequence =
						// node_set.get(node_index).sequence.substring(0,
						// start_pos + start) + str;
						// extend_flag = true;
						// }
						node_set.get(node_index).setSequence(
								node_set.get(node_index).getSequence().substring(0, start_pos + start) + str);
						extend_flag = true;
					} else {
						if (((int) str.length() > kh.pair_gap_length - 80) && (str.length() > extend_str.length()))
							extend_str = str.toString();
					}
					if (extend_flag) {
						add_used_kmers(kh, str.toString());
						break;
					}
				}

			}

			if (extend_flag) { // 那么将得到的str再向前扩展
				System.out.println("forward_extend_use_pairInfo！");
				String kmer = node_set.get(node_index).getSequence()
						.substring(node_set.get(node_index).sequence.length() - kh.kmer_length);
				long kmer_intval = baseOptions.kmerToIntval(kmer);
				String str = forward_extend(kmer_intval, bifurcation, kh);
				node_set.get(node_index).setSequence(node_set.get(node_index).sequence + str.substring(kh.kmer_length));
				// node_set.get(node_index).sequence = ;
				if (str.length() < 2 * kh.kmer_length)
					return;
			} else {

				// if (((int) extend_str.length() > kh.pair_gap_length - 80)
				// && ((int) extend_str.length() < kh.max_pair_gap_length)) {
				// node_set.get(node_index).sequence =
				// node_set.get(node_index).sequence + extend_str;
				// add_used_kmers(kh, extend_str);
				// }

				return;
			}
		}
	}

	public void reverse_extend_use_pairInfo(kmerHash kh, Vector<String> data, int node_index,
			Vector<Long> bifurcation) {
		int middle_read_id = data.size() / 2;
		while (true) {
			if (node_set.get(node_index).getSequence().length() < 2 * kh.kmer_length) {
				break;
			}
			String check = node_set.get(node_index).sequence.substring(0, 2 * kh.kmer_length);
			Set<Integer> reads = new HashSet<Integer>();
			set_reads(kh, check, reads);
			if (reads.size() == 0) {
				return;
			}
			long extend_val = 0l;
			boolean extend_flag = false;
			Iterator<Integer> it = reads.iterator();
			int paired_end_read_id = 0;
			String extend_str = "";
			while (it.hasNext()) {
				int read_id = it.next();
				if (kh.fr_strand == 1) { // --2--> ..........
											// <--1--//反向扩展，应该根据1找2，连接2
					if (read_id >= middle_read_id || (!compatible(check, data.get(read_id), kh))) {
						continue;
					} else {
						paired_end_read_id = read_id + middle_read_id;
					}
				} else if (kh.fr_strand == 2) { // --1--> .......... <--2--
					if (read_id < middle_read_id || (!compatible(check, data.get(read_id), kh))) {
						continue;
					} else {
						paired_end_read_id = read_id - middle_read_id;
					}
				}
				String paired_end_read = load_Read.read_vector.get(paired_end_read_id);
				long max = 0l;
				long max_read_kmer = 0l;
				for (int i = 0; i < paired_end_read.length() - kh.kmer_length; i++) {
					String kmer = paired_end_read.substring(i, i + kh.kmer_length);
					long intval = baseOptions.kmerToIntval(kmer);
					if (has_been_used(intval)) {
						max = 0;
						break;
					} else if (kh.get_readset_count(kh.kmer_hash, intval) > max) {
						max = kh.get_readset_count(kh.kmer_hash, intval);
						max_read_kmer = intval;
					}
				}

				if (max >= kh.min_seed_cov && baseOptions
						.computeEntropy(baseOptions.intvalToKmer(max_read_kmer, kh.kmer_length)) > kh.min_seed_entry) {
					long stop_kmer = baseOptions
							.kmerToIntval(node_set.get(node_index).sequence.substring(0, kh.kmer_length));
					StringBuffer str = new StringBuffer("");
					extend_val = find_tail_kmer(kh, max_read_kmer, stop_kmer, str);
					// System.out.println("str:" + str);
					String extend_kmer = baseOptions.intvalToKmer(extend_val, kh.kmer_length);
					String anchor = extend_kmer;
					int start = check.indexOf(anchor);
					if (start != -1) { // 如果找到了
						//if (is_similar(check.substring(0, start + 5), extend_kmer, 'R')) {
							node_set.get(node_index).sequence = str.toString()
									+ node_set.get(node_index).sequence.substring(start + 5);
							extend_flag = true;
					//	}
					} else {
						if (((int) str.length() > kh.pair_gap_length - 80) && (str.length() > extend_str.length()))
							extend_str = str.toString();
					}
					if (extend_flag) {
						add_used_kmers(kh, str.toString());
						break;
					}
				}
			}
			if (extend_flag) { // 那么将得到的str再逆向扩展
				System.out.println("reverse_extend_use_pairInfo！");
				String kmer = node_set.get(node_index).sequence.substring(0, kh.kmer_length);
				long kmer_intval = baseOptions.kmerToIntval(kmer);
				String str = reverse_extend(kmer_intval, bifurcation, kh);
				node_set.get(node_index).sequence = str.substring(0, str.length() - kh.kmer_length)
						+ node_set.get(node_index).sequence;
				if (str.length() < 2 * kh.kmer_length)
					return;
			} else {
				// 找不到就找不到，不加勉强
				// if (((int) extend_str.length() > kh.pair_gap_length - 80)
				// && ((int) extend_str.length() < kh.max_pair_gap_length)) {
				// node_set.get(node_index).sequence = extend_str +
				// node_set.get(node_index).sequence;
				// add_used_kmers(kh, extend_str);
				// }

				return;
			}

		}
	}

	public long find_tail_kmer(kmerHash kh, long seed_intval, long stop_kmer, StringBuffer str) {
		List<Map.Entry<Long, Long>> candidates = new ArrayList<Map.Entry<Long, Long>>();
		long intval = seed_intval;
		str.append(baseOptions.intvalToKmer(intval, kh.kmer_length));
		Map<Long, Boolean> use_kmers = new HashMap<Long, Boolean>();
		while (true) {
			candidates = kh.get_forward_candidates(intval, kh.kmer_hash);
			if (candidates.size() == 0) {
				break;
			}
			long candidate = 0l;
			boolean flag = false;
			for (Map.Entry<Long, Long> mapping : candidates) {
				candidate = mapping.getKey();
				if (!use_kmers.containsKey(candidate)) {
					flag = true;
					break;
				}
			}
			// 如果所有的candidate都在use_kmers中，或者candidate为stop_kmer 退出while循环
			if ((!flag) || candidate == stop_kmer) {
				break;
			}
			use_kmers.put(candidate, true);
			intval = candidate;
			int base_int = (int) (candidate & 3l);
			char base_char = baseOptions.intToBase(base_int);
			str.append(base_char);
		}
		// System.out.println("real_str:" + str);
		return intval;
	}

	public boolean compatible(String seq1, String seq2, kmerHash kh) {
		for (int i = 0; i < seq2.length() - kh.kmer_length; i++) {
			String kmer = seq2.substring(i, i + kh.kmer_length);
			int index = seq1.indexOf(kmer);
			if (index != -1) {
				if (index > i) {
					return is_aligned(seq1.substring(index - i), seq2);
				} else {
					return is_aligned(seq1, seq2.substring(i - index));
				}
			}

		}
		return false;

	}

	public boolean is_aligned(String seq1, String seq2) {
		int mismatch = 0;
		int length;
		if (seq1.length() >= seq2.length()) {
			length = seq2.length();
		} else {
			length = seq1.length();
		}
		for (int i = 0; i < length; ++i) {
			if (seq1.charAt(i) != seq2.charAt(i))
				mismatch++;
		}
		return (mismatch <= 2);
	}

	public boolean add_used_kmers(kmerHash kh, String str) {
		int length = str.length();
		if (length >= kh.kmer_length) {
			for (int i = 0; i <= length - kh.kmer_length; ++i) {
				String kmer = str.substring(i, i + kh.kmer_length);
				long intval = baseOptions.kmerToIntval(kmer);
				int cov = kh.get_readset_count(kh.kmer_hash, intval);
				used_kmers.put(intval, cov);
			}
		}

		return true;
	}

	public boolean is_similar(String seq1, String seq2, char direction) {
		int mismatch = 0;
		int length;
		if (seq1.length() >= seq2.length()) {
			length = seq2.length();
		} else {
			length = seq1.length();
		}
		if (length == 0) {
			if (seq1.length() == 0 && seq2.length() == 0) {
				return true;
			} else {
				return false;
			}
		}
		if (direction == 'F') { // 从前往后检查
			for (int i = 0; i < length; ++i) {
				if (seq1.charAt(i) != seq2.charAt(i))
					mismatch++;
			}
		} else {
			for (int i = 0; i < length; ++i) {
				if (seq1.charAt(seq1.length() - i - 1) != seq2.charAt(seq2.length() - i - 1))
					mismatch++;
			}
		}
		if ((float) mismatch / length < 0.35) {
			return true;
		} else {
			return false;
		}
	}

	public long find_head_kmer(kmerHash kh, long seed_intval, long stop_kmer, StringBuffer str) {
		List<Map.Entry<Long, Long>> candidates = new ArrayList<Map.Entry<Long, Long>>();
		long intval = seed_intval;
		str.append(baseOptions.intvalToKmer(intval, kh.kmer_length));
		Map<Long, Boolean> use_kmers = new HashMap<Long, Boolean>();
		while (true) {
			candidates = kh.get_reverse_candidates(intval, kh.kmer_hash);
			if (candidates.size() == 0) {
				break;
			}
			long candidate = 0l;
			boolean flag = false;
			for (Map.Entry<Long, Long> mapping : candidates) {
				candidate = mapping.getKey();
				if (!use_kmers.containsKey(candidate)) {
					flag = true;
					break;
				}
			}
			// 如果所有的candidate都在use_kmers中，或者candidate为stop_kmer 退出while循环
			if ((!flag) || candidate == stop_kmer) {
				break;
			}
			use_kmers.put(candidate, true);
			intval = candidate;
			int base_int = (int) ((candidate >> (kh.kmer_length * 2 - 2)) & 3l);
			char base_char = baseOptions.intToBase(base_int);
			str.insert(0, base_char);
		}
		// System.out.println("real_str:" + str);
		return intval;

	}

	public void set_reads(kmerHash kh, String seq, Set reads) {
		// 返回seq覆盖到的所有reads的编号 不重复
		if (seq.length() < kh.kmer_length) {
			return;
		}
		for (int i = 0; i < seq.length() - kh.kmer_length; i++) {
			String kmer = seq.substring(i, i + kh.kmer_length);
			long intval = baseOptions.kmerToIntval(kmer);
			Vector read_set = kh.get_readset(intval, kh.kmer_hash);
			if (read_set.size() == 0) {
				continue;
			}
			for (int j = 0; j < read_set.size(); j++) {
				reads.add(read_set.get(j));
			}
		}
	}

	public void rewrite_nodeSet(Set node_jihe) {
		// TODO Auto-generated method stub
		node_set.clear();
		Iterator<String> it = node_jihe.iterator();
		while (it.hasNext()) {
			String node = it.next();
			Node node1 = new Node();
			node1.sequence = node;
			node_set.add(node1);
		}
	}

	public void forward_check_and_extend(kmerHash kh, int node_index) {
		// System.out.println("处理第"+node_index+"个顶点"+"
		// "+node_set.get(node_index).sequence);
		int length = node_set.get(node_index).sequence.length() - kh.kmer_length;
		Vector<Long> bifurcations = new Vector<Long>();
		List<Map.Entry<Long, Long>> candidates = new ArrayList<Map.Entry<Long, Long>>();
		if (length < 0) {
			return;
		}
		for (int i = 0; i < length; i++) {
			String kmer = node_set.get(node_index).sequence.substring(i, i + kh.kmer_length);
			long intval = baseOptions.kmerToIntval(kmer);
			candidates = kh.get_forward_candidates(intval, kh.kmer_hash);
			if (candidates == null) {
				continue;
			}
			long candidate = 0l;
			for (Map.Entry<Long, Long> mapping : candidates) {
				candidate = mapping.getKey();
				if (!has_been_used(candidate)) {// 如果有未被使用的candidate
					bifurcations.add(intval);
					// System.out.println(bifurcations);
					break;
				}
			}

		}
		// System.out.println("初始的顶点序列："+node_set.get(0).getSequence());
		// for(int i=0;i<bifurcations.size();i++){
		// System.out.print(baseOptions.intvalToKmer(bifurcations.get(i),
		// kh.kmer_length)+" ");
		// }
		// System.out.println();
		// System.out.println("此时的分叉点"+bifurcations);
		grow_and_branch(kh, node_index, bifurcations);
	}

	public void reverse_check_and_extend(kmerHash kh, int node_index) {
		int length = node_set.get(node_index).sequence.length() - kh.kmer_length;
		Vector<Long> bifurcations = new Vector<Long>();
		List<Map.Entry<Long, Long>> candidates = new ArrayList<Map.Entry<Long, Long>>();
		if (length < 0) {
			return;
		}
		for (int i = 0; i < length; i++) {
			String kmer = node_set.get(node_index).sequence.substring(i, i + kh.kmer_length);
			long intval = baseOptions.kmerToIntval(kmer);
			candidates = kh.get_reverse_candidates(intval, kh.kmer_hash);
			if (candidates == null || candidates.size() <= 1) {
				continue;
			}
			long candidate = 0l;
			for (Map.Entry<Long, Long> mapping : candidates) {
				candidate = mapping.getKey();
				if (!has_been_used(candidate)) {// 如果有未被使用的candidate
					bifurcations.add(intval);
					// System.out.println(bifurcations);
					break;
				}
			}

		}
		// grow_and_branch(kh, node_index, bifurcations);
		for (int i = bifurcations.size(); i > 0; i--) {
			// System.out.println("顶点" + node_index +
			// "进行了reverse_check_and_extend!");
			long intval = bifurcations.get(i - 1);
			// System.out.println("分叉点为："+baseOptions.intvalToKmer(intval,
			// kh.kmer_length));
			Vector<Long> bifurcations_more = new Vector<Long>();
			String str = reverse_extend(intval, bifurcations_more, kh);
			String endkmer = str.substring(0, kh.kmer_length);
			long end_val = baseOptions.kmerToIntval(endkmer);
			List<Map.Entry<Long, Long>> candidates1 = new ArrayList<Map.Entry<Long, Long>>();
			candidates1 = kh.get_reverse_candidates(end_val, kh.kmer_hash);

			if (candidates1.size() > 0) { // add bubble (possible)

				int bubble_val = add_reverse_bubble(node_index, str, kh);
				if (bubble_val > 0) {
					reverse_branches.add(bubble_val);
					// System.out.println("reverse_add_bubble成功！");
				}
			} else { // add branch

				int bubble_val = add_reverse_branch(node_index, str, kh);
				if (bubble_val > 0)
					// System.out.println("reverse_add_branch成功！");

					reverse_branches.add(bubble_val);
			} // else
		} // for
	}

	public void reverse_check_and_extend(kmerHash kh) {
		for (int i = 0; i < node_set.size(); i++) {
			// if(forward_branches.contains(i)){
			// continue;
			// }
			reverse_check_and_extend(kh, i);
		}
	}

	int add_reverse_branch(int node_p, String str, kmerHash kh) {

		int str_length = str.length() - kh.kmer_length;
		if (str_length >= kh.min_exon_length) {

			String end_kmer = str.substring(str.length() - kh.kmer_length, str.length());
			int start = node_set.get(node_p).sequence.indexOf(end_kmer);
			if (start == -1) {
				restore_kmers(str.substring(0, str.length() - 1), kh);
				return -1;
			}
			// if (kh.is_paired_end) {
			// long count = kh.get_readset_count(kh.kmer_hash,
			// baseOptions.kmerToIntval(end_kmer));
			// boolean is_branch = check_reverse_branch_with_pair_info(kh, str,
			// count);
			// if (!is_branch) {
			// restore_kmers(str.substring(0, str.length() - 1), kh);
			// return -1;
			// }
			// }
			// Node node1,node2;
			Node node1 = new Node();
			Node node2 = new Node();
			node1.sequence = node_set.get(node_p).sequence.substring(start);
			node2.sequence = str.substring(0, str_length);
			if (start == 0) { // possible if not node 0
				int node2_index = add_node(node2);
				node_set.get(node2_index).addChildren(node_p); // node1 is p
				reverse_branches.add(node2_index);
				return node2_index;
			}

			if (is_similar(node_set.get(node_p).sequence.substring(0, start), node2.sequence, 'R')) {
				if (str_length > (int) start && node_p == 0) {
					node_set.get(node_p)
							.setSequence(str.substring(0, str_length) + node_set.get(node_p).sequence.substring(start));
					return -2;
				} else if (start > kh.kmer_length) {
					return -1;
				}
			}

			if (node_p == 0 && start < 2 * kh.kmer_length) {
				node_set.get(node_p)
						.setSequence(str.substring(0, str_length) + node_set.get(node_p).sequence.substring(start));
				return -2;
			}

			node_set.get(node_p).setSequence(node_set.get(node_p).sequence.substring(0, start));
			int node1_index = add_node(node1);
			for (int i = 0; i < node_set.get(node_p).children.size(); i++) {
				node_set.get(node1_index).addChildren(node_set.get(node_p).children.get(i));
			}
			// node_set.get(node1_index).children =
			// node_set.get(node_p).children;
			node_set.get(node_p).children.clear();
			node_set.get(node_p).addChildren(node1_index);
			int node2_index = add_node(node2);
			node_set.get(node2_index).addChildren(node1_index);
			reverse_branches.add(node2_index);

			return node2_index;

		} else {
			return -1;
		}
	}

	public int add_reverse_bubble(int node_p, String str, kmerHash kh) {

		if (str.length() < 2 * kh.min_anchor_length) {
			return -1;
		}
		// 可能报错！
		String anchor_right = str.substring(str.length() - kh.kmer_length,
				str.length() - kh.kmer_length + kh.kmer_length);
		int start = node_set.get(node_p).sequence.indexOf(anchor_right);
		if (start == -1) // possible if overlap
		{
			return -1;
		}

		String anchor_left = str.substring(0, kh.kmer_length - 1);
		int node_q = -1;
		while (anchor_left.length() > kh.min_anchor_length) {
			node_q = find_node_index(anchor_left);
			if (node_q >= 0)
				break;
			anchor_left = anchor_left.substring(1);
		}
		Set<Integer> checked = new HashSet<Integer>();
		if (node_q < 0 || node_p == node_q || has_path(node_p, node_q, checked)) {
			return -1;
		}

		int anchor_length = anchor_left.length();

		int end = node_set.get(node_q).sequence.indexOf(anchor_left);
		if (anchor_length + end == node_set.get(node_q).sequence.length()) {
			int length = str.length() - anchor_length - kh.kmer_length;
			if (length > 0 && start > 0 && length + start < 10) {
				return -1;
			}

			Node node1 = new Node();
			int node1_index = -1;
			// node_idx_t q1 = -1;
			if (length > 0) {
				node1.sequence = str.substring(anchor_length, length + anchor_length);
				node1_index = add_node(node1);
				reverse_branches.add(node1_index);
				node_set.get(node_q).addChildren(node1_index);
			} else {
				node1_index = node_q;
			}

			if (start == 0) {
				// node_set_[q1].add_child(p);
				node_set.get(node1_index).addChildren(node_p);
			} else {
				Node node2 = new Node();
				node2.sequence = node_set.get(node_p).sequence.substring(start);
				int node2_index = add_node(node2);
				node_set.get(node1_index).addChildren(node2_index);
				for (int i = 0; i < node_set.get(node_p).children.size(); i++) {
					node_set.get(node2_index).addChildren(node_set.get(node_p).children.get(i));
				}
				// node_set.get(node2_index).children =
				// node_set.get(node_p).children;
				node_set.get(node_p).children.clear();
				node_set.get(node_p).addChildren(node2_index);
				node_set.get(node_p).setSequence(node_set.get(node_p).sequence.substring(0, start));

			}

			if (node1_index != node_q)
				return node1_index;
		}

		return -1;
	}

	// long intvalg;
	public void grow_and_branch(kmerHash kh, int node_index, Vector<Long> bifurcations) {

		while (bifurcations.size() > 0) {
			long intval = bifurcations.lastElement();
			// System.out.println("处理分叉点："+intval);
			bifurcations.remove(bifurcations.size() - 1);
			Vector<Long> bifurcations_more = new Vector<Long>();
			String str = forward_extend(intval, bifurcations_more, kh);
			// System.out.println("分叉点："+intval+"
			// "+baseOptions.intvalToKmer(intval, kh.kmer_length));
			if (bifurcations_more.size() > 0 && bifurcations_more.get(0) == intval) {
				// 如果intval还有两个或以上没用的candidate
				bifurcations.add(intval);
				bifurcations_more.remove(0);
				// bifurcations.addAll(bifurcations_more);
			}
			// System.out.println(bifurcations);
			String end_kmer = str.substring(str.length() - kh.kmer_length); // str的最后一个kmer
			long end_intval = baseOptions.kmerToIntval(end_kmer);
			List<Map.Entry<Long, Long>> candidates = new ArrayList<Map.Entry<Long, Long>>();
			candidates = kh.get_forward_candidates(end_intval, kh.kmer_hash);
			if (candidates.size() > 0) {
				// intvalg=intval;
				// System.out.println("此时的顶点1："+node_set.get(node_index).getSequence());
				int bubble_val = add_bubble(node_index, str, kh);
				// System.out.println(bubble_val);
				// if(bubble_val==2){
				// System.out.println(str);
				// }
				// if(bubble_val==-6){
				// System.out.println("此时的顶点序列："+node_set.get(node_index).getSequence());
				// for(int i=0;i<bifurcations.size();i++){
				// System.out.print(baseOptions.intvalToKmer(bifurcations.get(i),
				// kh.kmer_length)+" ");
				// }
				// System.out.println();
				// }
				if (bubble_val > 0) {
					// System.out.println("forward_add_bubble成功！");
					forward_branches.add(bubble_val);
				}
				if (bubble_val == -2) {
					bubble_val = add_branch(node_index, str, kh);
					if (bubble_val > 0) {
						// System.out.println("forward_add_branch成功！");
						forward_branches.add(bubble_val);
						node_set.get(bubble_val)
								.setSequence(node_set.get(bubble_val).sequence.substring(kh.kmer_length));
					}
				}
			} else {
				int bubble_val = add_branch(node_index, str, kh);
				if (bubble_val > 0) {
					// System.out.println("forward_add_branch成功！");

					forward_branches.add(bubble_val);
					node_set.get(bubble_val).setSequence(node_set.get(bubble_val).sequence.substring(kh.kmer_length));
				}
			}

		} // while
	}

	public int add_bubble(int node_p, String str, kmerHash kh) {

		if (str.length() < 2 * kh.min_anchor_length) {
			return -1;
		}
		String anchor_left = str.substring(0, kh.kmer_length);
		int start = node_set.get(node_p).sequence.indexOf(anchor_left);
		// 如果找不到则不断缩小anchor_left
		if (start == -1) {
			while (anchor_left.length() > kh.min_anchor_length) {
				anchor_left = anchor_left.substring(0, anchor_left.length() - 1);
				start = node_set.get(node_p).sequence.indexOf(anchor_left);
				if (start != -1) {
					break; // 直至找到为止
				}
			}
			if (start == -1) {
				// System.out.println("node_p:"+node_set.get(node_p).getSequence());
				// System.out.println("str:"+str);
				// System.out.println("anchor_left:"+anchor_left);
				// System.out.println(baseOptions.intvalToKmer(intvalg,
				// kh.kmer_length));
				// System.out.println(node_set.get(node_p).getSequence());
				System.out.println("无法找到起始anchor！");
				return -6;
			}
		}

		String anchor_right = str.substring(str.length() - kh.kmer_length + 1);
		int node_q = -1;
		while (anchor_right.length() > kh.min_anchor_length) {
			node_q = find_node_index(anchor_right);
			if (node_q >= 0) {
				break;
			}
			anchor_right = anchor_right.substring(1);
		}
		if (node_q == -1) {
			restore_kmers(str, kh);
			return -2;
		}

		int anchor_left_length = anchor_left.length();
		int anchor_right_length = anchor_right.length();
		int end = node_set.get(node_q).sequence.indexOf(anchor_right);

		// 根据连接处的cov判断是否可以连接
		// System.out.println(0.5 * anchor_left_length+kh.kmer_length);
		// System.out.println(str.length());
		// if((0.5 * anchor_left_length + kh.kmer_length)>kh.kmer_length){
		// return -1;
		// }
		String jun1 = str.substring((int) (0.5 * anchor_left_length),
				(int) (0.5 * anchor_left_length + kh.kmer_length));
		long jun1_cov = kh.get_readset_count(kh.kmer_hash, baseOptions.kmerToIntval(jun1));
		// System.out.println("arl:"+anchor_right_length);
		// System.out.println(str.length() - kh.kmer_length );
		String jun2 = str.substring(str.length() - kh.kmer_length - (int) (0.5 * anchor_right_length),
				str.length() - kh.kmer_length - (int) (0.5 * anchor_right_length) + kh.kmer_length);
		long jun2_cov = kh.get_readset_count(kh.kmer_hash, baseOptions.kmerToIntval(jun2));
		// System.out.println(jun2_cov+" "+jun2);
		if (jun1_cov < kh.min_function_count || jun2_cov < kh.min_function_count) {
			System.out.println("连接点cov太小，无法连接！");
			return -1;
		}

		int length = str.length() - anchor_left_length - anchor_right_length; // str有多少未被使用

		if (node_p == node_q) { // 顶点本身回归
			// TODO Auto-generated method stub 可修改
			if (end <= start) {
				return -1;
			}
			int distance = end - start - anchor_left_length; // 表示end与start之间是否比anchor_left大
			if (length + distance < 4) {
				return -1;
			}
			if ((distance == 0 && length < kh.kmer_length) || (length == 0 && distance < kh.kmer_length)) {
				return -1;
			}
			if (distance > 0 && length == distance
					&& is_similar(
							node_set.get(node_p).sequence.substring(start + anchor_left_length,
									start + anchor_left_length + distance),
					str.substring(anchor_left_length, anchor_left_length + length), 'F')) {
				// 重合
				return -1;
			}
			if (length <= 0 && distance <= 0) {
				return -1;
			} else if (length < 0) {
				anchor_left_length = anchor_left_length + length;
			} else if (distance < 0) {
				anchor_left_length = anchor_left_length + distance; // 设置
																	// anchor_left与anchor_right刚好完全占据str
			}
			if ((int) start + anchor_left_length <= kh.kmer_length) // necessary
			{
				return -1;
			}
			Node node1 = new Node();
			Node node2 = new Node();
			// System.out.println("&&&&&&&:"+node_set.get(node_p).sequence);
			node1.sequence = node_set.get(node_p).sequence.substring(end);

			int node1_index = add_node(node1);
			// System.out.println("&&&&&&&:"+node_set.get(node1_index).sequence);
			// System.out.println("p的孩子："+node_set.get(node_p).children);
			for (int i = 0; i < node_set.get(node_p).children.size(); i++) {
				node_set.get(node1_index).addChildren(node_set.get(node_p).children.get(i));
			}
			// System.out.println("node1的孩子："+node_set.get(node1_index).children);
			node_set.get(node_p).children.clear();
			int node2_index = -1;
			if (distance <= 0) { // 此时length>0
				node_set.get(node_p)
						.setSequence(node_set.get(node_p).sequence.substring(0, start + anchor_left_length));
				node_set.get(node_p).addChildren(node1_index);
				// System.out.println(node_set.get(node_p).sequence);
			} else {
				Node node3 = new Node();
				if (length < 0) {
					node3.sequence = node_set.get(node_p).sequence.substring(start + anchor_left_length,
							start + anchor_left_length + distance - length);
				} else {
					node3.sequence = node_set.get(node_p).sequence.substring(start + anchor_left_length,
							start + anchor_left_length + distance);
				}
				int node3_index = add_node(node3);
				node_set.get(node3_index).addChildren(node1_index);
				node_set.get(node_p)
						.setSequence(node_set.get(node_p).sequence.substring(0, start + anchor_left_length));
				// node_set.get(node_p).sequence.substring(0, start +
				// anchor_left_length);
				// System.out.println(node_set.get(node_p).sequence);
				node_set.get(node_p).addChildren(node3_index);
			}
			if (length > 0) {
				if (distance < 0) {
					node2.sequence = str.substring(anchor_left_length, anchor_left_length + length - distance);
				} else {
					node2.sequence = str.substring(anchor_left_length, anchor_left_length + length);
				}
				node2_index = add_node(node2);
				node_set.get(node2_index).addChildren(node1_index);
				node_set.get(node_p).addChildren(node2_index);
			} else {
				node_set.get(node_p).addChildren(node1_index);
			}
			// System.out.println("node1的孩子："+node_set.get(node1_index).children);

			return node2_index; // 新增顶点标号
		} else {
			Set<Integer> checked = new HashSet<Integer>();
			if (node_q == 0 || has_path(node_p, node_q, checked)) {
				return -1;
			}
			if (node_set.get(node_p).isChild(node_q) && length > 0
					&& node_set.get(node_p).sequence.length() - start + end < kh.kmer_length) {
				return -1;
			}
			Node node1 = new Node();
			Node node2 = new Node();
			int node1_index = -1;
			int node2_index = -1;
			if (length < 0) {
				anchor_left_length = anchor_left_length + length;
			}
			if (start + anchor_left_length <= kh.kmer_length) {
				return -1;
			}
			if (start + anchor_left_length < node_set.get(node_p).sequence.length()) {
				// 只能小于或等于 如果等于，无需再进行操作
				node1.sequence = node_set.get(node_p).sequence.substring(start + anchor_left_length);
				node_set.get(node_p)
						.setSequence(node_set.get(node_p).sequence.substring(0, start + anchor_left_length));
				// node_set.get(node_p).sequence.substring(0, start +
				// anchor_left_length);
				node1_index = add_node(node1);
				for (int i = 0; i < node_set.get(node_p).children.size(); i++) {
					node_set.get(node1_index).addChildren(node_set.get(node_p).children.get(i));
				}
				// node_set.get(node1_index).children =
				// node_set.get(node_p).children;
				node_set.get(node_p).children.clear();
				node_set.get(node_p).addChildren(node1_index);
			}
			if (length > 0) {
				node2.sequence = str.substring(anchor_left_length, anchor_left_length + length);
				node2_index = add_node(node2);
				node_set.get(node_p).addChildren(node2_index);
			} else {
				node2_index = node_p; // ?
			}
			if (end == 0) {
				node_set.get(node2_index).addChildren(node_q);
			} else {
				Node node3 = new Node();
				if (node_set.get(node_q).sequence.length() >= end) {
					node3.sequence = node_set.get(node_q).sequence.substring(end);
					int node3_index = add_node(node3);
					for (int i = 0; i < node_set.get(node_q).children.size(); i++) {
						node_set.get(node3_index).addChildren(node_set.get(node_q).children.get(i));
					}
					// node_set.get(node3_index).children =
					// node_set.get(node_q).children;
					node_set.get(node_q).children.clear();
					node_set.get(node_q).addChildren(node3_index);
					node_set.get(node2_index).addChildren(node3_index);
					node_set.get(node_q).sequence = node_set.get(node_q).sequence.substring(0, end);
				}
			}
			// System.out.println("add_bubble之后的节点为：");
			// System.out.println(node2_index+" "+"node2"+"
			// "+node_set.get(node2_index).sequence);
			// System.out.println(node_p+" "+"p"+"
			// "+node_set.get(node_p).sequence);
			// System.out.println(node1_index+" "+"node1"+"
			// "+node_set.get(node1_index).sequence);

			if (node2_index == node_p) {
				return -1;
			} else {
				return node2_index;
			}
		}

	}

	public boolean has_path(int node_p, int node_q, Set checked) {
		// 判断q到p是否有路径
		boolean flag = false;
		if (node_set.get(node_q).children.size() == 0 || checked.contains(node_q)) {
			return false;
		} else {
			checked.add(node_q);
		}
		for (int i = 0; i < node_set.get(node_q).children.size(); i++) {
			if (node_set.get(node_q).children.get(i) == node_p) {
				flag = true;
				break;
			} else {
				flag = has_path(node_p, node_set.get(node_q).children.get(i), checked);
				if (flag) {
					break;
				}
			}
		}
		return flag;
	}

	public int find_node_index(String anchor) {
		int idx = -1;
		for (int i = 0; i < node_set.size(); ++i) {
			if (node_set.get(i).sequence.indexOf(anchor) != -1) {
				idx = i;
				break;
			}
		}
		return idx;
	}

	public void restore_kmers(String str, kmerHash kh) {

		if (str.length() < kh.kmer_length) {
			return;
		}

		for (int i = 0; i <= str.length() - kh.kmer_length; i++) {
			String kmer = str.substring(i, i + kh.kmer_length);
			long intval = baseOptions.kmerToIntval(kmer);
			restore_kmers.add(intval);
		}

	}

	// 如果在扩展过程中找不到对应的q 则回不到主干，只能作为分支
	public int add_branch(int node_p, String branch, kmerHash kh) {
		if (node_set.get(node_p).sequence.length() <= kh.kmer_length) {
			return -1;
		}
		if (branch.length() > kh.kmer_length + kh.min_exon_length) {
			String start_kmer = branch.substring(0, kh.kmer_length);
			int start = node_set.get(node_p).sequence.indexOf(start_kmer);
			if (start == -1) {
				restore_kmers(branch.substring(1), kh);
				return -1;
			}
			if (node_p == 0 && start < kh.pair_gap_length) {
				restore_kmers(branch.substring(1), kh);
				return -1;
			}
			if (kh.is_paired_end) {

				long count = kh.get_readset_count(kh.kmer_hash, baseOptions.kmerToIntval(start_kmer));
				boolean is_branch = check_forward_branch_with_pair_info(kh, branch, count);

				if (!is_branch) {
					restore_kmers(branch.substring(1), kh);
					return -1;
				}
			}
			Node node1 = new Node();
			Node node2 = new Node();
			if (start + kh.kmer_length < node_set.get(node_p).sequence.length()) {
				// 需要拆分
				node1.sequence = node_set.get(node_p).sequence.substring(start + kh.kmer_length);
			}
			node2.sequence = branch;
			if (is_similar(branch.substring(kh.kmer_length), node1.sequence, 'F')) {
				if (node_set.get(node_p).children.size() == 0 && node1.sequence.length() < kh.min_exon_length) {
					return -2;
				} else if (node1.sequence.length() + kh.kmer_length > branch.length()) {
					return -1;
				}
			}
			if (node_set.size() == 1 && node1.sequence.length() < 2 * kh.kmer_length) {
				node_set.get(node_p).setSequence(node_set.get(node_p).sequence.substring(0, start) + branch);
				return -2;
			}
			if (node1.sequence.length() != 0) {
				int node1_index = add_node(node1);
				for (int i = 0; i < node_set.get(node_p).children.size(); i++) {
					node_set.get(node1_index).addChildren(node_set.get(node_p).children.get(i));
				}
				// node_set.get(node1_index).children =
				// node_set.get(node_p).children;
				node_set.get(node_p).children.clear();
				node_set.get(node_p).addChildren(node1_index);
			}
			int node2_index = add_node(node2);
			node_set.get(node_p).addChildren(node2_index);
			node_set.get(node_p).setSequence(node_set.get(node_p).sequence.substring(0, start + kh.kmer_length));
			return node2_index;
		} else {
			return -1;
		}

	}

	public boolean check_forward_branch_with_pair_info(kmerHash kh, String branch, long count) {
		// 计算左边的paired_end_read
		boolean is_branch = false;
		if (3 * kh.kmer_length > branch.length()) {
			return false;
		}
		String check = branch.substring(kh.kmer_length, 3 * kh.kmer_length);
		int middle_read_id = load_Read.read_vector.size() / 2;
		Set<Integer> reads = new HashSet<Integer>();
		set_reads(kh, check, reads);
		Set<Long> kmers_in_branch = new HashSet<Long>();
		int support = 0;
		int unsupport = 0;
		Iterator<Integer> it = reads.iterator();
		int paired_end_read_id = 0;
		String extend_str = "";
		while (it.hasNext()) {
			int read_id = it.next();
			if (kh.fr_strand == 1) { // 2-> <-1 此时应该把1反过来 根据2找1
				if (read_id < middle_read_id) {
					String mate_left = load_Read.read_vector.get(read_id);
					if (compatible(branch, mate_left, kh)) {
						String mate_right = load_Read.read_vector.get(read_id + middle_read_id);
						int used = 0;
						for (int i = 0; i <= mate_right.length() - kh.kmer_length; i++) {
							long intval = baseOptions.kmerToIntval(mate_right.substring(i, i + kh.kmer_length));
							if (has_been_used(intval))
								used++;
						}
						if (used == 0) {
							unsupport++;
							if ((unsupport >= kh.min_ratio_welds * count) && (unsupport >= kh.min_reads_span_junction))
								break;
						} else if (used + 5 >= mate_right.length() - kh.kmer_length) {
							support++;
							if ((support >= kh.min_ratio_welds * count) && (support >= kh.min_reads_span_junction)) {
								is_branch = true;
								break;
							}
						}
					}
				}
			} else if (kh.fr_strand == 2) {// ->1 <-2
				if (read_id >= middle_read_id) {
					String mate_right = load_Read.read_vector.get(read_id);
					if (compatible(branch, mate_right, kh)) {
						String mate_left = load_Read.read_vector.get(read_id - middle_read_id);
						int used = 0;
						for (int i = 0; i <= (int) mate_left.length() - kh.kmer_length; ++i) {
							long intval = baseOptions.kmerToIntval(mate_left.substring(i, i + kh.kmer_length));
							if (has_been_used(intval))
								used++;
						}
						if (used == 0) {
							unsupport++;
							if ((unsupport >= kh.min_ratio_welds * count) && (unsupport >= kh.min_reads_span_junction))
								break;
						} else if (used + 5 >= (int) mate_left.length() - kh.kmer_length) {
							support++;
							if ((support >= kh.min_ratio_welds * count) && (support >= kh.min_reads_span_junction)) {
								is_branch = true;
								break;
							}
						}
					}
				}
			}

		} // else
		return is_branch;
	}

	public boolean check_reverse_branch_with_pair_info(kmerHash kh, String branch, long count) {
		// 计算右边的paired_end_read
		boolean is_branch = false;
		int middle_read_id = load_Read.read_vector.size() / 2;
		int str_length = branch.length() - kh.kmer_length;
		Set<Integer> reads = new HashSet<Integer>();
		int pos = str_length > 2 * kh.kmer_length ? (str_length - 2 * kh.kmer_length) : 0;
		set_reads(kh, branch.substring(pos, 2 * kh.kmer_length + pos), reads);

		Set<Long> kmers_in_branch = new HashSet<Long>();

		int support = 0;
		int unsupport = 0;

		Iterator<Integer> it = reads.iterator();
		int paired_end_read_id = 0;
		String extend_str = "";
		while (it.hasNext()) {
			int read_id = it.next();
			if (kh.fr_strand == 1) {
				if (read_id >= middle_read_id) {
					String mate_right = load_Read.read_vector.get(read_id);
					if (compatible(branch, mate_right, kh)) {
						String mate_left = load_Read.read_vector.get(read_id - middle_read_id);
						int used = 0;
						for (int i = 0; i < mate_left.length() - kh.kmer_length; i++) {
							String kmer = mate_left.substring(i, i + kh.kmer_length);
							long intval = baseOptions.kmerToIntval(kmer);
							if (has_been_used(intval)) {
								used++;
							}
						}
						if (used == 0) {
							unsupport++;
							if ((unsupport > kh.min_ratio_welds * count) && (unsupport >= kh.min_reads_span_junction)) {
								break;
							}
						} else if (used + 5 >= mate_left.length() - kh.kmer_length) {
							support++;
							if ((support >= kh.min_ratio_welds) && (support >= kh.min_reads_span_junction)) {
								is_branch = true;
								break;
							}
						}
					}
				}
			} else if (kh.fr_strand == 2) {
				if (read_id < middle_read_id) {
					String mate_left = load_Read.read_vector.get(read_id);
					if (compatible(branch, mate_left, kh)) {
						String mate_right = load_Read.read_vector.get(read_id + middle_read_id);
						int used = 0;
						for (int i = 0; i < mate_right.length() - kh.kmer_length; i++) {
							String kmer = mate_left.substring(i, i + kh.kmer_length);
							long intval = baseOptions.kmerToIntval(kmer);
							if (has_been_used(intval)) {
								used++;
							}
						}
						if (used == 0) {
							unsupport++;
							if ((unsupport > kh.min_ratio_welds * count) && (unsupport >= kh.min_reads_span_junction)) {
								break;
							}
						} else if (used + 5 >= mate_left.length() - kh.kmer_length) {
							support++;
							if ((support >= kh.min_ratio_welds) && (support >= kh.min_reads_span_junction)) {
								is_branch = true;
								break;
							}
						}
					}
				}
			}

		}

		return is_branch;
	}

	public boolean if_can_extend(kmerHash kh) {
		List<Map.Entry<Long, Long>> candidates = new ArrayList<Map.Entry<Long, Long>>();
		// boolean flag=false;
		for (int i = 0; i < node_set.size(); i++) {
			// 对于每一个顶点，检查其是否还有未使用的candidates
			// forward
			if (forward_branches.contains(node_set.get(i)) || reverse_branches.contains(node_set.get(i))) {
				for (int j = 0; j < node_set.get(i).getSequence().length(); j++) {
					String kmer = node_set.get(i).getSequence().substring(j, j + kh.kmer_length);
					long intval = baseOptions.kmerToIntval(kmer);
					candidates = kh.get_forward_candidates(intval, kh.kmer_hash);
					if (candidates.size() == 0) {
						break;
					}
					long candidate = 0l;
					for (Map.Entry<Long, Long> mapping : candidates) {
						candidate = mapping.getKey();
						if (!has_been_used(candidate)) {// 如果有未被使用的candidate
							// System.out.println("未被使用");
							// flag = true;
							// cov = mapping.getValue();
							// break;
							System.out.println("将进行再扩展！");
							return true;
						}
						// count++;
						// System.out.println();
					}
				}
				// reverse
				for (int j = 0; j < node_set.get(i).getSequence().length(); j++) {
					String kmer = node_set.get(i).getSequence().substring(j, j + kh.kmer_length);
					long intval = baseOptions.kmerToIntval(kmer);
					candidates = kh.get_reverse_candidates(intval, kh.kmer_hash);
					if (candidates.size() == 0) {
						break;
					}
					long candidate = 0l;
					for (Map.Entry<Long, Long> mapping : candidates) {
						candidate = mapping.getKey();
						if (!has_been_used(candidate)) {// 如果有未被使用的candidate
							// System.out.println("未被使用");
							// flag = true;
							// cov = mapping.getValue();
							// break;
							System.out.println("将进行再扩展！");
							return true;
						}
						// count++;
						// System.out.println();
					}
				}
			}
		}

		return false;
	}

	public boolean init_trunk(kmerHash kh, long seed_val, Set node_jihe, SplicingGraph splicing_graph) {
		// TODO Auto-generated method stub
		used_kmers.put(seed_val, (int) kh.kmer_map.get(seed_val));
		Vector<Long> bifurcation = new Vector<Long>();
		// String right = forward_extend(kh.list.get(0).getKey(), bifurcation,
		// kh);
		String right = forward_extend(seed_val, bifurcation, kh);
		String left = reverse_extend(seed_val, bifurcation, kh);
		String trunk = left + right.substring(kh.kmer_length);
		Node trunkN = new Node();
		trunkN.sequence = trunk;
		if (trunkN.sequence.length() <= 200) {
			return false;
		}
		// System.out.println("trunk："+trunk);
		// SplicingGraph splicing_graph = new SplicingGraph();
		int p = splicing_graph.add_node(trunkN);
		// System.out.println("******************" + p);
		// // System.out.println("原序列："+splicing_graph.node_set.size());
		// System.out.println("原序列" + " size:" +
		// splicing_graph.node_set.get(0).sequence.length() + " "
		// + splicing_graph.node_set.get(0).sequence);
		splicing_graph.extend_again(kh, bifurcation);
		// System.out.println("after
		// extend_again:"+node_set.get(0).getSequence());
		// System.out.println("第一次扩展:" + " size:" +
		// splicing_graph.node_set.get(0).sequence.length() + " "
		// + splicing_graph.node_set.get(0).sequence);
		// System.out.println();
		splicing_graph.forward_extend_use_pairInfo(kh, load_Read.read_vector, p, bifurcation);

		// splicing_graph.node_set.get(0).sequence.length() + " "
		// + splicing_graph.node_set.get(0).sequence);
		splicing_graph.reverse_extend_use_pairInfo(kh, load_Read.read_vector, p, bifurcation);
		// System.out.println("paired——end扩展："+node_set.get(0).getSequence());
		node_jihe.add(splicing_graph.node_set.get(p).sequence);
		// System.out.println(bifurcation);
		// System.out.println();
		// System.out.println(bifurcation);
		return true;
	}

	public void init_parents() {

		for (int i = 0; i < node_set.size(); i++) {
			Vector<Integer> children = new Vector<Integer>();
			children = node_set.get(i).getChildren();
			for (int j = 0; j < children.size(); j++) {
				node_set.get(children.get(j)).addParents(i);
			}

		}

	}

	/*
	 * public void compute_node_cov(kmerHash kh) {
	 * 
	 * for (int i = 0; i < node_set.size(); i++) { int cov = 0; for (int j = 0;
	 * j <= node_set.get(i).getSequence().length() - kh.kmer_length; j++) {
	 * String kmer = node_set.get(i).getSequence().substring(j, j +
	 * kh.kmer_length); long intval = baseOptions.kmerToIntval(kmer); if
	 * (kh.kmer_map.containsKey(intval)) { cov += kh.kmer_map.get(intval); }
	 * else { cov += 0; } } node_set.get(i).setcov(cov); } }
	 * 
	 * public void compute_edge_cov(kmerHash kh, int[][] edges) { int father =
	 * 0; int child = 0; for (int i = 0; i < node_set.size(); i++) { father = i;
	 * if (node_set.get(i).getChildren().size() == 0) { continue; } else {
	 * Vector<Integer> edge_kmer=new Vector<Integer>(); for (int j = 0; j <
	 * node_set.get(i).getChildren().size(); j++) { child = (int)
	 * node_set.get(i).getChildren().get(j); int edge_cov = 0; if
	 * (node_set.get(i).getSequence().length() >= kh.kmer_length &&
	 * node_set.get(child).getSequence().length() >= kh.kmer_length) { // int
	 * kmer_length=2*(kh.kmer_length); String edge_str =
	 * node_set.get(i).getSequence()
	 * .substring(node_set.get(i).getSequence().length() - kh.kmer_length + 1) +
	 * node_set.get(child).getSequence().substring(0, kh.kmer_length - 1);
	 * for(int k=0;k<=edge_str.length()-kh.kmer_length;k++){ String
	 * kmer=edge_str.substring(k,k+kh.kmer_length);
	 * if(!kh.kmer_map.containsKey(baseOptions.kmerToIntval(kmer))){
	 * edge_kmer.add(0); }else{ int
	 * cov1=kh.kmer_map.get(baseOptions.kmerToIntval(kmer));
	 * edge_kmer.add(cov1);} } edges[i][child] = compute_edge_bv(edge_kmer,
	 * edge_cov); } else if (node_set.get(i).getSequence().length() <
	 * kh.kmer_length && node_set.get(child).getSequence().length() >=
	 * kh.kmer_length) { String edge_str = node_set.get(i).getSequence() +
	 * node_set.get(child).getSequence().substring(0, kh.kmer_length - 1);
	 * for(int k=0;k<=edge_str.length()-kh.kmer_length;k++){ String
	 * kmer=edge_str.substring(k,k+kh.kmer_length); int
	 * cov1=kh.kmer_map.get(baseOptions.kmerToIntval(kmer));
	 * edge_kmer.add(cov1); } edges[i][child] = compute_edge_bv(edge_kmer,
	 * edge_cov); } else if (node_set.get(i).getSequence().length() >=
	 * kh.kmer_length && node_set.get(child).getSequence().length() <
	 * kh.kmer_length) { String edge_str = node_set.get(i).getSequence()
	 * .substring(node_set.get(i).getSequence().length() - kh.kmer_length + 1) +
	 * node_set.get(child).getSequence(); for(int
	 * k=0;k<=edge_str.length()-kh.kmer_length;k++){ String
	 * kmer=edge_str.substring(k,k+kh.kmer_length); int
	 * cov1=kh.kmer_map.get(baseOptions.kmerToIntval(kmer));
	 * edge_kmer.add(cov1); } edges[i][child] = compute_edge_bv(edge_kmer,
	 * edge_cov); } else { String edge_str = node_set.get(i).getSequence() +
	 * node_set.get(child).getSequence(); for(int
	 * k=0;k<=edge_str.length()-kh.kmer_length;k++){ String
	 * kmer=edge_str.substring(k,k+kh.kmer_length); int
	 * cov1=kh.kmer_map.get(baseOptions.kmerToIntval(kmer));
	 * edge_kmer.add(cov1); } edges[i][child] = compute_edge_bv(edge_kmer,
	 * edge_cov); }
	 * 
	 * } } }
	 * 
	 * }
	 * 
	 * public void compute_node(int node_index, kmerHash kh, int[][] edges) {
	 * Vector<Integer> node_kmer_count = new Vector<Integer>(); //
	 * Vector<Integer> node_kmer_count_asc = new Vector<Integer>(); int length =
	 * node_set.get(node_index).getSequence().length(); // 将顶点的所有kmerput到map中
	 * for (int j = 0; j <= length - kh.kmer_length; j++) { String kmer =
	 * node_set.get(node_index).getSequence().substring(j, j + kh.kmer_length);
	 * long intval = baseOptions.kmerToIntval(kmer); int cov =
	 * kh.kmer_map.get(intval); node_kmer_count.add(cov); } Vector<Integer>
	 * last=new Vector<Integer>(); // 将 顶点覆盖的kmer的read——count从小到大排序 int node_cov
	 * = compute_node_bv(0, node_kmer_count,last); int kmer_cov=0; for(int
	 * i=0;i<node_set.get(node_index).getChildren().size();i++){ int child=(int)
	 * node_set.get(node_index).getChildren().get(i); String
	 * kmer=node_set.get(node_index).getSequence().substring(length-kh.
	 * kmer_length+1)+node_set.get(child).getSequence().substring(0,1);
	 * kmer_cov+=kh.kmer_map.get(baseOptions.kmerToIntval(kmer)); }
	 * if(last.lastElement()>kmer_cov){ node_cov=node_cov-kmer_cov; } else{
	 * node_cov=node_cov-last.lastElement(); } System.out.println(node_cov);
	 * 
	 * }
	 * 
	 * 
	 * public static int compute_edge_bv(Vector<Integer> edge_kmer_count,int
	 * edge_cov){ if (edge_kmer_count.size() == 1) { edge_cov +=
	 * edge_kmer_count.get(0); return edge_cov; } Vector<Integer>
	 * node_kmer_count_asc = new Vector<Integer>(); for (int i = 0; i <
	 * edge_kmer_count.size(); i++) {
	 * node_kmer_count_asc.add(edge_kmer_count.get(i)); } //
	 * node_kmer_count_asc=node_kmer_count;
	 * Collections.sort(node_kmer_count_asc); if(node_kmer_count_asc.size()==0){
	 * // System.out.println("!!!!!!!!!!!!!!!!!!!!!!"); return 0; } int temp =
	 * node_kmer_count_asc.get(0);
	 * 
	 * edge_cov += temp; // 重置node_kmer_count for (int i = 0; i <
	 * edge_kmer_count.size(); i++) { edge_kmer_count.set(i,
	 * edge_kmer_count.get(i) - temp); } while (true) { //
	 * System.out.println(node_kmer_count.size()); if (edge_kmer_count.size() ==
	 * 0) { return edge_cov; } boolean flag = false; boolean flag1 = false;
	 * while (edge_kmer_count.size() != 0 && edge_kmer_count.get(0) == 0) {
	 * edge_kmer_count.remove(0); flag1 = true; } if (flag1) { continue; } for
	 * (int i = 0; i < edge_kmer_count.size(); i++) { if (edge_kmer_count.get(i)
	 * == 0) { Vector<Integer> node_kmer_c = new Vector<Integer>(); for (int j =
	 * 0; j < i; j++) { node_kmer_c.add(edge_kmer_count.get(j)); }
	 * edge_cov=compute_edge_bv(node_kmer_c,edge_cov); return edge_cov;
	 * 
	 * } } if (flag == false) { edge_cov = compute_edge_bv(edge_kmer_count,
	 * edge_cov); return edge_cov; } } }
	 * 
	 * 
	 * public static int compute_node_bv(int node_cov, Vector<Integer>
	 * node_kmer_count,Vector<Integer> last) { if (node_kmer_count.size() == 1)
	 * { node_cov += node_kmer_count.get(0); last.add(node_kmer_count.get(0));
	 * return node_cov; } Vector<Integer> node_kmer_count_asc = new
	 * Vector<Integer>(); for (int i = 0; i < node_kmer_count.size(); i++) {
	 * node_kmer_count_asc.add(node_kmer_count.get(i)); } //
	 * node_kmer_count_asc=node_kmer_count;
	 * Collections.sort(node_kmer_count_asc); int temp =
	 * node_kmer_count_asc.get(0); node_cov += temp; last.add(temp); //
	 * 重置node_kmer_count for (int i = 0; i < node_kmer_count.size(); i++) {
	 * node_kmer_count.set(i, node_kmer_count.get(i) - temp); } while (true) {
	 * // System.out.println(node_kmer_count.size()); if (node_kmer_count.size()
	 * == 0) { return node_cov; } boolean flag = false; boolean flag1 = false;
	 * while (node_kmer_count.size() != 0 && node_kmer_count.get(0) == 0) {
	 * node_kmer_count.remove(0); flag1 = true; } if (flag1) { continue; } for
	 * (int i = 0; i < node_kmer_count.size(); i++) { if (node_kmer_count.get(i)
	 * == 0) { Vector<Integer> node_kmer_c = new Vector<Integer>(); for (int j =
	 * 0; j < i; j++) { node_kmer_c.add(node_kmer_count.get(j)); } node_cov =
	 * compute_node_bv(node_cov, node_kmer_c,last); for (int j = 0; j <= i; j++)
	 * { node_kmer_count.remove(0); } flag = true; break; } } if (flag == false)
	 * { node_cov = compute_node_bv(node_cov, node_kmer_count,last); return
	 * node_cov; } } }
	 * 
	 * public List sort_kmer_asc(Map<Long, Integer> node_kmer) {
	 * 
	 * List list = new ArrayList<Map.Entry<Long,
	 * Integer>>(node_kmer.entrySet());
	 * 
	 * // 通过比较器实现比较排序 Collections.sort(list, new Comparator<Map.Entry<Long,
	 * Integer>>() { public int compare(Map.Entry<Long, Integer> o1,
	 * Map.Entry<Long, Integer> o2) { return
	 * o1.getValue().compareTo(o2.getValue()); // 倒序 } }); return list; }
	 */
	public void compute_edge_cov(int[][] edges, kmerHash kh) {
		for (int i = 0; i < node_set.size(); i++) {
			int child;
			if (node_set.get(i).getChildren().size() == 0) {
				continue;
			} else {
				Vector<Integer> edge_kmer = new Vector<Integer>();
				for (int j = 0; j < node_set.get(i).getChildren().size(); j++) {
					child = (int) node_set.get(i).getChildren().get(j);
					int edge_cov = 0;
					if (node_set.get(i).getSequence().length() >= kh.kmer_length
							&& node_set.get(child).getSequence().length() >= kh.kmer_length) {
						String edge_str = node_set.get(i).getSequence()
								.substring(node_set.get(i).getSequence().length() - kh.kmer_length + 1)
								+ node_set.get(child).getSequence().substring(0, kh.kmer_length - 1);
						for (int k = 0; k <= edge_str.length() - kh.kmer_length; k++) {
							String kmer = edge_str.substring(k, k + kh.kmer_length);
							int cov1;
							if(kh.kmer_map.containsKey(baseOptions.kmerToIntval(kmer))){
								cov1 = kh.kmer_map.get(baseOptions.kmerToIntval(kmer));
							}
							else{
								cov1=0;
							}
							edges[i][child] += cov1;
						}
					} else if (node_set.get(i).getSequence().length() < kh.kmer_length
							&& node_set.get(child).getSequence().length() >= kh.kmer_length) {
						String edge_str = node_set.get(i).getSequence()
								+ node_set.get(child).getSequence().substring(0, kh.kmer_length - 1);
						for (int k = 0; k <= edge_str.length() - kh.kmer_length; k++) {
							String kmer = edge_str.substring(k, k + kh.kmer_length);
							int cov1;
							if(kh.kmer_map.containsKey(baseOptions.kmerToIntval(kmer))){
								cov1 = kh.kmer_map.get(baseOptions.kmerToIntval(kmer));
							}
							else{
								cov1=0;
							}
							edges[i][child] = cov1;
						}
					} else if (node_set.get(i).getSequence().length() >= kh.kmer_length
							&& node_set.get(child).getSequence().length() < kh.kmer_length) {
						String edge_str = node_set.get(i).getSequence()
								.substring(node_set.get(i).getSequence().length() - kh.kmer_length + 1)
								+ node_set.get(child).getSequence();
						for (int k = 0; k <= edge_str.length() - kh.kmer_length; k++) {
							String kmer = edge_str.substring(k, k + kh.kmer_length);
							int cov1;
							if(kh.kmer_map.containsKey(baseOptions.kmerToIntval(kmer))){
								cov1 = kh.kmer_map.get(baseOptions.kmerToIntval(kmer));
							}
							else{
								cov1=0;
							}
							edges[i][child] = cov1;
						}

					} else {
						String edge_str = node_set.get(i).getSequence() + node_set.get(child).getSequence();
						for (int k = 0; k <= edge_str.length() - kh.kmer_length; k++) {
							String kmer = edge_str.substring(k, k + kh.kmer_length);
							int cov1;
							if(kh.kmer_map.containsKey(baseOptions.kmerToIntval(kmer))){
								cov1 = kh.kmer_map.get(baseOptions.kmerToIntval(kmer));
							}
							else{
								cov1=0;
							}
							edges[i][child] = cov1;
						}

					}

				}
			}
		}
	}

	
	

	public void compute_fragment_count(Vector<Integer> path, kmerHash kh,float[] bv,int[][] f) {
		int length = 0;// length代表了路径的总长度
		Vector<Integer> terminal = new Vector<Integer>();// 记录个顶点之间的端点
		String path_str = "";
		//System.out.println(node_set);
		for (int i = 0; i < path.size(); i++) {
			length += node_set.get(path.get(i)).getSequence().length();
			path_str += node_set.get(path.get(i)).getSequence();
			terminal.add(length - 1);
		}

		Vector<Integer> kmer_counts = new Vector<Integer>();
		for (int i = 0; i <= length - kh.kmer_length; i++) {
			String kmer = path_str.substring(i, i + kh.kmer_length);
			if(kh.kmer_map.containsKey(baseOptions.kmerToIntval(kmer))){
			kmer_counts.add(kh.kmer_map.get(baseOptions.kmerToIntval(kmer)));}
			else{
				kmer_counts.add(0);
			}
		}
		int start = 0;
		int end = length - kh.kmer_length;
		
		System.out.println("fragment_count:"+node_set.size());
		float[][] fragment_count = new float[node_set.size()][node_set.size()];

		compute_fragment(kmer_counts,start,end,path,terminal,fragment_count,kh);
		System.out.println("fragment_count!");
		for (int i = 0; i < node_set.size(); i++) {
			for (int j = 0; j < node_set.size(); j++) {
				System.out.print(fragment_count[i][j] + "  ");
			}
			System.out.println();
		}
		System.out.println("fragment_count!结束");
		for(int i=0;i<node_set.size();i++){
			//最终的顶点bv
			//减去其公共边
			for(int j=0;j<node_set.size();j++){
				if(fragment_count[i][j]!=0){
					int start_p=0;
					int end_p=0;
					for(int k=0;k<path.size();k++){
						if(path.get(k)==i){
							start_p=k;
							break;
						}
					}
					for(int k=0;k<path.size();k++){
						if(path.get(k)==j){
							end_p=k;
							break;
						}
					}
					for(int k=start_p;k<=end_p;k++){
						bv[path.get(k)]-=fragment_count[i][j]; //fxx
					}
				}
			}
		
		}
		for(int i=0;i<node_set.size();i++){
			fragment_count[i][i]=bv[i];
		}
		for(int i=0;i<node_set.size();i++){
			int fenzi=0;
			int fenmu=0;
			for(int j=0;j<node_set.size();j++){
				if(fragment_count[i][j]!=0){
					fenmu+=fragment_count[i][j];
				}
				if(fragment_count[j][i]!=0){
					fenzi+=fragment_count[j][i];
				}
			}
			bv[i]=(float)fenzi/fenmu;
		}
		
	}

	int start_biao=0;int end_biao=0;
	public void compute_fragment(Vector<Integer> kmer_counts, int start, int end, Vector<Integer> path,
			Vector<Integer> terminal, float[][] fragment_count,kmerHash kh) {
		
		int start_index = 0;
		int end_index = 0;
		for (int i = 0; i < terminal.size(); i++) {
			if (start < terminal.get(i)) {
				// 确定从哪个顶点开始
				start_index = i;
				break;
			}
		}		
		if(start+kh.kmer_length-1>terminal.get(start_index)){
			return;
		}
		if(kmer_counts.size()==1){
			if(start_biao+kh.kmer_length-1<terminal.get(start_index)){
				return;
			}
			else{
				for (int j = terminal.size() - 1; j >= 0; j--) {
					if (start_biao+kh.kmer_length-1 > terminal.get(j)) {
						//确定从哪个顶点结束
						end_index = j+1;
						break;
					}
				}
				fragment_count[path.get(start_index)][path.get(end_index)]+=kmer_counts.get(0);
				return;
			}
		}

		for (int j = terminal.size() - 1; j >= 0; j--) {
			if (end+kh.kmer_length-1 > terminal.get(j)) {
				//确定从哪个顶点结束
				end_index = j+1;
				break;
			}
		}
		if(end_index==start_index&&end+kh.kmer_length-1<=terminal.get(end_index)){
			return;
		}
		Vector<Integer> kmer_counts_asc = new Vector<Integer>();
		for (int i = 0; i < kmer_counts.size(); i++) {
			kmer_counts_asc.add(kmer_counts.get(i));
		}
	//	System.out.println("kmer_counts_asc:"+kmer_counts_asc);
		if(kmer_counts_asc.size()==0||kmer_counts_asc==null){
			return;
		}
		
		Collections.sort(kmer_counts_asc);
		
		
		int temp = kmer_counts_asc.get(0);// 取最小的cov
        kmer_counts.add(temp);//保证每次末尾都有0
        //System.out.println("start_index:"+start_index+"   end_index:"+end_index+"   path:"+path.get(start_index)+"->"+path.get(end_index));
		//System.out.println(path.size());
        fragment_count[path.get(start_index)][path.get(end_index)] += temp;
		// 重置node_kmer_count,并记录kmer——count为0的标号
		Vector<Integer> count_is_zero=new Vector<Integer>();
		for (int i = 0; i < kmer_counts.size(); i++) {
			kmer_counts.set(i, kmer_counts.get(i) - temp);
			if(kmer_counts.get(i)==0){
				count_is_zero.add(i);
			}
			
		}
 		int cer=start;
		//int flag=0;
		for(int i=0;i<count_is_zero.size();i++){
			if(i==count_is_zero.size()-1){
				//如果是最后一个0，则为人为添加，不做处理
				end_biao=cer+count_is_zero.get(i)-1;
				Vector<Integer> kmer_counts_c=new Vector<Integer>();
				for(int j=start_biao-cer;j<=end_biao-cer;j++){
					kmer_counts_c.add(kmer_counts.get(j));
				}
				compute_fragment(kmer_counts_c, start_biao, end_biao, path, terminal, fragment_count,kh);
				
				break;
			}
			if(count_is_zero.get(i)==0){
				//如果第一个点是0
				start_biao++;
				continue;
			}
			int index=count_is_zero.get(i);
			if(index+cer==start_biao){
				start_biao++;
				continue;
			}
			int index_index=0;
			//得到index所在的顶点标号
			//flag=index;
			int start_pos=0;
			for (int j = 0; j < terminal.size(); j++) {
				if (start_biao < terminal.get(j)) {
					// 确定从哪个顶点开始
					start_pos =j;
					break;
				}
			}
			boolean flag1=false;
			for (int j = terminal.size() - 1; j >= 0; j--) {
				if (index+cer-1+kh.kmer_length-1 > terminal.get(j)) {
					//确定从哪个顶点结束
					index_index = j+1;
					flag1=true;
					break;
				}
			}
			if(flag1==false){
				index_index=0;
			}
			//如果该片段在同一个顶点且不能与其他顶点相连时，跳过处理
			if(start_biao+index-1+kh.kmer_length-1<terminal.get(index_index)&&index_index==start_pos){
				start_biao=start_biao+index+1;
				continue;
			}
			if(start_biao+index-1+kh.kmer_length-1==terminal.get(index_index)&&index_index==start_pos){
				start_biao=start_biao+index+1;
				continue;
			}
			end_biao=cer+index-1;
			Vector<Integer> kmer_counts_c=new Vector<Integer>();
			for(int j=start_biao-cer;j<=end_biao-cer;j++){
				kmer_counts_c.add(kmer_counts.get(j));
			}
			compute_fragment(kmer_counts_c, start_biao, end_biao, path, terminal, fragment_count,kh);
			start_biao=cer+index+1;
		}
		
		
	}

	public  float compute_node_bv(float node_cov, Vector<Integer> node_kmer_count) {
		if (node_kmer_count.size() == 1) {
			node_cov += node_kmer_count.get(0);
		//	last.add(node_kmer_count.get(0));
			return node_cov;
		}
		Vector<Integer> node_kmer_count_asc = new Vector<Integer>();
		for (int i = 0; i < node_kmer_count.size(); i++) {
			node_kmer_count_asc.add(node_kmer_count.get(i));
		}
		Collections.sort(node_kmer_count_asc);
		if(node_kmer_count_asc.size()==0){
			return node_cov;
		}
		int temp = node_kmer_count_asc.get(0);
		node_cov += temp;
	//	last.add(temp);
		// 重置node_kmer_count
		for (int i = 0; i < node_kmer_count.size(); i++) {
			node_kmer_count.set(i, node_kmer_count.get(i) - temp);
		}
		while (true) {
			if (node_kmer_count.size() == 0) {
				return node_cov;
			}
			boolean flag = false;
			boolean flag1 = false;
			while (node_kmer_count.size() != 0 && node_kmer_count.get(0) == 0) {
				node_kmer_count.remove(0);
				flag1 = true;
			}
			if (flag1) {
				continue;
			}
			for (int i = 0; i < node_kmer_count.size(); i++) {
				if (node_kmer_count.get(i) == 0) {
					Vector<Integer> node_kmer_c = new Vector<Integer>();
					for (int j = 0; j < i; j++) {
						node_kmer_c.add(node_kmer_count.get(j));
					}
					node_cov = compute_node_bv(node_cov, node_kmer_c);
					for (int j = 0; j <= i; j++) {
						node_kmer_count.remove(0);
					}
					flag = true;
					break;
				}
			}
			if (flag == false) {
				node_cov = compute_node_bv(node_cov, node_kmer_count);
				return node_cov;
			}
		}
	}

	public List sort_kmer_asc(Map<Long, Integer> node_kmer) {

		List list = new ArrayList<Map.Entry<Long, Integer>>(node_kmer.entrySet());

		// 通过比较器实现比较排序
		Collections.sort(list, new Comparator<Map.Entry<Long, Integer>>() {
			public int compare(Map.Entry<Long, Integer> o1, Map.Entry<Long, Integer> o2) {
				return o1.getValue().compareTo(o2.getValue()); // 倒序
			}
		});
		return list;
	}

	Stack<Integer> stack = new Stack<Integer>();
	Vector<Integer> path = new Vector<Integer>();
	public void dfs(int[][] edges, int src, int des,kmerHash kh,float[] bv,int[][] f) {
		// TODO Auto-generated method stub
		//bv是每个顶点最开始的fragment大小，现在要减去公共边的数量
		stack.add(src);
		path.add(src);
		if (src == des) {
			System.out.println(path);
			compute_fragment_count(path, kh,bv,f);
			return;
		}
		for (int i = 0; i < edges.length; i++) {
			if (edges[src][i] != 0) {
				dfs(edges, i, des,kh,bv,f);
				path.remove(path.size() - 1);
			}
		}
	}
	
	public void compute_node_cov(int[] nodes, kmerHash kh,float[] bv) {
		for (int i = 0; i < node_set.size(); i++) {
			Vector<Integer> node_kmer_count=new Vector<Integer>();
			int length = node_set.get(i).getSequence().length();
			int node_cov = 0;
			for (int j = 0; j <= length - kh.kmer_length; j++) {
				String kmer = node_set.get(i).getSequence().substring(j, j + kh.kmer_length);
				int cov ;
				if(kh.kmer_map.containsKey(baseOptions.kmerToIntval(kmer))){
					
				cov = kh.kmer_map.get(baseOptions.kmerToIntval(kmer));
				node_cov += cov;
				node_kmer_count.add(cov);}
				else{
					node_kmer_count.add(0);
				}
				
			}
			nodes[i] = node_cov;
			bv[i]=compute_node_bv(0, node_kmer_count);
		}
	}


}
