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

public class test {
	int start_biao=0;int end_biao=0;
	public void compute_fragment(Vector<Integer> kmer_counts, int start, int end, Vector<Integer> path,
			Vector<Integer> terminal, int[][] fragment_count) {
		
		int start_index = 0;
		int end_index = 0;
		for (int i = 0; i < terminal.size(); i++) {
			if (start < terminal.get(i)) {
				// 确定从哪个顶点开始
				start_index = path.get(i);
				break;
			}
		}
		
		if(kmer_counts.size()==1){
			if(start_biao+5<terminal.get(start_index)){
				return;
			}
			else{
				for (int j = terminal.size() - 1; j >= 0; j--) {
					if (start_biao+5 > terminal.get(j)) {
						//确定从哪个顶点结束
						end_index = path.get(j + 1);
						break;
					}
				}
				fragment_count[start_index][end_index]+=kmer_counts.get(0);
				return;
			}
		}

		for (int j = terminal.size() - 1; j >= 0; j--) {
			if (end > terminal.get(j)) {
				//确定从哪个顶点结束
				end_index = path.get(j + 1);
				break;
			}
		}
		Vector<Integer> kmer_counts_asc = new Vector<Integer>();
		for (int i = 0; i < kmer_counts.size(); i++) {
			kmer_counts_asc.add(kmer_counts.get(i));
		}
		Collections.sort(kmer_counts_asc);
		int temp = kmer_counts_asc.get(0);// 取最小的cov
        kmer_counts.add(temp);//保证每次末尾都有0
		fragment_count[start_index][end_index] += temp;
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
				return;
			}
			if(count_is_zero.get(i)==0){
				//如果第一个点是0
				start_biao++;
				continue;
			}
			int index=count_is_zero.get(i);
			if(index==start_biao){
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
					start_pos = path.get(j);
					break;
				}
			}
			boolean flag1=false;
			for (int j = terminal.size() - 1; j >= 0; j--) {
				if (index+cer-1 > terminal.get(j)) {
					//确定从哪个顶点结束
					index_index = path.get(j + 1);
					flag1=true;
					break;
				}
			}
			if(flag1==false){
				index_index=0;
			}
			//如果该片段在同一个顶点且不能与其他顶点相连时，跳过处理
			if(start_biao+index-1+5<terminal.get(index_index)&&index_index==start_pos){
				start_biao=start_biao+index+1;
				continue;
			}
			if(start_biao+index-1+5==terminal.get(index_index)&&index_index==start_pos){
				start_biao=start_biao+index+1;
				continue;
			}
			end_biao+=index-1;
			Vector<Integer> kmer_counts_c=new Vector<Integer>();
			for(int j=start_biao;j<=end_biao;j++){
				kmer_counts_c.add(kmer_counts.get(j));
			}
			compute_fragment(kmer_counts_c, start_biao, end_biao, path, terminal, fragment_count);
		}
		
		
	}
/*
	public void compute_fragment(Vector<Integer> kmer_counts, int start, int end, Vector<Integer> path,
			Vector<Integer> terminal, int[][] fragment_count) {
		int start_index = 0;
		int end_index = 0;
		for (int i = 0; i < terminal.size(); i++) {
			if (start < terminal.get(i)) {
				// 确定从哪个顶点开始
				start_index = path.get(i);
				break;
			}
		}
		for (int j = terminal.size() - 1; j >= 0; j--) {
			if (end > terminal.get(j)) {
				end_index = path.get(j + 1);
				break;
			}
		}
		if (kmer_counts.size() == 1) {
			fragment_count[start_index][end_index] += kmer_counts.get(0);
			return;
		}
		Vector<Integer> kmer_counts_asc = new Vector<Integer>();
		for (int i = 0; i < kmer_counts.size(); i++) {
			kmer_counts_asc.add(kmer_counts.get(i));
		}
		Collections.sort(kmer_counts_asc);
		int temp = kmer_counts_asc.get(0);// 取最小的cov

		fragment_count[start_index][end_index] += temp;
		// 重置node_kmer_count
		for (int i = 0; i < kmer_counts.size(); i++) {
			kmer_counts.set(i, kmer_counts.get(i) - temp);
		}
		while (true) {
			while (kmer_counts.get(start) == 0) {
				kmer_counts.set(start, Integer.MAX_VALUE);
				start++;
				
			}
			int index = 0;
			for (int i = start; i < kmer_counts.size(); i++) {
				if (kmer_counts.get(i) == 0) {
					//kmer_counts.set(i, Integer.MAX_VALUE);
					index = i - 1;
					break;
				}
			}
			if (index + 5 < terminal.get(start_index)) {
				index++;
				start = index;
				continue;
			} else {
				Vector<Integer> kmer_counts_c = new Vector<Integer>();
				for (int k = start; k <= index; k++) {
					kmer_counts_c.add(kmer_counts.get(k));
				}
				compute_fragment(kmer_counts_c, start, index, path, terminal, fragment_count);
			}
		}
	}*/
	public static void main(String args[]) {
		Vector<Integer> kmer_counts = new Vector<Integer>();
		kmer_counts.add(5);
		kmer_counts.add(4);
		kmer_counts.add(6);
		kmer_counts.add(10);
		kmer_counts.add(3);
		kmer_counts.add(3);
		kmer_counts.add(10);
		kmer_counts.add(5);
		kmer_counts.add(6);
		kmer_counts.add(7);
		kmer_counts.add(9);
		kmer_counts.add(6);
		kmer_counts.add(6);
		kmer_counts.add(7);
		kmer_counts.add(9);
		kmer_counts.add(3);
		Vector<Integer> path = new Vector<Integer>();
		path.add(0);
		path.add(1);
		Vector<Integer> terminal = new Vector<Integer>();
		terminal.add(11);
		terminal.add(20);
		int[][] fragment_count = new int[path.size()][path.size()];
		test ts = new test();
		ts.compute_fragment(kmer_counts, 0, 15, path, terminal, fragment_count);
		for (int i = 0; i < path.size(); i++) {
			for (int j = 0; j < path.size(); j++) {
				System.out.print(fragment_count[i][j] + "  ");
			}
			System.out.println();
		}
	}
}
