package test;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.Vector;

import org.junit.Before;
import org.junit.Test;

public class testkmerHash {
	/**
	 * 对kmerHash类进行测试
	 */
	kmerHash testKh = new kmerHash();
	Map<Long, Vector<Integer>> testKmerHash = new HashMap<Long, Vector<Integer>>();
	Vector<String> testRead_vector = new Vector<String>();

	@Before
	public void setUp() throws Exception {
		 System.out.println("请输入kmer长度（kmer长度需在15-31之间，默认25）.........");
			Scanner scan = new Scanner(System.in);
			int read = scan.nextInt();
			while (read > 31 || read < 15) {
				System.out.println("您输入的数据不合法，请重新输入！");
				read = scan.nextInt();
			}
			System.out.println("您输入的kmer长度是：" + read);
			testKh.kmer_length = read;
			testRead_vector.add("GGGCCCCTGCCTGGGGGCTTGTCACCTCCCCCACCTTCTTCCTGAGTCATTCCTGCAGCC");
			testRead_vector.add("GAGCCATGCCTAGAGTGGGATGGGCCATTGTTCATCTTCTGGCCCCTGTTGTCTGCATGT");
			testRead_vector.add("CCTTTTCTTTGACCTCTTCTTTCTGTTCATGTGTATTTGCTGTCTCTTAGCCCAGACTTC");
			testRead_vector.add("AGGAGCCATGCCTAGAGTGGGATGGGCCATTGTTCATCTTCTGGCCCCTGTTGTCTGCAT");
			testRead_vector.add("GGGATGGGCCATTGTTCATCTTCTGGCCCCTGTTGTCTGCATGTAACTTAATACCACAAC");
			testKh.readsToKmer(testKmerHash, testRead_vector);
	}

//	@Test
	public void testReadsToKmer() {
//		System.out.println("请输入kmer长度（kmer长度需在15-31之间，默认25）.........");
//		Scanner scan = new Scanner(System.in);
//		int read = scan.nextInt();
//		while (read > 31 || read < 15) {
//			System.out.println("您输入的数据不合法，请重新输入！");
//			read = scan.nextInt();
//		}
//		System.out.println("您输入的kmer长度是：" + read);
//		testKh.kmer_length = read;
//		// testRead_vector.add("CAGCATCAGGTCTCCAGAGCTGCAGAAGACGACGGCCGACTTGGATCACACTCTTGTGAG");
//		// testRead_vector.add("TCCGTGAGATCTTCCCAGGGCAGCTCCCCTCTGTGGAATCCAATCTGTCTTCCATCCTGC");
//		testRead_vector.add("GGGCCCCTGCCTGGGGGCTTGTCACCTCCCCCACCTTCTTCCTGAGTCATTCCTGCAGCC");
//		testRead_vector.add("GAGCCATGCCTAGAGTGGGATGGGCCATTGTTCATCTTCTGGCCCCTGTTGTCTGCATGT");
//		testRead_vector.add("CCTTTTCTTTGACCTCTTCTTTCTGTTCATGTGTATTTGCTGTCTCTTAGCCCAGACTTC");
//		testRead_vector.add("AGGAGCCATGCCTAGAGTGGGATGGGCCATTGTTCATCTTCTGGCCCCTGTTGTCTGCAT");
//		testRead_vector.add("GGGATGGGCCATTGTTCATCTTCTGGCCCCTGTTGTCTGCATGTAACTTAATACCACAAC");
//		// testRead_vector.add("CTGGCCATCCGTGAGATCTTCCCAGGGCAGCTCCCCTCTGTGGAATCCAATCTGTCTTCC");
//		// testRead_vector.add("AGCATGACTATTTTTAGAGACCCCGTGTCTGTCACTGAAACCTTTTTTGTGGGAGACTAT");
//		// testRead_vector.add("CCGTCAGCCTTTTCTTTGACCTCTTCTTTCTGTTCATGTGTATTTGCTGTCTCTTAGCCC");
//		testKh.readsToKmer(testKmerHash, testRead_vector);
		// for (long key : testKmerHash.keySet()) {
		// System.out.println(key+" "+testKmerHash.get(key));
		// }

		Map<Long, Vector<Integer>> TempMap = new HashMap<Long, Vector<Integer>>();
		Vector<Integer> read_set = new Vector<Integer>();
		// 1,3
		read_set.add(1);
		read_set.add(3);
		TempMap.put(163714636091045l, read_set);
		TempMap.put(1030149314137591l, read_set);
		TempMap.put(1010601144521979l, read_set);
		TempMap.put(502178559377140l, read_set);
		TempMap.put(654858544364180l, read_set);
		TempMap.put(614230408543543l, read_set);
		TempMap.put(820487281955709l, read_set);
		TempMap.put(603878612444073l, read_set);
		TempMap.put(153557602135885l, read_set);
		TempMap.put(882814330665939l, read_set);
		TempMap.put(367634363771475l, read_set);
		TempMap.put(344637548243279l, read_set);
		TempMap.put(664704857560047l, read_set);
		TempMap.put(252650286130494l, read_set);
		TempMap.put(205121820488927l, read_set);
		TempMap.put(407019616554941l, read_set);
		// TempMap.put(key, read_set);
		// TempMap.put(key, read_set);
		// TempMap.put(key, read_set);
		// TempMap.put(key, read_set);
		// TempMap.put(key, read_set);

		// read_set.clear();
		// 1,3,4
		Vector<Integer> read_set1 = new Vector<Integer>();
		read_set1.add(1);
		read_set1.add(3);
		read_set1.add(4);
		TempMap.put(1030584872197781l, read_set1);
		TempMap.put(1076962669346526l, read_set1);
		TempMap.put(930150956858233l, read_set1);
		TempMap.put(726759259367767l, read_set1);
		TempMap.put(742897536022494l, read_set1);
		TempMap.put(257646218049445l, read_set1);
		TempMap.put(744639768263253l, read_set1);
		TempMap.put(627361507933673l, read_set1);
		TempMap.put(276885769115582l, read_set1);
		TempMap.put(350696418989551l, read_set1);
		TempMap.put(1052472585321453l, read_set1);
		TempMap.put(832190620757943l, read_set1);
		TempMap.put(342904106905060l, read_set1);
		TempMap.put(655237223785822l, read_set1);
		TempMap.put(245716520777619l, read_set1);
		TempMap.put(1107543076462331l, read_set1);
		TempMap.put(719790330404730l, read_set1);
		TempMap.put(369149081458043l, read_set1);
		// TempMap.put(key, read_set);

		// read_set.clear();
		// 0
		Vector<Integer> read_set2 = new Vector<Integer>();
		read_set2.add(0);
		TempMap.put(702388035149693l, read_set2);
		TempMap.put(539481934878021l, read_set2);
		TempMap.put(938183929334607l, read_set2);
		TempMap.put(557752326913525l, read_set2);
		TempMap.put(1032027832669463l, read_set2);
		TempMap.put(374244080399605l, read_set2);
		TempMap.put(553980103974473l, read_set2);
		TempMap.put(416345460430161l, read_set2);
		TempMap.put(419970002704274l, read_set2);
		TempMap.put(102620054117933l, read_set2);
		TempMap.put(410480216471732l, read_set2);
		TempMap.put(727298934963013l, read_set2);
		TempMap.put(307723101879673l, read_set2);
		TempMap.put(1105109400811479l, read_set2);
		TempMap.put(516020959044307l, read_set2);
		TempMap.put(375035996810557l, read_set2);
		TempMap.put(657395926166807l, read_set2);
		TempMap.put(750411610149983l, read_set2);
		TempMap.put(415042721498581l, read_set2);
		TempMap.put(385235657085301l, read_set2);
		TempMap.put(377783890981981l, read_set2);
		TempMap.put(534270979151701l, read_set2);
		TempMap.put(749846626914685l, read_set2);
		TempMap.put(793251810344312l, read_set2);
		TempMap.put(747586693973495l, read_set2);
		TempMap.put(1090020509055269l, read_set2);
		TempMap.put(1011184009764181l, read_set2);
		TempMap.put(921207427692002l, read_set2);
		TempMap.put(371076414755799l, read_set2);
		TempMap.put(1042737882718046l, read_set2);
		TempMap.put(358405752180574l, read_set2);
		TempMap.put(744774687162065l, read_set2);
		TempMap.put(307129990240139l, read_set2);
		TempMap.put(667036318528852l, read_set2);
		TempMap.put(104992500676068l, read_set2);
		TempMap.put(738546962208735l, read_set2);
		// TempMap.put(key, read_set);
		// TempMap.put(key, read_set);

		// read_set.clear();
		// 2
		Vector<Integer> read_set3 = new Vector<Integer>();
		read_set3.add(2);
		TempMap.put(103319594679219l, read_set3);
		TempMap.put(262576713202837l, read_set3);
		TempMap.put(562398148657019l, read_set3);
		TempMap.put(1116817052662253l, read_set3);
		TempMap.put(413278378716879l, read_set3);
		TempMap.put(914376491160711l, read_set3);
		TempMap.put(823527690717522l, read_set3);
		TempMap.put(1077226109853170l, read_set3);
		TempMap.put(1090573615157223l, read_set3);
		TempMap.put(1090584401927924l, read_set3);
		TempMap.put(544597239298935l, read_set3);
		TempMap.put(1119224976459901l, read_set3);
		TempMap.put(791544076211489l, read_set3);
		TempMap.put(1050306852811348l, read_set3);
		TempMap.put(1117071030613949l, read_set3);
		TempMap.put(1117071030613949l, read_set3);
		TempMap.put(279806244114975l, read_set3);
		TempMap.put(560851828207438l, read_set3);
		TempMap.put(527213608024895l, read_set3);
		TempMap.put(1117507405987131l, read_set3);
		TempMap.put(1089568490121143l, read_set3);
		TempMap.put(588779852091116l, read_set3);
		TempMap.put(347119155011365l, read_set3);
		TempMap.put(560679239876219l, read_set3);
		TempMap.put(832256480884604l, read_set3);
		TempMap.put(422074513874910l, read_set3);
		TempMap.put(1123692687785455l, read_set3);
		TempMap.put(1042310949184840l, read_set3);
		TempMap.put(931204718884809l, read_set3);
		TempMap.put(984594740101022l, read_set3);
		TempMap.put(982954525256958l, read_set3);
		TempMap.put(991619893154747l, read_set3);
		TempMap.put(1052489050353119l, read_set3);
		TempMap.put(984637887183827l, read_set3);
		TempMap.put(980574239956701l, read_set3);
		TempMap.put(554118380499961l, read_set3);
		TempMap.put(1092329903420654l, read_set3);
		// TempMap.put(key, read_set);
		// TempMap.put(key, read_set);

		// read_set.clear();
		// 4
		Vector<Integer> read_set4 = new Vector<Integer>();
		read_set4.add(4);
		System.out.println(read_set);
		TempMap.put(1054286943482053l, read_set4);
		TempMap.put(977719370388545l, read_set4);
		TempMap.put(807379796018448l, read_set4);
		TempMap.put(1105992399916113l, read_set4);
		TempMap.put(1046269879136580l, read_set4);
		TempMap.put(1089158540813548l, read_set4);
		TempMap.put(727309458140191l, read_set4);
		TempMap.put(839448053400340l, read_set4);
		TempMap.put(978934442726320l, read_set4);
		TempMap.put(538038050377409l, read_set4);
		TempMap.put(377952261816816l, read_set4);
		System.out.println("***" + TempMap.get(377952261816816l));
		TempMap.put(657438018875516l, read_set4);
		TempMap.put(385909140424643l, read_set4);
		TempMap.put(545046712581169l, read_set4);
		TempMap.put(417736654855948l, read_set4);
		TempMap.put(1026252294667015l, read_set4);

		// read_set.clear();
		// 1,4
		Vector<Integer> read_set5 = new Vector<Integer>();
		read_set5.add(1);
		read_set5.add(4);
		TempMap.put(553764611914043l, read_set5);
		TempMap.put(982866083110478l, read_set5);
		// TempMap.put(key, read_set);
		// TempMap.put(key, read_set);
		// TempMap.put(key, read_set);

		// read_set.clear();
		// 3
		Vector<Integer> read_set6 = new Vector<Integer>();
		read_set6.add(3);
		TempMap.put(178479901633082l, read_set6);
		TempMap.put(713919606532330l, read_set6);
		// TempMap.put(key, read_set);
		// TempMap.put(key, read_set);

		System.out.println(TempMap.size() + "   " + testKmerHash.size());
		if (TempMap.equals(testKmerHash)) {
			System.out.println("********");
		}
		assertEquals(TempMap, testKmerHash);
		for (long key : testKmerHash.keySet()) {
			System.out.print(key + "   " + testKmerHash.get(key) + "   ");
		}
		System.out.println();
		System.out.println(TempMap.get(377952261816816l));
		for (long key : TempMap.keySet()) {
			System.out.print(key + "   " + TempMap.get(key) + "   ");
		}
		System.out.println();
		System.out.println(TempMap.get(377952261816816l));
	}

//	@Test
	public void testGet_readset() {
//		System.out.println("请输入kmer长度（kmer长度需在15-31之间，默认25）.........");
//		Scanner scan = new Scanner(System.in);
//		int read = scan.nextInt();
//		while (read > 31 || read < 15) {
//			System.out.println("您输入的数据不合法，请重新输入！");
//			read = scan.nextInt();
//		}
//		System.out.println("您输入的kmer长度是：" + read);
//		testKh.kmer_length = read;
//		// testRead_vector.add("CAGCATCAGGTCTCCAGAGCTGCAGAAGACGACGGCCGACTTGGATCACACTCTTGTGAG");
//		// testRead_vector.add("TCCGTGAGATCTTCCCAGGGCAGCTCCCCTCTGTGGAATCCAATCTGTCTTCCATCCTGC");
//		testRead_vector.add("GGGCCCCTGCCTGGGGGCTTGTCACCTCCCCCACCTTCTTCCTGAGTCATTCCTGCAGCC");
//		testRead_vector.add("GAGCCATGCCTAGAGTGGGATGGGCCATTGTTCATCTTCTGGCCCCTGTTGTCTGCATGT");
//		testRead_vector.add("CCTTTTCTTTGACCTCTTCTTTCTGTTCATGTGTATTTGCTGTCTCTTAGCCCAGACTTC");
//		testRead_vector.add("AGGAGCCATGCCTAGAGTGGGATGGGCCATTGTTCATCTTCTGGCCCCTGTTGTCTGCAT");
//		testRead_vector.add("GGGATGGGCCATTGTTCATCTTCTGGCCCCTGTTGTCTGCATGTAACTTAATACCACAAC");
//		// testRead_vector.add("CTGGCCATCCGTGAGATCTTCCCAGGGCAGCTCCCCTCTGTGGAATCCAATCTGTCTTCC");
//		// testRead_vector.add("AGCATGACTATTTTTAGAGACCCCGTGTCTGTCACTGAAACCTTTTTTGTGGGAGACTAT");
//		// testRead_vector.add("CCGTCAGCCTTTTCTTTGACCTCTTCTTTCTGTTCATGTGTATTTGCTGTCTCTTAGCCC");
//		testKh.readsToKmer(testKmerHash, testRead_vector);
		Vector<Integer> read_temp = new Vector<Integer>();
		read_temp.add(1);
		read_temp.add(3);
		read_temp.add(4);
		System.out.println(testKmerHash);
		Vector<Integer> read_act = testKh.get_readset(1030584872197781l, testKmerHash);
		System.out.println(read_temp);
		System.out.println(read_act);
		assertEquals(read_temp, read_act);

	}

//	 @Test
	public void testGet_readset_count() {
//		 System.out.println("请输入kmer长度（kmer长度需在15-31之间，默认25）.........");
//			Scanner scan = new Scanner(System.in);
//			int read = scan.nextInt();
//			while (read > 31 || read < 15) {
//				System.out.println("您输入的数据不合法，请重新输入！");
//				read = scan.nextInt();
//			}
//			System.out.println("您输入的kmer长度是：" + read);
//			testKh.kmer_length = read;
//			testRead_vector.add("GGGCCCCTGCCTGGGGGCTTGTCACCTCCCCCACCTTCTTCCTGAGTCATTCCTGCAGCC");
//			testRead_vector.add("GAGCCATGCCTAGAGTGGGATGGGCCATTGTTCATCTTCTGGCCCCTGTTGTCTGCATGT");
//			testRead_vector.add("CCTTTTCTTTGACCTCTTCTTTCTGTTCATGTGTATTTGCTGTCTCTTAGCCCAGACTTC");
//			testRead_vector.add("AGGAGCCATGCCTAGAGTGGGATGGGCCATTGTTCATCTTCTGGCCCCTGTTGTCTGCAT");
//			testRead_vector.add("GGGATGGGCCATTGTTCATCTTCTGGCCCCTGTTGTCTGCATGTAACTTAATACCACAAC");
//			testKh.readsToKmer(testKmerHash, testRead_vector);
			assertEquals(3, testKh.get_readset_count(testKmerHash, 1030584872197781l));
		//fail("Not yet implemented");
	}

	// @Test
	public void testGet_forward_candidates() {
		 System.out.println("testGet_forward_candidates"+testKh.get_forward_candidates(377952261816816l,testKmerHash ));
		 List<Map.Entry<Long, Integer>> candidates =testKh.get_forward_candidates(377952261816816l,testKmerHash );
		long candidate=0l;
		 for (Map.Entry<Long, Integer> mapping : candidates) {
				candidate = mapping.getKey();
		 }
		 assertEquals(385909140424643l,candidate);
		 //fail("Not yet implemented");
	}

	 @Test
	public void testGet_reverse_candidates() {
		 System.out.println("testGet_reverse_candidates"+testKh.get_reverse_candidates(377952261816816l,testKmerHash ));
		 List<Map.Entry<Long, Integer>> candidates =testKh.get_reverse_candidates(377952261816816l,testKmerHash );
		 long candidate=0l;
		 for (Map.Entry<Long, Integer> mapping : candidates) {
				candidate = mapping.getKey();
				System.out.println(candidate);
		 }
		 assertEquals(657438018875516l,candidate);
		//fail("Not yet implemented");
	}

	 @Test
	public void testSort_kmer() {
		 List<Map.Entry<Long, Integer>> list=testKh.sort_kmer(testKmerHash);
		 System.out.println(testKh.sort_kmer(testKmerHash));
		 Vector<Long> sort_kmer=new Vector<Long>();
		 long kmer=0l;
		 for (Map.Entry<Long, Integer> mapping : list) {
				kmer = mapping.getKey();
				sort_kmer.add(kmer);
				//System.out.println(candidate);
		 }
		 Vector<Long> kmer_act=new Vector<Long>();
		 kmer_act.add(245716520777619l);
		 kmer_act.add(1107543076462331l);
		 kmer_act.add(369149081458043l);
		 kmer_act.add(719790330404730l);
		 kmer_act.add(832190620757943l);
		 kmer_act.add(342904106905060l);
		 kmer_act.add(655237223785822l);
		 kmer_act.add(1052472585321453l);
		 kmer_act.add(350696418989551l);
		 kmer_act.add(276885769115582l);
		 kmer_act.add(627361507933673l);
		 kmer_act.add(744639768263253l);
		 kmer_act.add(257646218049445l);
		 kmer_act.add(742897536022494l);
		 kmer_act.add(726759259367767l);
		 kmer_act.add(930150956858233l);
		 kmer_act.add(1030584872197781l);
		 kmer_act.add(1076962669346526l);
		 kmer_act.add(982866083110478l);
		 kmer_act.add(407019616554941l);
		 kmer_act.add(553764611914043l);
		 kmer_act.add(205121820488927l);
		 kmer_act.add(252650286130494l);
		 kmer_act.add(664704857560047l);
		 kmer_act.add(344637548243279l);
		 kmer_act.add(367634363771475l);
		 kmer_act.add(882814330665939l);
		 kmer_act.add(153557602135885l);
		 kmer_act.add(603878612444073l);
		 kmer_act.add(820487281955709l);
		 kmer_act.add(614230408543543l);
		 kmer_act.add(1010601144521979l);
		 kmer_act.add(502178559377140l);
		 kmer_act.add(654858544364180l);
		 kmer_act.add(1030149314137591l);
		 kmer_act.add(163714636091045l);
		 kmer_act.add(377952261816816l);
		 kmer_act.add(657438018875516l);
		 kmer_act.add(385909140424643l);
		 kmer_act.add(417736654855948l);
		 kmer_act.add(1042737882718046l);
		 kmer_act.add(991619893154747l);
		 kmer_act.add(1052489050353119l);
		 kmer_act.add(358405752180574l);
		 kmer_act.add(1026252294667015l);
		 kmer_act.add(744774687162065l);
		 kmer_act.add(307129990240139l);
		 kmer_act.add(984637887183827l);
		 kmer_act.add(667036318528852l);
		 kmer_act.add(104992500676068l);
		 kmer_act.add(713919606532330l);
		 kmer_act.add(980574239956701l);
		 kmer_act.add(554118380499961l);
		 kmer_act.add(545046712581169l);
		 kmer_act.add(738546962208735l);
		 kmer_act.add(1092329903420654l);
		 kmer_act.add(1090020509055269l);
		 kmer_act.add(1011184009764181l);
		 kmer_act.add(921207427692002l);
		 kmer_act.add(371076414755799l);
		 kmer_act.add(538038050377409l);
		 kmer_act.add(178479901633082l);
		 kmer_act.add(747586693973495l);
		 kmer_act.add(982954525256958l);
		 kmer_act.add(793251810344312l);
		 kmer_act.add(749846626914685l);
		 kmer_act.add(978934442726320l);
		 kmer_act.add(534270979151701l);
		 kmer_act.add(377783890981981l);
		 kmer_act.add(385235657085301l);
		 kmer_act.add(415042721498581l);
		 kmer_act.add(750411610149983l);
		 kmer_act.add(657395926166807l);
		 kmer_act.add(984594740101022l);
		 kmer_act.add(931204718884809l);
		 kmer_act.add(839448053400340l);
		 kmer_act.add(1042310949184840l);
		 kmer_act.add(1123692687785455l);
		 kmer_act.add(375035996810557l);
		 kmer_act.add(516020959044307l);
		 kmer_act.add(422074513874910l);
		 kmer_act.add(832256480884604l);
		 kmer_act.add(560679239876219l);
		 kmer_act.add(347119155011365l);
		 kmer_act.add(588779852091116l);
		 kmer_act.add(727309458140191l);
		 kmer_act.add(1105109400811479l);
		 kmer_act.add(307723101879673l);
		 kmer_act.add(1089568490121143l);
		 kmer_act.add(1117507405987131l);
		 kmer_act.add(527213608024895l);
		 kmer_act.add(727298934963013l);
		 kmer_act.add(1089158540813548l);
		 kmer_act.add(410480216471732l);
		 kmer_act.add(560851828207438l);
		 kmer_act.add(102620054117933l);
		 kmer_act.add(1046269879136580l);
		 kmer_act.add(1105992399916113l);
		 kmer_act.add(807379796018448l);
		 kmer_act.add(279806244114975l);
		 kmer_act.add(1117071030613949l);
		 kmer_act.add(419970002704274l);
		 kmer_act.add(1050306852811348l);
		 kmer_act.add(1119224976459901l);
		 kmer_act.add(791544076211489l);
		 kmer_act.add(977719370388545l);
		 kmer_act.add(416345460430161l);
		 kmer_act.add(544597239298935l);
		 kmer_act.add(1090584401927924l);
		 kmer_act.add(553980103974473l);
		 kmer_act.add(374244080399605l);
		 kmer_act.add(1090573615157223l);
		 kmer_act.add(1077226109853170l);
		 kmer_act.add(1032027832669463l);
		 kmer_act.add(823527690717522l);
		 kmer_act.add(1054286943482053l);
		 kmer_act.add(562398148657019l);
		 kmer_act.add(938183929334607l);
		 kmer_act.add(557752326913525l);
		 kmer_act.add(1116817052662253l);
		 kmer_act.add(413278378716879l);
		 kmer_act.add(914376491160711l);
		 kmer_act.add(262576713202837l);
		 kmer_act.add(539481934878021l);
		 kmer_act.add(103319594679219l);
		 kmer_act.add(702388035149693l);
		 assertEquals(true,sort_kmer.equals(kmer_act));
		//fail("Not yet implemented");
	}

	public void main(String args[]) {
		testReadsToKmer();
	}
}
